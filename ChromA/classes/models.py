from ..util.HMM import message_passing_multi, message_passing_posterior_state, message_passing_incremental
from ..util.ParamStorage import ParamStorage
from ..classes import data_handle

import multiprocessing
import numpy as np
import logging
import copy
import ray as ray
import os


class BayesianHsmmExperimentMultiProcessing:
    def __init__(self, states, pi_prior, tmat_prior, data=None, length=None, start=None, chrom=None,
                 blacklisted=True, save=False, top_states=None, compute_regions=False, datatype='atac'):

        self.logger = logging.getLogger()
        self.datatype = datatype

        # Data Containers, can be None and updated during Runtime
        self.data = data
        self.length = length
        self.start = start
        self.chrom = chrom
        self.blacklisted = blacklisted
        self.save = save
        self.compute_regions = compute_regions
        self.regions_test = []
        self.peaks = []

        # Data output Containers
        self.annotations = []
        self.annotations_chr = []
        self.annotations_start = []
        self.annotations_length = []

        # Elbo Parameters
        self.elbo = []

        # Specify model parameters (List of State Objects, Initial Prior, Compressed Tmat Prior)
        self.states = states
        self.pi_prior = pi_prior
        self.tmat_prior = tmat_prior
        self.k = np.size(tmat_prior, axis=0)

        # Loop through the states and assign them corresponding index in the transition
        self.s = sum(s.r for s in self.states)
        idx = 0
        if top_states is None:
            self.top_states = None
            for s in self.states:
                r = s.r
                s.idx = idx
                idx += r
                s.save = self.save
        else:
            self.top_states = top_states
            for i_, s in enumerate(self.states):
                s.idx = idx
                idx += s.r
                s.e_log_a = top_states[0].posterior.e_log_a[:, i_]
            idx = 0
            for s in self.top_states:
                s.idx = idx
                idx += s.r

        # Prior Allocation
        self.prior = ParamStorage(K=self.k)
        self.prior.setField('tmat', tmat_prior, dims=('K', 'K'))
        self.prior.setField('pi', pi_prior, dims='K')

        # Posterior Allocation
        self.posterior = ParamStorage(K=self.k)
        self.posterior.setField('tmat', tmat_prior, dims=('K', 'K'))
        self.posterior.setField('pi', pi_prior, dims='K')

    def train(self, filename, iterations, species=None, speciesfile=None, single_chr=None, opt="mo"):

        # ######################################################################################################
        # Get Logger Info
        if hasattr(self.logger.handlers[0], "baseFilename"):
            name = self.logger.handlers[0].baseFilename
        else:
            name = None

        # ######################################################################################################
        # Defining Ray Environment
        processors = multiprocessing.cpu_count()
        if processors > 22:
            if ray.utils.get_system_memory() < 80e9:
                memo = ray.utils.get_system_memory()
                self.logger.info("Recommended Memory > 80GB".format(memo))
            else:
                memo = int(150e9)
        else:
            self.logger.info("Number of Recommended Processors is > 22".format(int(processors)))
            memo = ray.utils.get_system_memory()

        self.logger.info("Running with {0} processors. Size of Plasma Storage {1}".format(int(processors), memo))
        if not ray.is_initialized():
            ray.init(num_cpus=int(processors), include_webui=False)

        # ######################################################################################################
        # Running Regions
        self.logger.info("Training on Regions")
        results = []
        chromosome = [Trainer.remote(-1, filename, species, speciesfile, self.blacklisted, self.states, self.prior,
                                     self.top_states, logger=logging.getLogger().getEffectiveLevel(),
                                     log_file=name, datatype=self.datatype)]
        results.append(chromosome[0].train.remote(iterations=iterations, msg="Th17 Regions: "))

        # Collect Results
        res, states = ray.get(results[0])
        for l_ in np.arange(len(res[0])):
            self.annotations.append(res[0][l_])
            self.annotations_chr.append(res[1][l_])
            self.annotations_start.append(res[2][l_])
            self.annotations_length.append(res[3][l_])
        posterior = ray.get(chromosome[0].get_posterior.remote())
        self.elbo = ray.get(chromosome[0].get_elbo.remote())

        # Validate Results
        self.validate_regions()

        # ######################################################################################################
        # Running Chromosomes
        if not self.compute_regions:
            # Prepare data structures and carry over states from Regions
            self.annotations = []
            self.annotations_chr = []
            self.annotations_start = []
            for i_ in np.arange(len(self.states)):
                if isinstance(self.states[0], type(states[0])):
                    self.states[i_] = copy.deepcopy(states[i_])
                else:
                    self.states[i_] = copy.deepcopy(states[0][i_])
                self.states[i_].prior = self.states[i_].posterior

            # Prune chromosomes
            chr_list = data_handle.validate_chr(filename, species, speciesfile, chr_list=single_chr)

            # Run Training in parallel
            while len(chr_list) > 0:
                results = []
                chromosome = []
                num_task = np.min([processors, len(chr_list)])
                for i_ in np.arange(num_task):
                    chr_ = chr_list[0]
                    self.logger.info("chr{}: Submitting job to Queue".format(chr_))
                    chromosome.append(Trainer.remote(chr_, filename, species, self.blacklisted, self.states, self.prior,
                                                     self.top_states, pi=posterior.pi, tmat=posterior.tmat,
                                                     logger=logging.getLogger().getEffectiveLevel(), log_file=name,
                                                     datatype=self.datatype))
                    results.append(chromosome[i_].train.remote(iterations=iterations, msg="chr{}: ".format(chr_)))
                    chr_list.remove(chr_)

                unfinished = results
                while len(unfinished) > 0:
                    finished, unfinished = ray.wait(unfinished)

                    for r_ in finished:
                        res, _ = ray.get(r_)
                        for l_ in np.arange(len(res[0])):
                            self.annotations.append(res[0][l_])
                            self.annotations_chr.append(res[1][l_])
                            self.annotations_start.append(res[2][l_])
                            self.annotations_length.append(res[3][l_])

        # Clean Ray
        ray.shutdown()

    def save_bedfile(self, path, name=None, thres=0.05, ext=100, merge=500, filterpeaks=0):
        # Format Filename
        if name is None:
            name = "region_"
        if path is "":
            path = os.getcwd()

        # Format Regions
        chromm = list()
        for i_ in np.arange(len(self.annotations_chr)):
            chromm.append(self.annotations_chr[i_][3:])

        # Save Bed-File All States
        regs = []
        for l_ in np.arange(len(self.annotations)):
            if self.annotations[0].shape[1] > 2:
                regs.append(self.annotations[l_][:, 1:].sum(axis=1))
            else:
                regs.append(self.annotations[l_][:, 1])
        peaks = data_handle.bed_result(os.path.join(path, name) + '_allpeaks.bed',
                                       regs, self.annotations_start, chromm, threshold=thres,
                                       bedext=ext, bedmerge=merge, filterpeaks=filterpeaks)
        self.peaks = peaks

        self.logger.info("Saved Bed File. ")

    def posterior_state(self, fname=os.getcwd()):
        chr_list = np.unique(self.annotations_chr)
        out = []
        for i_, c_ in enumerate(chr_list):
            out_s = []
            out_st = []
            out_l = []
            idx_c = np.where(np.array(self.annotations_chr) == c_)[0]
            for idx_reg in idx_c:
                out_s.append(np.float16(self.annotations[idx_reg][:, 1]))
                out_st.append(self.annotations_start[idx_reg])
                out_l.append(self.annotations_length[idx_reg])
            out.append([c_, out_st, out_l, out_s])
        path, file = os.path.split(fname)
        np.save(os.path.join(path, file + '.posterior_state.npy'), np.array(out))

    def validate_regions(self):
        # Get Logger
        logger = logging.getLogger('metrics')
        logger.info("METRICS ON REGIONS.")

        # Compute Metrics
        total_length = np.sum(self.annotations_length)
        n_states = self.annotations[0].shape[1]
        if n_states > 2:
            state1 = np.zeros(n_states - 1)
            for s_ in np.arange(n_states - 1):
                for r_ in self.annotations:
                    state1[s_] += r_[:, s_ + 1].sum()
            perc = state1 / total_length
            if perc[1] > 0.25:
                logger.info("ChromA annotated a high fraction of the regions. Verify Signal to Noise.")
                self.regions_test = False
            elif perc[1] < 0.01:
                logger.info("ChromA did not use High Signal state. Verify Signal to Noise.")
                self.regions_test = False
            else:
                logger.info("ChromA passed regions annotation test.")
                self.regions_test = True

        else:
            state1 = 0
            for r_ in self.annotations:
                state1 += r_[:, 1].sum()
            perc = state1 / total_length
            if perc > 0.25:
                logger.info("ChromA annotated a high fraction of the regions. Verify Signal to Noise.")
                self.regions_test = False
            else:
                logger.info("ChromA passed regions annotation test.")
                self.regions_test = True
    # #############################################################################################


@ray.remote(num_cpus=1)
class Trainer(object):
    def __init__(self, chr_, filename, species, speciesfile, blacklisted, states, prior,
                 top_states=None, pi=None, tmat=None, logger=None, log_file=None, datatype='atac'):
        # Init Logging Module
        if logger is None:
            data_handle.build_logger('0', filename=log_file, supress=True)
        elif logger == 0:
            data_handle.build_logger('0', filename=log_file, supress=True)
        elif logger == 10:
            data_handle.build_logger('2', filename=log_file, supress=True)
        elif logger == 20:
            data_handle.build_logger('1', filename=log_file, supress=True)
        self.logger = logging.getLogger()

        # Formatting States
        # State0 is used to fit individual Models.
        # State is used to fit Multi-experiment Models.
        self.states0 = copy.deepcopy(states)
        self.states = copy.deepcopy(states)
        for s_ in self.states:
            s_.mo = list()
        for s_ in self.states0:
            s_.mo = list()
        self.k = prior.tmat.shape[0]
        self.s = sum(s.r for s in self.states0)
        if top_states is not None:
            self.top_states = top_states

        # Getting Next Chromosome Data
        self.n_exp = len(filename)
        if datatype == 'atac':
            dnase = False
        elif datatype == 'dnase':
            dnase = True
        else:
            dnase = False

        if chr_ == -1:
            self.logger.info("Regions: Fetching Data")
            data, length, start, chrom = data_handle.regions_th17(filename=filename,
                                                                  species=species, dnase=dnase)
        else:
            chrom_str = "chr" + chr_.__str__()
            self.logger.info(chrom_str + ": Fetching Data")
            data, length, start, chrom = data_handle.regions_chr(filename=filename, chromosome=chrom_str,
                                                                 species=species, specfile=speciesfile,
                                                                 blacklisted=blacklisted, dnase=dnase)

        self.data = data
        self.length = length
        self.start = start
        self.chrom = chrom

        # Formatting Prior, Posterior
        self.prior = prior
        if self.n_exp == 1:
            self.posterior = ParamStorage(K=self.k)
            if tmat is None:
                tmat = self.prior.tmat
            self.posterior.setField('tmat', tmat, dims=('K', 'K'))
            if pi is None:
                pi = self.prior.pi
            self.posterior.setField('pi', pi, dims='K')
        elif self.n_exp > 1:
            self.posterior = ParamStorage(K=self.k, N=self.data.shape[0])
            self.posterior.setField('tmat', self.prior.tmat, dims=('K', 'K'))
            self.posterior.setField('pi', self.prior.pi, dims='K')
            self.posterior.setField('s_s', np.zeros((self.data.shape[0], self.k)), dims=('N', 'K'))
        self.elbo = []
        self.elbo_interrupted = 0

    def get_posterior(self):
        return self.posterior

    def get_elbo(self):
        return self.elbo[:self.elbo_interrupted, :]

    def train(self, iterations=20, msg=""):
        if len(self.data) > 0:
            if self.n_exp == 1:
                self.logger.info(msg + "Running Single File Routine")
                return self.train_single(iterations=iterations, msg=msg)
            elif self.n_exp > 1:
                self.logger.info(msg + "Running Multiple Files Routine")
                return self.train_multiple(iterations=iterations, msg=msg)
        else:
            return [[], [], [], []], []

    def train_single(self, iterations=20, msg=""):

        self.elbo = np.zeros([iterations, 2])
        for s in self.states0:
            s.it = 1

        # Iterations
        self.logger.info(msg + "Fitting model")
        self.elbo_interrupted = iterations
        for i_ in np.arange(0, iterations):
            self.vb_update()
            self.elbo[i_] = self.calc_elbo()
            self.print_iteration(elbo=self.elbo[i_], iteration=i_, message=msg)
            if np.abs(self.elbo[i_].sum() - self.elbo[i_ - 1].sum()) < 1:
                self.logger.debug(msg + "Finished by ELBO criteria.")
                self.elbo_interrupted = i_
                break

        # Formatting the Output
        self.logger.info(msg + "Calculate S")
        s_s = message_passing_posterior_state(self.posterior.pi, self.posterior.tmat, self.states0, self.s, self.k,
                                              self.length, data=self.data)
        self.logger.info(msg + "Formatting Output, Saving {} Regions.".format(len(self.length)))
        regions = []
        regions_start = []
        regions_chr = []
        regions_length = []
        count = 0
        for i_ in np.arange(len(self.length)):
            regions.append(s_s[count: count + self.length[i_]])
            regions_chr.append(self.chrom[i_])
            regions_start.append(self.start[i_])
            regions_length.append(self.length[i_])
            count += self.length[i_]
        self.logger.info(msg + "Output Stored.")

        return [regions, regions_chr, regions_start, regions_length], self.states0

    def train_multiple(self, iterations=5, msg=""):
        iterations = 10
        self.elbo = np.zeros([iterations, 2])
        for s in self.states:
            s.it = 1

        # Fit each Datasets Preliminary and Update Parameters
        states = list()
        for i_ in np.arange(self.n_exp):
            self.logger.info(msg + "Fitting Individual models. Model {}.".format(i_))
            [self.vb_update(exp=i_) for _ in np.arange(5)]
            self.posterior.s_s += message_passing_posterior_state(self.posterior.pi, self.posterior.tmat, self.states0,
                                                                  self.s, self.k, self.length,
                                                                  data=self.data[:, i_][:, None]) / self.n_exp
            states.append(copy.deepcopy(self.states0))
        self.states = states

        # Iterations
        self.logger.info(msg + "Fitting Consensus model")
        for i_ in np.arange(0, iterations):
            self.vb_update_multi()
            self.elbo[i_] = self.calc_elbo(state_flag=1)
            self.print_iteration(elbo=self.elbo[i_], iteration=i_, message=msg)
            # if np.abs(self.elbo[i_].sum() - self.elbo[i_ - 1].sum()) < 1:
            #     print(msg + "Finished by ELBO criteria.")
            #     break

        # Formatting the Output
        s_s = self.posterior.s_s
        self.logger.info(msg + "Formatting Output, Saving {} Regions.".format(len(self.length)))
        regions = []
        regions_start = []
        regions_chr = []
        regions_length = []
        count = 0
        for i_ in np.arange(len(self.length)):
            regions.append(s_s[count: count + self.length[i_]])
            regions_chr.append(self.chrom[i_])
            regions_start.append(self.start[i_])
            regions_length.append(self.length[i_])
            count += self.length[i_]
        self.logger.info(msg + "Output Stored.")

        return [regions, regions_chr, regions_start, regions_length], self.states

    def vb_update(self, exp=0):

        # MESSAGE PASSING: THIS VERSION ONLY TRAINS MO
        lmarg_pr_obs, ss_pi, ss_tmat = message_passing_incremental(self.posterior.pi, self.posterior.tmat,
                                                                   self.states0, self.s, self.k, self.length,
                                                                   data=self.data[:, exp], opt="mo")

        # PARAMETER UPDATES
        self.posterior.setField('lmarg_pr_obs', np.squeeze(lmarg_pr_obs), dims=None)
        self.posterior.setField('pi', self.prior.pi + ss_pi, dims='K')
        self.posterior.setField('tmat', self.prior.tmat + ss_tmat, dims=('K', 'K'))

    def vb_update_multi(self):

        p = self.posterior
        lmarg_pr_obs, ss_pi, ss_tmat = \
            message_passing_multi(p.pi, p.tmat, self.states, self.top_states, p.s_s, self.n_exp,
                                  self.s, self.k, self.length, data=self.data)

        # PARAMETER UPDATES
        self.posterior.setField('lmarg_pr_obs', lmarg_pr_obs[0, 0], dims=None)
        self.posterior.setField('pi', self.prior.pi + ss_pi, dims='K')
        self.posterior.setField('tmat', self.prior.tmat + ss_tmat, dims=('K', 'K'))

        # Update Parameters
        ts = self.top_states
        [s.update_parameters_ss() for s in self.top_states]
        [[s.update_parameters_ss(ts[i].posterior.e_log_a[i, :]) for i, s in enumerate(t)] for t in self.states]

    def calc_elbo(self, state_flag=0):
        kl_term = 0
        for s_ in self.states0:
            kl_term += s_.kl_term()

        if state_flag == 1:
            for exp in self.states:
                for s_ in exp:
                    kl_term += s_.kl_term()

        return self.posterior.lmarg_pr_obs, kl_term

    @staticmethod
    def print_iteration(message="", elbo=None, iteration=None):
        logger = logging.getLogger()
        logger.debug(message + 'iteration:' + iteration.__str__() + '    ELBO:' + elbo.__str__())
