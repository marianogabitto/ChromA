from ..util.HMM import message_passing_multi, message_passing_posterior_state, message_passing_incremental
from ..util.ParamStorage import ParamStorage
from ..classes import data_handle, states

import multiprocessing
import numpy as np
import logging
import copy
import ray
import os

import matplotlib
gui_env = ['Agg', 'TKAgg', 'GTKAgg', 'Qt4Agg', 'WXAgg']
for gui in gui_env:
    try:
        matplotlib.use(gui, warn=False, force=True)
        from matplotlib import pyplot as plt
        break
    except:
        continue
from matplotlib.backends.backend_pdf import PdfPages
eps = 1e-9


class BayesianHsmmExperimentMultiProcessing:
    def __init__(self, states, pi_prior, tmat_prior, data=None, length=None, start=None, chrom=None,
                 blacklisted=True, save=False, top_states=None, compute_regions=False):

        self.logger = logging.getLogger()
        self.species = None

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

    def train(self, filename, iterations, species=None, single_chr=None, opt="mo"):

        # ##################################
        # Get Logger Info
        if hasattr(self.logger.handlers[0], "baseFilename"):
            name = self.logger.handlers[0].baseFilename
        else:
            name = None

        # ##################################
        # Defining Ray Environment
        processors = multiprocessing.cpu_count()
        if processors > 22:
            if ray.utils.get_system_memory() < 80e9:
                memo = ray.utils.get_system_memory()
                self.logger.info("Recommended Memory > 80GB".format(memo))
            else:
                memo = int(100e9)
        else:
            self.logger.info("Number of Recommended Processors is > 22".format(int(processors)))
            memo = ray.utils.get_system_memory()

        self.logger.info("Running with {0} processors. Size of Plasma Storage {1}".format(int(processors), memo))
        if not ray.is_initialized():
            ray.init(num_cpus=int(processors), object_store_memory=memo, include_webui=False)

        # ##################################
        # Running Regions
        self.logger.info("Training on Regions")
        results = []
        chromosome = [Trainer.remote(-1, filename, species, self.blacklisted, self.states, self.prior,
                                     self.top_states, logger=logging.getLogger().getEffectiveLevel(), log_file=name)]
        results.append(chromosome[0].train.remote(iterations=50, msg="Th17 Regions: "))

        # Collect Results
        res, stts = ray.get(results[0])
        for l_ in np.arange(len(res[0])):
            self.annotations.append(res[0][l_])
            self.annotations_chr.append(res[1][l_])
            self.annotations_start.append(res[2][l_])
            self.annotations_length.append(res[3][l_])

        posterior = ray.get(chromosome[0].get_posterior.remote())
        self.elbo = ray.get(chromosome[0].get_elbo.remote())

        # Validate Results
        self.validate_regions()

        # ##################################
        # Running Chromosomes
        if not self.compute_regions:
            self.annotations = []
            self.annotations_chr = []
            self.annotations_start = []
            for i_ in np.arange(len(self.states)):
                # if single File. states is list of states
                if type(self.states[0]) == type(stts[0]):
                    self.states[i_] = copy.deepcopy(stts[i_])
                # if multiple Files. states is list of list of states
                else:
                    self.states[i_] = copy.deepcopy(stts[0][i_])
                self.states[i_].prior = self.states[i_].posterior

            chr_list = []
            if species is None or single_chr is not None:
                if single_chr is not None:
                    chr_list = single_chr
                else:
                    self.logger.error('Species and single_chr cannot be None at the same time.')
            elif species == 'mouse':
                self.species = 'mouse'
                chr_list = list(np.arange(1, 20))
                chr_list.append('X')
                chr_list.append('Y')
                self.logger.info('Running on mouse genome. 19 Chroms')
            elif species == 'human':
                self.species = 'human'
                chr_list = list(np.arange(1, 22))
                chr_list.append('X')
                chr_list.append('Y')
                self.logger.info('Running on human genome. 22 Chroms')
            elif species == 'fly':
                self.species = 'fly'
                chr_list = ['2L', '2R', '3L', '3R', '4', 'X', 'Y']
                self.logger.info('Running on fly genome. 7 Chroms')

            # Run Training in parallel
            while len(chr_list) > 0:
                results = []
                chromosome = []
                for i_ in np.arange(np.min([processors, len(chr_list)])):
                    chr_ = chr_list[0]
                    self.logger.info("chr{}: Submitting job to Queue".format(chr_))
                    chromosome.append(Trainer.remote(chr_, filename, species, self.blacklisted, self.states, self.prior,
                                                     self.top_states, pi=posterior.pi, tmat=posterior.tmat,
                                                     logger=logging.getLogger().getEffectiveLevel(), log_file=name))
                    results.append(chromosome[i_].train.remote(iterations=iterations, msg="chr{}: ".format(chr_)))
                    chr_list.remove(chr_)

                # Collect Results
                for r_ in reversed(results):
                    res, _ = ray.get(r_)
                    for l_ in np.arange(len(res[0])):
                        self.annotations.append(res[0][l_])
                        self.annotations_chr.append(res[1][l_])
                        self.annotations_start.append(res[2][l_])
                        self.annotations_length.append(res[3][l_])

    def save_bed(self, path, name=None):
        # Check that we have already computed the annotations and return the depth of the result.
        if len(self.annotations) == 0:
            print("Compute Annotations before writing to file.")
            return
        else:
            depth = lambda L: isinstance(L, list) and max(map(depth, L)) + 1
            d_annot = depth(self.annotations)

        # Format Filename
        if name is None:
            name = "region_"
        if path is "":
            path = os.getcwd()

        # Identify Single or Multiple Exp. Submit Bed File to be written.
        # Format Chromosome Names
        chromm = list()
        for i_ in np.arange(len(self.annotations_chr)):
            chromm.append(self.annotations_chr[i_][3:])

        if d_annot == 1:
            # Do single trial
            fname = os.path.join(path, name) + '_peaks.bed'
            self.save_bedfile(filename=fname,  annotations=self.annotations,
                              chrom=chromm, start=self.annotations_start)
        else:
            # Do multiple trials
            # First write Consensus
            peaks = list()
            fname = os.path.join(path, name) + '_cons_peaks.bed'
            peaks.append(self.save_bedfile(filename=fname, annotations=self.annotations, chrom=chromm,
                                           start=self.annotations_start, exp=0))
            # Then Experiments
            for i_ in np.arange(1, len(self.annotations[0])):
                fname = os.path.join(path, name) + '_exp_' + i_.__str__() + '_peaks.bed'
                peaks.append(self.save_bedfile(filename=fname, annotations=self.annotations, chrom=chromm,
                                               start=self.annotations_start, exp=i_))

            # Finally, Consensus Plus Experiments
            fname = os.path.join(path, name) + '_agg_peaks.bed'
            self.aggregate_peaks(peaks, filename=fname)

        self.logger.info("Saved Bed File. ")

    @staticmethod
    def save_bedfile(filename=None, annotations=None, chrom=None, start=None, exp=None):
        if (filename is None) or (annotations is None) or (chrom is None) or (start is None):
            print("Empty Annotations.")
            return

        regs = []
        if exp is None:
            for l_ in np.arange(len(annotations)):
                regs.append(annotations[l_][:, 1])
        else:
            for l_ in np.arange(len(annotations)):
                regs.append(annotations[l_][exp][:, 1])
        peaks = data_handle.bed_result(filename, regs, start, chrom, threshold=0.05)

        return peaks

    @staticmethod
    def aggregate_peaks(peaks, filename):
        # Parse Information
        # peaks[cons, exp1, exp2, ...][chr, regs]
        chr_l = [peaks[0][0]]
        exps = [peaks[0][1]]
        for e_ in np.arange(1, len(peaks)):
            chr_l.append(peaks[e_][0])
            exps.append(peaks[e_][1])

        # Intersect peaks
        reg_out = exps[0]
        chr_out = chr_l[0]

        def get_overlap(a, b):
            # [a] = [2], [b] = [n, 2]
            lower = np.min([a[1] * np.ones_like(b[:, 1]), b[:, 1]], axis=0)
            upper = np.max([a[0] * np.ones_like(b[:, 0]), b[:, 0]], axis=0)
            return np.max([np.zeros_like(lower), lower - upper], axis=0)

        for ne_, e_ in enumerate(exps[1:]):
            r = reg_out.tolist()
            for npeak_, peak_ in enumerate(e_):
                if np.all(get_overlap(peak_, reg_out) == 0):
                    r.append(peak_)
                    chr_out.append(chr_l[ne_ + 1][npeak_])
            reg_out = np.array(r)

        # Sort
        index = np.argsort(reg_out[:, 0])
        idx_out = []
        chr_out = np.array(chr_out)
        for c_ in np.sort(np.unique(chr_out)):
            idx_out.append(index[chr_out == c_])
        idx_out = np.concatenate(idx_out)
        chr_out = chr_out[idx_out]
        reg_out = reg_out[idx_out]

        # Write Output
        data_handle.write_bed(filename, data=np.array(chr_out), start=reg_out[:, 0], end=reg_out[:, 1])

    def validate_regions(self):
        # Get Logger
        logger = logging.getLogger('metrics')
        logger.info("METRICS ON REGIONS.")

        # Compute Metric on Single Trial or Consensus State
        depth = lambda l: isinstance(l, list) and max(map(depth, l)) + 1
        d_annot = depth(self.annotations)
        if d_annot == 1:
            # Annotations Single Trial
            # self.annotations[# regions][length x states]
            annot = self.annotations
        else:
            # Annotations Multiple Trial
            # self.annotations[# regions][Cons, exp1, exp2..][length x states]
            annot = list()
            for r_ in self.annotations:
                annot.append(r_[0])

        # Compute Metrics
        total_length = np.sum(self.annotations_length)
        n_states = annot[0].shape[1]
        if n_states > 2:
            state1 = np.zeros(n_states - 1)
            for s_ in np.arange(n_states - 1):
                for r_ in annot:
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
            for r_ in annot:
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
    def __init__(self, chr_, filename, species, blacklisted, states_, prior,
                 top_states=None, pi=None, tmat=None, logger=None, log_file=None):
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
        self.states0 = copy.deepcopy(states_)
        self.states = copy.deepcopy(states_)
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
        if chr_ == -1:
            self.logger.info("Regions: Fetching Data")
            data, length, start, chrom = data_handle.regions_th17(filename=filename, species=species)
        else:
            chrom_str = "chr" + chr_.__str__()
            self.logger.info(chrom_str + ": Fetching Data")
            data, length, start, chrom = data_handle.regions_chr(filename=filename, chromosome=chrom_str,
                                                                 species=species, blacklisted=blacklisted)

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

    def out(self, s_s, msg):
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

        return [regions, regions_chr, regions_start, regions_length]

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
        output = self.out(s_s=message_passing_posterior_state(self.posterior.pi, self.posterior.tmat, self.states0,
                                                              self.s, self.k, self.length, data=self.data), msg=msg)

        return output, self.states0

    def train_multiple(self, iterations=5, msg=""):
        iterations = 10
        self.elbo = np.zeros([iterations, 2])
        for s in self.states:
            s.it = 1

        # Fit each Datasets Preliminary and Update Parameters
        states = list()
        post_s = list()
        for i_ in np.arange(self.n_exp):
            self.logger.info(msg + "Fitting Individual models. Model {}.".format(i_))
            [self.vb_update(exp=i_) for _ in np.arange(5)]
            temp_s = message_passing_posterior_state(self.posterior.pi, self.posterior.tmat, self.states0,
                                                     self.s, self.k, self.length,
                                                     data=self.data[:, i_][:, None]) / self.n_exp
            self.posterior.s_s += temp_s
            states.append(copy.deepcopy(self.states0))
            post_s.append(self.out(s_s=temp_s, msg='Exp:' + i_.__str__() + ' ' + msg))
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
        # [consensus] = [annotations, chr, st, length]
        consensus_s = self.out(s_s=self.posterior.s_s, msg=msg)
        # [output] = [# regions][consensus, exp1, exp2, ...]
        output = list()
        for reg_ in np.arange(len(consensus_s[0])):
            temp = list()
            temp.append(consensus_s[0][reg_])
            for exp_ in np.arange(len(post_s)):
                temp.append(post_s[exp_][0][reg_])
            output.append(temp)

        # Returning
        # [annotations[# regs][consensus, exp1, exp2, ...], chr[# regs], st[# regs], length[# regs]]
        return [output, consensus_s[1], consensus_s[2], consensus_s[3]], self.states

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
