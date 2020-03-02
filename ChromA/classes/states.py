from ..util.ParamStorage import ParamStorage

from scipy.special import psi, gammaln
from scipy.stats import binom
import numpy as np
import logging


def build_states(typ='atac', filename=None, r=None):
    logger = logging.getLogger()

    class BedOptions(object):
        def __init__(self, thres=0.05, ext=100, merge=500, filterpeaks=50):
            self.thres = thres
            self.ext = ext
            self.merge = merge
            self.filterpeaks = filterpeaks

    if r is not None:
        logger.info("Building States r={}".format(r))

        pi_prior = np.array([1000, 1])
        tmat_prior = np.array([[0, 1],
                               [1, 0]])
        state_list = []
        closed_state = NegativeBinomialGEO(r=int(r[0]), p=r[2], cut0=1, cut1=50)
        open_state = NegativeBinomialGEO(r=int(r[1]), p=r[3], cut0=20, cut1=10)
        state_list.append(closed_state)
        state_list.append(open_state)

        top_states = None
        if len(filename) > 1:
            top_tmat0 = np.array([[100, 1],
                                  [1, 100]])
            top_states = list()
            top_closed_state = TopStateNegativeBinomial(r=int(r[0]), p=r[2], toptmat0=top_tmat0, order=0)
            top_open_state = TopStateNegativeBinomial(r=int(r[1]), p=r[3], toptmat0=top_tmat0, order=1)
            top_states.append(top_closed_state)
            top_states.append(top_open_state)

        bedopts = BedOptions(thres=r[4], ext=r[5], merge=r[6], filterpeaks=r[7])

    elif typ == "atac":
        r1 = 3
        r2 = 2
        logger.info("Running with 2 States: r={}, r={}".format(r1, r2))
        pi_prior = np.array([1000, 1])
        tmat_prior = np.array([[0, 1],
                               [1, 0]])
        state_list = []
        closed_state = NegativeBinomialGEO(r=r1, p=1e-4, cut0=1, cut1=50)
        open_state = NegativeBinomialGEO(r=r2, p=1e-4, cut0=20, cut1=10)
        state_list.append(closed_state)
        state_list.append(open_state)

        top_states = None
        if len(filename) > 1:
            top_tmat0 = np.array([[100, 1],
                                  [1, 100]])

            top_states = list()
            top_closed_state = TopStateNegativeBinomial(r=r1, p=1e-4, toptmat0=top_tmat0, order=0)
            top_open_state = TopStateNegativeBinomial(r=r2, p=1e-4, toptmat0=top_tmat0, order=1)
            top_states.append(top_closed_state)
            top_states.append(top_open_state)

        bedopts = BedOptions(thres=0.05, ext=100, merge=500, filterpeaks=0)

    elif typ == "dnase":
        r1 = 3
        r2 = 2
        logger.info("Running with 2 States: r={}, r={}".format(r1, r2))
        pi_prior = np.array([1000, 1])
        tmat_prior = np.array([[0, 1],
                               [1, 0]])
        state_list = []
        closed_state = NegativeBinomialGEO(r=r1, p=1.5e-4, cut0=1, cut1=50)
        open_state = NegativeBinomialGEO(r=r2, p=1.5e-4, cut0=20, cut1=10)
        state_list.append(closed_state)
        state_list.append(open_state)

        top_states = None
        if len(filename) > 1:
            top_tmat0 = np.array([[100, 1],
                                  [1, 100]])

            top_states = list()
            top_closed_state = TopStateNegativeBinomial(r=r1, p=1.5e-4, toptmat0=top_tmat0, order=0)
            top_open_state = TopStateNegativeBinomial(r=r2, p=1.5e-4, toptmat0=top_tmat0, order=1)
            top_states.append(top_closed_state)
            top_states.append(top_open_state)

        bedopts = BedOptions(thres=0.05, ext=0, merge=500, filterpeaks=75)

    else:
        raise SystemError

    return pi_prior, tmat_prior, state_list, top_states, bedopts


class NegativeBinomialGEO:
    def __init__(self, r, cut0, cut1, p=None, s_count=None, f_count=None):

        # Saving Variables
        self.save = False
        self.e_log_a = []

        # Empty container sufficient statistics
        self.mo = list()
        self.ss = list()
        self.ss.append(np.zeros(2))
        self.ss.append(np.zeros(2))

        # SO
        self.rho = 1.
        self.it = 1.
        self.rhodelay = 20.
        self.rhoexp = 0.55

        # Number of expanded states
        self.n = r
        self.r = r

        # Negative Binomial Parameter Count (successes and failures)
        if p is None:
            self.s_count = s_count
            self.f_count = f_count
            self.update = True
        else:
            self.update = False

        # Build the block prior matrix
        self.block_prior = np.zeros((r, r))

        # Define the prior bag of parameters
        self.prior = ParamStorage()
        if self.update is True:
            self.prior.setField('p', s_count, dims=None)
            self.prior.setField('f', f_count, dims=None)
        self.prior.setField('emit0', cut0, dims=None)
        self.prior.setField('emit1', cut1, dims=None)

        # Define the posterior bag of parameters
        self.posterior = ParamStorage()
        if self.update is True:
            self.posterior.setField('p', s_count, dims=None)
            self.posterior.setField('f', f_count, dims=None)
        self.posterior.setField('emit0', cut0, dims=None)
        self.posterior.setField('emit1', cut1, dims=None)

        # self.idx remembers this states location in an expanded transition matrix, it is set
        # by the HsMM model object when it builds an expanded transition matrix
        self.idx = None

        # Data Storage
        if self.update is True:
            self.E_log_p = psi(self.posterior.p) - psi(self.posterior.p + self.posterior.f)
            self.E_log_f = psi(self.posterior.f) - psi(self.posterior.p + self.posterior.f)
        else:
            self.E_log_p = np.log(p)
            self.E_log_f = np.log(1. - p)
        self.leaving_prob = self.E_log_p
        self.enter_prob = np.ones(self.n)
        e_f = np.exp(self.E_log_f) / (np.exp(self.E_log_p) + np.exp(self.E_log_f))
        self.enter_prob = np.log(binom.pmf(np.arange(self.n)[::-1], self.n - 1, e_f))

    def log_likelihood(self, obs=None, s=None):

        if obs is not None:
            e_success = self.posterior.emit0
            e_fail = self.posterior.emit1

            l_e = np.zeros(2)
            psi_norm = psi(e_success + e_fail)
            l_e[0] = psi(e_success) - psi_norm
            l_e[1] = psi(e_fail) - psi_norm

            cut = obs[:, 0]
            log_likelihood_vec = (cut * l_e[0] + l_e[1])
        else:
            log_likelihood_vec = np.zeros((s.shape[0]))

        if s is not None:
            if self.e_log_a is not None:
                log_likelihood_vec += np.dot(s, self.e_log_a)

        return log_likelihood_vec

    def kl_term(self):

        psi_norm = psi(self.posterior.emit0 + self.posterior.emit1)
        expt_terms = (self.prior.emit1 - self.posterior.emit1) * (psi(self.posterior.emit1) - psi_norm) + (
                self.prior.emit0 - self.posterior.emit0) * (psi(self.posterior.emit0) - psi_norm)

        post_kl = gammaln(self.posterior.emit1 + self.posterior.emit0) - gammaln(self.posterior.emit1) - gammaln(
            self.posterior.emit0)
        prior_kl = gammaln(self.prior.emit1 + self.prior.emit0) - gammaln(self.prior.emit1) - gammaln(self.prior.emit0)

        return expt_terms + prior_kl - post_kl

    def update_parameters_ss(self, e_log_a=None):

        # Re-casts Parameters
        ss_transitions = self.ss[0]
        ss_emissions = self.ss[1]

        if self.update is True:
            self.posterior.setField('f', self.prior.f + ss_transitions[0], dims='')
            self.posterior.setField('p', self.prior.p + ss_transitions[1], dims='')
            self.E_log_p = psi(self.posterior.p) - psi(self.posterior.p + self.posterior.f)
            self.E_log_f = psi(self.posterior.f) - psi(self.posterior.p + self.posterior.f)
            self.leaving_prob = self.E_log_p
            self.enter_prob = np.log(binom.pmf(np.arange(self.r)[::-1], self.r - 1, self.posterior.f))
            e_f = self.posterior.f / (self.posterior.p + self.posterior.f)
            self.enter_prob = np.log(binom.pmf(np.arange(self.n)[::-1], self.n - 1, e_f))

        # Update Parameters emissions
        self.posterior.setField('emit0', self.prior.emit0 + ss_emissions[0], dims=None)
        self.posterior.setField('emit1', self.prior.emit1 + ss_emissions[1], dims=None)
        if e_log_a is not None:
            self.e_log_a = e_log_a

        if self.save:
            emit = [self.posterior.emit0, self.posterior.emit1]
            np.save('emit_params{0}_{1}.npy'.format(self.idx, self.it), emit)
            self.it += 1

    def update_parameters_so(self, e_log_a=None):

        if self.it < 500:
            # Re-casts Parameters
            ss_transitions = self.ss[0]
            ss_emissions = self.ss[1]

            if self.update is True:
                p = self.prior.p + ss_transitions[0]
                f = self.prior.f + ss_transitions[1]
                new_p = p * self.rho + (1. - self.rho) * self.posterior.p
                new_f = f * self.rho + (1. - self.rho) * self.posterior.f
                self.posterior.setField('p', new_p, dims='')
                self.posterior.setField('f', new_f, dims='')

                self.E_log_p = psi(self.posterior.p) - psi(self.posterior.p + self.posterior.f)
                self.E_log_f = psi(self.posterior.f) - psi(self.posterior.p + self.posterior.f)
                self.leaving_prob = self.E_log_p
                self.enter_prob = np.log(binom.pmf(np.arange(self.r)[::-1], self.r - 1, self.posterior.f))

            # Update Parameters emissions
            emit0 = self.prior.emit0 + ss_emissions[0]
            emit1 = self.prior.emit1 + ss_emissions[1]
            emit0_new = emit0 * self.rho + (1. - self.rho) * self.posterior.emit0
            emit1_new = emit1 * self.rho + (1. - self.rho) * self.posterior.emit1
            self.posterior.setField('emit0', emit0_new, dims=None)
            self.posterior.setField('emit1', emit1_new, dims=None)
            if e_log_a is not None:
                self.e_log_a = e_log_a

            # Update Rho
            self.rho = (1 + self.it + self.rhodelay) ** (-1.0 * self.rhoexp)
            if self.save:
                emit = [self.posterior.emit0, self.posterior.emit1]
                np.save('emit_paramsSO{0}_{1}.npy'.format(self.idx, self.it), emit)
            self.it += 1

    def clear_ss(self):
        self.ss = list()
        self.ss.append(np.zeros(2))
        self.ss.append(np.zeros(2))

    def prepare_ss(self, opt, number):

        if opt == "batch":
            self.ss = list()
            self.ss.append(np.zeros(2))
            self.ss.append(np.zeros(2))
        elif opt == "so":
            self.ss = list()
            self.ss.append(np.zeros(2))
            self.ss.append(np.zeros(2))
            self.rho = (1. + self.it + self.rhodelay) ** (- self.rhoexp)
        elif opt == "mo":
            if len(self.mo) == 0:
                self.mo.append(np.zeros((number, 2)))
                self.mo.append(np.zeros((number, 2)))

    def increase_ss(self, ss=None, s=None, obs=None, suf_stat=None, amp=1.):

        if self.update:
            if ss is not None:
                # Update parameters p and f (using epsilons)
                p_ss = f_ss = 0
                for i in range(0, self.n):
                    f_ss += np.sum(ss[self.idx + i, self.idx + i])
                    p_ss += np.sum(ss[self.idx + i, :]) - np.sum(ss[self.idx + i, self.idx + i])

                index = np.setdiff1d(np.arange(ss.shape[0]), np.arange(self.idx, self.idx + self.r))
                for i in np.arange(self.r):
                    f_ss += np.sum((self.r - i - 1) * ss[index, self.idx + i])
                    p_ss += np.sum(i * ss[index, self.idx + i])

                # p_ss /= self.n

                self.ss[0][0] += amp * p_ss
                self.ss[0][1] += amp * f_ss

        if suf_stat is None:
            if (s is not None) and (obs is not None):
                # Update SS Emmissions
                self.ss[1][0] += amp * np.dot(obs[:, 0].T, s)
                self.ss[1][1] += amp * np.sum(s)
        else:
            self.ss[1][0] += amp * suf_stat[0]
            self.ss[1][1] += amp * suf_stat[1]

    def replace_ss(self, ss=None, s=None, obs=None, number=None):

        # Update SS Transitions
        if self.update:
            if ss is not None:
                p_ss = f_ss = 0
                for i in range(0, self.n):
                    f_ss += np.sum(ss[self.idx + i, self.idx + i])
                    p_ss += np.sum(ss[self.idx + i, :]) - np.sum(ss[self.idx + i, self.idx + i])

                index = np.setdiff1d(np.arange(ss.shape[0]), np.arange(self.idx, self.idx + self.r))
                for i in np.arange(self.r):
                    f_ss += np.sum((self.r - i - 1) * ss[index, self.idx + i])
                    p_ss += np.sum(i * ss[index, self.idx + i])

                # p_ss /= self.n

                self.mo[0][number, 0] = f_ss
                self.mo[0][number, 1] = p_ss

                self.ss[0][0] = np.sum(self.mo[0][:, 0])
                self.ss[0][1] = np.sum(self.mo[0][:, 1])

        # Update SS Emissions
        if (s is not None) and (obs is not None):
            self.mo[1][number, 0] = np.dot(obs[:, 0].T, s)
            self.mo[1][number, 1] = np.sum(s)

            self.ss[1][0] = np.sum(self.mo[1][:, 0])
            self.ss[1][1] = np.sum(self.mo[1][:, 1])

    def log_block(self):

        block = np.ones((self.n, self.n)) * -np.inf
        if self.update is True:
            self.E_log_p = psi(self.posterior.p) - psi(self.posterior.p + self.posterior.f)
            self.E_log_f = psi(self.posterior.f) - psi(self.posterior.p + self.posterior.f)

        r = self.n
        for i in range(0, r - 1):
            block[i, i] = self.E_log_f
            block[i, i + 1] = self.E_log_p
        block[r - 1, r - 1] = self.E_log_f

        return block

    def mf_enter_prob(self):
        r = self.n
        ep, e1mp = np.exp(self.E_log_p), np.exp(self.E_log_f)

        return self._mf_binom(np.arange(r)[::-1], r - 1, ep, e1mp)

    @staticmethod
    def _mf_binom(k, n, p1, p2):
        return np.exp(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
                      + k * p1 + (n - k) * p2)


class TopStateNegativeBinomial:
    def __init__(self, r, p=None, s_count=None, f_count=None, toptmat0=None, order=None, mu0=None, sigma02=None):

        # Parameter Dimensions
        # self.n_obs, self.k = s_s.shape
        self.k = toptmat0.shape[0]

        # Empty container sufficient statistics
        self.ss = list()
        self.ss.append(np.zeros(2))
        self.ss.append(np.zeros((self.k, self.k)))

        # Number of expanded states
        self.n = r
        self.r = r

        # Order in the array
        self.order = order

        # Negative Binomial Parameter Count (successes and failures)
        if p is None:
            self.mu0 = mu0
            self.sigma0 = sigma02
            self.update = True
        else:
            self.update = False

        # Build the block prior matrix
        self.block_prior = np.zeros((r, r))

        # Define the prior bag of parameters States
        self.prior = ParamStorage(K=self.k)
        if self.update is True:
            self.prior.setField('p', s_count, dims=None)
            self.prior.setField('f', f_count, dims=None)

        # Define the posterior bag of parameters
        self.posterior = ParamStorage(K=self.k)
        if self.update is True:
            self.posterior.setField('p', s_count, dims=None)
            self.posterior.setField('f', f_count, dims=None)
        else:
            self.posterior.setField('p', p, dims=None)
            self.posterior.setField('f', 1. - p, dims=None)
        # self.idx remembers this states location in an expanded transition matrix, it is set
        # by the HsMM model object when it builds an expanded transition matrix
        self.idx = None

        # Data Storage
        if self.update is True:
            self.E_log_p = psi(self.posterior.p) - psi(self.posterior.p + self.posterior.f)
            self.E_log_f = psi(self.posterior.f) - psi(self.posterior.p + self.posterior.f)
        else:
            self.E_log_p = np.log(p)
            self.E_log_f = np.log(1. - p)
        self.leaving_prob = self.E_log_p
        self.enter_prob = np.ones(self.r)
        e_f = self.posterior.f
        self.enter_prob = np.log(binom.pmf(np.arange(self.r)[::-1], self.r - 1, e_f))

        # Update Parameters emissions
        e_log_a0 = psi(toptmat0) - psi(np.sum(toptmat0, axis=1))
        self.prior.setField('e_log_a', e_log_a0, dims=('K', 'K'))
        self.posterior.setField('e_log_a', e_log_a0, dims=('K', 'K'))

    def log_likelihood(self, s):

        # Expand s into: s[Time, State, Experiment]
        e_log_a = self.posterior.e_log_a

        # Compute Likelihood
        log_likelihood_vec = np.dot(s, e_log_a.T)[:, self.order]

        return log_likelihood_vec

    def update_parameters_ss(self):

        # Re-casts Parameters
        ss_transitions = self.ss[0]
        ss_emissions = self.ss[1]

        if self.update:
            self.posterior.setField('f', self.prior.f + ss_transitions[0], dims='')
            self.posterior.setField('p', self.prior.p + ss_transitions[1], dims='')
            self.E_log_p = psi(self.posterior.p) - psi(self.posterior.p + self.posterior.f)
            self.E_log_f = psi(self.posterior.f) - psi(self.posterior.p + self.posterior.f)
            self.leaving_prob = self.E_log_p
            self.enter_prob = np.log(binom.pmf(np.arange(self.r)[::-1], self.r - 1, self.posterior.f))
            e_f = self.posterior.f / (self.posterior.p + self.posterior.f)
            self.enter_prob = np.log(binom.pmf(np.arange(self.n)[::-1], self.n - 1, e_f))

        # Update Parameters emissions
        toptmat = self.prior.e_log_a + ss_emissions
        e_log_a = psi(toptmat) - psi(np.sum(toptmat, axis=1))
        self.posterior.setField('e_log_a', e_log_a, dims=('K', 'K'))

    def clear_ss(self):
        self.ss = list()
        self.ss.append(np.zeros(2))
        self.ss.append(np.zeros((self.k, self.k)))

    def increase_ss(self, ss=None, aik=None, _=None):

        if self.update:
            if ss is not None:
                # Update SS Transitions
                a = b = 0.
                for i in np.arange(self.r):
                    a += np.sum(ss[self.idx + i, self.idx + i])
                    b += np.sum(ss[self.idx + i, :]) - np.sum(ss[self.idx + i, self.idx + i])

                index = np.setdiff1d(np.arange(ss.shape[0]), np.arange(self.idx, self.idx + self.r))
                for i in np.arange(self.r):
                    a += np.sum((self.r - i - 1) * ss[index, self.idx + i])
                    b += np.sum(i * ss[index, self.idx + i])

                self.ss[0][0] += a
                self.ss[0][1] += b

        if (aik is not None) and (_ is None):
            self.ss[1][:, :] += aik

    def log_block(self):

        block = np.ones((self.r, self.r)) * -np.inf
        if self.update is True:
            self.E_log_p = psi(self.posterior.p) - psi(self.posterior.p + self.posterior.f)
            self.E_log_f = psi(self.posterior.f) - psi(self.posterior.p + self.posterior.f)

        r = self.r
        for i in range(0, r - 1):
            block[i, i] = self.E_log_f
            block[i, i + 1] = self.E_log_p
        block[r - 1, r - 1] = self.E_log_f

        return block

    def mf_enter_prob(self):
        r = self.r
        ep, e1mp = np.exp(self.E_log_p), np.exp(self.E_log_f)
        return self._mf_binom(np.arange(r)[::-1], r - 1, ep, e1mp)

    @staticmethod
    def _mf_binom(k, n, p1, p2):
        return np.exp(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1) + k * p1 + (n - k) * p2)


class NegativeBinomialGEOBin:
    def __init__(self, r, cut0, cut1, p=None, s_count=None, f_count=None):

        # Saving Variables
        self.save = False
        self.e_log_a = []

        # Empty container sufficient statistics
        self.mo = list()
        self.ss = list()
        self.ss.append(np.zeros(2))
        self.ss.append(np.zeros(2))

        # SO
        self.rho = 1.
        self.it = 1.
        self.rhodelay = 20.
        self.rhoexp = 0.55

        # Number of expanded states
        self.n = r
        self.r = r

        # Negative Binomial Parameter Count (successes and failures)
        if p is None:
            self.s_count = s_count
            self.f_count = f_count
            self.update = True
        else:
            self.update = False

        # Build the block prior matrix
        self.block_prior = np.zeros((r, r))

        # Define the prior bag of parameters
        self.prior = ParamStorage()
        if self.update is True:
            self.prior.setField('p', s_count, dims=None)
            self.prior.setField('f', f_count, dims=None)
        self.prior.setField('emit0', cut0, dims=None)
        self.prior.setField('emit1', cut1, dims=None)

        # Define the posterior bag of parameters
        self.posterior = ParamStorage()
        if self.update is True:
            self.posterior.setField('p', s_count, dims=None)
            self.posterior.setField('f', f_count, dims=None)
        self.posterior.setField('emit0', cut0, dims=None)
        self.posterior.setField('emit1', cut1, dims=None)

        # self.idx remembers this states location in an expanded transition matrix, it is set
        # by the HsMM model object when it builds an expanded transition matrix
        self.idx = None

        # Data Storage
        if self.update is True:
            self.E_log_p = psi(self.posterior.p) - psi(self.posterior.p + self.posterior.f)
            self.E_log_f = psi(self.posterior.f) - psi(self.posterior.p + self.posterior.f)
        else:
            self.E_log_p = np.log(p)
            self.E_log_f = np.log(1. - p)
        self.leaving_prob = self.E_log_p
        self.enter_prob = np.ones(self.n)
        e_f = np.exp(self.E_log_f) / (np.exp(self.E_log_p) + np.exp(self.E_log_f))
        self.enter_prob = np.log(binom.pmf(np.arange(self.n)[::-1], self.n - 1, e_f))

    def log_likelihood(self, obs=None, s=None):

        if obs is not None:
            e_success = self.posterior.emit0
            e_fail = self.posterior.emit1

            l_e = np.zeros(2)
            psi_norm = psi(e_success + e_fail)
            l_e[0] = psi(e_success) - psi_norm
            l_e[1] = psi(e_fail) - psi_norm

            cut = obs[:, 0]
            log_likelihood_vec = (cut * l_e[0] + l_e[1])
        else:
            log_likelihood_vec = np.zeros((s.shape[0]))

        if s is not None:
            if self.e_log_a is not None:
                log_likelihood_vec += np.dot(s, self.e_log_a)

        return log_likelihood_vec

    def kl_term(self):

        psi_norm = psi(self.posterior.emit0 + self.posterior.emit1)
        expt_terms = (self.prior.emit1 - self.posterior.emit1) * (psi(self.posterior.emit1) - psi_norm) + (
                self.prior.emit0 - self.posterior.emit0) * (psi(self.posterior.emit0) - psi_norm)

        post_kl = gammaln(self.posterior.emit1 + self.posterior.emit0) - gammaln(self.posterior.emit1) - gammaln(
            self.posterior.emit0)
        prior_kl = gammaln(self.prior.emit1 + self.prior.emit0) - gammaln(self.prior.emit1) - gammaln(self.prior.emit0)

        return expt_terms + prior_kl - post_kl

    def update_parameters_ss(self, e_log_a=None):

        # Re-casts Parameters
        ss_transitions = self.ss[0]
        ss_emissions = self.ss[1]

        if self.update is True:
            self.posterior.setField('f', self.prior.f + ss_transitions[0], dims='')
            self.posterior.setField('p', self.prior.p + ss_transitions[1], dims='')
            self.E_log_p = psi(self.posterior.p) - psi(self.posterior.p + self.posterior.f)
            self.E_log_f = psi(self.posterior.f) - psi(self.posterior.p + self.posterior.f)
            self.leaving_prob = self.E_log_p
            self.enter_prob = np.log(binom.pmf(np.arange(self.r)[::-1], self.r - 1, self.posterior.f))
            e_f = self.posterior.f / (self.posterior.p + self.posterior.f)
            self.enter_prob = np.log(binom.pmf(np.arange(self.n)[::-1], self.n - 1, e_f))

        # Update Parameters emissions
        self.posterior.setField('emit0', self.prior.emit0 + ss_emissions[0], dims=None)
        self.posterior.setField('emit1', self.prior.emit1 + ss_emissions[1], dims=None)
        if e_log_a is not None:
            self.e_log_a = e_log_a

        if self.save:
            emit = [self.posterior.emit0, self.posterior.emit1]
            np.save('emit_params{0}_{1}.npy'.format(self.idx, self.it), emit)
            self.it += 1

    def update_parameters_so(self, e_log_a=None):

        if self.it < 500:
            # Re-casts Parameters
            ss_transitions = self.ss[0]
            ss_emissions = self.ss[1]

            if self.update is True:
                p = self.prior.p + ss_transitions[0]
                f = self.prior.f + ss_transitions[1]
                new_p = p * self.rho + (1. - self.rho) * self.posterior.p
                new_f = f * self.rho + (1. - self.rho) * self.posterior.f
                self.posterior.setField('p', new_p, dims='')
                self.posterior.setField('f', new_f, dims='')

                self.E_log_p = psi(self.posterior.p) - psi(self.posterior.p + self.posterior.f)
                self.E_log_f = psi(self.posterior.f) - psi(self.posterior.p + self.posterior.f)
                self.leaving_prob = self.E_log_p
                self.enter_prob = np.log(binom.pmf(np.arange(self.r)[::-1], self.r - 1, self.posterior.f))

            # Update Parameters emissions
            emit0 = self.prior.emit0 + ss_emissions[0]
            emit1 = self.prior.emit1 + ss_emissions[1]
            emit0_new = emit0 * self.rho + (1. - self.rho) * self.posterior.emit0
            emit1_new = emit1 * self.rho + (1. - self.rho) * self.posterior.emit1
            self.posterior.setField('emit0', emit0_new, dims=None)
            self.posterior.setField('emit1', emit1_new, dims=None)
            if e_log_a is not None:
                self.e_log_a = e_log_a

            # Update Rho
            self.rho = (1 + self.it + self.rhodelay) ** (-1.0 * self.rhoexp)
            if self.save:
                emit = [self.posterior.emit0, self.posterior.emit1]
                np.save('emit_paramsSO{0}_{1}.npy'.format(self.idx, self.it), emit)
            self.it += 1

    def clear_ss(self):
        self.ss = list()
        self.ss.append(np.zeros(2))
        self.ss.append(np.zeros(2))

    def prepare_ss(self, opt, number):

        if opt == "batch":
            self.ss = list()
            self.ss.append(np.zeros(2))
            self.ss.append(np.zeros(2))
        elif opt == "so":
            self.ss = list()
            self.ss.append(np.zeros(2))
            self.ss.append(np.zeros(2))
            self.rho = (1. + self.it + self.rhodelay) ** (- self.rhoexp)
        elif opt == "mo":
            if len(self.mo) == 0:
                self.mo.append(np.zeros((number, 2)))
                self.mo.append(np.zeros((number, 2)))

    def increase_ss(self, ss=None, s=None, obs=None, suf_stat=None, amp=1.):

        if self.update:
            if ss is not None:
                # Update parameters p and f (using epsilons)
                p_ss = f_ss = 0
                for i in range(0, self.n):
                    f_ss += np.sum(ss[self.idx + i, self.idx + i])
                    p_ss += np.sum(ss[self.idx + i, :]) - np.sum(ss[self.idx + i, self.idx + i])

                index = np.setdiff1d(np.arange(ss.shape[0]), np.arange(self.idx, self.idx + self.r))
                for i in np.arange(self.r):
                    f_ss += np.sum((self.r - i - 1) * ss[index, self.idx + i])
                    p_ss += np.sum(i * ss[index, self.idx + i])

                # p_ss /= self.n

                self.ss[0][0] += amp * p_ss
                self.ss[0][1] += amp * f_ss

        if suf_stat is None:
            if (s is not None) and (obs is not None):
                # Update SS Emmissions
                self.ss[1][0] += amp * np.dot(obs[:, 0].T, s)
                self.ss[1][1] += amp * np.sum(s)
        else:
            self.ss[1][0] += amp * suf_stat[0]
            self.ss[1][1] += amp * suf_stat[1]

    def replace_ss(self, ss=None, s=None, obs=None, number=None):

        # Update SS Transitions
        if self.update:
            if ss is not None:
                p_ss = f_ss = 0
                for i in range(0, self.n):
                    f_ss += np.sum(ss[self.idx + i, self.idx + i])
                    p_ss += np.sum(ss[self.idx + i, :]) - np.sum(ss[self.idx + i, self.idx + i])

                index = np.setdiff1d(np.arange(ss.shape[0]), np.arange(self.idx, self.idx + self.r))
                for i in np.arange(self.r):
                    f_ss += np.sum((self.r - i - 1) * ss[index, self.idx + i])
                    p_ss += np.sum(i * ss[index, self.idx + i])

                # p_ss /= self.n

                self.mo[0][number, 0] = f_ss
                self.mo[0][number, 1] = p_ss

                self.ss[0][0] = np.sum(self.mo[0][:, 0])
                self.ss[0][1] = np.sum(self.mo[0][:, 1])

        # Update SS Emissions
        if (s is not None) and (obs is not None):
            self.mo[1][number, 0] = np.dot(obs[:, 0].T, s)
            self.mo[1][number, 1] = np.sum(s)

            self.ss[1][0] = np.sum(self.mo[1][:, 0])
            self.ss[1][1] = np.sum(self.mo[1][:, 1])

    def log_block(self):

        block = np.ones((self.n, self.n)) * -np.inf
        if self.update is True:
            self.E_log_p = psi(self.posterior.p) - psi(self.posterior.p + self.posterior.f)
            self.E_log_f = psi(self.posterior.f) - psi(self.posterior.p + self.posterior.f)

        r = self.n
        for i in range(0, r - 1):
            block[i, i] = self.E_log_f
            block[i, i + 1] = self.E_log_p
        block[r - 1, r - 1] = self.E_log_f

        return block

    def mf_enter_prob(self):
        r = self.n
        ep, e1mp = np.exp(self.E_log_p), np.exp(self.E_log_f)

        return self._mf_binom(np.arange(r)[::-1], r - 1, ep, e1mp)

    @staticmethod
    def _mf_binom(k, n, p1, p2):
        return np.exp(gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
                      + k * p1 + (n - k) * p2)
