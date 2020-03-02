from .LibFwdBwd import FwdBwdAlg_cpp
from scipy.special import psi
import numpy as np


# ################################################################
# Compress - Expand Methods
def log_expand(pi, tmat, states, n, k, likelihood=None, obs1=None, obs2=None):

    # ######################################################################################################
    # Log( expand TMAT )
    # Calculate the compressed tmat transitions from their posteriors.
    if likelihood is not None:
        log_comp_tmat = np.log(np.clip(tmat, 1e-15, 1))
    else:
        log_comp_tmat = psi(tmat) - psi(np.sum(tmat, axis=1))

    log_comp_tmat[tmat == 0] = - np.inf
    for i_ in np.arange(tmat.shape[0]):
        if np.sum(tmat[i_, :] > 0) == 1:
            log_comp_tmat[i_, tmat[i_, :] > 0] = 0.

    # Build out the expanded transition matrix to feed to the forward-backward algorithm
    log_expand_tmat = - np.inf * np.ones((n, n))

    # Populate the expanded transition with blocks returned by state.block() method
    for s_ in states:
        log_expand_tmat[s_.idx:s_.idx + s_.r, s_.idx:s_.idx + s_.r] = s_.log_block()

    # Add transitions between compressed states to the expanded transition matrix
    for i in range(0, k):
        in_index = states[i].idx + states[i].r - 1
        leaving_prob = states[i].leaving_prob
        for j in range(0, k):
            if tmat[i, j] > 0:
                out_index = states[j].idx
                out_entering_lprob = states[j].enter_prob
                log_expand_tmat[in_index, out_index:out_index + states[j].r] = \
                    log_comp_tmat[i, j] + leaving_prob + out_entering_lprob
    # ######################################################################################################

    # ######################################################################################################
    # Make the compressed log pi vector from posterior distributions
    if likelihood is not None:
        log_comp_pi = np.log(pi)
    else:
        log_comp_pi = psi(pi) - psi(np.sum(pi))

    log_expand_pi = np.zeros(n)
    for i_, s_ in enumerate(states):
        log_expand_pi[s_.idx:s_.idx + s_.r] = log_comp_pi[i_] + s_.enter_prob
    # ######################################################################################################

    # ######################################################################################################
    # Compute log likelihood
    if likelihood is not None:  # Likelihood with/without obs2
        log_expand_likelihood = np.zeros((likelihood.shape[0], n))

        for i in range(0, k):
            if obs2 is not None:
                log_expand_likelihood[:, states[i].idx:states[i].idx + states[i].r] = \
                    np.log(likelihood[:, i])[:, None] + states[i].log_likelihood(s=obs2)[:, None]
            else:
                log_expand_likelihood[:, states[i].idx:states[i].idx + states[i].r] = \
                    np.log(likelihood[:, i])[:, None]

    elif obs1 is not None:  # Obs1 with/without obs2
        log_expand_likelihood = np.zeros((obs1.shape[0], n))

        for i in range(0, k):
            if obs2 is None:
                log_expand_likelihood[:, states[i].idx:states[i].idx + states[i].r] = \
                    np.squeeze(states[i].log_likelihood(obs1))[:, None]
            else:
                log_expand_likelihood[:, states[i].idx:states[i].idx + states[i].r] = \
                    states[i].log_likelihood(obs1, obs2)[:, None]

    else:  # Only Obs2
        log_expand_likelihood = np.zeros((obs2.shape[0], n))

        for i in range(0, k):
            log_expand_likelihood[:, states[i].idx:states[i].idx + states[i].r] += \
                np.squeeze(states[i].log_likelihood(s=obs2))[:, None]
    # ######################################################################################################

    return log_expand_pi, log_expand_tmat, log_expand_likelihood


def compress(resp, resp_pair, states):
    # Read Parameters
    t = resp.shape[0]
    n_out = len(states)

    # Compress Resp
    out_s_s = np.zeros((t, n_out))
    for i_, s_ in enumerate(states):
        out_s_s[:, i_] = np.sum(resp[:, s_.idx:s_.idx + s_.r], axis=1)

    # Compress RespPair
    out_s_ss = np.zeros((n_out, n_out))
    # Diagonal terms
    for i_ in np.arange(n_out):
        index1 = states[i_].idx + states[i_].r - 1
        index2 = states[i_].idx
        out_s_ss[i_, i_] = resp_pair[index1, index2]

    # Off Diagonal terms
    for i_ in np.arange(n_out):
        for j_ in np.arange(n_out):
            if not (i_ == j_):
                index1 = states[i_].idx + states[i_].r - 1
                index2 = states[j_].idx
                out_s_ss[i_, j_] = resp_pair[index1, index2]

    return out_s_s, out_s_ss
# ################################################################


# ################################################################
def message_passing_incremental(pi, tmat, states, n_expanded, n_compress, length, data=None, likelihood=None, opt='mo'):

    # Initialize Storage
    lmarg_pr_obs = 0
    ss_pi = np.zeros(n_compress)
    ss_tmat = np.zeros((n_compress, n_compress))

    # Clear sufficient statistics before pass
    [s.prepare_ss(opt, len(length)) for s in states]

    # PREPARE COMPUTATIONAL BATCHES
    position = np.insert(np.cumsum(length), 0, 0)
    np.random.seed(0)
    index = np.random.permutation(np.arange(position.shape[0] - 1))

    for i_ in index:
        st = position[i_]
        end = position[i_ + 1]

        # Log expected parameters of batch
        if likelihood is None:
            log_xpd_pi, log_xpd_tmat, log_xpd_likelihood = log_expand(pi, tmat, states, n_expanded, n_compress,
                                                                      obs1=data[st:end, None], likelihood=likelihood)
        else:
            log_xpd_pi, log_xpd_tmat, log_xpd_likelihood = log_expand(pi, tmat, states, n_expanded, n_compress,
                                                                      obs1=data, likelihood=likelihood[st:end])

        # Collect Forward - Backward
        resp, resp_pair, lmarg_pr_obs_temp = fw_bw(log_xpd_pi, log_xpd_tmat, log_xpd_likelihood)

        # Compress Matrices
        s, ss_tmat_acc, = compress(resp, resp_pair, states)
        lmarg_pr_obs += lmarg_pr_obs_temp

        if opt == 'so':
            amp_factor = data.shape[0] / np.abs(end - st)
            [t.increase_ss(resp_pair, s[:, e_], data[st:end], amp=amp_factor) for e_, t in enumerate(states)]
        elif opt == 'mo':
            [t.replace_ss(resp_pair, s[:, e_], data[st:end, None], i_) for e_, t in enumerate(states)]

        ss_tmat += ss_tmat_acc
        if st == 0:
            ss_pi += s[0]

        if opt == 'so':
            [t.update_parameters_so() for t in states]
            [s.clear_ss() for s in states]
        elif opt == 'mo':
            [t.update_parameters_ss() for t in states]

    return lmarg_pr_obs, ss_pi, ss_tmat


def message_passing_posterior_state(pi, tmat, states, n_expanded, n_compress, length, data=None, likelihood=None):

    # PREPARE COMPUTATIONAL BATCHES
    position = np.insert(np.cumsum(length), 0, 0)
    np.random.seed(0)
    index = np.random.permutation(np.arange(position.shape[0] - 1))
    s_out = np.zeros((np.sum(length), n_compress))

    # FORWARD - BACKWARD ON EACH DATASET
    for i_ in index:
        st = position[i_]
        end = position[i_ + 1]

        # Log expected parameters of batch
        if likelihood is not None:
            log_xpd_pi, log_xpd_tmat, log_xpd_likelihood = log_expand(pi, tmat, states, n_expanded, n_compress,
                                                                      likelihood=likelihood[st:end])
        else:
            log_xpd_pi, log_xpd_tmat, log_xpd_likelihood = log_expand(pi, tmat, states, n_expanded, n_compress,
                                                                      obs1=data[st:end])

        # Collect Forward - Backward
        resp, resp_pair, _ = fw_bw(log_xpd_pi, log_xpd_tmat, log_xpd_likelihood)

        s_out[st:end], _ = compress(resp, resp_pair, states)

    return s_out


def message_passing_multi(pi, tmat, states, top_states, s_s, n_exp, n_expanded, n_compress, length,
                          data=None, likelihood=None):

    # PRE-COMPUTING ORDERING OF DATASETS
    position = np.insert(np.cumsum(length), 0, 0)
    index = np.random.permutation(np.arange(position.shape[0] - 1))

    # Clear sufficient statistics before pass
    [t.clear_ss() for t in top_states]
    [[t.clear_ss() for t in s] for s in states]

    # Initialize Storage
    l_marg_pr = 0
    ss_pi = np.zeros(n_compress)
    ss_toptmat = ss_tmat = np.zeros((n_compress, n_compress))

    for i_ in index:

        # Compute for Bottom Chains
        cum_s_e = np.array([])
        st = position[i_]
        end = position[i_ + 1]

        for e_ in np.arange(n_exp):
            # Compute Forward - Backward
            if likelihood is not None:
                log_xpd_pi, log_xpd_tmat, log_xpd_likelihood = log_expand(pi, tmat, states[e_], n_expanded, n_compress,
                                                                          likelihood=likelihood[e_][st:end, :],
                                                                          obs2=s_s[st:end, :])
            else:
                log_xpd_pi, log_xpd_tmat, log_xpd_likelihood = log_expand(pi, tmat, states[e_], n_expanded, n_compress,
                                                                          obs1=data[st:end, e_][:, None],
                                                                          obs2=s_s[st:end, :])
            resp, resp_pair, l_marg_prob_temp = fw_bw(log_xpd_pi, log_xpd_tmat, log_xpd_likelihood)

            # Compress Matrices
            s_e, ss_tmat_acc, = compress(resp, resp_pair, states[e_])
            l_marg_pr += l_marg_prob_temp
            # l_marg_pr += - np.sum(np.einsum('bi, ij-> bij', s_s[st:end, :],
            #                                 top_states[0].posterior.e_log_a) * s_e[:, :, None])

            if likelihood is None:
                # Increase SS
                ss_tmat += ss_tmat_acc
                [t.increase_ss(resp_pair) for t in top_states]
                [t.increase_ss(None, s_e[:, l], data[st:end, e_][:, None]) for l, t in enumerate(states[e_])]
                [[t.increase_ss(resp_pair) for t in s] for s in states]

            if cum_s_e.shape[0] == s_e.shape[0]:
                cum_s_e += s_e
            else:
                cum_s_e = s_e
        cum_s_e = cum_s_e / n_exp

        # Compute for Top Chain
        # Log expected parameters
        log_xpd_pi, log_xpd_tmat, log_xpd_likelihood = log_expand(pi, tmat, top_states, n_expanded, n_compress,
                                                                  obs2=cum_s_e)
        resp, resp_pair, l_marg_prob_temp = fw_bw(log_xpd_pi, log_xpd_tmat, log_xpd_likelihood)
        s_s[st:end, :], ss_tmat_acc, = compress(resp, resp_pair, top_states)
        l_marg_pr += l_marg_prob_temp

        if st == 0:
            ss_pi += s_s[0] + cum_s_e[0]

        # lmarg_pr_obs += lmarg_pr_obs_temp
        if likelihood is None:
            ss_tmat += ss_tmat_acc
            [[t.increase_ss(resp_pair) for t in s] for s in states]

        ss_toptmat += np.dot(cum_s_e.T, s_s[st:end, :])
        [t.increase_ss(resp_pair, ss_toptmat) for t in top_states]

    return l_marg_pr, ss_pi, ss_tmat
# ################################################################


# ################################################################
# Fwd-Bwd Methods
def fw_bw(lpi_init, lpi_mat, log_soft_ev):
    init = np.exp(lpi_init)
    mat = np.exp(lpi_mat)

    log_soft_ev = np.asarray(log_soft_ev, dtype=np.float64)
    lognorm_c = np.max(log_soft_ev, axis=1)
    log_soft_ev = log_soft_ev - lognorm_c[:, np.newaxis]
    soft_ev = np.exp(log_soft_ev)

    resp, resp_pair, lmarg_pr_seq = FwdBwdAlg_cpp(init, mat, soft_ev)

    return resp, resp_pair, lmarg_pr_seq + lognorm_c.sum()
# ################################################################
