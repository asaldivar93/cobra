
import numpy as np
import pickle
import pystan


def run_sampling(
    interaction_matrix,
    model = [],
    positives = True,
    zero_cutoff = 0.01,
    a_mu_prior = [],
    a_phi_prior = [],
    b_mu_prior = 1,
    b_phi_prior = 1,
    a_beta_prior = 0.9,
    b_beta_prior = 0.9,
    chains = 4,
    iter = 2000
):
    IC = interaction_matrix.copy()
    if positives:
        IC[IC >= zero_cutoff] = 1
        IC[IC < zero_cutoff] = 0
    else:
        IC[IC > -zero_cutoff] = 0
        IC[IC <= -zero_cutoff] = 1

    n_models = len(IC.columns)
    non_zero_counts = np.zeros(n_models, dtype = np.int8)
    for j in range(0, n_models):
        non_zero_counts[j] = int(IC.iloc[:, j].sum())

    sample_mean = np.mean(non_zero_counts)
    sample_variance = np.var(non_zero_counts)
    a_H0 = np.median(np.linspace(1, n_models, n_models))

    if not a_mu_prior:
        a_mu_prior = sample_mean * b_mu_prior
    if not a_phi_prior:
        a_phi_prior = (sample_mean**2 / abs(sample_variance - sample_mean)) * b_phi_prior

    warmup = int(np.floor(iter / 2))
    if positives:
        data = dict(
            N = n_models,
            Y = non_zero_counts,
            a_m = a_mu_prior,
            b_m = b_mu_prior,
            a_p = a_phi_prior,
            b_p = b_phi_prior,
            a_b = a_beta_prior,
            b_b = b_beta_prior,
            a_0 = a_H0
        )
        init = [get_init(data) for b in range(chains)]
        model = pickle.load(open('stan/zinb_model.pkl', 'rb'))
        pars = ['mu', 'phi', 'theta', 'T_sum', 'T_max', 'p_sim']
    else:
        data = dict(
            N = n_models,
            Y = non_zero_counts,
            a_m = a_mu_prior,
            b_m = b_mu_prior,
            a_p = a_phi_prior,
            b_p = b_phi_prior,
            a_0 = a_H0
        )
        init = [dict(mu = np.mean(data['Y']), phi = np.var(data['Y']),) for b in range(chains)]
        model = pickle.load(open('stan/nb_model.pkl', 'rb'))
        pars = ['mu', 'phi', 'T_sum', 'T_max', 'p_sim']
    fit = model.sampling(
        data=data, pars=pars, chains=chains, init=init, iter=iter, warmup=warmup
    )
    return fit


def get_init(data):
    init = dict(
        mu = np.mean(data['Y']),
        phi = np.var(data['Y']),
        theta = sum(data['Y'] > 0) / data['N']
    )

    return init
