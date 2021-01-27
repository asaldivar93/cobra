import numpy as np
import pickle
import pystan


def run_zinb(
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
        pars = ['mu', 'phi', 'theta', 'T_sum', 'T_max', 'p_sum', 'p_max']
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
        pars = ['mu', 'phi', 'T_sum', 'T_max', 'p_sum', 'p_max']
    fit = model.sampling(
        data=data, pars=pars, chains=chains, init=init, iter=iter, warmup=warmup, seed=194838
    )
    return fit


def run_nb(
    Y,
    a_mu_prior = [],
    a_phi_prior = [],
    b_mu_prior = 1,
    b_phi_prior = 1,
    a_beta_prior = 0.9,
    b_beta_prior = 0.9,
    zero_inflated = True,
    chains = 4,
    iter = 2000
):
    N = len(Y)
    sample_mean = np.mean(Y)
    sample_variance = np.var(Y)
    a_H0 = np.median(np.linspace(1, N, N))

    if not a_mu_prior:
        a_mu_prior = sample_mean * b_mu_prior
    if not a_phi_prior:
        a_phi_prior = (sample_mean**2 / abs(sample_variance - sample_mean)) * b_phi_prior

    warmup = int(np.floor(iter / 2))
    if zero_inflated:
        data = dict(
            N = N,
            Y = Y,
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
            N = N,
            Y = Y,
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
        data=data, pars=pars, chains=chains, init=init, iter=iter, warmup=warmup, seed=194838
    )
    return fit


def get_init(data):
    init = dict(
        mu = np.mean(data['Y']),
        phi = np.var(data['Y']),
        theta = sum(data['Y'] > 0) / data['N']
    )

    return init


def run_st(
    Y,
    L,
    U,
    mu_0=0,
    sigma_0=2,
    a_n=2,
    b_n=0.1,
    a_s=2,
    b_s=0.1,
    chains=4,
    iter=2000
):

    data = dict(
        N=len(Y),
        Y=Y,
        L=L,
        U=U,
        mu_0=mu_0,
        sigma_0=sigma_0,
        a_n=a_n,
        b_n=b_n,
        a_s=a_s,
        b_s=b_s
    )
    warmup = int(np.floor(iter / 2))
    model = pickle.load(open('stan/st_model.pkl', 'rb'))
    fit = model.sampling(
        data=data, iter=iter, chains=chains, warmup=warmup, seed=194838
    )

    return fit


def run_bernoulli(
    Y,
    a,
    b,
    a0,
    b0,
    chains=4,
    iter=2000,
):
    warmup = int(np.floor(iter / 2))
    data = dict(
        N=len(Y),
        Y=Y,
        a=a,
        b=b,
        a0=a0,
        b0=b0,
    )
    model = pickle.load(open('stan/bernoulli.pkl', 'rb'))
    pars = ['theta', 'theta_0', 'T_sum', 'p_g', 'p_l']
    fit = model.sampling(
        data=data, pars=pars, iter=iter, chains=chains, warmup=warmup, seed=194838
    )

    return fit
