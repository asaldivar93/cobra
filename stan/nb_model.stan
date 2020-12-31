data {
  int N; // Number of features in sample
  int Y[N]; // Counts of non-zero interactions for each feature

  real a_m; // Shape parameter for gamma prior on mu
  real b_m; // Rate parameter for gamma prior on mu
  real a_p; // Shape parameter for gamma prior on phi
  real b_p; // Rate parameter for gamma prior on phi

  real a_0; // Value to use as the mean for null NB
}

parameters {
  real<lower=0> mu; // mean of negative binomial distribution
  real<lower=0> phi; // dispersion of negative binomial distribution
}

model {
  for(i in 1:N) {
    Y[i] ~ neg_binomial_2(mu, phi);
  }
  mu ~ gamma(a_m, b_m);
  phi ~ gamma(a_p, b_p);
}

generated quantities {
  real<lower=0> mu_0; // Mean of null NB distribution
  real Y_sim[N]; // Vector to store the simulated values
  real Y_0[N]; // Vectore to store H0 values
  real<lower=0> T_sum; // Test statistic indicating the probability that sum(Y_sim) >= sum(Y)
  real<lower=0> T_max; // Test statistic indicating the probability that max(Y_sim) > max(Y)
  real<lower=0> p_sim; // p value indicating the probability that sum(Y_sim) >= sum(Y_0)

  for(i in 1:N) {
    Y_sim[i] = neg_binomial_2_rng(mu, phi);
    mu_0 = gamma_rng(a_0, 1);
    Y_0[i] = neg_binomial_2_rng(mu_0, 4);
  }

  T_sum = sum(Y_sim) >= sum(Y);
  T_max = max(Y_sim) >= max(Y);
  p_sim = sum(Y_sim) >= sum(Y_0);
}
