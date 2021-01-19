data {
  int N; // Number of features in sample
  int Y[N]; // Counts of non-zero interactions for each feature

  real a_m; // Shape parameter for gamma prior on mu
  real b_m; // Rate parameter for gamma prior on mu
  real a_p; // Shape parameter for gamma prior on phi
  real b_p; // Rate parameter for gamma prior on phi

  real a_b; // Bias towards 1 parameter for beta prior
  real b_b; // Bias towards 0 parameter for beta prior

  real a_0; // Value to use as the mean for null NB
}

parameters {
  real<lower=0> mu; // mean of negative binomial distribution
  real<lower=0> phi; // dispersion of negative binomial distribution
  real<lower=0, upper=1> theta; // mean of bernoulli distribution
}

model {
  for(i in 1:N) {
    if (Y[i] == 0)
      target += log_sum_exp(bernoulli_lpmf(1 | theta), bernoulli_lpmf(0 | theta) + neg_binomial_2_lpmf( Y[i] | mu, phi));
    else
      target += bernoulli_lpmf(0 | theta) + neg_binomial_2_lpmf( Y[i] | mu, phi);
  }
  mu ~ gamma(a_m, b_m);
  phi ~ gamma(a_p, b_p);
  theta ~ beta(a_b, b_b);
}

generated quantities {
  real<lower=0, upper=1> binary; // Indicator variable for zero-inflated NB distribution
  real<lower=0, upper=1> binary_0; // Indicator variable for zero-inflated NB distribution
  real<lower=0> mu_0; // Mean of null NB distribution
  real Y_sim[N]; // Vector to store the simulated values
  real Y_0[N]; // Vectore to store H0 values
  real<lower=0> T_sum; // Test statistic indicating the probability that sum(Y_sim) >= sum(Y)
  real<lower=0> T_max; // Test statistic indicating the probability that max(Y_sim) > max(Y)
  real<lower=0> p_sum; // p value indicating the probability that sum(Y) >= sum(Y_0)
  real<lower=0> p_max; // p value indicating the probability that max(Y) >= max(Y_0)

  binary = bernoulli_rng(theta);
  binary_0 = bernoulli_rng(0.5);
  mu_0 = gamma_rng(a_0, 1);
  for(i in 1:N) {
    Y_sim[i] = (0) ^ binary * neg_binomial_2_rng(mu, phi) ^ (1 - binary);
    Y_0[i] = (0) ^ binary_0 * neg_binomial_2_rng(mu_0, 4) ^ (1 - binary_0);
  }

  T_sum = sum(Y_sim) >= sum(Y);
  T_max = max(Y_sim) >= max(Y);
  p_sum = sum(Y) >= sum(Y_0);
  p_max = max(Y) >= max(Y_0);
}
