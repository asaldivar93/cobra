data{
  int N; // Number of samples
  real Y[N]; // Interaction coefficients
  real L;
  real U;
  real mu_0; // mean prior for mu
  real sigma_0; // sigma prior for mu
  real a_n; // Shape parameter for gamma prior on nu
  real b_n; // Rate parameter for gamma prior on nu
  real a_s; // Shape parameter for gamma prior on sigma
  real b_s; // Rate parameter for gamma prior on sigma
}

parameters {
  real<lower=L, upper=U> mu; // mean of data
  real<lower=0> sigma; // variance
  real<lower=0> nu; // degrees of fredom
}

model {
  for(i in 1:N) {
    Y[i] ~ student_t(nu, mu, sigma);
  }
  mu ~ normal(mu_0, sigma_0);
  nu ~ gamma(a_n, b_n);
  sigma ~ gamma(a_s, b_s);
}
