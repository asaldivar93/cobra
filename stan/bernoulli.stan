data {
  int N; // number of observations
  int<lower=0, upper=1> Y[N]; // Vector of N observations
  real a; //
  real b; //
  real a0; //
  real b0; //
}

parameters {
  real<lower=0, upper=1> theta; //
}

model {
  Y ~ bernoulli(theta);
  theta ~ beta(a, b);
}

generated quantities {
  int<lower=0, upper=1> Y_sim[N];
  int<lower=0, upper=1> Y_0[N];
  real T_sum;
  real p_g;
  real p_l;
  real<lower=0, upper=1> theta_0; //

  theta_0 = beta_rng(a0, b0);
  for (i in 1:N){
    Y_0[i] = bernoulli_rng(theta_0);
    Y_sim[i] = bernoulli_rng(theta);
  }

  T_sum = sum(Y_sim) >= sum(Y);
  p_g = sum(Y) >= sum (Y_0);
  p_l = sum(Y) <= sum (Y_0);
}
