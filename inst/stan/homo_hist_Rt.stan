//
// Stan implementation of Bertozzi et al. 2020 within a Bayesian framework
// Adapted to discrete data.
//

data{
  int<lower=1> D_seed; // number of days to seed
  int<lower=D_seed> D; // number of days
  int<lower=0> y[D]; // daily count
  vector<lower=0>[D] mu; // proposed daily import rate
  int<lower=1> S; // max generation interval length
  real<lower=0> mean_omega; // generation interval mean
  real<lower=0> sd_omega; // generation interval standard deviation
  int<lower=1> delta; // bin width (days) for histogram estimator
  real<lower=0> mean_log_r; // prior mean for log reproduction numbers
  real<lower=0> sd_log_r; // prior standard deviation for log reproduction numbers
  int<lower = 0> y_ahead; // forward observations
  real<lower=0> mu_ahead; // forward daily import rate
  int<lower=0> c; // case count regularising constant
}
transformed data{
  vector<lower=0>[D] y_reg; // regularised counts
  int <lower=1> B; // number of bins for histogram estimator
  real<lower=0> a; // generation interval shape
  real<lower=0> b; // generation interval rate
  simplex[S] omega; // generation interval pmf
  vector<upper=0>[S] log_omega; // generation interval pmf
  //  regularise low counts
  for (i in 1:D) {
    y_reg[i] = max(y[i], c);
  }
  // bins
  B = (D / delta) + 1; // default integer rounding behaviour desired
  // set generation interval
  b = mean_omega / (sd_omega^2);
  a = mean_omega * b;
  omega[1] = gamma_cdf(1, a, b);
  for (i in 2:S) omega[i] = gamma_cdf(i, a, b) - gamma_cdf(i-1, a, b);
  omega = omega / sum(omega);
  log_omega = log(omega);
}
parameters{
  vector<lower=0>[B] r; // reproduction number bin values
}
transformed parameters{
  vector<lower=0>[D] R; // time-dependent reproductive numbers
  vector<lower=0>[D] psi; // local infection rate
  // reproduction number
  for (i in 1:D) {
    R[D - i + 1] = r[((i -  1) / delta) + 1]; // default integer rounding behaviour desired
  }
  // burden of infection
  psi = rep_vector(0, D);
  for (i in 2:D) {
    if (i <= S) {
      for (j in 1:(i - 1)) {
        psi[i] = psi[i] + (R[i-j] * y_reg[i-j] * omega[j]);
      }
    } else {
      for (j in 1:S) {
        psi[i] = psi[i] + (R[i-j] * y_reg[i-j] * omega[j]);
      }
    }
  }
}
model{
  y[(D_seed+1):D] ~ poisson(mu[(D_seed+1):D] + psi[(D_seed+1):D]); // likelihood omitting seed days
  r ~ lognormal(mean_log_r, sd_log_r); // reproduction number prior
}
generated quantities{
  vector[D-D_seed] y_rep; // replicated daily case counts
  vector[D-D_seed] log_lik; // log likelihood
  real<lower=0> psi_ahead; // forward burden of infection
  real log_lik_ahead; // forward likelihood

  for (i in (D_seed+1):D) {
    log_lik[i-D_seed] = poisson_lpmf(y[i] | mu[i] + psi[i]);
    y_rep[i-D_seed] = poisson_rng(mu[i] + psi[i]);
  }

  psi_ahead = 0;
  for (i in 1:S) psi_ahead = psi_ahead + (R[D - i + 1] * y_reg[D - i + 1] * omega[i]);
  log_lik_ahead = poisson_lpmf(y_ahead | mu_ahead + psi_ahead);
}

