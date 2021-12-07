//
// Estimation of case R_t.
// Gaussian process prior for R
//
data{
  int<lower=1> D_seed; // number of days to seed
  int<lower=2> D; // number of days
  int<lower=0> y[D]; // daily count
  real x[D]; // days observed
  vector<lower=0>[D] mu; // proposed daily import rate
  int<lower=1> S; // max generation interval length
  real<lower=0> mean_omega; // generation interval mean
  real<lower=0> sd_omega; // generation interval standard deviation
  real<lower=0> alpha; // GP prior amplitude
  real<lower=0> ell; // GP prior length scale
  real<lower=0> k; // case dispersion parameter
  int<lower=0, upper=1> M; // step ahead prediction
  int<lower=0> y_ahead; // forward observations
  real<lower=0> mu_ahead; // forward daily import rate
  real<lower=0> nugget; // gp regularising constant
  int<lower=0> c; // case count regularising constant
}
transformed data{
  vector<lower=0>[D] y_reg; // regularised counts
  cov_matrix[D] C; // f correlation
  matrix[D, D] L; // cholesky factorisation
  real<lower=0> a; // generation interval shape
  real<lower=0> b; // generation interval rate
  simplex[S] omega; // generation interval pmf
  vector<upper=0>[S] log_omega; // generation interval pmf
  int<lower=0> D_eta; // Momentum days
  //  regularise low counts
  for (i in 1:D) {
    y_reg[i] = max(y[i], c);
  }
  // Define Brownian prior correlation
  C = cov_exp_quad(x, 1, ell);
  for (i in 1:D) C[i, i] = C[i, i] + nugget;
  L = cholesky_decompose(C);
  // set generation interval
  b = mean_omega / (sd_omega^2);
  a = mean_omega * b;
  omega[1] = gamma_cdf(1.5, a, b);
  for (i in 2:S) omega[i] = gamma_cdf(i+0.5, a, b) - gamma_cdf(i-0.5, a, b);
  omega = omega / sum(omega);
  log_omega = log(omega);
  // Allow for homogeneous case
  if (k == 0) {
    D_eta = 0;
  } else {
    D_eta = D;
  }
}
parameters{
  vector[D_eta] log_eta; // total daily reproductive rate
  vector[D] epsilon; // underlying Gaussian random variables
}
transformed parameters{
  vector[D] f; // time-dependent log reproductive numbers
  vector<lower=0>[D] R; // time-dependent reproductive numbers
  vector<lower=0>[D] eta; // Disease momentum
  vector<lower=0>[D] psi; // local infection rate

  // reproduction numbers
  f = alpha * L * epsilon;
  R = exp(f);
  // disease momentum
  if (k == 0) {
    eta = R .* y_reg;
  } else {
    eta = exp(log_eta);
  }
  // burden of infection
  psi = rep_vector(0, D);
  psi[1] = 0;
  for (i in 2:D) {
    if (i <= S) {
      for (j in 1:(i - 1)) {
        psi[i] = psi[i] + (eta[i-j] * omega[j]);
      }
    } else {
      for (j in 1:S) {
        psi[i] = psi[i] + (eta[i-j] * omega[j]);
      }
    }
  }
}
model{
  if (k > 0) {
    target += poisson_lpmf(y[(D_seed+1):D] | mu[(D_seed+1):D] + psi[(D_seed+1):D]); // likelihood
    // y[(D_seed+1):D] ~ poisson(mu[(D_seed+1):D] + psi[(D_seed+1):D]); // likelihood
    target += gamma_lpdf( exp(log_eta)  | y_reg * k, k ./ R ) + sum(log_eta);
    // target += gamma_lpdf(eta  | y_reg * k, k ./ R );
    // eta ~ gamma(y_reg * k, k ./ R);
    // target += sum(log_eta); // jacobian correction
    target += std_normal_lpdf(epsilon); // GP random variates
    // epsilon ~ std_normal(); // GP random variates
  } else {
    target += poisson_lpmf(y[(D_seed+1):D] | mu[(D_seed+1):D] + psi[(D_seed+1):D]); // likelihood
    // y[(D_seed+1):D] ~ poisson(mu[(D_seed+1):D] + psi[(D_seed+1):D]); // likelihood
    target += std_normal_lpdf(epsilon); // GP random variates
    // epsilon ~ std_normal(); // GP random variates
  }
}
generated quantities{
  vector[D-D_seed] log_lik; // log likelihood
  real<lower=0> psi_ahead[M]; // forward burden of infection
  real log_lik_ahead[M]; // forward likelihood
  real<lower=0> y_rep_ahead[M]; // forward cases
  // log likelihood
  for (i in (D_seed+1):D) {
    log_lik[i-D_seed] = poisson_lpmf(y[i] | mu[i] + psi[i]);
  }
  // 1 step ahead validation
  if (M == 1){
    psi_ahead[M] = 0;
    for (i in 1:S) psi_ahead[M] = psi_ahead[M] + (eta[D - i + 1] * omega[i]);
    log_lik_ahead[M] = poisson_lpmf(y_ahead | mu_ahead + psi_ahead[M]);
    y_rep_ahead[M] = poisson_rng(mu_ahead + psi_ahead[M]);
  }
}
