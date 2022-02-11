//
// Estimation of case R_t.
// Gaussian process prior for R
//
data{
  int<lower=1> N0; // number of days to seed
  int<lower=2> N; // number of days
  int<lower=0> y[N]; // daily count
  real x[N]; // days observed
  vector<lower=0>[N] mu; // proposed daily import rate
  int<lower=1> S; // max generation interval length
  real<lower=0> generation_interval_mean; // generation interval mean
  real<lower=0> generation_interval_sd; // generation interval standard deviation
  real<lower=0> sigma; // GP prior amplitude
  real<lower=0> ls_prior_mean; // GP prior length scale
  real<lower=0> ls_prior_sd; // case dispersion parameter uncertainty
  real<lower=0> k_inv; // inverse case dispersion parameter
  real<lower=0> log_k_prior_sd; // case dispersion parameter uncertainty
  int<lower=0, upper=1> M; // step ahead prediction
  int<lower=0> y_ahead; // forward observations
  real<lower=0> mu_ahead; // forward daily import rate
  real<lower=0> nugget; // non-zero regularising constant
  int<lower=0> c; // case count regularising constant
}
transformed data{
  vector<lower=0>[N] y_reg; // regularised counts
  cov_matrix[N] fixed_C; // f correlation
  matrix[N, N] fixed_L; // cholesky factorisation
  real<lower=0> a; // generation interval shape
  real<lower=0> b; // generation interval rate
  simplex[S] omega; // generation interval pmf
  vector<upper=0>[S] log_omega; // generation interval pmf
  int<lower=0> N_eta; // Momentum days
  int<lower=0, upper=1> uncertain_ls; // Is there uncertainty on ls?
  int<lower=0, upper=1> uncertain_k; // Is there uncertainty on k?
  //  regularise low counts
  for (i in 1:N) {
    y_reg[i] = max(y[i], c);
  }
  // define GP prior correlation
  fixed_C = cov_exp_quad(x, 1, ls_prior_mean);
  for (i in 1:N) fixed_C[i, i] = fixed_C[i, i] + nugget;
  fixed_L = cholesky_decompose(fixed_C);
  // set generation interval
  b = generation_interval_mean / (generation_interval_sd^2);
  a = generation_interval_mean * b;
  omega[1] = gamma_cdf(1.5, a, b);
  for (i in 2:S) omega[i] = gamma_cdf(i+0.5, a, b) - gamma_cdf(i-0.5, a, b);
  omega = omega / sum(omega);
  log_omega = log(omega);
  // Allow for homogeneous case
  if (k_inv == 0) {
    N_eta = 0;
  } else {
    N_eta = N;
  }
  // Allow for prior on ls
  if (ls_prior_sd == 0) {
    uncertain_ls = 0;
  } else {
    uncertain_ls = 1;
  }
  // Allow for prior on k
  if (log_k_prior_sd == 0) {
    uncertain_k = 0;
  } else {
    uncertain_k = 1;
  }
}
parameters{
  vector[N_eta] log_eta; // total daily reproductive rate
  vector[N] epsilon; // underlying Gaussian random variables
  real z_ls[uncertain_ls]; // dispersion parameter random variable
  real z_k[uncertain_k]; // dispersion parameter random variable
}
transformed parameters{
  cov_matrix[N] C; // f correlation
  matrix[N, N] L; // cholesky factorisation
  real<lower=0> ls[uncertain_ls]; // random length-scale
  vector[N] f; // time-dependent log reproductive numbers
  vector<lower=0>[N] R; // time-dependent reproductive numbers
  vector<lower=0>[N] eta; // disease momentum
  matrix[N, S] tmp_psi_matrix; // temporary matrix
  vector<lower=0>[N] psi; // local infection rate
  real<lower=0> k[uncertain_k]; // reported dispersion parameter
  real<lower= 0> kk; // within-algorithm dispersion parameter
  // length-scale
  if (uncertain_ls == 1) {
    ls[1] = ls_prior_mean + (ls_prior_sd * z_ls[1]);
    C = cov_exp_quad(x, 1, ls[1]);
    for (i in 1:N) C[i, i] = C[i, i] + nugget;
    L = cholesky_decompose(C);
  } else {
    C = fixed_C;
    L = fixed_L;
  }
  // reproduction numbers
  f = sigma * L * epsilon;
  R = exp(f);
  // dispersion parameter
  if (uncertain_k == 1) {
    k[1] = exp(log(1 / k_inv) + log_k_prior_sd .* z_k[1]);
    kk = k[1];
  } else {
    kk = 1 / k_inv;
  }
  // disease momentum & burden of infection
  tmp_psi_matrix = rep_matrix(negative_infinity(), N, S);
  tmp_psi_matrix[1, 1] = log(nugget);
  if (k_inv == 0) {
    eta = R .* y_reg;
    for (i in 2:N) {
      if (i <= S) {
        for (j in 1:(i - 1)) {
          tmp_psi_matrix[i, j] = log(eta[i-j]) + log_omega[j];
        }
      } else {
        for (j in 1:S) {
          tmp_psi_matrix[i, j] = log(eta[i-j]) + log_omega[j];
        }
      }
    }
  } else {
    eta = exp(log_eta);
    for (i in 2:N) {
      if (i <= S) {
        for (j in 1:(i - 1)) {
          tmp_psi_matrix[i, j] = log_eta[i-j] + log_omega[j];
        }
      } else {
        for (j in 1:S) {
          tmp_psi_matrix[i, j] = log_eta[i-j] + log_omega[j];
        }
      }
    }
  }
  for(i in 1:N) {
    psi[i] = exp(log_sum_exp(tmp_psi_matrix[i, ]));
  }
}
model{
  target += poisson_lpmf(y[(N0+1):N] | mu[(N0+1):N] + psi[(N0+1):N]); // likelihood
  target += std_normal_lpdf(epsilon); // GP random variates
  if (k_inv != 0) {
    target += gamma_lpdf(exp(log_eta)  | y_reg * kk, kk ./ R);
    target += sum(log_eta); // jacobian correction
      if (uncertain_k == 1) {
        target += std_normal_lpdf(z_k);
      }
  }
  if (uncertain_ls == 1) target += std_normal_lpdf(z_ls);
}
generated quantities{
  vector[N-N0] y_rep;
  vector[N-N0] log_lik; // log likelihood
  real<lower=0> psi_ahead[M]; // forward burden of infection
  real log_lik_ahead[M]; // forward likelihood
  real<lower=0> y_rep_ahead[M]; // forward cases
  // posterior predictive & log likelihood
  for (i in (N0+1):N) {
    y_rep[i-N0] = poisson_rng(mu[i] + psi[i]);
    log_lik[i-N0] = poisson_lpmf(y[i] | mu[i] + psi[i]);
  }
  // 1 step ahead validation
  if (M == 1){
    psi_ahead[M] = 0;
    for (i in 1:S) psi_ahead[M] = psi_ahead[M] + (eta[N - i + 1] * omega[i]);
    log_lik_ahead[M] = poisson_lpmf(y_ahead | mu_ahead + psi_ahead[M]);
    y_rep_ahead[M] = poisson_rng(mu_ahead + psi_ahead[M]);
  }
}
