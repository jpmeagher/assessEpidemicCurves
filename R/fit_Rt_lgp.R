#' Bayesian Estimation of time-varying reproduction numbers
#'
#' @inheritParams fit_Rt_hist
#' @inheritParams convert_gamma_moments
#' @param gp_amplitude A positive scalar. The amplitude hyper-parameter for the
#'   underlying Gaussian process.
#' @param ls_prior_mean A positive scalar. The length-scale
#'   hyper-parameter for the underlying Gaussian process.
#' @param ls_prior_sd A non-negative scalar. The standard deviation of a
#'   Normal prior on the length-scale hyper-parameter for the underlying
#'   Gaussian process. Treats the length-scale as fixed when set to 0. Care
#'   should be taken when allowing for an unknown length-scale as the sampling
#'   algorithm can become unstable.
#' @param k A positive scalar. The case dispersion parameter. Setting to Inf
#'   allows for homogeneous disease reproduction.
#'
#' @return An object of class `stanfit` returned by `rstan::sampling`
#' @export
#'
#' @seealso homo_hist_Rt_stan
#'
#' @examples
#' D <- 30
#' df <- covid_incidence_roi_epidemiological_date[1:D, ]
#' # fit <- fit_Rt_lgp(
#' #   epidemic_curve = df$count, seed_days = 5,
#' #   import_rate = rep(1, D),
#' #   generation_interval_mean = 5, generation_interval_sd = 2.5,
#' #   generation_interval_length = 21,
#' #   gp_amplitude = 1, ls_prior_mean = 10,
#' #   log_k_prior_mean = log(0.1),
#' #   ahead = TRUE,
#' #   next_day_cases =  covid_incidence_roi_epidemiological_date$count[D+1],
#' #   next_day_import_rate = 1, iter = 1000
#' #   )
fit_Rt_lgp <- function(
  epidemic_curve, seed_days,
  import_rate,
  generation_interval_mean,
  generation_interval_sd,
  generation_interval_length,
  gp_amplitude = 1,
  ls_prior_mean = 17.5,
  ls_prior_sd = 0,
  log_k_prior_mean = 0, log_k_prior_sd = 0,
  ahead = FALSE,
  next_day_cases = 1, next_day_import_rate = 1,
  c = 1, nugget = 1e-6,
  perform_checks = TRUE,
  pars = c("R", "k", "ls", "y_rep", "log_lik", "y_rep_ahead", "log_lik_ahead"),
  ...) {
  if (perform_checks) {
    checkmate::assert_numeric(epidemic_curve, lower = 0, any.missing = FALSE)
    checkmate::assert_numeric(import_rate, lower = 0, any.missing = FALSE)
    checkmate::assert_number(generation_interval_mean, lower = 0)
    checkmate::assert_number(generation_interval_length, lower = 0)
    checkmate::assert_number(gp_amplitude, lower = 0)
    checkmate::assert_number(ls_prior_mean, lower = 0)
    checkmate::assert_number(ls_prior_sd, lower = 0)
    checkmate::assert_number(log_k_prior_sd, lower = 0)
    checkmate::assert_integerish(generation_interval_length, lower = 1, any.missing = FALSE, len = 1)
    checkmate::assert_integerish(seed_days, lower = 1, upper = generation_interval_length - 1, any.missing = FALSE, len = 1)
    checkmate::assert_integerish(length(epidemic_curve), lower = generation_interval_length, any.missing = FALSE, len = 1)
    checkmate::assert_number(log_k_prior_mean)
    checkmate::assert_false(is.infinite(log_k_prior_mean) & log_k_prior_mean < 0)
    checkmate::assert_number(log_k_prior_sd, lower=0)
    checkmate::assert_logical(ahead)
    checkmate::assert_number(next_day_cases, lower=0)
    checkmate::assert_number(next_day_import_rate, lower=0)
    checkmate::assert_number(c, lower=0)
    checkmate::assert_number(nugget, lower=0)
  }
  if (is.infinite(log_k_prior_mean)) {
    k_inv <- 0 # this is a hack to persuade stan to behave
    log_k_prior_sd <- 0 # no uncertainty on k allowed with homogeneous transmission
  } else {
    k_inv <-  1 / exp(log_k_prior_mean)
  }
  standata <- list(
    N0 = seed_days, N = length(epidemic_curve),
    y = epidemic_curve,
    x = seq_along(epidemic_curve),
    mu = import_rate,
    S = generation_interval_length,
    generation_interval_mean = generation_interval_mean,
    generation_interval_sd = generation_interval_sd,
    sigma = gp_amplitude,
    ls_prior_mean = ls_prior_mean,
    ls_prior_sd = ls_prior_sd,
    k_inv = k_inv, log_k_prior_sd = log_k_prior_sd,
    M = as.numeric(ahead),
    y_ahead = next_day_cases, mu_ahead = next_day_import_rate,
    c = c, nugget = nugget
  )
  out <- rstan::sampling(
    stanmodels$lgp_Rt, data = standata,
    pars = pars, ...
  )
  out
}
