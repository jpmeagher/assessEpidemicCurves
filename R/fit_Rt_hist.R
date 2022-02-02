#' Bayesian Estimation of time-varying reproduction numbers
#'
#' @export
#'
#' @param epidemic_curve A vector of integers. The epidemic curve to be
#'   assessed.
#' @param seed_days An integer scalar. The number of days included to seed the
#'   epidemic
#' @param import_rate A vector of integers. The same length as epidemic_curve.
#'   The rate at which cases are imported from outwith the population.
#' @param expected_generation_interval_mean A scalar. The mean of generation
#'   intervals between infector-infectee pairs.
#' @param generation_interval_mean_sd A scalar. The standard deviation of the
#'   mean generation interval between infector-infectee pairs. When set to 0,
#'   assumes that the mean generation interval is a known parameter of the
#'   model. Care should be taken when allowing for an unknown generation
#'   interval as the sampling algorithm can become unstable.
#' @param generation_interval_sd A scalar. The standard deviation of  generation
#'   intervals between infector-infectee pairs.
#' @param generation_interval_length  An integer scalar. The maximum generation
#'   interval.
#' @param bin_width An integer scalar. The bin width defining the histogram
#'   estimator for time-varying reproduction numbers.
#' @param log_Rt_prior_mean A scalar. The mean for a log-normal distribution.
#' @param log_Rt_prior_sd A positive scalar. The standard deviation for a
#'   log-normal distribution.
#' @param expected_k A positive scalar. The expected case dispersion parameter.
#'   Setting to Inf allows for homogeneous disease reproduction.
#' @param log_k_prior_sd A non-negative scalar. The standard deviation of a
#'   log-normal prior on  the case dispersion parameter. The model treats k as a
#'   fixed hyper-parameter when this variable is set to 0. Care should be taken
#'   when allowing for an unknown dispersion parameter as the sampling algorithm
#'   can become unstable.
#' @param ahead Logical. Include a 1 step ahead prediction.
#' @param next_day_cases An integer scalar. One day ahead case count. Allows
#'   model validation.
#' @param next_day_import_rate A positice scalar. Import rate for the day ahead.
#' @param c A small positive integer. A regularising constant.
#' @param nugget A small positive scalar. A regularising constant.
#' @param pars A character vector. Specifies the parameters returned by the
#'   underlying Stan program.
#' @param init Initial values for the sampling scheme. See `rstan::stan` for
#'   details.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @inheritParams fit_Rt_lgp
#' @inheritParams convert_gamma_moments
#'
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#' @examples
#' D <- 30
#' df <- covid_incidence_roi_epidemiological_date[1:D, ]
#' # fit <- fit_Rt_hist(
#' #   epidemic_curve = df$count, seed_days = 5,
#' #   import_rate = rep(1, D),
#' #   generation_interval_mean = 5, generation_interval_sd = 2.5,
#' #   generation_interval_length = 21,
#' #   gp_amplitude = 1, gp_length_scale = 10,
#' #   k = 0.1,
#' #   ahead = TRUE,
#' #   next_day_cases =  covid_incidence_roi_epidemiological_date$count[D+1],
#' #   next_day_import_rate = 1, iter = 1000
#' #   )
fit_Rt_hist <- function(
  epidemic_curve, seed_days,
  import_rate,
  expected_generation_interval_mean,
  generation_interval_mean_sd = 0,
  generation_interval_sd,
  generation_interval_length,
  bin_width = 7,
  log_Rt_prior_mean = 0, log_Rt_prior_sd = 1,
  expected_k = 1, log_k_prior_sd = 0,
  ahead = FALSE,
  next_day_cases = 1, next_day_import_rate = 1,
  c = 1, nugget = 1e-6,
  perform_checks = TRUE,
  pars = c("R", "k", "gi_mean", "y_rep_ahead", "log_lik_ahead"),
  init = function() ifelse(!is.infinite(expected_k), initialise_momentum(epidemic_curve = epidemic_curve, k = expected_k), list()),
  ...) {
  if (perform_checks) {
    checkmate::assert_numeric(epidemic_curve, lower = 0, any.missing = FALSE)
    checkmate::assert_numeric(import_rate, lower = 0, any.missing = FALSE)
    checkmate::assert_number(expected_generation_interval_mean, lower = 0)
    checkmate::assert_number(generation_interval_mean_sd, lower = 0)
    checkmate::assert_number(generation_interval_length, lower = 0)
    checkmate::assert_integerish(bin_width, lower = 1)
    checkmate::assert_number(log_Rt_prior_mean)
    checkmate::assert_number(log_k_prior_sd, lower = 0)
    checkmate::assert_integerish(generation_interval_length, lower = 1, any.missing = FALSE, len = 1)
    checkmate::assert_integerish(seed_days, lower = 1, upper = generation_interval_length - 1, any.missing = FALSE, len = 1)
    checkmate::assert_integerish(length(epidemic_curve), lower = generation_interval_length, any.missing = FALSE, len = 1)
    checkmate::assert_number(expected_k, lower = 0)
    checkmate::assert_true(expected_k != 0)
    checkmate::assert_number(log_k_prior_sd, lower=0)
    checkmate::assert_logical(ahead)
    checkmate::assert_number(next_day_cases, lower=0)
    checkmate::assert_number(next_day_import_rate, lower=0)
    checkmate::assert_number(c, lower=0)
    checkmate::assert_number(nugget, lower=0)
  }
  if (is.infinite(expected_k)) {
    expected_k_inv <- 0 # this is a hack to persuade stan to behave
    log_k_prior_sd <- 0 # no uncertainty on k allowed with homogeneous transmission
  } else {
    expected_k_inv <-  1 / expected_k
  }
  standata <- list(
    N0 = seed_days, N = length(epidemic_curve),
    y = epidemic_curve,
    mu = import_rate,
    S = generation_interval_length,
    expected_generation_interval_mean = expected_generation_interval_mean,
    generation_interval_mean_sd = generation_interval_mean_sd,
    generation_interval_sd = generation_interval_sd,
    delta = bin_width,
    log_Rt_prior_mean = log_Rt_prior_mean,
    log_Rt_prior_mean_sd = log_Rt_prior_sd,
    expected_k_inv = expected_k_inv, log_k_prior_sd = log_k_prior_sd,
    M = as.numeric(ahead),
    y_ahead = next_day_cases, mu_ahead = next_day_import_rate,
    c = c, nugget = nugget
  )
  out <- rstan::sampling(
    stanmodels$hist_Rt, data = standata,
    pars = pars,
    ...
  )
  out
}
