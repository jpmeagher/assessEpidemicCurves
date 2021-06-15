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
#' @param generation_interval_mean A scalar. The mean of generation intervals
#'   between infector-infectee pairs.
#' @param generation_interval_sd A scalar. The standard deviation of  generation
#'   intervals between infector-infectee pairs.
#' @param generation_interval_length  An integer scalar. The maximum generation
#'   interval.
#' @param bin_width An integer scalar. The bin width defining the histogram
#'   estimator for time-varying reproduction numbers.
#' @param r_prior_mean A scalar. The mean for a log-normal distribution.
#' @param r_prior_sd A positive scalar. The standard deviation for a log-normal
#'   distribution.
#' @param next_day_cases An integer scalar. One day ahead case count. Allows
#'   model validation.
#' @param next_day_import_rate A positice scalar. Import rate for the day ahead.
#' @param c A small positive integer A regularising constant.
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
  generation_interval_mean, generation_interval_sd,
  generation_interval_length,
  bin_width = 7,
  r_prior_mean = 0, r_prior_sd = 1,
  k = 1,
  ahead = FALSE,
  next_day_cases = 1, next_day_import_rate = 1,
  c = 1,
  perform_checks = TRUE,
  ...) {
  if (perform_checks) {
    checkmate::assert_integerish(generation_interval_length, lower = 1, any.missing = FALSE, len = 1)
    checkmate::assert_integerish(seed_days, lower = 1, upper = generation_interval_length - 1, any.missing = FALSE, len = 1)
    checkmate::assert_integerish(length(epidemic_curve), lower = generation_interval_length, any.missing = FALSE, len = 1)
    checkmate::assert_number(k, lower = 0)
    checkmate::assert_true(k != 0)
    checkmate::assert_logical(ahead)
  }
  if (is.infinite(k)) k <- 0 # this is a hack to persuade stan to behave
  standata <- list(
    D = length(epidemic_curve), D_seed = seed_days,
    y = epidemic_curve,
    mu = import_rate,
    S = generation_interval_length,
    mean_omega = generation_interval_mean, sd_omega = generation_interval_sd,
    delta = bin_width,
    mean_log_r = r_prior_mean, sd_log_r = r_prior_sd,
    k = k,
    M = as.numeric(ahead),
    y_ahead = next_day_cases, mu_ahead = next_day_import_rate,
    c = c
  )
  out <- rstan::sampling(
    stanmodels$hist_Rt, data = standata,
    ...
    )
  return(out)
}