#' Bayesian Estimation of time-varying reproduction numbers
#'
#' @export
#'
#' @param epidemic_curve A vector of integers. The epidemic curve to be assessed.
#' @param seed_days An integer scalar. The number of days included to seed the epidemic
#' @param import_rate A vector of integers. The same length as epidemic_curve. The rate at which cases are imported from outwith the population.
#' @param generation_interval_mean A scalar. The mean of generation intervals between infector-infectee pairs.
#' @param generation_interval_sd A scalar. The standard deviation of  generation intervals between infector-infectee pairs.
#' @param generation_interval_length  An integer scalar. The maximum generation interval.
#' @param bin_width An integer scalar. The bin width defining the histogram estimator for time-varying reproduction numbers.
#' @param r_prior_mean A scalar. The mean for a log-normal distribution.
#' @param r_prior_sd A positive scalar. The standard deviation for a log-normal distribution.
#' @param nugget A small positive scalar. A regularising constant.
#' @param c A small positive integer A regularising constant.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#'
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
homo_reproduction_hist_Rt_stan <- function(
  epidemic_curve, seed_days,
  import_rate,
  generation_interval_mean = 5, generation_interval_sd = 2.5,
  generation_interval_length = 21,
  bin_width = 7,
  r_prior_mean = 0, r_prior_sd = 1,
  nugget = 1e-6, c = 1,
  ...) {
  standata <- list(
    D = length(epidemic_curve), D_seed = seed_days,
    mu = import_rate,
    S = generation_interval_length,
    mean_omega = generation_interval_mean, sd_omega = generation_interval_sd,
    delta = bin_width,
    mean_log_r = r_prior_mean, sd_log_r = r_prior_sd,
    nugget = nugget, c = c
  )
  out <- rstan::sampling(stanmodels$homo_reproduction_hist_Rt, data = standata, ...)
  return(out)
}
