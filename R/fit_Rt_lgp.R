#' Bayesian Estimation of time-varying reproduction numbers
#'
#' @inheritParams homo_hist_Rt_stan
#' @inheritParams convert_gamma_moments
#' @param gp_amplitude A positive scalar. The amplitude hyper-parameter for the
#'   underlying Gaussian process.
#' @param gp_length_scale A positive scalar. The length-scale hyper-parameter
#'   for the underlying Gaussian process.
#' @param k A positive scalar. The case dispersion parameter. Setting to Inf
#'   allows for homogeneous disease reproduction.
#' @param nugget A small positive scalar. A regularising constant.
#'
#' @return An object of class `stanfit` returned by `rstan::sampling`
#' @export
#'
#' @seealso homo_hist_Rt_stan
#'
#' @examples
#' D <- 30
#' df <- covid_incidence_roi_epidemiological_date[1:D, ]
#' fit <- fit_Rt_lgp(
#'   epidemic_curve = df$count, seed_days = 5,
#'   import_rate = rep(1, D),
#'   generation_interval_mean = 5, generation_interval_sd = 2.5,
#'   generation_interval_length = 21,
#'   gp_amplitude = 1, gp_length_scale = 10,
#'   k = 0.1,
#'   next_day_cases =  covid_incidence_roi_epidemiological_date$count[D+1],
#'   next_day_import_rate = 1, iter = 1000
#'   )
fit_Rt_lgp <- function(
  epidemic_curve, seed_days,
  import_rate,
  generation_interval_mean = 5, generation_interval_sd = 2.5,
  generation_interval_length = 21,
  gp_amplitude = 1, gp_length_scale = 10,
  k = 0.1,
  next_day_cases = NA, next_day_import_rate = NA,
  nugget = 1e-6, c = 1,
  perform_checks = TRUE,
  ...) {
  if (perform_checks) {
    checkmate::assert_integerish(generation_interval_length, lower = 1, any.missing = FALSE, len = 1)
    checkmate::assert_integerish(seed_days, lower = 1, upper = generation_interval_length - 1, any.missing = FALSE, len = 1)
    checkmate::assert_integerish(length(epidemic_curve), lower = generation_interval_length, any.missing = FALSE, len = 1)
    checkmate::assert_number(k, lower = 0)
  }
  standata <- list(
    D = length(epidemic_curve), D_seed = seed_days,
    y = epidemic_curve,
    x = seq_along(epidemic_curve),
    mu = import_rate,
    S = generation_interval_length,
    mean_omega = generation_interval_mean, sd_omega = generation_interval_sd,
    alpha = gp_amplitude, ell = gp_length_scale,
    k = k,
    y_ahead = next_day_cases, mu_ahead = next_day_import_rate,
    nugget = nugget, c = c
  )
  if (is.na(next_day_cases)) {
    if (is.infinite(k)) {
      out <- rstan::sampling(stanmodels$lgp_Rt_homo, data = standata, ...)
    } else {
      out <- rstan::sampling(stanmodels$lgp_Rt_fixed_k, data = standata, ...)
    }
  } else {
    if (is.infinite(k)) {
      out <- rstan::sampling(stanmodels$lgp_Rt_homo_ahead, data = standata, ...)
    } else {
      out <- rstan::sampling(stanmodels$lgp_Rt_fixed_k_ahead, data = standata, ...)
    }
  }

  return(out)
}
