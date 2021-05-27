#' Bayesian Estimation of time-varying reproduction numbers
#'
#' @inheritParams homo_hist_Rt_stan
#' @param gp_amplitude A positive scalar. The amplitude hyper-parameter for the
#'   underlying Gaussian process.
#' @param gp_length_scale A positive scalar. The length-scale hyper-parameter
#'   for the underlying Gaussian process.
#' @param k A positive scalar. The case dispersion parameter.
#' @param nugget A small positive scalar. A regularising constant.
#'
#' @return
#' @export
#'
#' @seealso homo_hist_Rt_stan
#'
#' @examples
#' D <- 30
#' df <- covid_incidence_roi_epidemiological_date[1:D, ]
#' fit <- hetero_fixed_k_lgp_Rt_stan(
#'   epidemic_curve = df$count, seed_days = 5,
#'   import_rate = rep(1, D),
#'   next_day_cases =  covid_incidence_roi_epidemiological_date$count[D+1],
#'   next_day_import_rate = 1
#'   )
hetero_fixed_k_lgp_Rt_stan <- function(
  epidemic_curve, seed_days,
  import_rate,
  generation_interval_mean = 5, generation_interval_sd = 2.5,
  generation_interval_length = 21,
  gp_amplitude = 1, gp_length_scale = 10,
  k = 0.1,
  next_day_cases, next_day_import_rate,
  nugget = 1e-6, c = 1,
  ...) {
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
  out <- rstan::sampling(stanmodels$hetero_fixed_k_lgp_Rt, data = standata, ...)
  return(out)
}
