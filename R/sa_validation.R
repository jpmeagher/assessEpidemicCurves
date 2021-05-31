#' 1-step ahead validiation
#'
#' Compute the 1-step ahead elpd for and sample the posterior predictive
#' distribution for the specified epidemic model.
#'
#' @inheritParams fit_Rt_lgp
#' @inheritParams fit_Rt_hist
#' @inheritParams convert_gamma_moments
#' @param fit_fun Function call to fit an epidemic model to data. See
#'   `fit_Rt_lgp` or `fit_Rt_hist`.
#' @param ... Arguments passed to `fit_fun` including those passed to
#'   `rstan::sampling` (e.g. iter, chains).
#'
#' @seealso fit_Rt_lgp fit_Rt_hist
#'
#' @return A list with element `elpd`, the estimated elpd, and `y_rep`, a
#'   sampele from the posterior predictive distribution.
#' @export
#'
#' @examples
#'   # D <- 33
#'   # validation_day <- dmy(01102020)
#'   #
#'   # df <- covid_incidence_roi_epidemiological_date %>%
#'   #   mutate(ma_count = stats::filter(count, rep(1/7, 7)) %>% round)
#'   #
#'   # sub_df <- df %>%
#'   #   filter(date >= validation_day - D & date < validation_day)
#'
#'   # fit <- fit_Rt_lgp(
#'   #   epidemic_curve = sub_df$ma_count, seed_days = 5,
#'   #   import_rate = rep(1, D),
#'   #   generation_interval_mean = 5,
#'   #   generation_interval_sd = 2.3,
#'   #   generation_interval_length = 21,
#'   #   ahead = TRUE,
#'   #   next_day_cases = df$ma_count[df$date == validation_day],
#'   #   next_day_import_rate = 1,
#'   #   seed = 101,
#'   #   iter = 1000
#'   # )
sa_validation <- function(
  epidemic_curve, seed_days,
  import_rate,
  generation_interval_mean, generation_interval_sd,
  generation_interval_length,
  next_day_cases = NA, next_day_import_rate = NA,
  perform_checks = TRUE,
  fit_fun,
  ...
){
  fit <- fit_fun(
    epidemic_curve = epidemic_curve,
    seed_days = seed_days, import_rate = import_rate,
    generation_interval_mean = generation_interval_mean,
    generation_interval_sd = generation_interval_sd,
    generation_interval_length = generation_interval_length,
    ahead = TRUE,
    next_day_cases = next_day_cases,
    next_day_import_rate = next_day_import_rate,
    perform_checks = perform_checks,
    ...
  )
  log_lik <- unlist(rstan::extract(fit, "log_lik_ahead"))
  elpd <- matrixStats::logSumExp(log_lik) - log(length(log_lik))
  y_rep <- unlist(rstan::extract(fit, "y_rep_ahead"))
  list(elpd = elpd, y_rep = y_rep)
}
