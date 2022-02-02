#' Initialise sampling scheme.
#'
#' Initialise the sampling scheme for heterogeneous disease reproduction with a
#' log-Gaussian process prior on the time-varying reproduction number.
#'
#'
#' @inheritParams fit_Rt_lgp
#' @param seed A scalar. The random seed.
#'
#' @return A real valued vector. Log disease momentum for each day.
#'
#' @export
initialise_lgp_Rt <- function(
  epidemic_curve, gp_amplitude, k,
  seed = NULL
){
  set.seed(seed)
  D <- length(epidemic_curve)
  R <- exp(stats::rnorm(D, sd = gp_amplitude))
  list(
    log_eta = log(stats::rgamma(D, shape = epidemic_curve * k, rate = k / R))
    )
}

#' Initialise Momentum.
#'
#' Initialise the momentum in the sampling scheme for heterogeneous disease
#' reproduction assuming a reproduction number of 1.
#'
#'
#' @inheritParams fit_Rt_lgp
#' @param seed A scalar. The random seed.
#'
#' @return A real valued vector. Log disease momentum for each day.
#'
#' @export
initialise_momentum <- function(
  epidemic_curve, k,
  seed = NULL
){
  set.seed(seed)
  D <- length(epidemic_curve)
  list(
    log_eta = log(stats::rgamma(D, shape = epidemic_curve * k, rate = k))
  )
}
