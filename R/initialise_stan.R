#' Initialise sampling scheme.
#'
#' Initialise the sampling scheme for heterogeneous disease reproduction with a
#' log-Gaussian process prior on the time-varying reproduction number.
#'
#'
#' @inheritParams fit_Rt_lgp
#' @param k A positive scalar. A value for the dispersion parameter k that
#'   allows the momentum parameters to be initialised.
#' @param uncertain_k Logical. Is k treated as a random quantity in this model.
#' @param uncertain_length_scale Logical. Is the Gaussian process length_-cale treated as
#'   a random quantity in this model.
#' @param seed A scalar. The random seed.
#'
#' @return A real valued vector. Log disease momentum for each day.
initialise_lgp_Rt <- function(
  epidemic_curve, gp_amplitude = 0.1, k = 0.1,
  uncertain_k = FALSE, uncertain_length_scale = FALSE,
  seed = NULL
){
  set.seed(seed)
  D <- length(epidemic_curve)
  R <- exp(stats::rnorm(D, sd = gp_amplitude))
  init <- list(
    log_eta = log(stats::rgamma(D, shape = epidemic_curve * k, rate = k / R)),
    epsilon = stats::runif(D, min = -2, max = 2),
    z_ls = as.array(stats::runif(uncertain_length_scale, min = -2, max = 2)),
    z_k = as.array(stats::runif(uncertain_k, min = -2, max = 2))
    )
  init
}

#' Initialise Momentum.
#'
#' Initialise the momentum in the sampling scheme for heterogeneous disease
#' reproduction assuming a reproduction number of 1.
#'
#' @inheritParams fit_Rt_lgp
#' @param log_R_sd A positive scalar. The standard deviation of the zero-mean
#'   log-normal distribution providing initial time-varying reproduction
#'   numbers.
#' @param k A positive scalar. The value for the dispersion parameter k from
#'   which initial momentum parameters are sampled.
#' @param seed A scalar. The random seed.
#'
#' @return A real valued vector. Log disease momentum for each day.
#'
#' @export
initialise_momentum <- function(
  epidemic_curve, log_R_sd = 0.1, k = 0.1,
  seed = NULL
){
  set.seed(seed)
  D <- length(epidemic_curve)
  R <- exp(stats::rnorm(D, sd = log_R_sd))
  list(
    log_eta = log(stats::rgamma(D, shape = epidemic_curve * k, rate = k / R))
  )
}
