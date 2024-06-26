#' Convert gamma moments
#'
#' Convert the first two moments of a gamma distribution to shape and rate
#' parameters.
#'
#' @param expectation A scalar. The expected value of the gamma random variable.
#' @param standard_deviation A scalar. The standard deviation of the gamma
#'   random variable.
#' @param perform_checks Logical. Should function arguments be checked to ensure
#'   they satisfy any constraints.
#'
#' @return A named vector of length 2. The shape and rate for the gamma random
#'   variable.
#' @export
#' @examples
#' convert_gamma_moments(5, 1)
convert_gamma_moments <- function(expectation, standard_deviation, perform_checks = TRUE){
  if (perform_checks) {
    checkmate::assert_number(expectation, lower = 0, finite = TRUE)
    checkmate::assert_number(standard_deviation, lower = 0, finite = TRUE)
  }
  shape <- (expectation / standard_deviation)^2
  rate <- expectation / (standard_deviation^2)
  c("shape" = shape, "rate" = rate)
}

#' Convert inverse-gamma moments
#'
#' Convert the first two moments of an inverse-gamma distribution to shape and rate
#' parameters.
#'
#' @param expectation A scalar. The expected value of the inverse-gamma random variable.
#' @param standard_deviation A scalar. The standard deviation of the inverse-gamma
#'   random variable.
#' @param perform_checks Logical. Should function arguments be checked to ensure
#'   they satisfy any constraints.
#'
#' @return A named vector of length 2. The shape and rate for the inverse-gamma random
#'   variable.
#' @export
#'
#' @examples
#' convert_inverse_gamma_moments(5, 1)
convert_inverse_gamma_moments <- function(expectation, standard_deviation, perform_checks = TRUE){
  if (perform_checks) {
    checkmate::assert_number(expectation, lower = 0, finite = TRUE)
    checkmate::assert_number(standard_deviation, lower = 0, finite = TRUE)
  }
  shape <- (expectation / standard_deviation)^2 + 2
  rate <- expectation * (shape - 1)
  if (perform_checks) {
    checkmate::assert_number(shape, lower = 0, finite = TRUE)
    checkmate::assert_number(rate, lower = 0, finite = TRUE)
  }
  c("shape" = shape, "rate" = rate)
}
