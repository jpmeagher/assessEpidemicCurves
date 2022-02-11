#' Compute Transmission Quantiles
#'
#' Compute the expected proportion of transmission associated with the q most infectious index cases
#' for a given value of the dispersion parameter k. Assumes Negative-Binomial distributed secondary infections.
#'
#' @param q A an array of quantiles. Must take values on the unit interval. The proportion of most infectious index cases.
#' @param k A positive, real-values scalar. The dispersion parameter for the distribution of secondary infections.
#' @inheritParams convert_gamma_moments
#'
#' @return An array of values on the unit interval. The proportion  of transmission associated with each elemnt of q.
#' @export
#'
#' @examples compute_transmission_quantile(q = 0.2, k = 0.1)
compute_transmission_quantile <- function(
  q, k, perform_checks = TRUE
){
  if (perform_checks) {
    checkmate::assert_numeric(q, lower = 0, upper = 1, any.missing = FALSE)
    checkmate::assert_number(k, lower = 0)
  }
  x_q <- stats::qgamma(1 - q, shape = k, rate = k)
  t_q <- 1 - stats::pgamma(x_q, shape = k + 1, rate = k)
  t_q
}
