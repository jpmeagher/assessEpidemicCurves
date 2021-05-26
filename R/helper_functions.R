convert_gamma_moments <- function(expectation, standard_deviation){
  shape <- (expectation / standard_deviation)^2
  rate <- expectation / standard_deviation^2

  c("shape" = shape, "rate" = rate)
}
