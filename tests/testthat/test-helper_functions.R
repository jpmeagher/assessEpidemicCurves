test_that("Gamma converter works", {
  alpha <- 5
  beta <- 1
  mu <- alpha / beta
  sigma <- sqrt(alpha / (beta^2))

  expect_error(convert_gamma_moments(expectation = -1, standard_deviation = 1))
  expect_error(convert_gamma_moments(expectation = 1, standard_deviation = -1))
  expect_named(convert_gamma_moments(expectation = mu, standard_deviation = sigma), expected = c("shape", "rate"))
  expect_equal(convert_gamma_moments(expectation = mu, standard_deviation = sigma)[1], alpha, ignore_attr = TRUE)
  expect_equal(convert_gamma_moments(expectation = mu, standard_deviation = sigma)[2], beta, ignore_attr = TRUE)
})
