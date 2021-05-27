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

test_that("Inverse-gamma converter works", {
  alpha <- 20
  beta <- 10
  mu <- beta / (alpha - 1)
  sigma <- sqrt((beta^2) / ((alpha - 1)^2 * (alpha - 2)))

  expect_error(convert_inverse_gamma_moments(expectation = -1, standard_deviation = 1))
  expect_error(convert_inverse_gamma_moments(expectation = 1, standard_deviation = -1))
  expect_named(convert_inverse_gamma_moments(expectation = mu, standard_deviation = sigma), expected = c("shape", "rate"))
  expect_equal(convert_inverse_gamma_moments(expectation = mu, standard_deviation = sigma)[1], alpha, ignore_attr = TRUE)
  expect_equal(convert_inverse_gamma_moments(expectation = mu, standard_deviation = sigma)[2], beta, ignore_attr = TRUE)

  convert_inverse_gamma_moments(expectation = mu, standard_deviation = sigma)
})
