library(checkmate)

test_that("transmission quantiles", {
  q <- 0.2
  k <- 1
  expect_equal(
    compute_transmission_quantile(q, k),
    1 - pgamma(qgamma(1-q, k, k), k+1, k)
  )
  expect_numeric(
    compute_transmission_quantile(runif(10), 1)
  )

  expect_matrix(
    compute_transmission_quantile(matrix(runif(25), 5, 5), k / 10)
  )

})
