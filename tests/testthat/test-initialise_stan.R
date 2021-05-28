test_that("log-Gaussian process initialisation", {

  checkmate::expect_numeric(
    initialise_lgp_Rt(
      epidemic_curve = covid_incidence_roi_epidemiological_date$count,
      gp_amplitude = 1, k = 1
    )$log_eta,
    any.missing = FALSE, len = nrow(covid_incidence_roi_epidemiological_date)
  )

  checkmate::expect_list(
    initialise_lgp_Rt(
      epidemic_curve = covid_incidence_roi_epidemiological_date$count,
      gp_amplitude = 1, k = 1
    ),
    len = 1, names = "named"
  )

  checkmate::expect_list(
    lapply(
      1:4, function(i) {
        initialise_lgp_Rt(
          epidemic_curve = covid_incidence_roi_epidemiological_date$count,
          gp_amplitude = 1, k = 0.1)
      }),
    len = 4
  )
})
