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
    len = 4, names = "named"
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

  D <- 30
  df <- covid_incidence_roi_epidemiological_date[1:D, ]

  init_list <- lapply(
    1,
    function(i) {
      initialise_lgp_Rt(epidemic_curve = df$count)
    })

  fit <- fit_Rt_lgp(
    epidemic_curve = df$count, seed_days = 5,
    import_rate = rep(1, D),
    generation_interval_mean = 5,
    generation_interval_sd = 2.5,
    generation_interval_length = 21,
    ahead = TRUE,
    next_day_cases =  covid_incidence_roi_epidemiological_date$count[D+1],
    next_day_import_rate = 1,
    iter = 100,
    chains = 1, refresh = 0, show_messages = FALSE,
    init = init_list
  ) %>%
    suppressWarnings()

  init_list <- lapply(
    1,
    function(i) {
      initialise_lgp_Rt(
        epidemic_curve = df$count,
        uncertain_k = TRUE
        )
    })

  fit <- fit_Rt_lgp(
    epidemic_curve = df$count, seed_days = 5,
    import_rate = rep(1, D),
    generation_interval_mean = 5,
    generation_interval_sd = 2.5,
    generation_interval_length = 21,
    log_k_prior_sd = 1,
    ahead = TRUE,
    next_day_cases =  covid_incidence_roi_epidemiological_date$count[D+1],
    next_day_import_rate = 1,
    iter = 100,
    chains = 1, refresh = 0, show_messages = FALSE,
    init = init_list
  ) %>%
    suppressWarnings()

  init_list <- lapply(
    1,
    function(i) {
      initialise_lgp_Rt(
        epidemic_curve = df$count,
        uncertain_k = TRUE, uncertain_length_scale = TRUE
      )
    })

  fit <- fit_Rt_lgp(
    epidemic_curve = df$count, seed_days = 5,
    import_rate = rep(1, D),
    generation_interval_mean = 5,
    generation_interval_sd = 2.5,
    generation_interval_length = 21,
    log_k_prior_sd = 1,
    ls_prior_sd = 1,
    ahead = TRUE,
    next_day_cases =  covid_incidence_roi_epidemiological_date$count[D+1],
    next_day_import_rate = 1,
    iter = 100,
    chains = 1, refresh = 0, show_messages = FALSE,
    init = init_list
  ) %>%
    suppressWarnings()
})

test_that("Momentum initialisation", {

  checkmate::expect_numeric(
    initialise_momentum(
      epidemic_curve = covid_incidence_roi_epidemiological_date$count,k = 1
    )$log_eta,
    any.missing = FALSE, len = nrow(covid_incidence_roi_epidemiological_date)
  )

  checkmate::expect_list(
    initialise_momentum(
      epidemic_curve = covid_incidence_roi_epidemiological_date$count, k = 1
    ),
    len = 1, names = "named"
  )

  checkmate::expect_list(
    lapply(
      1:4, function(i) {
        initialise_momentum(
          epidemic_curve = covid_incidence_roi_epidemiological_date$count, k = 0.1)
      }),
    len = 4
  )
})
