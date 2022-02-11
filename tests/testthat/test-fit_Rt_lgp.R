library(ggplot2)

test_that("log-Gaussian process Rt with heterogeneous reproduction fits to data", {
  D <- 30
  df <- covid_incidence_roi_epidemiological_date[1:D, ]
  fit <- fit_Rt_lgp(
    epidemic_curve = df$count, seed_days = 5,
    import_rate = rep(1, D),
    generation_interval_mean = 5,
    generation_interval_sd = 2.5,
    generation_interval_length = 21,
    ahead = TRUE,
    next_day_cases =  covid_incidence_roi_epidemiological_date$count[D+1],
    next_day_import_rate = 1,
    iter = 1000,
    chains = 1, refresh = 0, show_messages = FALSE
  )


  expect_equal(
    rstan::extract(fit, "R") %>%
      as.data.frame() %>%
      ncol(),
    D
  )

  expect_error(
    fit_Rt_lgp(
      epidemic_curve = df$count, seed_days = 21,
      import_rate = rep(1, D), generation_interval_length = 21,
      generation_interval_mean = 5,
      generation_interval_sd = 2.5,
      generation_interval_length = 21,
      next_day_cases =  covid_incidence_roi_epidemiological_date$count[D+1],
      next_day_import_rate = 1
    )
  )

  expect_error(
    fit_Rt_lgp(
      epidemic_curve = df$count, seed_days = 5,
      import_rate = rep(1, D), generation_interval_length = D+1,
      generation_interval_mean = 5,
      generation_interval_sd = 2.5,
      next_day_cases =  covid_incidence_roi_epidemiological_date$count[D+1],
      next_day_import_rate = 1
    )
  )

  # y_scalar <- 5
  # rstan::extract(fit, "R") %>%
  #   as.data.frame() %>%
  #   magrittr::extract(sample.int(2000, 500), ) %>%
  #   as.matrix() %>%
  #   t() %>%
  #   data.frame(date = df$date, .) %>%
  #   reshape2::melt(id = "date") %>%
  #   ggplot() +
  #   geom_bar(
  #     data = df,
  #     aes(x = date, y = y_scalar * count /  max(count)  ), stat = "identity",
  #     alpha = 0.5
  #   ) +
  #   geom_line(aes(x = date, y = value, group = variable), color = "navy", alpha = 0.05) +
  #   geom_hline(yintercept = 1, lty = 3) +
  #   scale_y_continuous(
  #     bquote("R"["t"]),
  #     sec.axis = sec_axis(~ . * max(df$count) / y_scalar, name = "7 Day Moving Average Incidence")
  #   ) +
  #   theme_classic()

  checkmate::expect_number(rstan::extract(fit, "log_lik_ahead") %>% unlist() %>% mean(), na.ok = FALSE)

  fit <- fit_Rt_lgp(
    epidemic_curve = df$count, seed_days = 5,
    import_rate = rep(1, D),
    generation_interval_mean = 5,
    generation_interval_sd = 2.5,
    generation_interval_length = 21,
    iter = 1000,
    chains = 1, refresh = 0, show_messages = FALSE
  )

  expect_equal(
    rstan::extract(fit, "R") %>%
      as.data.frame() %>%
      ncol(),
    D
  )

  # y_scalar <- 5
  # rstan::extract(fit, "R") %>%
  #   as.data.frame() %>%
  #   magrittr::extract(sample.int(2000, 500), ) %>%
  #   as.matrix() %>%
  #   t() %>%
  #   data.frame(date = df$date, .) %>%
  #   reshape2::melt(id = "date") %>%
  #   ggplot() +
  #   geom_bar(
  #     data = df,
  #     aes(x = date, y = y_scalar * count /  max(count)  ), stat = "identity",
  #     alpha = 0.5
  #   ) +
  #   geom_line(aes(x = date, y = value, group = variable), color = "navy", alpha = 0.1) +
  #   geom_hline(yintercept = 1, lty = 3) +
  #   scale_y_continuous(
  #     bquote("R"["t"]),
  #     sec.axis = sec_axis(~ . * max(df$count) / y_scalar, name = "7 Day Moving Average Incidence")
  #   ) +
  #   theme_classic()

  expect_true(is.null(rstan::extract(fit, "log_lik_ahead") %>% unlist()))

  fit <- fit_Rt_lgp(
    epidemic_curve = df$count, seed_days = 5,
    import_rate = rep(1, D), log_k_prior_mean = Inf,
    generation_interval_mean = 5,
    generation_interval_sd = 2.5,
    generation_interval_length = 21,
    iter = 1000, chains = 1, refresh = 0,
    show_messages = FALSE
  )

  expect_equal(
    rstan::extract(fit, "R") %>%
      as.data.frame() %>%
      ncol(),
    D
  )

  # y_scalar <- 5
  # rstan::extract(fit, "R") %>%
  #   as.data.frame() %>%
  #   magrittr::extract(sample.int(2000, 100), ) %>%
  #   as.matrix() %>%
  #   t() %>%
  #   data.frame(date = df$date, .) %>%
  #   reshape2::melt(id = "date") %>%
  #   ggplot() +
  #   geom_bar(
  #     data = df,
  #     aes(x = date, y = y_scalar * count /  max(count)  ), stat = "identity",
  #     alpha = 0.5
  #   ) +
  #   geom_line(aes(x = date, y = value, group = variable), color = "navy", alpha = 0.1) +
  #   geom_smooth(aes(x = date, y = value), color = "navy", lwd = 2) +
  #   geom_hline(yintercept = 1, lty = 3) +
  #   scale_y_continuous(
  #     bquote("R"["t"]),
  #     sec.axis = sec_axis(~ . * max(df$count) / y_scalar, name = "7 Day Moving Average Incidence")
  #   ) +
  #   theme_classic()

  expect_error(
    fit_Rt_lgp(
      epidemic_curve = df$count, seed_days = 5,
      import_rate = rep(1, D), log_k_prior_mean = -Inf,
      generation_interval_mean = 5,
      generation_interval_sd = 2.5,
      generation_interval_length = 21,
      next_day_cases =  covid_incidence_roi_epidemiological_date$count[D+1],
      next_day_import_rate = 1,
      iter = 1000
    )
  )

  # Test stichastic k
  fit <- fit_Rt_lgp(
    epidemic_curve = df$count, seed_days = 5,
    import_rate = rep(1, D), log_k_prior_mean = 0,
    log_k_prior_sd = 1,
    generation_interval_mean = 5,
    generation_interval_sd = 2.5,
    generation_interval_length = 21,
    iter = 1000, chains = 1, refresh = 0,
    show_messages = FALSE
  ) %>%
    suppressWarnings()

  expect_equal(
    rstan::extract(fit, "R") %>%
      as.data.frame() %>%
      ncol(),
    D
  )

  expect_equal(
    rstan::extract(fit, "k") %>%
      unlist %>%
      length(),
    500
  )

  # Test stichastic ls
  fit <- fit_Rt_lgp(
    epidemic_curve = df$count, seed_days = 5,
    import_rate = rep(1, D), log_k_prior_mean = 0,
    generation_interval_mean = 5,
    generation_interval_sd = 2.5,
    generation_interval_length = 21,
    ls_prior_mean = 17.5,
    ls_prior_sd = 1,
    iter = 1000, chains = 1, refresh = 0,
    show_messages = FALSE
  ) %>%
    suppressWarnings()

  expect_equal(
    rstan::extract(fit, "R") %>%
      as.data.frame() %>%
      ncol(),
    D
  )

  expect_equal(
    rstan::extract(fit, "ls") %>%
      unlist %>%
      length(),
    500
  )

  # Test homogeneous
  fit <- fit_Rt_lgp(
    epidemic_curve = df$count, seed_days = 5,
    import_rate = rep(1, D),
    generation_interval_mean = 5,
    generation_interval_sd = 2.5,
    generation_interval_length = 21,
    log_k_prior_mean = Inf,
    ls_prior_mean = 17.5,
    ls_prior_sd = 1,
    iter = 1000, chains = 1, refresh = 0,
    show_messages = FALSE
  ) %>%
    suppressWarnings()

  expect_equal(
    rstan::extract(fit, "R") %>%
      as.data.frame() %>%
      ncol(),
    D
  )

  expect_equal(
    rstan::extract(fit, "ls") %>%
      unlist %>%
      length(),
    500
  )

})
