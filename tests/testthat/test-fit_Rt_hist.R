library(ggplot2)

test_that("Histogram estimator for Rt with homogeneous reproduction fits to data", {
  N <- 30
  df <- covid_incidence_roi_epidemiological_date[1:N, ]
  fit <- fit_Rt_hist(
    epidemic_curve = df$count, seed_days = 5,
    import_rate = rep(1, N),
    expected_generation_interval_mean = 5,
    generation_interval_sd = 2.5,
    generation_interval_length = 21,
    ahead = TRUE,
    next_day_cases =  covid_incidence_roi_epidemiological_date$count[N+1],
    next_day_import_rate = 1,
    chains = 1, refresh = 0
  )

  expect_true(
    rstan::extract(fit, "R") %>%
      as.data.frame() %>%
      magrittr::extract(, (N - 6): N) %>%
      as.matrix() %>%
      t() %>%
      apply(2, function(x) length(unique(x)) == 1) %>%
      all()
  )

  expect_error(
    fit_Rt_hist(
      epidemic_curve = df$count, seed_days = 21,
      import_rate = rep(1, N),
      expected_generation_interval_mean = 5,
      generation_interval_sd = 2.5,
      generation_interval_length = 21,
      next_day_cases =  covid_incidence_roi_epidemiological_date$count[N+1],
      next_day_import_rate = 1
    )
  )

  expect_error(
    fit_Rt_hist(
      epidemic_curve = df$count, seed_days = 5,
      expected_generation_interval_mean = 5,
      generation_interval_sd = 2.5,
      import_rate = rep(1, N), generation_interval_length = N+1,
      next_day_cases =  covid_incidence_roi_epidemiological_date$count[N+1],
      next_day_import_rate = 1
    )
  )

  checkmate::expect_number(rstan::extract(fit, "log_lik_ahead") %>% unlist() %>% mean(), na.ok = FALSE)

  # y_scalar <- 5
  # rstan::extract(fit, "R") %>%
  #   as.data.frame() %>%
  #   magrittr::extract(sample.int(4000, 500), ) %>%
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

  fit <- fit_Rt_hist(
    epidemic_curve = df$count, seed_days = 5,
    import_rate = rep(1, N),
    expected_generation_interval_mean = 5,
    generation_interval_sd = 2.5,
    generation_interval_length = 21,
    ahead = FALSE,
    next_day_cases =  covid_incidence_roi_epidemiological_date$count[N+1],
    next_day_import_rate = 1,
    chains = 1, refresh = 0
  )

  expect_true(
    rstan::extract(fit, "R") %>%
      as.data.frame() %>%
      magrittr::extract(, (N - 6): N) %>%
      as.matrix() %>%
      t() %>%
      apply(2, function(x) length(unique(x)) == 1) %>%
      all()
  )

  expect_true(is.null(rstan::extract(fit, "log_lik_ahead") %>% unlist()))

  # y_scalar <- 5
  # rstan::extract(fit, "R") %>%
  #   as.data.frame() %>%
  #   magrittr::extract(sample.int(4000, 500), ) %>%
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

  fit <- fit_Rt_hist(
    epidemic_curve = df$count, seed_days = 5,
    import_rate = rep(1, N),
    expected_generation_interval_mean = 5,
    generation_interval_sd = 2.5,
    generation_interval_length = 21,
    expected_k = Inf,
    ahead = FALSE,
    next_day_cases =  covid_incidence_roi_epidemiological_date$count[N+1],
    next_day_import_rate = 1,
    chains = 1, refresh = 0
  )

  expect_true(
    rstan::extract(fit, "R") %>%
      as.data.frame() %>%
      magrittr::extract(, (N - 6): N) %>%
      as.matrix() %>%
      t() %>%
      apply(2, function(x) length(unique(x)) == 1) %>%
      all()
  )

  expect_error(
    fit_Rt_hist(
      epidemic_curve = df$count, seed_days = 5,
      import_rate = rep(1, N), k = -1,
      expected_generation_interval_mean = 5,
      generation_interval_sd = 2.5,
      generation_interval_length = 21,
      next_day_cases =  covid_incidence_roi_epidemiological_date$count[N+1],
      next_day_import_rate = 1,
      iter = 1000
    )
  )

  expect_error(
    fit_Rt_hist(
      epidemic_curve = df$count, seed_days = 5,
      import_rate = rep(1, N), k = 0,
      expected_generation_interval_mean = 5,
      generation_interval_sd = 2.5,
      generation_interval_length = 21,
      next_day_cases =  covid_incidence_roi_epidemiological_date$count[N+1],
      next_day_import_rate = 1,
      iter = 1000
    )
  )

  # y_scalar <- 5
  # rstan::extract(fit, "R") %>%
  #   as.data.frame() %>%
  #   magrittr::extract(sample.int(4000, 500), ) %>%
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

  # fit <- fit_Rt_hist(
  #   epidemic_curve = df$count, seed_days = 5,
  #   import_rate = rep(1, N),
  #   expected_generation_interval_mean = 5,
  #   generation_interval_sd = 2.5,
  #   generation_interval_length = 21,
  #   expected_k = 1, log_k_prior_sd = 1,
  #   ahead = FALSE,
  #   next_day_cases =  covid_incidence_roi_epidemiological_date$count[N+1],
  #   next_day_import_rate = 1,
  #   control = list(adapt_delta = 0.99, max_treedepth = 15),
  #   cores = 4
  # )
  #
  # expect_true(
  #   rstan::extract(fit, "R") %>%
  #     as.data.frame() %>%
  #     magrittr::extract(, (N - 6): N) %>%
  #     as.matrix() %>%
  #     t() %>%
  #     apply(2, function(x) length(unique(x)) == 1) %>%
  #     all()
  # )
  #
  # expect_equal(
  #   rstan::extract(fit, "k") %>%
  #     unlist() %>%
  #     unname() %>%
  #     length(),
  #   4000
  # )

  # y_scalar <- 5
  # rstan::extract(fit, "R") %>%
  #   as.data.frame() %>%
  #   magrittr::extract(sample.int(4000, 500), ) %>%
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

  fit <- fit_Rt_hist(
    epidemic_curve = df$count, seed_days = 5,
    import_rate = rep(1, N),
    expected_generation_interval_mean = 5,
    generation_interval_mean_sd = 0.1,
    generation_interval_sd = 2.5,
    generation_interval_length = 21,
    ahead = TRUE,
    next_day_cases =  covid_incidence_roi_epidemiological_date$count[N+1],
    next_day_import_rate = 1,
    chains = 1, refresh = 0
  )

  expect_true(
    rstan::extract(fit, "R") %>%
      as.data.frame() %>%
      magrittr::extract(, (N - 6): N) %>%
      as.matrix() %>%
      t() %>%
      apply(2, function(x) length(unique(x)) == 1) %>%
      all()
  )

  expect_equal(
    rstan::extract(fit, "gi_mean") %>%
      unlist() %>%
      unname() %>%
      length(),
    1000
  )

  expect_error(
    fit_Rt_hist(
      epidemic_curve = df$count, seed_days = 5,
      import_rate = rep(1, N), k = -1,
      expected_generation_interval_mean = 5,
      generation_interval_sd = 2.5,
      generation_interval_length = 21,
      next_day_cases =  covid_incidence_roi_epidemiological_date$count[N+1],
      next_day_import_rate = 1,
      iter = 1000,
      chains = 1, refresh = 0
    )
  )

  expect_error(
    fit_Rt_hist(
      epidemic_curve = df$count, seed_days = 5,
      import_rate = rep(1, N), k = 0,
      expected_generation_interval_mean = 5,
      generation_interval_sd = 2.5,
      generation_interval_length = 21,
      next_day_cases =  covid_incidence_roi_epidemiological_date$count[N+1],
      next_day_import_rate = 1,
      iter = 1000
    )
  )

})
