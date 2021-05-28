library(ggplot2)

test_that("log-Gaussian process Rt with heterogeneous reproduction fits to data", {
  D <- 30
  df <- covid_incidence_roi_epidemiological_date[1:D, ]
  fit <- fit_Rt_lgp(
    epidemic_curve = df$count, seed_days = 5,
    import_rate = rep(1, D),
    ahead = TRUE,
    next_day_cases =  covid_incidence_roi_epidemiological_date$count[D+1],
    next_day_import_rate = 1,
    iter = 1000
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
      next_day_cases =  covid_incidence_roi_epidemiological_date$count[D+1],
      next_day_import_rate = 1
    )
  )

  expect_error(
    fit_Rt_lgp(
      epidemic_curve = df$count, seed_days = 5,
      import_rate = rep(1, D), generation_interval_length = D+1,
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
    iter = 1000
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
    import_rate = rep(1, D), k = Inf,
    iter = 1000
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
      import_rate = rep(1, D), k = -1,
      next_day_cases =  covid_incidence_roi_epidemiological_date$count[D+1],
      next_day_import_rate = 1,
      iter = 1000
    )
  )

  expect_error(
    fit_Rt_lgp(
      epidemic_curve = df$count, seed_days = 5,
      import_rate = rep(1, D), k = 0,
      next_day_cases =  covid_incidence_roi_epidemiological_date$count[D+1],
      next_day_import_rate = 1,
      iter = 1000
    )
  )
})
