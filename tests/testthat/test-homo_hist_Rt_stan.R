library(ggplot2)

test_that("Histogram estimator for Rt with homogeneous reproduction fits to data", {
  D <- 30
  df <- covid_incidence_roi_epidemiological_date[1:D, ]
  fit <- homo_hist_Rt_stan(
    epidemic_curve = df$count, seed_days = 5,
    import_rate = rep(1, D),
    next_day_cases =  covid_incidence_roi_epidemiological_date$count[D+1],
    next_day_import_rate = 1
  )

  expect_true(
    rstan::extract(fit, "R") %>%
      as.data.frame() %>%
      magrittr::extract(, (D - 6): D) %>%
      as.matrix() %>%
      t() %>%
      apply(2, function(x) length(unique(x)) == 1) %>%
      all()
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
  #     aes(x = date, y = y_scalar * count /  max(df$count)  ), stat = "identity",
  #     alpha = 0.5
  #   ) +
  #   geom_line(aes(x = date, y = value, group = variable), color = "navy", alpha = 0.1) +
  #   geom_hline(yintercept = 1, lty = 3) +
  #   scale_y_continuous(
  #     bquote("R"["t"]),
  #     sec.axis = sec_axis(~ . * max(df$count) / y_scalar, name = "7 Day Moving Average Incidence")
  #   ) +
  #   theme_classic()
})
