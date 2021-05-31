library(lubridate)
library(dplyr)
library(magrittr)
library(ggplot2)
# library(loo)

test_that("1 step ahead elpd computation", {
  SEED <- 101
  D <- 33
  validation_day <- dmy(01112020)
  df <- covid_incidence_roi_epidemiological_date %>%
    mutate(ma_count = stats::filter(count, rep(1/7, 7)) %>% round)

  sub_df <- df %>%
    filter(date >= validation_day - D & date < validation_day)

  set.seed(SEED)
  fit <- fit_Rt_hist(
    epidemic_curve = sub_df$ma_count, seed_days = 5,
    import_rate = rep(1, D),
    generation_interval_mean = 5,
    generation_interval_sd = 2.5,
    generation_interval_length = 21,
    ahead = TRUE,
    next_day_cases = df$ma_count[df$date == validation_day],
    next_day_import_rate = 1,
    iter = 2000, seed = SEED
  )

  set.seed(SEED)
  sa <- sa_validation(
    epidemic_curve = sub_df$ma_count, seed_days = 5,
    import_rate = rep(1, D),
    generation_interval_mean = 5,
    generation_interval_sd = 2.5,
    generation_interval_length = 21,
    fit_fun = fit_Rt_hist,
    next_day_cases = df$ma_count[df$date == validation_day],
    next_day_import_rate = 1,
    iter = 2000, seed = SEED
  )

  expect_equal(
    unlist(rstan::extract(fit, "y_rep_ahead")),
    sa$y_rep
  )

  expect_equal(
    unlist(rstan::extract(fit, "log_lik_ahead")) %>%
             matrixStats::logSumExp() %>%
             subtract(log(4000)),
    sa$elpd
  )

  # y_scalar <- 3
  # rstan::extract(full_fit, "R") %>%
  #   as.data.frame() %>%
  #   magrittr::extract(sample.int(2000, 250), ) %>%
  #   as.matrix() %>%
  #   t() %>%
  #   data.frame(date = full_df$date, .) %>%
  #   reshape2::melt(id = "date") %>%
  #   ggplot() +
  #   geom_bar(
  #     data = full_df,
  #     aes(x = date, y = y_scalar * ma_count /  max(ma_count)  ), stat = "identity",
  #     alpha = 0.5
  #   ) +
  #   geom_line(aes(x = date, y = value, group = variable), color = "navy", alpha = 0.1) +
  #   geom_hline(yintercept = 1, lty = 3) +
  #   scale_y_continuous(
  #     bquote("R"["t"]),
  #     sec.axis = sec_axis(~ . * max(full_df$ma_count) / y_scalar, name = "7 Day Moving Average Incidence")
  #   ) +
  #   theme_classic()

  # ## Set up LOO-CV for comparison
  # L <- 33
  # first_validation_day <- dmy(01102020)
  # last_validation_day <- dmy(01112020)
  # M <- as.numeric(last_validation_day - first_validation_day)
  #
  # df <- covid_incidence_roi_epidemiological_date %>%
  #   mutate(ma_count = stats::filter(count, rep(1/7, 7)) %>% round)
  #
  # ## Set up LOO-CV for comparison
  # full_df <- df %>%
  #   filter(date >= first_validation_day - L & date < last_validation_day)
  #
  # D <- nrow(full_df)
  # full_fit <- fit_Rt_hist(
  #   epidemic_curve = full_df$ma_count, seed_days = 5,
  #   import_rate = rep(1, D),
  #   generation_interval_mean = 5,
  #   generation_interval_sd = 2.5,
  #   generation_interval_length = 21,
  #   ahead = FALSE,
  #   seed = 101,
  #   iter = 2000
  # )
  #
  # full_loo <- loo(extract_log_lik(full_fit)[, (L + 1 - 5):(D - 5)])
  # full_loo
  #
  # n_samples <- prod(dim(full_fit)[1:2])
  #
  # ## 1 step ahead validation
  # timestamp()
  # i <- 1
  # day <- first_validation_day
  # posterior_predictive <- exact_log_lik <- array(dim = c(n_samples, M))
  # pb <- txtProgressBar(min = 0, max = M, style = 3)
  # while (day < last_validation_day ) {
  #   df_i <- df %>%
  #     filter(date >= day - L & date < day)
  #   fit_i <- fit_Rt_hist(
  #     epidemic_curve = df_i$ma_count, seed_days = 5,
  #     import_rate = rep(1, L),
  #     generation_interval_mean = 5,
  #     generation_interval_sd = 2.5,
  #     generation_interval_length = 21,
  #     ahead = TRUE,
  #     next_day_cases = full_df$ma_count[full_df$date == day],
  #     next_day_import_rate = 1,
  #     seed = 101,
  #     iter = 2000, refresh = 0
  #   )
  #   posterior_predictive[, i] <- unlist(rstan::extract(fit_i, "y_rep_ahead"))
  #   exact_log_lik[, i] <- unlist(rstan::extract(fit_i, "log_lik_ahead"))
  #   day <- day + 1
  #   setTxtProgressBar(pb, i)
  #   i <- i + 1
  # }
  # close(pb)
  # timestamp()
  #
  # elpd_i <- exact_log_lik %>%
  #   apply(2, matrixStats::logSumExp) %>%
  #   subtract(log(n_samples))
  #
  # log_lik_mcse <- exact_log_lik %>%
  #   apply(2, sd) %>%
  #   divide_by(sqrt(n_samples))
  #
  # exact_elpd <- sum(elpd_i)
  # se_elpd <- sqrt( M * var(elpd_i))
  #
  # c(
  #   elpd_loo = full_loo$estimates["elpd_loo", 1], elpd_loo_se = full_loo$estimates["elpd_loo", 2],
  #   elpd_lfo = exact_elpd, elpd_lfo_se = se_elpd
  # )
  #
  # lfo <- list(
  #   elpd_lfo = exact_elpd, se_elpd_lfo = se_elpd,
  #   pointwise = cbind(elpd_loo = elpd_i, se_log_lik = log_lik_mcse),
  # )
})

