## code to prepare `lab_testing` dataset goes here

library(readr)
lab_testing <- read_csv("~/R/Projects/assessing_epidemic_curves/read_only_data/lab_testing.csv")
usethis::use_data(lab_testing, overwrite = TRUE)

roi_covid_test_positivity_rate <- lab_testing %>%
  select(date = Date_HPSC, rate = PosR7) %>%
  mutate(date = strsplit(date, split = " ") %>%
           sapply( function(x) extract2(x, 1))
  ) %>%
  mutate(date = ymd(date))
use_data(roi_covid_test_positivity_rate, overwrite = TRUE)
