% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sa_validation.R
\name{sa_validation}
\alias{sa_validation}
\title{1-step ahead validiation}
\usage{
sa_validation(
  epidemic_curve,
  seed_days,
  import_rate,
  generation_interval_mean,
  generation_interval_sd,
  generation_interval_length,
  next_day_cases = NA,
  next_day_import_rate = NA,
  perform_checks = TRUE,
  fit_fun,
  ...
)
}
\arguments{
\item{epidemic_curve}{A vector of integers. The epidemic curve to be
assessed.}

\item{seed_days}{An integer scalar. The number of days included to seed the
epidemic}

\item{import_rate}{A vector of integers. The same length as epidemic_curve.
The rate at which cases are imported from outwith the population.}

\item{generation_interval_mean}{A scalar. The mean of generation intervals
between infector-infectee pairs.}

\item{generation_interval_sd}{A scalar. The standard deviation of  generation
intervals between infector-infectee pairs.}

\item{generation_interval_length}{An integer scalar. The maximum generation
interval.}

\item{next_day_cases}{An integer scalar. One day ahead case count. Allows
model validation.}

\item{next_day_import_rate}{A positice scalar. Import rate for the day ahead.}

\item{perform_checks}{Logical. Should function arguments be checked to ensure
they satisfy any constraints.}

\item{fit_fun}{Function call to fit an epidemic model to data. See
`fit_Rt_lgp` or `fit_Rt_hist`.}

\item{...}{Arguments passed to `fit_fun` including those passed to
`rstan::sampling` (e.g. iter, chains).}
}
\value{
A list with element `elpd`, the estimated elpd, and `y_rep`, a
  sampele from the posterior predictive distribution.
}
\description{
Compute the 1-step ahead elpd for and sample the posterior predictive
distribution for the specified epidemic model.
}
\examples{
  # D <- 33
  # validation_day <- dmy(01102020)
  #
  # df <- covid_incidence_roi_epidemiological_date \%>\%
  #   mutate(ma_count = stats::filter(count, rep(1/7, 7)) \%>\% round)
  #
  # sub_df <- df \%>\%
  #   filter(date >= validation_day - D & date < validation_day)

  # fit <- fit_Rt_lgp(
  #   epidemic_curve = sub_df$ma_count, seed_days = 5,
  #   import_rate = rep(1, D),
  #   generation_interval_mean = 5,
  #   generation_interval_sd = 2.3,
  #   generation_interval_length = 21,
  #   ahead = TRUE,
  #   next_day_cases = df$ma_count[df$date == validation_day],
  #   next_day_import_rate = 1,
  #   seed = 101,
  #   iter = 1000
  # )
}
\seealso{
fit_Rt_lgp fit_Rt_hist
}
