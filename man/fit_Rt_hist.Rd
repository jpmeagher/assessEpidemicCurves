% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_Rt_hist.R
\name{fit_Rt_hist}
\alias{fit_Rt_hist}
\title{Bayesian Estimation of time-varying reproduction numbers}
\usage{
fit_Rt_hist(
  epidemic_curve,
  seed_days,
  import_rate,
  generation_interval_mean,
  generation_interval_sd,
  generation_interval_length,
  bin_width = 7,
  log_Rt_prior_mean = 0,
  log_Rt_prior_sd = 1,
  log_k_prior_mean = 0,
  log_k_prior_sd = 0,
  ahead = FALSE,
  next_day_cases = 1,
  next_day_import_rate = 1,
  c = 1,
  nugget = 1e-06,
  perform_checks = TRUE,
  pars = c("R", "k", "y_rep", "log_lik", "y_rep_ahead", "log_lik_ahead"),
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

\item{bin_width}{An integer scalar. The bin width defining the histogram
estimator for time-varying reproduction numbers.}

\item{log_Rt_prior_mean}{A scalar. The mean for a log-normal distribution.}

\item{log_Rt_prior_sd}{A positive scalar. The standard deviation for a
log-normal distribution.}

\item{log_k_prior_mean}{A real-valued scalar. The prior mean of the log case
dispersion parameter. Setting to Inf allows for homogeneous disease
reproduction.}

\item{log_k_prior_sd}{A non-negative scalar. The standard deviation of a
log-normal prior on  the case dispersion parameter. The model treats k as a
fixed hyper-parameter when this variable is set to 0. Care should be taken
when allowing for an unknown dispersion parameter as the sampling algorithm
can become unstable.}

\item{ahead}{Logical. Include a 1 step ahead prediction.}

\item{next_day_cases}{An integer scalar. One day ahead case count. Allows
model validation.}

\item{next_day_import_rate}{A positice scalar. Import rate for the day ahead.}

\item{c}{A small positive integer. A regularising constant.}

\item{nugget}{A small positive scalar. A regularising constant.}

\item{perform_checks}{Logical. Should function arguments be checked to ensure
they satisfy any constraints.}

\item{pars}{A character vector. Specifies the parameters returned by the
underlying Stan program.}

\item{...}{Arguments passed to `rstan::sampling` (e.g. iter, chains).}
}
\value{
An object of class `stanfit` returned by `rstan::sampling`
}
\description{
Bayesian Estimation of time-varying reproduction numbers
}
\examples{
D <- 30
df <- covid_incidence_roi_epidemiological_date[1:D, ]
# fit <- fit_Rt_hist(
#   epidemic_curve = df$count, seed_days = 5,
#   import_rate = rep(1, D),
#   generation_interval_mean = 5, generation_interval_sd = 2.5,
#   generation_interval_length = 21,
#   gp_amplitude = 1, gp_length_scale = 10,
#   log_k_prior_mean = log(1),
#   ahead = TRUE,
#   next_day_cases =  covid_incidence_roi_epidemiological_date$count[D+1],
#   next_day_import_rate = 1, iter = 1000
#   )
}
