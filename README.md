
<!-- README.md is generated from README.Rmd. Please edit that file -->

# assessEpidemicCurves

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/jpmeagher/assessEpidemicCurves.svg?branch=master)](https://travis-ci.com/jpmeagher/assessEpidemicCurves)
<!-- badges: end -->

The goal of assessEpidemicCurves is to model epidemic curves under
heterogeneous disease reproduction, providing estimates for the
time-varying reproduction number and assessing epidemic curves for
evidence of superspreading. This package accompanies \[Reference our
paper here\]

## Installation

You can install the development version of assessEpidemicCurves from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jpmeagher/assessEpidemicCurves")
```

## Case Study: COVID-19 in the Republic of Ireland

``` r
library(assessEpidemicCurves)
library(rstan)
library(loo)
library(ggplot2)
library(dplyr)
library(magrittr)
library(lubridate)
library(EpiEstim)
library(knitr)
library(reshape2)
library(HDInterval)
library(latex2exp)
library(matrixStats)
```

As a case study we include the epidemic curve of confirmed COVID-19
cases in the Republic of Ireland from March 1, 2020 to 28 February,
2021. Cases are ordered by epidemiological date. This is either the date
of onset of symptoms, date of diagnosis, laboratory specimen collection
date, laboratory received date, laboratory reported date or the
notification date.

<img src="man/figures/README-raw_epidemic_curve-1.png" width="100%" />

## The generative model for epidemic curves

Analyses are based on a hierarchical model for the daily incidence count
of COVID-19

![
\\begin{aligned}
y\_t \\mid \\mu\_t, \\boldsymbol \\eta\_t, \\boldsymbol \\omega &\\sim \\operatorname{Pois} \\left( \\mu\_t + \\sum\_{s = 1}^t \\omega\_s \\eta\_{t-s} \\right), \\\\
\\eta\_t \\mid k, R\_t &\\sim \\operatorname{Gamma} \\left( y\_t k, \\frac{k}{R\_t}\\right),
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0Ay_t%20%5Cmid%20%5Cmu_t%2C%20%5Cboldsymbol%20%5Ceta_t%2C%20%5Cboldsymbol%20%5Comega%20%26%5Csim%20%5Coperatorname%7BPois%7D%20%5Cleft%28%20%5Cmu_t%20%2B%20%5Csum_%7Bs%20%3D%201%7D%5Et%20%5Comega_s%20%5Ceta_%7Bt-s%7D%20%5Cright%29%2C%20%5C%5C%0A%5Ceta_t%20%5Cmid%20k%2C%20R_t%20%26%5Csim%20%5Coperatorname%7BGamma%7D%20%5Cleft%28%20y_t%20k%2C%20%5Cfrac%7Bk%7D%7BR_t%7D%5Cright%29%2C%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
y_t \mid \mu_t, \boldsymbol \eta_t, \boldsymbol \omega &\sim \operatorname{Pois} \left( \mu_t + \sum_{s = 1}^t \omega_s \eta_{t-s} \right), \\
\eta_t \mid k, R_t &\sim \operatorname{Gamma} \left( y_t k, \frac{k}{R_t}\right),
\end{aligned}
")

where
![\\boldsymbol y = \\left(y\_0, y\_1, \\dots, y\_N \\right)^\\top](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20y%20%3D%20%5Cleft%28y_0%2C%20y_1%2C%20%5Cdots%2C%20y_N%20%5Cright%29%5E%5Ctop "\boldsymbol y = \left(y_0, y_1, \dots, y_N \right)^\top")
is the epidemic curve seeded by
![y\_0](https://latex.codecogs.com/png.latex?y_0 "y_0"),
![\\boldsymbol \\mu = \\left(\\mu\_0, \\mu\_1, \\dots, \\mu\_N \\right)^\\top](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Cmu%20%3D%20%5Cleft%28%5Cmu_0%2C%20%5Cmu_1%2C%20%5Cdots%2C%20%5Cmu_N%20%5Cright%29%5E%5Ctop "\boldsymbol \mu = \left(\mu_0, \mu_1, \dots, \mu_N \right)^\top")
is the rate at which cases are imported,
![\\boldsymbol \\omega = \\left(\\omega\_1, \\omega\_2, \\dots \\right)^\\top](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Comega%20%3D%20%5Cleft%28%5Comega_1%2C%20%5Comega_2%2C%20%5Cdots%20%5Cright%29%5E%5Ctop "\boldsymbol \omega = \left(\omega_1, \omega_2, \dots \right)^\top")
is the generation interval pmf,
![\\boldsymbol \\eta\_t = \\left( \\eta\_0, \\eta\_1, \\dots, \\eta\_{t - 1} \\right)](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Ceta_t%20%3D%20%5Cleft%28%20%5Ceta_0%2C%20%5Ceta_1%2C%20%5Cdots%2C%20%5Ceta_%7Bt%20-%201%7D%20%5Cright%29 "\boldsymbol \eta_t = \left( \eta_0, \eta_1, \dots, \eta_{t - 1} \right)")
is the disease momentum up to day
![t](https://latex.codecogs.com/png.latex?t "t"),
![k](https://latex.codecogs.com/png.latex?k "k") is the case dispersion
parameter and
![\\boldsymbol R = \\left(R\_0, R\_1, \\dots, R\_N \\right)^\\top](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20R%20%3D%20%5Cleft%28R_0%2C%20R_1%2C%20%5Cdots%2C%20R_N%20%5Cright%29%5E%5Ctop "\boldsymbol R = \left(R_0, R_1, \dots, R_N \right)^\top")
are the cohort reproduction numbers.

This model is underpinned by a branching process model for secondary
infections such that

![
\\begin{aligned}
z\_{t, i} \\mid \\nu\_{t\_i} &\\sim \\operatorname{Pois} \\left( \\nu\_{t, i} \\right), \\\\
\\nu\_{t, i} \\mid k, R\_t &\\sim \\operatorname{Gamma} \\left( k, \\frac{k}{R\_t}\\right),
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0Az_%7Bt%2C%20i%7D%20%5Cmid%20%5Cnu_%7Bt_i%7D%20%26%5Csim%20%5Coperatorname%7BPois%7D%20%5Cleft%28%20%5Cnu_%7Bt%2C%20i%7D%20%5Cright%29%2C%20%5C%5C%0A%5Cnu_%7Bt%2C%20i%7D%20%5Cmid%20k%2C%20R_t%20%26%5Csim%20%5Coperatorname%7BGamma%7D%20%5Cleft%28%20k%2C%20%5Cfrac%7Bk%7D%7BR_t%7D%5Cright%29%2C%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
z_{t, i} \mid \nu_{t_i} &\sim \operatorname{Pois} \left( \nu_{t, i} \right), \\
\nu_{t, i} \mid k, R_t &\sim \operatorname{Gamma} \left( k, \frac{k}{R_t}\right),
\end{aligned}
")

where
![z\_{t, i}](https://latex.codecogs.com/png.latex?z_%7Bt%2C%20i%7D "z_{t, i}")
is the number of secondary infections arising from the
![i^{th}](https://latex.codecogs.com/png.latex?i%5E%7Bth%7D "i^{th}")
index case on day ![t](https://latex.codecogs.com/png.latex?t "t") with
an individual reproduction number
![\\nu\_{t, i}](https://latex.codecogs.com/png.latex?%5Cnu_%7Bt%2C%20i%7D "\nu_{t, i}").
Within this model,
![\\mathbb E \\left\[ z\_{t, i}\\right\] = R\_t](https://latex.codecogs.com/png.latex?%5Cmathbb%20E%20%5Cleft%5B%20z_%7Bt%2C%20i%7D%5Cright%5D%20%3D%20R_t "\mathbb E \left[ z_{t, i}\right] = R_t")
and
![\\operatorname{Var} \\left( z\_{t, i}\\right) = R\_t + R\_t^2 / k](https://latex.codecogs.com/png.latex?%5Coperatorname%7BVar%7D%20%5Cleft%28%20z_%7Bt%2C%20i%7D%5Cright%29%20%3D%20R_t%20%2B%20R_t%5E2%20%2F%20k "\operatorname{Var} \left( z_{t, i}\right) = R_t + R_t^2 / k").
This links to the model for
![y\_t](https://latex.codecogs.com/png.latex?y_t "y_t") via

![
\\eta\_t = \\sum\_{i = 1}^{y\_t} \\nu\_{t, i}.
](https://latex.codecogs.com/png.latex?%0A%5Ceta_t%20%3D%20%5Csum_%7Bi%20%3D%201%7D%5E%7By_t%7D%20%5Cnu_%7Bt%2C%20i%7D.%0A "
\eta_t = \sum_{i = 1}^{y_t} \nu_{t, i}.
")

Thus, the hierarchical model allows for heterogeneous disease
reproduction resulting in over dispersed distributions of secondary
infections, that is, superspreading may be a feature of the epidemic.
Note that
![\\lim\_{k \\to \\infty} \\operatorname{Var} \\left( \\nu\_{t, i}\\right) = 0](https://latex.codecogs.com/png.latex?%5Clim_%7Bk%20%5Cto%20%5Cinfty%7D%20%5Coperatorname%7BVar%7D%20%5Cleft%28%20%5Cnu_%7Bt%2C%20i%7D%5Cright%29%20%3D%200 "\lim_{k \to \infty} \operatorname{Var} \left( \nu_{t, i}\right) = 0")
and so this model admits homogeneous disease reproduction as a special
case.

## Prior specification

In order to fit this generative model to data, we fix import rates
![\\boldsymbol \\mu = \\left(\\mu\_0, \\mu\_1, \\dots, \\mu\_N \\right)^\\top](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Cmu%20%3D%20%5Cleft%28%5Cmu_0%2C%20%5Cmu_1%2C%20%5Cdots%2C%20%5Cmu_N%20%5Cright%29%5E%5Ctop "\boldsymbol \mu = \left(\mu_0, \mu_1, \dots, \mu_N \right)^\top"),
the generation interval pmf
![\\boldsymbol \\omega = \\left(\\omega\_1, \\omega\_2, \\dots \\right)^\\top](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Comega%20%3D%20%5Cleft%28%5Comega_1%2C%20%5Comega_2%2C%20%5Cdots%20%5Cright%29%5E%5Ctop "\boldsymbol \omega = \left(\omega_1, \omega_2, \dots \right)^\top"),
and case dispersion parameter
![k](https://latex.codecogs.com/png.latex?k "k"). We then propose a
Gaussian process prior for
![\\boldsymbol R](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20R "\boldsymbol R")
such that

![
\\begin{aligned}
\\log R\_t &=  f \\left( t \\right), \\\\
f \\left( t \\right) &\\sim \\mathcal{GP} \\left( 0, k \\left(t, t' \\right) \\right), \\\\
k \\left(t, t' \\right) & = \\alpha^2 \\exp \\left( - \\frac{\\left(t - t'\\right)^2}{2 \\ell^2}\\right),
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0A%5Clog%20R_t%20%26%3D%20%20f%20%5Cleft%28%20t%20%5Cright%29%2C%20%5C%5C%0Af%20%5Cleft%28%20t%20%5Cright%29%20%26%5Csim%20%5Cmathcal%7BGP%7D%20%5Cleft%28%200%2C%20k%20%5Cleft%28t%2C%20t%27%20%5Cright%29%20%5Cright%29%2C%20%5C%5C%0Ak%20%5Cleft%28t%2C%20t%27%20%5Cright%29%20%26%20%3D%20%5Calpha%5E2%20%5Cexp%20%5Cleft%28%20-%20%5Cfrac%7B%5Cleft%28t%20-%20t%27%5Cright%29%5E2%7D%7B2%20%5Cell%5E2%7D%5Cright%29%2C%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
\log R_t &=  f \left( t \right), \\
f \left( t \right) &\sim \mathcal{GP} \left( 0, k \left(t, t' \right) \right), \\
k \left(t, t' \right) & = \alpha^2 \exp \left( - \frac{\left(t - t'\right)^2}{2 \ell^2}\right),
\end{aligned}
")

where amplitude
![\\alpha &gt; 0](https://latex.codecogs.com/png.latex?%5Calpha%20%3E%200 "\alpha > 0")
and length-scale
![\\ell &gt; 0](https://latex.codecogs.com/png.latex?%5Cell%20%3E%200 "\ell > 0")
are specified a priori.

## Model fitting

The model described here can be computationally expensive to implement
and so we restrict our analysis to a subset of the data. We examine
reported COVID-19 cases from December 10, 2020, to January 31 2020,
allowing the 5 days from December 5 to December 9 inclusive to seed the
epidemic.

We assume that
![\\mu\_t = 1](https://latex.codecogs.com/png.latex?%5Cmu_t%20%3D%201 "\mu_t = 1")
for all ![t](https://latex.codecogs.com/png.latex?t "t") and that
![\\boldsymbol \\omega](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Comega "\boldsymbol \omega")
is proportional to a discretised gamma distribution such that generation
intervals have a mean of 5 days, standard deviation of 2.5, and a
maximum of 21 days. The Gaussian process prior for
![\\boldsymbol R](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20R "\boldsymbol R")
is specified by
![\\alpha = 1](https://latex.codecogs.com/png.latex?%5Calpha%20%3D%201 "\alpha = 1")
and
![\\ell = 10](https://latex.codecogs.com/png.latex?%5Cell%20%3D%2010 "\ell = 10").
We fit two models,
![\\mathcal M\_{0.1}](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%7B0.1%7D "\mathcal M_{0.1}")
where
![k = 0.1](https://latex.codecogs.com/png.latex?k%20%3D%200.1 "k = 0.1")
and
![\\mathcal M\_\\infty](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%5Cinfty "\mathcal M_\infty")
where
![k \\to \\infty](https://latex.codecogs.com/png.latex?k%20%5Cto%20%5Cinfty "k \to \infty").

We smooth over day-of-the-week effects in the epidemic curve by taking
![y\_t](https://latex.codecogs.com/png.latex?y_t "y_t") to be the 7-day
moving average of reported cases on day
![t](https://latex.codecogs.com/png.latex?t "t").

``` r
first_day <- dmy(05122020)
last_day <- dmy(31012021)
df <- covid_incidence_roi_epidemiological_date %>% 
  mutate(ma_count = stats::filter(count, rep(1/7, 7)) %>% round) %>% 
  filter(date >= first_day & date <= last_day) 
D <- nrow(df)
# initialise heterogeneous disease reproduction
init_list <- lapply(
  1:4, function(i) {
    initialise_lgp_Rt(
      epidemic_curve = df$ma_count, 
      gp_amplitude = 1, k = 0.1
      )
    }
  )
# fit heterogeneous disease reproduction
M_0.1 <- fit_Rt_lgp(
  epidemic_curve = df$ma_count, seed_days = 5,
  import_rate = rep(1, D), 
  generation_interval_mean = 5, generation_interval_sd = 2.5,
  generation_interval_length = 21,
  gp_amplitude = 1, gp_length_scale = 10,
  k = 0.1,
  cores = 4, refresh = 500, 
  init = init_list
)
# fit homogeneous disease reproduction
M_inf <- fit_Rt_lgp(
  epidemic_curve = df$ma_count, seed_days = 5,
  import_rate = rep(1, D), 
  generation_interval_mean = 5, generation_interval_sd = 2.5,
  generation_interval_length = 21,
  gp_amplitude = 1, gp_length_scale = 10,
  k = Inf,
  cores = 4, refresh = 500,
  control = list(max_treedepth = 15)
)

wt <- wallinga_teunis(
  incid = df$ma_count, 
  method = "parametric_si",
  config = list(
    t_start = (5:(D-2)),
    t_end = 7:D,
    mean_si = 5,
    std_si = 2.5,
    n_sim = 100
  )
)
```

Having fit these models to the epidemic curve we can explore the fitted
posterior for
![\\boldsymbol R](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20R "\boldsymbol R").
Note that heterogeneous disease reproduction results in greater
uncertainty on estimates for
![R\_t](https://latex.codecogs.com/png.latex?R_t "R_t") on each day. We
include Wallinga & Teunis’ estimate for
![\\boldsymbol R](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20R "\boldsymbol R"),
implemented in the `EpiEstem` package, for comparison. All three
estimates are in broad agreement. Estimates for
![R\_t](https://latex.codecogs.com/png.latex?R_t "R_t") on the final ten
days of the period are omitted as these index cases remain infectious
and so their reproduction number depends on secondary infections that
may or may not occur in the future.

``` r
R_0.1 <- rstan::extract(M_0.1, "R") %>% 
  as.data.frame() %>% 
  colMeans()
ci_0.1 <- rstan::extract(M_0.1, "R") %>% 
  as.data.frame() %>% 
  hdi()

R_inf <- rstan::extract(M_inf, "R") %>% 
  as.data.frame() %>% 
  colMeans() 
ci_inf <- rstan::extract(M_inf, "R") %>% 
  as.data.frame() %>% 
  hdi()

y_scalar <- ci_0.1[, -(1:5)] %>% 
  max()

data.frame(
  date = df$date,
  a = R_0.1,
  b = R_inf,
  wt = c(rep(NA, 5), wt$R$`Mean(R)`, rep(NA, 1))
  ) %>%
  filter(date <= last_day -10) %>% 
  filter(date >= first_day + 5) %>% 
  reshape2::melt(id.vars = "date") %>%
  ggplot() +
  theme_classic(base_size = 10) +
  theme(legend.position = c(0.875, 0.75)) +
  geom_bar(
    data = df,
    aes(x = date, y =  y_scalar * ma_count / max(ma_count)), stat = "identity",
    alpha = 0.5
  ) +
  geom_ribbon(
    data = data.frame(
      date = df$date[6:(D-1)], 
      hdi = cbind(
        lower = wt$R$`Quantile.0.025(R)`, 
        upper = wt$R$`Quantile.0.975(R)`), 
      model =  "wt"
      ) %>% 
      rbind(data.frame(date = df$date, hdi = t(ci_0.1), model =  "a")) %>% 
      rbind(data.frame(date = df$date, hdi = t(ci_inf), model =  "b")) %>% 
      filter(date <= last_day - 10) %>% 
      filter(date >= first_day + 5),
    aes(x = date, ymin = hdi.lower, ymax = hdi.upper, fill = model),
    alpha= 0.25
  ) +
  geom_line(aes(x = date, y = value, color = variable)) +
  scale_color_viridis_d(labels = c(bquote("M"[0.1]), bquote("M"[infinity]), "W&T")) +
  scale_fill_viridis_d(labels = c(bquote("M"[0.1]), bquote("M"[infinity]), "W&T")) +
  geom_hline(yintercept = 1, lty = 3) +
  scale_y_continuous(
    TeX("Reproduction Number $(R_t)$"),
    sec.axis = sec_axis(
      ~ . * max(df$ma_count) / y_scalar, 
      name = TeX("7 Day Moving Average Incidence $(Y_t)$")
      )
  ) +
  labs(
    x = "Epidemiological Date",
    color = NULL,
    fill = NULL
  ) 
```

<img src="man/figures/README-time_varying_reproduction-1.png" width="100%" />

## Model comparison

We adopt a Leave-Future-Out (LFO) Cross-Validation (CV) approach to
model comparison, estimating the posterior predictive density for
![y\_{t+1}](https://latex.codecogs.com/png.latex?y_%7Bt%2B1%7D "y_{t+1}")
under
![\\mathcal M\_k](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_k "\mathcal M_k")

![
    p \\left( y\_{t+1} \\mid \\boldsymbol y\_{0:t}, \\mathcal M\_k \\right) = \\int p \\left( y\_{t+1} \\mid \\boldsymbol \\theta\_{0:t}, \\boldsymbol y\_{0:t}, \\mathcal M\_k \\right) p \\left( \\boldsymbol \\theta\_{0:t} \\mid \\boldsymbol y\_{0:t}, \\mathcal M\_k \\right) d \\boldsymbol \\theta\_{0:t},
](https://latex.codecogs.com/png.latex?%0A%20%20%20%20p%20%5Cleft%28%20y_%7Bt%2B1%7D%20%5Cmid%20%5Cboldsymbol%20y_%7B0%3At%7D%2C%20%5Cmathcal%20M_k%20%5Cright%29%20%3D%20%5Cint%20p%20%5Cleft%28%20y_%7Bt%2B1%7D%20%5Cmid%20%5Cboldsymbol%20%5Ctheta_%7B0%3At%7D%2C%20%5Cboldsymbol%20y_%7B0%3At%7D%2C%20%5Cmathcal%20M_k%20%5Cright%29%20p%20%5Cleft%28%20%5Cboldsymbol%20%5Ctheta_%7B0%3At%7D%20%5Cmid%20%5Cboldsymbol%20y_%7B0%3At%7D%2C%20%5Cmathcal%20M_k%20%5Cright%29%20d%20%5Cboldsymbol%20%5Ctheta_%7B0%3At%7D%2C%0A "
    p \left( y_{t+1} \mid \boldsymbol y_{0:t}, \mathcal M_k \right) = \int p \left( y_{t+1} \mid \boldsymbol \theta_{0:t}, \boldsymbol y_{0:t}, \mathcal M_k \right) p \left( \boldsymbol \theta_{0:t} \mid \boldsymbol y_{0:t}, \mathcal M_k \right) d \boldsymbol \theta_{0:t},
")

where
![\\boldsymbol y\_{0:t} = \\left( y\_0, y\_1, \\dots, y\_t \\right)](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20y_%7B0%3At%7D%20%3D%20%5Cleft%28%20y_0%2C%20y_1%2C%20%5Cdots%2C%20y_t%20%5Cright%29 "\boldsymbol y_{0:t} = \left( y_0, y_1, \dots, y_t \right)")
and
![\\boldsymbol \\theta\_{0:t}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Ctheta_%7B0%3At%7D "\boldsymbol \theta_{0:t}")
is the set of model parameters and latent variables up to time
![t](https://latex.codecogs.com/png.latex?t "t"). This amounts to
one-step-ahead prediction (1-SAP) for the epidemic curve.

Given
![\\boldsymbol \\theta\_{0:t}^{\\left( m \\right)}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Ctheta_%7B0%3At%7D%5E%7B%5Cleft%28%20m%20%5Cright%29%7D "\boldsymbol \theta_{0:t}^{\left( m \right)}")
drawn from
![p \\left( \\boldsymbol \\theta\_{0:t} \\mid \\boldsymbol y\_{0:t}, \\mathcal M\_k \\right)](https://latex.codecogs.com/png.latex?p%20%5Cleft%28%20%5Cboldsymbol%20%5Ctheta_%7B0%3At%7D%20%5Cmid%20%5Cboldsymbol%20y_%7B0%3At%7D%2C%20%5Cmathcal%20M_k%20%5Cright%29 "p \left( \boldsymbol \theta_{0:t} \mid \boldsymbol y_{0:t}, \mathcal M_k \right)")
for
![m = 1, \\dots, M](https://latex.codecogs.com/png.latex?m%20%3D%201%2C%20%5Cdots%2C%20M "m = 1, \dots, M")
we have that

![
    p \\left( y\_{t+1} \\mid \\boldsymbol y\_{0:t}, \\mathcal M\_k \\right) \\approx \\sum\_{m = 1}^M p \\left( y\_{t+1} \\mid \\boldsymbol \\theta\_{0:t}^{\\left( m \\right)}, \\boldsymbol y\_{0:t}, \\mathcal M\_k \\right)
](https://latex.codecogs.com/png.latex?%0A%20%20%20%20p%20%5Cleft%28%20y_%7Bt%2B1%7D%20%5Cmid%20%5Cboldsymbol%20y_%7B0%3At%7D%2C%20%5Cmathcal%20M_k%20%5Cright%29%20%5Capprox%20%5Csum_%7Bm%20%3D%201%7D%5EM%20p%20%5Cleft%28%20y_%7Bt%2B1%7D%20%5Cmid%20%5Cboldsymbol%20%5Ctheta_%7B0%3At%7D%5E%7B%5Cleft%28%20m%20%5Cright%29%7D%2C%20%5Cboldsymbol%20y_%7B0%3At%7D%2C%20%5Cmathcal%20M_k%20%5Cright%29%0A "
    p \left( y_{t+1} \mid \boldsymbol y_{0:t}, \mathcal M_k \right) \approx \sum_{m = 1}^M p \left( y_{t+1} \mid \boldsymbol \theta_{0:t}^{\left( m \right)}, \boldsymbol y_{0:t}, \mathcal M_k \right)
")

If ![L + 1](https://latex.codecogs.com/png.latex?L%20%2B%201 "L + 1")
observations are required for reliable 1-SAP, then we assess the
performance of
![\\mathcal M\_k](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_k "\mathcal M_k")
using its LFO expected log predictive density (ELPD)

![
    \\operatorname{ELPD}\_{\\operatorname{LFO}}^k = \\sum\_{t = L}^{N-1} \\log p \\left( y\_{t+1} \\mid \\boldsymbol y\_{0:t}, \\mathcal M\_k \\right).
](https://latex.codecogs.com/png.latex?%0A%20%20%20%20%5Coperatorname%7BELPD%7D_%7B%5Coperatorname%7BLFO%7D%7D%5Ek%20%3D%20%5Csum_%7Bt%20%3D%20L%7D%5E%7BN-1%7D%20%5Clog%20p%20%5Cleft%28%20y_%7Bt%2B1%7D%20%5Cmid%20%5Cboldsymbol%20y_%7B0%3At%7D%2C%20%5Cmathcal%20M_k%20%5Cright%29.%0A "
    \operatorname{ELPD}_{\operatorname{LFO}}^k = \sum_{t = L}^{N-1} \log p \left( y_{t+1} \mid \boldsymbol y_{0:t}, \mathcal M_k \right).
")

We refit
![\\mathcal M\_k](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_k "\mathcal M_k")
for each estimate of
![p \\left( y\_{t+1} \\mid \\boldsymbol y\_{0:t}, \\mathcal M\_k \\right)](https://latex.codecogs.com/png.latex?p%20%5Cleft%28%20y_%7Bt%2B1%7D%20%5Cmid%20%5Cboldsymbol%20y_%7B0%3At%7D%2C%20%5Cmathcal%20M_k%20%5Cright%29 "p \left( y_{t+1} \mid \boldsymbol y_{0:t}, \mathcal M_k \right)")
assuming that

![
    p \\left( y\_{t+1} \\mid \\boldsymbol y\_{0:t}, \\mathcal M\_k \\right) \\approx p \\left( y\_{t+1} \\mid \\boldsymbol y\_{(t - L):t}, \\mathcal M\_k \\right),
](https://latex.codecogs.com/png.latex?%0A%20%20%20%20p%20%5Cleft%28%20y_%7Bt%2B1%7D%20%5Cmid%20%5Cboldsymbol%20y_%7B0%3At%7D%2C%20%5Cmathcal%20M_k%20%5Cright%29%20%5Capprox%20p%20%5Cleft%28%20y_%7Bt%2B1%7D%20%5Cmid%20%5Cboldsymbol%20y_%7B%28t%20-%20L%29%3At%7D%2C%20%5Cmathcal%20M_k%20%5Cright%29%2C%0A "
    p \left( y_{t+1} \mid \boldsymbol y_{0:t}, \mathcal M_k \right) \approx p \left( y_{t+1} \mid \boldsymbol y_{(t - L):t}, \mathcal M_k \right),
")

such that we only fit to
![L + 1](https://latex.codecogs.com/png.latex?L%20%2B%201 "L + 1")
observations up to time ![t](https://latex.codecogs.com/png.latex?t "t")
to obtain the estimate (this becomes
![L + N\_0](https://latex.codecogs.com/png.latex?L%20%2B%20N_0 "L + N_0")
observations when
![N\_0](https://latex.codecogs.com/png.latex?N_0 "N_0") days seed the
epidemic).

For illustrative purposes, we present LFO-CV for the 15 days from 10-24
December, 2020 with
![L = 21](https://latex.codecogs.com/png.latex?L%20%3D%2021 "L = 21")
and
![N\_0 = 5](https://latex.codecogs.com/png.latex?N_0%20%3D%205 "N_0 = 5").
This comparison requires a few minutes to complete.

``` r
n_samples <- 4000
N0 <- 5
L <- 21 # lead in days to prediction
first_validation_day <- dmy(10122020)
last_validation_day <- dmy(25122020)
M <- as.numeric(last_validation_day - first_validation_day)
# Get 7-day moving average
df <- covid_incidence_roi_epidemiological_date %>%
  mutate(ma_count = stats::filter(count, rep(1/7, 7)) %>% round)
# candidate models
candidate_k <- c(0.1, Inf)

posterior_predictive <- log_lik <- list()
for (j in seq_along(candidate_k)) {
  i <- 1
  day <- first_validation_day
  tmp_posterior_predictive <- tmp_log_lik <- array(dim = c(n_samples, M))
  # pb <- txtProgressBar(min = 0, max = M, style = 3)
  while (day < last_validation_day) {
    df_i <- df %>%
      filter(date >= day - (L + N0) & date < day)
    fit_i <- fit_Rt_lgp(
      epidemic_curve = df_i$ma_count, seed_days = 5,
      import_rate = rep(1, (L + N0)),
      generation_interval_mean = 5,
      generation_interval_sd = 2.5,
      generation_interval_length = 21,
      gp_amplitude = 1, gp_length_scale = 10,
      k = candidate_k[j],
      ahead = TRUE,
      next_day_cases = df$ma_count[df$date == day],
      next_day_import_rate = 1, refresh = 0
    )
    tmp_posterior_predictive[, i] <- unlist(rstan::extract(fit_i, "y_rep_ahead"))
    tmp_log_lik[, i] <- unlist(rstan::extract(fit_i, "log_lik_ahead"))
    day <- day + 1
    # setTxtProgressBar(pb, i)
    i <- i + 1
  }
  posterior_predictive[[j]] <- tmp_posterior_predictive
  log_lik[[j]] <- tmp_log_lik
  # close(pb)
}
```

``` r
pointwise_elpd <- sapply(
  log_lik, 
  function(x) colLogSumExps(x) - log(n_samples)
  )
elpd_lfo <- apply(pointwise_elpd, 2, sum)
se_elpd_lfo <- apply(
  pointwise_elpd, 2, 
  function(x) sqrt(length(x) * var(x))
  )

best_fit <- which.max(elpd_lfo)
diff_elpd <- sweep(pointwise_elpd, 1, pointwise_elpd[, best_fit]) %>% 
  apply(2, sum)
se_diff_elpd <- sweep(pointwise_elpd, 1, pointwise_elpd[, best_fit]) %>% 
  apply(2, function(x) sqrt(length(x) * var(x)))
```

|        |   elpd | se\_elpd | elpd\_diff | se\_elpd\_diff |
|:-------|-------:|---------:|-----------:|---------------:|
| model1 | -80.40 |     3.06 |        0.0 |           0.00 |
| model2 | -82.19 |     3.54 |       -1.8 |           2.29 |

LFO model comparison of
![\\mathcal M\_{0.1}](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%7B0.1%7D "\mathcal M_{0.1}")
(model1) and
![\\mathcal M\_{\\infty}](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%7B%5Cinfty%7D "\mathcal M_{\infty}")
(model2).

Model comparison for the validation period provides some support for
![\\mathcal M\_{0.1}](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%7B0.1%7D "\mathcal M_{0.1}")
over
![\\mathcal M\_\\infty](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%5Cinfty "\mathcal M_\infty"),
indicating that heterogeneity is a feature of the COVID-19 epidemic in
the Republic of Ireland. This suggests that estimates for
![\\boldsymbol R](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20R "\boldsymbol R")
by
![\\mathcal M\_{0.1}](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%7B0.1%7D "\mathcal M_{0.1}")
provide more appropriate uncertainty quantification and should be
preferred to those of
![\\mathcal M\_{\\infty}](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%7B%5Cinfty%7D "\mathcal M_{\infty}").

## Alternative model comparison

Pareteo-smoothed importance sampling (PSIS) for approximate
leave-one-out cross-validation (LOO-CV) provides an easy to implement
measure of model fit. Although time series data are not exchangeable and
violate the assumptions underpinning LOO-CV, PSIS LOO-CV offers an
efficient approach to model comparison. Applying this technique to
![\\mathcal M\_{0.1}](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%7B0.1%7D "\mathcal M_{0.1}")
and
![\\mathcal M\_\\infty](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%5Cinfty "\mathcal M_\infty")
fit to COVID-19 cases from December 10, 2020, to January 31 2020, we
find that all Pareto-![k](https://latex.codecogs.com/png.latex?k "k")
diagnostic values are less than 0.7, indicating that both models fit the
data well.

``` r
loo_0.1 <- loo(M_0.1, moment_match = TRUE)
#> Warning: Some Pareto k diagnostic values are slightly high. See help('pareto-k-diagnostic') for details.
loo_inf <- loo(M_inf, moment_match = TRUE)
#> Warning: Some Pareto k diagnostic values are slightly high. See help('pareto-k-diagnostic') for details.
```

|        | elpd\_diff | se\_diff | elpd\_loo | se\_elpd\_loo | p\_loo | se\_p\_loo |  looic | se\_looic |
|:-------|-----------:|---------:|----------:|--------------:|-------:|-----------:|-------:|----------:|
| model1 |       0.00 |      0.0 |   -265.93 |          3.61 |   6.27 |       1.14 | 531.86 |      7.21 |
| model2 |     -16.07 |      5.8 |   -282.00 |          8.42 |   8.69 |       2.04 | 564.01 |     16.85 |

PSIS-LOO model comparison favours
![\\mathcal M\_{0.1}](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%7B0.1%7D "\mathcal M_{0.1}")
(model1) over
![\\mathcal M\_{\\infty}](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%7B%5Cinfty%7D "\mathcal M_{\infty}")
(model2).

As before, this model comparison supports
![\\mathcal M\_{0.1}](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%7B0.1%7D "\mathcal M_{0.1}")
over
![\\mathcal M\_\\infty](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%5Cinfty "\mathcal M_\infty"),
although it will over-estimate the 1-SAP accuracy.

## Alternative priors

The package also includes a simpler prior for
![\\boldsymbol R](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20R "\boldsymbol R")
which fits data more efficiently. The histogram estimator for
![\\boldsymbol R](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20R "\boldsymbol R")
is defined by

![
\\begin{aligned}
    R\_t &= \\sum\_{k = 1}^B r\_k \\mathbf 1 \\left\\{ t \\in I\_k \\right\\}, \\\\
    \\log r\_k  &\\sim \\mathcal N \\left( 0, 1\\right),
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0A%20%20%20%20R_t%20%26%3D%20%5Csum_%7Bk%20%3D%201%7D%5EB%20r_k%20%5Cmathbf%201%20%5Cleft%5C%7B%20t%20%5Cin%20I_k%20%5Cright%5C%7D%2C%20%5C%5C%0A%20%20%20%20%5Clog%20r_k%20%20%26%5Csim%20%5Cmathcal%20N%20%5Cleft%28%200%2C%201%5Cright%29%2C%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
    R_t &= \sum_{k = 1}^B r_k \mathbf 1 \left\{ t \in I_k \right\}, \\
    \log r_k  &\sim \mathcal N \left( 0, 1\right),
\end{aligned}
")

for the set of disjoint intervals
![I\_1, \\dots, I\_B](https://latex.codecogs.com/png.latex?I_1%2C%20%5Cdots%2C%20I_B "I_1, \dots, I_B")
discretising time. The interval
![I\_j](https://latex.codecogs.com/png.latex?I_j "I_j") defines the
![j^{th}](https://latex.codecogs.com/png.latex?j%5E%7Bth%7D "j^{th}")
bin for the estimator where

![
I\_j = \\left( N - j \\delta, N - (j - 1) \\delta \\right\],
](https://latex.codecogs.com/png.latex?%0AI_j%20%3D%20%5Cleft%28%20N%20-%20j%20%5Cdelta%2C%20N%20-%20%28j%20-%201%29%20%5Cdelta%20%5Cright%5D%2C%0A "
I_j = \left( N - j \delta, N - (j - 1) \delta \right],
")

covers a bin width of
![\\delta](https://latex.codecogs.com/png.latex?%5Cdelta "\delta") days,
then the histogram estimator for
![\\boldsymbol R](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20R "\boldsymbol R")
is defined by

![
B = \\left\\lceil \\frac{N+1}{\\delta} \\right\\rceil.
](https://latex.codecogs.com/png.latex?%0AB%20%3D%20%5Cleft%5Clceil%20%5Cfrac%7BN%2B1%7D%7B%5Cdelta%7D%20%5Cright%5Crceil.%0A "
B = \left\lceil \frac{N+1}{\delta} \right\rceil.
")

Assuming that
![\\delta = 7](https://latex.codecogs.com/png.latex?%5Cdelta%20%3D%207 "\delta = 7")
and
![k \\to \\infty](https://latex.codecogs.com/png.latex?k%20%5Cto%20%5Cinfty "k \to \infty"),
this model fits the epidemic curve from April 01, 2020 to February 14,
2021 in seconds.

``` r
first_day <- dmy(01042020)
last_day <- dmy(14022021)
df <- covid_incidence_roi_epidemiological_date %>% 
  mutate(ma_count = stats::filter(count, rep(1/7, 7)) %>% round) %>% 
  filter(date >= first_day & date <= last_day) 
D <- nrow(df)
# fit homogeneous disease reproduction
H_inf <- fit_Rt_hist(
  epidemic_curve = df$ma_count, seed_days = 5,
  import_rate = rep(1, D), 
  generation_interval_mean = 5, generation_interval_sd = 2.5,
  generation_interval_length = 21,
  bin_width = 7,
  k = Inf,
  cores = 4, refresh = 500
)
```

<img src="man/figures/README-histogram_estimator-1.png" width="100%" />

The LFO framework for model comparison can also be applied to this model
to objectively tune the bin width
![\\delta](https://latex.codecogs.com/png.latex?%5Cdelta "\delta"). The
model also accommodates heterogeneous disease reproduction given some
fixed value for ![k](https://latex.codecogs.com/png.latex?k "k").
