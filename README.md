
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
```

We consider the COVID-19 epidemic in the Republic of Ireland from March
1, 2020 to 28 February, 2021.

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

## The generative model for epidemic curves

Analyses are based on the hierarchical model for the daily incidence
count of COVID-19

![
\\begin{aligned}
p \\left( y\_t \\mid \\mu\_t, \\boldsymbol \\eta\_t, \\boldsymbol \\omega \\right) &= \\operatorname{Pois} \\left( y\_t \\mid \\mu\_t + \\sum\_{s = 1}^t \\omega\_s \\eta\_{t-s} \\right), \\\\
p \\left( \\eta\_t \\mid k, R\_t \\right) &= \\operatorname{Gamma} \\left( \\eta\_t \\mid y\_t k, \\frac{k}{R\_t}\\right),
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0Ap%20%5Cleft%28%20y_t%20%5Cmid%20%5Cmu_t%2C%20%5Cboldsymbol%20%5Ceta_t%2C%20%5Cboldsymbol%20%5Comega%20%5Cright%29%20%26%3D%20%5Coperatorname%7BPois%7D%20%5Cleft%28%20y_t%20%5Cmid%20%5Cmu_t%20%2B%20%5Csum_%7Bs%20%3D%201%7D%5Et%20%5Comega_s%20%5Ceta_%7Bt-s%7D%20%5Cright%29%2C%20%5C%5C%0Ap%20%5Cleft%28%20%5Ceta_t%20%5Cmid%20k%2C%20R_t%20%5Cright%29%20%26%3D%20%5Coperatorname%7BGamma%7D%20%5Cleft%28%20%5Ceta_t%20%5Cmid%20y_t%20k%2C%20%5Cfrac%7Bk%7D%7BR_t%7D%5Cright%29%2C%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
p \left( y_t \mid \mu_t, \boldsymbol \eta_t, \boldsymbol \omega \right) &= \operatorname{Pois} \left( y_t \mid \mu_t + \sum_{s = 1}^t \omega_s \eta_{t-s} \right), \\
p \left( \eta_t \mid k, R_t \right) &= \operatorname{Gamma} \left( \eta_t \mid y_t k, \frac{k}{R_t}\right),
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
p \\left( z\_{t, i} \\mid \\nu\_{t\_i} \\right) &= \\operatorname{Pois} \\left( z\_{t, i} \\mid \\nu\_{t, i} \\right), \\\\
p \\left( \\nu\_{t, i} \\mid k, R\_t \\right) &= \\operatorname{Gamma} \\left( k, \\frac{k}{R\_t}\\right),
\\end{aligned}
](https://latex.codecogs.com/png.latex?%0A%5Cbegin%7Baligned%7D%0Ap%20%5Cleft%28%20z_%7Bt%2C%20i%7D%20%5Cmid%20%5Cnu_%7Bt_i%7D%20%5Cright%29%20%26%3D%20%5Coperatorname%7BPois%7D%20%5Cleft%28%20z_%7Bt%2C%20i%7D%20%5Cmid%20%5Cnu_%7Bt%2C%20i%7D%20%5Cright%29%2C%20%5C%5C%0Ap%20%5Cleft%28%20%5Cnu_%7Bt%2C%20i%7D%20%5Cmid%20k%2C%20R_t%20%5Cright%29%20%26%3D%20%5Coperatorname%7BGamma%7D%20%5Cleft%28%20k%2C%20%5Cfrac%7Bk%7D%7BR_t%7D%5Cright%29%2C%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
p \left( z_{t, i} \mid \nu_{t_i} \right) &= \operatorname{Pois} \left( z_{t, i} \mid \nu_{t, i} \right), \\
p \left( \nu_{t, i} \mid k, R_t \right) &= \operatorname{Gamma} \left( k, \frac{k}{R_t}\right),
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
log-Gaussian process prior for
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
and so we only analyse a subset of the data here. We examine reported
COVID-19 cases from December 10, 2020, to January 31 2020, allowing the
5 days from December 5 to December 9 inclusive to seed the epidemic.

We assume that
![\\mu\_t = 1](https://latex.codecogs.com/png.latex?%5Cmu_t%20%3D%201 "\mu_t = 1")
for all ![t](https://latex.codecogs.com/png.latex?t "t") and that
![\\boldsymbol \\omega](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Comega "\boldsymbol \omega")
is a discretised gamma distribution such that generation intervals have
a mean of 5 days, standard deviation of 2.5, and a maximum of 21 days.
The log-Gaussian process prior for
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

We adjust the epidemic curve for day of the week effects by taking
![y\_t](https://latex.codecogs.com/png.latex?y_t "y_t") to be the 7-day
moving average of reported cases on day
![t](https://latex.codecogs.com/png.latex?t "t").

``` r
df <- covid_incidence_roi_epidemiological_date %>% 
  mutate(ma_count = stats::filter(count, rep(1/7, 7)) %>% round) %>% 
  filter(date >= dmy(05122020) & date <= dmy(31012021)) 

D <- nrow(df)
# initialise heterogeneous disease reproduction
init_list <- lapply(
  1:4, function(i) {
    initialise_lgp_Rt(
      epidemic_curve = df$ma_count, 
      gp_amplitude = 1, k = 0.1
)
  } )
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
    t_start = (6:D),
    t_end = 6:D,
    mean_si = 5,
    std_si = 2.5,
    n_sim = 3
  )
)
```

Having fit the models to the epidemic curve we can explore the fitted
posterior for
![\\boldsymbol R](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20R "\boldsymbol R").
Note that heterogeneous disease reproduction results in greater
uncertainty on estimates for
![R\_t](https://latex.codecogs.com/png.latex?R_t "R_t") on each day. We
include Wallinga & Teunis’ estimate for
![\\boldsymbol R](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20R "\boldsymbol R"),
implemented in the `EpiEstem` package, for comparison. All three
estimates are in broad agreement up until the final days of the assessed
period. The disagreement at this point is consequence of the differing
estimation procedures.

``` r
R_0.1 <- rstan::extract(M_0.1, "R") %>% 
  as.data.frame() %>% 
  magrittr::extract(sample.int(4000, 100), ) %>%
  t() %>% 
  unname() %>% 
  data.frame(date = df$date, R_0.1 = .) %>% 
  reshape2::melt(id = "date")

R_inf <- rstan::extract(M_inf, "R") %>% 
  as.data.frame() %>% 
  magrittr::extract(sample.int(4000, 100), ) %>%
  t() %>% 
  unname() %>% 
  data.frame(date = df$date, R_inf = .) %>% 
  reshape2::melt(id = "date") 



R <- rbind(
  cbind(R_0.1, model = "0.1"),
  cbind(R_inf, model = "inf")
) 

last_day <-dmy(20012021)
y_scalar <- 2.5
R %>%
  dplyr::filter(date <= last_day) %>% 
  ggplot() +
  geom_bar(
    data = df,
    aes(x = date, y = y_scalar * ma_count /  max(df$ma_count) ), stat = "identity",
    alpha = 0.5
  ) +
  geom_line(aes(x = date, y = value, group = variable, color = model), alpha = 0.25) +
  geom_hline(yintercept = 1, lty = 3) +
  scale_color_viridis_d(labels = c("0.1", bquote(infinity))) +
  scale_y_continuous(
    bquote("R"["t"]),
    sec.axis = sec_axis(~ . * max(df$ma_count) / y_scalar, name = "7 Day Moving Average Incidence")
  ) +
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  geom_line(
    data = data.frame(
      date = df$date[6:D], R = wt$R$`Mean(R)`
    ) %>% 
      dplyr::filter(date <= last_day),
    aes(x = date, y = R), lwd = 1
  ) +
  labs(
    title = "Fitted reproduction numbers",
    x = "Epidemiological Date",
    color = "k"
  )
#> Warning: Use of `df$ma_count` is discouraged. Use `ma_count` instead.
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

## Model comparison

PSIS-LOO provides an easy to implement measure of model fit. All
Pareto-![k](https://latex.codecogs.com/png.latex?k "k") diagnostic
values are less than 0.7, indicating that both models fit the data well
and estimates for the expected log point-wise predictive density (elpd)
for a new dataset are reliable.

``` r
loo_0.1 <- loo(M_0.1, moment_match = TRUE)
#> Warning: Some Pareto k diagnostic values are slightly high. See help('pareto-k-diagnostic') for details.
loo_inf <- loo(M_inf, moment_match = TRUE)
#> Warning: Some Pareto k diagnostic values are slightly high. See help('pareto-k-diagnostic') for details.

loo_compare(loo_0.1, loo_inf) %>% 
  kable(digits = 2, caption = "PSIS-LOO model selection favours $\\mathcal M_{0.1}$ (model1) over $\\mathcal M_{\\infty}$ (model2).")
```

|        | elpd\_diff | se\_diff | elpd\_loo | se\_elpd\_loo | p\_loo | se\_p\_loo |  looic | se\_looic |
|:-------|-----------:|---------:|----------:|--------------:|-------:|-----------:|-------:|----------:|
| model1 |       0.00 |     0.00 |   -265.86 |          3.60 |   6.22 |       1.11 | 531.72 |      7.20 |
| model2 |     -16.31 |     5.95 |   -282.17 |          8.52 |   8.36 |       1.90 | 564.34 |     17.04 |

PSIS-LOO model selection favours
![\\mathcal M\_{0.1}](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%7B0.1%7D "\mathcal M_{0.1}")
(model1) over
![\\mathcal M\_{\\infty}](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%7B%5Cinfty%7D "\mathcal M_{\infty}")
(model2).

This model comparison supports
![\\mathcal M\_{0.1}](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%7B0.1%7D "\mathcal M_{0.1}")
over
![\\mathcal M\_\\infty](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%5Cinfty "\mathcal M_\infty"),
indicating that superspreading is a feature of the COVID-19 epidemic in
the Republic of Ireland. This suggests that estimates for
![\\boldsymbol R](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20R "\boldsymbol R")
offered by
![\\mathcal M\_{0.1}](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%7B0.1%7D "\mathcal M_{0.1}")
provide a more appropriate quantification of uncertainty and should be
preferred to those of
![\\mathcal M\_{\\infty}](https://latex.codecogs.com/png.latex?%5Cmathcal%20M_%7B%5Cinfty%7D "\mathcal M_{\infty}").
