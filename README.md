
<!-- README.md is generated from README.Rmd. Please edit that file -->

# assessEpidemicCurves

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/jpmeagher/assessEpidemicCurves.svg?branch=master)](https://travis-ci.com/jpmeagher/assessEpidemicCurves)
<!-- badges: end -->

The goal of assessEpidemicCurves is to model epidemic curves under
heterogeneous disease reproduction, providing estimates for the
time-varying reproduction number and assessing epidemic curves for
evidence of superspreading. This package accompanies \[Reference\]

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
library(ggplot2)
```

We consider the COVID-19 epidemic in the Republic of Ireland from March
1, 2020 to 28 February, 2021.

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

## The model for epidemic curves

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
infections, that is, superspreading is a feature of the epidemic. Note
that
![\\lim\_{k \\to \\infty} \\operatorname{Var} \\left( \\nu\_{t, i}\\right) = 0](https://latex.codecogs.com/png.latex?%5Clim_%7Bk%20%5Cto%20%5Cinfty%7D%20%5Coperatorname%7BVar%7D%20%5Cleft%28%20%5Cnu_%7Bt%2C%20i%7D%5Cright%29%20%3D%200 "\lim_{k \to \infty} \operatorname{Var} \left( \nu_{t, i}\right) = 0")
and so this model admids homogeneous disease reproduction as a special
case.
