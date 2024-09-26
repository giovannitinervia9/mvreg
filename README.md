
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mvreg: Mean-Variance Regression

<!-- badges: start -->
<!-- badges: end -->

The mvreg package provides tools for estimating linear regression models
where both the mean and variance components are modeled as linear
functions of predictors. This package extends traditional linear
regression by allowing the variance of the response variable to depend
on a set of covariates, making it useful for scenarios where
heteroscedasticity (non-constant variance) is present.

## Installation

You can install the development version of mvreg from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("giovannitinervia9/mvreg")
```

## Model specification

The heteroskedastic linear model is specified as
$$Y_i|\mathbf{x}_i, \mathbf{z}_i \sim \mathcal{N}\left(\mu_i = \mathbf{x}_i'\boldsymbol{\beta}, \sigma^2_i = \exp\left\{\mathbf{z}_i'\boldsymbol{\tau}\right\}\right)$$.
$\mathbf{x}_i$ and $\mathbf{z}_i$ are both row vectors of covariates,
the former for the mean component and the latter for the variance
component, and they need not be equal. $\boldsymbol{\beta}$ and
$\boldsymbol{\tau}$ are the vectors of parameters for the mean and
variance components, respectively.

Variance component is modelled via a logarithmic link function

$$\log(\sigma^2_i) = \mathbf{z}_i'\boldsymbol{\tau}$$

$\hat{\boldsymbol{\beta}}$ and $\hat{\boldsymbol{\tau}}$ are estimated
via maximum likelihood estimation. Newton-Raphson method is implemented
to estimate the parameters.

Initial guess $\hat{\boldsymbol{\beta}}^0$ is found by regressing $Y$ on
covariates $\mathbf{X}$. The squared log of residuals
$r_i = \log((y_i - \mathbf{x}_i'\hat{\beta})^2)$ are then regressed on
covariates $\mathbf{Z}$ to find a first guess of
$\hat{\boldsymbol{\tau}}^0$.
