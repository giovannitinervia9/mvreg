
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mvreg <img src="man/figures/logo_01.png" align="right" width="118"/>

<!-- badges: start -->
<!-- badges: end -->

The `mvreg` package provides tools for estimating linear regression
models where both the mean and variance components are modeled as linear
functions of predictors. This package extends traditional linear
regression by allowing the variance of the response variable to depend
on a set of covariates, making it useful for scenarios where
heteroscedasticity is present.

## Installation

You can install the development version of mvreg from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("giovannitinervia9/mvreg")
```

## Model specification

The heteroscedastic linear model is specified as
$$Y_i|\mathbf{x}_i, \mathbf{z}_i \sim \mathcal{N}(\mu_i = \mathbf{x}_i'\boldsymbol{\beta}, \sigma^2_i = \exp\{\mathbf{z}_i'\boldsymbol{\tau}\})$$
where $\mathbf{x}_i$ and $\mathbf{z}_i$ are both row vectors of
covariates, the former for the mean component and the latter for the
variance component, and they need not be equal. $\boldsymbol{\beta}$ and
$\boldsymbol{\tau}$ are the vectors of parameters for the mean and
variance components, respectively.

Variance component is modelled via a logarithmic link function

$$\log(\sigma^2_i) = \mathbf{z}_i'\boldsymbol{\tau}$$

## How to fit a heteroscedastic linear model

To fit a heteroscedastic linear model, you can use `mvreg()` function,
in which you can specify a `formula.mu` for the mean component and a
`formula.s2` for the variance component.

Take, for example, the dataset `iris` and let’s fit a model in which
`Sepal.Length` depends on `Species` for both the mean and variance
components

``` r
mod <- mvreg(
  formula.mu = Sepal.Length ~ Species,
  formula.s2 = Sepal.Length ~ Species,
  data = iris
)
summary(mod)
#> 
#> Call:
#> mvreg(formula.mu = Sepal.Length ~ Species, formula.s2 = Sepal.Length ~ 
#>     Species, data = iris)
#> 
#> 
#> Residuals:
#>    Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#> -1.6880   -0.3285   -0.0060    0.0000    0.3120    1.3120
#> 
#> 
#>    logLik       AIC       BIC 
#> -103.4857  218.9715  237.0353 
#> 
#> 
#> Coefficients for mean component:
#>                      estimate      se z.value Pr(>|z|)    
#> mu.const              5.00600 0.04935  101.44   <2e-16 ***
#> mu.Speciesversicolor  0.93000 0.08751   10.63   <2e-16 ***
#> mu.Speciesvirginica   1.58200 0.10179   15.54   <2e-16 ***
#> 
#> 
#> Coefficients for log(variance) component:
#>                      estimate      se z.value Pr(>|z|)    
#> s2.const              -2.1057  0.2000 -10.528  < 2e-16 ***
#> s2.Speciesversicolor   0.7628  0.2828   2.697    0.007 ** 
#> s2.Speciesvirginica    1.1800  0.2828   4.172 3.02e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> 
#> LRT for comparison with a lm model
#>       -2logLik df   LRT Pr(>Chi)    
#> lm       223.5                      
#> mvreg    207.0  2 16.48 0.000264 ***
```
