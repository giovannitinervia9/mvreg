% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loglik_derivatives.R
\name{ll}
\alias{ll}
\title{Loglikelihood of heteroskedastic linear model}
\usage{
ll(y, x, z, b, t)
}
\arguments{
\item{y}{Vector of response variable.}

\item{x}{Matrix of explanatory variables for mean component.}

\item{z}{Matrix of explanatory variables for variance component.}

\item{b}{Vector of parameters for mean component.}

\item{t}{Vector of parameters for variance component.}
}
\value{
Value of loglikelihood for a given sample and parameters' values.
}
\description{
Loglikelihood of heteroskedastic linear model
}
\examples{
mvreg_mod <- mvreg(Sepal.Length ~ Species, data = iris) # same formula for mean and variance
y <- mvreg_mod$y # response variable
x <- mvreg_mod$x # model.matrix for mean component
z <- mvreg_mod$z # model.matrix for variance component
b <- coef(mvreg_mod, "mu") # coefficients of mean component
t <- coef(mvreg_mod, "s2") # coefficients of variance component

# loglikelihood
ll(y, x, z, b, t)

}
