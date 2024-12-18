% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvreg_simul.R
\name{print.simul_mvreg}
\alias{print.simul_mvreg}
\title{Print method for simul_mvreg}
\usage{
\method{print}{simul_mvreg}(x, digits = max(3L, getOption("digits") - 3L), ...)
}
\arguments{
\item{x}{A simul_mvreg object, containing the results of the simulation study.}

\item{digits}{The minimum number of significant digits to be used. Default is
set to max(3L, getOption("digits") - 3L).}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
Prints results of simul_mvreg.
}
\description{
This method prints the results of the simulation study for a \code{mvreg} model.
It includes the number of observations, the number of simulations, the number
of converged simulations, and the total time taken for the simulations.
}
\examples{
n <- 100
x <- cbind(1, rnorm(n))
z <- cbind(1, x[, 2], x[, 2]^2)
b <- rnorm(ncol(x))
t <- rnorm(ncol(z))
result <- mvreg_simul(x, z, b, t, nsim = 100, seed = 43)
print(result)
}
