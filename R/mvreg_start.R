#-------------------------------------------------------------------------------


#' Starting values for mvreg
#'
#' This function computes starting values for the parameters of a heteroskedastic linear model.
#' The initial estimates for the mean component parameters are calculated using ordinary least squares (OLS),
#' while the initial estimates for the variance component parameters can be derived from three different methods specified by `start.s2` argument.
#'
#' @param y A numeric vector representing the response variable.
#' @param x A numeric matrix representing the design matrix for the mean component of the model.
#' @param z A numeric matrix representing the design matrix for the variance component of the model.
#' @param start.s2 A character string indicating the method for selecting initial values for the variance component parameters.
#' Options are:
#' - "residuals": uses the logarithm of the squared residuals from the mean component model to estimate starting values for the variance component (default choice).
#' - "gamma": uses a Gamma regression on the squared residuals to estimate starting values for the variance component.
#' - "zero": sets all starting values for the variance component to zero, with the first parameter set to the logarithm of the sample variance of the response variable.
#'
#' @return A list containing:
#' - `start`: A numeric vector of initial parameter estimates, including both mean and variance component parameters.
#' - `b0`: A numeric vector of initial estimates for the mean component parameters.
#' - `t0`: A numeric vector of initial estimates for the variance component parameters.
#'
#' @details
#' A first guess for \eqn{\hat{\boldsymbol{\beta}}_0} is obtained using the standard OLS estimator
#' \eqn{\hat{\boldsymbol{\beta}}_0 = (\mathbf{X}'\mathbf{X})^{-1}\mathbf{X}'\mathbf{y}}.
#'
#' A first guess for \eqn{\hat{\boldsymbol{\tau}}_0} is obtained in different ways depending on the specification of the `start.s2` argument:
#' - If `start.s2` is set to "residuals", the first guess for \eqn{\hat{\boldsymbol{\tau}}_0} is obtained by regressing the logarithm of the squared empirical residuals \eqn{\hat{\varepsilon}_i = y_i - \mathbf{x}_i'\hat{\boldsymbol{\beta}}_0} on the covariates \eqn{\mathbf{Z}}.
#' - If `start.s2` is set to "gamma", a Gamma regression model is fit with the squared empirical residuals as the response variable and \eqn{\mathbf{Z}} as the covariates.
#' - If `start.s2` is set to "zero", the first guess is simply \eqn{\hat{\boldsymbol{\tau}}_0 = (\log\hat{S}^2_y, \boldsymbol{0}_{p-1})'} where \eqn{p} is the number of columns of \eqn{\mathbf{Z}}.
#'
#' @export
#'
#' @importFrom stats Gamma glm.fit
#' @examples
#' n <- 100
#' x1 <- rnorm(n)
#' x2 <- factor(sample(letters[1:3], n, TRUE))
#' x <- model.matrix(~ x1 + x2)
#' z1 <- factor(sample(letters[1:3], n, TRUE))
#' z <- model.matrix(~ z1)
#'
#' b <- rnorm(ncol(x))
#' t <- rnorm(ncol(z))
#'
#' y <- rnorm(n, mean = x %*% b, sd = sqrt(exp(z %*% t)))
#'
#' mvreg_start(y, x, z, start.s2 = "residuals")
#' mvreg_start(y, x, z, start.s2 = "gamma")
#' mvreg_start(y, x, z, start.s2 = "zero")
mvreg_start <- function(y, x, z, start.s2 = c("residuals", "gamma", "zero")) {
  start.s2 <- match.arg(start.s2)
  k <- ncol(x)
  p <- ncol(z)

  b0 <- as.vector(solve(crossprod(x)) %*% crossprod(x, y))

  if (start.s2 == "residuals") {
    r <- log((y - x %*% b0)^2)
    t0 <- as.vector(solve(crossprod(z)) %*% crossprod(z, r))
  } else if (start.s2 == "gamma") {
    r <- (y - x %*% b0)^2
    t0 <- glm.fit(z, r, family = Gamma(link = "log"))$coefficients
  } else if (start.s2 == "zero") {
    t0 <- rep(0, p)
    t0[1] <- log(var(y))
  }

  start <- c(b0, t0)
  names(start) <- c(colnames(x), colnames(z))
  b0 <- start[1:k]
  t0 <- start[(k + 1):length(start)]
  list(start = start, b0 = b0, t0 = t0)
}
