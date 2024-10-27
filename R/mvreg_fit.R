#-------------------------------------------------------------------------------


#' Fitter function for mvreg
#'
#' This function estimates the parameters of a heteroskedastic linear model using
#' either weighted least squares (WLS) or a full Newton-Raphson method.
#'
#' @param y A numeric vector representing the response variable.
#' @param x A numeric matrix representing the design matrix for the mean component of the model.
#' @param z A numeric matrix representing the design matrix for the variance component of the model.
#' @param b0 A numeric vector of starting values for the parameters of the mean component.
#' @param t0 A numeric vector of starting values for the parameters of the variance component.
#' @param tol A positive numeric value indicating the minimum difference between parameter estimates
#'            between two iterations required to stop the algorithm. Default is `1e-10`.
#' @param maxit An integer value indicating the maximum number of iterations allowed. Default is `100`.
#' @param method A character string specifying the method for estimation of parameters of the mean component. Options include:
#' - `"wls"`: uses weighted least squares for the mean parameters and Newton-Raphson for variance parameters.
#' - `"full_nr"`: uses the full Newton-Raphson method.
#'
#' @param vcov.type A character string specifying whether to use the observed or expected Fisher information matrix
#'                  to compute the variance-covariance matrix of the estimates. Options include:
#' - `"expected"`: uses the expected Fisher information.
#' - `"observed"`: uses the observed Fisher information.

#'
#' @return A list containing:
#' - `theta`: A numeric vector of estimated parameters, including both mean and variance component parameters.
#' - `b`: A numeric vector of estimated parameters for the mean component.
#' - `t`: A numeric vector of estimated parameters for the variance component.
#' - `vtheta`: A numeric matrix representing the variance-covariance matrix of the parameter estimates.
#' - `it`: An integer indicating the number of iterations performed.
#'
#' @details
#' The fitting process iteratively updates the estimates for the mean and variance parameters
#' until convergence is reached, as determined by the specified tolerance (`tol`) or until
#' the maximum number of iterations (`maxit`) is reached.
#'
#' If the `method` is set to `"wls"`, the fitting process updates the variance parameters
#' using a Newton-Raphson step based on the current estimates \eqn{\hat{\boldsymbol{\tau}}_{k-1}}.
#' The weights \eqn{w_{i}} are then calculated as \eqn{w_{i} = \dfrac{1}{\exp\left(\mathbf{z}_i'\hat{\boldsymbol{\tau}}_{k}\right)}},
#' which are subsequently used to compute the updated mean parameter estimates:
#' \deqn{\hat{\boldsymbol{\beta}}_k = (\mathbf{X}'\mathbf{W}\mathbf{X})^{-1}\mathbf{X}'\mathbf{W}\mathbf{y}}
#'
#' If the `method` is set to `"full_nr"`, the algorithm simultaneously updates both the mean and variance parameters.
#'
#' The variance-covariance matrix of the estimates can be computed using either the observed or
#' expected Fisher information, depending on the `vcov.type` specified.
#'
#' @export
#'
#' @examples
#' n <- 100
#' x1 <- rnorm(n)
#' x2 <- factor(sample(letters[1:3], n, TRUE))
#' x <- model.matrix(~ x1 + x2)
#' z1 <- factor(sample(letters[1:3], n, TRUE))
#' z <- model.matrix(~z1)
#'
#' b <- rnorm(ncol(x))
#' t <- rnorm(ncol(z))
#'
#' y <- rnorm(n, mean = x %*% b, sd = sqrt(exp(z %*% t)))
#'
#' start.list <- mvreg_start(y, x, z, start.s2 = "gamma")
#' b0 <- start.list$b0
#' t0 <- start.list$t0
#'
#' mvreg_fit(y, x, z, b0, t0)
mvreg_fit <- function(y, x, z, b0, t0, tol = 1e-10, maxit = 100, method = c("wls", "full_nr"), vcov.type = c("expected", "observed")) {

  if (!(maxit %% 1 == 0) | maxit < 0) {
    new.maxit <- round(abs(maxit))
    warning(paste0("maxit must be a positive integer, ", maxit, " taken as ", new.maxit))
    maxit <- new.maxit
  }
  if (maxit == 0) {
    new.maxit <- 100
    warning("maxit cannot be 0, set to default value of 100")
    maxit <- new.maxit
  }

  if(tol < 0) {
    new.tol <- abs(tol)
    warning(paste0("tol must be strictly positive, ", tol, " taken as ", new.tol))
    tol <- new.tol
  }
  if (tol == 0){
    new.tol <- 1e-10
    warning(paste0("tol must be strictly positive, ", tol, " set to default value of ", new.tol))
    tol <- new.tol
  }



  method <- match.arg(method)
  vcov.type <- match.arg(vcov.type)

  dev <- 2
  it <- 0L

  if (method == "wls") {
    while (any(abs(dev) > tol) && it < maxit) {
      t1 <- as.vector(t0 - solve(mvreg_hessian_s2(y, x, z, b0, t0, type = "observed")) %*% mvreg_gradient_s2(y, x, z, b0, t0))
      w <- as.vector(1 / exp(z %*% t1))
      b1 <- as.vector(solve(crossprod(x * w, x), crossprod(x * w, y)))
      dev <- c(b1, t1) - c(b0, t0)
      it <- it + 1L
      t0 <- t1
      b0 <- b1
    }
  } else {
    while (any(abs(dev) > tol) && it < maxit) {
      t1 <- as.vector(t0 - solve(mvreg_hessian_s2(y, x, z, b0, t0, type = "observed")) %*% mvreg_gradient_s2(y, x, z, b0, t0))
      b1 <- as.vector(b0 - solve(mvreg_hessian_mu(y, x, z, b0, t1, type = "observed")) %*% mvreg_gradient_mu(y, x, z, b0, t1))
      dev <- c(b1, t1) - c(b0, t0)
      it <- it + 1L
      t0 <- t1
      b0 <- b1
    }
  }

  theta0 <- c(b0, t0)
  names(theta0) <- c(colnames(x), colnames(z))
  vtheta <- solve(-mvreg_hessian(y, x, z, b0, t0, type = vcov.type))
  colnames(vtheta) <- rownames(vtheta) <- names(theta0)

  list(theta = theta0, b = b0, t = t0, vtheta = vtheta, it = it)
}
