#-------------------------------------------------------------------------------


#' Fitter function for mvreg
#'
#' @param y Response variable vector.
#' @param x Design matrix for mean component.
#' @param z Design matrix for variance component.
#' @param b0 Starting values for parameters of mean component.
#' @param t0 Starting values for parameters of variance component.
#' @param tol Positive value indicating what is the minimum difference between parameter estimates between two iterations to stop the algorithm.
#' @param maxit Integer value indicating the maximum number of iteration.
#' @param method Method chosen for estimation of parameters of mean component
#' @param vcov.type A string to specify whether to use observed or expected Fisher information matrix in order to compute variance-covariance matrix of estimates
#'
#' @return Estimate of coefficients and variance-covariance matrix
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

  method <- match.arg(method)
  vcov.type <- match.arg(vcov.type)

  p <- ncol(z)
  k <- ncol(x)

  dev <- 2
  it <- 0L

  if (method == "wls") {
    while (any(abs(dev) > tol) && it < maxit) {
      t1 <- as.vector(t0 - solve(d2ldt(y, x, z, b0, t0, type = "observed")) %*% dldt(y, x, z, b0, t0))
      w <- as.vector(1/exp(z%*%t1))
      b1 <- as.vector(solve(crossprod(x*w, x), crossprod(x*w, y)))
      dev <- c(b1, t1) - c(b0, t0)
      it <- it + 1L
      t0 <- t1
      b0 <- b1
    }
  } else {
    while (any(abs(dev) > tol) && it < maxit) {
      t1 <- as.vector(t0 - solve(d2ldt(y, x, z, b0, t0, type = "observed")) %*% dldt(y, x, z, b0, t0))
      b1 <- as.vector(b0 - solve(d2ldb(y, x, z, b0, t1, type = "observed")) %*% dldb(y, x, z, b0, t1))
      dev <- c(b1, t1) - c(b0, t0)
      it <- it + 1L
      t0 <- t1
      b0 <- b1
    }
  }

  theta0 <- c(b0, t0)
  names(theta0) <- c(colnames(x), colnames(z))
  vtheta <- solve(-d2l(y, x, z, b0, t0, type = vcov.type))
  colnames(vtheta) <- rownames(vtheta) <- names(theta0)

  list(theta = theta0, b = b0, t = t0, vtheta = vtheta, it = it)
}

