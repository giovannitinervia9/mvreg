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
#'
#' @return Estimate of coefficients and variance-covariance matrix
#' @export
#'
#' @examples
mvreg_fit <- function(y, x, z, b0, t0, tol = 1e-10, maxit = 100) {
  p <- ncol(z)
  k <- ncol(x)

  gb <- function(b, t) {
    dldb(y, x, z, b, t)
  }

  gt <- function(b, t) {
    dldt(y, x, z, b, t)
  }


  hb <- function(b, t) {
    d2ldb(y, x, z, b, t)
  }


  ht <- function(b, t) {
    d2ldt(y, x, z, b, t)
  }


  h <- function(theta) {
    b <- theta[1:k]
    t <- theta[(k + 1):length(theta)]
    d2l(y, x, z, b, t)
  }


  dev <- 2
  it <- 0L

  while (any(abs(dev) > tol) && it < maxit) {
    t1 <- as.vector(t0 - solve(ht(b0, t0)) %*% gt(b0, t0))
    b1 <- as.vector(b0 - solve(hb(b0, t1)) %*% gb(b0, t1))
    dev <- c(b1, t1) - c(b0, t0)
    it <- it + 1L
    t0 <- t1
    b0 <- b1
  }

  theta0 <- c(b0, t0)
  names(theta0) <- c(colnames(x), colnames(z))
  vtheta <- solve(-h(theta0))
  colnames(vtheta) <- rownames(vtheta) <- names(theta0)

  list(theta = theta0, b = b0, t = t0, vtheta = vtheta, it = it)
}
