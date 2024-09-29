#-------------------------------------------------------------------------------


#' Loglikelihood of heteroskedastic linear model
#'
#' @param y Vector of response variable.
#' @param x Matrix of explanatory variables for mean component.
#' @param z Matrix of explanatory variables for variance component.
#' @param b Vector of parameters for mean component.
#' @param t Vector of parameters for variance component.
#'
#' @return Value of loglikelihood for a given sample and parameters' values.
#' @export
#'
#' @examples
ll <- function(y, x, z, b, t) {
  eta.mu <- x %*% b
  eta.s2 <- z %*% t
  r <- -0.5 * sum(log(2 * pi) + eta.s2 + (y - eta.mu)^2 / exp(eta.s2))
  names(r) <- NULL
  r
}


# ------------------------------------------------------------------------------


#' Gradient of model loglikelihood with respect to beta
#'
#' @param y Vector of response variable.
#' @param x Matrix of explanatory variables for mean component.
#' @param z Matrix of explanatory variables for variance component.
#' @param b Vector of parameters for mean component.
#' @param t Vector of parameters for variance component.
#'
#' @return Gradient with respect to beta.
#' @export
#'
#' @examples
dldb <- function(y, x, z, b, t) {
  eta.mu <- x %*% b
  eta.s2 <- z %*% t
  const <- (y - eta.mu) / exp(eta.s2)
  r <- apply(x, 2, function(x) sum(const * x))
  names(r) <- colnames(x)
  r
}


# ------------------------------------------------------------------------------


#' Hessian of model loglikelihood with respect to beta
#'
#' @param y Vector of response variable.
#' @param x Matrix of explanatory variables for mean component.
#' @param z Matrix of explanatory variables for variance component.
#' @param b Vector of parameters for mean component.
#' @param t Vector of parameters for variance component.
#'
#' @return Hessian of loglikelihood with respect to beta.
#' @export
#'
#' @examples
d2ldb <- function(y, x, z, b, t) {
  h <- matrix(NA, ncol(x), ncol(x))
  eta.s2 <- z %*% t
  const <- exp(eta.s2)

  for (j in 1:ncol(x)) {
    for (l in 1:ncol(x)) {
      h[j, l] <- -sum(x[, j] * x[, l] / const)
    }
  }

  colnames(h) <- rownames(h) <- colnames(x)
  h
}


# ------------------------------------------------------------------------------


#' Gradient of model loglikelihood with respect to tau
#'
#' @param y Vector of response variable.
#' @param x Matrix of explanatory variables for mean component.
#' @param z Matrix of explanatory variables for variance component.
#' @param b Vector of parameters for mean component.
#' @param t Vector of parameters for variance component.
#'
#' @return Gradient with respect to tau.
#' @export
#'
#' @examples
dldt <- function(y, x, z, b, t) {
  eta.mu <- x %*% b
  eta.s2 <- z %*% t
  const <- 1 - (y - eta.mu)^2 * exp(-eta.s2)

  r <- apply(z, 2, function(z) {
    -0.5 * sum(const * z)
  })

  names(r) <- colnames(z)
  r
}


# ------------------------------------------------------------------------------


#' Hessian of model loglikelihood with respect to tau
#'
#' @param y Vector of response variable.
#' @param x Matrix of explanatory variables for mean component.
#' @param z Matrix of explanatory variables for variance component.
#' @param b Vector of parameters for mean component.
#' @param t Vector of parameters for variance component.
#'
#' @return Hessian of loglikelihood with respect to tau.
#' @export
#'
#' @examples
d2ldt <- function(y, x, z, b, t) {
  h <- matrix(NA, ncol(z), ncol(z))
  eta.mu <- x %*% b
  eta.s2 <- z %*% t
  const <- (y - eta.mu)^2 * exp(-eta.s2)

  for (j in 1:ncol(z)) {
    for (l in 1:ncol(z)) {
      h[j, l] <- -0.5 * sum(const * z[, j] * z[, l])
    }
  }

  colnames(h) <- rownames(h) <- colnames(z)
  h
}


# ------------------------------------------------------------------------------


#' Hessian of model loglikelihood with respect to beta and tau
#'
#' @param y Vector of response variable.
#' @param x Matrix of explanatory variables for mean component.
#' @param z Matrix of explanatory variables for variance component.
#' @param b Vector of parameters for mean component.
#' @param t Vector of parameters for variance component.
#'
#' @return Hessian of loglikelihood with respect to beta and tau.
#' @export
#'
#' @examples
d2ldbdt <- function(y, x, z, b, t) {
  h <- matrix(NA, ncol(x), ncol(z))

  eta.mu <- x %*% b
  eta.s2 <- z %*% t
  const <- (y - eta.mu) * exp(-eta.s2)

  for (j in 1:ncol(x)) {
    for (l in 1:ncol(z)) {
      h[j, l] <- -sum(const * x[, j] * z[, l])
    }
  }

  rownames(h) <- colnames(x)
  colnames(h) <- colnames(z)
  h
}


# ------------------------------------------------------------------------------


#' Full hessian of model loglikelihood
#'
#' @param y Vector of response variable.
#' @param x Matrix of explanatory variables for mean component.
#' @param z Matrix of explanatory variables for variance component.
#' @param b Vector of parameters for mean component.
#' @param t Vector of parameters for variance component.
#'
#' @return Full hessian of model loglikelihood.
#' @export
#'
#' @examples
d2l <- function(y, x, z, b, t) {
  hbb <- d2ldb(y, x, z, b, t)
  htt <- d2ldt(y, x, z, b, t)
  hbt <- d2ldbdt(y, x, z, b, t)
  rbind(cbind(hbb, hbt), cbind(t(hbt), htt))
}
