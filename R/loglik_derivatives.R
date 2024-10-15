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
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, data = iris) # same formula for mean and variance
#' y <- mvreg_mod$y # response variable
#' x <- mvreg_mod$x # model.matrix for mean component
#' z <- mvreg_mod$z # model.matrix for variance component
#' b <- coef(mvreg_mod, "mu") # coefficients of mean component
#' t <- coef(mvreg_mod, "s2") # coefficients of variance component
#'
#' # loglikelihood
#' mvreg_loglik(y, x, z, b, t)
#'
mvreg_loglik <- function(y, x, z, b, t) {
  eta.mu <- as.vector(x %*% b)
  eta.s2 <- as.vector(z %*% t)
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
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, data = iris) # same formula for mean and variance
#' y <- mvreg_mod$y # response variable
#' x <- mvreg_mod$x # model.matrix for mean component
#' z <- mvreg_mod$z # model.matrix for variance component
#' b <- coef(mvreg_mod, "mu") # coefficients of mean component
#' t <- coef(mvreg_mod, "s2") # coefficients of variance component
#'
#' # first derivatives
#' mvreg_gradient_mu(y, x, z, b, t) # w.r.t. mean coefficients
#' mvreg_gradient_s2(y, x, z, b, t) # w.r.t. variance coefficients
#'
mvreg_gradient_mu <- function(y, x, z, b, t) {
  eta.mu <- as.vector(x %*% b)
  eta.s2 <- as.vector(z %*% t)
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
#' @param type A string to specify whether to output the observed or expected matrix
#'
#' @return Hessian of loglikelihood with respect to beta.
#' @export
#'
#' @examples
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, data = iris) # same formula for mean and variance
#' y <- mvreg_mod$y # response variable
#' x <- mvreg_mod$x # model.matrix for mean component
#' z <- mvreg_mod$z # model.matrix for variance component
#' b <- coef(mvreg_mod, "mu") # coefficients of mean component
#' t <- coef(mvreg_mod, "s2") # coefficients of variance component
#'
#' # second derivatives
#' mvreg_hessian_mu(y, x, z, b, t) # w.r.t. mean coefficients
#' mvreg_hessian_s2(y, x, z, b, t) # w.r.t. variance coefficients
#' mvreg_hessian_mus2(y, x, z, b, t) # w.r.t. mean coefficients and variance coefficients
#' mvreg_hessian(y, x, z, b, t) # full hessian
#'
mvreg_hessian_mu <- function(y, x, z, b, t, type = c("observed", "expected")) {
  type <- match.arg(type)
  h <- matrix(NA, ncol(x), ncol(x))
  eta.s2 <- as.vector(z %*% t)
  const <- exp(eta.s2)

  for (j in seq_len(ncol(x))) {
    for (l in seq_len(ncol(x))) {
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
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, data = iris) # same formula for mean and variance
#' y <- mvreg_mod$y # response variable
#' x <- mvreg_mod$x # model.matrix for mean component
#' z <- mvreg_mod$z # model.matrix for variance component
#' b <- coef(mvreg_mod, "mu") # coefficients of mean component
#' t <- coef(mvreg_mod, "s2") # coefficients of variance component
#'
#' # first derivatives
#' mvreg_gradient_mu(y, x, z, b, t) # w.r.t. mean coefficients
#' mvreg_gradient_s2(y, x, z, b, t) # w.r.t. variance coefficients
#'
mvreg_gradient_s2 <- function(y, x, z, b, t) {
  eta.mu <- as.vector(x %*% b)
  eta.s2 <- as.vector(z %*% t)
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
#' @param type A string to specify whether to output the observed or expected matrix
#'
#' @return Hessian of loglikelihood with respect to tau.
#' @export
#'
#' @examples
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, data = iris) # same formula for mean and variance
#' y <- mvreg_mod$y # response variable
#' x <- mvreg_mod$x # model.matrix for mean component
#' z <- mvreg_mod$z # model.matrix for variance component
#' b <- coef(mvreg_mod, "mu") # coefficients of mean component
#' t <- coef(mvreg_mod, "s2") # coefficients of variance component
#'
#' # second derivatives
#' mvreg_hessian_mu(y, x, z, b, t) # w.r.t. mean coefficients
#' mvreg_hessian_s2(y, x, z, b, t) # w.r.t. variance coefficients
#' mvreg_hessian_mus2(y, x, z, b, t) # w.r.t. mean coefficients and variance coefficients
#' mvreg_hessian(y, x, z, b, t) # full hessian
#'
mvreg_hessian_s2 <- function(y, x, z, b, t, type = c("observed", "expected")) {
  type <- match.arg(type)
  h <- matrix(NA, ncol(z), ncol(z))

  if (type == "observed") {
    eta.mu <- as.vector(x %*% b)
    eta.s2 <- as.vector(z %*% t)
    const <- (y - eta.mu)^2 * exp(-eta.s2)
  } else if (type == "expected") {
    const <- 1
  }


  for (j in seq_len(ncol(z))) {
    for (l in seq_len(ncol(z))) {
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
#' @param type A string to specify whether to output the observed or expected matrix
#'
#' @return Hessian of loglikelihood with respect to beta and tau.
#' @export
#'
#' @examples
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, data = iris) # same formula for mean and variance
#' y <- mvreg_mod$y # response variable
#' x <- mvreg_mod$x # model.matrix for mean component
#' z <- mvreg_mod$z # model.matrix for variance component
#' b <- coef(mvreg_mod, "mu") # coefficients of mean component
#' t <- coef(mvreg_mod, "s2") # coefficients of variance component
#'
#' # second derivatives
#' mvreg_hessian_mu(y, x, z, b, t) # w.r.t. mean coefficients
#' mvreg_hessian_s2(y, x, z, b, t) # w.r.t. variance coefficients
#' mvreg_hessian_mus2(y, x, z, b, t) # w.r.t. mean coefficients and variance coefficients
#' mvreg_hessian(y, x, z, b, t) # full hessian
#'
mvreg_hessian_mus2 <- function(y, x, z, b, t, type = c("observed", "expected")) {
  type <- match.arg(type)


  if (type == "observed") {
    h <- matrix(NA, ncol(x), ncol(z))
    eta.mu <- as.vector(x %*% b)
    eta.s2 <- as.vector(z %*% t)
    const <- (y - eta.mu) * exp(-eta.s2)

    for (j in seq_len(ncol(x))) {
      for (l in seq_len(ncol(z))) {
        h[j, l] <- -sum(const * x[, j] * z[, l])
      }
    }
  } else if (type == "expected") {
    h <- matrix(0, ncol(x), ncol(z))
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
#' @param type A string to specify whether to output the observed or expected matrix
#'
#' @return Full hessian of model loglikelihood.
#' @export
#'
#' @examples
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, data = iris) # same formula for mean and variance
#' y <- mvreg_mod$y # response variable
#' x <- mvreg_mod$x # model.matrix for mean component
#' z <- mvreg_mod$z # model.matrix for variance component
#' b <- coef(mvreg_mod, "mu") # coefficients of mean component
#' t <- coef(mvreg_mod, "s2") # coefficients of variance component
#'
#' # second derivatives
#' mvreg_hessian_mu(y, x, z, b, t) # w.r.t. mean coefficients
#' mvreg_hessian_s2(y, x, z, b, t) # w.r.t. variance coefficients
#' mvreg_hessian_mus2(y, x, z, b, t) # w.r.t. mean coefficients and variance coefficients
#' mvreg_hessian(y, x, z, b, t) # full hessian
#'
mvreg_hessian <- function(y, x, z, b, t, type = c("observed", "expected")) {
  type <- match.arg(type)
  hbb <- mvreg_hessian_mu(y, x, z, b, t, type)
  htt <- mvreg_hessian_s2(y, x, z, b, t, type)
  hbt <- mvreg_hessian_mus2(y, x, z, b, t, type)
  rbind(cbind(hbb, hbt), cbind(t(hbt), htt))
}
