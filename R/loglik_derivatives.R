#-------------------------------------------------------------------------------


#' Log-likelihood of a Heteroskedastic Linear Model
#'
#' This function computes the log-likelihood of a heteroskedastic linear model, where the variance component is modeled as a linear function of a set of explanatory variables, possibly different from those of the mean component.
#'
#' @param y A numeric vector representing the response variable.
#' @param x A numeric matrix of explanatory variables for the mean component of the model.
#' @param z A numeric matrix of explanatory variables for the variance component of the model.
#' @param b A numeric vector of parameters for the mean component.
#' @param t A numeric vector of parameters for the variance component.
#'
#' @return A numeric value representing the log-likelihood of the model for the given sample and parameter values.
#'
#' @details
#' The log-likelihood is computed under the assumption that the response variable follows a normal distribution with mean \eqn{\mathbf{x}_i' \boldsymbol{\beta}} and variance \eqn{\exp\left\{\mathbf{z}_i' \boldsymbol{\tau}\right\}}, where \eqn{\boldsymbol{\beta}} and \eqn{\boldsymbol{\tau}} are the parameters vectors for the mean and variance components, respectively.
#' The log-likelihood is computed as \deqn{\ell(\boldsymbol{\beta},\boldsymbol{\tau}) = -\dfrac{1}{2}\sum_{i=1}^n \left\{\log(2\pi) + \mathbf{z}_i'\boldsymbol{\tau} + \dfrac{\left(y_i - \mathbf{x}_i'\boldsymbol{\beta}\right)^2}{\exp\left\{\mathbf{z}_i'\boldsymbol{\tau}\right\}} \right\}}
#' @export
#'
#' @examples
#' # Fit a heteroskedastic linear model to the iris dataset
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, data = iris) # same formula for mean and variance
#'
#' # Extract model components
#' y <- mvreg_mod$y # Response variable
#' x <- mvreg_mod$x # Model matrix for the mean component
#' z <- mvreg_mod$z # Model matrix for the variance component
#' b <- coef(mvreg_mod, "mu") # Coefficients of the mean component
#' t <- coef(mvreg_mod, "s2") # Coefficients of the variance component
#'
#' # Compute log-likelihood
#' mvreg_loglik(y, x, z, b, t)
mvreg_loglik <- function(y, x, z, b, t) {
  eta.mu <- as.vector(x %*% b)
  eta.s2 <- as.vector(z %*% t)
  r <- -0.5 * sum(log(2 * pi) + eta.s2 + (y - eta.mu)^2 / exp(eta.s2))
  names(r) <- NULL
  r
}


# ------------------------------------------------------------------------------


#' Gradient of model log-likelihood with respect to mean coefficients (beta)
#'
#' This function computes the gradient of the log-likelihood of a heteroskedastic linear model with respect to the coefficients of the mean component (denoted as \eqn{\boldsymbol{\beta}}).
#'
#' @param y A numeric vector representing the response variable.
#' @param x A numeric matrix of explanatory variables for the mean component of the model.
#' @param z A numeric matrix of explanatory variables for the variance component of the model.
#' @param b A numeric vector of parameters for the mean component.
#' @param t A numeric vector of parameters for the variance component.
#'
#' @return A numeric vector representing the gradient of the log-likelihood with respect to the mean component coefficients \eqn{\boldsymbol{\beta}}.
#'
#' @details
#' The gradient is calculated based on the partial derivatives of the log-likelihood function with respect to the mean coefficients \eqn{\boldsymbol{\beta}}.
#' The \eqn{j}-th element is computed as \deqn{\dfrac{\partial \ell}{\partial \beta_j} = \sum_{i=1}^n \dfrac{y_i - \mathbf{x}_i'\boldsymbol{\beta}}{\exp\left\{\mathbf{z}_i'\boldsymbol{\tau}\right\}}x_{ij}}
#'
#' @export
#'
#' @examples
#' # Fit a heteroskedastic linear model to the iris dataset
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, data = iris) # same formula for mean and variance
#'
#' # Extract model components
#' y <- mvreg_mod$y # Response variable
#' x <- mvreg_mod$x # Model matrix for the mean component
#' z <- mvreg_mod$z # Model matrix for the variance component
#' b <- coef(mvreg_mod, "mu") # Coefficients of the mean component
#' t <- coef(mvreg_mod, "s2") # Coefficients of the variance component
#'
#' # Compute gradient with respect to the mean coefficients
#' mvreg_gradient_mu(y, x, z, b, t)
mvreg_gradient_mu <- function(y, x, z, b, t) {
  eta.mu <- as.vector(x %*% b)
  eta.s2 <- as.vector(z %*% t)
  const <- (y - eta.mu) / exp(eta.s2)
  r <- apply(x, 2, function(x) sum(const * x))
  names(r) <- colnames(x)
  r
}


# ------------------------------------------------------------------------------


#' Hessian of model log-likelihood with Respect to mean coefficients (beta)
#'
#' This function computes the Hessian matrix of the log-likelihood of a heteroskedastic linear model with respect to the coefficients of the mean component.
#' The Hessian can be calculated as either the observed or expected matrix, depending on the `type` argument.
#'
#' @param y A numeric vector representing the response variable.
#' @param x A numeric matrix of explanatory variables for the mean component of the model.
#' @param z A numeric matrix of explanatory variables for the variance component of the model.
#' @param b A numeric vector of parameters for the mean component.
#' @param t A numeric vector of parameters for the variance component.
#' @param type A character string specifying whether to return the observed or expected Hessian matrix.
#' The default is "observed".
#'
#' @return A numeric matrix representing the Hessian of the log-likelihood with respect to the mean component coefficients \eqn{\boldsymbol{\beta}}.
#'
#'
#' @details
#' The Hessian is calculated as the matrix of second-order partial derivatives of the log-likelihood function with respect to \eqn{\boldsymbol{\beta}}.
#' The \eqn{(j, l)}-th element of the Hessian is computed using the formula:
#' \deqn{\dfrac{\partial^2 \ell}{\partial \beta_j \partial \beta_l} = -\sum_{i=1}^n \dfrac{x_{ij}x_{il}}{\exp\left\{\mathbf{z}_i'\boldsymbol{\tau}\right\}}}
#' In this case, the observed hessian and the expected hessian coincide. The argument `type` is maintained just for consistency.
#'
#'
#' @export
#'
#' @examples
#' # Fit a heteroskedastic linear model to the iris dataset
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, data = iris) # same formula for mean and variance
#'
#' # Extract model components
#' y <- mvreg_mod$y # Response variable
#' x <- mvreg_mod$x # Model matrix for the mean component
#' z <- mvreg_mod$z # Model matrix for the variance component
#' b <- coef(mvreg_mod, "mu") # Coefficients of the mean component
#' t <- coef(mvreg_mod, "s2") # Coefficients of the variance component
#'
#' # Compute Hessian with respect to the mean coefficients
#' mvreg_hessian_mu(y, x, z, b, t, type = "observed")
#' mvreg_hessian_mu(y, x, z, b, t, type = "expected")
mvreg_hessian_mu <- function(y, x, z, b, t, type = c("observed", "expected")) {
  type <- match.arg(type)
  k <- ncol(x)
  w <- as.vector(1/exp(z%*%t))
  h <- -crossprod(x, x*w)
  colnames(h) <- rownames(h) <- colnames(x)
  h
}


# ------------------------------------------------------------------------------


#' Gradient of model log-likelihood with respect to variance coefficients (tau)
#'
#' This function computes the gradient of the log-likelihood of a heteroskedastic linear model with respect to the coefficients of the variance component (denoted as \eqn{\boldsymbol{\tau}}).
#'
#' @param y A numeric vector representing the response variable.
#' @param x A numeric matrix of explanatory variables for the mean component of the model.
#' @param z A numeric matrix of explanatory variables for the variance component of the model.
#' @param b A numeric vector of parameters for the mean component.
#' @param t A numeric vector of parameters for the variance component.
#'
#' @return A numeric vector representing the gradient of the log-likelihood with respect to the variance component coefficients \eqn{\boldsymbol{\tau}}.
#'
#' @details
#' The gradient is calculated based on the partial derivatives of the log-likelihood function with respect to the variance coefficients \eqn{\boldsymbol{\tau}}.
#' The \eqn{j}-th element of the gradient is computed as:
#' \deqn{\dfrac{\partial \ell}{\partial \tau_j} = -\dfrac{1}{2}\sum_{i=1}^n\left\{ 1 - \dfrac{\left(y_i - \mathbf{x}_i'\boldsymbol{\beta}\right)^2}{\exp\left\{\mathbf{z}_i'\boldsymbol{\tau}\right\}}z_{ij}\right\}}
#'
#' @export
#'
#' @examples
#' # Fit a heteroskedastic linear model to the iris dataset
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, data = iris) # same formula for mean and variance
#'
#' # Extract model components
#' y <- mvreg_mod$y # Response variable
#' x <- mvreg_mod$x # Model matrix for the mean component
#' z <- mvreg_mod$z # Model matrix for the variance component
#' b <- coef(mvreg_mod, "mu") # Coefficients of the mean component
#' t <- coef(mvreg_mod, "s2") # Coefficients of the variance component
#'
#' # Compute gradient with respect to the variance coefficients
#' mvreg_gradient_s2(y, x, z, b, t)
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


#' Hessian of model log-likelihood with respect to variance coefficients (tau)
#'
#' This function computes the Hessian matrix of the log-likelihood of a heteroskedastic linear model with respect to the coefficients of the variance component.
#' The Hessian can be calculated as either the observed or expected matrix, depending on the `type` argument.
#'
#' @param y A numeric vector representing the response variable.
#' @param x A numeric matrix of explanatory variables for the mean component of the model.
#' @param z A numeric matrix of explanatory variables for the variance component of the model.
#' @param b A numeric vector of parameters for the mean component.
#' @param t A numeric vector of parameters for the variance component.
#' @param type A character string specifying whether to return the observed or expected Hessian matrix.
#' The default is "observed".
#'
#' @return A numeric matrix representing the Hessian of the log-likelihood with respect to the variance component coefficients \eqn{\boldsymbol{\tau}}.
#'
#' @details
#' The observed Hessian is calculated as the matrix of second-order partial derivatives of the log-likelihood function with respect to \eqn{\boldsymbol{\tau}}.
#' The \eqn{(j, l)}-th element of the observed Hessian is computed using the formula:
#' \deqn{\dfrac{\partial^2 \ell}{\partial \tau_j \partial \tau_l} = -\dfrac{1}{2}\sum_{i=1}^n \dfrac{\left(y_i - \mathbf{x}_i'\boldsymbol{\beta}\right)^2}{\exp\left\{\mathbf{z}_i'\boldsymbol{\tau}\right\}}z_{ij}z_{il} }
#'
#' The \eqn{(j, l)}-th element of the expected Hessian is computed calculating the expectation of the \eqn{(j, l)}-th element of the observed Hessian:
#' \deqn{\mathbb{E}\left(\dfrac{\partial^2 \ell}{\partial \tau_j \partial \tau_l}\right) = -\dfrac{1}{2}\sum_{i=1}^n z_{ij}z_{il}}
#'
#' @export
#'
#' @examples
#' # Fit a heteroskedastic linear model to the iris dataset
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, data = iris) # same formula for mean and variance
#'
#' # Extract model components
#' y <- mvreg_mod$y # Response variable
#' x <- mvreg_mod$x # Model matrix for the mean component
#' z <- mvreg_mod$z # Model matrix for the variance component
#' b <- coef(mvreg_mod, "mu") # Coefficients of the mean component
#' t <- coef(mvreg_mod, "s2") # Coefficients of the variance component
#'
#' # Compute Hessian with respect to the variance coefficients
#' mvreg_hessian_s2(y, x, z, b, t, type = "observed")
#' mvreg_hessian_s2(y, x, z, b, t, type = "expected")
mvreg_hessian_s2 <- function(y, x, z, b, t, type = c("observed", "expected")) {
  type <- match.arg(type)
  p <- ncol(z)
  if (type == "observed") {
    eta.mu <- as.vector(x %*% b)
    eta.s2 <- as.vector(z %*% t)
    const <- (y - eta.mu)^2 * exp(-eta.s2)
  }
  else if (type == "expected") {
    const <- 1
  }
  h <- -0.5*crossprod(z, z*const)
  colnames(h) <- rownames(h) <- colnames(z)
  h
}


# ------------------------------------------------------------------------------


#' Hessian of model log-likelihood with respect to mean and variance coefficients (beta and tau)
#'
#' This function computes the Hessian matrix of the log-likelihood of a heteroskedastic linear model
#' with respect to the coefficients of the mean component (denoted as \eqn{\boldsymbol{\beta}})
#' and the coefficients of the variance component (denoted as \eqn{\boldsymbol{\tau}}).
#' The Hessian can be calculated as either the observed or expected matrix, depending on the `type` argument.
#'
#' @param y A numeric vector representing the response variable.
#' @param x A numeric matrix of explanatory variables for the mean component of the model.
#' @param z A numeric matrix of explanatory variables for the variance component of the model.
#' @param b A numeric vector of parameters for the mean component.
#' @param t A numeric vector of parameters for the variance component.
#' @param type A character string specifying whether to return the observed or expected Hessian matrix.
#' The default is "observed".
#'
#' @return A numeric matrix representing the Hessian of the log-likelihood with respect to the mean and variance component coefficients.
#'
#' @details
#' The observed Hessian is calculated as the matrix of mixed second-order partial derivatives of the log-likelihood function
#' with respect to \eqn{\boldsymbol{\beta}} and \eqn{\boldsymbol{\tau}}.
#' The \eqn{(j, l)}-th element of the observed Hessian is computed using the formula:
#' \deqn{\dfrac{\partial^2 \ell}{\partial \beta_j \partial \tau_j} = -\sum_{i=1}^n \dfrac{y_i - \mathbf{x}_i'\boldsymbol{\beta}}{\exp\left\{\mathbf{z}_i'\boldsymbol{\tau}\right\}}x_{ij}z_{il}}
#'
#' The expected Hessian is a zero matrix.
#'
#' @export
#'
#' @examples
#' # Fit a heteroskedastic linear model to the iris dataset
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, data = iris) # same formula for mean and variance
#'
#' # Extract model components
#' y <- mvreg_mod$y # Response variable
#' x <- mvreg_mod$x # Model matrix for the mean component
#' z <- mvreg_mod$z # Model matrix for the variance component
#' b <- coef(mvreg_mod, "mu") # Coefficients of the mean component
#' t <- coef(mvreg_mod, "s2") # Coefficients of the variance component
#'
#' # Compute Hessian with respect to mean and variance coefficients
#' mvreg_hessian_mus2(y, x, z, b, t, type = "observed")
#' mvreg_hessian_mus2(y, x, z, b, t, type = "expected")
mvreg_hessian_mus2 <- function(y, x, z, b, t, type = c("observed", "expected")) {
  type <- match.arg(type)
  if (type == "observed") {
    eta.mu <- as.vector(x %*% b)
    eta.s2 <- as.vector(z %*% t)
    const <- (y - eta.mu) * exp(-eta.s2)
    h <- -crossprod(x, const*z)
  }
  else if (type == "expected") {
    h <- matrix(0, ncol(x), ncol(z))
  }
  rownames(h) <- colnames(x)
  colnames(h) <- colnames(z)
  h
}


# ------------------------------------------------------------------------------


#' Full hessian of model log-likelihood
#'
#' This function computes the full hessian matrix of the log-likelihood of a heteroskedastic linear model.
#' The hessian can be calculated as either the observed or expected matrix, depending on the `type` argument.
#'
#' @param y A numeric vector representing the response variable.
#' @param x A numeric matrix of explanatory variables for the mean component of the model.
#' @param z A numeric matrix of explanatory variables for the variance component of the model.
#' @param b A numeric vector of parameters for the mean component.
#' @param t A numeric vector of parameters for the variance component.
#' @param type A character string specifying whether to return the observed or expected Hessian matrix.
#' The default is "observed".
#'
#' @return A numeric matrix representing the full hessian of the log-likelihood with respect to the mean and variance component coefficients.
#'
#' @details
#' The full Hessian is constructed by combining the hessians of the mean and variance components,
#' as well as the mixed partial derivatives between them. The structure of the returned matrix is:
#' \deqn{\mathcal{H}(\boldsymbol{\beta}, \boldsymbol{\tau}) =
#' \begin{bmatrix}
#' \mathcal{H}_{\boldsymbol{\beta}}(\boldsymbol{\beta}, \boldsymbol{\tau}) & \mathcal{H}_{\boldsymbol{\beta} \boldsymbol{\tau}}(\boldsymbol{\beta}, \boldsymbol{\tau}) \\
#' \mathcal{H}_{\boldsymbol{\beta} \boldsymbol{\tau}}'(\boldsymbol{\beta}, \boldsymbol{\tau}) & \mathcal{H}_{\boldsymbol{\tau}}(\boldsymbol{\beta}, \boldsymbol{\tau})
#' \end{bmatrix}}
#'
#' @export
#'
#' @examples
#' # Fit a heteroskedastic linear model to the iris dataset
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, data = iris) # same formula for mean and variance
#'
#' # Extract model components
#' y <- mvreg_mod$y # Response variable
#' x <- mvreg_mod$x # Model matrix for the mean component
#' z <- mvreg_mod$z # Model matrix for the variance component
#' b <- coef(mvreg_mod, "mu") # Coefficients of the mean component
#' t <- coef(mvreg_mod, "s2") # Coefficients of the variance component
#'
#' # Compute the full Hessian
#' mvreg_hessian(y, x, z, b, t, type = "observed")
#' mvreg_hessian(y, x, z, b, t, type = "expected")
mvreg_hessian <- function(y, x, z, b, t, type = c("observed", "expected")) {
  type <- match.arg(type)
  hbb <- mvreg_hessian_mu(y, x, z, b, t, type)
  htt <- mvreg_hessian_s2(y, x, z, b, t, type)
  hbt <- mvreg_hessian_mus2(y, x, z, b, t, type)
  rbind(cbind(hbb, hbt), cbind(t(hbt), htt))
}
