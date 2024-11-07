#-------------------------------------------------------------------------------


#' Heteroscedastic linear model with variance as a function of linear predictor.
#'
#' @param formula.mu A two-sided formula describing the model for mean component.
#' @param formula.s2 A one-sided formula describing the model for variance component.
#' @param data An optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param tol Positive value indicating what is the minimum difference between parameter estimates between two iterations to stop the algorithm.
#' @param maxit Integer value indicating the maximum number of iteration.
#' @param start.s2 A character vector indicating how to select initial values for variance component parameters. See the documentation of \code{\link{mvreg_start}} for details.
#' @param method A character vector indicating the method chosen for estimation of parameters of mean component. See the documentation of \code{\link{mvreg_fit}} for details.
#' @param vcov.type A character vector to specify whether to use observed or expected Fisher information matrix in order to compute variance-covariance matrix of estimates. See the documentation of \code{\link{mvreg_fit}}, \code{\link{mvreg_hessian_mu}}, \code{\link{mvreg_hessian_s2}}, \code{\link{mvreg_hessian_mus2}} and \code{\link{mvreg_hessian}} for details.
#'
#' @return An object of class mvreg containing:
#' \item{coefficients}{The estimated coefficients for both the mean and variance components.}
#' \item{coefficients.mu}{The estimated coefficients for the mean component.}
#' \item{coefficients.s2}{The estimated coefficients for the variance component.}
#' \item{vcov}{The variance-covariance matrix of all the estimated coefficients.}
#' \item{vcov.mu}{The variance-covariance matrix of the coefficients for the mean component.}
#' \item{vcov.s2}{The variance-covariance matrix of the coefficients for the variance component.}
#' \item{logLik}{The log-likelihood of the fitted model.}
#' \item{fit.mu}{The fitted values for the mean component.}
#' \item{fit.log.s2}{The fitted values for the log of the variance component.}
#' \item{fit.s2}{The fitted values for the variance component.}
#' \item{residuals}{The residuals from the mean component.}
#' \item{it}{The number of iterations the algorithm took to converge.}
#' \item{start}{Starting values of the estimates for both mean and variance components.}
#' \item{y}{The response vector.}
#' \item{xd}{The design matrix for the mean component.}
#' \item{zd}{The design matrix for the variance component.}
#' \item{nobs}{The number of observations in the dataset.}
#' \item{df.residual}{The residual degrees of freedom.}
#' \item{call}{The matched function call.}
#' \item{response}{The name of the response variable.}
#' \item{colx}{The column names of the explanatory variables for the mean component.}
#' \item{colz}{The column names of the explanatory variables for the variance component.}
#' \item{formula.mu}{The formula used for the mean component.}
#' \item{formula.s2}{The formula used for the variance component.}
#' \item{terms.mu}{Terms attribute of model.frame() of mean component}
#' \item{terms.s2}{Terms attribute of model.frame() of variance component}
#'
#' @export
#'
#' @importFrom stats terms model.frame model.response var as.formula
#' @examples
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, data = iris) # same formula for mean and variance
#' mvreg_mod1 <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris) # different formulas
#' summary(mvreg_mod)
#' summary(mvreg_mod1)
#'
#' # coef
#' coef(mvreg_mod)
#' coef(mvreg_mod, partition = "mu")
#'
#' # vcov
#' vcov(mvreg_mod)
#' vcov(mvreg_mod, partition = "s2")
#'
#' # fitted
#' fitted(mvreg_mod)
#' fitted(mvreg_mod, type = "mu")
#' fitted(mvreg_mod, type = "log.s2")
#' fitted(mvreg_mod, type = "s2")
#'
#' # logLik
#' logLik(mvreg_mod)
#' logLik(mvreg_mod1)
#'
#' # predict without newdata
#' predict(mvreg_mod1)
#' predict(mvreg_mod1, type = "mu")
#' predict(mvreg_mod1, type = "log.s2", se.fit = TRUE)
#' predict(mvreg_mod1, type = "s2", interval = TRUE)
#' predict(mvreg_mod1, se.fit = TRUE, interval = TRUE, level = 0.99)
#'
#' # predict with newdata
#' unique(c(mvreg_mod1$colx, mvreg_mod1$colz)) # getting names of explanatory variables
#'
#' newdata <- data.frame(
#'   Species = levels(iris$Species),
#'   Sepal.Width = c(min(iris$Sepal.Width), mean(iris$Sepal.Width), max(iris$Sepal.Width))
#' )
#'
#' newdata <- expand.grid(newdata)
#'
#' predict(mvreg_mod1, newdata = newdata, se.fit = TRUE, interval = TRUE)
mvreg <- function(formula.mu,
                  formula.s2 = NULL,
                  data = NULL,
                  tol = 1e-10,
                  maxit = 100L,
                  start.s2 = c("residuals", "gamma", "zero"),
                  method = c("wls", "full_nr"),
                  vcov.type = c("expected", "observed")) {
  cl <- match.call()
  method <- match.arg(method)
  vcov.type <- match.arg(vcov.type)
  start.s2 <- match.arg(start.s2)

  if (is.null(data)) {
    data <- environment()
  }

  mf.mu <- model.frame(formula.mu, data)
  terms.mu <- attr(mf.mu, "terms")
  removed.formula.mu <- as.vector(attr(mf.mu, "na.action"))
  if (!is.null(removed.formula.mu)) {
    stop("mvreg() can't handle NAs")
  }
  colx <- colnames(mf.mu)[-1]
  response <- colnames(mf.mu)[1]
  y <- model.response(mf.mu)
  xd <- model.matrix(attr(mf.mu, "terms"), mf.mu)

  if (is.null(formula.s2)) {
    formula.s2 <- formula.mu
  } else {
    if (length(formula.s2) == 2) {
      if (sum(grepl(response, formula.s2)) == 1) {
        stop("response variable cannot appear on right-hand side of the formula")
      } else {
        formula.s2 <- as.formula(paste0(response, " ", paste0(formula.s2, collapse = " ")))
      }
    } else {
      formula.s2 <- formula.s2
    }

  }

  mf.s2 <- model.frame(formula.s2, data)
  terms.s2 <- attr(mf.s2, "terms")
  removed.formula.s2 <- as.vector(attr(mf.s2, "na.action"))
  if (!is.null(removed.formula.s2)) {
    stop("mvreg() can't handle NAs")
  }
  colz <- colnames(mf.s2)[-1]
  zd <- model.matrix(attr(mf.s2, "terms"), mf.s2)


  cl$formula.s2 <- formula.s2


  colnames(xd)[which(colnames(xd) == "(Intercept)")] <- "const"
  colnames(zd)[which(colnames(zd) == "(Intercept)")] <- "const"

  colnames(xd) <- paste0("mu.", colnames(xd))
  colnames(zd) <- paste0("s2.", colnames(zd))
  coefnames <- c(colnames(xd), colnames(zd))

  k <- ncol(xd)
  p <- ncol(zd)
  nobs <- nrow(xd)


  start.list <- mvreg_start(y, xd, zd, start.s2 = start.s2)

  start <- start.list$start
  names(start) <- coefnames

  b0 <- start[1:k]
  t0 <- start[(k + 1):length(start)]

  fit.list <- mvreg_fit(y, xd, zd, b0, t0, tol = tol, maxit = maxit, method = method, vcov.type = vcov.type)

  it <- fit.list$it

  theta0 <- fit.list$theta
  names(theta0) <- coefnames
  b0 <- theta0[1L:k]
  t0 <- theta0[(k + 1L):length(theta0)]
  vtheta <- fit.list$vtheta

  colnames(vtheta) <- rownames(vtheta) <- names(theta0)

  fit.mu <- as.vector(xd %*% b0)
  fit.log.s2 <- as.vector(zd %*% t0)
  fit.s2 <- exp(fit.log.s2)
  residuals <- as.vector(y - fit.mu)


  results <- list(
    coefficients = theta0,
    coefficients.mu = theta0[1L:k],
    coefficients.s2 = theta0[(k + 1L):length(theta0)],
    vcov = vtheta,
    vcov.mu = vtheta[1L:k, 1L:k],
    vcov.s2 = vtheta[(k + 1L):nrow(vtheta), (k + 1L):ncol(vtheta)],
    logLik = mvreg_loglik(y, xd, zd, b0, t0),
    fit.mu = fit.mu,
    fit.log.s2 = fit.log.s2,
    fit.s2 = fit.s2,
    it = it,
    start = start,
    y = y,
    xd = xd,
    zd = zd,
    nobs = nobs,
    call = cl,
    df.residual = length(y) - (p + k),
    residuals = residuals,
    response = response,
    colx = colx,
    colz = colz,
    formula.mu = formula.mu,
    formula.s2 = formula.s2,
    terms.mu = terms.mu,
    terms.s2 = terms.s2
  )

  class(results) <- "mvreg"

  results
}
