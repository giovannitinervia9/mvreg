#-------------------------------------------------------------------------------


#' Heteroscedastic linear model with variance as a function of linear predictor.
#'
#' @param formula.mu A two-sided formula describing the model for mean component.
#' @param formula.s2 A one-sided formula describing the model for variance component.
#' @param data An optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param tol Positive value indicating what is the minimum difference between parameter estimates between two iterations to stop the algorithm.
#' @param maxit Integer value indicating the maximum number of iteration.
#' @param start.s2 A character vector indicating how to select initial values for variance component parameters.
#' @param method Method chosen for estimation of parameters of mean component
#'
#' @return An object of class mvreg
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
#' predict(mvreg_mod1, se.fit = TRUE, interval = TRUE, sig.level = 0.99)
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
                  method = c("wls", "full_nr")) {
  cl <- match.call()
  method <- match.arg(method)

  if (is.null(data)) {
    data <- environment()
  }

  response <- all.vars(formula.mu)[attr(terms(formula.mu), "response")]
  colx <- all.vars(formula.mu)[-which(all.vars(formula.mu) == response)]


  start.s2 <- match.arg(start.s2)


  mf.mu <- model.frame(formula.mu, data)
  y <- model.response(mf.mu)
  x <- model.matrix(attr(mf.mu, "terms"), mf.mu)




  if (is.null(formula.s2)) {
    z <- x
    colz <- colx
  } else {
    mf.s2 <- model.frame(formula.s2, data)
    z <- model.matrix(attr(mf.s2, "terms"), mf.s2)
    if (any(all.vars(formula.s2) == response)) {
      colz <- all.vars(formula.s2)[-which(all.vars(formula.s2) == response)]
    } else {
      colz <- all.vars(formula.s2)
    }
  }


  colnames(x)[which(colnames(x) == "(Intercept)")] <- "const"
  colnames(z)[which(colnames(z) == "(Intercept)")] <- "const"

  colnames(x) <- paste0("mu.", colnames(x))
  colnames(z) <- paste0("s2.", colnames(z))

  k <- ncol(x)
  p <- ncol(z)
  nobs <- nrow(x)


  start.list <- mvreg_start(y, x, z, start.s2 = start.s2)

  start <- start.list$start
  names(start) <- c(colnames(x), colnames(z))

  b0 <- start[1:k]
  t0 <- start[(k + 1):length(start)]

  fit.list <- mvreg_fit(y, x, z, b0, t0, tol = tol, maxit = maxit, method = method)

  it <- fit.list$it

  theta0 <- fit.list$theta
  names(theta0) <- c(colnames(x), colnames(z))
  b0 <- theta0[1L:k]
  t0 <- theta0[(k + 1L):length(theta0)]
  vtheta <- fit.list$vtheta

  colnames(vtheta) <- rownames(vtheta) <- names(theta0)

  fit.mu <- as.vector(x %*% b0)
  fit.log.s2 <- as.vector(z %*% t0)
  fit.s2 <- exp(fit.log.s2)
  residuals <- as.vector(y - fit.mu)


  results <- list(
    coefficients = theta0,
    coefficients.mu = theta0[1L:k],
    coefficients.s2 = theta0[(k + 1L):length(theta0)],
    vcov = vtheta,
    vcov.mu = vtheta[1L:k, 1L:k],
    vcov.s2 = vtheta[(k + 1L):nrow(vtheta), (k + 1L):ncol(vtheta)],
    logLik = ll(y, x, z, b0, t0),
    fit.mu = fit.mu,
    fit.log.s2 = fit.log.s2,
    fit.s2 = fit.s2,
    it = it,
    start = start,
    y = y,
    x = x,
    z = z,
    nobs = nobs,
    call = cl,
    df.residual = length(y) - (p + k),
    residuals = residuals,
    response = response,
    colx = colx,
    colz = colz,
    formula.mu = as.formula(formula.mu),
    formula.s2 = as.formula(formula.s2)
  )

  class(results) <- "mvreg"

  results
}
