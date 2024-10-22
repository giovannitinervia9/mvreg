#-------------------------------------------------------------------------------


#' Print method for mvreg objects
#'
#' @param x A `mvreg` object.
#' @param digits Minimal number of significant digits to be printed, see `print.default`.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns the `mvreg` object. Prints the call, coefficients,
#'         and any additional relevant model information to the console.
#' @export
#'
#' @importFrom stats coef
#'
#' @examples
#' mod <- mvreg(Sepal.Length ~ Species, data = iris)
#' print(mod)
print.mvreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits),
      print.gap = 2L,
      quote = FALSE
    )
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(x)
}


#-------------------------------------------------------------------------------


#' Variance-Covariance matrix for mvreg objects
#'
#' @param object A `mvreg` object.
#' @param partition A character vector that specifies which partition of the matrix to return.
#'                  Options are `"all"`, `"mu"` for mean component, or `"s2"` for variance component.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A matrix representing the specified partition of the Variance-Covariance matrix
#'         of the model's parameters. If `partition` is `"all"`, returns the full matrix;
#'         if `"mu"`, returns the matrix for mean parameters;
#'         if `"s2"`, returns the matrix for variance parameters.
#' @export
#'
#' @examples
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, ~Species, data = iris)
#' summary(mvreg_mod)
#' vcov(mvreg_mod)           # Returns the full variance-covariance matrix
#' vcov(mvreg_mod, partition = "mu")  # Returns the variance-covariance matrix for mean parameters
#' vcov(mvreg_mod, partition = "s2")   # Returns the variance-covariance matrix for variance parameters
vcov.mvreg <- function(object, partition = c("all", "mu", "s2"), ...) {
  partition <- match.arg(partition)
  if (partition == "all") {
    r <- as.matrix(object$vcov)
    colnames(r) <- rownames(r) <- names(coef(object))
  } else if (partition == "mu") {
    r <- as.matrix(object$vcov.mu)
    colnames(r) <- rownames(r) <- names(coef(object, "mu"))
  } else {
    r <- as.matrix(object$vcov.s2)
    colnames(r) <- rownames(r) <- names(coef(object, "s2"))
  }
  r
}


#-------------------------------------------------------------------------------


#' Extract coefficients from mvreg objects
#'
#' @param object A `mvreg` object.
#' @param partition A character vector that specifies which partition of coefficients to return.
#'                  Options are `"all"`, `"mu"` for mean component, or `"s2"` for variance component.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A named vector of coefficients from the `mvreg` model. If `partition` is `"all"`,
#'         returns all coefficients; if `"mu"`, returns coefficients for the mean component;
#'         if `"s2"`, returns coefficients for the variance component.
#' @export
#'
#' @examples
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, ~Species, data = iris)
#' summary(mvreg_mod)
#' coef(mvreg_mod)             # Returns all coefficients
#' coef(mvreg_mod, partition = "mu")  # Returns coefficients for mean component
#' coef(mvreg_mod, partition = "s2")  # Returns coefficients for variance component
coef.mvreg <- function(object, partition = c("all", "mu", "s2"), ...) {
  partition <- match.arg(partition)
  if (partition == "all") {
    r <- object$coefficients
  } else if (partition == "mu") {
    r <- object$coefficients.mu
  } else {
    r <- object$coefficients.s2
  }
  r
}


#-------------------------------------------------------------------------------


#' Summary method for mvreg objects
#'
#' @param object A `mvreg` object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A list containing:
#'         \item{call}{the call used to fit the model.}
#'         \item{residuals}{residuals from the fitted model.}
#'         \item{coefficients}{a combined data frame of coefficients for both mean and variance components.}
#'         \item{coefficients.mu}{a data frame of coefficients for the mean component.}
#'         \item{coefficients.s2}{a data frame of coefficients for the variance component.}
#'         \item{df}{degrees of freedom.}
#'         \item{vcov}{variance-covariance matrix of all coefficients.}
#'         \item{vcov.mu}{variance-covariance matrix of mean component coefficients.}
#'         \item{vcov.s2}{variance-covariance matrix of variance component coefficients.}
#'         \item{loglik}{log-likelihood of the fitted model.}
#'         \item{AIC}{Akaike Information Criterion.}
#'         \item{BIC}{Bayesian Information Criterion.}
#' @export
#'
#' @importFrom stats vcov pnorm logLik AIC BIC
#'
#' @examples
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, data = iris) # same formula for mean and variance
#' mvreg_mod1 <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris) # different formulas
#' summary(mvreg_mod1)
#' summary(mvreg_mod)
summary.mvreg <- function(object, ...) {
  object.call <- object$call

  mu.tab <- data.frame(
    estimate = coef(object, "mu"),
    se = sqrt(diag(vcov(object, "mu")))
  )
  mu.tab$z.value <- mu.tab$estimate / mu.tab$se
  mu.tab$`Pr(>|z|)` <- 2 * pmin(
    pnorm(mu.tab$z.value, lower.tail = T),
    pnorm(mu.tab$z.value, lower.tail = F)
  )

  s2.tab <- data.frame(
    estimate = coef(object, "s2"),
    se = sqrt(diag(vcov(object, "s2")))
  )
  s2.tab$z.value <- s2.tab$estimate / s2.tab$se
  s2.tab$`Pr(>|z|)` <- 2 * pmin(
    pnorm(s2.tab$z.value, lower.tail = T),
    pnorm(s2.tab$z.value, lower.tail = F)
  )

  coefficients <- rbind(mu.tab, s2.tab)
  df <- object$df.residual


  res <- list(
    call = object.call,
    residuals = object$residuals,
    coefficients = coefficients,
    coefficients.mu = mu.tab,
    coefficients.s2 = s2.tab,
    df = df,
    vcov = vcov(object),
    vcov.mu = vcov(object, "mu"),
    vcov.s2 = vcov(object, "s2"),
    loglik = logLik(object),
    AIC = AIC(object),
    BIC = BIC(object)
  )

  class(res) <- "summary.mvreg"
  res
}


#-------------------------------------------------------------------------------


#' Print method for summary.mvreg objects
#'
#' @param x A `summary.mvreg` object.
#' @param digits Minimal number of significant digits, see `print.default`.
#' @param signif.stars logical. If TRUE, ‘significance stars’ are printed for each coefficient.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A printed summary of the `mvreg` model, including:
#' \item{Call}{The call used to fit the model.}
#' \item{Residuals}{Summary statistics of the residuals from the fitted model.}
#' \item{Coefficients}{Coefficients for the mean and log(variance) components, including standard errors and significance.}
#' \item{Model Fit Statistics}{Log-likelihood, AIC, and BIC of the fitted model.}
#' \item{Likelihood Ratio Test (LRT)}{Comparison with a linear model (lm), including likelihood ratio, degrees of freedom, and p-value.}
#' @export
#'
#' @importFrom stats printCoefmat AIC BIC pchisq
#' @examples
#' mod <- mvreg(Sepal.Length ~ Species, data = iris)
#' s <- summary(mod)
#' print(s)
print.summary.mvreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                                signif.stars = getOption("show.signif.stars"), ...) {
  residuals.summary <- summary(x$residuals)
  muMat <- x$coefficients.mu
  s2Mat <- x$coefficients.s2

  # CALL
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
    "\n\n",
    sep = ""
  )


  # RESIDUALS
  rounded.residuals <- round(residuals.summary, digits)
  residuals.value.width <- max(nchar(rounded.residuals))
  residuals.formatted.labels <- sprintf(paste0("%", residuals.value.width, "s"), names(residuals.summary))
  residuals.formatted.values <- sprintf(paste0("%", residuals.value.width, ".", digits, "f"), rounded.residuals)
  cat("\nResiduals:\n",
    paste(residuals.formatted.labels, collapse = "   "), "\n",
    paste(residuals.formatted.values, collapse = "   "), "\n\n",
    sep = ""
  )

  cat("\n")

  print(c(logLik = x$loglik, AIC = x$AIC, BIC = x$BIC))

  cat("\n")


  cat("\nCoefficients for mean component:\n")

  printCoefmat(as.matrix(muMat),
    digits = digits, cs.ind = 1:2, tst.ind = 3,
    signif.stars = signif.stars, signif.legend = F, P.values = NULL
  )

  cat("\n")


  cat("\nCoefficients for log(variance) component:\n")

  printCoefmat(as.matrix(s2Mat),
    digits = digits, cs.ind = 1:2, tst.ind = 3,
    signif.stars = signif.stars, P.values = NULL
  )

  cat("\n")

  cat("\nLRT for comparison with a lm model\n")

  mod_lm <- mvreg_to_lm(x)
  dim_h0 <- length(coef(mod_lm)) + 1
  dim_h1 <- nrow(x$coefficients)
  loglik_lm <- logLik(mod_lm)
  lrt <- as.vector(-2*(loglik_lm - x$loglik))
  df <- dim_h1 - dim_h0

  pval <- 1 - pchisq(lrt, df)

  dd <- rbind(c(-2*loglik_lm, NA, NA, NA),
              c(-2*x$loglik, df, lrt, pval))


  colnames(dd) = c("-2logLik", "df", "LRT", "Pr(>Chi)")
  rownames(dd) = c("mvreg", "lm")

  printCoefmat(dd,
               digits = digits, cs.ind = 1, tst.ind = 3,
               signif.stars = signif.stars, signif.legend = F, P.values = TRUE,
               has.Pvalue = TRUE, na.print = ""
  )


}


#-------------------------------------------------------------------------------


#' Fitted values of mvreg model
#'
#' @param x A `mvreg` object.
#' @param type A character vector specifying the component for which to return fitted values.
#'             Options include:
#'             \itemize{
#'               \item{"all"}: Returns fitted values for both mean and variance components.
#'               \item{"mu"}: Returns fitted values for the mean component.
#'               \item{"log.s2"}: Returns fitted values for the log of the variance component.
#'               \item{"s2"}: Returns fitted values for the variance component.
#'             }
#'
#' @return A vector or data frame of fitted values for the specified component of the `mvreg` model.
#' @export
#' @exportS3Method confint mvreg
#'
#' @examples
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris) # different formulas
#' # fitted values
#' fitted(mvreg_mod)            # returns all fitted values
#' fitted(mvreg_mod, type = "mu")  # returns fitted values for the mean component
#' fitted(mvreg_mod, type = "log.s2")  # returns fitted values for the log of variance
#' fitted(mvreg_mod, type = "s2")  # returns fitted values for the variance component
fitted.mvreg <- function(x, type = c("all", "mu", "log.s2", "s2")) {
  type <- match.arg(type)
  if (type == "all") {
    fit.mu <- x$fit.mu
    fit.log.s2 <- x$fit.log.s2
    fit.s2 <- x$fit.s2
    r <- data.frame(fit.mu, fit.log.s2, fit.s2)
  }
  if (type == "mu") {
    r <- x$fit.mu
  } else if (type == "log.s2") {
    r <- x$fit.log.s2
  } else if (type == "s2") {
    r <- x$fit.s2
  }
  r
}


#-------------------------------------------------------------------------------


#' logLik method for mvreg
#'
#' @param object A `mvreg` object.
#' @param ... Additional optional arguments (not currently used).
#'
#' @return A log-likelihood value of the fitted `mvreg` model, with attributes:
#'         \item{df}{Degrees of freedom associated with the model.}
#'         \item{nobs}{Number of observations used in the fitting.}
#' @export
#'
#' @examples
#' mod <- mvreg(Sepal.Length ~ Species, data = iris)
#' logLik(mod)  # returns the log-likelihood of the fitted model
logLik.mvreg <- function(object, ...) {
  val <- object$logLik
  df <- ncol(object$x) + ncol(object$z)
  attr(val, "df") <- df
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  val
}


#-------------------------------------------------------------------------------


#' Predict method for mvreg
#'
#' @param object A `mvreg` object.
#' @param type A character vector specifying the component for which to return predicted values. Options include "all", "mu", "log.s2", or "s2".
#' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param se.fit A logical value indicating if standard errors are to be returned.
#' @param interval A logical value indicating if confidence intervals are to be returned.
#' @param level Confidence level for confidence intervals (default is `0.95`).
#'
#' @return A list containing predicted values for the different components of the `mvreg` model. If specified, standard error estimates and confidence intervals are also returned.
#' @export
#'
#' @importFrom stats model.matrix update qnorm formula
#'
#' @examples
#' mvreg_mod1 <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris) # different formulas
#'
#' # Predict without newdata
#' predict(mvreg_mod1)
#' predict(mvreg_mod1, type = "mu")
#' predict(mvreg_mod1, type = "log.s2", se.fit = TRUE)
#' predict(mvreg_mod1, type = "s2", interval = TRUE)
#' predict(mvreg_mod1, se.fit = TRUE, interval = TRUE, level = 0.99)
#'
#' # Predict with newdata
#' newdata <- data.frame(
#'   Species = levels(iris$Species),
#'   Sepal.Width = c(min(iris$Sepal.Width), mean(iris$Sepal.Width), max(iris$Sepal.Width))
#' )
#'
#' predict(mvreg_mod1, newdata = newdata, se.fit = TRUE, interval = TRUE)
predict.mvreg <- function(object, type = c("all", "mu", "log.s2", "s2"), newdata, se.fit = F, interval = F, level = 0.95) {
  type <- match.arg(type)
  coln <- unique(c(object$colx, object$colz))

  if (type == "all") {
    select <- c("mu", "log.s2", "s2")
  } else {
    select <- type
  }

  noData <- (missing(newdata) || is.null(newdata))

  if (noData) {
    x <- object$x
    z <- object$z
  } else if (setequal(coln, colnames(newdata))) {
    formula.mu <- formula(paste("~ ", deparse(object$formula.mu[[3]])))
    formula.s2 <- formula(paste("~ ", deparse(object$formula.s2[[3]])))
    x <- model.matrix(formula.mu, newdata)
    z <- model.matrix(formula.s2, newdata)
  } else {
    stop("newdata must be a data.frame whose column names must be the
         same as the names of the variables in the model")
  }

  b <- object$coefficients.mu
  t <- object$coefficients.s2

  pred.mu <- as.vector(x %*% b)
  pred.log.s2 <- as.vector(z %*% t)
  pred.s2 <- exp(pred.log.s2)

  if (se.fit == FALSE && interval == FALSE) {
    return(list(mu = pred.mu, log.s2 = pred.log.s2, s2 = pred.s2)[select])
  } else if (se.fit == TRUE && interval == FALSE) {
    vcov.mu <- object$vcov.mu
    vcov.log.s2 <- object$vcov.s2

    se.mu <- sqrt(rowSums((x %*% vcov.mu) * x))
    se.log.s2 <- sqrt(rowSums((z %*% vcov.log.s2) * z))
    se.s2 <- sqrt(exp(2 * pred.log.s2) * se.log.s2^2) # delta method

    dmu <- data.frame(pred.mu, se.mu)
    dlogs2 <- data.frame(pred.log.s2, se.log.s2)
    ds2 <- data.frame(pred.s2, se.s2)

    return(list(mu = dmu, log.s2 = dlogs2, s2 = ds2)[select])
  } else if (se.fit == TRUE && interval == TRUE) {
    vcov.mu <- object$vcov.mu
    vcov.log.s2 <- object$vcov.s2

    se.mu <- sqrt(rowSums((x %*% vcov.mu) * x))
    se.log.s2 <- sqrt(rowSums((z %*% vcov.log.s2) * z))
    se.s2 <- sqrt(exp(2 * pred.log.s2) * se.log.s2^2) # delta method

    quant <- qnorm(1 - (1 - level) / 2)
    confint.mu <- cbind(pred.mu - quant * se.mu, pred.mu + quant * se.mu)
    confint.log.s2 <- cbind(pred.log.s2 - quant * se.log.s2, pred.log.s2 + quant * se.log.s2)
    confint.s2 <- exp(confint.log.s2)

    colnames(confint.mu) <- colnames(confint.log.s2) <- colnames(confint.s2) <- paste0(c("lwr", "upr"), level)

    dmu <- data.frame(pred.mu, se.mu, confint.mu)
    dlogs2 <- data.frame(pred.log.s2, se.log.s2, confint.log.s2)
    ds2 <- data.frame(pred.s2, se.s2, confint.s2)

    return(list(mu = dmu, log.s2 = dlogs2, s2 = ds2)[select])
  } else if (se.fit == FALSE && interval == TRUE) {
    vcov.mu <- object$vcov.mu
    vcov.log.s2 <- object$vcov.s2

    se.mu <- sqrt(rowSums((x %*% vcov.mu) * x))
    se.log.s2 <- sqrt(rowSums((z %*% vcov.log.s2) * z))
    se.s2 <- sqrt(exp(2 * pred.log.s2) * se.log.s2^2) # delta method

    quant <- qnorm(1 - (1 - level) / 2)
    confint.mu <- cbind(pred.mu - quant * se.mu, pred.mu + quant * se.mu)
    confint.log.s2 <- cbind(pred.log.s2 - quant * se.log.s2, pred.log.s2 + quant * se.log.s2)
    confint.s2 <- exp(confint.log.s2)
    colnames(confint.mu) <- colnames(confint.log.s2) <- colnames(confint.s2) <- paste0(c("lwr", "upr"), level)

    dmu <- data.frame(pred.mu, confint.mu)
    dlogs2 <- data.frame(pred.log.s2, confint.log.s2)
    ds2 <- data.frame(pred.s2, confint.s2)

    return(list(mu = dmu, log.s2 = dlogs2, s2 = ds2)[select])
  }
}



#-------------------------------------------------------------------------------



#' Simulate Responses from a mvreg object
#'
#' @param object A `mvreg` object.
#' @param nsim Number of response vectors to simulate.
#' @param seed NULL or an integer that will be used in a call to `set.seed` before simulating the response vectors. If set, the value is saved as the "seed" attribute of the returned value. The default, NULL will not change the random generator state, and return .Random.seed as the "seed" attribute.
#' @param ... Additional optional arguments.
#'
#' @details
#' The `simulate` method generates simulated responses based on the fitted `mvreg` model.
#' The responses are drawn from a normal distribution with a mean specified by the fitted values of the model and a standard deviation that is derived from the estimated variance.
#' The number of simulated response vectors can be specified using the `nsim` parameter.
#'
#' If a seed is provided, the random number generator will be set to this seed to ensure reproducibility of the simulation. The original random seed is restored upon completion.
#'
#' @return A `data.frame` where each column corresponds to a simulated response vector.
#' The attribute "seed" contains the random seed used for the simulation.
#' @export
#'
#' @importFrom stats rnorm runif simulate
#'
#' @examples
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, data = iris)
#' simulate(mvreg_mod, nsim = 100, seed = 43)
simulate.mvreg <- function(object, nsim = 1, seed = NULL, ...) {
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    runif(1)
  }
  if (is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }

  x <- object$x
  z <- object$z
  b <- object$coefficients.mu
  t <- object$coefficients.s2
  n <- object$nobs

  mu <- as.vector(x %*% b)
  s <- as.vector(sqrt(exp(z %*% t)))


  sim <- do.call(
    cbind,
    lapply(
      seq_len(nsim),
      function(i) {
        rnorm(n, mu, s)
      }
    )
  )
  colnames(sim) <- paste0("sim_", seq_len(nsim))

  sim <- as.data.frame(sim)

  attr(sim, "seed") <- RNGstate

  sim
}


#-------------------------------------------------------------------------------


#' Update and Re-fit a mvreg model
#'
#' @param object A `mvreg` object.
#' @param new.formula.mu Changes to `formula.mu` – see `update.formula` for details.
#' @param new.formula.s2 Changes to `formula.s2` – see `update.formula` for details.
#' @param ... Additional arguments to the call, or arguments with changed values. Use `name = NULL` to remove the argument name.
#' @param evaluate If TRUE, evaluate the new call; else return the call.
#'
#' @details
#' The `update` method allows for modifying the existing `mvreg` model formulas (`formula.mu` and `formula.s2`)
#' by specifying new formulas. The method can also accept additional arguments to change other components of the
#' model fitting process (e.g., `method`, `vcov.type`, `start.s2`).
#'
#' If `evaluate` is set to TRUE, the updated model is refitted and returned; otherwise, the updated call is returned.
#'
#' @return A new `mvreg` object with the updated formulas or the updated call.
#' @export
#'
#' @importFrom stats update getCall
#'
#' @examples
#' mod <- mvreg(Sepal.Length ~ Species, data = iris)
#'
#' # changing formulas
#' update(mod, new.formula.mu = . ~ . + Sepal.Width)
#' update(mod, new.formula.mu = . ~ . + Sepal.Width, new.formula.s2 = . ~ . + Petal.Length)
#'
#' # changing method, vcov.type, start.s2
#' update(mod, method = "full_nr", vcov.type = "observed", start.s2 = "gamma")
update.mvreg <- function(object, new.formula.mu, new.formula.s2, ..., evaluate = TRUE) {
  if (is.null(call <- getCall(object))) {
    stop("need an object with call component")
  }

  extras <- match.call(expand.dots = FALSE)$...

  if (missing(new.formula.mu)) {
    call$formula.mu <- object$formula.mu
  } else {
    call$formula.mu <- update(object$formula.mu, new.formula.mu)
  }

  if (missing(new.formula.s2)) {
    call$formula.s2 <- object$formula.s2
  } else {
    call$formula.s2 <- update(object$formula.s2, new.formula.s2)
  }

  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }

  if (evaluate) {
    eval(call, parent.frame())
  } else {
    call
  }
}


#-------------------------------------------------------------------------------


#' Confidence intervals for mvreg model parameters
#'
#' This function computes the confidence intervals for the parameters of
#' an `mvreg` model.
#'
#' @param object A `mvreg` object from which to compute the confidence intervals.
#' @param parm A specification of which parameters to compute confidence intervals for.
#'   This can be either a vector of indices (numbers) or a vector of parameter names (character strings).
#'   If missing, all parameters will be included.
#' @param level A numeric value specifying the confidence level for the intervals (default is `0.95` for 95% confidence intervals).
#' @param ... Additional argument(s) for methods.
#'
#' @return A matrix with columns giving the lower and upper confidence limits for each parameter specified.
#'
#'
#' @details
#' The confidence intervals are calculated using the normal approximation
#' based on the estimated coefficients and their standard errors. The
#' intervals are of the form:
#' \deqn{ \hat{\theta} \pm z_{1 - \frac{\alpha}{2}} \cdot \text{SE}(\hat{\theta}) }
#' where \eqn{\hat{\theta}} is the estimated parameter, \eqn{\text{SE}(\hat{\theta})} is the standard error,
#' and \eqn{z_{1 - \frac{\alpha}{2}}} is the critical value from the standard normal distribution
#' corresponding to the specified confidence level.
#'
#' @export
#'
#' @examples
#' mod <- mvreg(Sepal.Length ~ Species, data = iris)
#' confint(mod)
confint.mvreg <- function(object, parm, level = 0.95, ...){
  cf <- coef(object)
  ses <- sqrt(diag(vcov(object)))
  pnames <- names(cf)
  if (missing(parm)) {
    parm <- pnames
  } else if (is.numeric(parm)) {
    parm <- pnames[parm]
  }

  a <- (1 - level)/2
  a <- c(a, 1 - a)
  z <- qnorm(a)

  ci <- apply(as.matrix(z), 1, function(z) cf[parm] + z*ses[parm])

  if(length(parm) == 1){
    ci <- t(ci)
    rownames(ci) <- parm}

  colnames(ci) <- paste0(c("lwr", "upr"), level)
  ci
}
