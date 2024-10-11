#-------------------------------------------------------------------------------


#' Print method for mvreg
#'
#' @param x A mvreg object.
#' @param digits Minimal number of significant digits, see print.default.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Print informations about the model.
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


#' Variance-Covariance matrix for mvreg
#'
#' @param object A mvreg object.
#' @param partition A character vector that specifies which partition of the matrix to return.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Specified partition of Variance-Covariance matrix of model's parameters.
#' @export
#'
#' @examples
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, ~Species, data = iris)
#' summary(mvreg_mod)
#' vcov(mvreg_mod)
#' vcov(mvreg_mod, partition = "mu")
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


#' Extract mvreg coefficients
#'
#' @param object A mvreg object.
#' @param partition A character vector that specifies which partition of coefficients to return.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Specified coefficients of mvreg model.
#' @export
#'
#' @examples
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, ~Species, data = iris)
#' summary(mvreg_mod)
#' coef(mvreg_mod)
#' coef(mvreg_mod, partition = "s2")
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


#' Summary method for mvreg
#'
#' @param object A mvreg object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A summary.mvreg object.
#' @export
#'
#' @importFrom stats vcov pnorm
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
    vcov.s2 = vcov(object, "s2")
  )

  class(res) <- "summary.mvreg"
  res
}


#-------------------------------------------------------------------------------


#' Print method for summary.mvreg
#'
#' @param x A mvreg object.
#' @param digits Minimal number of significant digits, see print.default.
#' @param signif.stars logical. If TRUE, ‘significance stars’ are printed for each coefficient.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Print summary for mvreg.
#' @export
#'
#' @importFrom stats printCoefmat
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

  cat("\n\n")
}


#-------------------------------------------------------------------------------


#' Fitted values of mvreg model
#'
#' @param x A mvreg object.
#' @param type A character vector specifying the component for which to return fitted values.
#'
#' @return Fitted values for a given component of a mvreg model.
#' @export
#'
#' @examples
#' mvreg_mod <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris) # different formulas
#' # fitted
#' fitted(mvreg_mod)
#' fitted(mvreg_mod, type = "mu")
#' fitted(mvreg_mod, type = "log.s2")
#' fitted(mvreg_mod, type = "s2")
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


#' logLik for mvreg
#'
#' @param x A mvreg object.
#'
#' @return Loglikelihood value of a fitted mvreg model.
#' @export
#'
#' @examples
#' mod <- mvreg(Sepal.Length ~ Species, data = iris)
#' logLik(mod)
logLik.mvreg <- function(x) {
  val <- x$logLik
  df <- ncol(x$x) + ncol(x$z)
  attr(val, "df") <- df
  attr(val, "nobs") <- x$nobs
  class(val) <- "logLik"
  val
}


#-------------------------------------------------------------------------------


#' Predict method for mvreg
#'
#' @param object A mvreg object.
#' @param type A character vector specifying the component for which to return predicted values.
#' @param newdata An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param se.fit A logical value indicating if standard errors are to be returned.
#' @param interval A logical value indicating if confidence intervals are to be returned.
#' @param sig.level Confidence level for confidence intervals.
#'
#' @return Predicted values for the different components of mvreg model. If specified, standard error estimates and confidence interval are returned.
#' @export
#'
#' @importFrom stats model.matrix update qnorm
#'
#' @examples
#' mvreg_mod1 <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris) # different formulas
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
predict.mvreg <- function(object, type = c("all", "mu", "log.s2", "s2"), newdata, se.fit = F, interval = F, sig.level = 0.95) {
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
    x <- model.matrix(update(object$formula.mu, NULL ~ .), newdata)
    z <- model.matrix(update(object$formula.s2, NULL ~ .), newdata)
  } else {
    stop("newdata must be a data.frame whose column names must be the
         same as the names of the variables in the model")
  }

  b <- object$coefficients.mu
  t <- object$coefficients.s2

  pred.mu <- as.vector(x %*% b)
  pred.log.s2 <- as.vector(z %*% t)
  pred.s2 <- exp(pred.log.s2)

  if (se.fit == F & interval == F) {
    return(list(mu = pred.mu, log.s2 = pred.log.s2, s2 = pred.s2)[select])
  } else if (se.fit == T & interval == F) {
    vcov.mu <- object$vcov.mu
    vcov.log.s2 <- object$vcov.s2

    se.mu <- sqrt(rowSums((x %*% vcov.mu) * x))
    se.log.s2 <- sqrt(rowSums((z %*% vcov.log.s2) * z))
    se.s2 <- sqrt(exp(2 * pred.log.s2) * se.log.s2) # delta method

    dmu <- data.frame(pred.mu, se.mu)
    dlogs2 <- data.frame(pred.log.s2, se.log.s2)
    ds2 <- data.frame(pred.s2, se.s2)

    return(list(mu = dmu, log.s2 = dlogs2, s2 = ds2)[select])
  } else if (se.fit == T & interval == T) {
    vcov.mu <- object$vcov.mu
    vcov.log.s2 <- object$vcov.s2

    se.mu <- sqrt(rowSums((x %*% vcov.mu) * x))
    se.log.s2 <- sqrt(rowSums((z %*% vcov.log.s2) * z))
    se.s2 <- sqrt(exp(2 * pred.log.s2) * se.log.s2) # delta method

    quant <- qnorm(1 - (1 - sig.level) / 2)
    confint.mu <- cbind(pred.mu - quant * se.mu, pred.mu + quant * se.mu)
    confint.log.s2 <- cbind(pred.log.s2 - quant * se.log.s2, pred.log.s2 + quant * se.log.s2)
    confint.s2 <- exp(confint.log.s2)

    colnames(confint.mu) <- colnames(confint.log.s2) <- colnames(confint.s2) <- paste0(c("lwr", "upr"), sig.level)

    dmu <- data.frame(pred.mu, se.mu, confint.mu)
    dlogs2 <- data.frame(pred.log.s2, se.log.s2, confint.log.s2)
    ds2 <- data.frame(pred.s2, se.s2, confint.s2)

    return(list(mu = dmu, log.s2 = dlogs2, s2 = ds2)[select])
  } else if (se.fit == F & interval == T) {
    vcov.mu <- object$vcov.mu
    vcov.log.s2 <- object$vcov.s2

    se.mu <- sqrt(rowSums((x %*% vcov.mu) * x))
    se.log.s2 <- sqrt(rowSums((z %*% vcov.log.s2) * z))
    se.s2 <- sqrt(exp(2 * pred.log.s2) * se.log.s2) # delta method

    quant <- qnorm(1 - (1 - sig.level) / 2)
    confint.mu <- cbind(pred.mu - quant * se.mu, pred.mu + quant * se.mu)
    confint.log.s2 <- cbind(pred.log.s2 - quant * se.log.s2, pred.log.s2 + quant * se.log.s2)
    confint.s2 <- exp(confint.log.s2)
    colnames(confint.mu) <- colnames(confint.log.s2) <- colnames(confint.s2) <- paste0(c("lwr", "upr"), sig.level)

    dmu <- data.frame(pred.mu, confint.mu)
    dlogs2 <- data.frame(pred.log.s2, confint.log.s2)
    ds2 <- data.frame(pred.s2, confint.s2)

    return(list(mu = dmu, log.s2 = dlogs2, s2 = ds2)[select])
  }
}



#-------------------------------------------------------------------------------



#' Simulate Responses from a mvreg object
#'
#' @param object A mvreg object.
#' @param nsim Number of response vectors to simulate
#' @param seed NULL or an integer that will be used in a call to `set.seed` before simulating the response vectors. If set, the value is saved as the "seed" attribute of the returned value. The default, NULL will not change the random generator state, and return .Random.seed as the "seed" attribute.
#' @param ... Additional optional arguments.
#'
#' @return A data.frame which columns are simulated response vectors.
#' @export
#'
#' @importFrom stats rnorm runif simulate
#'
#' @examples
#' mvreg.mod <- mvreg(Sepal.Length ~ Species, data = iris)
#' simulate(mvreg.mod, nsim = 100, seed = 43)
#'
simulate.mvreg <- function(object, nsim = 1, seed = NULL, ...){

  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
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

  mu <- as.vector(x%*%b)
  s <- as.vector(sqrt(exp(z%*%t)))


  sim <- do.call(cbind,
                 lapply(seq_len(nsim),
                        function(i){
                          y <- rnorm(n, mu, s)
                        }))
  colnames(sim) <- paste0("sim_", seq_len(nsim))

  sim <- as.data.frame(sim)

  attr(sim, "seed") <- RNGstate

  sim

}

