#-------------------------------------------------------------------------------


#' ANOVA for mvreg models
#'
#' This function performs an ANOVA-like analysis for `mvreg` objects to compare models
#' with different formulas for the mean and variance components. It calculates the likelihood ratio test
#' statistics and their associated p-values.
#'
#' @param object A `mvreg` object. If only one `mvreg` object is specified, the model is compared with all the reduced models obtained by sequentially dropping one term at a time from the original formula.
#' @param ... Additional `mvreg` objects to be compared with the first model.
#'
#'
#' @return A list of class `anova.mvreg` containing:
#'   - If one model is provided (i.e., `...` is empty):
#'     - `mu.tests`: A data frame with results of the ANOVA for the mean component, including:
#'       - `-2logLik`: The negative two times the log-likelihood of each model.
#'       - `n.param`: The number of parameters in each model.
#'       - `LRT`: The likelihood ratio test statistic.
#'       - `df`: The degrees of freedom for the test.
#'       - `Pr(>Chi)`: The p-value associated with the likelihood ratio test.
#'     - `s2.tests`: A similar data frame for the variance component.
#'     - `new.formula.mu`: A character vector of formulas used for the mean component.
#'     - `new.formula.s2`: A character vector of formulas used for the variance component.
#'     - `n_models`: The number of models compared.
#'
#'   - If more than one model is provided:
#'     - `tests`: A data frame with results of the likelihood ratio tests, including:
#'       - `-2logLik`: The negative two times the log-likelihood of each model.
#'       - `n.param`: The number of parameters in each model.
#'       - `LRT`: The likelihood ratio test statistic.
#'       - `Pr(>Chi)`: The p-value associated with the likelihood ratio test.
#'     - `formulas.mu`: A character vector of formulas used for the mean component.
#'     - `formulas.s2`: A character vector of formulas used for the variance component.
#'     - `n_models`: The number of models compared.
#'
#' @export
#'
#' @importFrom stats df.residual
#'
#' @examples
#' # Fit two models
#' mod1 <- mvreg(Sepal.Length ~ Species + Sepal.Width, data = iris)
#' mod2 <- mvreg(Sepal.Length ~ Species, ~1, data = iris)
#' anova(mod1, mod2)
#'
#' # Fit a single model
#' mod3 <- mvreg(Sepal.Length ~ Species, ~ Sepal.Width, data = iris)
#' anova(mod3)
anova.mvreg <- function(object, ...){

  models <- list(object, ...)

  if (!all(sapply(models, inherits, "mvreg"))) {
    stop("All objects must be of class 'mvreg'.")
  }

  n_models <- length(models)
  n <- models[[1]]$nobs


  if(n_models == 1) {

    model <- models[[1]]

    new.formula.mu <- get_reduced_formulas(model$response, model$formula.mu)
    new.formula.mu <- new.formula.mu[length(new.formula.mu):1]
    new.formula.s2 <- get_reduced_formulas(model$response, model$formula.s2)
    new.formula.s2 <- new.formula.s2[length(new.formula.s2):1]

    new.mod.mu <- lapply(new.formula.mu, function(new.formula.mu) {
      update(model, new.formula.mu = new.formula.mu)
    })
    names(new.mod.mu) <- as.character(new.formula.mu)

    new.mod.s2 <- lapply(new.formula.s2, function(new.formula.s2) {
      update(model, new.formula.s2 = new.formula.s2)
    })
    names(new.mod.s2) <- as.character(new.formula.s2)

    model.mu <- paste0("model", 1:length(new.mod.mu))
    model.s2 <- paste0("model", 1:length(new.mod.s2))

    mu.tests <- lapply(new.mod.mu, function(mod) {
      c(-2*as.numeric(logLik(mod)),
        n - df.residual(mod),
        AIC(mod),
        df.residual(mod))})

    mu.tests <- as.data.frame(do.call(rbind, mu.tests))
    colnames(mu.tests) <- c("-2logLik", "n.param", "AIC", "df.residuals")
    rownames(mu.tests) <- model.mu
    mu.tests$LRT <- NA
    mu.tests$df <- NA
    mu.tests$`Pr(>Chi)` <- NA
    if(nrow(mu.tests) > 1) {
      for(i in 2:nrow(mu.tests)){
        lrt <- mu.tests$`-2logLik`[i - 1] - mu.tests$`-2logLik`[i]
        df <- mu.tests$n.param[i] - mu.tests$n.param[i - 1]
        pval <- 1 - pchisq(lrt, df = df)
        mu.tests$LRT[i] <- lrt
        mu.tests$df[i] <- df
        mu.tests$`Pr(>Chi)`[i] <- pval
      }
    }



    s2.tests <- lapply(new.mod.s2, function(mod) {
      c(-2*as.numeric(logLik(mod)),
        n - df.residual(mod),
        AIC(mod),
        df.residual(mod))})

    s2.tests <- as.data.frame(do.call(rbind, s2.tests))
    colnames(s2.tests) <- c("-2logLik", "n.param", "AIC", "df.residuals")
    rownames(s2.tests) <- model.s2
    s2.tests$LRT <- NA
    s2.tests$df <- NA
    s2.tests$`Pr(>Chi)` <- NA

    if(nrow(s2.tests) > 1) {
      for(i in 2:nrow(s2.tests)){
        lrt <- s2.tests$`-2logLik`[i - 1] - s2.tests$`-2logLik`[i]
        df <-s2.tests$n.param[i] - s2.tests$n.param[i - 1]
        pval <- 1 - pchisq(lrt, df = df)
        s2.tests$LRT[i] <- lrt
        s2.tests$df[i] <- df
        s2.tests$`Pr(>Chi)`[i] <- pval
      }
    }


    res <- list(mu.tests = mu.tests,
                s2.tests = s2.tests,
                new.formula.mu = as.character(new.formula.mu),
                new.formula.s2 = as.character(new.formula.s2),
                n_models = n_models)


  } else {

    tests <- lapply(models, function(mod) {
      c(-2*as.numeric(logLik(mod)),
        n - df.residual(mod),
        AIC(mod),
        df.residual(mod))
      })

    tests <- as.data.frame(do.call(rbind, tests))
    colnames(tests) <- c("-2logLik", "n.param", "AIC", "df.residuals")
    tests$LRT <- NA
    tests$df <- NA
    tests$`Pr(>Chi)` <- NA

    reorder.ind <- order(tests$n.param, decreasing = F)
    tests <- tests[reorder.ind,]
    models <- models[reorder.ind]
    models.names <- paste0("model", 1:length(models))
    rownames(tests) <- models.names

    n.param.equal <- rep(NA, length(2:nrow(tests)))
    for(i in 2:nrow(tests)){
      lrt <- tests$`-2logLik`[i - 1] - tests$`-2logLik`[i]
      df <- tests$n.param[i] - tests$n.param[i - 1]

      if(df > 0){
        pval <- 1 - pchisq(lrt, df = df)
        } else {pval <- NA}

      tests$LRT[i] <- lrt
      tests$df[i] <- df
      tests$`Pr(>Chi)`[i] <- pval
      n.param.equal[i - 1] <- tests$n.param[i] == tests$n.param[i - 1]
    }

    if(any(n.param.equal == T)){warning("LRT test with models with the same number of parameters is meaningless")}

    formulas.mu <- unlist(lapply(models, function(mod) deparse(mod$formula.mu)))
    formulas.s2 <- unlist(lapply(models, function(mod) deparse(mod$formula.s2)))

    res <- list(tests = tests,
                formulas.mu = formulas.mu,
                formulas.s2 = formulas.s2,
                n_models = n_models)

  }
  class(res) <- "anova.mvreg"
  res

}


#------------------------------------------------------------------------------


#' Print Method for ANOVA Results of mvreg models
#'
#' This function formats and prints the results of an ANOVA-like analysis
#' for `mvreg` objects. It provides a summary of the tests conducted on
#' the mean and variance components of the models.
#'
#' @param x An object of class `anova.mvreg`, which contains the
#' results of the ANOVA analysis.
#' @param digits An integer specifying the number of significant digits
#' to print. Default is the maximum of 3 or the `digits` option minus 3.
#' @param ... Additional arguments (not used).
#'
#' @details
#' If only one model was provided to the `anova.mvreg` function, the output
#' includes tests for both the mean and variance components, along with their
#' respective formulas. If multiple models were provided, it summarizes the
#' likelihood ratio tests for the models compared, along with their formulas.
#'
#' @return
#' This function prints the results directly to the console and does not
#' return a value.
#'
#' @export
#'
#' @examples
#' # Fit two models
#' mod1 <- mvreg(Sepal.Length ~ Species + Sepal.Width, data = iris)
#' mod2 <- mvreg(Sepal.Length ~ Species, ~1, data = iris)
#' anova(mod1, mod2)
#'
#' # Fit a single model
#' mod3 <- mvreg(Sepal.Length ~ Species, ~ Sepal.Width, data = iris)
#' anova(mod3)
print.anova.mvreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  if(x$n_models == 1) {

    model.mu <- paste0("model", 1:nrow(x$mu.tests))
    model.s2 <- paste0("model", 1:nrow(x$s2.tests))


    print.mu <- paste0(model.mu, " = ", as.character(x$new.formula.mu))
    print.s2 <- paste0(model.s2, " = ", as.character(x$new.formula.s2))

    cat("Mean Model Comparison (variance component taken as in the full model)\n")
    cat(print.mu, sep = "\n")
    cat("\n")
    printCoefmat(x$mu.tests, signif.legend = F, na.print = "",
                 cs.ind = 1, tst.ind = 3)
    cat("\n")
    cat("---\n\n")

    cat("Variance Model Comparison (mean component taken as in the full model):\n")
    cat(print.s2, sep = "\n")
    cat("\n")
    printCoefmat(x$s2.tests, signif.legend = F, na.print = "",
                 cs.ind = 1, tst.ind = 5)
    cat("\n")

  } else {
    cat("Model comparison: mean and variance components\n\n")
    cat("Mean component formulas:\n")
    cat(paste0("model", 1:nrow(x$tests), ": ", x$formulas.mu, sep = "\n", collapse = ""))
    cat("\n")
    cat("Variance component formulas:\n")
    cat(paste0("model", 1:nrow(x$tests), ": ", x$formulas.s2, sep = "\n", collapse = ""))
    cat("\n")
    printCoefmat(x$tests, signif.legend = F, na.print = "",
                 cs.ind = 1, tst.ind = 5)
    cat("\n")

  }


}


