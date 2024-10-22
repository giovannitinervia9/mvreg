#-------------------------------------------------------------------------------


#' ANOVA-like Comparison for mvreg Models
#'
#' This function performs an ANOVA-like analysis for comparing one or more `mvreg` models. It uses likelihood ratio tests (LRT)
#' to compare nested models, providing test statistics and p-values. If only one `mvreg` model is provided, the function
#' compares the model with reduced versions obtained by sequentially removing terms from the
#' mean and variance formulas.
#'
#' @param object A `mvreg` object, representing the main model to compare.
#'        If this is the only model provided, the function compares this model with reduced versions of itself.
#' @param order.models A logical value indicating if the models should be ordered.
#'        If `TRUE` (default), models are ranked by the number of parameters in increasing order,
#'        with ties broken by the value of `-2logLik` in decreasing order.
#' @param ... Additional `mvreg` objects to compare with the first model.
#'
#' @return A list of class `anova.mvreg`:
#'   - If only one model is provided:
#'     - `mu.tests`: A data frame with results of the ANOVA for the mean component, including:
#'       - `-2logLik`: Negative two times the log-likelihood of each model.
#'       - `n.param`: Number of parameters in each model.
#'       - `AIC`: Akaike Information Criterion for each model.
#'       - `LRT`: Likelihood ratio test statistic.
#'       - `df`: Degrees of freedom for the test.
#'       - `Pr(>Chi)`: P-value associated with the likelihood ratio test.
#'     - `s2.tests`: A similar data frame for the variance component.
#'     - `new.formula.mu`: A character vector of formulas used for the mean component.
#'     - `new.formula.s2`: A character vector of formulas used for the variance component.
#'     - `n_models`: Number of models compared.
#'   - If more than one model is provided:
#'     - `tests`: A data frame with results of the likelihood ratio tests, including:
#'       - `-2logLik`: Negative two times the log-likelihood of each model.
#'       - `n.param`: Number of parameters in each model.
#'       - `AIC`: Akaike Information Criterion for each model.
#'       - `LRT`: Likelihood ratio test statistic.
#'       - `Pr(>Chi)`: P-value associated with the likelihood ratio test.
#'     - `formulas.mu`: A character vector of formulas used for the mean component.
#'     - `formulas.s2`: A character vector of formulas used for the variance component.
#'     - `n_models`: Number of models compared.
#'
#' @export
#'
#' @importFrom stats df.residual logLik pchisq na.omit
#'
#' @details
#' The `anova.mvreg` function compares `mvreg` models using likelihood ratio tests
#' (LRT). If multiple models are provided, they are compared based on their likelihoods
#' and parameter counts. If only one model is provided, the function automatically generates
#' reduced models by sequentially dropping terms from the mean and variance formulas.
#' This allows for comparison of nested models and the identification of significant terms.
#'
#' The function checks whether the models being compared are nested using the `are_models_nested()` function.
#' If non-nested models are provided, a warning is issued, as the likelihood ratio test is only valid for nested models.
#' Additionally, models with the same number of parameters cannot be meaningfully compared using the LRT.
#'
#' If `order.models = TRUE`, models are ranked by the number of parameters (increasing order) and further sorted by
#' the `-2logLik` value (decreasing order) if the number of parameters is equal. This ensures that models with the
#' same number of parameters are presented in order of likelihood.
#'
#' The returned object contains detailed information about the comparison, including AIC values,
#' degrees of freedom, likelihood ratio statistics, and their associated p-values.
#'
#' @note If two models have the same number of parameters and/or are nested, the likelihood ratio test is meaningless,
#' and the function will issue a warning. Negative LRT or degrees of freedom values are also flagged as errors.
#'
#'
#' @examples
#' # Fit two models
#' mod1 <- mvreg(Sepal.Length ~ Species + Sepal.Width, data = iris)
#' mod2 <- mvreg(Sepal.Length ~ Species, ~1, data = iris)
#' anova(mod1, mod2)
#'
#' # Compare a single model with its reduced forms
#' mod3 <- mvreg(Sepal.Length ~ Species, ~ Sepal.Width, data = iris)
#' anova(mod3)
anova.mvreg <- function(object, ..., order.models = T){

  models <- list(object, ...)

  if (!all(sapply(models, inherits, "mvreg"))) {
    stop("All objects must be of class 'mvreg'.")
  }

  n_models <- length(models)


  if(n_models == 1) {

    model <- models[[1]]
    n <- model$nobs

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
        NA,
        df.residual(mod))})

    mu.tests <- as.data.frame(do.call(rbind, mu.tests))
    colnames(mu.tests) <- c("-2logLik", "n.param", "AIC", "df.residuals")
    rownames(mu.tests) <- model.mu
    mu.tests$AIC <- mu.tests$`-2logLik` + 2*mu.tests$n.param
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
        NA,
        df.residual(mod))})

    s2.tests <- as.data.frame(do.call(rbind, s2.tests))
    colnames(s2.tests) <- c("-2logLik", "n.param", "AIC", "df.residuals")
    rownames(s2.tests) <- model.s2
    s2.tests$AIC <- s2.tests$`-2logLik` + 2*s2.tests$n.param
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
                n_models = n_models,
                response = model$response)


  } else {

    response <- unique(unlist(lapply(models, function(models) models$response)))
    if(length(response) > 1){
      stop(paste0("Models have different response variables: ",
                  paste(response, collapse = ", "),
                  ". Please provide models fitted with the same response variable."))
    }

    n <- unlist(lapply(models, function(models) models$nobs))
    if(length(unique(n)) > 1){
      stop(paste0("Models are fitted with different numbers of observations: ",
                  paste(unique(n), collapse = ", "),
                  ". Please provide models fitted with the same number of observations."))
    }



    tests <- lapply(models, function(mod) {
      c(-2*as.numeric(logLik(mod)),
        n - df.residual(mod),
        NA,
        df.residual(mod))
    })

    tests <- as.data.frame(do.call(rbind, tests))
    colnames(tests) <- c("-2logLik", "n.param", "AIC", "df.residuals")
    tests$AIC <- tests$`-2logLik` + 2*tests$n.param
    tests$LRT <- NA
    tests$df <- NA
    tests$`Pr(>Chi)` <- NA

    if (order.models == T) {
      reorder.ind <- order(tests$n.param, -tests$`-2logLik`, decreasing = F)
    } else {
      reorder.ind <- 1:length(models)
    }

    tests <- tests[reorder.ind,]
    models <- models[reorder.ind]
    models.names <- paste0("model", 1:length(models))
    rownames(tests) <- models.names

    n.param.equal <- rep(NA, length(2:nrow(tests)))
    model.nested <- n.param.equal

    for(i in 2:nrow(tests)){
      lrt <- tests$`-2logLik`[i - 1] - tests$`-2logLik`[i]
      df <- tests$n.param[i] - tests$n.param[i - 1]
      model.nested[i - 1] <- are_models_nested(models[[i]], models[[i - 1]])

      if(df > 0 & model.nested[i - 1]){
        pval <- 1 - pchisq(lrt, df = df)
      } else {pval <- NA}

      tests$LRT[i] <- lrt
      tests$df[i] <- df
      tests$`Pr(>Chi)`[i] <- pval
      n.param.equal[i - 1] <- tests$n.param[i] == tests$n.param[i - 1]
    }

    if (any(n.param.equal == TRUE)) {
      warning("LRT test with models with the same number of parameters is meaningless")
    }

    if (any(model.nested == FALSE)) {
      warning("LRT test with non nested models is meaningless")
    }

    if (any(na.omit(tests$LRT) < 0) | any(na.omit(tests$df) < 0)) {
      warning("Negatives LRT tests and degrees of freedom are not possible, check your models to see if they are correctly ordered and nested")
    }

    formulas.mu <- unlist(lapply(models, function(mod) deparse(mod$formula.mu)))
    formulas.s2 <- unlist(lapply(models, function(mod) deparse(mod$formula.s2)))

    res <- list(tests = tests,
                formulas.mu = formulas.mu,
                formulas.s2 = formulas.s2,
                n_models = n_models,
                response = response)

  }
  class(res) <- "anova.mvreg"
  res

}


#------------------------------------------------------------------------------


#' Print Method for ANOVA Results of mvreg Models
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
#'
print.anova.mvreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  if(x$n_models == 1) {

    model.mu <- paste0("model", 1:nrow(x$mu.tests))
    model.s2 <- paste0("model", 1:nrow(x$s2.tests))
    new.formula.mu <- sub(".*?(~)", "\\1", x$new.formula.mu)
    new.formula.s2 <- sub(".*?(~)", "\\1", x$new.formula.s2)

    print.mu <- paste0(model.mu, " = ", as.character(new.formula.mu))
    print.s2 <- paste0(model.s2, " = ", as.character(new.formula.s2))
    cat("\n")
    cat(paste0("Response: ", x$response))
    cat("\n\n")

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
    formulas.mu <- sub(".*?(~)", "\\1", x$formulas.mu)
    formulas.s2 <- sub(".*?(~)", "\\1", x$formulas.s2)

    cat("\n")
    cat("Model comparison: mean and variance components\n\n")
    cat(paste0("Response: ", x$response))
    cat("\n\n")
    cat("Mean component formulas:\n")
    cat(paste0("model", 1:nrow(x$tests), ": ", formulas.mu, sep = "\n", collapse = ""))
    cat("\n")
    cat("Variance component formulas:\n")
    cat(paste0("model", 1:nrow(x$tests), ": ", formulas.s2, sep = "\n", collapse = ""))
    cat("\n")
    printCoefmat(x$tests, signif.legend = F, na.print = "",
                 cs.ind = 1, tst.ind = 5)
    cat("\n")

  }


}


