#' Convert a mvreg model to lm
#'
#' This function converts an object of class `mvreg` into a standard linear model (`lm` object).
#' It extracts the formula for the mean component (i.e., `formula.mu`) from the `mvreg` model and fits a linear model
#' to the same data with `lm()` function.
#'
#' @param object A `mvreg` or `summary.mvreg` object.
#' @return A `lm` object corresponding to the mean component of the `mvreg` model.
#'         The output is a fitted linear model using the same formula as `formula.mu` from the original `mvreg` object.
#'
#' @importFrom stats lm
#' @export
#'
#'
#' @seealso \code{\link{mvreg}} for fitting heteroscedastic regression models
#' and \code{\link[stats]{lm}} for standard linear regression.
#'
#' @note This function only captures the mean component of the `mvreg` model.
#' The variance component (which models heteroscedasticity) is not included in the resulting `lm` object.
#'
#' @details The \code{mvreg_to_lm} function extracts the formula for the mean component from the `mvreg` model, which
#'         is stored in `formula.mu`. It then fits an ordinary least squares regression model using this formula.
#'         If the `mvreg` model was fit with a dataset, the same dataset is used for fitting the `lm` model.
#'         This enables a direct comparison between the two models.
#'
#' @param object A fitted model object of class `mvreg` or `summary.mvreg`.
#'
#' @return A fitted model object of class `lm`, which is the linear model fitted to the mean component of the original
#'         `mvreg` model.
#'
#' @examples
#' # Fit a mvreg model to the cars dataset
#' mod.mvreg <- mvreg(dist ~ speed, data = cars)
#'
#' # Convert the mvreg model to an lm object
#' mod.lm <- mvreg_to_lm(mod.mvreg)
mvreg_to_lm <- function(object) {
  if (!(class(object) %in% c("mvreg", "summary.mvreg"))) {
    stop("Object is not of class 'mvreg' or 'summary.mvreg'.")
  }

  formula.mu <- object$call$formula.mu


  if (!is.null(object$call$data)) {
    data_name <- as.character(object$call$data)
    data <- eval(object$call$data, envir = parent.frame())

    lm_model <- lm(formula.mu, data = data)
    lm_model$call$data <- as.name(data_name)
  } else {
    lm_model <- lm(formula.mu)
  }

  lm_model
}
