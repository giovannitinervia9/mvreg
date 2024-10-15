#' Convert mvreg model to lm
#'
#' @param object A mvreg object
#'
#' @return A lm object with same formula as the formula.mu of the initial mvreg object
#' @export
#'
#' @importFrom stats lm
#'
#' @examples
#' mod.mvreg <- mvreg(dist ~ speed, data = cars)
#' mvreg_to_lm(mod.mvreg)
mvreg_to_lm <- function(object) {
  if(!(class(object) %in% c("mvreg", "summary.mvreg"))){
    stop("Object is not of class 'mvreg' or 'summary.mvreg'.")
  }
  formula.mu <- object$call$formula.mu
  data <- eval(object$call$data, envir = parent.frame())
  lm_model <- lm(formula.mu, data = data)
  lm_model
}
