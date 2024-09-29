

#-------------------------------------------------------------------------------


#' Starting values for mvreg
#'
#' @param y Response variable vector.
#' @param x Design matrix for mean component.
#' @param z Design matrix for variance component.
#' @param start.s2 A character vector indicating how to select initial values for variance component parameters.
#'
#' @return Starting values of estimates of parameters
#' @export
#'
#' @importFrom stats Gamma glm.fit
#' @examples
mvreg_start <- function(y, x, z, start.s2 = c("residuals", "gamma", "zero")){
  start.s2 <- match.arg(start.s2)
  p <- ncol(z)
  k <- ncol(x)

  b0 <- as.vector(solve(crossprod(x)) %*% crossprod(x, y))

  if (start.s2 == "residuals") {
    r <- log((y - x %*% b0) ^ 2)
    t0 <- as.vector(solve(crossprod(z)) %*% crossprod(z, r))
  }
  else if (start.s2 == "gamma"){
    r <- (y - x %*% b0)^2
    t0 <- glm.fit(z, r, family = Gamma(link = "log"))$coefficients
  }
  else if (start.s2 == "zero") {
    t0 <- rep(0, p)
    t0[1] <- log(var(y))
  }

  start <- c(b0, t0)
  names(start) <- c(colnames(x), colnames(z))
  list(start = start, b0 = b0, t0 = t0)
}
