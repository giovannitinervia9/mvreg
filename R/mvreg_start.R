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
#' n <- 100
#' x1 <- rnorm(n)
#' x2 <- factor(sample(letters[1:3], n, TRUE))
#' x <- model.matrix(~ x1 + x2)
#' z1 <- factor(sample(letters[1:3], n, TRUE))
#' z <- model.matrix(~z1)
#'
#' b <- rnorm(ncol(x))
#' t <- rnorm(ncol(z))
#'
#' y <- rnorm(n, mean = x %*% b, sd = sqrt(exp(z %*% t)))
#'
#' mvreg_start(y, x, z, start.s2 = "residuals")
#' mvreg_start(y, x, z, start.s2 = "gamma")
#' mvreg_start(y, x, z, start.s2 = "zero")
mvreg_start <- function(y, x, z, start.s2 = c("residuals", "gamma", "zero")) {
  start.s2 <- match.arg(start.s2)
  k <- ncol(x)
  p <- ncol(z)

  b0 <- as.vector(solve(crossprod(x)) %*% crossprod(x, y))

  if (start.s2 == "residuals") {
    r <- log((y - x %*% b0)^2)
    t0 <- as.vector(solve(crossprod(z)) %*% crossprod(z, r))
  } else if (start.s2 == "gamma") {
    r <- (y - x %*% b0)^2
    t0 <- glm.fit(z, r, family = Gamma(link = "log"))$coefficients
  } else if (start.s2 == "zero") {
    t0 <- rep(0, p)
    t0[1] <- log(var(y))
  }

  start <- c(b0, t0)
  names(start) <- c(colnames(x), colnames(z))
  b0 <- start[1:k]
  t0 <- start[(k + 1):length(start)]
  list(start = start, b0 = b0, t0 = t0)
}
