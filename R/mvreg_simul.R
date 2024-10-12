
#' Simulation study for mvreg models
#'
#' @param x Design matrix for mean component.
#' @param z Design matrix for variance component.
#' @param b True values for parameters of mean component.
#' @param t True values for parameters of variance component.
#' @param n Sample size.
#' @param nsim Number of samples to simulate.
#' @param sig.level Significance level for Wald tests.
#' @param seed Seed of the simulation.
#' @param method Method chosen for estimation of parameters of mean component.
#' @param vcov.type A string to specify whether to use observed or expected Fisher information matrix in order to compute variance-covariance matrix of estimates.
#' @param start.s2 A character vector indicating how to select initial values for variance component parameters.
#' @param tol Positive value indicating what is the minimum difference between parameter estimates between two iterations to stop the algorithm.
#' @param maxit Integer value indicating the maximum number of iteration.
#'
#' @return Results of a simulation study to check properties of estimators and statistical tests.
#' @export
#' @importFrom stats rnorm pnorm
#'
#' @examples
#' n <- 100
#' x <- cbind(1, rnorm(n))
#' z <- cbind(1, x[,2], x[,2]^2)
#' b <- rnorm(ncol(x))
#' t <- rnorm(ncol(z))
#' mvreg_simul(x, z, b, t, n = n, nsim = 100, seed = 43)
#'
#'
mvreg_simul <- function(x, z, b, t,
                        n = 100, nsim = 100, sig.level = 0.05,
                        seed = NULL,
                        method = c("wls", "full_nr"),
                        vcov.type = c("expected", "observed"),
                        start.s2 = c("residuals", "gamma", "zero"),
                        tol = 1e-10,
                        maxit = 100) {

  method <- match.arg(method)
  vcov.type <- match.arg(vcov.type)
  start.s2 <- match.arg(start.s2)
  k <- ncol(as.matrix(x))
  p <- ncol(as.matrix(z))

  if(is.null(colnames(x))){
    if(k == 1) {colnames(x)[1] <- "mu.const"}
    else if(k > 1) {
      colnames(x) <- c("mu.const", paste0("mu.x", 1:(k-1)))
    }
  } else {
    if(colnames(x)[1] %in% c("", "(Intercept)")){
      colnames(x)[1] <- "mu.const"
    }
  }

  if(is.null(colnames(z))){
    if(p == 1) {colnames(z)[1] <- "s2.const"}
    else if(k > 1) {
      colnames(z) <- c("s2.const", paste0("s2.z", 1:(p-1)))
    }
  } else {
    if(colnames(z)[1] %in% c("", "(Intercept)")){
      colnames(z)[1] <- "s2.const"
    }
  }


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

  mu <- as.vector(x%*%b)
  s <- as.vector(sqrt(exp(z%*%t)))

  y <- lapply(seq_len(nsim),
              function(i){
                rnorm(n, mean = mu, sd = s)
              })

  names(y) <- paste0("sim_", seq_len(nsim))

  start <- lapply(y, function(y){
    mvreg_start(y, x, z, start.s2 = start.s2)$start
  })

  sim_res <- mapply(function(y, start){
    mvreg_fit(y, x, z, start[1:k], start[(k+1):length(start)],
              tol = tol, maxit = maxit,
              method = method, vcov.type = vcov.type)
  }, y, start, SIMPLIFY  = F)

  sim_theta <- do.call(rbind, lapply(sim_res, function(res) res$theta))

  sim_vtheta <- lapply(sim_res, function(res) res$vtheta)

  sim_se <- do.call(rbind, lapply(sim_vtheta, function(res) sqrt(diag(res))))
  sim_t <- sim_theta/sim_se
  sim_p <- 2*(1 - pnorm(abs(sim_t)))
  sim_it <- unlist(lapply(sim_res, function(res) res$it))


  true_value <- c(b, t)
  prop_rej_h0 <- colSums(sim_p < sig.level)/nsim
  distortion <- colMeans(sim_theta) - true_value
  variance <- apply(sim_theta, 2, var)
  SE <- sqrt(variance)
  MSE <- distortion^2 + variance
  RMSE <- sqrt(MSE)


  tab <- data.frame(true_value = true_value,
                    prop_rej_h0 = prop_rej_h0,
                    distortion = distortion,
                    variance = variance,
                    SE = SE,
                    MSE = MSE,
                    RMSE = RMSE)

  rownames(tab) <- c(colnames(x), colnames(z))

  res <- list(tab = tab,
              theta = sim_theta,
              vtheta = sim_vtheta,
              it = sim_it,
              y = data.frame(y),
              n = n,
              nsim = nsim,
              seed = RNGstate,
              k = ncol(x),
              p = ncol(z))
  class(res) = "simul_mvreg"

  res

}



#' Print method for simul_mvreg
#'
#' @param x A simul_mvreg object.
#' @param digits The minimum number of significant digits to be used.
#' @param ... Further arguments passed to or from other methods
#'
#' @return Prints results of simul_mvreg
#' @export
#' @exportS3Method print simul_mvreg
#'
#' @examples
#' n <- 100
#' x <- cbind(1, rnorm(n))
#' z <- cbind(1, x[,2], x[,2]^2)
#' b <- rnorm(ncol(x))
#' t <- rnorm(ncol(z))
#' mvreg_simul(x, z, b, t, n = n, nsim = 100, seed = 43)
print.simul_mvreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  cat("\nSimulation study for a mvreg model\n")
  cat(paste0("\nn = ", x$n))
  cat(paste0("\nnumber of simulations = ", x$nsim, "\n"))
  cat(paste0("\nnumber of parameters for mean component = ", x$k))
  cat(paste0("\nnumber of parameters for variance component = ", x$p, "\n"))
  cat("\n")

  print(x$tab, digits = digits)

  }


