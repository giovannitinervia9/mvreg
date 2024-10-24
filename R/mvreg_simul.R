#-------------------------------------------------------------------------------



#' Simulation Study for Multivariate Regression Models
#'
#' This function performs a simulation study to evaluate the properties of
#' estimators for the parameters of `mvreg` models. It generates
#' simulated response variables based on the specified design matrices for the
#' mean and variance components and true parameter values, and then fits the
#' model to the simulated data.
#'
#' @param x A numeric matrix representing the design matrix for the mean component.
#'          Each column corresponds to a predictor variable. The first column should be
#'          an intercept term.
#' @param z A numeric matrix representing the design matrix for the variance component.
#'          Each column corresponds to a predictor variable for the variance. The first
#'          column should also be an intercept term.
#' @param b A numeric vector of true values for the parameters of the mean component.
#'          Its length must match the number of columns in `x`.
#' @param t A numeric vector of true values for the parameters of the variance component.
#'          Its length must match the number of columns in `z`.
#' @param nsim An integer indicating the number of samples to simulate. Default is `100`.
#' @param sig.level A numeric value representing the significance level for Wald tests.
#'                  Default is `0.05`.
#' @param seed An optional integer to set the seed for random number generation. This ensures
#'              reproducibility of the results. Default is NULL.
#' @param method A character vector indicating the method chosen for estimation of parameters of mean component. See the documentation of \code{\link{mvreg_fit}} for details.
#' @param vcov.type A character vector to specify whether to use observed or expected Fisher information matrix in order to compute variance-covariance matrix of estimates. See the documentation of \code{\link{mvreg_fit}}, \code{\link{mvreg_hessian_mu}}, \code{\link{mvreg_hessian_s2}}, \code{\link{mvreg_hessian_mus2}} and \code{\link{mvreg_hessian}} for details.
#' @param start.s2 A character vector indicating how to select initial values for variance component parameters. See the documentation of \code{\link{mvreg_start}} for details.
#' @param tol A positive numeric value indicating the minimum difference between parameter
#'            estimates between two iterations to stop the algorithm. Default is `1e-10`.
#' @param maxit An integer value specifying the maximum number of iterations for the fitting
#'               algorithm. Default is `100`.
#'
#' @return A list containing results of the simulation study, including:
#' \describe{
#'   \item{tab}{A data frame summarizing the true values, mean estimates, distortion,
#'              variance, standard errors (SE), mean squared error (MSE), root mean squared
#'              error (RMSE), and proportion of rejected null hypotheses.}
#'   \item{theta}{A data frame of parameter estimates from the simulated models.}
#'   \item{mean_vtheta}{A matrix containing the average variance-covariance matrices of the
#'                      estimates across simulations.}
#'   \item{vtheta}{A list of variance-covariance matrices for each simulated model.}
#'   \item{it}{A vector of the number of iterations taken to converge for each simulation.}
#'   \item{y}{A data frame of the simulated response variables.}
#'   \item{n}{The number of observations used in the simulation.}
#'   \item{nsim}{The number of simulations performed.}
#'   \item{seed}{The seed used for random number generation.}
#'   \item{k}{The number of predictors in the mean component.}
#'   \item{p}{The number of predictors in the variance component.}
#'   \item{converged}{The number of simulations that converged successfully.}
#'   \item{non_converged}{The number of simulations that did not converge.}
#'   \item{total_time}{The total time taken to perform the simulations.}
#' }
#'
#' @export
#' @importFrom stats rnorm pnorm
#'
#' @examples
#' n <- 100
#' x <- cbind(1, rnorm(n)) # Design matrix for mean component
#' z <- cbind(1, x[, 2], x[, 2]^2) # Design matrix for variance component
#' b <- rnorm(ncol(x)) # True parameters for mean component
#' t <- rnorm(ncol(z)) # True parameters for variance component
#' mvreg_simul(x, z, b, t, nsim = 100, seed = 43)
mvreg_simul <- function(x, z, b, t, nsim = 100, sig.level = 0.05,
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
  n <- nrow(x)

  if (n != nrow(z)) {
    stop("x and z must have the same number of rows")
  }

  if (is.null(colnames(x))) {
    if (k == 1) {
      colnames(x)[1] <- "mu.const"
    } else if (k > 1) {
      colnames(x) <- c("mu.const", paste0("mu.x", 1:(k - 1)))
    }
  } else {
    if (colnames(x)[1] %in% c("", "(Intercept)")) {
      colnames(x)[1] <- "mu.const"
      if (k > 1) {
        colnames(x)[-1] <- paste0("mu.", colnames(x)[-1])
      }
    }
  }

  if (is.null(colnames(z))) {
    if (p == 1) {
      colnames(z)[1] <- "s2.const"
    } else if (p > 1) {
      colnames(z) <- c("s2.const", paste0("s2.z", 1:(p - 1)))
    }
  } else {
    if (colnames(z)[1] %in% c("", "(Intercept)")) {
      colnames(z)[1] <- "s2.const"
      if (p > 1) {
        colnames(z)[-1] <- paste0("s2.", colnames(z)[-1])
      }
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

  start_time <- proc.time()

  mu <- as.vector(x %*% b)
  s <- as.vector(sqrt(exp(z %*% t)))

  y <- lapply(
    seq_len(nsim),
    function(i) {
      rnorm(n, mean = mu, sd = s)
    }
  )

  names(y) <- paste0("sim_", seq_len(nsim))

  start <- lapply(y, function(y) {
    mvreg_start(y, x, z, start.s2 = start.s2)$start
  })

  sim_res <- mapply(
    function(y, start) {
      tryCatch(
        {
          mvreg_fit(y, x, z, start[1:k], start[(k + 1):length(start)],
            tol = tol, maxit = maxit,
            method = method, vcov.type = vcov.type
          )
        },
        error = function(e) conditionMessage(e)
      )
    },
    y, start,
    SIMPLIFY = F
  )

  ok <- which(unlist(lapply(sim_res, class)) != "character")

  converged <- length(ok)
  non_converged <- nsim - length(ok)

  if (non_converged == nsim) {
    stop("None of the simulations converged")
  }

  sim_res <- sim_res[ok]

  sim_theta <- do.call(rbind, lapply(sim_res, function(res) res$theta))


  sim_vtheta <- lapply(sim_res, function(res) res$vtheta)

  sim_se <- do.call(rbind, lapply(sim_vtheta, function(res) sqrt(diag(res))))
  sim_t <- sim_theta / sim_se
  sim_p <- 2 * (1 - pnorm(abs(sim_t)))
  sim_it <- unlist(lapply(sim_res, function(res) res$it))

  sim_theta <- as.data.frame(sim_theta)
  colnames(sim_theta) <- c(colnames(x), colnames(z))


  true_value <- c(b, t)
  mean_value <- colMeans(sim_theta)
  prop_rej_h0 <- colSums(sim_p < sig.level) / converged
  distortion <- mean_value - true_value
  variance <- apply(sim_theta, 2, var)
  SE <- sqrt(variance)
  MSE <- distortion^2 + variance
  RMSE <- sqrt(MSE)

  mean_vtheta <- matrix(NA, nrow = k + p, ncol = k + p)

  for (i in 1:(k + p)) {
    for (j in 1:(k + p)) {
      mean_vtheta[i, j] <- mean(unlist(
        lapply(sim_vtheta, function(x) x[i, j])
      ))
    }
  }

  colnames(mean_vtheta) <- rownames(mean_vtheta) <- c(colnames(x), colnames(z))


  end_time <- proc.time()
  total_time <- end_time - start_time

  tab <- data.frame(
    true_value = true_value,
    mean_value = mean_value,
    distortion = distortion,
    variance = variance,
    SE = SE,
    MSE = MSE,
    RMSE = RMSE,
    prop_rej_h0 = prop_rej_h0
  )

  rownames(tab) <- c(colnames(x), colnames(z))

  res <- list(
    tab = tab,
    theta = sim_theta,
    mean_vtheta = mean_vtheta,
    vtheta = sim_vtheta,
    it = sim_it,
    y = data.frame(y),
    n = n,
    nsim = nsim,
    seed = RNGstate,
    k = ncol(x),
    p = ncol(z),
    converged = converged,
    non_converged = non_converged,
    total_time = total_time
  )
  class(res) <- "simul_mvreg"

  res
}



#' Print method for simul_mvreg
#'
#' This method prints the results of the simulation study for a `mvreg` model.
#' It includes the number of observations, the number of simulations, the number
#' of converged simulations, and the total time taken for the simulations.
#'
#' @param x A simul_mvreg object, containing the results of the simulation study.
#' @param digits The minimum number of significant digits to be used. Default is
#'               set to max(3L, getOption("digits") - 3L).
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Prints results of simul_mvreg.
#' @export
#' @exportS3Method print simul_mvreg
#'
#' @examples
#' n <- 100
#' x <- cbind(1, rnorm(n))
#' z <- cbind(1, x[, 2], x[, 2]^2)
#' b <- rnorm(ncol(x))
#' t <- rnorm(ncol(z))
#' result <- mvreg_simul(x, z, b, t, nsim = 100, seed = 43)
#' print(result)
print.simul_mvreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nSimulation study for a mvreg model\n")
  cat(paste0("\nn = ", x$n))
  cat(paste0("\nnumber of simulations = ", x$nsim))
  cat(paste0("\nconverged simulations = ", x$converged, " (", round((x$converged / x$nsim) * 100, 3), "%)"))
  cat(paste0("\ntime needed = ", round(x$total_time[3], 8), " seconds\n"))
  cat("\n")

  print(x$tab, digits = digits)
}
