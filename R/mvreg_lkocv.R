#-------------------------------------------------------------------------------


#' Leave-k-out Cross-validation for mvreg models
#'
#' @description
#' Performs leave-k-out cross-validation for models fitted by `mvreg()`. It supports both exact computation of all possible k-fold
#' combinations and approximation using random sampling when the number of possible combinations is large.
#'
#' @param object A mvreg object.
#' @param k Number of observations to leave out in each fold. Must be a positive integer
#'   and less than the total number of observations. Default is 1.
#' @param num_samples Number of random combinations to sample when using
#'   approximation. Only used when the total number of possible combinations exceeds
#'   this value. Default is 1000.
#' @param seed Integer or NULL. Random seed for reproducibility. If NULL, the current
#'   random seed state is used. Default is NULL.
#' @param mc.cores Number of cores to use for parallel processing. Default is
#'   half of the available cores.
#'
#' @return The average mean squared prediction error divided by k, indicating the model's
#' predictive performance, where lower values signify better predictions.
#'
#' @details
#' For each leave-k-out sample, the function:
#' \enumerate{
#'   \item Refits the `mvreg` model excluding k observations
#'   \item Uses the fitted model to predict the responses for the left-out observations
#'   \item Calculates the prediction error
#' }
#'
#' When the number of possible combinations (`choose(n,k)`) exceeds `num_samples`,
#' the function uses random sampling to approximate the leave-k-out cross-validation.
#'
#' Parallel processing is implemented using either \code{mclapply} (Unix-like systems)
#' or \code{parLapply} (Windows) to improve computational efficiency.
#'
#' @examples
#' \dontrun{
#' # Fit a model with mvreg
#' fit <- mvreg(Sepal.Length ~ Species, data = iris)
#'
#' # Perform leave-1-out cross-validation
#' mvreg_lkocv(fit)
#'
#' # Perform leave-2-out cross-validation with approximation
#' mvreg_lkocv(fit, k = 2)
#'
#' # Perform leave-2-out cross-validation without approximation
#' mvreg_lkocv(object, k = 2, num_samples = choose(object$nobs, 2))
#'
#' # Use approximation with 500 samples
#' mvreg_lkocv(fit, k = 2, num_samples = 500, seed = 123)
#' }
#'
#'
#' @importFrom parallel mclapply detectCores makeCluster stopCluster clusterExport parLapply
#' @importFrom utils combn
#'
#' @export
mvreg_lkocv <- function(object, k = 1, num_samples = 1000, seed = NULL, mc.cores = round(parallel::detectCores() / 2)) {
  # Input validation
  if (!inherits(object, "mvreg")) {
    stop("object must be a mvreg model")
  }
  n <- object$nobs
  if (!is.numeric(k) || k <= 0 || k >= n) {
    stop("k must be a positive integer less than the number of observations")
  }
  if (!is.numeric(num_samples) || num_samples <= 0) {
    stop("num_samples must be a positive integer")
  }

  # Manage random seed state
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    runif(1) # Initialize RNG state if not set
  }
  if (!is.null(seed)) {
    old_seed <- .Random.seed
    set.seed(seed)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }


  y <- object$y
  x <- object$xd
  z <- object$zd
  b <- object$coefficients.mu
  t <- object$coefficients.s2


  # Check if approximation is needed
  max_combinations <- choose(n, k)
  approximate <- num_samples < max_combinations
  if (approximate) {
    warning("Using approximation due to large number of combinations")
  }

  # Generate indices for Leave-k-out cross-validation
  if (approximate) {
    ny <- seq_len(n)
    indices <- vector("list", length = num_samples)
    generated_combinations <- new.env(hash = TRUE, parent = emptyenv())
    it <- 1
    max_attempts <- num_samples * 100  # Prevent infinite loop
    attempts <- 0

    while (it <= num_samples && attempts < max_attempts) {
      # Generate a unique combination of indices
      ind <- sort(sample(ny, k, replace = FALSE))
      ind_key <- paste(ind, collapse = ",")
      if (!exists(ind_key, envir = generated_combinations)) {
        indices[[it]] <- ind
        assign(ind_key, TRUE, envir = generated_combinations)
        it <- it + 1
      }
      attempts <- attempts + 1
    }

    if (it <= num_samples) {
      warning("Could not generate requested number of unique combinations")
    }
  } else {
    indices <- utils::combn(n, k, simplify = FALSE)
  }

  # Function to process each fold
  process_fold <- function(leave_out) {
    tryCatch({
      # Fit model without the k observations in leave_out
      fit <- mvreg_fit(y[-leave_out], x[-leave_out, , drop = FALSE],
                       z[-leave_out, , drop = FALSE], b, t)

      # Calculate prediction errors for the left-out observations
      fit_errors <- sapply(leave_out, function(i) {
        fit_y_i <- as.vector(x[i, , drop = FALSE] %*% fit$b)
        (y[i] - fit_y_i)^2
      })

      sum(fit_errors) # Return sum of squared errors for this combination
    }, error = function(e) {
      warning("Error in fold calculation: ", e$message)
      NA
    })
  }

  # Choose parallel processing method based on OS
  if (.Platform$OS.type == "unix") {
    cv_errors <- parallel::mclapply(indices, process_fold, mc.cores = mc.cores)
  } else {
    # For Windows, use parLapply
    cl <- parallel::makeCluster(mc.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export necessary objects to cluster
    parallel::clusterExport(cl, c("y", "x", "z", "b", "t", "mvreg_fit"),
                  envir = environment())

    cv_errors <- parallel::parLapply(cl, indices, process_fold)
  }

  # Calculate and return the average mean squared error
  cv_errors_clean <- unlist(cv_errors)
  if (all(is.na(cv_errors_clean))) {
    stop("All cross-validation calculations failed")
  }

  mean(cv_errors_clean, na.rm = TRUE) / k
}






