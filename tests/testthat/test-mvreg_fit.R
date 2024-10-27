# Setup example data for tests
set.seed(123)
n <- 100
x1 <- rnorm(n)
x2 <- factor(sample(letters[1:3], n, TRUE))
x <- model.matrix(~ x1 + x2)

z1 <- factor(sample(letters[1:3], n, TRUE))
z <- model.matrix(~z1)

b <- rnorm(ncol(x))
t <- rnorm(ncol(z))

y <- rnorm(n, mean = x %*% b, sd = sqrt(exp(z %*% t)))

# Define initial values for parameters using mvreg_start
start.list <- mvreg_start(y, x, z, start.s2 = "gamma")
b0 <- start.list$b0
t0 <- start.list$t0

# Test Weighted Least Squares (WLS) method
test_that("mvreg_fit() with WLS method converges and returns correct structure", {
  fit <- mvreg_fit(y, x, z, b0, t0, method = "wls")

  # Check output structure
  expect_type(fit, "list")
  expect_named(fit, c("theta", "b", "t", "vtheta", "it", "maxit", "tol"))

  # Check dimensions of estimates
  expect_length(fit$theta, ncol(x) + ncol(z)) # Full parameter vector
  expect_length(fit$b, ncol(x)) # Mean component parameters
  expect_length(fit$t, ncol(z)) # Variance component parameters

  # Check that the number of iterations is within the limit
  expect_true(fit$it <= 100)

  # Check variance-covariance matrix dimensions
  expect_equal(dim(fit$vtheta), c(ncol(x) + ncol(z), ncol(x) + ncol(z)))
})

# Test Full Newton-Raphson (full_nr) method
test_that("mvreg_fit() with full NR method converges and returns correct structure", {
  fit <- mvreg_fit(y, x, z, b0, t0, method = "full_nr")

  # Check output structure
  expect_type(fit, "list")
  expect_named(fit, c("theta", "b", "t", "vtheta", "it", "maxit", "tol"))

  # Check dimensions of estimates
  expect_length(fit$theta, ncol(x) + ncol(z)) # Full parameter vector
  expect_length(fit$b, ncol(x)) # Mean component parameters
  expect_length(fit$t, ncol(z)) # Variance component parameters

  # Check that the number of iterations is within the limit
  expect_true(fit$it <= 100)

  # Check variance-covariance matrix dimensions
  expect_equal(dim(fit$vtheta), c(ncol(x) + ncol(z), ncol(x) + ncol(z)))
})

# Test expected Fisher information for vcov.type
test_that("mvreg_fit() with vcov.type = 'expected' returns correct vtheta", {
  fit <- mvreg_fit(y, x, z, b0, t0, vcov.type = "expected")

  # Check that vtheta matrix is positive definite
  eigenvalues <- eigen(fit$vtheta)$values
  expect_true(all(eigenvalues > 0))
})

# Test observed Fisher information for vcov.type
test_that("mvreg_fit() with vcov.type = 'observed' returns correct vtheta", {
  fit <- mvreg_fit(y, x, z, b0, t0, vcov.type = "observed")
  vtheta <- -solve(mvreg_hessian(y, x, z, fit$b, fit$t, type = "observed"))

  # Check that vtheta matrix is positive definite
  eigenvalues <- eigen(fit$vtheta)$values
  expect_true(all(eigenvalues > 0))

  expect_equal(fit$vtheta, vtheta)
})

test_that("mvreg_fit() with vcov.type = 'expected' returns correct vtheta", {
  fit <- mvreg_fit(y, x, z, b0, t0, vcov.type = "expected")
  vtheta <- -solve(mvreg_hessian(y, x, z, fit$b, fit$t, type = "expected"))

  # Check that vtheta matrix is positive definite
  eigenvalues <- eigen(fit$vtheta)$values
  expect_true(all(eigenvalues > 0))

  expect_equal(fit$vtheta, vtheta)
})


# Test handling of edge cases with tolerance and max iterations
test_that("mvreg_fit() correctly handles tolerance and maxit parameters", {
  # Test with high tolerance (should converge in fewer iterations)
  fit_tol <- mvreg_fit(y, x, z, b0, t0, tol = 1e-2)
  expect_true(fit_tol$it < 100) # Should converge faster

  # Test with low maxit (may not fully converge)
  fit_maxit <- mvreg_fit(y, x, z, b0, t0, maxit = 1)
  expect_equal(fit_maxit$it, 1) # Should stop after 1 iteration

  # Check that fit_tol parameters are of correct dimensions even if convergence is quicker
  expect_length(fit_tol$theta, ncol(x) + ncol(z))
})

# Test handling of incorrect method and vcov.type parameters
test_that("mvreg_fit() warns with invalid method and vcov.type", {
  expect_error(
    mvreg_fit(y, x, z, b0, t0, method = "invalid_method"),
    "should be one of"
  )
  expect_error(
    mvreg_fit(y, x, z, b0, t0, vcov.type = "invalid_vcov"),
    "should be one of"
  )
})

# Test handling of misspecification of maxit and tol
test_that("mvreg_fit() correctly corrects misspecified maxit and tol", {
  expect_warning(mvreg_fit(y, x, z, b0, t0, tol = -0.0001), "tol must be strictly positive")
  fit1 <- suppressWarnings(mvreg_fit(y, x, z, b0, t0, tol = -0.0001))
  expect_equal(fit1$tol, 0.0001)

  expect_warning(mvreg_fit(y, x, z, b0, t0, tol = 0), "tol must be strictly positive")
  fit2 <- suppressWarnings(mvreg_fit(y, x, z, b0, t0, tol = 0))
  expect_equal(fit2$tol, 1e-10)


  expect_warning(mvreg_fit(y, x, z, b0, t0, maxit = -2), "maxit must be a positive integer")
  fit3 <- suppressWarnings(mvreg_fit(y, x, z, b0, t0, maxit = -2))
  expect_equal(fit3$maxit, 2)

  expect_warning(mvreg_fit(y, x, z, b0, t0, maxit = 0), "maxit cannot be 0, set to default value of 100")
  fit4 <- suppressWarnings(mvreg_fit(y, x, z, b0, t0, maxit = 0))
  expect_equal(fit4$maxit, 100)
})
