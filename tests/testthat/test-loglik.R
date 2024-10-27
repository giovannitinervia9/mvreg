test_that("mvreg_loglik computes correct log-likelihood value", {
  # Example data
  y <- c(2.5, 3.0, 3.5)
  x <- matrix(c(1, 1, 1, 2, 3, 4), nrow = 3, ncol = 2) # Intercept and predictor for mean component
  z <- matrix(c(1, 1, 1, 0.5, -0.5, 0.5), nrow = 3, ncol = 2) # Intercept and predictor for variance component

  # Parameters
  b <- c(0.5, 1.0) # Intercept and slope for mean component
  t <- c(0.1, 0.2) # Intercept and slope for variance component

  # Expected log-likelihood calculated manually
  expected_loglik <- -0.5 * sum(log(2 * pi) + z %*% t + (y - x %*% b)^2 / exp(z %*% t))

  # Log-likelihood computed with the function
  computed_loglik <- mvreg_loglik(y, x, z, b, t)

  # Comparison
  expect_equal(computed_loglik, expected_loglik, tolerance = 1e-6)
})

test_that("mvreg_loglik handles cases with constant variance", {
  # Case with constant variance
  y <- c(1, 2, 3)
  x <- matrix(c(1, 1, 1, 1, 2, 3), ncol = 2)
  z <- matrix(c(1, 1, 1, 0, 0, 0), ncol = 2)
  b <- c(0, 1)
  t <- c(5, 0) # Constant variance

  computed_loglik <- mvreg_loglik(y, x, z, b, t)
  expect_true(is.numeric(computed_loglik))
  expect_false(is.nan(computed_loglik)) # Check to avoid NaN
})

test_that("mvreg_loglik returns NA or error with invalid input sizes", {
  # Length of b or t does not match dimensions of x or z
  y <- c(1, 2, 3)
  x <- matrix(c(1, 1, 1, 1, 2, 3), ncol = 2)
  z <- matrix(c(1, 1, 1, 0, 0, 0), ncol = 2)

  b <- c(0, 1, 2) # Incorrect length of b
  t <- c(0, 0)

  expect_error(mvreg_loglik(y, x, z, b, t), "non-conformable arguments") # Error for incorrect length of b

  b <- c(0, 1)
  t <- c(0, 0, 0) # Incorrect length of t
  expect_error(mvreg_loglik(y, x, z, b, t), "non-conformable arguments") # Error for incorrect length of t
})

test_that("mvreg_loglik handles large input values correctly", {
  # Data with large values to test numerical stability
  y <- c(1e10, 2e10, 3e10)
  x <- matrix(c(1, 1, 1, 1e2, 1e2, 1e2), ncol = 2)
  z <- matrix(c(1, 1, 1, 1e2, -1e2, 1e2), ncol = 2)
  b <- c(1e10, 1e10)
  t <- c(1e-2, -1e-2)

  computed_loglik <- mvreg_loglik(y, x, z, b, t)
  expect_true(is.finite(computed_loglik)) # Check to avoid overflow or infinities
})

test_that("mvreg_loglik returns numeric output for valid inputs", {
  # Check that the output is numeric for valid data
  y <- rnorm(10)
  x <- matrix(rnorm(20), ncol = 2)
  z <- matrix(rnorm(20), ncol = 2)
  b <- rnorm(2)
  t <- rnorm(2)

  computed_loglik <- mvreg_loglik(y, x, z, b, t)
  expect_type(computed_loglik, "double")
  expect_length(computed_loglik, 1)
})
