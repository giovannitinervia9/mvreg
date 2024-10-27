# Tests for mvreg_gradient_mu
test_that("mvreg_gradient_mu computes correct gradient for mean coefficients", {
  # Sample data for testing
  y <- c(1.5, 2.0, 2.5)
  x <- matrix(c(1, 1, 1, 0.5, -0.5, 0.5), nrow = 3, ncol = 2) # intercept and predictor for mean
  z <- matrix(c(1, 1, 1, 0.3, -0.2, 0.4), nrow = 3, ncol = 2) # intercept and predictor for variance

  # Parameters for mean and variance components
  b <- c(0.2, 0.5)
  t <- c(0.1, -0.1)

  # Manual calculation of expected gradient
  eta.mu <- as.vector(x %*% b)
  eta.s2 <- as.vector(z %*% t)
  const <- (y - eta.mu) / exp(eta.s2)
  expected_gradient <- apply(x, 2, function(xi) sum(const * xi))

  # Calculated gradient
  computed_gradient <- mvreg_gradient_mu(y, x, z, b, t)

  # Check computed gradient against expected values
  expect_equal(computed_gradient, expected_gradient, tolerance = 1e-6)
})

test_that("mvreg_gradient_mu returns numeric vector of correct length", {
  y <- rnorm(5)
  x <- matrix(rnorm(10), ncol = 2)
  z <- matrix(rnorm(10), ncol = 2)
  b <- rnorm(2)
  t <- rnorm(2)

  # Compute gradient
  gradient_mu <- mvreg_gradient_mu(y, x, z, b, t)

  # Test output type and length
  expect_type(gradient_mu, "double")
  expect_length(gradient_mu, ncol(x))
})

# Tests for mvreg_gradient_s2
test_that("mvreg_gradient_s2 computes correct gradient for variance coefficients", {
  # Sample data for testing
  y <- c(1.5, 2.0, 2.5)
  x <- matrix(c(1, 1, 1, 0.5, -0.5, 0.5), nrow = 3, ncol = 2) # intercept and predictor for mean
  z <- matrix(c(1, 1, 1, 0.3, -0.2, 0.4), nrow = 3, ncol = 2) # intercept and predictor for variance

  # Parameters for mean and variance components
  b <- c(0.2, 0.5)
  t <- c(0.1, -0.1)

  # Manual calculation of expected gradient
  eta.mu <- as.vector(x %*% b)
  eta.s2 <- as.vector(z %*% t)
  const <- 1 - (y - eta.mu)^2 * exp(-eta.s2)
  expected_gradient <- apply(z, 2, function(zi) -0.5 * sum(const * zi))

  # Calculated gradient
  computed_gradient <- mvreg_gradient_s2(y, x, z, b, t)

  # Check computed gradient against expected values
  expect_equal(computed_gradient, expected_gradient, tolerance = 1e-6)
})

test_that("mvreg_gradient_s2 returns numeric vector of correct length", {
  y <- rnorm(5)
  x <- matrix(rnorm(10), ncol = 2)
  z <- matrix(rnorm(10), ncol = 2)
  b <- rnorm(2)
  t <- rnorm(2)

  # Compute gradient
  gradient_s2 <- mvreg_gradient_s2(y, x, z, b, t)

  # Test output type and length
  expect_type(gradient_s2, "double")
  expect_length(gradient_s2, ncol(z))
})

# Edge Case Tests
test_that("mvreg_gradient_mu and mvreg_gradient_s2 handle zero variance cases", {
  y <- c(1, 1, 1)
  x <- matrix(c(1, 1, 1, 0, 0, 0), ncol = 2)
  z <- matrix(c(1, 1, 1, 0, 0, 0), ncol = 2)
  b <- c(0, 0)
  t <- c(0, 0)

  # Ensure that the output is numeric and does not contain NaNs
  gradient_mu <- mvreg_gradient_mu(y, x, z, b, t)
  gradient_s2 <- mvreg_gradient_s2(y, x, z, b, t)

  expect_type(gradient_mu, "double")
  expect_type(gradient_s2, "double")
  expect_false(any(is.nan(gradient_mu)))
  expect_false(any(is.nan(gradient_s2)))
})

