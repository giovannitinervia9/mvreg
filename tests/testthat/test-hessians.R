# Sample data for testing
set.seed(123)
n <- 30
x <- cbind(1, matrix(rnorm(n * 2), n, 2)) # 2 explanatory variables
z <- cbind(1, matrix(rnorm(n), n, 1)) # 1 explanatory variable for variance
y <- rnorm(n) # Response variable
b <- c(0.5, -0.2, 0.2) # Coefficients for mean
t <- c(0.3, 0.2) # Coefficients for variance

test_that("Hessian dimensions match the number of mean predictors", {
  hessian_obs <- mvreg_hessian_mu(y, x, z, b, t, type = "observed")
  hessian_exp <- mvreg_hessian_mu(y, x, z, b, t, type = "expected")
  expect_equal(dim(hessian_obs), c(ncol(x), ncol(x)))
  expect_equal(dim(hessian_exp), c(ncol(x), ncol(x)))
})


test_that("Hessian dimensions match the number of variance predictors", {
  hessian_obs <- mvreg_hessian_s2(y, x, z, b, t, type = "observed")
  hessian_exp <- mvreg_hessian_s2(y, x, z, b, t, type = "expected")
  expect_equal(dim(hessian_obs), c(ncol(z), ncol(z)))
  expect_equal(dim(hessian_exp), c(ncol(z), ncol(z)))
})


test_that("Hessian dimensions match the number of mean (row) and variance (column) predictor", {
  hessian_obs <- mvreg_hessian_mus2(y, x, z, b, t, type = "observed")
  hessian_exp <- mvreg_hessian_mus2(y, x, z, b, t, type = "expected")
  expect_equal(dim(hessian_obs), c(ncol(x), ncol(z)))
  expect_equal(dim(hessian_exp), c(ncol(x), ncol(z)))
})


test_that("Hessian matrix for mean component should be symmetric", {
  hessian_obs <- mvreg_hessian_mu(y, x, z, b, t, type = "observed")
  hessian_exp <- mvreg_hessian_mu(y, x, z, b, t, type = "expected")
  expect_true(all(hessian_obs == t(hessian_obs)))
  expect_true(all(hessian_exp == t(hessian_exp)))
})


test_that("Hessian matrix for variance component should be symmetric", {
  hessian_obs <- mvreg_hessian_s2(y, x, z, b, t, type = "observed")
  hessian_exp <- mvreg_hessian_s2(y, x, z, b, t, type = "expected")
  expect_true(all(hessian_obs == t(hessian_obs)))
  expect_true(all(hessian_exp == t(hessian_exp)))
})



test_that("Observed and expected Hessians should be equal", {
  observed_hessian <- mvreg_hessian_mu(y, x, z, b, t, type = "observed")
  expected_hessian <- mvreg_hessian_mu(y, x, z, b, t, type = "expected")
  expect_equal(observed_hessian, expected_hessian)
})


test_that("Observed and expected Hessians for variance coefficients should differ in values", {
  hessian_obs <- mvreg_hessian_s2(y, x, z, b, t, type = "observed")
  hessian_exp <- mvreg_hessian_s2(y, x, z, b, t, type = "expected")
  expect_false(all(hessian_obs == hessian_exp))
})


test_that("mvreg_hessian_mus2 returns zero matrix for expected Hessian", {
  expected_hessian <- matrix(0, ncol(x), ncol(z))

  hessian_result <- mvreg_hessian_mus2(y, x, z, b, t, type = "expected")

  expect_equal(hessian_result, expected_hessian)
})


test_that("mvreg_hessian computes full observed Hessian correctly", {
  hbb <- mvreg_hessian_mu(y, x, z, b, t, type = "observed")
  htt <- mvreg_hessian_s2(y, x, z, b, t, type = "observed")
  hbt <- mvreg_hessian_mus2(y, x, z, b, t, type = "observed")

  # Manually construct the expected full Hessian matrix
  expected_hessian <- rbind(cbind(hbb, hbt), cbind(t(hbt), htt))

  # Run the function and check if the result matches the expected Hessian
  hessian_result <- mvreg_hessian(y, x, z, b, t, type = "observed")
  expect_equal(hessian_result, expected_hessian, tolerance = 1e-6)
})


test_that("mvreg_hessian computes full expected Hessian correctly", {
  # Compute individual Hessian blocks for validation
  hbb <- mvreg_hessian_mu(y, x, z, b, t, type = "expected")
  htt <- mvreg_hessian_s2(y, x, z, b, t, type = "expected")
  hbt <- mvreg_hessian_mus2(y, x, z, b, t, type = "expected")

  # Manually construct the expected full Hessian matrix
  expected_hessian <- rbind(cbind(hbb, hbt), cbind(t(hbt), htt))

  # Run the function and check if the result matches the expected Hessian
  hessian_result <- mvreg_hessian(y, x, z, b, t, type = "expected")
  expect_equal(hessian_result, expected_hessian, tolerance = 1e-6)
})
