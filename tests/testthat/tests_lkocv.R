fit <- mvreg(Sepal.Length ~ Species, data = iris)

test_that("Input validation works correctly", {
  expect_error(mvreg_lkocv("not_a_model"), "object must be a mvreg model")
  expect_error(mvreg_lkocv(fit, k = -1), "k must be a positive integer less than the number of observations")
  expect_error(mvreg_lkocv(fit, k = 151), "k must be a positive integer less than the number of observations")
  expect_error(mvreg_lkocv(fit, num_samples = -1), "num_samples must be a positive integer")
})

test_that("Leave-1-out cross-validation returns a numeric value", {
  result <- mvreg_lkocv(fit, k = 1, mc.cores = 1)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("Leave-2-out cross-validation works with approximation", {
  result <- suppressWarnings(mvreg_lkocv(fit, k = 2, mc.cores = 1))
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("Leave-2-out cross-validation works without approximation", {
  n <- nrow(iris)
  result <- mvreg_lkocv(fit, k = 2, num_samples = choose(n, 2), mc.cores = 1)
  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("Leave-k-out cross-validation handles different seeds", {
  set.seed(123)
  result1 <- suppressWarnings(mvreg_lkocv(fit, k = 2, num_samples = 500, mc.cores = 1))
  set.seed(123)
  result2 <- suppressWarnings(mvreg_lkocv(fit, k = 2, num_samples = 500, mc.cores = 1))
  expect_equal(result1, result2, tolerance = 1e-8)
})

test_that("Function warns when using approximation", {
  n <- nrow(iris)
  expect_warning(mvreg_lkocv(fit, k = 3, num_samples = 100, mc.cores = 1), "Using approximation due to large number of combinations")
})


rm("fit")

