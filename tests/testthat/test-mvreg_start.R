# Example data
n <- 100
x1 <- rnorm(n)
x2 <- factor(sample(letters[1:3], n, TRUE))
x <- model.matrix(~ x1 + x2)
z1 <- factor(sample(letters[1:3], n, TRUE))
z <- model.matrix(~z1)
b <- rnorm(ncol(x))
t <- rnorm(ncol(z))
y <- rnorm(n, mean = x %*% b, sd = sqrt(exp(z %*% t)))



test_that("mvreg_start returns a list with correct structure", {
  result <- mvreg_start(y, x, z, start.s2 = "residuals")
  expect_type(result, "list")
  expect_named(result, c("start", "b0", "t0"))
  expect_length(result$start, ncol(x) + ncol(z))
  expect_length(result$b0, ncol(x))
  expect_length(result$t0, ncol(z))
})

# Test for start.s2 = "residuals"
test_that("mvreg_start with start.s2 = 'residuals' returns reasonable starting values", {
  result <- mvreg_start(y, x, z, start.s2 = "residuals")
  expect_true(is.numeric(result$start))
  expect_true(all(result$start != 0))
})

# Test for start.s2 = "gamma"
test_that("mvreg_start with start.s2 = 'gamma' returns reasonable starting values", {
  result <- mvreg_start(y, x, z, start.s2 = "gamma")
  expect_true(is.numeric(result$start))
  expect_true(all(result$start != 0))
})

# Test for start.s2 = "zero"
test_that("mvreg_start with start.s2 = 'zero' returns correct starting values", {
  result <- mvreg_start(y, x, z, start.s2 = "zero")
  expect_true(is.numeric(result$start))
  expect_equal(c(result$start[ncol(x) + 1], use.names = F), log(var(y)))
  expect_true(all(c(result$start[-(1:(ncol(x) + 1))], use.names = F) == 0))
})

# # Test for invalid input types
# test_that("mvreg_start handles invalid input types", {
#   expect_error(mvreg_start("invalid", x, z, start.s2 = "residuals"), "non-numeric argument")
#   expect_error(mvreg_start(y, "invalid", z, start.s2 = "residuals"), "non-numeric argument")
#   expect_error(mvreg_start(y, x, "invalid", start.s2 = "residuals"), "non-numeric argument")
# })

# Test for different methods in start.s2
test_that("mvreg_start handles different start.s2 methods", {
  result_resid <- mvreg_start(y, x, z, start.s2 = "residuals")
  result_gamma <- mvreg_start(y, x, z, start.s2 = "gamma")
  result_zero <- mvreg_start(y, x, z, start.s2 = "zero")
  expect_true(length(result_resid$start) == length(result_gamma$start) &&
    length(result_gamma$start) == length(result_zero$start))
})


# Test handling of incorrect start.s2 parameter
test_that("mvreg_start() warns with invalid start.s2", {
  expect_error(
    mvreg_start(y, x, z, start.s2 = "invalid_start.s2"),
    "should be one of"
  )
})


rm(list = c("n", "x1", "x2", "x", "z1", "z", "b", "t", "y"))
