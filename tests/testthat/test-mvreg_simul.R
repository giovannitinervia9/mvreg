# Define some test data
n <- 100
x <- cbind(1, rnorm(n)) # Design matrix for mean component
z <- cbind(1, x[, 2], x[, 2]^2) # Design matrix for variance component
b <- rnorm(ncol(x)) # True parameters for mean component
t <- rnorm(ncol(z)) # True parameters for variance component

test_that("mvreg_simul runs without errors for valid inputs", {
  expect_no_error(
    result <- mvreg_simul(x, z, b, t, nsim = 100, seed = 42)
  )

  expect_true(is.list(result))
  expect_s3_class(result, "simul_mvreg")
})

test_that("mvreg_simul returns the expected components", {
  result <- mvreg_simul(x, z, b, t, nsim = 100, seed = 42)

  expect_true(all(c("tab", "theta", "mean_vtheta", "vtheta", "it", "y",
                    "n", "nsim", "seed", "k", "p", "converged",
                    "non_converged", "total_time") %in% names(result)))

  expect_equal(nrow(result$tab), ncol(x) + ncol(z)) # Should match number of parameters
  expect_equal(ncol(result$theta), ncol(x) + ncol(z)) # Should match number of parameters
  expect_equal(nrow(result$y), 100) # Check number of observations in y
  expect_equal(ncol(result$y), result$converged) # Check correct number of simulations
  expect_equal(nrow(result$theta), result$converged) # Check correct number of simulations
  expect_equal(length(result$vtheta), result$converged) # Check correct number of simulations
  expect_equal(length(result$it), result$converged) # Check correct number of simulations
})

test_that("mvreg_simul handles input errors correctly", {
  expect_error(
    mvreg_simul(x, z[1, , drop = FALSE], b, t),
    "x and z must have the same number of rows"
  )

  expect_error(
    mvreg_simul(x, z, b[-1], t),
    "b must have length equal to the number of columns in x"
  )

  expect_error(
    mvreg_simul(x, z, b, t[-1]),
    "t must have length equal to the number of columns in z"
  )

  expect_error(
    mvreg_simul(x, z, b, t, nsim = -10),
    "nsim must be a positive integer"
  )
})

test_that("mvreg_simul returns consistent results with fixed seed", {
  set.seed(42)
  result1 <- mvreg_simul(x, z, b, t, nsim = 100)

  set.seed(42)
  result2 <- mvreg_simul(x, z, b, t, nsim = 100)

  expect_equal(result1$tab, result2$tab)
  expect_equal(result1$theta, result2$theta)
  expect_equal(result1$y, result2$y)
  expect_equal(result1$mean_vtheta, result2$mean_vtheta)
  expect_equal(result1$vtheta, result2$vtheta)
})

test_that("mvreg_simul returns consistent names of coefficients that starts with mu and s2", {
  set.seed(42)
  result <- mvreg_simul(x, z, b, t, nsim = 100)
  expect_equal(colnames(result$theta), c("mu.const", "mu.x1", "s2.const", "s2.z1", "s2.z2"))
})

