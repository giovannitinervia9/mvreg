# Define some test data
set.seed(42)
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

  expect_true(all(c(
    "tab", "theta", "mean_vtheta", "vtheta", "it", "y",
    "n", "nsim", "seed", "k", "p", "converged",
    "non_converged", "total_time"
  ) %in% names(result)))

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
  colnames(x) <- c("(Intercept)", "x1")
  colnames(z) <- c("(Intercept)", "z1", "z2")
  result <- mvreg_simul(x, z, b, t, nsim = 100)
  expect_equal(colnames(result$theta), c("mu.const", "mu.x1", "s2.const", "s2.z1", "s2.z2"))
})

test_that("mvreg_simul calls first columns of x and z mu.const and s2.const", {
  set.seed(42)
  xd <- x[, 1]
  zd <- z[, 1]
  bd <- b[1]
  td <- t[1]
  result <- mvreg_simul(xd, zd, bd, td, nsim = 100)
  expect_equal(colnames(result$theta), c("mu.const", "s2.const"))
})

test_that("mvreg_simul use a random seed if it is not present in the environment", {
  rm(.Random.seed, envir = globalenv())
  result <- mvreg_simul(x, z, b, t, nsim = 100)
  expect_type(result$seed, "integer")
})

test_that("mvreg_simul gives error if no simulation converges", {
  set.seed(43)
  n <- 1000
  x <- cbind(1, rnorm(n, 0, 1000), rnorm(n, 0, 1000)) # Design matrix for mean component
  z <- cbind(1, rnorm(n, 0, 20), rnorm(n, 0, 20)^2) # Design matrix for variance component
  b <- rnorm(ncol(x)) # True parameters for mean component
  t <- rnorm(ncol(z)) # True parameters for variance component
  expect_error(mvreg_simul(x, z, b, t, nsim = 50), "None of the simulations converged")
})



#### print() ####

mock_simul_mvreg <- mvreg_simul(x, z, b, t, nsim = 100)

# Test printing functionality
test_that("print.simul_mvreg prints expected output", {
  # Capture the output of the print function
  output <- capture.output(print(mock_simul_mvreg))

  # Check if the output contains expected strings
  expect_true(any(grepl("Simulation study for a mvreg model", output)))
  expect_true(any(grepl("n = 100", output)))
  expect_true(any(grepl("number of simulations = 100", output)))
  expect_true(any(grepl("converged simulations = 100", output)))


  # Check if the table is printed correctly
  expect_true(any(grepl("mu.const", output)))
  expect_true(any(grepl("SE", output)))
  expect_true(any(grepl("s2.const", output)))
  expect_true(any(grepl("s2.z1", output)))
})

# Test digits argument functionality
test_that("print.simul_mvreg respects digits argument", {
  # Capture the output with different digits
  output <- capture.output(print(mock_simul_mvreg, digits = 5))

  # Check if the estimates are printed with the expected number of significant digits
  expect_true(any(grepl("1.20097", output)))
})


rm(list = c("n", "x", "z", "b", "t", "mock_simul_mvreg"))
