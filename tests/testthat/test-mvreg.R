test_that("mvreg() returns a mvreg class object", {
  expect_equal(class(mvreg(Sepal.Length ~ Species, data = iris)), "mvreg")
})


test_that("mvreg() returns correct structure and elements", {
  result <- list(mvreg(Sepal.Length ~ Species, data = iris),
   mvreg(Sepal.Length ~ Species, ~ Species, data = iris),
   mvreg(Sepal.Length ~ Species, Sepal.Length ~ Species, data = iris))

  # Check that result contains all expected elements
  lapply(result, function(result){
    expect_type(result$coefficients, "double")
    expect_type(result$coefficients.mu, "double")
    expect_type(result$coefficients.s2, "double")
    expect_type(result$vcov, "double")
    expect_type(result$vcov.mu, "double")
    expect_type(result$vcov.s2, "double")
    expect_type(result$logLik, "double")
    expect_type(result$fit.mu, "double")
    expect_type(result$fit.log.s2, "double")
    expect_type(result$fit.s2, "double")
    expect_type(result$residuals, "double")
    expect_type(result$it, "integer")
    expect_type(result$start, "double")
    expect_type(result$y, "double")
    expect_type(result$xd, "double")
    expect_type(result$zd, "double")
    expect_type(result$nobs, "integer")
    expect_type(result$df.residual, "integer")
    expect_type(result$call, "language")
    expect_type(result$response, "character")
    expect_type(result$colx, "character")
    expect_type(result$colz, "character")
    expect_type(result$formula.mu, "language")
    expect_type(result$formula.s2, "language")
    expect_type(result$terms.mu, "language")
    expect_type(result$terms.s2, "language")

  })

})



test_that("mvreg() works with a single formula (same model for mean and variance)", {
  result <- mvreg(Sepal.Length ~ Species, data = iris)

  # Check if the object is of class 'mvreg'
  expect_s3_class(result, "mvreg")

  # Check if coefficients are correctly returned
  expect_true(!is.null(result$coefficients))
  expect_equal(length(result$coefficients.mu), length(result$coefficients.s2))
  expect_equal(length(result$coefficients.mu) + length(result$coefficients.s2), length(result$coefficients))

  # Check if log-likelihood is computed
  expect_true(!is.null(result$logLik))

  # Check if residuals are computed
  expect_equal(length(result$residuals), nrow(iris))

  # Check if vcov matrices for mean and variance components are computed
  expect_true(is.matrix(result$vcov.mu))
  expect_true(is.matrix(result$vcov.s2))
  expect_equal(dim(result$vcov.mu), dim(result$vcov.s2))
  expect_equal(dim(result$vcov.mu) + dim(result$vcov.s2), dim(result$vcov))

  # Check if the two formulas are equal
  expect_equal(result$formula.mu, result$formula.s2)
})


test_that("mvreg() works with different formulas for mean and variance components", {
  result <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris)

  # Check if the object is of class 'mvreg'
  expect_s3_class(result, "mvreg")

  # Check if coefficients are correctly returned
  expect_true(!is.null(result$coefficients))
  expect_equal(length(result$coefficients.mu) + length(result$coefficients.s2), length(result$coefficients))

  # Check if log-likelihood is computed
  expect_true(!is.null(result$logLik))

  # Check if residuals are computed
  expect_equal(length(result$residuals), nrow(iris))

  # Check if vcov matrices for mean and variance components are computed
  expect_true(is.matrix(result$vcov.mu))
  expect_true(is.matrix(result$vcov.s2))
  expect_equal(dim(result$vcov.mu) + dim(result$vcov.s2), dim(result$vcov))

  # Check if the two formulas are different
  expect_true(result$formula.mu != result$formula.s2)
})


test_that("mvreg() gives errors when NAs are present in response variable", {
  iris_with_na <- iris
  iris_with_na$Sepal.Length[1] <- NA

  # Check if the function throws an error for missing values
  expect_error(mvreg(Sepal.Length ~ Species, Sepal.Length ~ Species, data = iris_with_na))
})

test_that("mvreg() gives errors when NAs are present in mean component", {
  iris_with_na <- iris
  iris_with_na$Sepal.Width[1] <- NA

  # Check if the function throws an error for missing values
  expect_error(mvreg(Sepal.Length ~ Sepal.Width, Sepal.Length ~ Species, data = iris_with_na))
})

test_that("mvreg() gives errors when NAs are present in variance component", {
  iris_with_na <- iris
  iris_with_na$Species[1] <- NA

  # Check if the function throws an error for missing values
  expect_error(mvreg(Sepal.Length ~ Sepal.Width, Sepal.Length ~ Species, data = iris_with_na))
})

test_that("If formula.s2 is not specified, formula.mu is used as formula.s2", {
  result <- mvreg(Sepal.Length ~ Species, data = iris)
  expect_equal(result$formula.mu, result$formula.s2)
})

test_that("If in formula.s2 the left-hand-side is not specified than it is taken from formula.mu", {
  result <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris)
  expect_equal(result$formula.s2[[2]], result$formula.mu[[2]])
})


test_that("If in formula.s2 the left-hand-side and it is the same as in formula.mu, formula.s2 is correctly built", {
  result <- mvreg(Sepal.Length ~ Species, Sepal.Length ~ Sepal.Width, data = iris)
  expect_equal(result$formula.s2[[2]], result$formula.mu[[2]])
})

test_that("mvreg() stops if response variable appears in the left-hand-side of formula,s2", {
  expect_error(
    mvreg(Sepal.Length ~ Species, ~Sepal.Length, data = iris),
    "response variable cannot appear on right-hand side of the formula"
  )
})


test_that("mvreg() convergence criteria and iteration count", {

  result <- mvreg(Sepal.Length ~ Species, data = iris, tol = 1e-8, maxit = 5)

  # Check if the iteration count is returned correctly
  expect_true(result$it <= 5)
  expect_true(result$it > 0)
})


test_that("mvreg() produces consistent output with different estimation methods", {
  # Test with method "wls"
  result_wls <- mvreg(Sepal.Length ~ Species, data = iris, method = "wls")

  # Test with method "full_nr"
  result_full_nr <- mvreg(Sepal.Length ~ Species, data = iris, method = "full_nr")

  # Check if both results are of class "mvreg"
  expect_s3_class(result_wls, "mvreg")
  expect_s3_class(result_full_nr, "mvreg")

  # Check if results are consistent in terms of structure
  expect_equal(names(result_wls), names(result_full_nr))
})



test_that("mvreg() handles different initial values for variance component", {
  # Using different initial values for the variance component
  result_residuals <- mvreg(Sepal.Length ~ Species, data = iris, start.s2 = "residuals")
  result_gamma <- mvreg(Sepal.Length ~ Species, data = iris, start.s2 = "gamma")
  result_zero <- mvreg(Sepal.Length ~ Species, data = iris, start.s2 = "zero")

  # Check if each result is of class "mvreg"
  expect_s3_class(result_residuals, "mvreg")
  expect_s3_class(result_gamma, "mvreg")
  expect_s3_class(result_zero, "mvreg")

  # Check if they all have consistent structure
  expect_equal(names(result_residuals), names(result_gamma))
  expect_equal(names(result_gamma), names(result_zero))
})


test_that("If data is not passed, then call doesn't contain data specification", {
  # Temporarily assign `y` and `x` to the global environment
  original_y <- if (exists("y", .GlobalEnv)) get("y", .GlobalEnv) else NULL
  original_x <- if (exists("x", .GlobalEnv)) get("x", .GlobalEnv) else NULL

  assign("y", iris$Sepal.Length, envir = .GlobalEnv)
  assign("x", iris$Species, envir = .GlobalEnv)

  # Run mvreg without the `data` argument
  result <- mvreg(y ~ x)

  # Test to ensure that "data" is not in the function call
  expect_true(!("data" %in% names(result$call)))

  # Restore the original `y` and `x` in the global environment (if they existed)
  if (!is.null(original_y)) {
    assign("y", original_y, .GlobalEnv)
  } else {
    rm("y", envir = .GlobalEnv)
  }

  if (!is.null(original_x)) {
    assign("x", original_x, .GlobalEnv)
  } else {
    rm("x", envir = .GlobalEnv)
  }
})



test_that("If data is passed, then call contains data specification", {
  result <- mvreg(Sepal.Length ~ Species, data = iris)
  expect_true(("data" %in% names(result$call)))
  expect_equal(as.character(result$call$data), "iris")
})


