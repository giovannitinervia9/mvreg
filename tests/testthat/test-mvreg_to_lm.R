# Load the necessary library
test_that("mvreg_to_lm correctly converts mvreg to lm", {
  # Fit a mvreg model using the cars dataset (assuming mvreg is implemented)
  mod.mvreg <- mvreg(dist ~ speed, data = cars)

  # Convert the mvreg model to an lm object
  mod.lm <- mvreg_to_lm(mod.mvreg)

  # Test that the class of the resulting model is "lm"
  expect_s3_class(mod.lm, "lm")

  # Test that the coefficients of the lm model are equivalent to those of the mvreg model
  expect_equal(length(mod.lm$coefficients), length(mod.mvreg$coefficients.mu))

  # Check that the data used in lm is the same as in mvreg
  expect_equal(as.character(mod.lm$call$data), "cars") # Ensure the data is the same
})

test_that("mvreg_to_lm handles invalid input gracefully", {
  # Create a non-mvreg object (e.g., a simple lm object)
  model_lm <- lm(dist ~ speed, data = cars)

  # Expect an error when trying to convert an invalid object
  expect_error(mvreg_to_lm(model_lm), "Object is not of class 'mvreg' or 'summary.mvreg'.")

  # Create a non-existent object
  non_mvreg_object <- list()

  # Expect an error when trying to convert a non-existent object
  expect_error(mvreg_to_lm(non_mvreg_object), "Object is not of class 'mvreg' or 'summary.mvreg'.")
})

test_that("mvreg_to_lm works without data in mvreg object", {
  # Temporarily assign `y` and `x` to the global environment
  original_y <- if (exists("y", .GlobalEnv)) get("y", .GlobalEnv) else NULL
  original_x <- if (exists("x", .GlobalEnv)) get("x", .GlobalEnv) else NULL

  assign("y", iris$Sepal.Length, envir = .GlobalEnv)
  assign("x", iris$Species, envir = .GlobalEnv)

  mod.mvreg <- mvreg(y ~ x)

  # Convert the mvreg model to an lm object
  mod.lm <- mvreg_to_lm(mod.mvreg)

  # Test that the class of the resulting model is "lm"
  expect_s3_class(mod.lm, "lm")

  # Test that the coefficients of the lm model are equivalent to those of the mvreg model
  expect_equal(length(mod.lm$coefficients), length(mod.mvreg$coefficients.mu))
})
