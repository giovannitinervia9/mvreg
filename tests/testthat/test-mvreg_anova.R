#### anova() ####

test_that("anova accepts only mvreg objects", {
  mod1 <- mvreg(Sepal.Length ~ Species + Sepal.Width, data = iris)
  non_mvreg <- lm(Sepal.Length ~ Species + Sepal.Width, data = iris)
  expect_error(anova.mvreg(non_mvreg), "All objects must be of class 'mvreg'")
})

test_that("anova performs single model comparison", {
  mod <- mvreg(Sepal.Length ~ Species + Sepal.Width, data = iris)
  result <- anova(mod)

  expect_true("mu.tests" %in% names(result))
  expect_true("s2.tests" %in% names(result))

  # Check the result structure
  expect_s3_class(result$mu.tests, "data.frame")
  expect_s3_class(result$s2.tests, "data.frame")

  # Verify at least two reduced models were compared
  expect_gt(nrow(result$mu.tests), 1)
  expect_gt(nrow(result$s2.tests), 1)
})

test_that("anova performs multiple model comparison", {
  mod1 <- mvreg(Sepal.Length ~ Species + Sepal.Width, data = iris)
  mod2 <- mvreg(Sepal.Length ~ Species, ~1, data = iris)

  result <- anova(mod1, mod2)

  # Check structure and presence of key fields
  expect_true("tests" %in% names(result))
  expect_s3_class(result$tests, "data.frame")

  # Verify the warning for non-nested models, if applicable
  expect_warning(anova(mod1, mod1), "LRT test with models with the same number of parameters is meaningless")
})

test_that("anova orders models correctly when order.models is TRUE", {
  mod1 <- mvreg(Sepal.Length ~ Species + Sepal.Width, data = iris)
  mod2 <- mvreg(Sepal.Length ~ Species, data = iris)

  result_ordered <- anova(mod1, mod2, order.models = TRUE)
  params <- result_ordered$tests$n.param

  # Check that parameters are in ascending order
  expect_true(all(diff(params) >= 0))

  result_unordered <- suppressWarnings(anova(mod1, mod2, order.models = FALSE))
  expect_warning(anova(mod1, mod2, order.models = FALSE))
  # Ensure no reordering occurs
  expect_equal(rownames(result_unordered$tests), c("model1", "model2"))
})

test_that("anova calculations are consistent", {
  mod1 <- mvreg(Sepal.Length ~ Species + Sepal.Width, data = iris)
  mod2 <- mvreg(Sepal.Length ~ Species, data = iris)

  result <- anova(mod1, mod2)

  # Test LRT consistency
  lrt_value <- result$tests$LRT[2]
  df_value <- result$tests$df[2]

  # Calculate p-value manually and check for consistency
  expect_equal(result$tests$`Pr(>Chi)`[2], 1 - pchisq(lrt_value, df = df_value), tolerance = 1e-5)
})


test_that("anova doesn't accept models with different responses", {
  mod1 <- mvreg(Sepal.Length ~ Species + Sepal.Width, data = iris)
  mod2 <- mvreg(Petal.Length ~ Species, data = iris)
  expect_error(anova(mod1, mod2), "Models have different response variables")
})


test_that("anova doesn't accept models with different number of observation", {
  mod1 <- mvreg(Sepal.Length ~ Species + Sepal.Width, data = iris)
  mod2 <- mvreg(Sepal.Length ~ Species, data = iris[1:130, ])
  expect_error(anova(mod1, mod2), "Models are fitted with different numbers of observations")
})

test_that("anova raises a warning when non nested models are passed", {
  mod1 <- mvreg(Sepal.Length ~ Species + Sepal.Width, data = iris)
  mod2 <- mvreg(Sepal.Length ~ Species, ~Petal.Width, data = iris)
  expect_warning(anova(mod1, mod2), "LRT test with non nested models is meaningless")
})



#### print.anova() ####
mod1 <- mvreg(Sepal.Length ~ Species + Sepal.Width, data = iris)
mod2 <- mvreg(Sepal.Length ~ Species, ~1, data = iris)
mod3 <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris)

# Test that the print.anova function produces expected output for a single model
# Test that the print.anova function produces expected output for a single model
test_that("print.anova displays correctly for single model input", {
  anova_result <- anova(mod3)

  # Use regular expressions to allow for flexibility in whitespace and special characters
  expect_output(print(anova_result), "Response:\\s*Sepal\\.Length")
  expect_output(print(anova_result), "Mean Model Comparison \\(variance component taken as in the full model\\)")
  expect_output(print(anova_result), "Variance Model Comparison \\(mean component taken as in the full model\\):")
  expect_output(print(anova_result), "Species")
  expect_output(print(anova_result), "Sepal\\.Width")
})


# Test that the function displays correctly for multiple model inputs
test_that("print.anova displays correctly for multiple models comparison", {
  anova_result_multi <- anova(mod1, mod2)

  expect_output(print(anova_result_multi), "Model comparison: mean and variance components")
  expect_output(print(anova_result_multi), "Response: Sepal.Length")
  expect_output(print(anova_result_multi), "Mean component formulas:")
  expect_output(print(anova_result_multi), "Variance component formulas:")
  expect_output(print(anova_result_multi), "model1")
  expect_output(print(anova_result_multi), "model2")
})

# Test that digits argument changes output precision by directly setting the digits option
test_that("print.anova respects digits argument for output precision", {
  # Temporarily store the current option for digits
  original_digits <- getOption("digits")

  # Test with lower precision
  options(digits = 4)
  anova_result <- anova(mod3)
  expect_output(print(anova_result, digits = 2), "\\d\\.\\d{2}")

  # Test with higher precision
  options(digits = 8)
  expect_output(print(anova_result, digits = 5), "\\d\\.\\d{5}")

  # Restore the original digits option
  options(digits = original_digits)
})

# Test that additional arguments (...) are accepted without errors
test_that("print.anova accepts additional arguments", {
  anova_result <- anova(mod3)
  expect_no_warning(print(anova_result, digits = 3, another_unused_arg = "test"))
})

rm(list = c("mod1", "mod2", "mod3"))
