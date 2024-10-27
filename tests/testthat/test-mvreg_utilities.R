#### get_reduced_formula() ####
test_that("get_reduced_formulas generates correct reduced formulas", {

  # Test 1: Basic formula with multiple terms
  formula1 <- y ~ x1 + x2 + x3
  reduced1 <- get_reduced_formulas("y", formula1)

  expect_equal(length(reduced1), 4)  # Expect 4 formulas (3 terms + intercept)
  expect_equal(deparse(reduced1[[1]]), deparse(as.formula("y ~ x1 + x2 + x3")))
  expect_equal(deparse(reduced1[[2]]), deparse(as.formula("y ~ x1 + x2")))
  expect_equal(deparse(reduced1[[3]]), deparse(as.formula("y ~ x1")))
  expect_equal(deparse(reduced1[[4]]), deparse(as.formula("y ~ 1")))

  # Test 2: Formula with interaction terms
  formula2 <- y ~ x1 * x2
  reduced2 <- get_reduced_formulas("y", formula2)

  expect_equal(length(reduced2), 4)  # Expect 3 formulas (2 terms + intercept)
  expect_equal(deparse(reduced2[[1]]), deparse(as.formula("y ~ x1 + x2 + x1:x2")))
  expect_equal(deparse(reduced2[[2]]), deparse(as.formula("y ~ x1 + x2")))
  expect_equal(deparse(reduced2[[3]]), deparse(as.formula("y ~ x1")))
  expect_equal(deparse(reduced2[[4]]), deparse(as.formula("y ~ 1")))

  # Test 3: Formula with a single term
  formula3 <- y ~ x1
  reduced3 <- get_reduced_formulas("y", formula3)

  expect_equal(length(reduced3), 2)  # Expect 2 formulas (1 term + intercept)
  expect_equal(deparse(reduced3[[1]]), deparse(as.formula("y ~ x1")))
  expect_equal(deparse(reduced3[[2]]), deparse(as.formula("y ~ 1")))

  # Test 4: Formula with no terms
  formula4 <- y ~ 0
  reduced4 <- get_reduced_formulas("y", formula4)

  expect_equal(length(reduced4), 1)  # Expect only intercept formula
  expect_equal(deparse(reduced4[[1]]), deparse(as.formula("y ~ 1")))

  # Test 5: Formula with multiple factors
  formula5 <- y ~ factor1 + factor2 + factor3
  reduced5 <- get_reduced_formulas("y", formula5)

  expect_equal(length(reduced5), 4)  # 3 factors + intercept
  expect_equal(deparse(reduced5[[1]]), deparse(as.formula("y ~ factor1 + factor2 + factor3")))
  expect_equal(deparse(reduced5[[2]]), deparse(as.formula("y ~ factor1 + factor2")))
  expect_equal(deparse(reduced5[[3]]), deparse(as.formula("y ~ factor1")))
  expect_equal(deparse(reduced5[[4]]), deparse(as.formula("y ~ 1")))

  # Test 6: Handling non-standard variable names
  formula6 <- y ~ `var with spaces` + `another-var`
  reduced6 <- get_reduced_formulas("y", formula6)

  expect_equal(length(reduced6), 3)  # 2 terms + intercept
  expect_equal(deparse(reduced6[[1]]), deparse(as.formula("y ~ `var with spaces` + `another-var`")))
  expect_equal(deparse(reduced6[[2]]), deparse(as.formula("y ~ `var with spaces`")))
  expect_equal(deparse(reduced6[[3]]), deparse(as.formula("y ~ 1")))

  # Test 8: Very complex formula with nested interactions
  formula7 <- y ~ x1 * (x2 + x3) + x4
  reduced7 <- get_reduced_formulas("y", formula7)

  expect_equal(length(reduced7), 7)  # 4 terms + intercept
  expect_equal(deparse(reduced7[[1]]), deparse(as.formula("y ~ x1 + x2 + x3 + x4 + x1:x2 + x1:x3")))
  expect_equal(deparse(reduced7[[2]]), deparse(as.formula("y ~ x1 + x2 + x3 + x4 + x1:x2")))
  expect_equal(deparse(reduced7[[3]]), deparse(as.formula("y ~ x1 + x2 + x3 + x4")))
  expect_equal(deparse(reduced7[[4]]), deparse(as.formula("y ~ x1 + x2 + x3")))
  expect_equal(deparse(reduced7[[5]]), deparse(as.formula("y ~ x1 + x2")))
  expect_equal(deparse(reduced7[[6]]), deparse(as.formula("y ~ x1")))
  expect_equal(deparse(reduced7[[7]]), deparse(as.formula("y ~ 1")))
})


#### are_models_nested() ####
test_that("are_models_nested correctly identifies nested models", {

  # Test 1: Basic nested models
  model1 <- lm(mpg ~ wt, data = mtcars) # model with one predictor
  model2 <- lm(mpg ~ wt + hp, data = mtcars) # model with two predictors
  expect_true(are_models_nested(model1, model2)) # model1 is nested in model2
  expect_true(are_models_nested(model2, model1)) # model2 contains model1

  # Test 2: Non-nested models
  model3 <- lm(mpg ~ qsec, data = mtcars) # unrelated predictor
  expect_false(are_models_nested(model1, model3)) # no nesting between model1 and model3
  expect_false(are_models_nested(model3, model1)) # no nesting between model3 and model1

  # Test 3: Both models with no common terms
  model4 <- lm(mpg ~ drat, data = mtcars) # another unrelated predictor
  expect_false(are_models_nested(model1, model4)) # no nesting
  expect_false(are_models_nested(model4, model1)) # no nesting

  # Test 4: Nested models with interaction terms
  model5 <- lm(mpg ~ wt * hp, data = mtcars) # interaction model
  model6 <- lm(mpg ~ wt + hp + wt:hp, data = mtcars) # interaction model
  expect_true(are_models_nested(model5, model6)) # model5 is nested in model6
  expect_true(are_models_nested(model6, model5)) # model6 cannot be nested in model5

  # Test 5: Models with only intercept
  model7 <- lm(mpg ~ 1, data = mtcars) # intercept only
  expect_true(are_models_nested(model7, model1)) # intercept model is nested in any other model
  expect_true(are_models_nested(model1, model7)) # any other model contains intercept

  # Test 6: Identical models
  model8 <- lm(mpg ~ wt + hp, data = mtcars) # same as model2
  expect_true(are_models_nested(model2, model8)) # identical models should be nested
  expect_true(are_models_nested(model8, model2)) # identical models should be nested

  # Test 7: Models with different sets of factors
  model9 <- lm(mpg ~ factor(cyl) + factor(gear), data = mtcars)
  model10 <- lm(mpg ~ factor(cyl), data = mtcars)
  expect_true(are_models_nested(model10, model9))
  expect_true(are_models_nested(model9, model10))

})



