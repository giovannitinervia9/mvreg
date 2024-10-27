#### vcov() ####
test_that("vcov matrix is correctly computed", {

  mod.obs <- mvreg(Sepal.Length ~ Species, data = iris, vcov.type = "observed")
  k <- ncol(mod.obs$xd)
  p <- ncol(mod.obs$zd)
  vtheta.obs <- -solve(mvreg_hessian(mod.obs$y, mod.obs$xd, mod.obs$zd, mod.obs$coefficients.mu, mod.obs$coefficients.s2, type = "observed"))
  expect_equal(vtheta.obs, vcov(mod.obs))
  expect_equal(vtheta.obs[1:k, 1:k], vcov(mod.obs, "mu"))
  expect_equal(vtheta.obs[(k+1):(k+p), (k+1):(k+p)], vcov(mod.obs, "s2"))

  mod.exp <- mvreg(Sepal.Length ~ Species, data = iris, vcov.type = "expected")
  k <- ncol(mod.exp$xd)
  p <- ncol(mod.exp$zd)
  vtheta.obs <- -solve(mvreg_hessian(mod.exp$y, mod.exp$xd, mod.exp$zd, mod.exp$coefficients.mu, mod.exp$coefficients.s2, type = "observed"))
  expect_equal(vtheta.obs, vcov(mod.exp))
  expect_equal(vtheta.obs[1:k, 1:k], vcov(mod.exp, "mu"))
  expect_equal(vtheta.obs[(k+1):(k+p), (k+1):(k+p)], vcov(mod.exp, "s2"))

})

test_that("vcov.mvreg() handles invalid partition argument gracefully", {
  mod <- mvreg(Sepal.Length ~ Species, data = iris)

  expect_error(vcov(mod, partition = "invalid"),
               "should be one") # Customize the error message as needed
})


#### summary() ####


test_that("summary.mvreg returns correct summary object",{
  mod <- mvreg(Sepal.Length ~ Species, data = iris)
  s <- summary(mod)
  expect_type(s, "list")
  expect_equal(class(s), "summary.mvreg")
  expect_equal(names(s), c("call", "residuals", "coefficients", "coefficients.mu",
                           "coefficients.s2", "df", "vcov",
                           "vcov.mu", "vcov.s2", "loglik", "AIC", "BIC"))
  expect_equal(nrow(s$coefficients), length(mod$coefficients))
  expect_equal(nrow(s$coefficients.mu), length(mod$coefficients.mu))
  expect_equal(nrow(s$coefficients.s2), length(mod$coefficients.s2))
  expect_type(s$call, "language")
  expect_type(s$residuals, "double")
  expect_type(s$coefficients, "list")
  expect_type(s$coefficients.mu, "list")
  expect_type(s$coefficients.s2, "list")
  expect_type(s$df, "integer")
  expect_type(s$vcov, "double")
  expect_type(s$vcov.mu, "double")
  expect_type(s$vcov.s2, "double")
  expect_type(s$loglik, "double")
  expect_type(s$AIC, "double")
  expect_type(s$BIC, "double")
  })




#### fitted() ####

test_that("fitted.mvreg returns correct fitted values", {
  mod <- mvreg(Sepal.Length ~ Species, data = iris)
  fit.all <- fitted(mod)
  fit.mu <- fitted(mod, "mu")
  fit.ls2 <- fitted(mod, "log.s2")
  fit.s2 <- fitted(mod, "s2")

  expect_type(fit.all, "list")
  expect_equal(ncol(fit.all), 3)
  expect_equal(nrow(fit.all), 150)
  expect_equal(colnames(fit.all), c("fit.mu", "fit.log.s2", "fit.s2"))

  expect_type(fit.mu, "double")
  expect_length(fit.mu, 150)

  expect_type(fit.ls2, "double")
  expect_length(fit.ls2, 150)

  expect_type(fit.s2, "double")
  expect_length(fit.s2, 150)

  })

test_that("fitted.mvreg handles invalid type argument gracefully", {
  mod <- mvreg(Sepal.Length ~ Species, data = iris)

  expect_error(fitted(mod, type = "invalid"),
               "should be one") # Customize the error message as needed
})



#### predict() ####

test_that("predict.mvreg returns correct predictions for different components and settings", {
  # Fit a model
  mvreg_mod <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris)

  # Test: Default prediction (all components, no newdata, no SE, no interval)
  pred_all <- predict(mvreg_mod)
  expect_type(pred_all, "list")
  expect_named(pred_all, c("mu", "log.s2", "s2"))

  # Test: Predict mean component only
  pred_mu <- predict(mvreg_mod, type = "mu")
  expect_type(pred_mu, "list")
  expect_type(pred_mu$mu, "double")
  expect_length(pred_mu$mu, nrow(iris))

  # Test: Predict log variance component only
  pred_log_s2 <- predict(mvreg_mod, type = "log.s2")
  expect_type(pred_log_s2, "list")
  expect_type(pred_log_s2$log.s2, "double")
  expect_length(pred_log_s2$log.s2, nrow(iris))

  # Test: Predict variance component only
  pred_s2 <- predict(mvreg_mod, type = "s2")
  expect_type(pred_s2, "list")
  expect_type(pred_s2$s2, "double")
  expect_length(pred_s2$s2, nrow(iris))

  # Test: Prediction with standard errors (all components)
  pred_se <- predict(mvreg_mod, se.fit = TRUE)
  expect_type(pred_se, "list")
  expect_named(pred_se, c("mu", "log.s2", "s2"))
  expect_s3_class(pred_se$mu, "data.frame")
  expect_named(pred_se$mu, c("pred.mu", "se.mu"))
  expect_s3_class(pred_se$log.s2, "data.frame")
  expect_named(pred_se$log.s2, c("pred.log.s2", "se.log.s2"))
  expect_s3_class(pred_se$s2, "data.frame")
  expect_named(pred_se$s2, c("pred.s2", "se.s2"))

  # Test: Prediction with confidence intervals (all components)
  pred_interval <- predict(mvreg_mod, interval = TRUE)
  expect_type(pred_interval, "list")
  expect_named(pred_interval, c("mu", "log.s2", "s2"))
  expect_s3_class(pred_interval$mu, "data.frame")
  expect_named(pred_interval$mu, c("pred.mu", "lwr0.95", "upr0.95"))
  expect_s3_class(pred_interval$log.s2, "data.frame")
  expect_named(pred_interval$log.s2, c("pred.log.s2", "lwr0.95", "upr0.95"))
  expect_s3_class(pred_interval$s2, "data.frame")
  expect_named(pred_interval$s2, c("pred.s2", "lwr0.95", "upr0.95"))

  # Test: Prediction with standard errors and confidence intervals (all components)
  pred_se_interval <- predict(mvreg_mod, se.fit = TRUE, interval = TRUE)
  expect_type(pred_se_interval, "list")
  expect_named(pred_se_interval, c("mu", "log.s2", "s2"))
  expect_s3_class(pred_se_interval$mu, "data.frame")
  expect_named(pred_se_interval$mu, c("pred.mu", "se.mu", "lwr0.95", "upr0.95"))
  expect_s3_class(pred_se_interval$log.s2, "data.frame")
  expect_named(pred_se_interval$log.s2, c("pred.log.s2", "se.log.s2", "lwr0.95", "upr0.95"))
  expect_s3_class(pred_se_interval$s2, "data.frame")
  expect_named(pred_se_interval$s2, c("pred.s2", "se.s2", "lwr0.95", "upr0.95"))

  # Test: Prediction with custom confidence level (e.g., 99%)
  pred_interval_99 <- predict(mvreg_mod, interval = TRUE, level = 0.99)
  expect_type(pred_interval_99, "list")
  expect_s3_class(pred_interval_99$mu, "data.frame")
  expect_named(pred_interval_99$mu, c("pred.mu", "lwr0.99", "upr0.99"))
  expect_s3_class(pred_interval_99$log.s2, "data.frame")
  expect_named(pred_interval_99$log.s2, c("pred.log.s2", "lwr0.99", "upr0.99"))
  expect_s3_class(pred_interval_99$s2, "data.frame")
  expect_named(pred_interval_99$s2, c("pred.s2", "lwr0.99", "upr0.99"))

  # Test: Prediction with newdata
  newdata <- data.frame(Species = levels(iris$Species),
                        Sepal.Width = c(min(iris$Sepal.Width),
                                        mean(iris$Sepal.Width),
                                        max(iris$Sepal.Width)))
  pred_newdata <- predict(mvreg_mod, newdata = newdata)
  expect_type(pred_newdata, "list")
  expect_length(pred_newdata$mu, nrow(newdata))
  expect_length(pred_newdata$log.s2, nrow(newdata))
  expect_length(pred_newdata$s2, nrow(newdata))

  # Test error handling for incorrect newdata
  incorrect_data <- data.frame(OtherColumn = 1:3)
  expect_error(suppressMessages(predict(mvreg_mod, newdata = incorrect_data)),
               "newdata must be a data.frame whose column names must be the
         same as the names of the variables in the model")
})




#### simulate() ####
test_that("simulate.mvreg generates correct response structure", {
  # Fit a model
  mvreg_mod <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris)

  # Test: Simulate with default parameters (nsim = 1, no seed)
  sim_default <- simulate(mvreg_mod)
  expect_s3_class(sim_default, "data.frame")
  expect_equal(nrow(sim_default), nrow(iris))
  expect_equal(ncol(sim_default), 1) # Default nsim is 1

  # Test: Simulate multiple response vectors (nsim = 10)
  sim_multiple <- simulate(mvreg_mod, nsim = 10)
  expect_equal(ncol(sim_multiple), 10)
  expect_equal(nrow(sim_multiple), nrow(iris))

  # Check that each column name is as expected
  expect_named(sim_multiple, paste0("sim_", seq_len(10)))
})

test_that("simulate.mvreg generates reproducible results with set seed", {
  # Fit a model
  mvreg_mod <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris)

  # Simulate responses with a fixed seed
  sim_seed_1 <- simulate(mvreg_mod, nsim = 5, seed = 123)
  sim_seed_2 <- simulate(mvreg_mod, nsim = 5, seed = 123)

  # Check if simulations with the same seed are identical
  expect_equal(sim_seed_1, sim_seed_2)

  # Check if the seed attribute is correctly set and retrieved
  expect_equal(attr(sim_seed_1, "seed")[1], 123)
  expect_equal(attr(sim_seed_2, "seed")[1], 123)
})

test_that("simulate.mvreg generates different results without setting a seed", {
  # Fit a model
  mvreg_mod <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris)

  # Simulate responses without a fixed seed
  sim_no_seed_1 <- simulate(mvreg_mod, nsim = 3)
  sim_no_seed_2 <- simulate(mvreg_mod, nsim = 3)

  # Check if simulations are different
  expect_false(identical(sim_no_seed_1, sim_no_seed_2))
})



#### update() ####
test_that("update.mvreg updates formula.mu correctly", {
  # Fit initial model
  mod <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris)

  # Update model with new predictor in formula.mu
  updated_mod <- update(mod, new.formula.mu = . ~ . + Sepal.Width)

  # Check that the formula.mu has been updated correctly
  expect_equal(deparse(updated_mod$formula.mu), deparse(Sepal.Length ~ Species + Sepal.Width))
})

test_that("update.mvreg updates formula.s2 correctly", {
  # Fit initial model
  mod <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris)

  # Update model with new predictor in formula.s2
  updated_mod <- update(mod, new.formula.s2 = . ~ . + Petal.Length)

  # Check that the formula.s2 has been updated correctly
  expect_equal(deparse(updated_mod$formula.s2), deparse(Sepal.Length ~ Sepal.Width + Petal.Length))
})

test_that("update.mvreg updates multiple formulas correctly", {
  # Fit initial model
  mod <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris)

  # Update both formula.mu and formula.s2
  updated_mod <- update(mod, new.formula.mu = . ~ . + Sepal.Width, new.formula.s2 = . ~ . + Petal.Length)

  # Check that both formulas have been updated correctly
  expect_equal(deparse(updated_mod$formula.mu), deparse(Sepal.Length ~ Species + Sepal.Width))
  expect_equal(deparse(updated_mod$formula.s2), deparse(Sepal.Length ~ Sepal.Width + Petal.Length))
})

test_that("update.mvreg allows additional parameters to be changed", {
  # Fit initial model
  mod <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris, method = "wls", vcov.type = "expected")

  # Update the method and vcov.type
  updated_mod <- update(mod, method = "full_nr", vcov.type = "observed")

  # Check that the additional parameters have been updated
  expect_equal(updated_mod$call$method, "full_nr")
  expect_equal(updated_mod$call$vcov.type, "observed")
})

test_that("update.mvreg returns the call when evaluate is FALSE", {
  # Fit initial model
  mod <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris)

  # Update model with evaluate = FALSE
  updated_call <- update(mod, new.formula.mu = . ~ . + Sepal.Width, evaluate = FALSE)

  # Check that the result is a call and not evaluated
  expect_type(updated_call, "language")
})

test_that("update.mvreg handles removing parameters with name = NULL", {
  # Fit initial model with additional parameter
  mod <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris, start.s2 = "gamma")

  # Update model and remove start.s2 parameter
  updated_mod <- update(mod, start.s2 = NULL)

  # Check that start.s2 has been removed from the updated model call
  expect_null(updated_mod$call$start.s2)
})

test_that("update.mvreg maintains consistency of object structure after update", {
  # Fit initial model
  mod <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris)

  # Update the model and check structure
  updated_mod <- update(mod, new.formula.mu = . ~ . + Sepal.Width)

  # Check that the updated model is a mvreg object with all required components
  expect_s3_class(updated_mod, "mvreg")
  expect_true(all(names(mod) == names(updated_mod)))
})


#### confint() ####

test_that("confint.mvreg computes confidence intervals for all parameters", {
  # Fit initial model
  mod <- mvreg(Sepal.Length ~ Species, ~ Sepal.Width, data = iris)

  # Compute confidence intervals
  ci <- confint(mod)

  # Check that the result is a matrix with the correct dimensions
  expect_type(ci, "double")
  expect_equal(nrow(ci), length(coef(mod)))
  expect_equal(ncol(ci), 2)

  # Check that column names are correct
  expect_equal(colnames(ci), c("lwr0.95", "upr0.95"))
})




test_that("confint.mvreg computes confidence intervals for specified parameters by name", {
  # Fit initial model
  mod <- mvreg(Sepal.Length ~ Species, ~ Sepal.Width, data = iris)

  # Compute confidence intervals for a specific parameter by name
  ci <- confint(mod, parm = "mu.Speciesversicolor")

  # Check that the result is a matrix with the correct dimensions
  expect_type(ci, "double")
  expect_equal(nrow(ci), 1)
  expect_equal(ncol(ci), 2)

  # Check that row and column names are correct
  expect_equal(rownames(ci), "mu.Speciesversicolor")
  expect_equal(colnames(ci), c("lwr0.95", "upr0.95"))
})



test_that("confint.mvreg computes confidence intervals for specified parameters by index", {
  # Fit initial model
  mod <- mvreg(Sepal.Length ~ Species, ~ Sepal.Width, data = iris)

  # Compute confidence intervals for a specific parameter by index
  ci <- confint(mod, parm = 2)

  # Check that the result is a matrix with the correct dimensions
  expect_type(ci, "double")
  expect_equal(nrow(ci), 1)
  expect_equal(ncol(ci), 2)

  # Check that row and column names are correct
  expect_equal(rownames(ci), names(coef(mod))[2])
  expect_equal(colnames(ci), c("lwr0.95", "upr0.95"))
})

test_that("confint.mvreg computes confidence intervals with different confidence levels", {
  # Fit initial model
  mod <- mvreg(Sepal.Length ~ Species, ~Sepal.Width, data = iris)

  # Compute confidence intervals with 90% confidence level
  ci_90 <- confint(mod, level = 0.90)

  # Compute confidence intervals with 99% confidence level
  ci_99 <- confint(mod, level = 0.99)

  # Check that column names reflect confidence levels
  expect_equal(colnames(ci_90), c("lwr0.9", "upr0.9"))
  expect_equal(colnames(ci_99), c("lwr0.99", "upr0.99"))

  # Ensure that intervals widen with higher confidence levels
  expect_true(all(ci_90[, "lwr0.9"] > ci_99[, "lwr0.99"]))
  expect_true(all(ci_90[, "upr0.9"] < ci_99[, "upr0.99"]))
})


#### model.matrix() ####
test_that("model.matrix.mvreg returns both design matrices when type is 'all'", {
  # Fit an mvreg model
  mod <- mvreg(Sepal.Length ~ Species + Sepal.Width, ~Species, data = iris)

  # Get both design matrices
  matrices <- model.matrix(mod, type = "all")

  # Check that it returns a list with two elements
  expect_type(matrices, "list")
  expect_named(matrices, c("xd", "zd"))

  # Check that xd and zd are matrices
  expect_true(is.matrix(matrices$xd))
  expect_true(is.matrix(matrices$zd))

  # Check dimensions of xd and zd (number of rows should match the number of observations)
  expect_equal(nrow(matrices$xd), nrow(iris))
  expect_equal(nrow(matrices$zd), nrow(iris))
})

test_that("model.matrix.mvreg returns only the mean component matrix when type is 'mu'", {
  # Fit an mvreg model
  mod <- mvreg(Sepal.Length ~ Species + Sepal.Width, ~Species, data = iris)

  # Get only the mean component design matrix
  mean_matrix <- model.matrix(mod, type = "mu")

  # Check that result is a matrix
  expect_true(is.matrix(mean_matrix))

  # Check dimensions of the mean component matrix
  expect_equal(nrow(mean_matrix), nrow(iris))
  expect_equal(ncol(mean_matrix), length(coef(mod, partition = "mu")))
})

test_that("model.matrix.mvreg returns only the variance component matrix when type is 's2'", {
  # Fit an mvreg model
  mod <- mvreg(Sepal.Length ~ Species + Sepal.Width, ~Species, data = iris)

  # Get only the variance component design matrix
  var_matrix <- model.matrix(mod, type = "s2")

  # Check that result is a matrix
  expect_true(is.matrix(var_matrix))

  # Check dimensions of the variance component matrix
  expect_equal(nrow(var_matrix), nrow(iris))
  expect_equal(ncol(var_matrix), length(coef(mod, partition = "s2")))
})

test_that("model.matrix.mvreg handles invalid type argument with an informative error", {
  # Fit an mvreg model
  mod <- mvreg(Sepal.Length ~ Species + Sepal.Width, ~Species, data = iris)

  # Check that invalid type argument triggers an error
  expect_error(model.matrix(mod, type = "invalid"), "should be one of")
})

test_that("model.matrix.mvreg defaults to type 'all' if no type specified", {
  # Fit an mvreg model
  mod <- mvreg(Sepal.Length ~ Species + Sepal.Width, ~Species, data = iris)

  # Call model.matrix without specifying type
  matrices <- model.matrix(mod)

  # Check that the default is equivalent to type = "all"
  expect_type(matrices, "list")
  expect_named(matrices, c("xd", "zd"))
  expect_true(is.matrix(matrices$xd))
  expect_true(is.matrix(matrices$zd))
})







