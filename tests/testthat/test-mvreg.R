test_that("mvreg() returns a mvreg class object", {
  expect_equal(class(mvreg(Sepal.Length ~ Species, data = iris)), "mvreg")
})

