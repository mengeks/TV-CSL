library(testthat)


trained_model <- lm(mpg ~ wt + factor(am), data = mtcars)
X_example <- mtcars[, c("wt", "am")] 
W_example <- 1  


test_that("get_counterfactual_estimation returns numeric vector of correct length", {
  result <- 
    get_counterfactual_estimation(trained_model, W_example, X_example)
  expect_type(result, "double") 
  expect_equal(length(result), nrow(X_example))  # The length should match the number of observations in X
})

# Test 2: Check if the function handles different binary inputs for W
test_that("get_counterfactual_estimation handles different binary values for W", {
  result_treated <- get_counterfactual_estimation(trained_model, W = 1, X = X_example)
  result_control <- get_counterfactual_estimation(trained_model, W = 0, X = X_example)
  
  expect_type(result_treated, "double")
  expect_type(result_control, "double")
  expect_false(identical(result_treated, result_control))  # The results for treated and control should differ
})