source("R/helpers.R")

testthat::test_that("smoothed approximation functions are correct", {
  x <- seq(-1, 1, 0.1)
  # Check boundaries
  testthat::expect_equal(s_gt(0, 0, 0.01), 0)
  testthat::expect_equal(s_gte(0, 0, 0.01), 1)
  testthat::expect_equal(s_lt(0, 0, 0.01), 0)
  testthat::expect_equal(s_lte(0, 0, 0.01), 1)
  
  # Check conditions
  testthat::expect_true(all(s_gt(x, 0, 1) <= ifelse(x > 0, 1, 0)))
  testthat::expect_true(all(s_lt(x, 0, 1) <= ifelse(x < 0, 1, 0)))
  testthat::expect_true(all(s_gte(x, 0, 1) >= ifelse(x >= 0, 1, 0)))
  testthat::expect_true(all(s_lte(x, 0, 1) >= ifelse(x <= 0, 1, 0)))
})
