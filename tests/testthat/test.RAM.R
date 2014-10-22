test_that("RAM passes 1D sanity test", {
  set.seed(0)
  data <- RAM(function(x) {dnorm(x, 100,10)}, 5000, 0, 1)
  
  expect_equal(object = mean(data), expected = 100, tolerance = 0.5)
  expect_equal(object = sd(data), expected = 10, tolerance = 0.5)
})

test_that("RAM passes 2D sanity test", {
  set.seed(0)
  data <- RAM(function(x) {dmvnorm(x, c(10, 5), matrix(c(10, 2, 2, 10), 2, 2))}, 5000, rep(0, 2), diag(2))
  
  expect_equal(object = colMeans(data), expected = c(10, 5), tolerance = 0.5)
  expect_equal(object = cov(data), expected = matrix(c(10, 2, 2, 10), 2, 2), tolerance = 0.5)
})
