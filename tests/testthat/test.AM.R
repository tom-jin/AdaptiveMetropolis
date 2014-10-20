test_that("AM can find the parameters of a normal distribution", {
  data <- AM(function(x) {dnorm(x, 100,10)}, 1100, 0, 1, 100)

  expect_equal(object = mean(data[101:1100]), expected = 100, tolerance = 0.5)
  expect_equal(object = sd(data[101:1100]), expected = 10, tolerance = 0.5)
})