context("CountSet")
test_that("CountSet works correctly", {
  data <- cbind(X=1, Y=2, C=3)
  cs <- CountSet$new(list(test=data), c("red"))
  expect_that(cs$getIDs(), equals("test"))
  expect_that(cs$getFeature("test"), equals(data))
  expect_error(CountSet$new(list(cbind(X=1, Y=2, C=3, D=4)), c("red")))
  expect_error(CountSet$new(list(cbind(X=1, Y=2)), c("red")))
})
