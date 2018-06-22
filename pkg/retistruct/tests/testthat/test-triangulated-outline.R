context("TriangulatedOutline")
test_that("TriangulatedOutlines work correctly", {
  ## 2x2 square - area 4
  P <- rbind(c(0,0),
             c(2,0),
             c(2,2),
             c(0,2))
  ## Check area OK
  o <- TriangulatedOutline$new(P)
  o$triangulate()
  expect_equal(sum(o$A), 4)
  expect_true(all(is.numeric(o$L)))
  expect_equal(nrow(o$getPoints()), length(o$gf))
  expect_equal(nrow(o$getPoints()), length(o$gb))
  expect_equal(sum(o$getOutlineLengths()), 8)
  
  ## Check area OK with xy scale
  o1 <- TriangulatedOutline$new(P, scale=2)
  o1$triangulate()
  expect_equal(sum(o1$A), 16)
  expect_equal(o1$A.tot, 16)
  expect_equal(sum(o1$L), 2*sum(o$L))
  expect_equal(nrow(o1$getPoints()), length(o1$gf))
  expect_equal(nrow(o1$getPoints()), length(o1$gb))
  expect_equal(sum(o1$getOutlineLengths()), 16)
})
