context("TriangulatedOutline")
test_that("TriangulatedOutlines work correctly", {
  ## 2x2 square - area 4
  P <- list(rbind(c(0,0),
                  c(2,0),
                  c(2,2),
                  c(0,2)))
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

    ## Check area OK with xy scale and z scale
  o2 <- TriangulatedOutline$new(P, scale=2, scalez=1)
  o2$triangulate()
  expect_equal(sum(o2$A), 16)
  expect_equal(o2$A.tot, 16)
  expect_equal(nrow(o2$getPoints()), length(o2$gf))
  expect_equal(nrow(o2$getPoints()), length(o2$gb))


  ## 3, 4, 5 triangle in y-z direction

  dm <- rbind(c(0, 0),
              c(4/3, 4/3),
              c(8/3, 8/3),
              c(4, 4))

  P <- list(rbind(c(0.5, 0.5),
                  c(1.5, 0.5),
                  c(1.5, 3.5),
                  c(0.5, 3.5)))

  ## Check area OK with xy scale and z scale
  o3 <- TriangulatedOutline$new(P, dm=dm, scale=1, scalez=1)
  o3$triangulate()
  expect_equal(sum(o3$A), 1*5)
  expect_equal(o3$A.tot, 1*5)
  expect_true(all(is.numeric(o3$L)))
  expect_equal(nrow(o3$getPoints()), length(o3$gf))
  expect_equal(nrow(o3$getPoints()), length(o3$gb))

  ## Double-check without depthmap
  o4 <- TriangulatedOutline$new(P, scale=1, scalez=1)
  o4$triangulate()
  expect_equal(sum(o4$A), 1*3)
  expect_equal(o4$A.tot, 1*3)
  expect_true(sum(o3$L) > sum(o4$L))
  expect_equal(nrow(o4$getPoints()), length(o4$gf))
  expect_equal(nrow(o4$getPoints()), length(o4$gb))
})
