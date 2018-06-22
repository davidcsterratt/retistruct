context("Outline")
test_that("Outlines work correctly", {
  ## Constructor without any arguments
  a <- Outline$new()

  ## Adding points - get pids of all points added in each case
  pids <- a$addPoints(matrix(1:2, 1, 2))
  expect_equal(pids, 1)
  pids <- a$addPoints(matrix(3:4, 1, 2))
  expect_equal(pids, 2)
  pids <- a$addPoints(rbind(1:2, 3:4))
  expect_equal(pids, 1:2)

  ## getPoint() should give us the points we have entered thus far
  expect_equal(a$getPoints(), rbind(c(X=1, Y=2), c(3, 4)))

  ## Full outline
  P <- rbind(c(X=1,Y=1),
             c(2,1),
             c(2,-1),
             c(1,-1),
             c(1,-2),
             c(-1,-2),
             c(-1,-1),
             c(-2,-1),
             c(-2,1),
             c(-1,1),
             c(-1,2),
             c(1,2))
  
  o1 <- Outline$new(P)
  expect_equal(o1$getPoints(), P)
  
  ## Construct Outline with X-Y scale set
  o2 <- Outline$new(P, scale=2)
  expect_equal(o2$getPointsScaled(), P*2)

  ## But area of Triangulated Outline will be 4x bigger 

  ## FIXME (See Issue https://github.com/davidcsterratt/retistruct/issues/6)
  ## There should be more meaningful errors for outlines with crossings:
  ## > Outline$new(list(rbind(c(-4,-1), c(-1,-1), c(-3,2), c(-1,-4))))
  ## Error in P[Pt[, 1], , drop = FALSE] : subscript out of bounds
  ## > TriangulatedFragment$new(rbind(c(-4,-1), c(-1,-1), c(-3,2), c(-1,-4)))
  ## Error in out$P[1:nrow(P), ] : subscript out of bounds

  
})
