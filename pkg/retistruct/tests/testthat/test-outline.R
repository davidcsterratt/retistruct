context("Outline")
test_that("Outlines work correctly", {
  ## Constructor without any arguments
  a <- Outline$new()

  ## Adding points - get pids of all points added in each case
  pids <- a$addPoints(matrix(1:2, 1, 2), 1)
  expect_equal(pids, 1)
  pids <- a$addPoints(matrix(3:4, 1, 2), 1)
  expect_equal(pids, 2)
  pids <- a$addPoints(rbind(1:2, 3:4), 1)
  expect_equal(pids, 1:2)

  ## getFragmementPointIDs() should give the IDs of points in fragment 1
  expect_equal(a$getFragmentPointIDs(1), 1:2)
  
  ## Get fragment should give us the points we have entered thus far
  expect_equal(a$getFragmentPoints(1), rbind(c(X=1, Y=2), c(3, 4)))

  f1 <- a$getFragment(1)
  expect_equal(f1$P, rbind(c(X=1, Y=2), c(3, 4)))
  
  ## Full outline
  P <- list(rbind(c(X=1,Y=1),
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
             c(1,2)))
  
  o1 <- Outline$new(P)
  expect_equal(o1$getPoints(), P[[1]])
  
  ## Construct Outline with X-Y scale set
  o2 <- Outline$new(P, scale=2)
  expect_equal(o2$getPointsScaled(), P[[1]]*2)

  ## Constructing multi-fragment outline
  P <- list(rbind(c(1,1),
                  c(2,1),
                  c(2.5,2),
                  c(3,1),
                  c(4,1),
                  c(1,4)),
            rbind(c(-1,1),
                  c(-1,4),
                  c(-4,1)),
            rbind(c(-4,-1),
                  c(-1,-1),
                  c(-1,-4)),
            rbind(c(1,-1),
                  c(2,-1),
                  c(2.5,-2),
                  c(3,-1),
                  c(4,-1),
                  c(1,-4)))

  o1 <- Outline$new(P)
  nP <- sum(sapply(P, nrow))
  expect_equal(nrow(o1$P), nP)
  expect_equal(o1$h, 1:nP)

  ## getFragmementPointIDs() should give the IDs of points in fragment 1
  expect_equal(o1$getFragmentPointIDs(4), 13:18)

  f4 <- o1$getFragment(4)
  expect_equal(f4$P, cbind(rbind(c(X=1,Y=-1),
                                 c(2,-1),
                                 c(2.5,-2),
                                 c(3,-1),
                                 c(4,-1),
                                 c(1,-4))))
  
  expect_equal(f4$gf, c(6, 1:5))
  expect_equal(f4$gb, c(2:6, 1))
  expect_equal(f4$h, 1:6)
  
  f1 <- TriangulatedFragment$new(o1$getFragment(1), n=NA)
  expect_equal(f1$A.tot, 4)

  ## Construct multi-fragment outline with X-Y scale set
  o2 <- Outline$new(P, scale=2)
  f2 <- TriangulatedFragment$new(o2$getFragment(1), n=NA)
  ## Area of Fragment should still be 4
  expect_equal(f2$A.tot, 4)

  ## But area of Triangulated Outline will be 4x bigger 

  ## FIXME (See Issue https://github.com/davidcsterratt/retistruct/issues/6)
  ## There should be more meaningful errors for outlines with crossings:
  ## > Outline$new(list(rbind(c(-4,-1), c(-1,-1), c(-3,2), c(-1,-4))))
  ## Error in P[Pt[, 1], , drop = FALSE] : subscript out of bounds
  ## > TriangulatedFragment$new(rbind(c(-4,-1), c(-1,-1), c(-3,2), c(-1,-4)))
  ## Error in out$P[1:nrow(P), ] : subscript out of bounds
  

  
})
