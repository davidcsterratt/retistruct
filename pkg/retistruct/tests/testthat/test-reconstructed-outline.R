context("ReconstructedOutline")
test_that("ReconstructededOutlines with a single fragment work correctly", {
  P <- rbind(c(1,1),
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

  ## Stitched outlines
  a <- StitchedOutline$new(P)

  ## Set a fixed point
  ## One that is in the rim should be fine
  a$setFixedPoint(5, "Nasal")
  expect_equal(a$i0, c(Nasal=5))

  ## One that is not in the rim should be moved
  a$addTear(c(3, 4, 5))
  a$addTear(c(6, 7, 8))
  a$addTear(c(9, 10, 11))
  a$addTear(c(12, 1, 2))

  ## Reconstruct
  r <- ReconstructedOutline$new()
  r$loadOutline(a)
  r$mergePointsEdges()
  expect_equal(length(r$Lt), nrow(r$Cut))
  expect_true(!any(is.na(r$Lt)))
  r$projectToSphere()
  expect_equal(r$ol$A.tot, 4 + 4*2)
  expect_true(is.numeric(r$R))
  r$getStrains()

  ## FIXME: Test of scaling required. i.e. suppose we scale the
  ## original points (P) by a factor, do we get the same mapping?
})

test_that("ReconstructededOutlines with mutliple fragments work correctly", {
  P <- list(rbind(c(1,1),
                  c(2,1),
                  c(2.5,2),
                  c(3,1),
                  c(4,1),
                  c(1,4)),
            rbind(c(-1,1),
                  c(-1,4),
                  c(-2,3),
                  c(-2,2),
                  c(-3,2),
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

  ## Stitched outlines
  a <- StitchedOutline$new(P)

  ## Set a fixed point
  ## One that is in the rim should be fine
  a$setFixedPoint(5, "Nasal")
  expect_equal(a$i0, c(Nasal=5))

  ## One that is not in the rim should be moved
  a$addTear(c(9, 10, 11))
  a$addTear(c(2, 3, 4))
  a$addTear(c(17, 18, 19))

  ## Add fullcuts
  a$addFullCut(c(1, 5, 16, 20))
  a$addFullCut(c(1, 7, 6, 8))
  a$addFullCut(c(7, 14, 12, 13))
  a$addFullCut(c(14, 15, 16, 21))
  r <- ReconstructedOutline$new()
  r$loadOutline(a)
  r$mergePointsEdges()
  expect_equal(length(r$Lt), nrow(r$Cut))
  expect_true(!any(is.na(r$Lt)))
  r$projectToSphere()
  expect_equal(r$ol$A.tot, 16.5)
  expect_true(is.numeric(r$R))
  r$getStrains()

  #r <- ReconstructedOutline$new(a)
  #r$reconstruct()

  ## FIXME: Test of scaling required. i.e. suppose we scale the
  ## original points (P) by a factor, do we get the same mapping?
})

test_that("ReconstructededOutlines with multiple fragments with a hole work correctly", {
  ## Constructing multi-fragment outline
  P <- list(rbind(c(1,1.5),
                  c(1.5,1),
                  c(2,1),
                  c(2.5,2),
                  c(3,1),
                  c(4,1),
                  c(1,4)),
            rbind(c(-1.5,1),
                  c(-1,1.5),
                  c(-1,4),
                  c(-2,3),
                  c(-2,2),
                  c(-3,2),
                  c(-4,1)),
            rbind(c(-4,-1),
                  c(-1.5,-1),
                  c(-1,-1.5),
                  c(-1,-4)),
            rbind(c(1,-1.5),
                  c(1.5,-1),
                  c(2,-1),
                  c(2.5,-2),
                  c(3,-1),
                  c(4,-1),
                  c(1,-4)))
  ## Stitched outlines
  a <- StitchedOutline$new(P)

  expect_false(a$isStitched())

  ## Set a fixed point
  ## One that is in the rim should be fine
  a$setFixedPoint(6, "Nasal")
  expect_equal(a$i0, c(Nasal=6))

  ## One that is not in the rim should be moved
  a$addTear(c(11, 12, 13))
  a$addTear(c(3, 4, 5))
  a$addTear(c(21, 22, 23))

  ## Add fullcuts
  a$addFullCut(c(2, 6, 20, 24))
  a$addFullCut(c(1, 7, 9, 10))
  a$addFullCut(c(8, 14, 15, 16))
  a$addFullCut(c(17, 18, 19, 25))

  r <- ReconstructedOutline$new()
  r$loadOutline(a)
  r$mergePointsEdges()
  expect_equal(length(r$Lt), nrow(r$Cut))
  expect_true(!any(is.na(r$Lt)))
  r$projectToSphere()
  expect_equal(r$ol$A.tot, 16.0)
  expect_true(is.numeric(r$R))
  r$getStrains()
  r$reconstruct()
})
