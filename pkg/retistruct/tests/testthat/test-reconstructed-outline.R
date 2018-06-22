context("ReconstructedOutline")
test_that("ReconstructededOutlines work correctly", {
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
  r <- ReconstructedOutline$new(a)
  r$mergePointsEdges()
  expect_equal(length(r$Lt), nrow(r$Cut))
  expect_true(!any(is.na(r$Lt)))
  r$projectToSphere()
  expect_equal(r$ol$A.tot, 4 + 4*2)
  expect_true(is.numeric(r$R))
  r$getStrains()

  #r <- ReconstructedOutline$new(a)
  #r$reconstruct()

  ## FIXME: Test of scaling required. i.e. suppose we scale the
  ## original points (P) by a factor, do we get the same mapping?
})

