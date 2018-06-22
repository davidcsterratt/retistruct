context("StitchedOutline")
test_that("StitchedOutlines work correctly", {
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

  ## Stitch tears without triangulation
  a$stitchTears()
  ## hf and hb points point accross tears
  expect_equal(which(a$hf != 1:length(a$hf)), c(2, 5, 8, 11))
  expect_equal(which(a$hb != 1:length(a$hb)), c(3, 6, 9, 12))

  ## Test a path around the rim
  expect_equal(path(5, 12, g=a$gf, h=a$hf), c(5, 3, 2,  12))
  expect_equal(nrow(a$P), length(a$gf))

  ## Stitch tears after triangulation
  b <- StitchedOutline$new(P)
  b$addTear(c(3, 4, 5))
  b$addTear(c(6, 7, 8))
  b$addTear(c(9, 10, 11))
  b$addTear(c(12, 1, 2))
  
  ## Triangulate
  b$triangulate()
  expect_equal(nrow(b$P), length(b$gf))
  b$stitchTears()
  b$triangulate(suppress.external.steiner=TRUE)
  expect_equal(length(a$L), nrow(a$Cu))
  b$stitchTears()
})
