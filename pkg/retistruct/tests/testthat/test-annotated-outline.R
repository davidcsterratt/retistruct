

context("AnnotatedOutline")
test_that("AnnotatedOutlines with tears work correctly", {
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

  ## Annotated outlines
  a <- AnnotatedOutline$new(P)
  expect_equal(sum(a$getRimLengths()), 16)
  
  ## Check triangulation works
  f <- TriangulatedFragment$new(a, n=NA)
  expect_equal(f$A.tot, 4 + 4*2)
  expect_equal(sum(a$getRimLengths()), 16)
  
  # If we ask for a non-existent tear, get NA
  expect_equal(NA, a$getTear(1))
  
  ## Add a tear
  a$addTear(c(3, 4, 5))
  tear <- a$getTear(1)
  expect_equal(tear, c(V0=4, VB=5, VF=3))

  ## Overlapping tears throw an error
  expect_error(a$addTear(c(4, 5, 6)))
  
  ## Remove a tear
  a$removeTear(1)
  expect_equal(NA, a$getTear(1))

  ## Check whichTear()
  a$addTear(c(3, 4, 5))
  expect_equal(a$whichTear(3), 1)

  ## Check checkTears()
  expect_equal(length(a$checkTears()), 0)
  
  ## Set a fixed point
  ## One that is in the rim should be fine
  a$setFixedPoint(5, "Nasal")
  expect_equal(a$i0, c(Nasal=5))
  
  ## A fixed point in a tear should be moved
  a$setFixedPoint(4, "Nasal")
  expect_false(a$i0 == 4)

  a$addTear(c(6, 7, 8))
  a$addTear(c(9, 10, 11))
  a$addTear(c(12, 1, 2))
  expect_true(setequal(a$getRimSet(), c(2, 3, 5, 6, 8, 9, 11, 12)))
})


