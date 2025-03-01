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

test_that("AnnotatedOutlines comprising fragments with tears work correctly", {
  ## Constructing multi-fragment outline
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

  ## Annotated outlines
  a <- AnnotatedOutline$new(P)

  ## Check triangulation works
  f <- TriangulatedFragment$new(a$getFragment(1), n=NA)
  expect_equal(f$A.tot, 4)

  # If we ask for a non-existent tear, get NA
  expect_equal(NA, a$getTear(1))

  ## Add a tear

  ## Points across two fragments won't work
  expect_error(a$addTear(c(14, 15, 16)))

  a$addTear(c(2, 3, 4))
  tear <- a$getTear(1)
  expect_equal(tear, c(V0=3, VB=2, VF=4))

  ## Overlapping tears throw an error
  expect_error(a$addTear(c(15, 16, 17)))

  ## Remove a tear
  a$removeTear(1)
  expect_equal(NA, a$getTear(1))

  ## Check whichTear()
  a$addTear(c(2, 3, 4))
  a$addTear(c(17, 18, 19))

  expect_equal(a$whichTear(2), 1)
  expect_equal(a$whichTear(18), 2)

  ## Check checkTears()
  expect_equal(length(a$checkTears()), 0)

  ## Check we can add a cut

  ## Wrong number of points won't work
  expect_error(a$addFullCut(c(76, 35, 44, 96, 67)))

  ## Non-existent points won't work
  expect_error(a$addFullCut(c(76, 35, 37, 44)))

  ## 3 points in one fragment won't work
  expect_error(a$addFullCut(a, c(1, 2, 3, 9)))

  ##  Points in more than two fragments won't work
  expect_error(a$addFullCut(c(1, 7, 10, 15)))

  ## Add a cut
  a$addFullCut(c(1, 5, 16, 20))
  expect_true(setequal(a$getFullCuts(), cbind(1, 5, 16, 20)))

  ## Check
  cr <- a$computeFullCutRelationships(a$fullcuts)
  expect_true(setequal(cr$Rset, c(1, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 20, 21)))

  ## Set a fixed point
  ## One that is in the rim should be fine
  a$setFixedPoint(5, "Nasal")
  expect_equal(a$i0, c(Nasal=5))

  ## One that is not in the rim should be moved
  a$setFixedPoint(2, "Nasal")
  expect_false(a$i0 == 2)

  ## Add a tear outside a cut
  a$addTear(c(9, 10, 11))
  tr <- a$computeTearRelationships(a$getTears())
  expect_true(setequal(tr$Rset, (1:21)[-c(3, 10, 18)]))

  ## A fixed point in a tear should be moved
  a$setFixedPoint(10, "Nasal")
  expect_false(a$i0 == 10)

  ## Add more fullcuts
  a$addFullCut(c(1, 7, 6, 8))
  a$addFullCut(c(7, 14, 12, 13))
  a$addFullCut(c(14, 15, 16, 21))

  expect_true(setequal(a$getRimSet(), c(5, 6, 8, 9, 11, 12, 13, 15, 21, 20)))
})

test_that("AnnotatedOutlines comprising fragments without tears work correctly", {
    P <- list(rbind(c(1,1),
                  c(4,1),
                  c(1,4)),
            rbind(c(-1,1),
                  c(-1,4),
                  c(-4,1)),
            rbind(c(-4,-1),
                  c(-1,-1),
                  c(-1,-4)),
            rbind(c(1,-1),
                  c(4,-1),
                  c(1,-4)))

    ## Annotated outlines
    a <- AnnotatedOutline$new(P)

    ## Add fullcuts
    a$addFullCut(c(1, 10, 11, 2))
    a$addFullCut(c(1, 4, 3, 5))
    a$addFullCut(c(4, 8, 6, 7))
    a$addFullCut(c(8, 10, 9, 12))
})
