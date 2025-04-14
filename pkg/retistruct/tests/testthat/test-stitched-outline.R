context("StitchedOutline")
test_that("StitchedOutlines with a single fragment work correctly", {
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
  expect_false(a$isStitched())
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

test_that("StitchedOutlines with multiple fragments work correctly", {
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

  ## hf and hb point across on rim
  cr <- a$computeFullCutRelationships(a$fullcuts)
  expect_equal(cr$hf[6], 8)
  expect_equal(cr$hf[12], 13)
  expect_equal(cr$hf[15], 21)
  expect_equal(cr$hf[20], 5)
  expect_equal(cr$hb[8],  6)
  expect_equal(cr$hb[13], 12)
  expect_equal(cr$hb[21], 15)
  expect_equal(cr$hb[5],  20)

  ## Stitch tears
  expect_equal(nrow(a$P), length(a$gf))
  a$triangulate()
  expect_equal(nrow(a$P), length(a$gf))
  a$stitchTears()
  expect_equal(nrow(a$P), length(a$gf))
  expect_equal(which(a$hf != 1:length(a$hf)), c(2, 9, 19))
  expect_equal(which(a$hb != 1:length(a$hb)), c(4, 11, 17))
  ## Points in tears should not be in rim set
  expect_false(any(c(3, 33, 28, 223, 18, 230, 98, 10, 83) %in% a$Rset))

  cr <- a$computeFullCutRelationships(a$fullcuts)
  expect_true(setequal(cr$CFset[[1]],
                       c(1, 45, 2, 4, 63, 5)))

  expect_true(setequal(cr$CBset[[1]],
                       c(16, 195, 17, 19, 205, 20)))

  ## hf and hb point across on rim
  expect_equal(cr$hf[6], 8)
  expect_equal(cr$hf[12], 13)
  expect_equal(cr$hf[15], 21)
  expect_equal(cr$hf[20], 5)
  expect_equal(cr$hb[8],  6)
  expect_equal(cr$hb[13], 12)
  expect_equal(cr$hb[21], 15)
  expect_equal(cr$hb[5],  20)

  ## hf and hb points point accross tears
  expect_equal(which(a$hf != 1:length(a$hf)), c(2, 9, 19))
  expect_equal(which(a$hb != 1:length(a$hb)), c(4, 11, 17))

  a$triangulate(suppress.external.steiner=TRUE)
  ## Points in tears should not be in rim set
  expect_false(any(c(3, 33, 28, 223, 18, 230, 98, 10, 83) %in% a$Rset))
  expect_equal(nrow(a$P), length(a$gf))
  expect_equal(nrow(a$P), length(a$hf))
  expect_equal(nrow(a$P), length(a$h))
  expect_equal(nrow(a$P), length(a$hb))

  ## hf and hb points point accross tears
  expect_equal(which(a$hf != 1:length(a$hf)), c(2, 9, 19))
  expect_equal(which(a$hb != 1:length(a$hb)), c(4, 11, 17))

  ## Stitch correspondences
  a$stitchFullCuts()
  ## Points in tears should not be in rim set
  expect_false(any(c(3, 33, 28, 223, 18, 230, 98, 10, 83) %in% a$Rset))

  trueRset <- c(21, 218, 190, 198, 194, 178, 183, 181, 179, 204, 186, 207, 20,
               5, 64, 66, 24, 26, 27, 22, 34, 32, 51, 6,
               8, 112, 114, 9, 11, 90, 89, 100, 12,
               13, 172, 148, 165, 123, 130, 125, 126, 15)
  expect_true(setequal(a$Rset, trueRset))
  expect_true(setequal(a$getRimSet(), trueRset))

  expect_equal(a$h[146], 213)
  ## Test a path around the rim
  expect_equal(path(204, 66, g=a$gf, h=a$h), c(204, 186, 207,  20, 5, 64, 66))
  path(a$Rset[1], a$Rset[2], g=a$gf, a$hf)
  path(a$Rset[2], a$Rset[1], g=a$gf, a$h)

  a$triangulate(suppress.external.steiner=TRUE)

  expect_equal(length(a$L), nrow(a$Cu))

  path(a$Rset[1], a$Rset[2], g=a$gf, a$hf)
  path(a$Rset[2], a$Rset[1], g=a$gf, a$h)
})

test_that("StitchedOutlines with multiple fragments with a hole work correctly", {
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

  ## hf and hb point across gaps rim
  ## hf points anti-clockwise and hb points clockwise
  cr <- a$computeFullCutRelationships(a$fullcuts)
  expect_equal(cr$hf[7], 10)
  expect_equal(cr$hf[14], 15)
  expect_equal(cr$hf[18], 25)
  expect_equal(cr$hf[24], 6)
  expect_equal(cr$hb[10],  7)
  expect_equal(cr$hb[15], 14)
  expect_equal(cr$hb[25], 18)
  expect_equal(cr$hb[6],  24)

  ## hf and hb point across gaps on the internal boundary
  ## hf points clockwise and hb points anti-clockwise
  expect_equal(cr$hf[2] , 20)
  expect_equal(cr$hf[19], 17)
  expect_equal(cr$hf[16],  8)
  expect_equal(cr$hf[9] ,  1)
  expect_equal(cr$hb[1],   9)
  expect_equal(cr$hb[8],  16)
  expect_equal(cr$hb[17], 19)
  expect_equal(cr$hb[20],  2)

  ## Stitch tears
  expect_equal(nrow(a$P), length(a$gf))
  a$triangulate()
  expect_equal(nrow(a$P), length(a$gf))
  a$stitchTears()
  expect_equal(nrow(a$P), length(a$gf))
  expect_equal(which(a$hf != 1:length(a$hf)), c(3, 11, 23))
  expect_equal(which(a$hb != 1:length(a$hb)), c(5, 13, 21))
  ## Points in tears should not be in rim set
  ## expect_false(any(c(3, 33, 28, 223, 18, 230, 98, 10, 83) %in% a$Rset))

  cr <- a$computeFullCutRelationships(a$fullcuts)
  expect_true(setequal(cr$CFset[[1]],
                       c(2, 3, 5, 53, 6)))

  expect_true(setequal(cr$CBset[[1]],
                       c(20, 21, 23, 222, 24)))

  ## hf and hb point across gaps rim
  ## hf points anti-clockwise and hb points clockwise
  cr <- a$computeFullCutRelationships(a$fullcuts)
  expect_equal(cr$hf[7], 10)
  expect_equal(cr$hf[14], 15)
  expect_equal(cr$hf[18], 25)
  expect_equal(cr$hf[24], 6)
  expect_equal(cr$hb[10],  7)
  expect_equal(cr$hb[15], 14)
  expect_equal(cr$hb[25], 18)
  expect_equal(cr$hb[6],  24)

  ## hf and hb point across gaps on the internal boundary
  ## hf points clockwise and hb points anti-clockwise
  expect_equal(cr$hf[2] , 20)
  expect_equal(cr$hf[19], 17)
  expect_equal(cr$hf[16],  8)
  expect_equal(cr$hf[9] ,  1)
  expect_equal(cr$hb[1],   9)
  expect_equal(cr$hb[8],  16)
  expect_equal(cr$hb[17], 19)
  expect_equal(cr$hb[20],  2)

  a$triangulate(suppress.external.steiner=TRUE)

  ## Points in tears should not be in rim set
  expect_false(any(c(4, 29, 42, 12, 90, 109, 22, 199, 233, 188, 186) %in% a$Rset))
  expect_equal(nrow(a$P), length(a$gf))
  expect_equal(nrow(a$P), length(a$hf))
  expect_equal(nrow(a$P), length(a$h))
  expect_equal(nrow(a$P), length(a$hb))

  ## hf and hb points point accross tears
  expect_equal(which(a$hf != 1:length(a$hf)), c(3, 11, 23))
  expect_equal(which(a$hb != 1:length(a$hb)), c(5, 13, 21))

  ## Stitch correspondences
  a$stitchFullCuts()
  expect_true(a$isStitched())

  bs <- a$getBoundarySets()
  expect_equal(length(bs), 2)
  expect_true("Rim" %in% names(bs))
  expect_true(is.vector(bs[["Rim"]]))
  expect_true(is.vector(bs[["n1"]]))

  ## The points on the hole should be ordered in the direction of the forward pointer gf
  expect_equal(bs[["n1"]], c(20, 196, 19, 17, 125, 16, 8, 77, 9, 1, 2))

  ## Points

  ## Points in tears or fullcuts should not be in rim set
  # expect_false(any(c(3, 33, 28, 223, 18, 230, 98, 10, 83) %in% a$Rset))

  trueRset <- c(56,28,41,26,32,46,30,69,7,10,103,93,104,11,13,115,113,117,14,
                15,172,147,165,123,151,139,161,150,159,18,25,208,184,202,187,181,198,183,
                221,219,224,24,6,54)

  ## FIXME?? a$Rset isn't the rim set any more...
  ## expect_true(setequal(a$Rset, trueRset))
  expect_true(setequal(a$getRimSet(), trueRset))

  expect_equal(a$h[24], 6)
  ## Test a path around the rim
  expect_equal(path(219, 54, g=a$gf, h=a$h), c(219, 224, 24, 6, 54))
  rs <- a$getRimSet()
  path(rs[1], rs[2], g=a$gf, a$hf)
  path(rs[2], rs[1], g=a$gf, a$h)

  a$triangulate(suppress.external.steiner=TRUE)

  expect_equal(length(a$L), nrow(a$Cu))

  rs <- a$getRimSet()
  path(rs[1], rs[2], g=a$gf, a$hf)
  path(rs[2], rs[1], g=a$gf, a$h)
})
