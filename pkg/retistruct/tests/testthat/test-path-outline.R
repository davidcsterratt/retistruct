context("PathOutline")
test_that("paths with no correspondances work correctly", {

  ## Simple example
  po <- PathOutline$new()
  P <- rbind(cbind(0, 0:2), cbind(1, seq(0, 2, by=0.5)))
  colnames(P) <- c("X", "Y")
  po$addPoints(P)
  po$gf <- c(2, 3, NA, NA, 4, 5, 6, 7)
  po$gb <- c(NA, 1, 2, 5, 6, 7, 8, NA)

  po$stitchSubpaths(VF0=1, VF1=3, VB0=4, VB1=8, epsilon=0.01)
  expect_equal(po$getPoints(), rbind(P, c(0, 0.5), c(0, 1.5)))
  expect_equal(po$h, c(1, 6, 3, 4, 5, 6, 7, 8, 5, 7))
  expect_equal(po$gf, c(9, 10, NA, NA, 4, 5, 6, 7, 2, 3))
  expect_equal(po$gb, c(NA, 9, 10,  5, 6, 7, 8, NA, 1, 2))

  ## Inverse of simple example; keep same points but make forward path
  ## backward path and vice versa
  po <- PathOutline$new()
  po$addPoints(P)
  po$gb <- c(2, 3, NA, NA, 4, 5, 6, 7)
  po$gf <- c(NA, 1, 2, 5, 6, 7, 8, NA)

  po$stitchSubpaths(VF0=4, VF1=8, VB0=1, VB1=3, epsilon=0.01)

  expect_equal(po$getPoints(), rbind(P, c(0, 0.5), c(0, 1.5)))
  expect_equal(po$h, c(1, 2, 3, 4, 5, 2, 7, 8, 5, 7))
  expect_equal(po$gf, c(NA, 9, 10,  5, 6, 7, 8, NA, 1, 2))
  expect_equal(po$gb, c(9, 10, NA, NA, 4, 5, 6, 7, 2, 3))
  
})

test_that("paths with a correspondance in one path work correctly", {

  P <- rbind(c(X=1,Y=1),
             c(2,1),
             c(2.5,2),
             c(3,1),
             c(4,1),
             c(1,-1),
             c(4,-1))
  po <- PathOutline$new()
  po$addPoints(P)
  po$gf <- c(2, 3, 4, 5, NA, NA, 6)
  po$gb <- c(NA, 1, 2, 3, 4, 7, NA)
  po$hf[2] <- 4
  po$hb[4] <- 2

  po$stitchSubpaths(VF0=1, VF1=5, VB0=6, VB1=7, epsilon=0.01)
  expect_equal(po$getPoints(), rbind(P, c(2.5, -1)))
  expect_equal(po$h, c(1, 2, 3, 2, 5, 6, 7, 2))
  expect_equal(po$gf, c(2, 3, 4, 5, NA, NA, 8, 6))
  expect_equal(po$gb, c(NA, 1, 2, 3, 4, 8, NA, 7))

  ## Add extra points in the backwards path
  P <- rbind(c(X=1,Y=1),
             c(2,1),
             c(2.5,2),
             c(3,1),
             c(4,1),
             c(1,-1),
             c(2,-1),
             c(3,-1),
             c(4,-1))

  po <- PathOutline$new()
  po$addPoints(P)
  po$gf <- c(2, 3, 4, 5, NA, NA, 6, 7, 8)
  po$gb <- c(NA, 1, 2, 3, 4, 7, 8, 9, NA)
  po$hf[2] <- 4
  po$hb[4] <- 2
  po$stitchSubpaths(VF0=1, VF1=5, VB0=6, VB1=9, epsilon=0.01)
  expect_equal(po$getPoints(), rbind(P, c(5/3, 1), c(2.5, -1), c(10/3, 1)))
  expect_equal(po$h, c(1, 2, 3, 2, 5, 6, 7, 8, 9, 7, 2, 8))
  expect_equal(po$gf, c(10, 3, 4, 12, NA, NA, 6, 11, 8, 2, 7, 5))
  expect_equal(po$gb, c(NA, 10, 2, 3, 12, 7, 11, 9, NA, 1, 8, 4))
  
})

test_that("paths with correspondances in each path work correctly", {

  P <- rbind(c(X=1,Y=1),
             c(2,1),
             c(2.5,2),
             c(3,1),
             c(4,1),
             c(1,-1),
             c(2,-1),
             c(2.5,-2),
             c(3,-1),
             c(4,-1))

  po <- PathOutline$new()
  po$addPoints(P)
  po$gf <- c(2, 3, 4, 5, NA, NA, 6, 7, 8, 9)
  po$gb <- c(NA, 1, 2, 3, 4, 7, 8, 9, 10, NA)
  po$hf[2] <- 4
  po$hb[4] <- 2
  po$hf[9] <- 7
  po$hb[7] <- 9
  
  po$stitchSubpaths(VF0=1, VF1=5, VB0=6, VB1=10, epsilon=0.01)
  expect_equal(po$getPoints(), P)
  expect_equal(po$h, c(1, 7, 3, 7, 5, 6, 7, 8, 7, 10))
  expect_equal(po$gf, c(2, 3, 4, 5, NA, NA, 6, 7, 8, 9))
  expect_equal(po$gb, c(NA, 1, 2, 3, 4, 7, 8, 9, 10, NA))
})

test_that("paths with a double correspondance in one path works correctly", {
  P <- rbind(c(X=1,Y=1),
             c(2,1),
             c(2.5,2),
             c(3,1),
             c(4,1),
             c(1,-1),
             c(2,-1),
             c(2.2,-2),
             c(2.5,-1),
             c(2.8,-2),
             c(3,-1),
             c(4,-1))
  
  po <- PathOutline$new()
  po$addPoints(P)
  gf <- c(2, 3, 4, 5, NA, NA, 6, 7, 8, 9, 10, 11)
  po$gf <- gf
  gb <- c(NA, 1, 2, 3, 4, 7, 8, 9, 10, 11, 12, NA)
  po$gb <- gb
  po$hf[2] <- 4
  po$hb[4] <- 2
  po$hf[9] <- 7
  po$hb[7] <- 9
  po$hf[11] <- 9
  po$hb[9] <- 11
  
  po$stitchSubpaths(VF0=1, VF1=5, VB0=6, VB1=12, epsilon=0.01)
  expect_equal(po$getPoints(), P)
  expect_equal(po$h, c(1, 7, 3, 7, 5, 6, 7, 8, 7, 10, 7, 12))
  expect_equal(po$gf, gf)
  expect_equal(po$gb, gb)
})
