context("paths")

## test_that("path length works correctly", {
## })
 
test_that("paths with no correspondances work correctly", {

  ## Simple example
  P <- rbind(cbind(0, 0:2), cbind(1, seq(0, 2, by=0.5)))
  gf <- c(2, 3, NA, NA, 4, 5, 6, 7)
  gb <- c(NA, 1, 2, 5, 6, 7, 8, NA)
  hf <- 1:nrow(P)
  hb <- 1:nrow(P)
  out <- stitch.paths(P, VF0=1, VF1=3, VB0=4, VB1=8, gf, gb, hf, hb, h=hf, epsilon=0.01)
  expect_equal(out$P, rbind(P, c(0, 0.5), c(0, 1.5)))
  expect_equal(out$h, c(1, 6, 3, 4, 5, 6, 7, 8, 5, 7))
  expect_equal(out$gf, c(9, 10, NA, NA, 4, 5, 6, 7, 2, 3))
  expect_equal(out$gb, c(NA, 9, 10,  5, 6, 7, 8, NA, 1, 2))

  ## Inverse of simple example; keep same points but make forward path
  ## backward path and vice versa
  gb <- c(2, 3, NA, NA, 4, 5, 6, 7)
  gf <- c(NA, 1, 2, 5, 6, 7, 8, NA)
  hf <- 1:nrow(P)
  hb <- 1:nrow(P)
  out <- stitch.paths(P, VF0=4, VF1=8, VB0=1, VB1=3, gf, gb, hf, hb, h=hf, epsilon=0.01)

  expect_equal(out$P, rbind(P, c(0, 0.5), c(0, 1.5)))
  expect_equal(out$h, c(1, 2, 3, 4, 5, 2, 7, 8, 5, 7))
  expect_equal(out$gf, c(NA, 9, 10,  5, 6, 7, 8, NA, 1, 2))
  expect_equal(out$gb, c(9, 10, NA, NA, 4, 5, 6, 7, 2, 3))
  
})

test_that("paths with a correspondance in one path work correctly", {

  P <- rbind(c(1,1),
             c(2,1),
             c(2.5,2),
             c(3,1),
             c(4,1),
             c(1,-1),
             c(4,-1))

  gf <- c(2, 3, 4, 5, NA, NA, 6)
  gb <- c(NA, 1, 2, 3, 4, 7, NA)
  h <- 1:nrow(P)
  hf <- h
  hb <- h
  hf[2] <- 4
  hb[4] <- 2

  out <- stitch.paths(P, VF0=1, VF1=5, VB0=6, VB1=7, gf, gb, hf, hb, h=hf, epsilon=0.01)
  expect_equal(out$P, rbind(P, c(2.5, -1)))
  expect_equal(out$h, c(1, 4, 3, 4, 5, 6, 7, 2))
  expect_equal(out$gf, c(2, 3, 4, 5, NA, NA, 8, 6))
  expect_equal(out$gb, c(NA, 1, 2, 3, 4, 8, NA, 7))

  ## Add extra points in the backwards path
  P <- rbind(c(1,1),
             c(2,1),
             c(2.5,2),
             c(3,1),
             c(4,1),
             c(1,-1),
             c(2,-1),
             c(3,-1),
             c(4,-1))

  gf <- c(2, 3, 4, 5, NA, NA, 6, 7, 8)
  gb <- c(NA, 1, 2, 3, 4, 7, 8, 9, NA)
  h <- 1:nrow(P)
  hf <- h
  hb <- h
  hf[2] <- 4
  hb[4] <- 2
  out <- stitch.paths(P, VF0=1, VF1=5, VB0=6, VB1=9, gf, gb, hf, hb, h=hf, epsilon=0.01)
  expect_equal(out$P, rbind(P, c(5/3, 1), c(2.5, -1), c(10/3, 1)))
  expect_equal(out$h, c(1, 4, 3, 4, 5, 6, 7, 8, 9, 7, 2, 8))
  expect_equal(out$gf, c(10, 3, 4, 12, NA, NA, 6, 11, 8, 2, 7, 5))
  expect_equal(out$gb, c(NA, 10, 2, 3, 12, 7, 11, 9, NA, 1, 8, 4))
  
})

test_that("paths with correspondances in each path work correctly", {

  P <- rbind(c(1,1),
             c(2,1),
             c(2.5,2),
             c(3,1),
             c(4,1),
             c(1,-1),
             c(2,-1),
             c(2.5,-2),
             c(3,-1),
             c(4,-1))

  gf <- c(2, 3, 4, 5, NA, NA, 6, 7, 8, 9)
  gb <- c(NA, 1, 2, 3, 4, 7, 8, 9, 10, NA)
  h <- 1:nrow(P)
  hf <- h
  hb <- h
  hf[2] <- 4
  hb[4] <- 2
  hf[9] <- 7
  hb[7] <- 9
  
  out <- stitch.paths(P, VF0=1, VF1=5, VB0=6, VB1=10, gf, gb, hf, hb, h=hf, epsilon=0.01)
  expect_equal(out$P, P)
  expect_equal(out$h, c(1, 7, 3, 7, 5, 6, 7, 8, 7, 10))
  expect_equal(out$gf, c(2, 3, 4, 5, NA, NA, 6, 7, 8, 9))
  expect_equal(out$gb, c(NA, 1, 2, 3, 4, 7, 8, 9, 10, NA))
})

test_that("paths with a double correspondance in one path works correctly", {
  P <- rbind(c(1,1),
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

  gf <- c(2, 3, 4, 5, NA, NA, 6, 7, 8, 9, 10, 11)
  gb <- c(NA, 1, 2, 3, 4, 7, 8, 9, 10, 11, 12, NA)
  h <- 1:nrow(P)
  hf <- h
  hb <- h
  hf[2] <- 4
  hb[4] <- 2
  hf[9] <- 7
  hb[7] <- 9
  hf[11] <- 9
  hb[9] <- 11
  
  out <- stitch.paths(P, VF0=1, VF1=5, VB0=6, VB1=12, gf, gb, hf, hb, h=hf, epsilon=0.01)
  expect_equal(out$P, P)
  expect_equal(out$h, c(1, 7, 3, 7, 5, 6, 7, 8, 7, 10, 7, 12))
  expect_equal(out$gf, gf)
  expect_equal(out$gb, gb)
})
