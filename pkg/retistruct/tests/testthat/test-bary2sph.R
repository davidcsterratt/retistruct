context("bary2sph")
test_that("bary2sph works correctly", {
  ## Construct a diamond round the orgin
  P <- rbind(c(0, 0, 1),
             c(1, 0, 0),
             c(0, 1, 0),
             c(-1, 0, 0),
             c(0, -1, 0),
             c(0, 0, -1))

  T <- rbind(c(1, 2, 3),
             c(1, 3, 4),
             c(1, 4, 5),
             c(1, 5, 2),
             c(6, 2, 3),
             c(6, 3, 4),
             c(6, 4, 5),
             c(6, 5, 2))
  Ib <- list(idx=c(1, 1, 1, 1, 2, 3, 5),
             p=rbind(c(1, 0, 0),
                     c(0, 1, 0),
                     c(0, 0, 1),
                     c(0, 1/2, 1/2),
                     c(0, 0 , 1),
                     c(0, 0, 1),
                     c(1, 0, 0)))
  expect_equal(bary2sph(Ib, T=T, P=P)/pi*2,
               rbind(c(phi=1, lambda=0),
                     c(0, 0),
                     c(0, 1),
                     c(0, 0.5),
                     c(0, 2.0),
                     c(0,-1.0),
                     c(-1, 0.0)))
  expect_equal(bary2sph(Ib=list(idx=c(1, NA), p=rbind(c(1, 0, 0), c(NA, NA, NA))), T=T, P=P)/pi*2,
               rbind(c(phi=1, lambda=0),
                     c(NA, NA)))
})
