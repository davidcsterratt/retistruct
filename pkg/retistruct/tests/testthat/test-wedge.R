grid.coords <- function(grid.int.f=0.5, grid.int.psi=45,
                                              phi0=135*pi/180) {
  ## Funny things happen when f=0
  fs <- seq(1E-6, 1-1E-6, by=grid.int.f - 1E-6)
  psis <- seq(-90, 90, by=grid.int.psi)*pi/180
  
  gfs   <- outer(fs, psis*0, "+")
  gpsis <- outer(fs*0, psis, "+")

  gc <- cbind(psi=as.vector(gpsis), f=as.vector(gfs))
  return(gc)
}

context("Checking wedge coordinates")

test_that("All points lie on sphere", {
  phi0 <- 135*pi/180
  gc <- grid.coords(phi0=phi0)
  P <- sphere.wedge.to.sphere.cart(gc[,"psi"], gc[,"f"], phi0=phi0)
  expect_that(rowSums(P^2), equals(rep(1, nrow(P))))
})

test_that("Points convert back", {
  phi0 <- 135*pi/180
  gc <- grid.coords(phi0=phi0)
  Pc <- sphere.wedge.to.sphere.cart(gc[,"psi"], gc[,"f"], phi0=phi0)
  Pt <- sphere.cart.to.sphere.wedge(Pc, phi0=phi0)
  ## print(Pt)
  ## print(gc)
  expect_that(Pt, equals(gc))
})
