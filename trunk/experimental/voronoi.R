## This code was used to produce a Voronoi tesselation in projection.ReconstructedDataset()

plot.mosaic <- !is.null(args$mosaic) && args$mosaic

## Voroni
  if (plot.mosaic) {
    Dss <- getDss(r)
    if (length(Dss)) {
      for (i in 1:length(Dss)) {
        if (nrow(Dss[[i]]) >= 2) {
          ## Convert to angle-preserving coordinates
          pos <- sphere.spherical.to.polar.cart(Dss[[i]], preserve="angle")
          ## Create Voronoi mosaic in area-preserving coordinates
          vm <- voronoi.mosaic(pos[,1], pos[,2])

          ## Convert back to polar coords
          vmxy.polar <- polar.cart.to.sphere.spherical(cbind(x=vm$x, y=vm$y),
                                                       preserve="angle")
          ## Convert into whichever representation we're using
          ## and put back on to voroini object for plotting
          pos <- rho.to.degrees(sphere.spherical.to.polar.cart(vmxy.polar, pa), r$phi, pa)
          vm$x <- pos[,1]
          vm$y <- pos[,2]
          ## Do the same for the dummy coordinates
          vmxy.polar <- polar.cart.to.sphere.spherical(cbind(x=vm$dummy.x,
                                                             y=vm$dummy.y),
                                                       preserve="angle")
          ## Convert into whichever representation we're using
          ## and put back on to voroini object for plotting
          pos <- rho.to.degrees(sphere.spherical.to.polar.cart(vmxy.polar, pa), r$phi, pa)
          vm$dummy.x <- pos[,1]
          vm$dummy.y <- pos[,2]

          plot.voronoi.circular(vm, R=r$phi0*180/pi + 90, col=r$cols[[names(Dss)[i]]])
        }
      }
    }
  }
