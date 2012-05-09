## This code was used to plot an intensity plot of the KDE in projection.ReconstructedDataset()
plot.kde <- !is.null(args$kde) && args$kde

## FIXME: This should either go or be fixed
  if (plot.kde) {
    KDE <- getKDE(r)
    image(rho.to.degrees(KDE$red$g$xs, r$phi0, pa),
          rho.to.degrees(KDE$red$g$ys, r$phi0, pa),
          KDE$red$g$f, col=gray((0:100)/100))
    Dss <- getDss(r)
    pos <- sphere.spherical.to.polar.cart(Dss[["red"]], pa)
        suppressWarnings(points(rho.to.degrees(pos, r$phi0, pa),
                                col="red",
                                pch=20, cex=0.2, ...))

  }
  
