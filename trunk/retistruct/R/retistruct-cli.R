retistruct.cli <- function(dataset, cpu.time.limit=Inf, outputdir=NA) {
  setTimeLimit(cpu=cpu.time.limit)
  dataset <<- dataset
  out <- try(retistruct.cli.process(outputdir=outputdir))
  mess <- geterrmessage()
  if (inherits(out, "try-error")) {
    if (grepl("reached CPU time limit", mess)) {
      quit(status=1)
    } else {
      ## Unknown error
      quit(status=2)
    }
  }
  ## Success
  quit(status=0)
}

retistruct.cli.process <- function(outputdir=NA) {
  retistruct.read.dataset()
  retistruct.reconstruct()
  retistruct.save()
  if (!is.na(outputdir)) {
    retistruct.cli.figures(outputdir)
  }
}

retistruct.cli.figures <- function(outputdir) {
  basepath <- gsub("/", "_", dataset)
  ## Flat plot
  pdf(file=paste(outputdir, "/", basepath, "-flat.pdf", sep=""),
      width=6, height=6)
  plot.outline.flat(P, gb, axt="s")
  if (!is.null(r$t) && !is.null(r$r)) {
    with(r, plot.gridlines.flat(t$P, t$T, r$phi, r$lambda, m$Tt, p$phi0*180/pi))
  }
  plot.datapoints.flat(Ds)
  plot.landmarks.flat(Ss, col="orange")
  title(dataset)
  dev.off()

  ## Polar plot
  pdf(file=paste(outputdir, "/", basepath, "-polar.pdf", sep=""),
      width=6, height=6)
  plot.polar(phi0)
  if (!is.null(r$Dss)) {
    plot.outline.polar(r)
    plot.datapoints.polar(r$Dss, cex=5)
  }
  if (!is.null(r$Sss)) {
    plot.landmarks.polar(r$Sss, col="orange")
  }
  title(dataset)
  dev.off()
}
