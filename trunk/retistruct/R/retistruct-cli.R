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
    retistruct.cli.figure(outputdir)
  }
}

## retistruct.cli.figure - Print a figure to file
##
## It requires the global variable dataset to be set
retistruct.cli.figure <- function(outputdir) {
  retistruct.read.recdata()
  if (!is.null(r)) {
    ## Determine the name of a figure
    basepath <- gsub("\\./", "", dataset)
    basepath <- gsub("/", "_", basepath)
    
    ## Flat plot
    pdf(file=file.path(outputdir, paste(basepath, "-flat.pdf", sep="")),
        width=6, height=6)
    plot.outline.flat(r$P, r$gb, axt="s")
    with(r, plot.gridlines.flat(P, T, phi, lambda, Tt, phi0*180/pi))
    plot.datapoints.flat(r$Ds)
    plot.landmarks.flat(r$Ss, col="orange")
    title(dataset)
    dev.off()

    ## Polar plot
    pdf(file=file.path(outputdir, paste(basepath, "-polar.pdf", sep="")),
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
}
