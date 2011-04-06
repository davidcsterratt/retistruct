retistruct.cli.revision <- function() {
  return(as.integer(gsub("Rev: ", "" ,gsub("\\$", "", "$Rev$"))))
}

retistruct.cli <- function(dataset, cpu.time.limit=Inf, outputdir=NA,
                           device="pdf") {
  setTimeLimit(cpu=cpu.time.limit)
  dataset <<- dataset
  out <- try(retistruct.cli.process(outputdir=outputdir, device=device))
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

retistruct.cli.process <- function(outputdir=NA, device="pdf") {
  ## Processing
  retistruct.read.dataset()
  retistruct.read.markup()
  out <- try(retistruct.reconstruct())

  ## Output
  if (!inherits(out, "try-error")) {
    retistruct.save.recdata()
    if (!is.na(outputdir)) {
      retistruct.cli.figure(outputdir, device=device)
    }
  }
}

## retistruct.cli.basepath - generate a path based on the elided directory name
##
retistruct.cli.basepath <- function(dataset) {
  basepath <- gsub("\\./", "", dataset)
  basepath <- gsub("/", "_", basepath)
  return(basepath)
}

## retistruct.cli.figure - Print a figure to file
##
## It requires the global variable dataset to be set
retistruct.cli.figure <- function(outputdir, device="pdf", width=6, height=6,
                                  res=100) {
  retistruct.read.recdata()
  units <- NULL
  if (device!="pdf") {
    height <- height*res
    width  <- width*res
  }
  suffix <- paste(".", device, sep="")
  dev <- switch(device,
                pdf=pdf,
                png=png,
                jpeg=jpeg,
                tiff=tiff)
  if (is.null(dev)) {
    stop(paste("Device", device, "is not supported"))
  }
  if (!is.null(r)) {
    ## Determine the name of a figure
    basepath <- retistruct.cli.basepath(dataset)
    
    ## Flat plot
    dev(file=file.path(outputdir, paste(basepath, "-flat", suffix, sep="")),
           width=width, height=height)
    plot.outline.flat(r$P, r$gb, axt="s")
    with(r, plot.gridlines.flat(P, T, phi, lambda, Tt, phi0*180/pi))
    plot.datapoints.flat(r$Ds)
    plot.landmarks.flat(r$Ss, col="orange")
    plot.stitch.flat(r, add=TRUE)
    title(dataset)
    dev.off()

    ## Polar plot
    dev(file=file.path(outputdir, paste(basepath, "-polar", suffix, sep="")),
           width=width, height=height)
    plot.polar(r$phi0*180/pi)
    if (!is.null(r$Dss)) {
      plot.outline.polar(r)
      plot.datapoints.polar(r$Dss, cex=5)
    }
    if (!is.null(r$Sss)) {
      plot.landmarks.polar(r$Sss, col="orange")
    }
    title(dataset)
    if (!is.null(r$EOD)) {
      text.polar(paste("OD displacement:", format(r$EOD, digits=3, nsmall=2), "deg"))
    }
    dev.off()

    ## Strain plot
    dev(file=file.path(outputdir, paste(basepath, "-strain", suffix, sep="")),
           width=width, height=height)
    plot.outline.flat(r$P, r$gb, axt="s")
    plot.strain.flat(r)
    dev.off()

    ## l.vs.L plot
    dev(file=file.path(outputdir, paste(basepath, "-strain-lvsL", suffix, sep="")),
           width=width, height=height)
    plot.l.vs.L(r)
    dev.off()
  }
}
