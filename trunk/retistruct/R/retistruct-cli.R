retistruct.cli <- function(dataset, cpu.time.limit=Inf, outputdir=NA,
                           device="pdf") {
  ## Return code
  status <- 0
  setTimeLimit(cpu=cpu.time.limit)
  syst <- system.time(out <- tryCatch(retistruct.cli.process(dataset,
                                                             outputdir=outputdir, device=device),
                                      error=function(e) {return(e)}))
  mess <- "Success"
  if (inherits(out, "error")) {
    mess <- as.character(out)
    if (grepl("reached CPU time limit", mess)) {
      status <- 1
    } else {
      ## Unknown error
      status <- 2
    }
  }
  ## Success
  return(list(status=status, time=syst["user.self"], mess=mess))
}

retistruct.cli.process <- function(dataset, outputdir=NA, device="pdf") {
  ## Processing
  r <- retistruct.read.dataset(dataset)
  r <- retistruct.read.markup(r)
  r <- retistruct.reconstruct(r)

  ## Output
  retistruct.save.recdata(r)
  
  if (!is.na(outputdir)) {
    retistruct.cli.figure(dataset, outputdir, device=device)
  }

  ## Export to matlab
  retistruct.export.matlab(r)
}

## retistruct.cli.basepath - generate a path based on the elided directory name
##
retistruct.cli.basepath <- function(dataset) {
  basepath <- gsub("\\./", "", dataset)
  basepath <- gsub("/", "_", basepath)
  basepath <- gsub(" ", "_", basepath)
  return(basepath)
}

## retistruct.cli.figure - Print a figure to file
##
## It requires the global variable dataset to be set
retistruct.cli.figure <- function(dataset,
                                  outputdir, device="pdf", width=6, height=6,
                                  res=100) {
  r <- retistruct.read.recdata(list(dataset=dataset))
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
    par(mar=c(1, 1, 1, 1))
    plot.flat(r, axt="n",
              datapoints=TRUE,
              landmarks=TRUE,
              markup=FALSE,
              stitch=TRUE,
              grid=TRUE,
              mesh=FALSE,
              strain=FALSE)
    title(dataset)
    dev.off()

    ## Polar plot
    dev(file=file.path(outputdir, paste(basepath, "-polar", suffix, sep="")),
           width=width, height=height)
    par(mar=c(2, 2, 2, 2))
    plot.polar(r)
    title(dataset)
    if (!is.null(r$EOD)) {
      text.polar(paste("OD displacement:", format(r$EOD, digits=3, nsmall=2), "deg"))
    }
    dev.off()

    ## Strain plot
    dev(file=file.path(outputdir, paste(basepath, "-strain", suffix, sep="")),
           width=width, height=height)
    par(mar=c(1, 1, 1, 1))
    plot.flat(r, axt="n",
              datapoints=FALSE,
              landmarks=FALSE,
              markup=FALSE,
              stitch=FALSE,
              grid=FALSE,
              mesh=FALSE,
              strain=TRUE)
    title(dataset)
    dev.off()

    ## l.vs.L plot
    dev(file=file.path(outputdir, paste(basepath, "-strain-lvsL", suffix, sep="")),
           width=width, height=height)
    par(mar=c(3.0, 3.0, 1.5, 0.5))
    par(mgp=c(1.5, 0.5, 0))
    par(tcl=-0.3)
    
    plot.l.vs.L(r)
    title(dataset)
    dev.off()
  }
}
