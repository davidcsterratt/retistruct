##' List valid datasets underneath a directory. This reports all
##' directories that appear to be valid.
##'
##' @title List datasets underneath a directory
##' @param path Directory path to start searching from
##' @param verbose If \code{TRUE} report on progress
##' @return A vector of directories containing datasets
##' @author David Sterratt
##' @export
list.datasets <- function(path='.', verbose=FALSE) {
  ## We have to define a function to do directory listings as
  ## (a) list.files doesn't work recursively on all platforms and (b)
  ## in any case, it doesn't list directories when running recursively
  list.dirs <- function(path='.') {
    files <- list.files(path, recursive=FALSE, full.name=TRUE)
    fi <- file.info(files)
    dirs <- row.names(fi[fi$isdir,])
    for (d in dirs) {
      dirs <- c(dirs,list.dirs(d))
    }
    return(dirs)
  }

  ## Now get the directories
  dirs <- list.dirs(path)

  ## Go through directories determining if datasets are valid
  datasets <- c()
  for (d in dirs) {
    ## Determine if directory is a valid data directory. If it's a
    ## faulty one, we will let it pass to be picked up later on
    
    id.data.dir <- TRUE
    ## Case of faulty directory
    tryCatch({
      is.data.dir <- is.character(check.datadir(d))
    }, error=function(e) {})

    if (!is.data.dir) {
      if (verbose) message(paste(d, "is not a data directory."))
      next
    }
    datasets <- c(datasets, d)
    if (verbose) message(paste(d, "is a data directory."))
  }
  return(datasets)
}

##' This function reconstructs a number of  datasets, using the R
##' \code{multicore} package to distribute the reconstruction of
##' multiple datasets across CPUs. If \code{datasets} is not specified
##' the function recurses through a directory tree starting at
##' \code{tldir}, determining whether the directory contains valid raw
##' data and markup, and performing the reconstruction if it does.
##'
##' @title Batch operation using the multicore package
##' @param tldir If datasets is not specified, the top level of the
##' directory tree through which to recurse in order to find datasets.
##' @param outputdir directory in which to dump a log file and images
##' @param datasets Vector of dataset directories to reconstruct
##' @param device string indicating what type of graphics output
##' required. Options are "pdf" and "png".
##' @param cpu.time.limit amount of CPU after which to terminate the
##' process
##' @param mc.cores The number of cores to use. Defaults to the total
##' number available.
##' @author David Sterratt
##' @export
retistruct.batch <- function(tldir='.', outputdir=tldir, datasets=NULL, 
                             device="pdf", titrate=FALSE,
                             cpu.time.limit=3600,
                             mc.cores=getOption("cores")) {
  ## Get datasets
  if (is.null(datasets)) {
    datasets <- list.datasets(tldir)
  }
  message(paste("About to reconstruct", length(datasets), "datasets."))

  ## Function to pass to mclapply
  multicore.call <- function(dataset) {
    logfile <- file.path(outputdir,
                         paste(retistruct.cli.basepath(dataset),
                               ".log", sep=""))
    message(paste("Reconstructing. Logging to", logfile))
    flog <- file(logfile, open="wt")
    sink(flog)
    sink(flog, type="message")
    return(retistruct.cli(dataset, cpu.time.limit, outputdir, device,
                          titrate=titrate))
  }

  ## Run the reconstructions
  ret <- mclapply(datasets, multicore.call, mc.preschedule=FALSE,
                  mc.cores=mc.cores)

  ## Extract data from the return structures
  ## Function to replace NULL with NA - needed for creating data frames
  n <- function(x) {
    return(ifelse(is.null(x), NA, x))
  }

  dat <- data.frame(cbind(dataset=datasets,
                          status=sapply(ret, function(x) {return(n(x$status))}),
                          mess=as.character(sapply(ret, function(x) {return(x$mess)})),
                          time=sapply(ret, function(x) {return(n(x$time))})))

  summ <- retistruct.batch.summary(tldir)
  dat <- merge(summ, dat, by="dataset", all=TRUE)
  write.csv(dat, file.path(outputdir, "retistruct-batch.csv"))
  return(dat)
}

##' Recurse through a directory tree, determining whether the
##' directory contains valid derived data and extracting summary data
##' if it does.
##'
##' @title Extract summary data for a batch of reconstructions
##' @param tldir The top level directory of the tree through which to
##' recurse.
##' @param cache If \code{TRUE} use the cached statistics rather than
##' generate on the fly (which is slower).
##' @return Data frame containing summary data
##' @author David Sterratt
##' @export
retistruct.batch.summary <- function(tldir=".", cache=TRUE) {
  datasets <- list.datasets(tldir)
  logdat <- data.frame()

  ## Function to replace NULL with NA - needed for creating data frames
  n <- function(x) {
    return(ifelse(is.null(x), NA, x))
  }

  ## Go through datasets
  for (dataset in datasets) {
    message(paste("Reading", dataset))
    r <- retistruct.read.recdata(list(dataset=dataset))
    if (!is.null(r)) {
      dat <- data.frame(dataset=dataset,
                        E=n(r$opt$value),
                        El=n(r$E.l),
                        nflip=n(r$nflip),
                        EOD=n(r$EOD),
                        sqrt.E=n(sqrt(r$E.l)),
                        mean.strain=n(r$mean.strain),
                        mean.logstrain=n(r$mean.logstrain),
                        OD.phi=n(r$Dss$OD[1,"phi"]),
                        OD.lambda=n(r$Dss$OD[1,"lambda"]),
                        mean.dtheta=n(r$titration$Dtheta.mean),
                        phi0d=n(r$phi0*180/pi),
                        phi0d.opt=n(r$titration$phi0d.opt))
      message(paste("Getting KDE"))
      KDE <- getKDE(r, cache)
      if (length(KDE) > 0) {
        KDEdat <- lapply(KDE, function(x) {x$h})
        names(KDEdat) <- paste("h.", names(KDEdat), sep="")
        dat <- cbind(dat, KDEdat)
      }
      logdat <- merge(logdat, dat, all=TRUE)
    }
  }
  return(logdat)
}

##' Recurse through a directory tree, determining whether the
##' directory contains valid derived data and plotting graphs if it
##' does.
##'
##' @title Plot figures for a batch of reconstructions
##' @param tldir The top level directory of the tree through which to
##' recurse.
##' @param outputdir Directory in which to dump a log file and images
##' @param ... Parameters passed to plotting functions
##' @author David Sterratt
##' @export
retistruct.batch.figures <- function(tldir=".", outputdir=tldir, ...) {
  datasets <- list.datasets(tldir)
  for (dataset in datasets) {
    message(paste("Attempting to produce figures from", dataset))
    try(retistruct.cli.figure(dataset, outputdir, ...))
  }
}

##' Recurse through a directory tree, determining whether the
##' directory contains valid derived data and converting r.rData files
##' to files in matlab format named r.mat
##'
##' @title Export data from reconstruction data files to matlab
##' @param tldir The top level of the directory tree through which to
##' recurse
##' @author David Sterratt
##' @export
retistruct.batch.export.matlab <- function(tldir=".") {
  datasets <- list.datasets(tldir)
  for (dataset in datasets) {
    r <- retistruct.read.recdata(list(dataset=dataset))
    retistruct.export.matlab(r)
  }
}

##' Extract statistics from the retistruct-batch.csv summary file
##'
##' @title Extract statistics from the retistruct-batch.csv summary file
##' @param file The path to the retistruct-batch.csv
##' @return list of various statistics
##' @author David Sterratt
##' @export
retistruct.batch.analyse.summary <- function(path) {
  dat <- read.csv(file.path(path, "retistruct-batch.csv"))
  par(mfcol=c(1, 3))
  par(mar=c(2.4, 2.3, 0.7, 0.2))
  par(mgp=c(1.3, 0.3, 0), tcl=-0.3)
  
  ## Detailed output codes
  etab <- sort(table(dat[,"mess"]), decreasing=TRUE)
  message("OUPUT CODES")
  message(rbind(format(etab, width=4), " ", paste(gsub("\n", "\n   ", gsub("\n$", "", names(etab))), "\n")))

  sqrt.E.fail <- 0.2
  
  ## Get number of failures due to lack of time
  N.outtime <- sum(na.omit(dat$status) == 1)

  ## Get successful reconstructions
  sdat <- subset(dat, !is.na(sqrt.E) & status==0)

  ## Get number of failures due to sqrt.E being too large
  N.fail <- nrow(subset(sdat, sqrt.E >= sqrt.E.fail))
  sdat <- subset(sdat, sqrt.E < sqrt.E.fail)
  
  message("\nSTATISTICS")
  
  message("sqrt.E")
  sqrt.E <- summary(sdat[,"sqrt.E"])
  print(sqrt.E)
  message("mean.strain")
  mean.strain <- summary(sdat[,"mean.strain"])
  print(mean.strain)
  message("mean.logstrain")
  mean.logstrain <- summary(sdat[,"mean.logstrain"])
  print(mean.logstrain)
  message("time")
  time <- summary(sdat[,"time"])
  print(time)
  message("nflip")
  nflip <- summary(sdat[,"nflip"])
  print(nflip)
  message("%with.flips")
  with.flips <- mean(sdat[,"nflip"] > 0) * 100
  
  hist(sdat[,"sqrt.E"], breaks=seq(0, max(sdat[,"sqrt.E"]), len=100),
       xlab=expression(sqrt(italic(E)[L])), main="")
  mtext("A", adj=-0.15, font=2, line=-0.7)
  
  message("\nOUTLIERS")
  outliers <- subset(sdat, sqrt.E > (mean(sqrt.E) + 2*sd(sqrt.E)))
  outliers <- outliers[order(outliers[,"sqrt.E"], decreasing=TRUE),]
  ## print(outliers)

  ## Plot of age versus goodness
  ## Find datasets containing ("Pxx")
  fage <- grepl("^.*(P\\d+).*$",  sdat$dataset)
  sdat$age <- sub("^.*P(\\d+).*$", "\\1", sdat$dataset)
  sdat$age <- sub(".*adult.*", "adult", sdat$age)
  sdat$age[grepl(".{6}", sdat$age)] <- NA
  sdat$age <- ordered(sdat$age, c(unique(sort(as.numeric(sdat$age))), "adult"))
  ## print(factor(sdat$age))
  ## print(sort(as.numeric(factor(sdat$age))))
  ## levels(sdat$age) <- sub("(\\d+)", "P\\1", levels(sdat$age))
  levels(sdat$age) <- sub("adult", "A", levels(sdat$age))
  ##  print(factor(sdat$age))
  with(sdat, boxplot(sqrt.E ~ age,
                     xaxt="n",
                     xlab="Postnatal day",
                     ylab=expression(sqrt(italic(E)[L]))))
  mtext("B", adj=-0.15, font=2, line=-0.7)
  axis(1, labels=NA, at=seq(1, len=length(levels(sdat$age))))
  mtext(levels(sdat$age), 1, at=seq(1, len=length(levels(sdat$age))), line=0.3, cex=0.8)

  par(mar=c(1,1,0.7,1))
  retistruct.batch.plot.ods(sdat)
  mtext("C", adj=-0.05, font=2, line=-0.7)
  
  ## with(sdat, table(sqrt.E ~ age))
  dev.copy2pdf(file=file.path(path, "retistruct-goodness.pdf"), width=6.83, height=6.83/3)
  return(invisible(list(N=nrow(sdat),
                        N.outtime=N.outtime,
                        N.fail=N.fail,
                        sqrt.E=sqrt.E, mean.strain=mean.strain,
                        mean.logstrain=mean.logstrain,
                        time=time, nflip=nflip,
                        with.flips=with.flips,
                        outliers=outliers,
                        sdat=sdat)))
}

##' Extract statistics from a directory containing
##' reconstruction directories. 
##'
##' @title Extract statistics from a directory containing
##' reconstruction directories. 
##' @param path Directory containing recontstruction directories
##' @return Data frame containg various statistics 
##' @author David Sterratt
##' @export
retistruct.batch.analyse.summaries <- function(path) {
  files <- list.files(path, recursive=FALSE, full.name=TRUE)
  fi <- file.info(files)
  dirs <- row.names(fi[fi$isdir,])
  out <- data.frame()
  for (d in dirs) {
    file <- file.path(d, "retistruct-batch.csv")
    if (file.exists(file)) {
      print(file)
      summ <- try(retistruct.batch.analyse.summary(d))
      try(print(summ$sqrt.E["Median"]))
      try(out <- rbind(out, data.frame(file=file,
                                       N=summ$N,
                                       N.outtime=summ$N.outtime,
                                       N.fail=summ$N.fail,
                                       sqrt.E.Median=summ$sqrt.E["Median"],
                                       sqrt.E.Mean=summ$sqrt.E["Mean"],
                                       nflip.Median=summ$nflip["Median"],
                                       nflip.Mean=summ$nflip["Mean"],
                                       nflip.Max=summ$nflip["Max."],
                                       pc.flipped=summ$with.flips,
                                       time.Mean=summ$time["Mean"])))
    }
  }
  return(out)
}

##' Polar plot of ODs of a group of retinae.
##'
##' @title Superposed plot of ODs on polar axes
##' @param summ Summary object returned by
##' \code{\link{retistruct.batch.summary}}
##' @return A pseudo retina, in which the optic disks are treated as
##' datapoints 
##' @author David Sterratt
##' @export
retistruct.batch.plot.ods <- function(summ) {
  ## Make a dummy retina
  o <- list()
  class(o) <- "reconstructedOutline"
  o$phi0 <- -60*pi/180
  r <- ReconstructedDataset(o)
  summ <- subset(summ, age=="A")
  r$Dss$OD <- na.omit(summ[,c("OD.phi","OD.lambda")])
  colnames(r$Dss$OD) <- c("phi", "lambda")
  r$cols["OD"] <- "blue"
  
  km <- karcher.mean.sphere(r$Dss$OD, na.rm=TRUE, var=TRUE)
  message(nrow(summ), " points")
  message("Mean: Lat ", format(km$mean["phi"]*180/pi, digits=3),
          " Long ", format(km$mean["lambda"]*180/pi, digits=3),
          " ; SD: ", format(sqrt(km$var)*180/pi, digits=3))
  message("Mean location is ", 180/pi*central.angle(km$mean["phi"],
                                      km$mean["lambda"],
                                      -pi/2,
                                      0),  " away from geometric centre")
  
  summ$OD.res <- 180/pi*central.angle(km$mean["phi"],
                                      km$mean["lambda"],
                                      r$Dss$OD[,"phi"],
                                      r$Dss$OD[,"lambda"])

  summlm <- lm(OD.res ~ sqrt.E, summ)
  print(summary(summlm))
  plot.polar(r, datapoint.contours=FALSE)
  ##with(summ, plot(sqrt.E, OD.res))
  ##abline(summlm)
  
  return(r)
}
