##' List valid datasets underneath a directory. This reports all
##' directories that appear to be valid.
##'
##' @title List datasets underneath a directory
##' @param path Directory path to start searching from
##' @param verbose If \code{TRUE} report on progress
##' @return A vector of directories containing datasets
##' @author David Sterratt
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
      is.data.dir <- check.datadir(d)
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

##' Recurse through a directory tree, determining whether the
##' directory contains valid raw data and markup, and performing the
##' reconstruction if it does
##'
##' @title Retistruct batch operation 
##' @param tldir the top level of the tree through which to recurse
##' @param outputdir directory in which to dump a log file and images
##' @param cpu.time.limit amount of CPU after which to terminate the process
##' @param device string indicating what type of graphics output
##' required. Options are "pdf" and "png".
##' @author David Sterratt
retistruct.batch <- function(tldir='.', outputdir=tldir, cpu.time.limit=3600,
                             device="pdf") {
  print(outputdir)
  datasets <- list.datasets(tldir)
  logdat <- data.frame()
  for (dataset in datasets) {
    Result <- ""
    ret <- -1
    logfile <- file.path(outputdir,
                         paste(retistruct.cli.basepath(dataset),
                               ".log", sep=""))
    message(paste("Trying to open", dataset, "..."))
    is.data.dir <- try(check.datadir(dataset))
    EOD <- NA
    nflip <- NA
    nflip.init <- NA
    E <- NA
    El <- NA
    sqrt.El <- NA
    mean.strain <- NA
    mean.logstrain <- NA
    
    if (inherits(is.data.dir, "try-error")) {
      ret <- 1
      Result <- gsub("\n$", "", geterrmessage())
    } else {
      if (!is.data.dir) {
        message("... not a data directory.")
        next
      }
      if (is.data.dir) {
        message(paste("... reconstructing. Logging to", logfile))
        ret <- system(paste("R --vanilla >", logfile, "2>&1"),
                      input=paste("library(retistruct)
status = retistruct.cli(\"",
                        dataset, "\",",
                        cpu.time.limit, ",\"",
                        outputdir, "\", \"",
                        device, "\")
quit(status=status)",
                        sep=""),
                      intern=FALSE, wait=TRUE)
        if (ret==0) {
          r <- retistruct.read.recdata(list(dataset=dataset))
          if (!is.null(r)) {
            if (!is.null(r$EOD))         { EOD <- r$EOD      }
            if (!is.null(r$nflip))       { nflip <- r$nflip  }
            if (!is.null(r$nflip0))      { nflip0 <- r$nflip0  }
            if (!is.null(r$opt$value))   { E  <- r$opt$value }
            if (!is.null(r$E.l))         { El <- r$E.l
                                           sqrt.El <- sqrt(r$E.l) }
            if (!is.null(r$mean.strain)) { mean.strain <- r$mean.strain }
            if (!is.null(r$mean.logstrain)) { mean.logstrain <- r$mean.logstrain }
          }
        }
        message(paste("Return value", ret, ". Result:"))
        out <- read.csv(logfile)
        Result <- out[nrow(out),1]
        print(as.vector(Result))
      }
    }
    logdat <- rbind(logdat, data.frame(Dataset=dataset,
                                       Return=ret,
                                       Result=Result,
                                       E=E,
                                       El=El,
                                       nflip=nflip,
                                       nflip0=nflip0,
                                       EOD=EOD,
                                       sqrt.E=sqrt.E,
                                       mean.strain=mean.strain,
                                       mean.logstrain=mean.logstrain))
    write.csv(logdat, file.path(outputdir, "retistruct-batch.csv"))
  }
}

##' This function reconstructs a number of  datasets, using the R
##' \code{multicore} package to distribute the reconstruction of
##' multiple datasets across CPUs. If \code{datasets} is not specified
##' the function recurses through a directory tree starting at
##' \code{tldir}, determining whether the directory contains valid raw
##' data and markup, and performing the reconstruction if it does.
##'
##' @title Batch operation using multicore
##' @param datasets Vector of dataset directories to reconstruct.
##' @param tldir If datasets is not specified, the top level of the
##' directory tree through which to recurse in order to find datasets.
##' @param outputdir directory in which to dump a log file and images
##' @param device string indicating what type of graphics output
##' required. Options are "pdf" and "png".
##' @param cpu.time.limit amount of CPU after which to terminate the process
##' @author David Sterratt
retistruct.multicore <- function(tldir='.', datasets=NULL, outputdir=tldir,
                                 device="pdf",
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
    return(retistruct.cli(dataset, cpu.time.limit, outputdir, device))
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
##' @title Plot figures for a batch of reconstructions
##' @param tldir The top level directory of the tree through which to
##' recurse.
##' @return Data frame containing summary data
##' @author David Sterratt
retistruct.batch.summary <- function(tldir=".") {
  datasets <- list.datasets(tldir)
  logdat <- data.frame()

  ## Function to replace NULL with NA - needed for creating data frames
  n <- function(x) {
    return(ifelse(is.null(x), NA, x))
  }

  ## Go through datasets
  for (dataset in datasets) {
    r <- retistruct.read.recdata(list(dataset=dataset))
    if (!is.null(r)) {
      logdat <- rbind(logdat, data.frame(dataset=dataset,
                                         E=n(r$opt$value),
                                         El=n(r$E.l),
                                         nflip=n(r$nflip),
                                         EOD=n(r$EOD),
                                         sqrt.E=n(sqrt(r$E.l)),
                                         mean.strain=n(r$mean.strain),
                                         mean.logstrain=n(r$mean.logstrain),
                                         OD.phi=n(r$Dss$OD[1,"phi"]),
                                         OD.lambda=n(r$Dss$OD[1,"lambda"])))
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
retistruct.batch.export.matlab <- function(tldir=".") {
  datasets <- list.datasets(tldir)
  for (dataset in datasets) {
    print(dataset)
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
retistruct.batch.analyse.summary <- function(file) {
  dat <- read.csv(file)

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
       xlab=expression(sqrt(E[L])), main="")

  message("\nOUTLIERS")
  outliers <- subset(sdat, sqrt.E > (mean(sqrt.E) + 2*sd(sqrt.E)))
  outliers <- outliers[order(outliers[,"sqrt.E"], decreasing=TRUE),]
  ## print(outliers)
  return(invisible(list(N=nrow(sdat),
                        N.outtime=N.outtime,
                        N.fail=N.fail,
                        sqrt.E=sqrt.E, mean.strain=mean.strain,
                        mean.logstrain=mean.logstrain,
                        time=time, nflip=nflip,
                        with.flips=with.flips,
                        outliers=outliers)))
}

##' Extract statistics from a directory containing
##' reconstruction directories. 
##'
##' @title Extract statistics from a directory containing
##' reconstruction directories. 
##' @param path Directory containing recontstruction directories
##' @return Data frame containg various statistics 
##' @author David Sterratt
retistruct.batch.analyse.summaries <- function(path) {
  files <- list.files(path, recursive=FALSE, full.name=TRUE)
  fi <- file.info(files)
  dirs <- row.names(fi[fi$isdir,])
  out <- data.frame()
  for (d in dirs) {
    file <- file.path(d, "retistruct-batch.csv")
    if (file.exists(file)) {
      print(file)
      summ <- try(retistruct.batch.analyse.summary(file))
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
##' @param summ Summary object returned by \link{\code{retistruct.batch.summary}} 
##' @return A pseudo retina, in which the optic disks are treated as
##' datapoints 
##' @author David Sterratt
retistruct.batch.plot.ods <- function(summ) {
  ## Make a dummy retina
  o <- list()
  class(o) <- "reconstructedOutline"
  o$phi0 <- 20*pi/180
  r <- ReconstructedDataset(o)
  r$Dss$OD <- na.omit(summ[,c("OD.phi","OD.lambda")])
  colnames(r$Dss$OD) <- c("phi", "lambda")
  r$cols["OD"] <- "blue"
  plot.polar(r)
  return(r)
}
