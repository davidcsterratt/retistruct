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
                                       EOD=EOD,
                                       sqrt.E=sqrt.E,
                                       mean.strain=mean.strain,
                                       mean.logstrain=mean.logstrain))
    write.csv(logdat, file.path(outputdir, "retistruct-batch.csv"))
  }
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

## retistruct.batch.rdata2hdf() - Recurse through a directory tree,
## determining whether the directory contains valid derived data and
## converting r.rData files to r.h5
##
## tldir     - the top level of the tree through which to recurse
##
retistruct.batch.rdata2hdf <- function(tldir=".", ...) {
  datasets <- list.datasets(tldir)
  for (dataset in datasets) {
    print(dataset)
    r <- retistruct.read.recdata(list(dataset=dataset))
    if (!is.null(r)) {
      f <- file.path(dataset, "r.h5")
      print(paste("Saving", f))
      hdf5save(f, "r")
      print(paste("Trying to load", f))
      try(hdf5load(f))
    }
  }
}

retistruct.batch.testhdf <- function(tldir=".", ...) {
  datasets <- list.datasets(tldir)
  failures <- c()
  successes <- c()
  for (dataset in datasets) {
    print(dataset)
    f <- file.path(dataset, "r.h5")
    if (file.exists(f)) {
      print(paste(f, "exists"))
      e <- try(hdf5load(f, verbosity=0))
      if (inherits(e, "try-error")) {
        failures <- c(failures, f)
      } else {
        successes <- c(successes, f)
      }
    }
  }
  print(failures)
  print(successes)
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
                                 cpu.time.limit=3600) {
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

  ret <- mclapply(datasets, multicore.call, mc.preschedule=FALSE)
  return(data.frame(cbind(dataset=datasets, ret=ret)))
}

