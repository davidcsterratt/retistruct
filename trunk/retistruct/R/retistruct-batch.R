list.dirs <- function(path='.') {
  files <- list.files(path, recursive=FALSE, full.name=TRUE)
  fi <- file.info(files)
  dirs <- row.names(fi[fi$isdir,])
  for (d in dirs) {
    dirs <- c(dirs,list.dirs(d))
  }
  return(dirs)
}

##' Recurse through a directory tree, determining whether the
##' directory contains valid raw data and markup, and performing the
##' reconstruction if it does
##'
##' @title Retistruct batch operation 
##' @param tldir the top level of the tree through which to recurse
##' @param outputdir directory in which to dump a log file and images
##' @param cpu.time.limit amount of CPU after which to terminate the process
##' @param device string indicating what type of graphics output required. Options are "pdf" and "png" 
##' @author David Sterratt
retistruct.batch <- function(tldir='.', outputdir=tldir, cpu.time.limit=1800,
                             device="pdf") {
  print(outputdir)
  datasets <- list.dirs(tldir)
  logdat <- data.frame()
  for (d in datasets) {
    dataset <<- d
    Result <- ""
    ret <- -1
    logfile <- file.path(outputdir,
                         paste(retistruct.cli.basepath(dataset),
                               ".log", sep=""))
    print(dataset)
    is.data.dir <- try(check.datadir(dataset))
    EOD <- NA
    nflip <- NA
    E <- NA
    El <- NA
    
    if (inherits(is.data.dir, "try-error")) {
      ret <- 1
      Result <- gsub("\n$", "", geterrmessage())
    } else {
      if (!is.data.dir) {
        next
      }
      if (is.data.dir) {
        ret <- system(paste("R --vanilla >", logfile, "2>&1"),
                      input=paste("library(retistruct)
retistruct.cli(\"",
                        dataset, "\",",
                        cpu.time.limit, ",\"",
                        outputdir, "\", \"",
                        device, "\")",
                        sep=""),
                      intern=FALSE, wait=TRUE)
        if (ret==0) {
          retistruct.read.recdata()
          if (!is.null(r)) {
            if (!is.null(r$EOD))       { EOD <- r$EOD      }
            if (!is.null(r$nflip))     { nflip <- r$nflip  }
            if (!is.null(r$opt$value)) { E  <- r$opt$value }
            if (!is.null(r$E.l))       { El <- r$E.l       }
          }
        }
        print(ret)
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
                                       EOD=EOD))
    write.csv(logdat, file.path(outputdir, "retistruct-batch.csv"))
  }
}

## retistruct.batch.figures() - Recurse through a directory tree,
## determining whether the directory contains valid derived data and
## plotting graphs if it does.
##
## tldir     - the top level of the tree through which to recurse
## outputdir - directory in which to dump a log file and images
##
retistruct.batch.figures <- function(tldir=".", outputdir=tldir, ...) {
  retistruct.initialise.userdata()
  datasets <- list.dirs(tldir)
  for (d in datasets) {
    print(d)
    dataset <<- d
    try(retistruct.cli.figure(outputdir, ...))
  }
}

## retistruct.batch.export.matlab() - Recurse through a directory tree,
## determining whether the directory contains valid derived data and
## converting r.rData files to files in matlab format named r.mat
##
## tldir     - the top level of the tree through which to recurse
##
retistruct.batch.export.matlab <- function(tldir=".") {
  retistruct.initialise.userdata()
  datasets <- list.dirs(tldir)
  for (d in datasets) {
    print(d)
    dataset <<- d
    r <<- NULL
    retistruct.read.recdata()
    retistruct.export.matlab()
  }
}

## retistruct.batch.rdata2hdf() - Recurse through a directory tree,
## determining whether the directory contains valid derived data and
## converting r.rData files to r.h5
##
## tldir     - the top level of the tree through which to recurse
##
retistruct.batch.rdata2hdf <- function(tldir=".", ...) {
  retistruct.initialise.userdata()
  datasets <- list.dirs(tldir)
  for (d in datasets) {
    print(d)
    dataset <<- d
    r <<- NULL
    retistruct.read.recdata()
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
  retistruct.initialise.userdata()
  datasets <- list.dirs(tldir)
  failures <- c()
  successes <- c()
  for (d in datasets) {
    print(d)
    dataset <<- d
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
