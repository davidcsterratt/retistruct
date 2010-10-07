retistruct.batch.revision <- function() {
  return(as.integer(gsub("Rev: ", "" ,gsub("\\$", "", "$Rev$"))))
}

list.dirs <- function(path='.') {
  files <- list.files(path, recursive=FALSE, full.name=TRUE)
  fi <- file.info(files)
  dirs <- row.names(fi[fi$isdir,])
  for (d in dirs) {
    dirs <- c(dirs,list.dirs(d))
  }
  return(dirs)
}

## retistruct.bactch() - Recurse through a directory tree, determining
## whether the directory contains valid raw data and markup, and
## performing the reconstruction if it does
##
## tldir     - the top level of the tree through which to recurse
## outputdir - directory in which to dump a log file and images
##
retistruct.batch <- function(tldir='.', outputdir=tldir, cpu.time.limit=1800,
                             device="pdf") {
  print(outputdir)
  datasets <- list.dirs(tldir)
  logdat <- data.frame()
  for (dataset in datasets) {
    Result <- ""
    ret <- -1
    logfile <- file.path(outputdir,
                         paste(retistruct.cli.basepath(dataset),
                               ".log", sep=""))
    print(dataset)
    is.data.dir <- try(check.datadir(dataset))

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
        print(ret)
        out <- read.csv(logfile)
        Result <- out[nrow(out),1]
        print(as.vector(Result))
      }
    }
    EOD <- NA
    nflip <- NA
    E <- NA
    if (ret==0) {
      retistruct.read.recdata()
      if (!is.null(r)) {
        EOD <- r$EOD
        nflip <- r$nflip
        E <- r$opt$value
      }
    }
    logdat <- rbind(logdat, data.frame(Dataset=dataset,
                                       Return=ret,
                                       Result=Result,
                                       E=E,
                                       nflip=nflip,
                                       EOD=EOD))
    write.csv(logdat, file.path(outputdir, "retistruct-batch.csv"))
  }
}

## retistruct.batch() - Recurse through a directory tree, determining
## whether the directory contains valid derived data and
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
