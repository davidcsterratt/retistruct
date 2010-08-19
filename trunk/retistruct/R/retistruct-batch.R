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
retistruct.batch <- function(tldir='.', outputdir=tldir) {
  print(outputdir)
  datasets <- list.dirs(tldir)
  logdat <- data.frame()
  for (dataset in datasets) {
    logfile <- paste(dataset, "/retistruct.log", sep="")
    print(dataset)
    ret <<- system(paste("R --vanilla >", logfile, "2>&1"),
                  input=paste("library(retistruct)
retistruct.cli(\"", dataset, "\", 600, \"", outputdir, "\")", sep=""),
                  intern=FALSE, wait=TRUE)
    stdout <- read.csv(logfile)
    logdat <- rbind(logdat, data.frame(Dataset=dataset,
                                       Return=ret,
                                       Result=stdout[nrow(stdout),1]))
    write.csv(logdat, paste(outputdir, "/retistruct-batch.csv", sep=""))
  }
}

## retistruct.bactch() - Recurse through a directory tree, determining
## whether the directory contains valid derived data and
## plotting graphs if it does.
##
## tldir     - the top level of the tree through which to recurse
## outputdir - directory in which to dump a log file and images
##
retistruct.batch.figures <- function(tldir=".", outputdir=tldir) {
  datasets <- list.dirs(tldir)
  for (d in datasets) {
    print(d)
    dataset <<- d
    retistruct.read.recdata()
    if (!is.null(r)) {
      retistruct.cli.figures(outputdir)
    }
  }
}
