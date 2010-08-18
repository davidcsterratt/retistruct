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
retistruct.batch <- function(tldir='.', outputdir='.') {
  datasets <- list.dirs(tldir)
  logdat <- data.frame()
  for (dataset in datasets) {
    ret <- system(paste("R --vanilla <<EOF
library(retistruct)
retistruct.cli(\"", dataset, "\", 600)
EOF", sep=""), intern=FALSE, wait=TRUE)
    logdat <- rbind(logdat, data.frame(Dataset=dataset, Return=ret))
    write.csv(logdat, paste(outputdir, "/retistruct.log", sep=""))
  }
}
