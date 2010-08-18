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
    print(dataset)
    ret <<- system("R --vanilla 2>&1",
                  input=paste("library(retistruct)
retistruct.cli(\"", dataset, "\", 600)", sep=""),
                  intern=TRUE, wait=TRUE)
    print(ret)
    write(ret, paste(dataset, "/retistruct.log", sep=""))
    logdat <- rbind(logdat, data.frame(Dataset=dataset, Result=ret[length(ret)]))
    write.csv(logdat, paste(outputdir, "/retistruct-batch.csv", sep=""))
  }
}
