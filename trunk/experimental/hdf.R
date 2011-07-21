
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

