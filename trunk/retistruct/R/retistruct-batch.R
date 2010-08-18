list.dirs <- function(path='.') {
  files <- list.files(path, recursive=FALSE, full.name=TRUE)
  fi <- file.info(files)
  dirs <- row.names(fi[fi$isdir,])
  for (d in dirs) {
    dirs <- c(dirs,list.dirs(d))
  }
  return(dirs)
}
  
retistruct.batch <- function(path='.') {
  datasets <- list.dirs(path)
  for (dataset in datasets) {
    ret <- system(paste("R --vanilla <<EOF
library(retistruct)
retistruct.cli(\"", dataset, "\", 600)
EOF", sep=""), intern=FALSE, wait=TRUE)
  }
}
