list.dirs <- function(path='.') {
  files <- list.files(path, recursive=FALSE, full.name=TRUE)
  fi <- file.info(files)
  dirs <- row.names(fi[fi$isdir,])
  for (d in dirs) {
    dirs <- c(dirs,list.dirs(d))
  }
  return(dirs)
}
  
