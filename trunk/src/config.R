libdir <- file.path(getwd(), "library")
if (is.na(file.info(libdir)$isdir)) {
   system(paste("mkdir", libdir))
}
.libPaths(new=libdir)
