libdir <- file.path(getwd(), "library")
system(paste("mkdir", libdir))
.libPaths(new=libdir)
