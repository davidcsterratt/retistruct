## Set up libPaths to ensure that the libraries are installed in
## the R_LIBS_USER directory. This can be specified using the 
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

## Need to make sure we have a valid repository
options(repos="http://cran.r-project.org")

## It ought to be possible to get these dependendcies from the
## retistruct DESCRIPTION, but I don't know how!!
if (!require("foreign"))       install.packages("foreign")
if (!require("geometry"))      install.packages("geometry")
if (!require("gWidgets"))      install.packages("gWidgets")
if (!require("gWidgetsRGtk2")) install.packages("gWidgetsRGtk2")
if (!require("cairoDevice"))   install.packages("cairoDevice")
if (!require("rgl"))           install.packages("rgl")
if (!require("plotrix"))       install.packages("plotrix")
update.packages(Sys.getenv("R_LIBS_USER"))
## Now install retistruct
install.packages("retistruct_0.1.tar.gz")
