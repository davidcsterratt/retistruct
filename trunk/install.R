## Set up libPaths to ensure that the libraries are installed in
## the R_LIBS_USER directory. This can be specified using the 
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

## Need to make sure we have a valid repository
options(repos="http://cran.r-project.org")

## It ought to be possible to get these dependendcies from the
## retistruct DESCRIPTION, but I don't know how!!
if (!require("foreign"))       install.packages("foreign")
if (!require("geometry"))      install.packages("geometry", repos="http://R-Forge.R-project.org")
if (!require("gWidgets"))      install.packages("gWidgets")
if (!require("gWidgetsRGtk2")) install.packages("gWidgetsRGtk2")
if (!require("cairoDevice"))   install.packages("cairoDevice")
if (!require("rgl"))           install.packages("rgl")
if (!require("R.matlab"))      install.packages("R.matlab")
if (!require("ttutils"))       install.packages("ttutils")
if (!require("multicore"))     install.packages("multicore")
if (!require("png"))           install.packages("png")
if (!require("akima"))         install.packages("akima")
if (!require("tripack"))       install.packages("tripack")
if (!require("Triangle"))      install.packages("Triangle", repos="http://R-Forge.R-project.org")
update.packages(Sys.getenv("R_LIBS_USER"))
## These lines needed while  Triangle and geometry only on R-Forge
install.packages("Triangle", repos="http://R-Forge.R-project.org")
## install.packages("geometry", repos="http://R-Forge.R-project.org")
## Now install retistruct
install.packages("retistruct_0.4.10.tar.gz")
install.packages("retistructgui_0.4.10.tar.gz")

library(retistructgui)
retistruct()
