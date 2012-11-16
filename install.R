## Set up libPaths to ensure that the libraries are installed in
## the R_LIBS_USER directory. This can be specified using the 
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

## Need to make sure we have a valid repository
options(repos=c("http://cran.r-project.org"))

## It ought to be possible to get these dependendcies from the
## retistruct DESCRIPTION, but I don't know how!!
if (!require("foreign"))       install.packages("foreign")
if (!require("geometry"))      install.packages("geometry")
if (!require("rgl"))           install.packages("rgl")
if (!require("R.matlab"))      install.packages("R.matlab")
if (!require("ttutils"))       install.packages("ttutils")
if (!require("png"))           install.packages("png")
if (!require("sp"))            install.packages("sp")
if (!require("RTriangle"))     install.packages("RTriangle")

## Updating packates may be necessary
## update.packages(Sys.getenv("R_LIBS_USER"))
## Now install retistruct
install.packages("retistruct_0.5.7.tar.gz")
install.packages("retistructdemos_0.5.7.tar.gz")

