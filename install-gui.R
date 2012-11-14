source("install.R")

## It ought to be possible to get these dependendcies from the
## retistruct DESCRIPTION, but I don't know how!!
if (!require("gWidgets"))      install.packages("gWidgets")
if (!require("gWidgetsRGtk2")) install.packages("gWidgetsRGtk2")
if (!require("cairoDevice"))   install.packages("cairoDevice")
## These lines needed while  Triangle and geometry only on R-Forge
## update.packages(Sys.getenv("R_LIBS_USER"))
## Now install retistruct
install.packages("retistructgui_0.5.5.tar.gz")

library(retistructgui)
retistruct()
