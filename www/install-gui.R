install.packages(c("gWidgets", "gWidgetsRGtk2", "cairoDevice"),
                 repos=c("http://cran.r-project.org"))
install.packages(c("retistructgui"),
                 repos=c("http://www.neuralmapformation.org/R/",
                     "http://cran.r-project.org", 
                     "http://R-Forge.R-project.org"),
                 type=ifelse(Sys.info()[["sysname"]] == "Windows", "win.binary", "source"))
library(retistructgui)
