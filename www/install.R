install.packages(c("geometry", "R.matlab", "ttutils", "rgl", "png", "sp", "RImageJROI", "RTriangle"), repos="http://cran.r-project.org")
install.packages(c("retistruct", "retistructdemos"),
                 repos=c("http://www.neuralmapformation.org/R/",
                     "http://cran.r-project.org",
                     "http://R-Forge.R-project.org"),
                 type=ifelse(Sys.info()[["sysname"]] == "Windows", "win.binary", "source"))


