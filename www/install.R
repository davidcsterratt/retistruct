install.packages(c("retistruct", "retistructdemos"),
                 repos=c("http://cran.r-project.org",
                   "http://R-Forge.R-project.org"),
                 type=ifelse(Sys.info()[["sysname"]] == "Windows", "win.binary", "source"))


