source("datafile-utils.R")

## install.packages('gWidgetsRGtk2') first if not installed
if (!require("gWidgetsRGtk2")) install.packages("gWidgetsRGtk2")
if (!require("cairoDevice")) install.packages("cairoDevice")
library(gWidgetsRGtk2)
options(guiToolkit = "RGtk2")

add.tear <- function(h, ...) {
  print("Tear added")
}

fileChoose <- function(action="print", text = "Select a file...",
                       type="open", ...) {
  gfile(text=text, type=type, ..., action = action, handler =
        function(h,...) {
          do.call(h$action, list(h$file))
        })
}

datadir <- NULL

tbl <- glayout(container = gwindow("Tear editor"), spacing = 0)
tbl[1, 1:2, anchor = c(0, 0), expand = TRUE] = g.f = ggraphics(container = tbl,
    expand = TRUE, ps = 11)
d1 <- dev.cur()
tbl[1, 3:4, anchor = c(0, 0), expand = TRUE] = g.f2 = ggraphics(container = tbl,
    expand = TRUE, ps = 11)
d2 <- dev.cur()
tbl[2, 1,   anchor = c(1, 0)] = "Dataset"
tbl[2, 2:4, anchor = c(0, 0), expand = TRUE] <- g.dfn <- gbutton("/afs/inf.ed.ac.uk/user/s/sterratt/projects/rettect/data/Anatomy/marked-up-retinae-2010-03-24/gm119-5-adult-C57BL6", handler = function(h, ...) {
  curdir <- getwd()
  setwd(svalue(h$obj))
  setwd("..")
  gfile(type="selectdir", text="Select a directory...",
        handler = function(h, ...) {
          print(h$file)
          svalue(g.dfn) <- h$file
        })
  setwd(curdir)
  map <<- read.map(svalue(h$obj))
  sys <<- read.sys(svalue(h$obj))
}) 
tbl[3, 1, anchor = c(1, 0)] = "denominator df"
tbl[3, 2, anchor = c(0, 0), expand = TRUE] <- g.dfd <- gbutton("Add tear", handler=add.tear)
tbl[3, 3, anchor = c(0, 0), expand = TRUE] <- g.b.plot <- gbutton("Plot sys and map",
                              handler=function(h, ...) {
                                dev.set(d1)
                                plot.sys.map(sys, map)
                                })

tbl[3, 4, anchor = c(0, 0), expand = TRUE] <- g.b.plot <- gbutton("Plot map",
                              handler=function(h, ...) {
                                dev.set(d2)
                                plot.map(map)
                                })

