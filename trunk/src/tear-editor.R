source("datafile-utils.R")

## install.packages('gWidgetsRGtk2') first if not installed
if (!require("gWidgetsRGtk2")) install.packages("gWidgetsRGtk2")
if (!require("cairoDevice")) install.packages("cairoDevice")
library(gWidgetsRGtk2)
options(guiToolkit = "RGtk2")

A <- c()
VB <- c()
VF <- c()

enable.group <- function(widgets, state=TRUE) {
  for (w in widgets) {
      enabled(w) <- state
    }
}

enable.widgets <- function(state) {
  enable.group(c(g.add, g.move, g.remove), state)
}

add.tear <- function(h, ...) {
  enable.widgets(FALSE)
  dev.set(d1)
  id <- identify(P[,1], P[,2], n=1)
  A  <<- c(A , id)
  VF <<- c(VF, id+1)
  VB <<- c(VB, id-1)
  do.plot()
  enable.widgets(TRUE)
}

remove.tear <- function(h, ...) {
  enable.widgets(FALSE)
  id <- identify(P[,1], P[,2], n=1)
  tid <- which(id==A)
  A  <<- A[-tid]
  VF <<- VF[-tid]
  VB <<- VB[-tid]
  do.plot()
  enable.widgets(TRUE)
}

move.point <- function(h, ...) {
  enable.widgets(FALSE)
  id <- identify(P[,1], P[,2], n=1)
  tid <- which(id==VF)
  if (length(tid)) {
    points(P[VF[tid],1], P[VF[tid],2], col="yellow")
    id <- identify(P[,1], P[,2], n=1)
    VF[tid] <<- id
  } 
  tid <- which(id==VB)
  if (length(tid)) {
    points(P[VB[tid],1], P[VB[tid],2], col="yellow")
    id <- identify(P[,1], P[,2], n=1)
    VB[tid] <<- id
  }
  tid <- which(id==A)
  if (length(tid)) {
    points(P[A[tid],1], P[A[tid],2], col="yellow")
    id <- identify(P[,1], P[,2], n=1)
    A[tid] <<- id
  }
  do.plot()
  enable.widgets(TRUE)
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
  segs <- map.to.segments(map)
  P <<- segments.to.outline(segs)
  do.plot()
}) 
tbl[3, 1, anchor = c(1, 0)] = "Mode"
tbl[3, 2, anchor = c(0, 0), expand = TRUE] <- g.add <- gbutton("Add tear",
                              handler=add.tear)
tbl[3, 3, anchor = c(0, 0), expand = TRUE] <- g.move <- gbutton("Move Point",
                              handler=move.point)
tbl[3, 4, anchor = c(0, 0), expand = TRUE] <- g.remove <- gbutton("Remove tear",
                              handler=remove.tear)

do.plot <- function() {
  dev.set(d1)
  if ("Sys" %in% svalue(g.show)) {
    plot.sys.map(sys, map)
  } else {
    plot.map(map)
  }

  if (length(A) > 0) {
    points(P[VF,], col="purple", pch="+")
    points(P[VB,], col="blue", pch="+")
    points(P[A, ], col="cyan", pch="+")
  }
}

tbl[4, 1, anchor = c(1, 0)] <- "Show"
tbl[4, 2, anchor = c(0, 0), expand = TRUE] <- g.show <- gcheckboxgroup(c("Sys"),
                              handler=function(h, ...) {
                                do.plot()
                              })
