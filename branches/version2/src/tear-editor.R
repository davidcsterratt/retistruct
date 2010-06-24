source("datafile-utils.R")
source("fold-sphere2.R")
require(geometry)
source("triangle/triangle.R", chdir=TRUE)


## install.packages('gWidgetsRGtk2') first if not installed
if (!require("gWidgetsRGtk2")) install.packages("gWidgetsRGtk2")
if (!require("cairoDevice")) install.packages("cairoDevice")
library(gWidgetsRGtk2)
options(guiToolkit = "RGtk2")

## Global variables
dataset <- NULL
initial.dir <- "/afs/inf.ed.ac.uk/user/s/sterratt/projects/rettect/data/Anatomy/marked-up-retinae-2010-03-24/"
A <- c()
VB <- c()
VF <- c()

## Convenience functions for handlers
enable.group <- function(widgets, state=TRUE) {
  for (w in widgets) {
      enabled(w) <- state
    }
}

enable.widgets <- function(state) {
  enable.group(c(g.add, g.move, g.remove, g.fold, g.save), state)
}

## Editting handlers
h.add <- function(h, ...) {
  enable.widgets(FALSE)
  dev.set(d1)
  id <- identify(P[,1], P[,2], n=1)
  A  <<- c(A , id)
  VF <<- c(VF, gf[id])
  VB <<- c(VB, gb[id])
  do.plot()
  enable.widgets(TRUE)
}

h.remove <- function(h, ...) {
  enable.widgets(FALSE)
  dev.set(d1)
  id <- identify(P[,1], P[,2], n=1)
  print(id)
  tid <- which(id==A)
  print(tid)
  A  <<- A[-tid]
  VF <<- VF[-tid]
  VB <<- VB[-tid]
  do.plot()
  enable.widgets(TRUE)
}

h.move <- function(h, ...) {
  enable.widgets(FALSE)
  dev.set(d1)
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

## Handler for saving state
h.save <- function(h, ...) {
  datadir <- dataset
  write.csv(cbind(A, VB, VF), file.path(datadir, "T.csv"),  row.names=FALSE)
  write.csv(P, file.path(datadir, "P.csv"), row.names=FALSE)
}

## Handler for brining up a file dialogue to open a dataset
## 
## Changes the following global variables:
##   dataset - directory in which data is contained
##   map     - the map data
##   sys     - the sys data
##   P       - the outline
##   gf      - forward pointers
##   gb      - backward pointers
##   A       - tear apices
##   VF      - tear forward verticies
##   VB      - tear backward verticies
##
## Produces a plot of the retina in device d1
## 
h.open <- function(h, ...) {
  curdir <- getwd()
  if (is.null(dataset)) {
    info = file.info(initial.dir)
    if (!is.na(info$isdir)) {
      setwd(initial.dir)
    }
  } else {
    setwd(dataset)
    setwd("..")
  } 
  gfile(type="selectdir", text="Select a directory...",
        handler = function(h, ...) {
          print(h$file)
          dataset <<- h$file
        })
  setwd(curdir)
  map <<- read.map(dataset)
  sys <<- read.sys(dataset)
  segs <- map.to.segments(map)
  P <<- segments.to.outline(segs)
  svalue(g.dataset) <- dataset 
  
  ## Create forward and backward pointers
  t <- make.triangulation(P, n=NA)
  gf <<- t$gf
  gb <<- t$gb

  ## Read in tearfile
  tearfile <- file.path(dataset, "T.csv")
  if (file.exists(tearfile)) {
    T <- read.csv(tearfile)
    A  <<- T[,1]                            # apicies of tears
    VB <<- T[,2]                           # forward verticies
    VF <<- T[,3]                           # backward verticies
  }
  do.plot()
}

## Handler to start folding the outline
h.fold <- function(h, ...) {
  enable.widgets(FALSE)
  dev.set(d2)
  fold.retina(P, cbind(A, VB, VF), graphical=TRUE)
  enable.widgets(TRUE)
}

## Unused handlers
h.stitch.outline <- function(h, ...) {
  s <- stitch.outline(P, cbind(A, VB, VF))
  dev.set(d2)
  plot.stitch(s)
}

h.triangulate.retina <- function(h, ...) {
  out <- make.triangulation(P, n=400)
  dev.set(d2)
  with(out, trimesh(T, P))
}

## Plot in edit pane
do.plot <- function() {
  dev.set(d1)
  if ("Sys" %in% svalue(g.show)) {
    plot.sys.map(sys, map)
  } else {
    plot.map(map)
  }

  if (length(A) > 0) {
    points(P[VF,,drop=FALSE], col="red", pch="+")
    segments(P[A,1], P[A,2], P[VF,1], P[VF,2], col="red")
    points(P[VB,,drop=FALSE], col="orange", pch="+")
    segments(P[A,1], P[A,2], P[VB,1], P[VB,2], col="orange")
    points(P[A,,drop=FALSE], col="cyan", pch="+")
    text(P[A,,drop=FALSE]+100, labels=1:length(A), col="cyan")
  }
}

##
## GUI Layout
## 
tbl <- glayout(container = gwindow("Tear editor"), spacing = 0)

## Toolbar in row 1
tbl[1, 1, anchor = c(0, 0), expand = TRUE] <- g.open <- gbutton("Open",
                              handler=h.open)
tbl[1, 2, anchor = c(0, 0), expand = TRUE] <- g.save <- gbutton("Save",
                              handler=h.save)
tbl[1, 3, anchor = c(0, 0), expand = TRUE] <- g.fold <- gbutton("Fold retina",
                              handler=h.fold)

## Name of dataset in row 2
tbl[2, 1:5, anchor = c(0, 0), expand = TRUE] <- g.dataset <- gbutton("No datasetselected")
enabled(g.dataset) <- FALSE

## Tear editor down left side
tbl[3, 1, anchor = c(0, 0), expand = TRUE] <- g.add <- gbutton("Add tear",
                              handler=h.add)
tbl[4, 1, anchor = c(0, 0), expand = TRUE] <- g.move <- gbutton("Move Point",
                              handler=h.move)
tbl[5, 1, anchor = c(0, 0), expand = TRUE] <- g.remove <- gbutton("Remove tear",
                              handler=h.remove)
tbl[6, 1, anchor = c(1, 0)] <- "Show"
tbl[7, 1, anchor = c(0, 0), expand = TRUE] <- g.show <- gcheckboxgroup(c("Sys"),
                              handler=function(h, ...) {
                                do.plot()
                              })

## Graphs at right
tbl[3:20, 2:3, anchor = c(0, 0), expand = TRUE] = g.f = ggraphics(container = tbl,
    expand = TRUE, ps = 11)
d1 <- dev.cur()
tbl[3:20, 4:5, anchor = c(0, 0), expand = TRUE] = g.f2 = ggraphics(container = tbl,
    expand = TRUE, ps = 11)
d2 <- dev.cur()

