source("datafile-utils.R")
source("fold-sphere2.R")
require(geometry)
source("triangle/triangle.R", chdir=TRUE)


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

move.point <- function(h, ...) {
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

save.state <- function(h, ...) {
  datadir <- svalue(g.dataset)
  write.csv(cbind(A, VB, VF), file.path(datadir, "T.csv"),  row.names=FALSE)
  write.csv(P, file.path(datadir, "P.csv"), row.names=FALSE)
}

h.fold.retina <- function(h, ...) {
  dev.set(d2)
  fold.retina(P, cbind(A, VB, VF), graphical=TRUE)
}

h.stitch.retina <- function(h, ...) {
  s <- stitch.retina(P, cbind(A, VB, VF))
  dev.set(d2)
  plot.stitch(s)
}

h.triangulate.retina <- function(h, ...) {
  out <- triangulate(P)
  dev.set(d2)
  with(out, trimesh(T, Q))
}

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

fileChoose <- function(action="print", text = "Select a file...",
                       type="open", ...) {
  gfile(text=text, type=type, ..., action = action, handler =
        function(h,...) {
          do.call(h$action, list(h$file))
        })
}

dataset <- NULL
initial.dir <- "/afs/inf.ed.ac.uk/user/s/sterratt/projects/rettect/data/Anatomy/marked-up-retinae-2010-03-24/"

tbl <- glayout(container = gwindow("Tear editor"), spacing = 0)
tbl[1, 1:2, anchor = c(0, 0), expand = TRUE] = g.f = ggraphics(container = tbl,
    expand = TRUE, ps = 11)
d1 <- dev.cur()
tbl[1, 3:4, anchor = c(0, 0), expand = TRUE] = g.f2 = ggraphics(container = tbl,
    expand = TRUE, ps = 11)
d2 <- dev.cur()
tbl[2, 1,   anchor = c(1, 0)] = "Dataset"
tbl[2, 2:4, anchor = c(0, 0), expand = TRUE] <- g.dataset <- gbutton("Select dataset", handler = function(h, ...) {
  curdir <- getwd()
  if (svalue(h$obj) == "Select dataset") {
    info = file.info(initial.dir)
    if (!is.na(info$isdir)) {
      setwd(initial.dir)
    }
  } else {
    setwd(svalue(h$obj))
    setwd("..")
  } 
  gfile(type="selectdir", text="Select a directory...",
        handler = function(h, ...) {
          print(h$file)
          svalue(g.dataset) <- h$file
        })
  setwd(curdir)
  dataset <<- svalue(h$obj)
  map <<- read.map(dataset)
  sys <<- read.sys(dataset)
  segs <- map.to.segments(map)
  P <<- segments.to.outline(segs)
  tearfile <- file.path(dataset, "T.csv")
  if (file.exists(tearfile)) {
    T <- read.csv(tearfile)
    A  <<- T[,1]                            # apicies of tears
    VB <<- T[,2]                           # forward verticies
    VF <<- T[,3]                           # backward verticies
  }
  do.plot()
}) 
tbl[3, 1, anchor = c(1, 0)] = "Actions"
tbl[3, 2, anchor = c(0, 0), expand = TRUE] <- g.add <- gbutton("Add tear",
                              handler=add.tear)
tbl[4, 2, anchor = c(0, 0), expand = TRUE] <- g.move <- gbutton("Move Point",
                              handler=move.point)
tbl[5, 2, anchor = c(0, 0), expand = TRUE] <- g.remove <- gbutton("Remove tear",
                              handler=remove.tear)
tbl[6, 2, anchor = c(0, 0), expand = TRUE] <- g.fold <- gbutton("Stitch retina",
                              handler=h.stitch.retina)
tbl[6, 4, anchor = c(0, 0), expand = TRUE] <- g.fold <- gbutton("Fold retina",
                              handler=h.fold.retina)

tbl[7, 2, anchor = c(0, 0), expand = TRUE] <- g.save <- gbutton("Save",
                              handler=save.state)
tbl[7, 4, anchor = c(0, 0), expand = TRUE] <- g.triangulate <- gbutton("Triangulate retina",
                              handler=h.triangulate.retina)

tbl[3, 3, anchor = c(1, 0)] <- "Show"
tbl[3, 4, anchor = c(0, 0), expand = TRUE] <- g.show <- gcheckboxgroup(c("Sys"),
                              handler=function(h, ...) {
                                do.plot()
                              })

