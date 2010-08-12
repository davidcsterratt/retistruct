## Return initialised userdata list
initialise.userdata <- function() {
  A <<- c()          # Indices of apices of tears
  VB <<- c()         # Indices of forward verticies of tears
  VF <<- c()         # Indices of backward verticies of tears
  phi0 <<- 50        # Height of rim of retina in degrees
  r <<- NULL         # Reconstruction object
  iN <<- NA          # Index of nasal point
  iD <<- NA          # Index of dorsal point
  recfile.version <<- 1   # Version of reconstruction file data format
}

## Convenience functions for handlers
enable.group <- function(widgets, state=TRUE) {
  for (w in widgets) {
      enabled(w) <- state
    }
}

enable.widgets <- function(state) {
  enable.group(c(g.add, g.move, g.remove, g.reconstruct,
                 g.mark.n, g.mark.d,
                 g.phi0, g.show), state)
  if (is.na(iD) && is.na(iN)) {
    enable.group(c(g.reconstruct), FALSE)
  }
}

unsaved.data <- function(state) {
  enable.group(c(g.save), state)
}

## Function to report to set status
set.status <- function(s) {
  svalue(g.status) <- s
}

## Utility function
identify.abort.text <- function() {
  if (.Platform$GUI == "X11") {
    return("Right-click to abort.")
  }
  if (.Platform$GUI == "AQUA") {
    return("Press ESC to abort.")
  }
}

## Editting handlers
h.add <- function(h, ...) {
  unsaved.data(TRUE)
  enable.widgets(FALSE)
  svalue(g.status) <- paste("Click on the three points of the tear in any order.",
                            identify.abort.text())
  dev.set(d1)
  id <- identify(P[,1], P[,2], n=3)
  M <- markers.to.apex.vertices(id, gf, gb, P)
  A  <<- c(A,  M["A"])
  VF <<- c(VF, M["VF"])
  VB <<- c(VB, M["VB"])
  do.plot()
  svalue(g.status) <- ""
  enable.widgets(TRUE)
}

h.remove <- function(h, ...) {
  unsaved.data(TRUE)
  enable.widgets(FALSE)
  svalue(g.status) <- paste("Click on the apex of the tear to remvoe.",
                            identify.abort.text())
  dev.set(d1)
  id <- identify(P[,1], P[,2], n=1, plot=FALSE)
  tid <- which(id==A)
  if (length(tid) == 1) {
    A  <<- A[-tid]
    VF <<- VF[-tid]
    VB <<- VB[-tid]
  }
  do.plot()
  svalue(g.status) <- ""
  enable.widgets(TRUE)
}

## Handler for moving points
h.move <- function(h, ...) {
  unsaved.data(TRUE)
  enable.widgets(FALSE)
  dev.set(d1)
  ## Find the intial point
  svalue(g.status) <- paste("Click on apex or vertex to move.",
                            identify.abort.text())
  id <- identify(P[,1], P[,2], n=1, plot=FALSE)

  ## Locate tear in which th point occurs
  T <- cbind(A, VF, VB)                 # Tear matrix
  tid <- which(apply(id==T, 1, any))[1]
  pid <- which(id==T[tid,])[1]

  ## If there is a tear in which it occurs, select a point to move it to
  if (length(tid)) {
    svalue(g.status) <- paste("Click on point to move it to.",
                            identify.abort.text())
    points(P[T[tid,pid],1], P[T[tid,pid],2], col="yellow")
    id <- identify(P[,1], P[,2], n=1)
    if (length(id)) T[tid,pid] <- id
    ## It is possible to get the apex and vertex mixed up when moving points.
    ## Fix any errors.
    M <- markers.to.apex.vertices(T[tid,], gf, gb, P)
    A[tid]  <<- M["A"]
    VF[tid] <<- M["VF"]
    VB[tid] <<- M["VB"]
  }

  ## Display and cleanup
  do.plot()
  svalue(g.status) <- ""
  enable.widgets(TRUE)
}

## Handler for marking nasal point
h.mark.n <- function(h, ...) {
  unsaved.data(TRUE)
  enable.widgets(FALSE)
  svalue(g.status) <- paste("Click on nasal point.",
                            identify.abort.text())
  dev.set(d1)
  id <- identify(P[,1], P[,2], n=1)
  iN <<- id
  iD <<- NA
  do.plot()
  svalue(g.status) <- ""
  enable.widgets(TRUE)
}

## Handler for marking dorsal point
h.mark.d <- function(h, ...) {
  unsaved.data(TRUE)
  enable.widgets(FALSE)
  svalue(g.status) <- paste("Click on dorsal point.",
                            identify.abort.text())
  dev.set(d1)
  id <- identify(P[,1], P[,2], n=1)
  iD <<- id
  iN <<- NA
  do.plot()
  svalue(g.status) <- ""
  enable.widgets(TRUE)
}

## Handler for setting phi0
h.phi0 <- function(h, ...) {
  unsaved.data(TRUE)
  v <- svalue(g.phi0)
  if (v < -80) {
    v <- -89
  }
  if (v > 89) {
    v <- 89
  }
  phi0 <<- v
}

## Handler for saving state
h.save <- function(h, ...) {
  if (!is.null(dataset)) {
    write.csv(cbind(A, VB, VF), file.path(dataset, "T.csv"),  row.names=FALSE)
    write.csv(P, file.path(dataset, "P.csv"), row.names=FALSE)

    markup <- data.frame(iD=iD, iN=iN, phi0=phi0)    
    write.csv(markup, file.path(dataset, "markup.csv"))
    r$version <- recfile.version        # Datafile version
    if (!is.null(r)) {
      save(r, file=file.path(dataset, "r.Rdata"))
    }
  }
  unsaved.data(FALSE)
}

## Handler for brining up a file dialogue to open a dataset
## 
## Changes the following global variables:
##   dataset - directory in which data is contained
##   map     - the map data
##   Ds      - list of datapoints
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
          dataset <<- h$file
        })
  setwd(curdir)
  initialise.userdata()
  map <<- read.map(dataset)
  sys <<- read.sys(dataset)
  Ds <<- list(green=cbind(na.omit(sys[,'XGREEN']), na.omit(sys[,'YGREEN'])),
              red  =cbind(na.omit(sys[,'XRED'])  , na.omit(sys[,'YRED'])))
  segs <- map.to.segments(map)
  P <<- segments.to.outline(segs)
  svalue(g.dataset) <- dataset 
  
  ## Create forward and backward pointers
  t <- triangulate.outline(P, n=NA)
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
  markupfile <- file.path(dataset, "markup.csv")
  if (file.exists(markupfile)) {
    M <<- read.csv(markupfile)
    iD <<- M[1, "iD"]
    iN <<- M[1, "iN"]
    phi0 <<- M[1, "phi0"]
    svalue(g.phi0) <- phi0
  }

  ## The format in f.Rdata is deprecated
  foldfile <- file.path(dataset, "f.Rdata")
  if (file.exists(foldfile)) {
    file.remove(foldfile)
    gmessage("The algorithm has changed significantly since this retina was last reconstructed, so the cached reconstruction data has been deleted.",
             title="Warning", icon="warning")
  }
  recfile <- file.path(dataset, "r.Rdata")
  if (file.exists(recfile)) {
    load(recfile, globalenv())
    if (r$version != recfile.version) {
      gmessage("The algorithm has changed significantly since this retina was last reconstructed, so the cached reconstruction data has been deleted.",
               title="Warning", icon="warning")
      r <<- NULL
    }
  } else {
    r <<- NULL
  }
  unsaved.data(FALSE)
  enable.widgets(TRUE)
  dev.set(d2)
  plot.new()
  do.plot()
}

## Handler to start reconstructing the retina
h.reconstruct <- function(h, ...) {
  unsaved.data(TRUE)
  enable.widgets(FALSE)
  dev.set(d1)
  ct <- check.tears(cbind(A, VF, VB), gf, gb, P)
  if (length(ct)) {
    gmessage(paste("Invalid tears marked up. Somehow tears",
                   toString(ct),
                   "got messed up.  The red lines should always point in the clockwise direction and the orange ones in the anticlockwise direction. To fix, select \"Move Point\", and double click on tear",
                   toString(ct), "."), title="Error", icon="error")

    enable.widgets(TRUE)
    return()
  }
  i0 <- 0
  lambda0 <- 0
  if (!is.na(iD)) {
    i0 <- iD
    lambda0 <- 90
  }
  if (!is.na(iN)) {
    i0 <- iN
    lambda0 <- 0
  }
  r <<- fold.outline(P, cbind(A, VB, VF), phi0, i0=i0, lambda0=lambda0,
                     Ds=Ds,
                     graphical=TRUE, report=set.status)
  enable.widgets(TRUE)
  do.plot()
}

## Unused handlers
h.stitch.outline <- function(h, ...) {
  s <- stitch.outline(P, cbind(A, VB, VF))
  dev.set(d2)
  plot.stitch.flat(s)
}

h.triangulate.retina <- function(h, ...) {
  out <- triangulate.outline(P, n=400)
  dev.set(d2)
  with(out, trimesh(T, P))
}

## Handler for showing data
h.show <- function(h, ...) {
  do.plot()
}

## Plot in edit pane
do.plot <- function() {
  dev.set(d1)
  plot.outline(P, gb, axt="s")
  if ("Datapoints" %in% svalue(g.show)) {
    plot.datapoints.flat(Ds)
  }

  if ("Strain" %in% svalue(g.show)) {
    if (!is.null(r)) {
      dev.set(d1)
      plot.strain.flat(r)
      dev.set(d2)
      plot.l.vs.L(r)
      dev.set(d1)
    }
  } else {
    if (!is.null(r$Dss)) {
      dev.set(d2)
      plot.polar(r$phi0 * 180/pi)
      plot.datapoints.polar(r$Dss, cex=5)
      dev.set(d1)
    }
  }
  
  if ("Landmarks" %in% svalue(g.show)) {
    if (length(A) > 0) {
      points(P[VF,,drop=FALSE], col="red", pch="+")
      segments(P[A,1], P[A,2], P[VF,1], P[VF,2], col="red")
      points(P[VB,,drop=FALSE], col="orange", pch="+")
      segments(P[A,1], P[A,2], P[VB,1], P[VB,2], col="orange")
      points(P[A,,drop=FALSE], col="cyan", pch="+")
      text(P[A,,drop=FALSE]+100, labels=1:length(A), col="cyan")
    }
    if (!is.na(iD)) {
      text(P[iD,1], P[iD,2], "D")
    }
    if (!is.na(iN)) {
      text(P[iN,1], P[iN,2], "N")
    }
  }
  if ("Stitch" %in% svalue(g.show)) {
    if (!is.null(r$s)) {
      plot.stitch.flat(r$s, add=TRUE)
    }  
  }
  if ("Grid" %in% svalue(g.show)) {
    
    if (!is.null(r$t) && !is.null(r$r)) {
      with(r, plot.gridlines.flat(t$P, t$T, r$phi, r$lambda, m$Tt, p$phi0*180/pi))
    }
  }
}

## It would be nice to have error messages displayed graphically.
## This function should work, but has the problem that it always gives an
## "Error in get(\"toolkit\", inherits = TRUE) : object 'toolkit' not found\n"
## error itself.
h.error <- function() {
  gmessage(geterrmessage(), title="Error", icon="error")
}

retistruct <- function() {
  options(guiToolkit = "RGtk2")

  ## Global variables
  dataset <<- NULL                         # Directory of dataset
  initial.dir <<- "/afs/inf.ed.ac.uk/user/s/sterratt/projects/rettect/data/Anatomy/marked-up-retinae-2010-03-24/"

  initialise.userdata()

  ##
  ## GUI Layout
  ## 
  g.win <<- gwindow("Retistruct")

  g.rows <<- ggroup(horizontal=FALSE, container=g.win)
  ## Toolbar in row 1
  g.open         <<- gaction("Open", icon="open", handler=h.open)
  g.save         <<- gaction("Save", icon="save", handler=h.save)
  g.reconstruct  <<- gaction("Reconstuct retina", icon="polar", handler=h.reconstruct)
  g.toolbar <<- gtoolbar(list(open=g.open, save=g.save, reconstruct=g.reconstruct),
                         container=g.rows, style="both")

  ## Name of dataset in row 2
  g.dataset.row <<- ggroup(container=g.rows)
  g.dataset <<- glabel("No dataset selected", anchor=0, container=g.dataset.row)
  addSpring(g.dataset.row)
  ## enabled(g.dataset) <<- FALSE

  ## Body of interface
  g.body <<- ggroup(container=g.rows)

  ## Tear editor down left side
  g.editor <<- ggroup(horizontal = FALSE, container=g.body)

  g.add  <<-   gbutton("Add tear",    handler=h.add,    container=g.editor)
  g.move <<-   gbutton("Move Point",  handler=h.move,   container=g.editor)
  g.remove <<- gbutton("Remove tear", handler=h.remove, container=g.editor)
  g.mark.n <<- gbutton("Mark nasal",  handler=h.mark.n, container=g.editor)
  g.mark.d <<- gbutton("Mark dorsal", handler=h.mark.d, container=g.editor)

  ## Editing of phi0
  g.phi0.frame <<- gframe("Phi0", container=g.editor)
  g.phi0 <<- gedit(phi0, handler=h.phi0, width=5, coerce.with=as.numeric,
                   container=g.phi0.frame)

  ## What to show
  g.show.frame <<- gframe("Show", container=g.editor)
  g.show <<- gcheckboxgroup(c("Landmarks", "Stitch", "Grid", "Datapoints",
                              "Strain"),
                            checked=c(TRUE, FALSE, FALSE, FALSE),
                            handler=h.show, container=g.show.frame)

  ## Graphs at right
  g.f <<-  ggraphics(expand=TRUE, ps=11, container=g.body)
  d1 <<- dev.cur()
  g.f2 <<- ggraphics(expand=TRUE, ps=11, container=g.body)
  d2 <<- dev.cur()

  ## Status bar
  ## g.statusbar <<- ggroup(container=g.rows)
  g.statusbar <<- gframe("", expand=TRUE, container=g.rows)
  g.status <<- glabel("", container=g.statusbar)
  addSpring(g.statusbar)

  ## Disable buttons initially
  unsaved.data(FALSE)
  enable.widgets(FALSE)
}
