retistruct.gui.revision <- function() {
  return(as.integer(gsub("Rev: ", "" ,gsub("\\$", "", "$Rev$"))))
}

## Convenience functions for handlers
enable.group <- function(widgets, state=TRUE) {
  for (w in widgets) {
      enabled(w) <- state
    }
}

enable.widgets <- function(state) {
  enable.group(c(g.add, g.move, g.remove, g.reconstruct,
                 g.mark.n, g.mark.d, g.mark.od,
                 g.phi0, g.show), state)
  if (!retistruct.potential.od()) {
    enable.group(c(g.mark.od), FALSE)
  }
  if (!retistruct.check.markup()) {
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
  V0 <<- c(V0, M["V0"])
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
  tid <- which(id==V0)
  if (length(tid) == 1) {
    V0 <<- V0[-tid]
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

  ## Locate tear in which the point occurs
  T <- cbind(V0, VF, VB)                 # Tear matrix
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
    V0[tid] <<- M["V0"]
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

## Handler for marking optic disc
h.mark.od <- function(h, ...) {
  unsaved.data(TRUE)
  enable.widgets(FALSE)
  svalue(g.status) <- paste("Click on a point on the optic disc.",
                            identify.abort.text())
  dev.set(d1)
  ## Convert list of segments to a matrix
  Sm <- NULL
  
  for (S in Ss) {
    Sm <- rbind(Sm, S)
  }
  id <- identify(Sm[,1], Sm[,2], n=1)
  N <- 0
  i <- 1
  while (id <= N && i<=length(Ss)) {
    N <- N + nrow(Ss[i])
    i <- i + 1
  }
  iOD <<- i
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
  retistruct.save.markup()
  retistruct.save.recdata()
  retistruct.export.matlab()
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
##   V0      - tear apices
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

  ## Read the raw data
  out <- try(retistruct.read.dataset(mess=gmessage))
  if (inherits(out, "try-error")) {
    gmessage(out, title="Error", icon="error")
  } else{
    ## Read the markup
    try(retistruct.read.markup(mess))
  
    ## Read the reconstruction data
    try(retistruct.read.recdata(mess=gmessage))

    svalue(g.dataset) <- dataset 
    svalue(g.phi0)    <- phi0
  
    unsaved.data(FALSE)
    enable.widgets(TRUE)
  }
  dev.set(d2)
  plot.new()
  do.plot()
}

## Handler to start reconstructing the retina
h.reconstruct <- function(h, ...) {
  unsaved.data(TRUE)
  enable.widgets(FALSE)
  out <- try(retistruct.reconstruct(mess=gmessage, report=set.status,
                        plot.3d=TRUE, dev.grid=d1, dev.polar=d2))
  if (inherits(out, "try-error")) {
    gmessage(out, title="Error", icon="error")
  }
  enable.widgets(TRUE)
  do.plot()
}

## Unused handlers
h.stitch.outline <- function(h, ...) {
  s <- stitch.outline(P, cbind(V0, VB, VF))
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
  plot.outline.flat(P, gb, axt="s")
  
  if ("Datapoints" %in% svalue(g.show)) {
    plot.datapoints.flat(Ds)
  }
  
  if ("Strain" %in% svalue(g.show)) {   # Strain plot
    if (!is.null(r)) {
      dev.set(d1)
      plot.strain.flat(r)
      dev.set(d2)
      plot.l.vs.L(r)
      dev.set(d1)
    }
  } else {                              # Polar plot
    dev.set(d2)
    plot.polar(phi0)
    if (!is.null(r$Dss)) {
      plot.outline.polar(r)
      if ("Datapoints" %in% svalue(g.show)) {
        plot.datapoints.polar(r$Dss, cex=5)
      }
    }
    if (!is.null(r$Sss) && ("Landmarks" %in% svalue(g.show))) {
      if (is.na(iOD)) {
         plot.landmarks.polar(r$Sss, col="orange")
      } else {
        plot.landmarks.polar(r$Sss[-iOD], col="orange")
        plot.landmarks.polar(r$Sss[iOD], col="blue")
      }
    }
    if (!is.null(r$EOD)) {
      text.polar(paste("OD displacement:", format(r$EOD, digits=3, nsmall=2), "deg"))
    }
    dev.set(d1)
  }
  
  if ("Markup" %in% svalue(g.show)) {
    if (length(V0) > 0) {
      points(P[VF,,drop=FALSE], col="red", pch="+")
      segments(P[V0,1], P[V0,2], P[VF,1], P[VF,2], col="red")
      points(P[VB,,drop=FALSE], col="orange", pch="+")
      segments(P[V0,1], P[V0,2], P[VB,1], P[VB,2], col="orange")
      points(P[V0,,drop=FALSE], col="cyan", pch="+")
      text(P[V0,,drop=FALSE]+100, labels=1:length(V0), col="cyan")
    }
    if (!is.na(iD)) {
      text(P[iD,1], P[iD,2], "D")
    }
    if (!is.na(iN)) {
      text(P[iN,1], P[iN,2], "N")
    }
  }

  if ("Stitch" %in% svalue(g.show)) {
    if (!is.null(r$gb)) {
      plot.stitch.flat(r, add=TRUE)
    }  
  }

  if ("Grid" %in% svalue(g.show)) {
    if (!is.null(r$phi)) {
      with(r, plot.gridlines.flat(P, T, phi, lambda, Tt, phi0*180/pi))
    }
  }

  if ("Landmarks" %in% svalue(g.show)) {
    if (is.na(iOD)) {
      plot.landmarks.flat(Ss, col="orange")
    } else {
      plot.landmarks.flat(Ss[-iOD], col="orange")
      plot.landmarks.flat(Ss[iOD], col="blue")
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

retistruct <- function(guiToolkit="RGtk2") {
  options(guiToolkit=guiToolkit)

  ## Global variables
  dataset <<- NULL                         # Directory of dataset
  initial.dir <<- "/afs/inf.ed.ac.uk/user/s/sterratt/projects/rettect/data/Anatomy/marked-up-retinae-2010-03-24/"

  retistruct.initialise.userdata()

  ##
  ## GUI Layout
  ## 
  g.win <<- gwindow(paste("Retistruct ",
                          packageDescription("retistruct", fields="Version"),
                          " (Revision", retistruct.global.revision(), ")",
                          sep=""))

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

  g.add     <<- gbutton("Add tear",    handler=h.add,     container=g.editor)
  g.move    <<- gbutton("Move Point",  handler=h.move,    container=g.editor)
  g.remove  <<- gbutton("Remove tear", handler=h.remove,  container=g.editor)
  g.mark.n  <<- gbutton("Mark nasal",  handler=h.mark.n,  container=g.editor)
  g.mark.d  <<- gbutton("Mark dorsal", handler=h.mark.d,  container=g.editor)
  g.mark.od <<- gbutton("Mark OD",     handler=h.mark.od, container=g.editor)
  ## Editing of phi0
  g.phi0.frame <<- gframe("Phi0", container=g.editor)
  g.phi0 <<- gedit(phi0, handler=h.phi0, width=5, coerce.with=as.numeric,
                   container=g.phi0.frame)

  ## What to show
  g.show.frame <<- gframe("Show", container=g.editor)
  g.show <<- gcheckboxgroup(c("Markup", "Stitch", "Grid", "Datapoints",
                              "Landmarks", "Strain"),
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
