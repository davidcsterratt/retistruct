## Convenience functions for handlers
enable.group <- function(widgets, state=TRUE) {
  for (w in widgets) {
      enabled(w) <- state
    }
}

enable.widgets <- function(state) {
  enable.group(c(g.add, g.move, g.remove, g.reconstruct,
                 g.mark.n, g.mark.d, g.mark.od,
                 g.phi0d, g.show, g.data, g.eye), state)
  enable.group(c(g.mark.od), retistruct.potential.od(r))
  if (!retistruct.check.markup(r)) {
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
  pids <- with(r, identify(P[,1], P[,2], n=3))
  tryCatch({
    r <<- addTear(r, pids)
  }, warning=h.warning, error=h.warning)
  do.plot()
  svalue(g.status) <- ""
  enable.widgets(TRUE)
}

## Handler for removing points
h.remove <- function(h, ...) {
  unsaved.data(TRUE)
  enable.widgets(FALSE)
  svalue(g.status) <- paste("Click on the apex of the tear to remvoe.",
                            identify.abort.text())
  dev.set(d1)
  id <- with(r, identify(P[,1], P[,2], n=1, plot=FALSE))
  r <<- removeTear(r, whichTear(r, id))
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
  id1 <- with(r, identify(P[,1], P[,2], n=1, plot=FALSE))
  
  ## Locate tear ID in which the point occurs
  tid <- whichTear(r, id1)

  ## If there is a tear in which it occurs, select a point to move it to
  if (!is.na(tid)) {
    svalue(g.status) <- paste("Click on point to move it to.",
                            identify.abort.text())

    ## Label first point
    with(r, points(P[id1,1], P[id1,2], col="yellow"))

    ## Select second point
    id2 <- with(r, identify(P[,1], P[,2], n=1))

    ## Get point ids of exsiting tear
    pids <- getTear(r, tid)

    ## Replace old point with desired new point
    if (length(id2)) pids[pids==id1] <- id2

    ## It is possible to get the apex and vertex mixed up when moving points.
    ## Fix any errors.
    pids <- labelTearPoints(r, pids)
    r <<- removeTear(r, tid)
    r <<- addTear(r, pids)
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
  id <- with(r, identify(P[,1], P[,2], n=1))
  tryCatch({
    r <<- setFixedPoint(r, id, "Nasal")
  }, warning=h.warning, error=h.warning)
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
  id <- with(r, identify(P[,1], P[,2], n=1))
  tryCatch({
    r <<- setFixedPoint(r, id, "Dorsal")
  }, warning=h.warning, error=h.warning)
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
  ## Convert list of segments to a matrix of points
  Sm <- NULL
  for (S in r$Ss) {
    Sm <- rbind(Sm, S)
  }

  ## Identify a point
  id <- identify(Sm[,1], Sm[,2], n=1)

  ## 
  N <- 0
  i <- 1
  while (id <= N && i<=length(Ss)) {
    N <- N + nrow(r$Ss[i])
    i <- i + 1
  }
  r <<- nameLandmark(r, i, "OD")
  do.plot()
  svalue(g.status) <- ""
  enable.widgets(TRUE)
}

## Handler for setting phi0d
h.phi0d <- function(h, ...) {
  unsaved.data(TRUE)
  v <- svalue(g.phi0d)
  if (v < -80) {
    v <- -89
  }
  if (v > 89) {
    v <- 89
  }
  r$phi0 <<- v*pi/180
}

## Handler for saving state
h.save <- function(h, ...) {
  retistruct.save.markup(r)
  retistruct.save.recdata(r)
  retistruct.export.matlab(r)
  unsaved.data(FALSE)
}

## Handler for brining up a file dialogue to open a dataset
##
## Changes the global objec r
##
## Produces a plot of the retina in device d1
## 
h.open <- function(h, ...) {
  curdir <- getwd()
  if (is.null(r$dataset)) {
    info = file.info(initial.dir)
    if (!is.na(info$isdir)) {
      setwd(initial.dir)
    }
  } else {
    setwd(r$dataset)
    setwd("..")
  } 
  gfile(type="selectdir", text="Select a directory...",
        handler = function(h, ...) {
          r$dataset <<- h$file
        })
  setwd(curdir)

  ## Read the raw data
  tryCatch({
    r <<- retistruct.read.dataset(r$dataset)
  }, warning=h.warning, error=h.error)
  
  ## Read the markup
  tryCatch({
    r <<- retistruct.read.markup(r, error=message)
  }, warning=h.warning, error=h.warning)
  
  ## Read the reconstruction data
  tryCatch({
    r <<- retistruct.read.recdata(r)
  }, warning=h.warning, error=h.error)

  svalue(g.dataset) <- r$dataset 
  svalue(g.phi0d)   <- r$phi0*180/pi
  
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
  tryCatch({
    r <<- retistruct.reconstruct(r, report=set.status,
                                 plot.3d=TRUE, dev.grid=d1, dev.polar=d2)
  }, warning=h.warning, error=h.error)  
  enable.widgets(TRUE)
  do.plot()
}

## Handler for showing data
h.show <- function(h, ...) {
  do.plot()
}

## Handler for flipping DV axis
h.flipdv <- function(h, ...) {
  r$DVflip <<- ("Flip DV" %in% svalue(g.data))
  do.plot()
}

## Handler for dealing with data
h.eye <- function(h, ...) {
  r$side <<- svalue(g.eye)
  do.plot()
}

## Plot in edit pane
do.plot <- function() {
  dev.set(d1)
  plot.flat(r, axt="s",
            datapoints=("Datapoints" %in% svalue(g.show)),
            landmarks=("Landmarks" %in% svalue(g.show)),
            markup=("Markup" %in% svalue(g.show)))
  with(r, {
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
      plot.polar(r, cex=5)
      ## if (!is.null(r$Dss)) {
      ##   plot.outline.polar(r)
      ##   if ("Datapoints" %in% svalue(g.show)) {
      ##     plot.datapoints.polar(r$Dss, r$D.cols, cex=5)
      ##   }
      ## }
      ## if (!is.null(r$Sss) && ("Landmarks" %in% svalue(g.show))) {
      ##   if (is.na(r$iOD)) {
      ##     plot.landmarks.polar(r$Sss, col="orange")
      ##   } else {
      ##     plot.landmarks.polar(r$Sss[-iOD], col="orange")
      ##     plot.landmarks.polar(r$Sss[iOD], col="blue")
      ##   }
      ## }
      if (!is.null(r$EOD)) {
        text.polar(paste("OD displacement:", format(r$EOD, digits=3, nsmall=2), "deg"))
      }
      dev.set(d1)
    }

    if ("Stitch" %in% svalue(g.show)) {
      if (!is.null(r$gb) && !is.null(r$TFset)) {
        plot.stitch.flat(r, add=TRUE)
      }  
    }

    if ("Grid" %in% svalue(g.show)) {
      if (!is.null(r$phi)) {
        with(r, plot.gridlines.flat(P, T, phi, lambda, Tt, phi0*180/pi))
      }
    }
  })
}

## It would be nice to have error messages displayed graphically.
## This function should work, but has the problem that it always gives an
## "Error in get(\"toolkit\", inherits = TRUE) : object 'toolkit' not found\n"
## error itself.
h.error <- function(e) {
  gmessage(e, title="Error", icon="error")
  stop(e)
}

## Warning message
h.warning <- function(e) {
  gmessage(e, title="Warning", icon="warning")
}

## The GUI itself
retistruct <- function(guiToolkit="RGtk2") {
  options(guiToolkit=guiToolkit)

  ## Global variables
  dataset <<- NULL                         # Directory of dataset
  initial.dir <<- "."

  r <<- NULL

  ##
  ## GUI Layout
  ## 
  g.win <<- gwindow(paste("Retistruct ",
                          packageDescription("retistruct", fields="Version"),
                          " (Revision ", retistruct.global.revision, ")",
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

  ## Editting of data
  g.data.frame <<- gframe("Data", container=g.editor, horizontal=FALSE)
  g.data <<- gcheckboxgroup(c("Flip DV"),
                            checked=c(FALSE),
                            handler=h.flipdv, container=g.data.frame)
  g.eye.frame <<- gframe("Eye", container=g.editor, horizontal=FALSE)
  g.eye <<- gradio(c("Right", "Left"),
                            checked=c(FALSE),
                            handler=h.eye, container=g.eye.frame)
  
  ## Editing of phi0
  g.phi0d.frame <<- gframe("Phi0", container=g.editor)
  g.phi0d <<- gedit(0, handler=h.phi0d, width=5, coerce.with=as.numeric,
                   container=g.phi0d.frame)

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
