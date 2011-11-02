## Convenience functions for handlers
enable.group <- function(widgets, state=TRUE) {
  for (w in widgets) {
      enabled(w) <- state
    }
}

enable.widgets <- function(state) {
  enable.group(c(g.add, g.move, g.remove, g.reconstruct,
                 g.mark.n, g.mark.d, g.mark.od,
                 g.phi0d, g.show, g.data, g.eye,
                 g.print1, g.print2), state)
  if (state) 
    enable.group(c(g.mark.od), retistruct.potential.od(a))
  if (!retistruct.check.markup(a)) {
    enable.group(c(g.reconstruct), FALSE)
  }
}

unsaved.data <- function(state) {
  if (state) {
    r <<- NULL
  }
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
  pids <- with(a, identify(P[,1], P[,2], n=3))
  withCallingHandlers({
    a <<- addTear(a, pids)
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
  id <- with(a, identify(P[,1], P[,2], n=1, plot=FALSE))
  a <<- removeTear(a, whichTear(a, id))
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
  id1 <- with(a, identify(P[,1], P[,2], n=1, plot=FALSE))
  
  ## Locate tear ID in which the point occurs
  tid <- whichTear(a, id1)

  ## If there is a tear in which it occurs, select a point to move it to
  if (!is.na(tid)) {
    svalue(g.status) <- paste("Click on point to move it to.",
                            identify.abort.text())

    ## Label first point
    with(a, points(P[id1,1], P[id1,2], col="yellow"))

    ## Select second point
    id2 <- with(a, identify(P[,1], P[,2], n=1))

    ## Get point ids of exsiting tear
    pids <- getTear(a, tid)

    ## Replace old point with desired new point
    if (length(id2)) pids[pids==id1] <- id2

    ## It is possible to get the apex and vertex mixed up when moving points.
    ## Fix any errors.
    pids <- labelTearPoints(a, pids)
    a <<- removeTear(a, tid)
    a <<- addTear(a, pids)
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
  id <- with(a, identify(P[,1], P[,2], n=1))
  withCallingHandlers({
    a <<- setFixedPoint(a, id, "Nasal")
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
  id <- with(a, identify(P[,1], P[,2], n=1))
  withCallingHandlers({
    a <<- setFixedPoint(a, id, "Dorsal")
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
  for (S in a$Ss) {
    Sm <- rbind(Sm, S)
  }

  ## Identify a point
  id <- identify(Sm[,1], Sm[,2], n=1)
  
  ## Idendify segment in which point appears
  i <- 0
  N <- 0
  while (id > N && i < length(a$Ss)) {
    i <- i + 1
    N <- N +  nrow(a$Ss[[i]])
  }
  ## Set "OD" landmark
  a <<- nameLandmark(a, i, "OD")
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
  a$phi0 <<- v*pi/180
}

## Handler for saving state
h.save <- function(h, ...) {
  retistruct.save.markup(a)
  ## If the reconstruction doesn't exist, remove the reconstruction
  ## file to ensure consistency
  if (is.null(r)) {
    unlink(file.path(a$dataset, "r.Rdata"))
  } else {
    retistruct.save.recdata(r)
  }
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
  if (is.null(a$dataset)) {
    info = file.info(initial.dir)
    if (!is.na(info$isdir)) {
      setwd(initial.dir)
    }
  } else {
    setwd(a$dataset)
    setwd("..")
  } 
  gfile(type="selectdir", text="Select a directory...",
        handler = function(h, ...) {
          a$dataset <<- h$file
        })
  setwd(curdir)

  ## Read the raw data
  withCallingHandlers({
    a <<- retistruct.read.dataset(a$dataset)
  }, warning=h.warning, error=h.error)
  
  ## Read the markup
  withCallingHandlers({
    a <<- retistruct.read.markup(a, error=message)
  }, warning=h.warning, error=h.warning)
  svalue(g.dataset) <- a$dataset 
  svalue(g.phi0d)   <- a$phi0*180/pi
  
  ## Read the reconstruction data
  withCallingHandlers({
    r <<- retistruct.read.recdata(a)
  }, warning=h.warning, error=h.error)
  ## If there is no reconstruction data, show the markup so that we
  ## don't think there is no markup.
  if (is.null(r)) {
    svalue(g.show) <- unique(c("Markup", svalue(g.show)))
  }

  unsaved.data(FALSE)
  enable.widgets(TRUE)
  do.plot()
}

## Handler to start reconstructing the retina
h.reconstruct <- function(h, ...) {
  unsaved.data(TRUE)
  enable.widgets(FALSE)
  withCallingHandlers({
    r <<- retistruct.reconstruct(a, report=set.status,
                                 plot.3d=getOption("show.sphere"),
                                 dev.grid=d1, dev.polar=d2)
  }, warning=h.warning, error=h.warning)  
  enable.widgets(TRUE)
  do.plot()
}

## Handler for showing data
h.show <- function(h, ...) {
  do.plot()
}

## Handler for flipping DV axis
h.flipdv <- function(h, ...) {
  unsaved.data(TRUE)
  a$DVflip <<- ("Flip DV" %in% svalue(g.data))
  do.plot()
}

## Handler for dealing with data
h.eye <- function(h, ...) {
  unsaved.data(TRUE)
  a$side <<- svalue(g.eye)
  do.plot()
}

## Print device d to file
print.image <- function(d, file) {
  dev <- NULL
  if (grepl(".png$", file, ignore.case=TRUE)) 
    dev <- png
  if (grepl(".jpeg$", file, ignore.case=TRUE) ||
      grepl(".jpg$", file, ignore.case=TRUE))
    dev <- jpeg
  if (grepl(".tif$", file, ignore.case=TRUE) ||
      grepl(".tiff$", file, ignore.case=TRUE))
    dev <- tiff
  if (is.null(dev)) {
    file <- paste(file, ".png")
    dev <- png
  }
  dev.set(d)
  dev.print(dev, file, width=1000, height=1000)
}

h.print.image <- function(d, initialfilename) {
  setwd(a$dataset)  
  gfile(type="save", text="Select a filename to save image to...",
        initialfilename=initialfilename,
        handler=function(h, ...) {
          print.image(d, h$file)
        })
}

h.print1 <- function(h, ...) {
  h.print.image(d1, initialfilename="image-flat.png")
}

h.print2 <- function(h, ...) {
  h.print.image(d2, initialfilename="image-polar.png")
}

## Plot in edit pane
do.plot <- function() {
  if (is.null(r)) {
    r <- a
  }
  if ("Strain" %in% svalue(g.show)) {   # Strain plot
    dev.set(d1)
    par(mar=c(0.5, 0.5, 0.5, 0.5))
    plot.flat(r, axt="n",
              datapoints=FALSE,
              landmarks=FALSE,
              markup=FALSE,
              stitch=FALSE,
              grid=FALSE,
              mesh=FALSE,
              strain=TRUE,
              scalebar=1)
    dev.set(d2)
    par(mar=c(4.5, 4.5, 0.5, 0.5))
    plot.l.vs.L(r)
  } else {
    dev.set(d1)
    par(mar=c(0.5, 0.5, 0.5, 0.5))
    plot.flat(r, axt="n",
              datapoints=("Datapoints" %in% svalue(g.show)),
              landmarks=("Landmarks" %in% svalue(g.show)),
              markup=("Markup" %in% svalue(g.show)),
              stitch=("Stitch" %in% svalue(g.show)),
              grid=("Grid" %in% svalue(g.show)),
              mesh=FALSE,
              scalebar=1)
    dev.set(d2)
    plot.polar(r,
               datapoints=("Datapoints" %in% svalue(g.show)),
               datapoint.means=("Means" %in% svalue(g.show)),
               landmarks=("Landmarks" %in% svalue(g.show)),
               preserve.area=("Preserve area" %in% svalue(g.show)),
               datapoint.contours=("Contours" %in% svalue(g.show)))
    ## FIXME: EOD not computed
    if (!is.null(r$EOD)) {
      text.polar(paste("OD displacement:",
                       format(r$EOD, digits=3, nsmall=2), "deg"))
    }
  }
  dev.set(d1)
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
  invokeRestart("muffleWarning")
}

## The GUI itself
retistruct <- function(guiToolkit="RGtk2") {
  options(guiToolkit=guiToolkit)

  ## Global variables
  dataset <<- NULL                         # Directory of dataset
  initial.dir <<- "."

  ## Annotation object
  a <<- NULL

  ## Reconstruction object
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
  g.show <<- gcheckboxgroup(c("Markup", "Stitch", "Grid", "Datapoints", "Means",
                              "Landmarks", "Strain", "Preserve area", "Contours"),
                            checked=c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
                            handler=h.show, container=g.show.frame)

  ## Graphs at right
  g.f1 <<- ggroup(horizontal = FALSE, container=g.body)
  g.fd1 <<- ggraphics(expand=TRUE, ps=11, container=g.f1)
  d1 <<- dev.cur()
  g.print1     <<- gbutton("Print", handler=h.print1, container=g.f1)

  g.f2 <<- ggroup(horizontal = FALSE, container=g.body)
  g.fd2 <<- ggraphics(expand=TRUE, ps=11, container=g.f2)
  d2 <<- dev.cur()
  g.print2     <<- gbutton("Print", handler=h.print2, container=g.f2)

  ## Status bar
  ## g.statusbar <<- ggroup(container=g.rows)
  g.statusbar <<- gframe("", expand=TRUE, container=g.rows)
  g.status <<- glabel("", container=g.statusbar)
  addSpring(g.statusbar)

  ## Disable buttons initially
  unsaved.data(FALSE)
  enable.widgets(FALSE)
}
