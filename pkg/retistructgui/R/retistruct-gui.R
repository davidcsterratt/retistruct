## Convenience functions for handlers
enable.group <- function(widgets, state=TRUE) {
  for (w in widgets) {
      enabled(w) <- state
    }
}

enable.widgets <- function(state) {
  enable.group(c(g.add, g.move, g.remove, g.reconstruct, g.properties,
                 g.mark.n, g.mark.d, g.mark.od,
                 g.phi0d, g.show, g.edit.show, g.data, g.eye,
                 g.print1, g.print2,
                 g.print.pdf1, g.print.pdf2), state)
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
  svalue(g.win)   <- paste(version.string(), "-" ,a$dataset)
  svalue(g.phi0d) <- a$phi0*180/pi
  svalue(g.eye)   <- a$side
  
  ## Read the reconstruction data
  withCallingHandlers({
    r <<- retistruct.read.recdata(a, check=TRUE)
  }, warning=h.warning, error=h.error)
  ## If there is no reconstruction data, show the markup so that we
  ## don't think there is no markup.
  if (is.null(r)) {
    svalue(g.show) <- unique(c("Markup", svalue(g.show)))
    svalue(g.nb) <- 1                   # Set "Edit" tab
  } else {
    svalue(g.nb) <- 2                   # Set "View" tab
  }
  delete(g.ids.frame, g.ids)
  ids <- getIDs(a)
  if (!is.null(ids)) {
    g.ids <<- gcheckboxgroup(ids, checked=rep(TRUE, length(ids)),
                             handler=h.show, container=g.ids.frame)
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
                                 dev.flat=d1, dev.polar=d2)
  }, warning=h.warning, error=h.warning)  
  enable.widgets(TRUE)
  do.plot()
}

## Handler for showing data
h.show <- function(h, ...) {
  if(!is.null(h$pageno)) {
    if (h$pageno == 1) {
      do.plot(markup=TRUE)
    } else {
      do.plot(markup=("Markup" %in% (svalue(g.show))))
    }
  } else {
    do.plot()
  }
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
print.bitmap <- function(d, file) {
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

## Print device d to file
print.pdf <- function(d, file) {
  dev.set(d)
  dev.print(pdf, file,
            width=getOption("retistruct.print.pdf.width"),
            height=getOption("retistruct.print.pdf.width"))
}

h.print.bitmap <- function(d, initialfilename) {
  curdir <- getwd()
  setwd(a$dataset)  
  gfile(type="save", text="Select a filename to save image to...",
        initialfilename=initialfilename,
        handler=function(h, ...) {
          print.bitmap(d, h$file)
        })
  setwd(curdir)
}

h.print.pdf <- function(d, initialfilename) {
  g.pdf <- gbasicdialog(title="PDF options", 
                        handler=function(h, ...) {
                          options(retistruct.print.pdf.width=svalue(g.width))
                        })

  g.pdf.group <<- ggroup(container=g.pdf, horizontal=TRUE)
  glabel("Width & height (inches)", container=g.pdf.group)
  g.width <<- gedit(getOption("retistruct.print.pdf.width"),
                    width=5, coerce.with=as.numeric,
                    container=g.pdf.group)
  visible(g.pdf, set=TRUE)

  curdir <- getwd()
  setwd(a$dataset)  
  gfile(type="save", text="Select a filename to save image to...",
        initialfilename=initialfilename,
        handler=function(h, ...) {
          print.pdf(d, h$file)
        })
  setwd(curdir)
}

h.print1 <- function(h, ...) {
  h.print.bitmap(d1, initialfilename="image-flat.png")
}

h.print.pdf1 <- function(h, ...) {
  h.print.pdf(d1, initialfilename="image-flat.pdf")
}

h.print2 <- function(h, ...) {
  h.print.bitmap(d2, initialfilename="image-polar.png")
}

h.print.pdf2 <- function(h, ...) {
  h.print.pdf(d2, initialfilename="image-polar.pdf")
}

getProjections <- function() {
  return(list("Azimuthal Equidistant"=azimuthal.equidistant,
              "Azimuthal Equal Area" =azimuthal.equalarea,
              "Azimuthal Conformal"  =azimuthal.conformal,
              "Sinusoidal"           =sinusoidal,
              "Orthographic"         =orthographic))
}

getTransforms <- function() {
    return(list("None"  =identity.transform,
                "Invert"=invert.sphere,
                "Invert to hemisphere"=invert.sphere.to.hemisphere))
}

## Plot in edit pane
do.plot <- function(markup=("Markup" %in% (svalue(g.show))) | (svalue(g.nb) == 1)) {
  
  if (is.null(r)) {
    r <- a
  }
  if ("Strain" %in% svalue(g.edit.show)) {   # Strain plot
    dev.set(d1)
    par(mar=c(0.5, 0.5, 0.5, 0.5))
    flatplot(r, axt="n",
             datapoints=FALSE,
             landmarks=FALSE,
             markup=FALSE,
             stitch=FALSE,
             grid=FALSE,
             mesh=FALSE,
             strain=TRUE,
             scalebar=1)
    dev.set(d2)
    par(mar=c(4.5, 4.5, 1, 0.5))
    lvsLplot(r)
  } else {
    dev.set(d1)
    par(mar=c(0.5, 0.5, 0.5, 0.5))
    flatplot(r, axt="n",
             datapoints=("Points" %in% svalue(g.show)),
             grouped=("Counts" %in% svalue(g.show)),
             landmarks=("Landmarks" %in% svalue(g.show)),
             markup=markup,
             stitch=("Stitch" %in% svalue(g.show)),
             grid=("Grid" %in% svalue(g.show)),
             ids=svalue(g.ids),
             mesh=FALSE,
             scalebar=1)
    dev.set(d2)
    par(mar=c(0.7, 0.7, 0.7, 0.7))
    projection(r,
               datapoints=("Points" %in% svalue(g.show)),
               datapoint.means=("Point means" %in% svalue(g.show)),
               landmarks=("Landmarks" %in% svalue(g.show)),
               transform=getTransforms()[[svalue(g.transform)]],
               projection=getProjections()[[svalue(g.projection)]],
               axisdir=cbind(phi=svalue(g.axis.el), lambda=svalue(g.axis.az)),
               proj.centre=cbind(phi=svalue(g.pc.el), lambda=svalue(g.pc.az)),
               datapoint.contours=("Point contours" %in% svalue(g.show)),
               grouped=("Counts" %in% svalue(g.show)),
               grouped.contours=("Count contours" %in% svalue(g.show)),
               ids=svalue(g.ids))
    ## FIXME: EOD not computed
    if (!is.null(r$EOD)) {
      polartext(paste("OD displacement:",
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


## Poperties dialogue
h.properties <- function(h, ...) {
  g.win <- gwindow("Properties",
                   parent=g.win)
  g.props <- ggroup(cont=g.win, horizontal=FALSE)
  g.colours <- gframe("Colours", container=g.props, horizontal=FALSE)
  
  g.prop.dl <- function(name, property, container) {
    g.prop.dl.group <- ggroup(container=container)
    glabel(name, container=g.prop.dl.group)
    g.dl <- gdroplist(palette(),
                      selected=which(palette() == options(property)),
                      container=g.prop.dl.group,
                      handler=function(h, ...) {
                        eval(parse(text=paste("options(", property, "=svalue(g.dl))")))
                        do.plot()})
  }

  g.prop.dl("Outline colour", "outline.col", g.colours)
  g.prop.dl("Stitch colour", "stitch.col", g.colours)
  g.prop.dl("Major gridline colour", "grid.maj.col", g.colours)
  g.prop.dl("Minor gridline colour", "grid.min.col", g.colours)

  gbutton("Close", container=g.props,
          handler = function(h,...) dispose(g.win))
}

version.string <- function() {
  return(paste("Retistruct ",
               packageDescription("retistruct", fields="Version"),
               " (Revision ", retistruct.global.revision, ")",
               sep=""))
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
  g.win <<- gwindow(version.string())

  g.rows <<- ggroup(horizontal=FALSE, container=g.win)
  ## Toolbar in row 1
  g.open         <<- gaction("Open", icon="open", handler=h.open)
  g.save         <<- gaction("Save", icon="save", handler=h.save)
  g.reconstruct  <<- gaction("Reconstuct retina", icon="polar", handler=h.reconstruct)
  g.properties   <<- gaction("Properties", icon="properties", handler=h.properties)
  g.toolbar <<- gtoolbar(list(open=g.open,
                              save=g.save,
                              reconstruct=g.reconstruct,
                              options=g.properties),
                         container=g.rows, style="both")

  ## Body of interface
  g.body <<- ggroup(container=g.rows)

  ## "Edit" and "View" tabs
  g.nb <<- gnotebook(container=g.body)

  ## Edit tab
  
  ## Tear editor
  g.editor <<- ggroup(horizontal = FALSE, container=g.nb, label="Edit")

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

  ## Whether to show strain
  g.edit.show.frame <<- gframe("Show", container=g.editor)
  g.edit.show <<- gcheckboxgroup(c("Strain"),
                                 checked=c(FALSE),
                                 handler=h.show, container=g.edit.show.frame)

  ## View Tab

  ## What to show
  g.view <<- ggroup(horizontal=FALSE, container=g.nb, label="View")
  g.show.frame <<- gframe("Show", container=g.view)
  g.show <<- gcheckboxgroup(c("Markup", "Stitch", "Grid", "Landmarks",
                              "Points", "Point means", "Point contours",
                              "Counts", "Count contours"),
                            checked=c(TRUE, FALSE, FALSE, FALSE,
                              FALSE, FALSE, FALSE,
                              FALSE, FALSE),
                            handler=h.show, container=g.show.frame)

  ## Group IDs
  g.ids.frame <<- gframe("IDs", container=g.view)
  g.ids <<- gcheckboxgroup("All", checked=TRUE,
                           handler=h.show, container=g.ids.frame)

  
  ## Projection type
  g.projection.frame <<- gframe("Projection", container=g.view)
  g.projection <<- gdroplist(names(getProjections()), selected=1,
                             handler=h.show, 
                             action=NULL, container=g.projection.frame)

  ## Projection centre
  g.pc.frame <<- gframe("Projection centre", container=g.view, horizontal=TRUE)
  glabel("El", container=g.pc.frame)
  g.pc.el <<- gedit("0", handler=h.show, width=5, coerce.with=as.numeric,
                      container=g.pc.frame)
  glabel("Az", container=g.pc.frame)
  g.pc.az <<- gedit("0", handler=h.show, width=5, coerce.with=as.numeric,
                      container=g.pc.frame)

  ## Transform
  g.transform.frame <<- gframe("Transform", container=g.view)
  g.transform <<- gdroplist(names(getTransforms()), selected = 1,  handler = h.show, 
                                 action = NULL, container = g.transform.frame)

  ## Axis direction
  g.axisdir.frame <<- gframe("Axis direction", container=g.view, horizontal=TRUE)
  glabel("El", container=g.axisdir.frame)
  g.axis.el <<- gedit("90", handler=h.show, width=5, coerce.with=as.numeric,
                      container=g.axisdir.frame)
  glabel("Az", container=g.axisdir.frame)
  g.axis.az <<- gedit("0", handler=h.show, width=5, coerce.with=as.numeric,
                      container=g.axisdir.frame)

  ## Graphs at right

  ## Flat plot
  g.f1 <<- ggroup(horizontal=FALSE, container=g.body)
  ## Buttons
  g.f1.buttons <<- ggroup(horizontal=TRUE, container=g.f1)
  g.print1     <<- gbutton("Bitmap", handler=h.print1,     container=g.f1.buttons)
  g.print.pdf1 <<- gbutton("PDF",    handler=h.print.pdf1, container=g.f1.buttons)
  ## Device itself
  g.fd1 <<- ggraphics(expand=TRUE, width=500, height=500, ps=11, container=g.f1)
  d1 <<- dev.cur()

  ## Projection
  g.f2 <<- ggroup(horizontal=FALSE, container=g.body)
  ## Buttons  
  g.f2.buttons <<- ggroup(horizontal=TRUE, container=g.f2)  
  g.print2     <<- gbutton("Bitmap", handler=h.print2,     container=g.f2.buttons)
  g.print.pdf2 <<- gbutton("PDF",    handler=h.print.pdf2, container=g.f2.buttons)
  ## Device itself
  g.fd2 <<- ggraphics(expand=TRUE, , width=500, height=500, ps=11, container=g.f2)
  d2 <<- dev.cur()
  
  ## Status bar
  ## g.statusbar <<- ggroup(container=g.rows)
  g.statusbar <<- gframe("", expand=TRUE, container=g.rows)
  g.status <<- glabel("", container=g.statusbar)
  addSpring(g.statusbar)

  ## Disable buttons initially
  unsaved.data(FALSE)
  enable.widgets(FALSE)

  ## Have to add the hander to the notebook at the end, otherwise
  ## there are complaints about various components not being defined.
  addHandlerChanged(g.nb, handler=h.show)
}

