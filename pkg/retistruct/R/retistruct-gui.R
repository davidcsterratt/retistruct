##' Start the Retistruct GUI
##' @seealso gWidgets2
##' @return Object with \code{getData()} method to return
##' reconstructed retina data and environment \code{this} which
##' contains variables in object.
##' @importFrom grDevices dev.cur dev.set dev.print
##' @importFrom graphics identify
##' @importFrom utils packageVersion
##' @export
retistruct <- function() {
  ## This function is essentially the constructor for a class. The
  ## environment 'this' contains all member variables and functions of
  ## the class. We do not need to refer to 'this' in the code, but we
  ## will return it to facilitate debugging. 
  this <- environment()

  ## Global variables
  dataset <- NULL                         # Directory of dataset
  initial.dir <- "."
  a <- NULL                             # Annotation object
  r <- NULL                             # Reconstruction object
  ## Path to extdata for demos
  extdata       <- file.path(system.file(package = "retistruct"), "extdata")
  extdata.demos <- file.path(system.file(package = "retistructdemos"), "extdata")
  
  ## Accessor functions
  getData <- function() {
    return(r)
  }

  ##
  ## Load and install graphics packages
  ##
  
  ## TODO: It would be good to make it possible to set
  ## guiToolkit=tcltk as an option to the function so that there is no
  ## need to install gWidgets2RGtk2, which is problematic on a Mac
  ##
  ## @param guiToolkit The toolkit \code{gWidgets2} toolkit that will
  ## be used for the GUI
  guiToolkit <- "RGtk2"
  require.package <- function(pkg) {
    suggests <- parse.dependencies(utils::installed.packages()["retistruct","Suggests"])
    suggests <- suggests[suggests[,1] == pkg]
    uptodate <- TRUE
    if (pkg %in% utils::installed.packages()[,"Package"]) {
      uptodate <- ifelse(is.na(suggests[2]),
                         TRUE,
                         eval(parse(text=paste("packageVersion(suggests[1])", suggests[2],  "suggests[3]"))))
    }
    if (!suppressWarnings(require(pkg, character.only=TRUE)) | 
        !uptodate) {
      message(paste("Trying to install required package", pkg))
      utils::install.packages(pkg)
      if (pkg=="RGtk2" & .Platform$OS.type=="windows") {
        require(pkg, character.only=TRUE)
        return(FALSE)
      }
      if (!suppressWarnings(require(pkg, character.only=TRUE))) {
        stop(paste("Could not install", pkg))
      }
    }
    return(TRUE)
  }
  if (!require.package(guiToolkit)) {
    return("To finish install, restart R and type retistruct::retistruct")
  }
  require.package(paste0("gWidgets2", guiToolkit))
  if ((guiToolkit == "RGtk2") &&
      (packageVersion("gWidgets2RGtk2") == "1.0.5")) {
    message("Error is likely with version 1.0.5 of gWidgets2RGtk2.")
    message("Try installing using install.packages(\"gWidgets2RGtk2\") .")
    message("You may quit and start R before doing this.")
  }
  require.package("cairoDevice")
  options(guiToolkit=guiToolkit)
  
  ##
  ## Convenience functions for handlers
  ## 

  enable.group <- function(widgets, state=TRUE) {
    for (w in widgets) {
      gWidgets2::enabled(w) <- state
    }
  }

  enable.widgets <- function(state) {
    enable.group(c(g.add, g.move, g.remove, g.reconstruct, g.properties,
                   g.mark.n, g.mark.d, g.mark.od,
                   g.phi0d, g.show, g.edit.show, g.data, g.eye,
                   g.print1, g.print2,
                   g.print.pdf1, g.print.pdf2), state)
    if (state) 
      enable.group(c(g.mark.od), length(a$getFeatureSet("LandmarkSet")$getIDs() > 0))
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
  set.status <- function(...) {
    ## FIXME #46: status is not work; retistruct::report() added to
    ## ensure output
    retistruct::report(...)
    gWidgets2::svalue(g.status) <- paste0(...)
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

  ##
  ## Editting handlers
  ##

  ## Handler for adding a tear
  h.add <- function(h, ...) {
    unsaved.data(TRUE)
    enable.widgets(FALSE)
    gWidgets2::svalue(g.status) <- paste("Click on the three points of the tear in any order.",
                              identify.abort.text())
    dev.set(d1)
    P <- a$getPoints()
    pids <- identify(P[,"X"], P[,"Y"], n=3, col=getOption("TF.col"))
    withCallingHandlers({
      a$addTear(pids)
    }, warning=h.warning, error=h.warning)
    do.plot()
    gWidgets2::svalue(g.status) <- ""
    enable.widgets(TRUE)
  }

  ## Handler for removing a tear
  h.remove <- function(h, ...) {
    unsaved.data(TRUE)
    enable.widgets(FALSE)
    gWidgets2::svalue(g.status) <- paste("Click on the apex of the tear to remvoe.",
                              identify.abort.text())
    dev.set(d1)
    P <- a$getPoints()
    id <- identify(P[,"X"], P[,"Y"], n=1, plot=FALSE)
    a$removeTear(a$whichTear(id))
    do.plot()
    gWidgets2::svalue(g.status) <- ""
    enable.widgets(TRUE)
  }

  ## Handler for moving a point in a tear
  h.move <- function(h, ...) {
    unsaved.data(TRUE)
    enable.widgets(FALSE)
    dev.set(d1)
    ## Find the intial point
    gWidgets2::svalue(g.status) <- paste("Click on apex or vertex to move.",
                                         identify.abort.text())
    P <- a$getPoints()
    id1 <- identify(P[,"X"], P[,"Y"], n=1, plot=FALSE)
    
    ## Locate tear ID in which the point occurs
    tid <- a$whichTear(id1)

    ## If there is a tear in which it occurs, select a point to move it to
    if (!is.na(tid)) {
      gWidgets2::svalue(g.status) <- paste("Click on point to move it to.",
                                           identify.abort.text())

      ## Label first point
      points(P[id1,1], P[id1,2], col="yellow")

      ## Select second point
      id2 <- identify(P[,"X"], P[,"Y"], n=1)

      ## Get point ids of exsiting tear
      pids <- a$getTear(tid)

      ## Replace old point with desired new point
      if (length(id2)) pids[pids==id1] <- id2

      ## It is possible to get the apex and vertex mixed up when moving points.
      ## Fix any errors.
      pids <- a$labelTearPoints(pids)
      a$removeTear(tid)
      a$addTear(pids)
    }

    ## Display and cleanup
    do.plot()
    gWidgets2::svalue(g.status) <- ""
    enable.widgets(TRUE)
  }

  ## Handler for marking nasal point
  h.mark.n <- function(h, ...) {
    unsaved.data(TRUE)
    enable.widgets(FALSE)
    gWidgets2::svalue(g.status) <- paste("Click on nasal point.",
                                         identify.abort.text())
    dev.set(d1)
    P <- a$getPoints()
    id <- identify(P[,"X"], P[,"Y"], n=1)
    withCallingHandlers({
      a$setFixedPoint(id, "Nasal")
    }, warning=h.warning, error=h.warning)
    do.plot()
    gWidgets2::svalue(g.status) <- ""
    enable.widgets(TRUE)
  }

  ## Handler for marking dorsal point
  h.mark.d <- function(h, ...) {
    unsaved.data(TRUE)
    enable.widgets(FALSE)
    gWidgets2::svalue(g.status) <- paste("Click on dorsal point.",
                              identify.abort.text())
    dev.set(d1)
    P <- a$getPoints()
    id <- identify(P[,"X"], P[,"Y"], n=1)
    withCallingHandlers({
      a$setFixedPoint(id, "Dorsal")
    }, warning=h.warning, error=h.warning)
    do.plot()
    gWidgets2::svalue(g.status) <- ""
    enable.widgets(TRUE)
  }

  ## Handler for marking optic disc
  h.mark.od <- function(h, ...) {
    unsaved.data(TRUE)
    enable.widgets(FALSE)
    gWidgets2::svalue(g.status) <- paste("Click on a point on the optic disc.",
                                         identify.abort.text())
    dev.set(d1)
    ## Convert list of segments to a matrix of points
    Sm <- NULL
    fs <- a$getFeatureSet("LandmarkSet")
    Ss <- fs$getFeatures()
    for (S in Ss) {
      Sm <- rbind(Sm, S)
    }

    ## Identify a point
    id <- identify(Sm[,1], Sm[,2], n=1)
    
    ## Identify segment in which point appears
    i <- 0
    N <- 0
    while (id > N && i < length(Ss)) {
      i <- i + 1
      N <- N +  nrow(Ss[[i]])
    }
    ## Set "OD" landmark
    fs$setID(i, "OD")

    ## Update IDs panel
    checked <- a$getIDs() %in% c(gWidgets2::svalue(g.ids), "OD")
    gWidgets2::delete(g.ids.frame, g.ids)
    ids <- a$getIDs()
    if (length(ids) > 0) {
      g.ids <<- gWidgets2::gcheckboxgroup(ids, checked=checked,
                                          handler=h.show,
                                          container=g.ids.frame)
    }
    
    do.plot()
    gWidgets2::svalue(g.status) <- ""
    enable.widgets(TRUE)
  }

  ## Handler for setting phi0d
  h.phi0d <- function(h, ...) {
    unsaved.data(TRUE)
    v <- gWidgets2::svalue(g.phi0d)
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
  ## Changes the object a
  ##
  ## Produces a plot of the retina in device d1
  ## 
  h.select <- function(h, ...) {
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
    dirname <- gWidgets2::gfile(type="selectdir", text="Select a directory...")
    if (length(dirname) > 0) {
      a$dataset <<- dirname
      h.open()
    }
    setwd(curdir)
  }

  ## Handler for opening a file
  h.open <- function(h, ...) {
    ## Read the raw data
    withCallingHandlers({
      a <<- retistruct.read.dataset(a$dataset, report=FALSE)
    }, warning=h.warning, error=h.error)
    
    ## Read the markup
    withCallingHandlers({
      a <<- retistruct.read.markup(a, error=message)
    }, warning=h.warning, error=h.warning)
    gWidgets2::svalue(g.win)   <- paste(version.string(), "-" ,a$dataset)
    gWidgets2::svalue(g.phi0d) <- a$phi0*180/pi
    gWidgets2::svalue(g.eye)   <- a$side
    gWidgets2::svalue(g.data) <-  ifelse (a$DVflip, "Flip DV", "")
    
    ## Read the reconstruction data
    withCallingHandlers({
      r <<- retistruct.read.recdata(a, check=TRUE)
    }, warning=h.warning, error=h.error)
    ## If there is no reconstruction data, show the markup so that we
    ## don't think there is no markup.
    if (is.null(r)) {
      gWidgets2::svalue(g.show) <- unique(c("Markup", gWidgets2::svalue(g.show)))
      gWidgets2::svalue(g.nb) <- 1                   # Set "Edit" tab
    } else {
      gWidgets2::svalue(g.nb) <- 2                   # Set "View" tab
    }
    gWidgets2::delete(g.ids.frame, g.ids)
    ids <- a$getIDs()
    if (length(ids) > 0) {
      g.ids <<- gWidgets2::gcheckboxgroup(ids, checked=rep(TRUE, length(ids)),
                                          handler=h.show, container=g.ids.frame)
    }
    unsaved.data(FALSE)
    enable.widgets(TRUE)
    do.plot()
  }

  h.retistructdemos <- function(dir1="Figure_6-data", dir2="left-ipsi", ...) {
    dataset <- file.path(extdata.demos, dir1, dir2)
    print(dataset)
    if (!file.exists(dataset)) {
      gWidgets2::gmessage("Install the retistructdemos package using devtools::install_github(\"davidcsterratt/retistruct/pkg/retistructdemos\")", title="Demo not installed", icon="error")
    } else {
      a$dataset <<- dataset
      h.open()
    }
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
    set.status("")
  }

  ## Handler for showing data
  h.show <- function(h, ...) {
    tab <<- "Edit"
    if (!is.null(h$page.no)) {
      ## h$page.no is set when the notebook handler is called
      ## print(paste("h.show: pageno:", h$page.no))
      if (h$page.no == 2) { tab <<- "View" }
    } else {
      ## Notebook not active
      ## print(paste("h.show: g.nb:", gWidgets2::svalue(g.nb)))
      if (gWidgets2::svalue(g.nb) == 2) { tab <<- "View" }
    }
    ## print(paste(tab, "Tab"))
    if (tab == "Edit") {
      do.plot(markup=TRUE)
    } else {
      do.plot(markup=("Markup" %in% (gWidgets2::svalue(g.show))))
    }
  }

  ## Handler for flipping DV axis
  h.flipdv <- function(h, ...) {
    unsaved.data(TRUE)
    a$DVflip <<- ("Flip DV" %in% gWidgets2::svalue(g.data))
    do.plot()
  }

  ## Handler for dealing with data
  h.eye <- function(h, ...) {
    unsaved.data(TRUE)
    a$side <<- gWidgets2::svalue(g.eye)
    do.plot()
  }

  ## Print device or function d to file
  print.bitmap <- function(d, file) {
    dev <- NULL
    if (grepl(".png$", file, ignore.case=TRUE)) 
      dev <- grDevices::png
    if (grepl(".jpeg$", file, ignore.case=TRUE) ||
        grepl(".jpg$", file, ignore.case=TRUE))
      dev <- grDevices::jpeg
    if (grepl(".tif$", file, ignore.case=TRUE) ||
        grepl(".tiff$", file, ignore.case=TRUE))
      dev <- grDevices::tiff
    if (is.null(dev)) {
      file <- paste(file, ".png")
      dev <- grDevices::png
    }
    if (is.function(d)) {
      dev(file, width=getOption("max.proj.dim"), height=getOption("max.proj.dim"))
      d(max.proj.dim=getOption("max.proj.dim"))
      dev.off()
    } else {
      dev.set(d)
      dev.print(dev, file, width=getOption("max.proj.dim"), height=getOption("max.proj.dim"))
    }
  }

  ## Print device d to file
  print.pdf <- function(d, file) {
    dev.set(d)
    dev.print(grDevices::pdf, file,
              width=getOption("retistruct.print.pdf.width"),
              height=getOption("retistruct.print.pdf.width"))
  }

  ## Handler for printing to a bitmap
  h.print.bitmap <- function(d, initial.filename) {
    curdir <- getwd()
    setwd(a$dataset)
    fname <- gWidgets2::gfile(type="save",
                              text="Select a filename to save image to...",
                              initial.filename=initial.filename)
    if (length(fname) > 0) {
      print.bitmap(d, fname)
    }
    setwd(curdir)
  }

  ## Handler for printing to a pdf file
  h.print.pdf <- function(d, initial.filename) {
    h.width <- function(h, ...) {
      options(retistruct.print.pdf.width=gWidgets2::svalue(g.width))
    }
    g.pdf <- gWidgets2::gbasicdialog(title="PDF options", 
                          handler=h.width)
    g.pdf.group <- gWidgets2::ggroup(container=g.pdf, horizontal=TRUE)
    gWidgets2::glabel("Width & height (inches)", container=g.pdf.group)
    g.width <- gWidgets2::gedit(getOption("retistruct.print.pdf.width"),
                     width=5, coerce.with=as.numeric,
                     container=g.pdf.group)


    gWidgets2::visible(g.pdf, set=TRUE)

    curdir <- getwd()
    setwd(a$dataset)  
    fname <- gWidgets2::gfile(type="save",
                              text="Select a filename to save image to...",
                              initial.filename=initial.filename)
    if (length(fname) > 0) {
      print.pdf(d, fname)
    }
    setwd(curdir)
  }

  ## Handlers for printing bitmaps and PDFs from the two graphics
  ## devices
  h.print1 <- function(h, ...) {
    h.print.bitmap(d1, initial.filename="image-flat.png")
  }

  h.print.pdf1 <- function(h, ...) {
    h.print.pdf(d1, initial.filename="image-flat.pdf")
  }

  h.print2 <- function(h, ...) {
    h.print.bitmap(plotProjection, initial.filename="image-polar.png")
  }

  h.print.pdf2 <- function(h, ...) {
    h.print.pdf(d2, initial.filename="image-polar.pdf")
  }

  ## Get and name available projections
  getProjections <- function() {
    return(list("Azimuthal Equidistant"=azimuthal.equidistant,
                "Azimuthal Equal Area" =azimuthal.equalarea,
                "Azimuthal Conformal"  =azimuthal.conformal,
                "Sinusoidal"           =sinusoidal,
                "Orthographic"         =orthographic))
  }

  ## Get and name available transforms
  getTransforms <- function() {
    return(list("None"  =identity.transform,
                "Invert"=invert.sphere,
                "Invert to hemisphere"=invert.sphere.to.hemisphere))
  }

  plotProjection <- function(max.proj.dim=getOption("max.proj.dim"),
                             markup=("Markup" %in% (gWidgets2::svalue(g.show)))) {
    projection(r,
               datapoints=("Points" %in% gWidgets2::svalue(g.show)),
               datapoint.means=("Point means" %in% gWidgets2::svalue(g.show)),
               landmarks=("Landmarks" %in% gWidgets2::svalue(g.show)),
               transform=getTransforms()[[gWidgets2::svalue(g.transform)]],
               projection=getProjections()[[gWidgets2::svalue(g.projection)]],
               axisdir=cbind(phi=gWidgets2::svalue(g.axis.el), lambda=gWidgets2::svalue(g.axis.az)),
               proj.centre=cbind(phi=gWidgets2::svalue(g.pc.el), lambda=gWidgets2::svalue(g.pc.az)),
               datapoint.contours=("Point contours" %in% gWidgets2::svalue(g.show)),
               grouped=("Counts" %in% gWidgets2::svalue(g.show)),
               grouped.contours=("Count contours" %in% gWidgets2::svalue(g.show)),
               grid=("Grid" %in% gWidgets2::svalue(g.show)),
               markup=markup,
               ids=gWidgets2::svalue(g.ids),
               max.proj.dim=max.proj.dim)
    
    if (!is.null(r$EOD)) {
      polartext(paste("OD displacement:",
                      format(r$EOD, digits=3, nsmall=2), "deg"))
    }
  }  
  
  ## Plot in edit pane
  do.plot <- function(markup=("Markup" %in% (gWidgets2::svalue(g.show)))) {
    
    if (is.null(r)) {
      r <- a
    }
    if (("Strain" %in% gWidgets2::svalue(g.edit.show)) & (tab == "Edit")) {   # Strain plot
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
      sphericalplot(r, strain=TRUE, datapoints=FALSE)
    } else {
      dev.set(d1)
      par(mar=c(0.5, 0.5, 0.5, 0.5))
      flatplot(r, axt="n",
               datapoints=("Points" %in% gWidgets2::svalue(g.show)),
               grouped=("Counts" %in% gWidgets2::svalue(g.show)),
               landmarks=("Landmarks" %in% gWidgets2::svalue(g.show)),
               markup=markup,
               stitch=("Stitch" %in% gWidgets2::svalue(g.show)),
               grid=("Grid" %in% gWidgets2::svalue(g.show)),
               ids=gWidgets2::svalue(g.ids),
               mesh=FALSE,
               scalebar=1)
      dev.set(d2)
      par(mar=c(0.7, 0.7, 0.7, 0.7))
      plotProjection(max.proj.dim=400, markup=markup)
      sphericalplot(r, datapoints=("Points" %in% gWidgets2::svalue(g.show)))
    }
    dev.set(d1)
    set.status("")
  }

  ## It would be nice to have error messages displayed graphically.
  ## This function should work, but has the problem that it always gives an
  ## "Error in get(\"toolkit\", inherits = TRUE) : object 'toolkit' not found\n"
  ## error itself.
  h.error <- function(e) {
    gWidgets2::gmessage(e, title="Error", icon="error")
    stop(e)
  }

  ## Warning message
  prior.warnings <- c()
  h.warning <- function(e) {
    if (!(any(e %in% prior.warnings))) {
      gWidgets2::gmessage(e, title="Warning", icon="warning")
      prior.warnings <<- c(prior.warnings, e)
    }
  }

  ## Poperties dialogue
  h.properties <- function(h, ...) {
    g.win <- gWidgets2::gwindow("Properties",
                     parent=g.win)
    g.props <- gWidgets2::ggroup(container=g.win, horizontal=FALSE)
    g.colours <- gWidgets2::gframe("Colours", container=g.props, horizontal=FALSE)
    cols <- c("black", "red", "green3", "blue", "cyan", "magenta", "yellow",
              "gray")
    g.prop.dl <- function(name, property, container) {
      g.prop.dl.group <- gWidgets2::ggroup(container=container)
      gWidgets2::glabel(name, container=g.prop.dl.group)
      g.dl <- gWidgets2::gcombobox(
        cols,
        selected=which(cols == options(property[1])),
        container=g.prop.dl.group,
        expand=TRUE,
        handler=function(h, ...) {
          for (p in property) {
            eval(parse(text=paste("options(", p, "=gWidgets2::svalue(g.dl))")))
          }
          do.plot()})
    }

    g.prop.dl("Outline colour", "outline.col", g.colours)
    g.prop.dl("Tear colour", c("TF.col", "TB.col", "V.col"), g.colours)
    g.prop.dl("Stitch colour", "stitch.col", g.colours)
    g.prop.dl("Major gridline colour", "grid.maj.col", g.colours)
    g.prop.dl("Minor gridline colour", "grid.min.col", g.colours)

    g.printing <- gWidgets2::gframe("Bitmap/PDF output", container=g.props, horizontal=FALSE)
    g.group <- gWidgets2::ggroup(container=g.printing)
    gWidgets2::glabel("Maximum width of projection", container=g.group)
    property <- "max.proj.dim"
    g.max.proj.dim <- gWidgets2::gedit(
                                   getOption(property),
                                   width=5, coerce.with=as.numeric,
                                   container=g.group)
    h.max.proj.dim <- function(h, ...) {
      eval(parse(text=paste("options(", property, "=gWidgets2::svalue(g.max.proj.dim))")))
    }
    gWidgets2::addHandlerKeystroke(g.max.proj.dim, handler=h.max.proj.dim)
    gWidgets2::addHandlerBlur(g.max.proj.dim, handler=h.max.proj.dim)

    gWidgets2::glabel("pixels", container=g.group)
    ## gWidgets2::svalue(g.max.proj.dim) <- getOption(property)
    
    gWidgets2::gbutton("Close", container=g.props,
            handler = function(h,...) gWidgets2::dispose(g.win))
  }

  ## Construct the version string
  version.string <- function() {
    return(paste0("Retistruct ",
                  utils::packageDescription("retistruct", fields="Version"),
                  " (", utils::packageDescription("retistruct", fields="Date"), ")"))
  }
  
  ##
  ## GUI Layout
  ##

  ## Top level window
  g.win <- gWidgets2::gwindow(version.string())

  ## Menu in row 0
  mbl <- list()
  mbl$File$Open <- gWidgets2::gaction(label="Open", handler=h.select)
  mbl$File$Save <- gWidgets2::gaction(label="Save", handler=h.save)
  mbl$Edit$Reconstruct <- gWidgets2::gaction(label="Reconstruct retina", handler=h.reconstruct)
  mbl$Edit$Properties <- gWidgets2::gaction(label="Properties", handler=h.properties)
  mbl$Demos$fig1 <-
    gWidgets2::gaction(label="Figure 1",
            handler=function(h, ...) {
              a$dataset <<- file.path(extdata, "GM509", "R-CONTRA")
              h.open()
            })
  mbl$Demos$fig2low <-
    gWidgets2::gaction(label="Figure 2A-D: Low deformation",
            handler=function(h, ...) {
              a$dataset <<- file.path(extdata, "GMB530", "R-CONTRA")
              h.open()
            })
  mbl$Demos$fig2high <-
    gWidgets2::gaction(label="Figure 2E-H: High deformation",
            handler=function(h, ...) {
              a$dataset <<- file.path(extdata, "GM182-4", "R-CONTRA")
              h.open()
            })
  mbl$Demos$smi32 <- gWidgets2::gaction(
    label="SMI32",
    handler=function(h, ...) {
      a$dataset <<- file.path(extdata, "smi32")
      h.open()
    })
  mbl$Demos$left.contra <-
    gWidgets2::gaction(label="Figure 6 Left Contra",
                       handler=function(h, ...) {
                         h.retistructdemos("Figure_6-data", "left-contra")
                       })
  mbl$Demos$left.ipsi <-
    gWidgets2::gaction(label="Figure 6 Left Ipsi",
                       handler=function(h, ...) {
                         h.retistructdemos("Figure_6-data", "left-ipsi")
                       })
  mbl$Demos$right.contra <-
    gWidgets2::gaction(label="Figure 6 Right Contra",
                       handler=function(h, ...) {
                         h.retistructdemos("Figure_6-data", "right-contra")
                       })
  mbl$Demos$right.ipsi <-
    gWidgets2::gaction(label="Figure 6 Right Ipsi",
                       handler=function(h, ...) {
                         h.retistructdemos("Figure_6-data", "right-ipsi")
                       })
  mbl$Help$About <- gWidgets2::gaction(
    label="About",
    handler=function(h, ...) {
      gWidgets2::gmessage(
        "Retistruct was written by David Sterratt at the University of Edinburgh, and tested by Daniel Lyngholm and Ian Thompson at the MRC Centre for Developmental Neurobiology, KCL.

This work was supported by a Programme Grant from the Wellcome Trust (G083305). ",
        title="About")
  })

  g.menu <- gWidgets2::gmenu(mbl, container=g.win)
  
  g.rows <- gWidgets2::ggroup(horizontal=FALSE, container=g.win)
  ## Toolbar in row 1
  g.open         <- gWidgets2::gaction("Open", icon="open", handler=h.select)
  g.save         <- gWidgets2::gaction("Save", icon="save", handler=h.save)
  g.reconstruct  <- gWidgets2::gaction("Reconstuct retina", icon="execute", handler=h.reconstruct)
  g.properties   <- gWidgets2::gaction("Properties", icon="properties", handler=h.properties)
  g.toolbar <- gWidgets2::gtoolbar(list(open=g.open,
                             save=g.save,
                             reconstruct=g.reconstruct,
                             options=g.properties),
                        container=g.win, style="both")

  ## Body of interface
  g.body <- gWidgets2::ggroup(container=g.rows, expand=TRUE)

  ## "Edit" and "View" tabs
  g.nb <- gWidgets2::gnotebook(container=g.body)
  ## Keep track of state in global variable, which can be "Edit" or "View"
  tab <- "Edit"

  ## Edit tab
  
  ## Tear editor
  g.editor <- gWidgets2::ggroup(horizontal = FALSE, container=g.nb, label="Edit")

  g.add     <- gWidgets2::gbutton("Add tear",    handler=h.add,     container=g.editor)
  g.move    <- gWidgets2::gbutton("Move Point",  handler=h.move,    container=g.editor)
  g.remove  <- gWidgets2::gbutton("Remove tear", handler=h.remove,  container=g.editor)
  g.mark.n  <- gWidgets2::gbutton("Mark nasal",  handler=h.mark.n,  container=g.editor)
  g.mark.d  <- gWidgets2::gbutton("Mark dorsal", handler=h.mark.d,  container=g.editor)
  g.mark.od <- gWidgets2::gbutton("Mark OD",     handler=h.mark.od, container=g.editor)
  
  ## Editting of data
  g.data.frame <- gWidgets2::gframe("Data", container=g.editor, horizontal=FALSE)
  g.data <- gWidgets2::gcheckboxgroup(c("Flip DV"),
                           checked=c(FALSE),
                           handler=h.flipdv, container=g.data.frame)
  g.eye.frame <- gWidgets2::gframe("Eye", container=g.editor, horizontal=FALSE)
  g.eye <- gWidgets2::gradio(c("Right", "Left"),
                  checked=c(FALSE),
                  handler=h.eye, container=g.eye.frame)
  
  ## Editing of phi0
  g.phi0d.frame <- gWidgets2::gframe("Phi0", container=g.editor)
  g.phi0d <- gWidgets2::gedit("0", handler=h.phi0d, width=5, coerce.with=as.numeric,
                              container=g.phi0d.frame)

  ## Whether to show strain
  g.edit.show.frame <- gWidgets2::gframe("Show", container=g.editor)
  g.edit.show <- gWidgets2::gcheckboxgroup(c("Strain"),
                                checked=c(FALSE),
                                handler=h.show, container=g.edit.show.frame)

  ## View Tab

  ## What to show
  g.view <- gWidgets2::ggroup(horizontal=TRUE, container=g.nb, label="View")
  g.view1 <- gWidgets2::ggroup(horizontal=FALSE, container=g.view)

  g.show.frame <- gWidgets2::gframe("Show", container=g.view1)
  g.show <- gWidgets2::gcheckboxgroup(c("Markup", "Stitch", "Grid", "Landmarks",
                             "Points", "Point means", "Point contours",
                             "Counts", "Count contours"),
                           checked=c(TRUE, FALSE, FALSE, TRUE,
                             TRUE, FALSE, FALSE,
                             TRUE, FALSE),
                           handler=h.show, container=g.show.frame)

  ## Group IDs
  g.ids.frame <- gWidgets2::gframe("IDs", container=g.view1)
  g.ids <- gWidgets2::gcheckboxgroup("All", checked=TRUE,
                          handler=h.show, container=g.ids.frame)


  g.view2 <- gWidgets2::ggroup(horizontal=FALSE, container=g.view)
  ## Projection type
  g.projection.frame <- gWidgets2::gframe("Projection", container=g.view2)
  g.projection <- gWidgets2::gcombobox(
    names(getProjections()), selected=1, editable=FALSE,
    handler=h.show, action=NULL,
    ellipsize="none",
    container=g.projection.frame)

  ## Projection centre
  g.pc.frame <- gWidgets2::gframe("Projection centre", container=g.view2, horizontal=TRUE)
  gWidgets2::glabel("El", container=g.pc.frame)
  g.pc.el <- gWidgets2::gedit("0", handler=h.show, width=5, coerce.with=as.numeric,
                   container=g.pc.frame)
  gWidgets2::glabel("Az", container=g.pc.frame)
  g.pc.az <- gWidgets2::gedit("0", handler=h.show, width=5, coerce.with=as.numeric,
                   container=g.pc.frame)

  ## Transform
  g.transform.frame <- gWidgets2::gframe("Transform", container=g.view2)
  g.transform <- gWidgets2::gcombobox(
    names(getTransforms()), selected = 1, editable=FALSE,
    ellipsize="none",
    handler = h.show, action = NULL, 
    container = g.transform.frame)

  ## Axis direction
  g.axisdir.frame <- gWidgets2::gframe("Axis direction", container=g.view2, horizontal=TRUE)
  gWidgets2::glabel("El", container=g.axisdir.frame)
  g.axis.el <- gWidgets2::gedit("90", handler=h.show, width=5, coerce.with=as.numeric,
                     container=g.axisdir.frame)
  gWidgets2::glabel("Az", container=g.axisdir.frame)
  g.axis.az <- gWidgets2::gedit("0", handler=h.show, width=5, coerce.with=as.numeric,
                     container=g.axisdir.frame)

  ## Graphs at right

  ## Flat plot
  g.f1 <- gWidgets2::ggroup(horizontal=FALSE, container=g.body, expand=TRUE)
  ## Buttons
  g.f1.buttons <- gWidgets2::ggroup(horizontal=TRUE, container=g.f1)
  g.print1     <- gWidgets2::gbutton("Bitmap", handler=h.print1,     container=g.f1.buttons)
  g.print.pdf1 <- gWidgets2::gbutton("PDF",    handler=h.print.pdf1, container=g.f1.buttons)
  ## Device itself
  g.fd1 <- gWidgets2::ggraphics(expand=TRUE, ps=11, container=g.f1)
  d1 <- dev.cur()

  ## Projection
  g.f2 <- gWidgets2::ggroup(horizontal=FALSE, container=g.body, expand=TRUE)
  ## Buttons  
  g.f2.buttons <- gWidgets2::ggroup(horizontal=TRUE, container=g.f2)  
  g.print2     <- gWidgets2::gbutton("Bitmap", handler=h.print2,     container=g.f2.buttons)
  g.print.pdf2 <- gWidgets2::gbutton("PDF",    handler=h.print.pdf2, container=g.f2.buttons)
  ## Device itself
  g.fd2 <- gWidgets2::ggraphics(expand=TRUE, ps=11, container=g.f2)
  d2 <- dev.cur()
  
  ## Status bar
  ## g.statusbar <- ggroup(container=g.rows)
  g.statusbar <- gWidgets2::gframe("", expand=TRUE, container=g.rows)
  g.status <- gWidgets2::glabel("", container=g.statusbar)
  gWidgets2::addSpring(g.statusbar)

  ## Disable buttons initially
  unsaved.data(FALSE)
  enable.widgets(FALSE)

  ## Have to add the hander to the notebook at the end, otherwise
  ## there are complaints about various components not being defined.
  gWidgets2::addHandlerChanged(g.nb, handler=h.show)

  return(invisible(list(getData=getData, env=this)))
}

