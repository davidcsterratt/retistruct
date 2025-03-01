## ---------- Translates strings from UI to function references  ----------

## Converts projections into a function, as values returned from ui are strings.
translateProjections <- function() {
  return(list("0" = azimuthal.equidistant,
              "1" = azimuthal.equalarea,
              "2" = azimuthal.conformal,
              "3" = sinusoidal,
              "4" = orthographic))
}

## Converts transforms into a function, as values returned from ui are strings.
translateTransforms <- function() {
  return(list("0" = identity.transform,
              "1" = invert.sphere,
              "2" = invert.sphere.to.hemisphere))
}

## ---------- Enables and disables UI elements to prevent user error  ----------

## Enables or disables a vector of UI elements, using shinyJS.
enable.group <- function(widgets, save = TRUE) {
  for (w in widgets) {
    if (save) {
      shinyjs::enable(w)
    } else {
      shinyjs::disable(w)
    }
  }
}

## Enables or disables all UI elements
enable.widgets <- function(save, state) {
  enable.group(c("add_tear", "remove_tear", "add_fullcut", "remove_fullcut",
                 "move_point", "reconstruct", "properties",
                 "mark_nasal", "mark_dorsal", "mark_od",
                 "phi0", "show", "strain", "flip_dv", "eye",
                 "bitmap1", "bitmap2", "pdf1", "pdf2", "projection",
                 "center.el", "center.az", "transform", "ax.el", "ax.az",
                 "ids"), save)
  if (save  && !is.null(state$a))
    enable.group(c("mark_od"), length(state$a$getFeatureSet("LandmarkSet")$getIDs() > 0))
  if (!retistruct.check.markup(state$a)) {
    enable.group(c("reconstruct"), FALSE)
  }

  ## Turns the cancel button on or off depending on the server mode, see
  ## state$mode in server.R to understand why.
  if (state$mode!=0) {
    enable.group("cancel", TRUE)
    shinyjs::addClass("cancel", "red")
  } else {
    enable.group("cancel", FALSE)
    shinyjs::removeClass("cancel", "red")
  }
}

## Enables or Disables the save button if there is work to save
unsaved.data <- function(save, state) {
  if (save) {
    state$r <- NULL
  }
  enable.group(c("save"), save)
}

## Function to report to set status
set.status <- function(output, ...) {
  ## FIXME #46: status is not work; report() added to
  ## ensure output
  report(...)
  output$status <- renderText(paste0(...))
}

## ---------- Plotting functions  ----------

plotProjection <- function(max.proj.dim=getOption("max.proj.dim"),
                           markup=NULL, state, input) {
  validate(
    need(is.numeric(input$center.el), "Center el is not a numeric. Check that the input is not an expression, or non-empty."),
    need(is.numeric(input$center.az), "Center az is not a numeric. Check that the input is not an expression, or non-empty."),
    need(is.numeric(input$ax.el), "Axis el is not a numeric. Check that the input is not an expression, or non-empty."),
    need(is.numeric(input$ax.az), "Axis az is not a numeric. Check that the input is not an expression, or non-empty."),
    need(el.low <= input$center.el && input$center.el <= el.high, "Center el must be between 0 and 90"),
    need(az.low <= input$center.az && input$center.az <= az.high, "Center az must be between 0 and 90"),
    need(el.low <= input$ax.el && input$ax.el <= el.high, "Ax el must be between 0 and 90"),
    need(az.low <= input$ax.az && input$ax.az <= az.high, "Ax az must be between 0 and 90")
    )

  if (is.null(markup)) {
    markup <- ("markup" %in% input$show)
  }

  projection(state$r,
             datapoints=("points" %in% input$show),
             datapoint.means=("point_mean" %in% input$show),
             landmarks=("landmarks" %in% input$show),
             transform=translateTransforms()[[input$transform]],
             projection=translateProjections()[[input$projection]],
             axisdir=cbind(phi=input$ax.el, lambda=input$ax.az),
             proj.centre=cbind(phi=input$center.el, lambda=input$center.az),
             datapoint.contours=("point_contour" %in% input$show),
             grouped=("counts" %in% input$show),
             grouped.contours=("point_contour" %in% input$show),
             grid=("grid" %in% input$show),
             markup=markup,
             ids=input$ids,
             max.proj.dim=max.proj.dim)

  if (!is.null(state$r$EOD)) {
    polartext(paste("OD displacement:",
                    format(state$r$EOD, digits=3, nsmall=2), "deg"))
  }
}

## Similar to do.plot, captures into a graphics device instead of a shiny output
## used for saving to pdf or bitmap
##' @importFrom grDevices dev.set
do.plot.silent <- function(markup=NULL, state, input, d1=NULL, d2=NULL) {
    if (is.null(markup)) {
      markup = ("markup" %in% input$show)
    }

    if (is.null(state$r)) {
      r <- state$a
    } else {
      r <- state$r
    }

    if (input$strain) {   # Strain plot
      ## Capture the plot for printing to file
      if (!is.null(d1)) {
        dev.set(d1)
        par(mar=c(0.5, 0.5, 0.5, 0.5), ps=11)
        flatplot(r, axt="n",
                 datapoints=FALSE,
                 landmarks=FALSE,
                 markup=FALSE,
                 stitch=FALSE,
                 grid=FALSE,
                 mesh=FALSE,
                 strain=TRUE,
                 scalebar=1)
        dev.off()
      }
      if (!is.null(d2)) {
        dev.set(d2)
        lvsLplot(r)
        dev.off()
      }
    } else {

      if (!is.null(d1)) {
        dev.set(d1)
        par(mar=c(0.5, 0.5, 0.5, 0.5), ps=11)
        flatplot(r, axt="n",
                 datapoints=("points" %in% input$show),
                 grouped=("counts" %in% input$show),
                 landmarks=("landmarks" %in% input$show),
                 markup=markup,
                 stitch=("stitch" %in% input$show),
                 grid=("grid" %in% input$show),
                 ids=input$ids,
                 mesh=FALSE,
                 scalebar=1)
        dev.off()
      }

      if (!is.null(d2)) {
        dev.set(d2)
        par(mar=c(0.7, 0.7, 0.7, 0.7), ps=11)
        plotProjection(max.proj.dim=400, markup=markup, state=state, input=input)
        dev.off()
      }
    }
  }

## Plot in edit pane
do.plot <- function(markup=NULL, state, input, output) {
  if (is.null(markup)) {
    markup = ("markup" %in% input$show)
  }

  if (is.null(state$r)) {
    r <- state$a
  } else {
    r <- state$r
  }

  if (input$strain) {   # Strain plot
    output$plot1 <- renderPlot({
      ## Capture the plot for printing to file
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
    })

    output$plot2 <- renderPlot({
      lvsLplot(r)
    })

    output$plot3 <- renderRglwidget({
      sphericalplot(r, strain=TRUE, datapoints=FALSE)
      rglwidget()
    })

  } else {
    output$plot1 <- renderPlot({
      par(mar=c(0.5, 0.5, 0.5, 0.5))
      flatplot(r, axt="n",
               datapoints=("points" %in% input$show),
               grouped=("counts" %in% input$show),
               landmarks=("landmarks" %in% input$show),
               markup=markup,
               stitch=("stitch" %in% input$show),
               grid=("grid" %in% input$show),
               ids=input$ids,
               mesh=FALSE,
               scalebar=1)
    })

    output$plot2 <- renderPlot({
      par(mar=c(0.7, 0.7, 0.7, 0.7))
      plotProjection(max.proj.dim=400, markup=markup, state=state, input=input)
    })

    output$plot3 <- renderRglwidget({
      sphericalplot(r, datapoints=("points" %in% input$show))
      rglwidget()
    })
  }
  set.status(output, "")
}

## ---------- Title bar button handlers ----------
h.open <- function(state, input, output, session) {
  req(state$dataset)
  ## Read the raw data
  catchErrorsRecordWarnings({
    state$a <- retistruct.read.dataset(state$dataset, report=FALSE)
  }, warning=function(w) h.warning(w, state, session), error=function(e) h.error(e, session))

  ## Read the markup
  catchErrorsRecordWarnings({
    state$a <- retistruct.read.markup(state$a, error=message)
  }, warning=function(w) h.warning(w, state, session), error=function(e) h.error(e, session))

  # Update UI based on loaded data
  updateNumericInput(session, "phi0", value = state$a$phi0 * 180 / pi)
  updateRadioButtons(session, "eye", selected = state$a$side)
  updateCheckboxInput(session, "flip_dv", value = state$a$DVflip)

  catchErrorsRecordWarnings({
    state$r <- retistruct.read.recdata(state$a, check=TRUE)
  }, warning=function(w) h.warning(w, state, session), error=function(e) h.error(e, session))

  ids <- state$a$getIDs()
  if (length(ids) > 0) {
    updateCheckboxGroupInput(session, "ids", choices=ids, selected=ids)
  }

  unsaved.data(FALSE, state=state)
  enable.widgets(TRUE, state)
  do.plot(state=state, input=input, output=output)
}

## Handler for saving state
h.save <- function(h, state, ...) {
  retistruct.save.markup(state$a)
  ## If the reconstruction doesn't exist, remove the reconstruction
  ## file to ensure consistency
  if (is.null(state$r)) {
    unlink(file.path(state$a$dataset, "r.Rdata"))
  } else {
    retistruct.save.recdata(state$r)
  }
  retistruct.export.matlab(state$r)
  unsaved.data(FALSE, state=state)
}

## Handler for reconstruction
h.reconstruct <- function(h, state, input, output, session, ...) {
  validate(
    need(is.numeric(input$center.el), "Center el is not a numeric. Check that the input is not an expression, or non-empty."),
    need(is.numeric(input$center.az), "Center az is not a numeric. Check that the input is not an expression, or non-empty."),
    need(is.numeric(input$ax.el), "Axis el is not a numeric. Check that the input is not an expression, or non-empty."),
    need(is.numeric(input$ax.az), "Axis az is not a numeric. Check that the input is not an expression, or non-empty."),
    need(el.low <= input$center.el && input$center.el <= el.high, "Center el must be between 0 and 90"),
    need(az.low <= input$center.az && input$center.az <= az.high, "Center az must be between 0 and 90"),
    need(el.low <= input$ax.el && input$ax.el <= el.high, "Ax el must be between 0 and 90"),
    need(az.low <= input$ax.az && input$ax.az <= az.high, "Ax az must be between 0 and 90")
  )

  unsaved.data(TRUE, state)
  enable.widgets(FALSE, state)
  catchErrorsRecordWarnings({
    state$r <- retistruct.reconstruct(state$a, report=function(m) set.status(output, m),
                                      plot.3d=getOption("show.sphere"),
                                      shinyOutput=output)
  }, warning=function(w) h.warning(w, state, session), error=function(e) h.error(e, session))
  enable.widgets(TRUE, state)
  do.plot(state=state, input=input, output=output)
  set.status(output, "")
}

## ----------  Handlers for loading demos ----------

h.demo1 <- function(state, input, output, session, extdata, directory1, directory2) {
    state$dataset <- file.path(extdata, directory1, directory2)
    h.open(state, input, output, session)
}

h.demo2 <- function(state, input, output, session, extdata.demos, directory1, directory2) {
  dataset <- file.path(extdata.demos, directory1, directory2)

  if (!file.exists(dataset)) {
    showNotification(
      "Install the retistructdemos package using
      devtools::install_github(\n\"davidcsterratt/retistruct/pkg/retistructdemos\n\")",
      duration=NULL, closeButton=TRUE, type="error", session=session)
    stop()
  } else {
    state$dataset <- dataset
    h.open(state, input, output, session)
  }
}


## ---------- Handlers for modifying shiny server state  ----------

## Convenience function for setting the server mode
set.state <- function(state, mode) {
  state$mode <- mode
  clear.points(state)
  unsaved.data(TRUE, state)
  enable.widgets(FALSE, state)
}

## Convenience function for resetting server mode
reset.state <- function(state) {
  state$mode <- 0
  clear.points(state)
  unsaved.data(TRUE, state)
  enable.widgets(TRUE, state)
}

## Convenience function for capturing a click to the server state
add.point <- function(state, x, y) {
  state$points_x <- unique(c(state$points_x, x))
  state$points_y <- unique(c(state$points_y, y))
}

## Convencience function for resetting captured clicks
clear.points <- function(state) {
  state$points_x <- c()
  state$points_y <- c()
}

## A simplified version of the built in identify function, that doesn't require
## a graphics device, which shiny does not support.
h.identify <- function(click_x, click_y, points_x, points_y) {
  selected <- c()
  distances <- sqrt((points_x - click_x)^2 + (points_y - click_y)^2)
  closest_index <- which.min(distances)
  selected <- c(selected, closest_index)
  return(selected)
}

## ---------- Handlers for modifying left plot  ----------

## Add tear handler
h.add <- function(state, input, output, session, xs, ys, ...) {
  P <- state$a$getPoints()
  pids <- c()
  for (i in 0:3) {
    pids <- c(pids, h.identify(xs[i], ys[i], P[,"X"], P[,"Y"]))
  }

  catchErrorsRecordWarnings({
      state$a$addTear(pids)
    }, error=function(e) h.error(e, session), warning=function(w) h.warning(w, state, session))
  do.plot(state=state, input=input, output=output)
}

## Handler for moving a point in a tear
h.move <- function(state, input, output, session, xs, ys, ...) {
  P <- state$a$getPoints()

  id1 <- h.identify(xs[1], ys[1], P[,"X"], P[,"Y"])

  ## Locate tear ID in which the point occurs
  tid <- state$a$whichTear(id1)

  ## If there is a tear in which it occurs, select a point to move it to
  if (!is.na(tid)) {
    ## Select second point
    id2 <- h.identify(xs[2], ys[2], P[,"X"], P[,"Y"])

    ## Get point ids of exsiting tear
    pids <- state$a$getTear(tid)

    ## Replace old point with desired new point
    if (length(id2)) pids[pids==id1] <- id2

    ## It is possible to get the apex and vertex mixed up when moving points.
    ## Fix any errors.
    pids <- state$a$labelTearPoints(pids)
    state$a$removeTear(tid)
    state$a$addTear(pids)
    set.status(output, "Moved points successfully")
  } else {
    set.status(output, "Error: no tear at inital selected point")
  }

  ## Display and cleanup
  do.plot(state=state, input=input, output=output)
}

## Remove tear handler
h.remove <- function(state, input, output, x, y, ...) {
  P <- state$a$getPoints()
  id <- c()
  id <- h.identify(x, y, P[,"X"], P[,"Y"])
  state$a$removeTear(state$a$whichTear(id))
  do.plot(state=state, input=input, output=output)
}

## Handler for adding a cut
h.add.fullcut <- function(state, input, output, session, xs, ys, ...) {
  P <- state$a$getPoints()
  pids <- c()
  for (i in 0:4) {
    pids <- c(pids, h.identify(xs[i], ys[i], P[,"X"], P[,"Y"]))
  }

  catchErrorsRecordWarnings({
      state$a$addFullCut(pids)
    }, error=function(e) h.error(e, session), warning=function(w) h.warning(w, state, session))
  do.plot(state=state, input=input, output=output)   #DOTHIS
}

## Remove tear handler
h.remove.fullcut <- function(state, input, output, x, y, ...) {
  P <- state$a$getPoints()
  id <- c()
  id <- h.identify(x, y, P[,"X"], P[,"Y"])
  state$a$removeFullCut(state$a$whichFullCut(id))
  do.plot(state=state, input=input, output=output)
}

## Mark nasal handler
h.mark.n <- function(state, input, output, session, x, y, ...) {
  P <- state$a$getPoints()
  id <- h.identify(x, y, P[,"X"], P[,"Y"])
  catchErrorsRecordWarnings({
    state$a$setFixedPoint(id, "Nasal")
  }, warning=function(w) h.warning(w, state, session), error=function(e) h.error(e, session))
  do.plot(state=state, input=input, output=output)
}

## Mark dorsal handler
h.mark.d <- function(state, input, output, session, x, y, ...) {
  P <- state$a$getPoints()
  id <- h.identify(x, y, P[,"X"], P[,"Y"])
  catchErrorsRecordWarnings({
    state$a$setFixedPoint(id, "Dorsal")
  }, warning=function(w) h.warning(w, state, session), error=function(e) h.error(e, session))
  do.plot(state=state, input=input, output=output)
}

## Handler for marking optic disc
h.mark.od <- function(state, input, output, session, x, y, ...) {

  ## Convert list of segments to a matrix of points
  Sm <- NULL
  fs <- state$a$getFeatureSet("LandmarkSet")
  Ss <- fs$getFeatures()
  for (S in Ss) {
    Sm <- rbind(Sm, S)
  }

  ## Identify a point
  id <- h.identify(x, y, Sm[,1], Sm[,2])

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
  checked <- state$a$getIDs() %in% c(input$ids, "OD")
  ids <- state$a$getIDs()
  if (length(ids) > 0) {
    updateCheckboxGroupInput(session, "ids", choices=ids, selected=ids)
  }

  do.plot(state=state, input=input, output=output)
}

##'@importFrom grDevices png jpeg tiff dev.cur
save.bitmap <- function(state, input, output, session, file, left) {
  x = getOption("max.proj.dim")

  if (grepl("\\.png$", file, ignore.case=TRUE)) {
    png(file, width=x, height=x)
  } else if (grepl("\\.jpeg$|\\.jpg$", file, ignore.case=TRUE)) {
    jpeg(file, width=x, height=x)
  } else if (grepl("\\.tif$|\\.tiff$", file, ignore.case=TRUE)) {
    tiff(file, width=x, height=x)
  } else {
    file <- paste0(file, ".png")
    png(file, width=x, height=x)
  }

  did <- dev.cur()

  # Call the plotting function
  if (left) {
    do.plot.silent(state=state, input=input, d1=did)
  } else {
    do.plot.silent(state=state, input=input, d2=did)
  }
}

##' @importFrom grDevices pdf dev.cur
save.pdf <- function(state, input, output, session, file, left) {
  width = getOption("retistruct.print.pdf.width")
  pdf(file, width=width, height=width)
  did <- dev.cur()
  # Call the plotting function
  if (left) {
    do.plot.silent(state=state, input=input, d1=did)
  } else {
    do.plot.silent(state=state, input=input, d2=did)
  }
}

## ---------- Warning and error handlers  ----------
# Error Message
h.error <- function(e, session) {
  showNotification(conditionMessage(e), duration=10, type="error", session=session)
}

## Warning message
h.warning <- function(w, state, session) {
  e <- conditionMessage(w)
  if (!(any(e %in% state$prior.warnings))) {
    showNotification(e, duration=10, type="warning", session=session)
    state$prior.warnings <- c(state$prior.warnings, e)
  }
}

## Allows shiny app to continue running, whilst preventing warnings from
## stopping functions completing
catchErrorsRecordWarnings <- function(expr, warning, error) {
  withCallingHandlers({
    tryCatch({
      expr
    }, error=error)
  }, warning=warning)
}
