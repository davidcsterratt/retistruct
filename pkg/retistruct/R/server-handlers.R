## Converts projections into a function, as values returned from ui are strings.
translateProjections <- function() {
  return(list('0' = azimuthal.equidistant,
              '1' = azimuthal.equalarea,
              '2' = azimuthal.conformal,
              '3' = sinusoidal,
              '4' = orthographic))
}

## Converts transforms into a function, as values returned from ui are strings.
translateTransforms <- function() {
  return(list("0" = identity.transform,
              "1" = invert.sphere,
              "2" = invert.sphere.to.hemisphere))
}


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
  enable.group(c("add_tear", "move_point", "remove_tear", "reconstruct", "properties",
                 "mark_nasal", "mark_dorsal", "mark_od",
                 "phi0", "show", "strain", "flip_dv", "eye",
                 "bitmap1", "bitmap2", "pdf1", "pdf2", "projection",
                 "center.el", "center.az", "transform", "ax.el", "ax.az",
                 "ids"), save)
  if (save) 
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

plotProjection <- function(max.proj.dim=getOption("max.proj.dim"),
                           markup=NULL, state, input) {
  
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


h.open <- function(state, input, output, session) {
  req(state$dataset)
  ## Read the raw data
  print(state$dataset)
  withCallingHandlers({
    state$a <- retistruct.read.dataset(state$dataset, report=FALSE)
  }, warning=function(w) h.warning(w, state, session), error=function(e) h.error(e, session))
  
  ## Read the markup
  withCallingHandlers({
    state$a <- retistruct.read.markup(state$a, error=message)
  }, warning=function(w) h.warning(w, state, session), error=function(e) h.error(e, session))
  
  # Update UI based on loaded data
  updateNumericInput(session, "phi0", value = state$a$phi0 * 180 / pi)
  updateRadioButtons(session, "eye", selected = state$a$side)
  updateCheckboxInput(session, "flip_dv", value = state$a$DVflip)
  
  withCallingHandlers({
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
  unsaved.data(TRUE, state)
  enable.widgets(FALSE, state)
  withCallingHandlers({
    state$r <- retistruct.reconstruct(state$a, report=function(m) set.status(output, m),
                                      plot.3d=getOption("show.sphere"), 
                                      shinyOutput=output)
  }, warning=function(w) h.warning(w, state, session), error=function(e) h.error(e, session))  
  enable.widgets(TRUE, state)
  do.plot(state=state, input=input, output=output)
  set.status(output, "")
}

# Error Message
h.error <- function(e, session) {
  showNotification(conditionMessage(e), duration=NULL, closeButton=TRUE, type="error", session=session)
  stop(e)
}

## Warning message
h.warning <- function(w, state, session) {
  e <- conditionMessage(w)
  if (!(any(e %in% state$prior.warnings))) {
    showNotification(e, duration=NULL, closeButton=TRUE, type="warning", session=session)
    state$prior.warnings <- c(state$prior.warnings, e)
  }
}

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
  state$points_x <- c(state$points_x, x)
  state$points_y <- c(state$points_y, y)
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

## Add tear handler
h.add <- function(state, input, output, session, xs, ys, ...) {
  P <- state$a$getPoints()
  pids <- c()
  for (i in 0:3) {
    pids <- c(pids, h.identify(xs[i], ys[i], P[,"X"], P[,"Y"]))
  }
  withCallingHandlers({
    state$a$addTear(pids)
  }, warning=function(w) h.warning(w, state, session), error=function(e) h.error(e, session))  
  do.plot(state=state, input=input, output=output)   #DOTHIS
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
  do.plot(state=state, input=input, output=output)   #DOTHIS
}

## Mark nasal handler
h.mark.n <- function(state, input, output, session, x, y, ...) {
  P <- state$a$getPoints()
  id <- h.identify(x, y, P[,"X"], P[,"Y"])
  withCallingHandlers({
    state$a$setFixedPoint(id, "Nasal")
  }, warning=function(w) h.warning(w, state, session), error=function(e) h.error(e, session))  
  do.plot(state=state, input=input, output=output)
}

## Mark dorsal handler
h.mark.d <- function(state, input, output, session, x, y, ...) {
  P <- state$a$getPoints()
  id <- h.identify(x, y, P[,"X"], P[,"Y"])
  withCallingHandlers({
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




