## --------- All the server handlers are located in R/server-handler.R ---------
time_out <- 500 # How long to leave a status text before clearing in ms
abort.text <- "Press cancel to abort adding the tear."
cols <- c("black", "red", "green3", "blue", "cyan", "magenta", "yellow", "gray")
el.low <- -90
el.high <- 90
az.low <- -180
az.high <- 180

##' File system directories used by shinyFiles
##' @importFrom fs path_home
directories <- function() {c(Home=fs::path_home())}


##' @title Retistruct Shiny Server
##' @author Jan Okul
##' @description The R shiny server responsible for storing a state for each
##' session, handling inputs from the UI to the server, and plotting outputs
##' to the UI. The arguments are all handled by the shiny package and this
##' function should not be instantiated manually.
##' @param input object that holds the UI state (Managed automatically by shiny)
##' @param output sends new outputs to the UI (Managed automatically by shiny)
##' @param session controls each open instance (Managed automatically by shiny)
##' @importFrom shiny observeEvent renderPlot renderText
##' @importFrom shinyjs useShinyjs enable disable delay
##' @import shinyFiles
server <- function(input, output, session) {
  ## The server state which remains seperate for each instance of shiny running
  state <- new.env()
  state$dataset <- NULL
  state$a <- NULL ## Annotation object
  state$r <- NULL ## Reconstruction object
  state$prior.warnings <- c() # #For the warning handler

  ## The "running mode" of the server, 0 is default, 1-6 is set to let the
  ## server know that incoming clicks should be captured by the click handler,
  ## and to recognize which handler to call when points have been captured.
  state$mode <- 0
  ## The x and y coordinates captured by the click handler
  state$points_x <- c()
  state$points_y <- c()
  extdata       <- file.path(system.file(package = "retistruct"), "extdata")
  extdata.demos <- file.path(system.file(package = "retistructdemos"), "extdata")

  ## shinyFiles handlers
  shinyDirChoose(input, "open", session=session, roots=directories())
  shinyFileSave(input, "bitmap1", session=session, roots=directories())
  shinyFileSave(input, "pdf1", session=session, roots=directories())
  shinyFileSave(input, "bitmap2", session=session, roots=directories())
  shinyFileSave(input, "pdf2", session=session, roots=directories())

  ## ------------------- Navbar handlers -------------------
  ## Open project handler
  observeEvent(input$open, {
    dirname <- parseDirPath(roots=directories(), input$open)
    if (length(dirname) > 0) {
      state$dataset <- dirname
      tryCatch({
        h.open(state=state, input=input, output=output, session=session)
      }, error=function(e) return())
    }
  })

  ## Save project handler
  observeEvent(input$save, {
    h.save(state=state)
  })

  ## Reconstruct handler
  observeEvent(input$reconstruct, {
    h.reconstruct(state=state, input=input, output=output, session=session)
  })

  ## Properties handler
  observeEvent(input$properties, showModal(properties.ui()))

  ## ---------- Properties option handlers ----------
  observeEvent(input$out_colour, {
    options("outline.col" = input$out_colour)
    do.plot(state=state, input=input, output=output)
  }, ignoreInit = TRUE, ignoreNULL = TRUE)

  observeEvent(input$tear_colour, {
    options("TF.col" = input$tear_colour)
    options("TB.col" = input$tear_colour)
    options("V.col" = input$tear_colour)
    do.plot(state=state, input=input, output=output)
  }, ignoreInit = TRUE, ignoreNULL = TRUE)

  observeEvent(input$stitch_colour, {
    options("stitch.col" = input$stitch_colour)
    do.plot(state=state, input=input, output=output)
  }, ignoreInit = TRUE, ignoreNULL = TRUE)

  observeEvent(input$major_colour, {
    options("grid.maj.col" = input$major_colour)
    do.plot(state=state, input=input, output=output)
  }, ignoreInit = TRUE, ignoreNULL = TRUE)

  observeEvent(input$minor_colour, {
    options("grid.min.col" = input$minor_colour)
    do.plot(state=state, input=input, output=output)
  }, ignoreInit = TRUE, ignoreNULL = TRUE)

  observeEvent(input$output_width, {
    options("max.proj.dim" = input$output_width)
    do.plot(state=state, input=input, output=output)
  }, ignoreInit = TRUE, ignoreNULL = TRUE)

  observeEvent(input$pdf_width, {
    options("retistruct.print.pdf.width" = input$pdf_width)
    do.plot(state=state, input=input, output=output)
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
  ## Demo button handler
  observeEvent(input$demo, showModal(demo.ui))

  ## ---------- Demo button option handlers ----------
  observeEvent(input$fig1, {
    enable.widgets(FALSE, state)
    removeModal() ## Closes the overlay automatically
    h.demo1(state, input, output, session, extdata, "GM509", "R-CONTRA")
    enable.widgets(TRUE, state)
  })

  observeEvent(input$fig2a, {
    enable.widgets(FALSE, state)
    removeModal()
    h.demo1(state, input, output, session, extdata, "GMB530", "R-CONTRA")
    enable.widgets(TRUE, state)
  })

  observeEvent(input$fig2e, {
    enable.widgets(FALSE, state)
    removeModal()
    h.demo1(state, input, output, session, extdata, "GM182-4", "R-CONTRA")
    enable.widgets(TRUE, state)
  })

  observeEvent(input$smi32, {
    enable.widgets(FALSE, state)
    removeModal()
    h.demo1(state, input, output, session, extdata, "smi32", "")
    enable.widgets(TRUE, state)
  })

  observeEvent(input$octants, {
    enable.widgets(FALSE, state)
    removeModal()
    h.demo1(state, input, output, session, extdata, "ijroimulti", "")
    enable.widgets(TRUE, state)
  })

  observeEvent(input$fig6lc, {
    enable.widgets(FALSE, state)
    removeModal()
    ## Stops enabling widgets if no directory
    tryCatch({
      h.demo2(state, input, output, session, extdata.demos, "Figure_6-data", "left-contra")
    }, error=function(e){return()})
    enable.widgets(TRUE, state)
  })

  observeEvent(input$fig6li, {
    enable.widgets(FALSE, state)
    removeModal()
    tryCatch({
      h.demo2(state, input, output, session, extdata.demos, "Figure_6-data", "left-ipsi")
    }, error=function(e){return()})
    enable.widgets(TRUE, state)
  })

  observeEvent(input$fig6rc, {
    enable.widgets(FALSE, state)
    removeModal()
    tryCatch({
      h.demo2(state, input, output, session, extdata.demos, "Figure_6-data", "right-contra")
    }, error=function(e){return()})
    enable.widgets(TRUE, state)
  })

  observeEvent(input$fig6ri, {
    enable.widgets(FALSE, state)
    removeModal()
    tryCatch({
      h.demo2(state, input, output, session, extdata.demos, "Figure_6-data", "right-ipsi")
    }, error=function(e){return()})
    enable.widgets(TRUE, state)
  })

  ## About button handler
  observeEvent(input$about, showModal(about.ui))


  ## Cancel button, this is to let the server know that the mode should be reset
  ## and the captured points should be cleared.
  observeEvent(input$cancel, {
    reset.state(state)
    set.status(output, "Operation Cancelled.")
    delay(time_out, {set.status(output, "")})
  })

  ## Click handler for any click made onto the plot, doesn't do anything if
  ## server mode is 0.
  observeEvent(input$plot1click, {
    ## Add tear handler
    if (state$mode==1) {
      add.point(state, input$plot1click$x, input$plot1click$y)
      ## Call the handler once 3 points are captured.
      if (length(state$points_x)==3) {
        h.add(state, input, output, session, state$points_x, state$points_y)
        reset.state(state)
        set.status(output, "Added tear successfully")
        delay(time_out, {set.status(output, "")})
      }
    }

    ## Remove tear handler
    if (state$mode==3) {
      add.point(state, input$plot1click$x, input$plot1click$y)
      ## Whilst if statement not necessary, it adds redundancy in case there is
      ## more than one point added.
      if (length(state$points_x)==1) {
        h.remove(state, input, output, state$points_x, state$points_y)
        reset.state(state)
        set.status(output, "Removed tear successfully")
        delay(time_out, {set.status(output, "")})
      }
    }

    ## Add full cut handler
    if (state$mode==7) {
      add.point(state, input$plot1click$x, input$plot1click$y)
      ## Whilst if statement not necessary, it adds redundancy in case there is
      ## more than one point added.
      if (length(state$points_x)==4) {
        h.add.fullcut(state, input, output, session, state$points_x, state$points_y)
        reset.state(state)
        set.status(output, "Added full cut successfully")
        delay(time_out, {set.status(output, "")})
      }
    }

    ## Remove full cut handler
    if (state$mode==8) {
      add.point(state, input$plot1click$x, input$plot1click$y)
      ## Whilst if statement not necessary, it adds redundancy in case there is
      ## more than one point added.
      if (length(state$points_x)==1) {
        h.remove.fullcut(state, input, output, state$points_x, state$points_y)
        reset.state(state)
        set.status(output, "Removed full cut successfully")
        delay(time_out, {set.status(output, "")})
      }
    }

    ## Move points handler
    if (state$mode==2) {
      add.point(state, input$plot1click$x, input$plot1click$y)
      if (length(state$points_x) == 1) {
        set.status(output, paste("Click on point to move it to.", abort.text))
      }
      ## Call the handler once 2 points are captured.
      if (length(state$points_x)==2) {
        h.move(state, input, output, session, state$points_x, state$points_y)
        reset.state(state)
        set.status(output, "Moved points successfully")
        delay(time_out, {set.status(output, "")})
      }
    }

    ## Mark nasal handler
    if (state$mode==4) {
      add.point(state, input$plot1click$x, input$plot1click$y)
      ## Whilst if statement not necessary, it adds redundancy in case there is
      ## more than one point added.
      if (length(state$points_x)==1) {
        h.mark.n(state, input, output, session, state$points_x, state$points_y)
        reset.state(state)
        set.status(output, "Marked nasal successfully")
        delay(time_out, {set.status(output, "")})
      }
    }

    ## Mark dorsal handler
    if (state$mode==5) {
      add.point(state, input$plot1click$x, input$plot1click$y)
      ## Whilst if statement not necessary, it adds redundancy in case there is
      ## more than one point added.
      if (length(state$points_x)==1) {
        h.mark.d(state, input, output, session, state$points_x, state$points_y)
        reset.state(state)
        set.status(output, "Marked dorasal successfully")
        delay(time_out, {set.status(output, "")})
      }
    }

    ## Mark OD handler
    if (state$mode==6) {
      add.point(state, input$plot1click$x, input$plot1click$y)
      ## Whilst if statement not necessary, it adds redundancy in case there is
      ## more than one point added.
      if (length(state$points_x)==1) {
        h.mark.od(state, input, output, session, state$points_x, state$points_y)
        reset.state(state)
        set.status(output, "Marked OD successfully")
        delay(time_out, {set.status(output, "")})
      }
    }
  })

  ## These button handlers only change the state of the server, so that the
  ## server is willing to store click inputs, without blocking the rest of the
  ## server, the actual click handling is in the plot click handler.
  observeEvent(input$add_tear, {
    set.state(state, 1)
    set.status(output, paste("Click on the three points of the tear in any order.", abort.text))
  })

  observeEvent(input$remove_tear, {
    set.state(state, 3)
    set.status(output, paste("Click on the apex of the tear to remove.", abort.text))
  })

  observeEvent(input$add_fullcut, {
    set.state(state, 7)
    set.status(output, paste("Click on the four points of the full cut in any order.", abort.text))
  })

  observeEvent(input$remove_fullcut, {
    set.state(state, 8)
    set.status(output, paste("Click on one of the points of the full cut to remove.", abort.text))
  })

  observeEvent(input$move_point, {
    set.state(state, 2)
    set.status(output, paste("Click on apex or vertex to move.", abort.text))
  })

  observeEvent(input$mark_nasal, {
    set.state(state, 4)
    set.status(output, paste("Click on nasal point.", abort.text))
  })

  observeEvent(input$mark_dorsal, {
    set.state(state, 5)
    set.status(output, paste("Click on dorsal point.", abort.text))
  })

  observeEvent(input$mark_od, {
    set.state(state, 6)
    set.status(output, paste("Click on a point on the optic disc.", abort.text))
  })

  # flip dv handler
  observeEvent(input$flip_dv, {
      unsaved.data(TRUE, state)
      state$a$DVflip <- input$flip_dv
      do.plot(state=state, input=input, output=output)
  }, ignoreInit = TRUE)

  # phi0 handler
  observeEvent(input$phi0, {
    if (!is.numeric(input$phi0)) {
      showNotification("Phi0 cannot be an expression, or non-empty. Phi0 will not be updated.", type="error")
      return()
    }
    if (0 > input$phi0 || input$phi0 > 100) {
      showNotification("Phi0 must be between 0 and 100. Phi0 will not be updated.", type="error")
      return()
    }

    unsaved.data(TRUE, state)
    v <- input$phi0
    if (v < -80) {
      v <- -89
    }
    if (v > 89) {
      v <- 89
    }
    state$a$phi0 <- v*pi/180
  }, ignoreInit = TRUE)

  # Eye handler
  observeEvent(input$eye, {
    unsaved.data(TRUE, state)
    state$a$side <- input$eye
    do.plot(state=state, input=input, output=output)
  }, ignoreInit = TRUE)

  # Strain handler
  observeEvent(input$strain, {
    do.plot(markup=TRUE, state=state, input=input, output=output)
  })

  # ------------------- Print to bitmap/PDF handlers -------------------
  observeEvent(input$bitmap1, {
    file <- parseSavePath(roots=directories(), input$bitmap1)
    if (length(file$datapath) > 0) {
      save.bitmap(state, input, output, session, file$datapath, TRUE)
    }
  }, ignoreInit=TRUE, ignoreNULL=TRUE)

  observeEvent(input$pdf1, {
      file <- parseSavePath(roots=directories(), input$pdf1)
      if (length(file$datapath) > 0) {
        save.pdf(state, input, output, session, file$datapath,
                  TRUE)
      }
  }, ignoreInit=TRUE, ignoreNULL=TRUE)

  observeEvent(input$bitmap2, {
    file <- parseSavePath(roots=directories(), input$bitmap2)
    if (length(file$datapath) > 0) {
      save.bitmap(state, input, output, session, file$datapath, FALSE)
    }
  }, ignoreInit=TRUE, ignoreNULL=TRUE)

  observeEvent(input$pdf2, {
    file <- parseSavePath(roots=directories(), input$pdf2)
    if (length(file$datapath) > 0) {
      save.pdf(state, input, output, session, file$datapath,
                FALSE)
    }
  }, ignoreInit=TRUE, ignoreNULL=TRUE)


  ## Needed to render text boxes
  output$projcentre <- renderText("Projection Centre")
  output$axdir <- renderText("Axis Direction")

  ## Initially disable widgets and disable save button.
  enable.widgets(FALSE, state)
  unsaved.data(FALSE, state)
}
