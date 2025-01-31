## --------- All the server handlers are located in R/server-handler.R ---------
time_out <- 1000 # How long to leave a status text before clearing in ms
abort.text <- "Press the red cancel button at the top of the screen to cancel"

##' File system directories used by shinyFiles
##' @importFrom fs path_home
directories <- c(Home=fs::path_home())
extdata       <- file.path(system.file(package = "retistruct"), "extdata")
extdata.demos <- file.path(system.file(package = "retistructdemos"), "extdata")

##' @title Retistruct Shiny Server
##' @description The R shiny server responsible for storing a state for each 
##' session, handling inputs from the UI to the server, and plotting outputs
##' to the UI. The arguments are all handled by the shiny package and this 
##' function should not be instantiated manually.
##' @param input object that holds the UI state (Managed automatically by shiny)
##' @param output sends new outputs to the UI (Managed automatically by shiny)
##' @param session controls each open instance (Managed automatically by shiny)
##' @importFrom shiny observeEvent renderPlot renderText
##' @importFrom shinyjs useShinyjs enable disable delay
##' @importFrom shinyFiles shinyDirChoose parseDirPath
server <- function(input, output, session) {
  useShinyjs()
  
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
  
  ## Open project handler
  observeEvent(input$open, {
    shinyDirChoose(input, "open", roots = directories, session = session)
    dirname <- parseDirPath(roots = directories, input$open)
    if (length(dirname) > 0) {
      state$dataset <- dirname
      h.open(state=state, input=input, output=output, session=session)
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
  
  ## Demo handler
  observe({
    showModal(
      modalDialog(
        title = "Choose one of the following demos to load",
        easy_close = TRUE,
        size="s",
        fluidPage(
          fluidRow(actionButton("fig1", "Figure 1")),
          fluidRow(actionButton("fig2a", "Figure 2A-D: Low deformation")),
          fluidRow(actionButton("fig2e", "Figure 2E-H: High deformation")),
          fluidRow(actionButton("smi32", "smi32")),
          fluidRow(actionButton("fig6lc", "Figure 6 Left Contra")),
          fluidRow(actionButton("fig6li", "Figure 6 Left Ipsi")),
          fluidRow(actionButton("fig6rc", "Figure 6 Right Contra")),
          fluidRow(actionButton("fig6ri", "Figure 6 Right Ipsi"))
        )
      )
    )
  }) |> bindEvent(input$demo)
  
  ## Demo button handlers
  observeEvent(input$fig1, {
    h.demo1(state, input, output, session, "GM509", "R-CONTRA")
    removeModal()
  })
  
  observeEvent(input$fig2a, {
    h.demo1(state, input, output, session, "GMB530", "R-CONTRA")
    removeModal()
  })
  
  observeEvent(input$fig2e, {
    h.demo1(state, input, output, session, "GM182-4", "R-CONTRA")
    removeModal()
  })
  
  observeEvent(input$smi32, {
    h.demo1(state, input, output, session, "smi32", "")
    removeModal()
  })
  
  observeEvent(input$fig6lc, {
    h.demo2(state, input, output, session, "Figure_6-data", "left-contra")
    removeModal()
  })
  
  observeEvent(input$fig6li, {
    h.demo2(state, input, output, session, "Figure_6-data", "left-ipsi")
    removeModal()
  })
  
  observeEvent(input$fig6rc, {
    h.demo2(state, input, output, session, "Figure_6-data", "right-contra")
    removeModal()
  })
  
  observeEvent(input$fig6ri, {
    h.demo2(state, input, output, session, "Figure_6-data", "right-ipsi")
    removeModal()
  })
  
  ## About handler
  observe({
    showModal(
      modalDialog(
        title = "About",
        easy_close = TRUE,
        "Retistruct was written by David Sterratt at the University of Edinburgh
        , and tested by Daniel Lyngholm and Ian Thompson at the MRC Centre for
        Developmental Neurobiology, KCL. This work was supported by a Programme
        Grant from the Wellcome Trust (G083305)."
      )
    )
  }) |> bindEvent(input$about)
  
  ## Cancel button, this is to let the server know that the mode should be reset 
  ## and the captured points should be cleared.
  observeEvent(input$cancel, {
    reset.state(state)
    set.status(output, "Operation Cancelled.")
    delay(time_out, {set.status(output, "")})
  })
  
  observeEvent(input$smi32, {
    removeModal()
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
    
    ## Move points handler
    if (state$mode==2) {
      add.point(state, input$plot1click$x, input$plot1click$y)
      print(length(state$points_x))
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
        print(c(state$points_x, state$points_y))
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
    set.status(output, "Click on the three points of the tear in any order. Double click to remove the latest point added and press cancel, to cancel.")
  })
  
  observeEvent(input$move_point, {
    set.state(state, 2)
    set.status(output, "Click on apex or vertex to move.")
    print(length(state$points_x))
  })
  
  observeEvent(input$remove_tear, {
    set.state(state, 3)
    set.status(output, "Click on any point of the tear to remove it and press cancel, to cancel.")
  })
  
  observeEvent(input$mark_nasal, {
    set.state(state, 4)
    set.status(output, "Click on any point of the tear to remove it and press cancel, to cancel.")
  })
  
  observeEvent(input$mark_od, {
    set.state(state, 5)
    set.status(output, "Click on any point of the tear to remove it and press cancel, to cancel.")
  })
  
  
  # flip dv handler
  observeEvent(input$flip_dv, {
      unsaved.data(TRUE, state)
      state$a$DVflip <- input$flip_dv
      do.plot(state=state, input=input, output=output)
  }, ignoreInit = TRUE)
  
  # phi0 handler
  observeEvent(input$phi0, {
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
  
  ## Needed to render text boxes
  output$projcentre <- renderText("Projection Centre")
  output$axdir <- renderText("Axis Direction")
  
  ## Initially disable widgets and disable save button.
  enable.widgets(FALSE, state)
  unsaved.data(FALSE, state)
}