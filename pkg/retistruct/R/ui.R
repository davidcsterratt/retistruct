library(shiny)
library(shinyjs)
library(shinyFiles)
library(bslib)

## Get and name available projections
getProjections <- function() {
  return(list("Azimuthal Equidistant" = '0',
              "Azimuthal Equal Area"  = '1',
              "Azimuthal Conformal"   = '2',
              "Sinusoidal"            = '3',
              "Orthographic"          = '4'))
}

## Get and name available transforms
getTransforms <- function() {
  return(list("None"   = "0",
              "Invert" = "1",
              "Invert to hemisphere" = "2"))
}

## Construct the version string
version.string <- function() {
  return(paste0("Retistruct ",
                utils::packageDescription("retistruct", fields="Version"),
                " (", utils::packageDescription("retistruct", fields="Date"), ")"))
}

##' @description
##' The Shiny UI element, runs on a browser and is similar to HTML, attempt to 
##' mimic the original Retistruct UI as closely as possible.
ui <- fluidPage(
  useShinyjs(),
  ## CSS class to turn cancel button red when active.
  tags$style(HTML("
      .red {
         background-color: red;
        color: white;
        font-weight: bold;
        padding: 5px;
      }
    ")),
  
  titlePanel(
    windowTitle=version.string(),
    fluidRow(
      "Retistruct",
      ## Top level Nav bar
      shinyDirButton("open", "Open", "Select a directory...", icon=icon("folder-open")),
      actionButton("save", "Save", icon=icon("floppy-disk")),
      actionButton("reconstruct", "Reconstruct", icon=icon("hammer")),
      actionButton("properties", "Properties", icon=icon("gear")),
      actionButton("demo", "Demo"),
      actionButton("about", "About"),
      actionButton("cancel", "Cancel")
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel(
          "Edit",
          # Edit panel options
          actionButton("add_tear", "Add Tear", style="width:100%"),
          br(),
          actionButton("move_point", "Move Points", style="width:100%"),
          br(),
          actionButton("remove_tear", "Remove Tear", style="width:100%"),
          br(),
          actionButton("mark_nasal", "Mark Nasal", style="width:100%"),
          br(),
          actionButton("mark_dorsal", "Mark Dorsal", style="width:100%"),
          br(),
          actionButton("mark_od", "Mark OD", style="width:100%"),
          br(), 
          checkboxInput("flip_dv", "Flip DV"),
          radioButtons(
            inputId = "eye",
            label = "Eye",
            choices = list(
              "Left" = "Left",
              "Right" = "Right"
            )
          ),
          numericInput("phi0", "Phi0", 0, 0, 100),
          checkboxInput("strain", "Strain"),
        ),
        # View panel
        tabPanel(
          "View",
          # UI check boxes for all the show options
          checkboxGroupInput(
            "show",
            "Show",
            c(
              "Markup" = "markup",
              "Stitch" = "stitch",
              "Grid" = "grid",
              "Landmarks" = "landmarks",
              "Points" = "points",
              "Point Means" = "point_mean",
              "Counts" = "counts",
              "Point Contours" = "point_contour"
            ),
            selected = c("markup", "landmarks", "points", "counts")
          ),
          selectInput("projection", "Projection", choices = getProjections()),
          tags$strong(textOutput("projcentre")),
          fluidRow(
            column(6, numericInput("center.el", "El", value = 0, min = 0, max = 90)),
            column(6, numericInput("center.az", "Az", value = 0, min = 0, max = 90))
          ),
          selectInput("transform", "Transform", choices = getTransforms()),
          tags$strong(textOutput("axdir")),
          fluidRow(
            # Width is 6
            column(6, numericInput("ax.el", "El", value = 90, min = 0, max = 90)),
            column(6, numericInput("ax.az", "Az", value = 0, min = 0, max = 90))
          ),
          checkboxGroupInput("ids", "IDs", choices=c("All"), selected =c("All")),
        )
      ),
      width = 3 # Sidebar width
    ),
    
    mainPanel(
      navset_tab(
        nav_panel(
          "Main Panel",
          fluidRow(
            column(
              width=6,
              actionButton("bitmap1", "BitMap"),
              actionButton("pdf1", "PDF"),
              plotOutput("plot1", height="60vh", click="plot1click")
            ),
            column(
              width=6,
              actionButton("bitmap2", "BitMap"),
              actionButton("pdf2", "PDF"),
              plotOutput("plot2", height="60vh")
            ),
          )
        ),
        nav_panel(
          "3D Panel",
          rglwidgetOutput("plot3"),
        )
      )
    )
  ),
  textOutput("status")
)
