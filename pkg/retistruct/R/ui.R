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

## Default values for properties
default.val <- c(
  outline = "black",
  tear = "black",
  stitch = "cyan",
  maj = "black",
  min = "grey",
  proj = 400
)

properties.ui <- modalDialog(
    title="Properties",
    easy_close=TRUE,
    tags$strong("Colours"),
    selectInput("out_colour", "Outline Colour", choices=cols,
                selected=getOption("outline.col", default.val["outline"])),
    selectInput("tear_colour", "Tear Colour", choices=cols,
                selected=getOption("TF.col", default=default.val["tear"])),
    selectInput("stitch_colour", "Stitch Colour", choices=cols,
                selected=getOption("stitch.col", default=default.val["stitch"])),
    selectInput("major_colour", "Major Gridline Colour", choices=cols,
                selected=getOption("grid.maj.col", default=default.val["maj"])),
    selectInput("minor_colour", "Minor Gridline Colour", choices=cols,
                selected=getOption("grid.min.col", default.val["min"])),
    tags$strong("Bitmap/PDF Output"),
    numericInput("output_width", "Maximum Width of Projection",
                 value=getOption("max.proj.dim", default=default.val["proj"])),
  ) 

demo.ui <- modalDialog(
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

about.ui <- modalDialog(
  title = "About",
  easy_close = TRUE,
  "Retistruct was written by David Sterratt at the University of Edinburgh
        , and tested by Daniel Lyngholm and Ian Thompson at the MRC Centre for
        Developmental Neurobiology, KCL. This work was supported by a Programme
        Grant from the Wellcome Trust (G083305)."
  )

##' @title Retistruct UI
##' @description
##' The Shiny UI element, runs on a browser and is similar to HTML, attempt to 
##' mimic the original Retistruct UI as closely as possible.
##' @importFrom shinyjs useShinyjs
##' @import shiny
##' @importFrom shinyFiles shinyDirButton
##' @importFrom bslib navset_tab nav_panel
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
