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
picture.types <- list(picture=c("png", "jpg", "jpeg", "tif", "tiff"))
pdf.type <- list(files=c("pdf"))

## Function re-renders the modal every time it's called
properties.ui <- function() {
  modalDialog(
    title="Properties",
    easy_close=TRUE,
    tags$strong("Colours"),
    selectInput("out_colour", "Outline Colour", choices=cols,
                selected=getOption("outline.col")),
    selectInput("tear_colour", "Tear Colour", choices=cols,
                selected=getOption("TF.col")),
    selectInput("stitch_colour", "Stitch Colour", choices=cols,
                selected=getOption("stitch.col")),
    selectInput("major_colour", "Major Gridline Colour", choices=cols,
                selected=getOption("grid.maj.col")),
    selectInput("minor_colour", "Minor Gridline Colour", choices=cols,
                selected=getOption("grid.min.col")),
    tags$strong("Bitmap/PDF Output"),
    numericInput("output_width", "Maximum Width of Projection",
                 value=getOption("max.proj.dim")),
    numericInput("pdf_width", "PDF Width",
                 value=getOption("retistruct.print.pdf.width"))
  )
}

demo.ui <- modalDialog(
    title = "Choose one of the following demos to load",
    easy_close = TRUE,
    size="s",
    fluidPage(
      fluidRow(actionButton("fig1", "Figure 1")),
      fluidRow(actionButton("fig2a", "Figure 2A-D: Low deformation")),
      fluidRow(actionButton("fig2e", "Figure 2E-H: High deformation")),
      fluidRow(actionButton("smi32", "SMI32")),
      fluidRow(actionButton("octants", "Octants")),
      fluidRow(actionButton("fig6lc", "Figure 6 Left Contra")),
      fluidRow(actionButton("fig6li", "Figure 6 Left Ipsi")),
      fluidRow(actionButton("fig6rc", "Figure 6 Right Contra")),
      fluidRow(actionButton("fig6ri", "Figure 6 Right Ipsi"))
    )
  )

about.ui <- modalDialog(
  title = "About",
  easy_close = TRUE,
  "Retistruct was written by David Sterratt at the University of Edinburgh and
  tested by Daniel Lyngholm and Ian Thompson at the MRC Centre for
  Developmental Neurobiology, KCL. This work was supported by a Programme
  Grant from the Wellcome Trust (G083305).",
  br(),
  "Improvements to image handing and refactoring the code (released in
  v0.6.0) were supported by The Jackson Laboratory (Bar Harbor, ME, USA)
  Scientific Services Innovation Fund from 2016-2017 and an NIH R21
  Laboratory.",
  br(),
  "The capabilities to reconstruct tissue comprised of separate fragments
  (released in v0.7.0) and to reconstruct 3D data comprising an overhead
  image and depth map (released in v0.7.2), and user interface
  improvements (released in v0.7.0) were supported by an NIH R21 grant
  (EY027894) from 2018-2020 to Dr. Mark P. Krebs, The Jackson
  Laboratory.",
  br(),
  "Jan Okul converted the user interface to Shiny as part of his
  undergraduate project in the University of Edinburgh in 2024-2025
  (released in v0.8.0).",
  br(),br(),
  paste("Version:", version.string())
 )

## UI for file management
title.bar <- fluidRow(
  "Retistruct",
  shinyDirButton("open", "Open", "Select a directory...", icon=icon("folder-open")),
  actionButton("save", "Save", icon=icon("floppy-disk")),
  actionButton("reconstruct", "Reconstruct", icon=icon("hammer")),
  actionButton("properties", "Properties", icon=icon("gear")),
  actionButton("demo", "Demo"),
  actionButton("about", "About"),
  actionButton("cancel", "Cancel")
)

## The "edit" tab of the sidebar
edit.panel <- tabPanel(
  "Edit",
  actionButton("add_tear", "Add tear", style="width:100%"),
  br(),
  actionButton("remove_tear", "Remove tear", style="width:100%"),
  br(),
  actionButton("add_fullcut", "Add full cut", style="width:100%"),
  br(),
  actionButton("remove_fullcut", "Remove full cut", style="width:100%"),
  br(),
  actionButton("move_point", "Move point", style="width:100%"),
  br(),
  actionButton("mark_nasal", "Mark nasal", style="width:100%"),
  br(),
  actionButton("mark_dorsal", "Mark dorsal", style="width:100%"),
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
)

## The "view" section of the sidebar
view.panel <- tabPanel(
  "View",
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
    column(6, numericInput("center.el", "El", value = 0, min = el.low, max = el.high)),
    column(6, numericInput("center.az", "Az", value = 0, min = az.low, max = az.high))
  ),
  selectInput("transform", "Transform", choices = getTransforms()),
  tags$strong(textOutput("axdir")),
  fluidRow(
    # Width is 6
    column(6, numericInput("ax.el", "El", value = 90, min = el.low, max = el.high)),
    column(6, numericInput("ax.az", "Az", value = 0, min = az.low, max = az.high))
  ),
  checkboxGroupInput("ids", "IDs", choices=c("All"), selected =c("All")),
)

sidebar.ui <- sidebarPanel(
  tabsetPanel(edit.panel, view.panel),
  width = 2 # Sidebar width
)

## UI for output panel for plots to be displayed on
main.ui <- mainPanel(
  navset_tab(
    nav_panel(
      "Main Panel",
      fluidRow(
        column(
          width=6,
          shinySaveButton("bitmap1", "BitMap", "Select file to save to...",
                          filetype=picture.types),
          shinySaveButton("pdf1", "PDF", "Select file to save to..." ,
                          filetype=pdf.type),
          plotOutput("plot1", width="100%", height="80vh", fill=TRUE,
                     click="plot1click")
        ),
        column(
          width=6,
          shinySaveButton("bitmap2", "BitMap", "Select file to save to...",
                          filetype=picture.types),
          shinySaveButton("pdf2", "PDF", "Select file to save to...",
                          filetype=pdf.type),
          plotOutput("plot2", width="100%", height="80vh", fill=TRUE)
        ),
      )
    ),
    nav_panel(
      "3D Panel",
      rglwidgetOutput("plot3"),
    )
  ),
  width=10
)

##' @title Retistruct UI
##' @author Jan Okul
##' @description
##' The Shiny UI element, runs on a browser and is similar to HTML, attempted to
##' mimic the original Retistruct UI as closely as possible.
##' @importFrom shinyjs useShinyjs
##' @import shiny
##' @importFrom shinyFiles shinyDirButton shinySaveButton
##' @importFrom bslib navset_tab nav_panel
ui <- function() {
  fluidPage(
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
  
    titlePanel(windowTitle=version.string(), title.bar),
    sidebarLayout(sidebar.ui, main.ui),
    textOutput("status")
  )
}
