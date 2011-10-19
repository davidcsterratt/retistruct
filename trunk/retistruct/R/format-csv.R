csv.check.datadir <- function(dir=NULL) {
  return(file.exists(file.path(dir, "outline.csv")))
}

##' Read a retinal dataset in CSV format. Each dataset is a folder
##' containing a file called outline.csv that specifies the outline in
##' X-Y coordinates.
##' 
##' @title Read a retinal dataset in CSV format
##' @param dataset Path to directory containing \code{outline.csv}
##' @return A \code{RetinalDataset} object
##' \item{dataset}{The path to the directory given as an argument}
##' \item{raw}{List containing\describe{
##'    \item{\code{outline}}{The raw outline.csv data}
##' }}
##' \item{P}{The points of the outline}
##' \item{gf}{Forward pointers along the outline}
##' \item{gb}{Backward pointers along the outline}
##' \item{Ds}{List of datapoints}
##' \item{Ss}{List of landmark lines}
##' }
##' @author David Sterratt
csv.read.dataset <- function(dataset) {
  ## Read the raw data
  out <- read.csv(file.path(dataset, "outline.csv"))

  ## If there is a scale file, read it
  scale <- 0
  scfile <- file.path(dataset, "scale.csv")
  if (file.exists(scfile)) {
    print("Reading scale file")
    sc <- read.csv(file.path(dataset, "scale.csv"))
    scale <- sc[1,1]
    if (!is.numeric(scale)) {
      stop("Scale file has not been read correctly. Check it is in the correct format.")
    }
  } else {
    warning("Scale file does not exist. Scale bar will not be set.")
  }
  
  ## If there is an image, read it
  im <- NULL
  imfile <- file.path(dataset, "image.png")
  if (file.exists(imfile)) {
    print("Reading image")
    im <- as.raster(readPNG(imfile))
  }
  
  ## Extract datapoints
  ##
  ## At present, for the plotting functions to work, the name of each
  ## group has to be a valid colour.
  Ds <- list()
  cols <- list()

  ## The outline (P) is the longest connected segment and the outline
  ## is removed from the list of segments
  P  <- out
  if (!is.null(im)) {
    message("Image present, so flipping P to align coordinates")
    P[,2] <- nrow(im) - P[,2] + 1
  }
  Ss <- list()
  
  ## Create forward and backward pointers
  o <- Outline(P, scale, im)
  o <- simplify.outline(o)
  
  ## Check that P is more-or-less closed
  ## if (vecnorm(P[1,] - P[nrow(P),]) > (d.close * diff(range(P[,1])))) {
  ##    stop("Unable to find a closed outline.")
  ## }

  d <- Dataset(o, dataset, Ds, Ss, cols=cols, raw=list(outline=out))
  a <- AnnotatedOutline(d)
  a <- RetinalDataset(a)
  return(a)
}

