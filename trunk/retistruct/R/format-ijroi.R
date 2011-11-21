ijroi.check.datadir <- function(dir=NULL) {
  return(file.exists(file.path(dir, "outline.roi")))
}

##' Read a retinal dataset in IJROI format. Each dataset is a folder
##' containing a file called outline.roi that specifies the outline in
##' X-Y coordinates.
##' 
##' @title Read a retinal dataset in IJROI format
##' @param dataset Path to directory containing \code{outline.roi}
##' @return A \code{RetinalDataset} object
##' \item{dataset}{The path to the directory given as an argument}
##' \item{raw}{List containing\describe{
##'    \item{\code{outline}}{The raw outline.roi data}
##' }}
##' \item{P}{The points of the outline}
##' \item{gf}{Forward pointers along the outline}
##' \item{gb}{Backward pointers along the outline}
##' \item{Ds}{List of datapoints}
##' \item{Ss}{List of landmark lines}
##' }
##' @author David Sterratt
ijroi.read.dataset <- function(dataset) {
  ## Read the raw data
  roi <- read.ijroi(file.path(dataset, "outline.roi"))
  out <- roi$coords

  ## Read scale
  scale <- read.scale(dataset)
  
  ## If there is an image, read it
  im <- read.image(dataset)
  
  ## Extract datapoints
  ##
  ## At present, for the plotting functions to work, the name of each
  ## group has to be a valid colour. There are no datapoints in this
  ## format, but we may have landmarks.
  Ds <- list()
  cols <- list(OD="blue",
               default="orange")

  ## ImageJ ROI format plots has the coordinate (0, 0) in the top
  ## left.  We have the coordinate (0, 0) in the bottom left. We need
  ## to transform P so that the outline appears in the correct
  ## orientation.
  P <- out
  offset <- ifelse(is.null(im), max(P[,2]), nrow(im))
  P[,2] <- offset - P[,2]

  Ss <- list()

  ## Read in an Optic Disc. FIXME: this should actually be marked as
  ## the Optic Disc
  od.file <- file.path(dataset, "od.roi")
  if (file.exists(od.file)) {
    roi <- read.ijroi(od.file)
    out <-  roi$coords
    offset <- ifelse(is.null(im), max(out[,2]), nrow(im))
    out[,2] <- offset - out[,2]
    Ss[["OD"]] <- out
  }
  
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

