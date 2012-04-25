csv.checkDatadir <- function(dir=NULL) {
  return(file.exists(file.path(dir, "outline.csv")))
}

##' Read a retinal dataset in CSV format. Each dataset is a folder
##' containing a file called outline.csv that specifies the outline in
##' X-Y coordinates. It may also contain a file \code{datapoints.csv},
##' containing the locations of data points; see
##' \code{\link{read.datapoints}} for the format of this file.
##' 
##' @title Read a retinal dataset in CSV format
##' @param dataset Path to directory containing \code{outline.csv}
##' @return A \code{RetinalDataset} object
##' @author David Sterratt
csv.read.dataset <- function(dataset) {
  ## Read the raw data
  out <- read.csv(file.path(dataset, "outline.csv"))

  ## Read scale
  scale <- read.scale(dataset)
  
  ## If there is an image, read it
  im <- read.image(dataset)

  ## Extract datapoints
  ##
  ## At present, for the plotting functions to work, the name of each
  ## group has to be a valid colour.
  Ds <- list()
  cols <- list()
  dat <- read.datapoints(dataset)
  Ds <- c(Ds, dat$Ds)
  cols <- c(cols, dat$cols)

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

