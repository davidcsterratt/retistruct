csv.check.datadir <- function(dir=NULL) {
  return(file.exists(file.path(dir, "outline.csv")))
}

##' Read a retinal dataset in CSV format. Each dataset is a folder
##' containing a file called outline.csv that specifies the outline in
##' X-Y coordinates.
##' 
##' @title Read a retinal dataset in CSV format
##' @param dataset Path to directory containing \code{outline.csv}
##' @return 
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

  ## Extract datapoints
  ##
  ## At present, for the plotting functions to work, the name of each
  ## group has to be a valid colour.
  Ds <- list()
  cols <- list()

  ## The outline (P) is the longest connected segment and the outline
  ## is removed from the list of segments
  P  <- out
  Ss <- list()
  
  ## Create forward and backward pointers
  o <- Outline(P)
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

