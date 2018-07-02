csv.checkDatadir <- function(dir=NULL) {
  return(file.exists(file.path(dir, "outline.csv")))
}

##' Read a retinal dataset in CSV format. Each dataset is a folder
##' containing a file called outline.csv that specifies the outline in
##' X-Y coordinates. It may also contain a file \code{datapoints.csv},
##' containing the locations of data points and a file
##' \code{datacounts.csv}, containing the locations of data counts;
##' see \code{\link{read.datapoints}} and
##' \code{\link{read.datacounts}} for the formats of these files. The
##' folder may also contain a file \code{od.csv} specifying the
##' coordinates of the optic disc.
##' 
##' @title Read a retinal dataset in CSV format
##' @param dataset Path to directory containing \code{outline.csv}
##' @param report Function to report progress
##' @return A \code{\link{RetinalOutline}} object
##' @author David Sterratt
csv.read.dataset <- function(dataset, report=message) {
  ## Read the raw data
  out <- as.matrix(read.csv(file.path(dataset, "outline.csv")))

  ## Read scale
  scale <- read.scale(dataset)
  
  ## If there is an image, read it
  im <- read.image(dataset, report=report)
  
  ## ImageJ ROI format plots has the coordinate (0, 0) in the top
  ## left.  We have the coordinate (0, 0) in the bottom left. We need
  ## to transform P so that the outline appears in the correct
  ## orientation.
  P <- out
  offset <- ifelse(is.null(im), max(P[,2]), nrow(im))
  P[,2] <- offset - P[,2]
  
  ## Extract datapoints
  ##
  ## At present, for the plotting functions to work, the name of each
  ## group has to be a valid colour. There are no datapoints in this
  ## format, but we may have landmarks.
  Ds <- list()

  cols <- c(OD="blue",
            default="orange")
  dat <- read.datapoints(dataset)
  Ds <- c(Ds, dat$Ds)
  cols <- c(cols, dat$cols)
  Ds <- lapply(Ds, function(P) {cbind(X=P[,1], Y=(offset - P[,2]))})
  
  ## Extract landmarks (currently optic disc)
  Ss <- list()

  ## Read in an Optic Disc. FIXME: this should actually be marked as
  ## the Optic Disc
  od.file <- file.path(dataset, "od.csv")
  if (file.exists(od.file)) {
    out <-  as.matrix(read.csv(od.file))
    out[,2] <- offset - out[,2]
    colnames(out) <- c("X", "Y")
    Ss[["OD"]] <- out
  }

  ## Extract datapoints
  ##
  ## At present, for the plotting functions to work, the name of each
  ## group has to be a valid colour. There are no datapoints in this
  ## format, but we may have landmarks.
  Gs <- list()
  dat <- read.datacounts(dataset)
  Gs <- c(Gs, dat$Gs)
  cols <- c(cols, dat$cols)
  Gs <- lapply(Gs, function(P) {cbind(P[,1], offset - P[,2], P[,3])})
  
  ## Create forward and backward pointers
  o <- RetinalOutline$new(P, scale=scale["Scale"], im=im,
                          units=scale["Units"],
                          dataset=dataset, report=report)
  
  ## Check that P is more-or-less closed
  ## if (vecnorm(P[1,] - P[nrow(P),]) > (d.close * diff(range(P[,1])))) {
  ##    stop("Unable to find a closed outline.")
  ## }

  o$addFeatureSet(PointSet$new(data=Ds, cols=cols))
  o$addFeatureSet(LandmarkSet$new(data=Ss, cols=cols))
  return(o)
}
