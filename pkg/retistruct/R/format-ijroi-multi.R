ijroimulti.checkDatadir <- function(dir=NULL) {
  if (file.exists(file.path(dir, "RoiSet.zip")) |
      file.exists(file.path(dir, "outline.roi"))) {
    return(TRUE)
  }
  return(FALSE)
}

##' Read a retinal dataset in IJROI format. Each dataset is a folder
##' containing a file called \code{outline.roi} that specifies the
##' outline in X-Y coordinates. It may also contain a file
##' \code{datapoints.csv}, containing the locations of data points;
##' see \code{\link{read.datapoints}} for the format of this file. The
##' folder may also contain a file \code{od.roi} specifying the
##' coordinates of the optic disc.
##'
##' @title Read a retinal dataset in IJROI format
##' @param dataset Path to directory containing \code{outline.roi}
##' @return A \code{\link{RetinalOutline}} object
##' @author David Sterratt
ijroimulti.read.dataset <- function(dataset) {
  ## Read the raw data
  if (file.exists(file.path(dataset, "RoiSet.zip"))) {
    rois <- RImageJROI::read.ijzip(file.path(dataset, "RoiSet.zip"))
  } else {
    roifile <- file.path(dataset, "outline.roi")
    if (file.exists(roifile)) {
      rois <- list(RImageJROI::read.ijroi(roifile))
    }
  }
  out <- lapply(rois, function(roi) {roi$coords})

  ## Read scale
  scale <- read.scale(dataset)

  ## If there is an image, read it
  im <- read.image(dataset)

  ## If there is a depth map read it
  dm <- read.depthmap(dataset)
  if (!is.null(dm)) {
    if (is.null(im)) {
      stop("There is a depthmap (depthmap.tif) but no image (image.png)")
    }
    if (!("Z" %in% names(scale))) {
      stop("\"Z\" must be specified in scale.csv when depthmap is present")
    }
  }

  ## ImageJ ROI format plots has the coordinate (0, 0) in the top
  ## left.  We have the coordinate (0, 0) in the bottom left. We need
  ## to transform P so that the outline appears in the correct
  ## orientation.
  P <- do.call(rbind, out)
  offset <- ifelse(is.null(im), max(P[,2]), nrow(im))
  Ps <- lapply(out, function(P) {
    P[,2] <- offset - P[,2]
    return(P)
  })

  ## Extract datapoints
  ##
  ## At present, for the plotting functions to work, the name of each
  ## group has to be a valid colour. There are no datapoints in this
  ## format, but we may have landmarks.
  Ds <- list()
  cols <- list(OD="blue",
               default="orange")
  dat <- read.datapoints(dataset)
  Ds <- c(Ds, dat$Ds)
  cols <- c(cols, dat$cols)
  Ds <- lapply(Ds, function(P) {cbind(P[,1], offset - P[,2])})

  Ss <- list()

  ## Read in an Optic Disc. FIXME: this should actually be marked as
  ## the Optic Disc
  od.file <- file.path(dataset, "od.roi")
  if (file.exists(od.file)) {
    roi <- RImageJROI::read.ijroi(od.file)
    out <-  roi$coords
    colnames(out) <- c("X", "Y")
    out[,2] <- offset - out[,2]
    Ss[["OD"]] <- out
  }

  ## Create forward and backward pointers
  o <- RetinalOutline$new(fragments=Ps, scale=scale[["XY"]], im=im,
                          scalez=scale[["Z"]], dm=dm,
                          dataset=dataset)
  # o <- simplify.outline(o)

  o$addFeatureSet(LandmarkSet$new(data=Ss, cols=cols))
  return(o)
}
