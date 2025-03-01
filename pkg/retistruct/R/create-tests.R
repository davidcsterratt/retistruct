##' Morph a flat dataset to a parabola for testing purposes
##' @param orig.dataset Directory of original dataset
##' @param morphed.dataset Directory to write morphed dataset to. If
##'   NA, a temporary directory will be created
##' @param f Focus of parabola
##' @return Path to \code{morphed.dataset}
##' @author David Sterratt
##' @export
morph.dataset.to.parabola <- function(orig.dataset=file.path(system.file(package = "retistruct"), "extdata", "smi32-csv"),
                                      morphed.dataset=NA,
                                      f=100) {

  ## Set morphed.dataset if necessary
  if (is.na(morphed.dataset)) {
    morphed.dataset <- tempfile()
  }

  ## Read outline
  ## P <- as.matrix(read.csv(file.path(orig.dataset, "P.csv")))
  o <- retistruct.read.dataset(orig.dataset)
  P <- o$getPoints()
  P <- P[,c("X", "Y")]

  ## Points should have 0,0 at bottom left
  im <- NULL
  ## offset <- ifelse(is.null(im), max(P[,2]), nrow(im))
  ## P[,2] <- offset - P[,2]

  ## Read optic disk
  ## od <- as.matrix(read.csv(file.path(orig.dataset, "od.csv")))
  ## od[,2] <- offset - od[,2]

  ## Make grid of points
  xs <- seq(from=min(P[,1]) + 2, to=max(P[,1]) - 2, length.out=20)
  ys <- seq(from=min(P[,2]) + 2, to=max(P[,2]) - 2, length.out=20)

  Ds <- cbind(as.vector(outer(ys*0, xs, FUN="+")),
              as.vector(outer(ys, xs*0, FUN="+")))
  colnames(Ds) <- c("Test", "red")

  ## Parabolic transform of outline. s will be the path length along
  ## the parabola, but is now the y-coordinate of each point
  s <- P[,2] - min(P[,2])

  ## The parabola should be centred, so limits should be determined
  ## from the y-coordinate corresponding to half the arc
  y2 <- parabola.invarclength(0, max(s)/2, f)
  y1 <- -y2
  ## Want centre (i.e. Z = 0) to be on 0.5, to be in register with image
  ycent <- ceiling(y2 - 0.5) + 1
  ## Transformed Pt
  Pt <- cbind(P[,1], parabola.invarclength(y1, s, f) + ycent)

  ## Transform optic disk
  ## s <- od[,2] - min(P[,2])
  ## od[,2] <- parabola.invarclength(y1, s, f) + y2 + 1

  ## Transform points
  s <- Ds[,2] - min(P[,2])
  Dst <- cbind(Ds[,1], parabola.invarclength(y1, s, f) + ycent)

  ## Create a parabolic depthmap image
  Mt <- ceiling(max(Pt[,2]))              # Height
  Nt <- ceiling(max(Pt[,1]))              # Width
  ## a$Z <- (a$P[,2]^2)/4/f

  dm <- matrix(((Mt:1 - 0.5) - ycent)^2/4/f, Mt, Nt)
  Zmax <- max(dm)
  scalez <- Zmax/255
  dm <- dm/max(dm)

  ## Transform points so that they (0, 0) is at top left
  Pt[,2] <- Mt - Pt[,2]
  Dst[,2] <- Mt - Dst[,2]
  colnames(Dst) <- colnames(Ds)

  ## Write the dataset
  dir.create(morphed.dataset, showWarnings=FALSE)
  imtiff <- tiff::writeTIFF(dm, file.path(morphed.dataset, "depthmap.tif"), reduce=TRUE)
  ## png::writePNG(im, file.path(morphed.dataset, "image.png"))

  write.csv(Pt, file.path(morphed.dataset, "outline.csv"), row.names=FALSE)
  ## write.csv(od, file.path(morphed.dataset, "od.csv"), row.names=FALSE)
  write.csv(Dst, file.path(morphed.dataset, "datapoints.csv"), row.names=FALSE)
  file.copy(file.path(orig.dataset, "T.csv"), morphed.dataset, overwrite=TRUE)
  markup <- read.csv(file.path(orig.dataset, "markup.csv"))
  markup[1,"iOD"] <- NA
  write.csv(markup, file.path(morphed.dataset, "markup.csv"), row.names=FALSE)
  scale <- data.frame(XY=1, Z=scalez, Units="um")
  write.csv(scale, file.path(morphed.dataset, "scale.csv"), row.names=FALSE)
  return(morphed.dataset)
}
