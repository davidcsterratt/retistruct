read.scale <- function(dataset) {
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
  return(scale)
}

read.image <- function(dataset) {
  im <- NULL
  imfile <- file.path(dataset, "image.png")
  if (file.exists(imfile)) {
    message("Reading image")
    im <- as.raster(readPNG(imfile))
  }
  return(im)
}
