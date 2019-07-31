read.scale <- function(dataset, report=message) {
  ## If there is a scale file, read it
  scale <- c(Scale=NA, Units=NA)
  scfile <- file.path(dataset, "scale.csv")
  if (file.exists(scfile)) {
    report("Reading scale file")
    sc <- read.csv(scfile)
    valid.colnames <- c("Scale", "Units")
    if (!all(colnames(sc) %in% valid.colnames)) {
      stop(paste("Unknown column names",
                 paste0("\"", setdiff(colnames(sc), valid.colnames), "\"",
                        collapse=", "),
                 "in", scfile, ". Valid column names:",
                 paste(valid.colnames, collapse=", ")))
    }
    scale <- as.matrix(sc)[1,]
    if (!("Scale" %in% names(scale)) | !is.numeric(scale["Scale"])) {
      stop("Scale file has not been read correctly. Check it is in the correct format.")
    }
    if (!("Units" %in% names(scale))) {
      scale["Units"] <- NA
    }
  } else {
    warning("Scale file \"scale.csv\" does not exist. Scale bar will not be set.")
  }
  return(scale)
}

read.image <- function(dataset, report=message) {
  im <- NULL
  imfile <- file.path(dataset, "image.png")
  if (file.exists(imfile)) {
    report("Reading image")
    im <- grDevices::as.raster(png::readPNG(imfile))
  }
  return(im)
}

## Copied from demo("colors")
## @title Comparing Colors
## @param col
## @param nrow
## @param ncol
## @param txt.col
## @return the grid layout, invisibly
## @author Marius Hofert, originally
plotCol <- function(col, nrow=1, ncol=ceiling(length(col) / nrow),
                    txt.col="black") {
    stopifnot(nrow >= 1, ncol >= 1)
    if(length(col) > nrow*ncol)
        warning("some colors will not be shown")
    grid::grid.newpage()
    gl <- grid::grid.layout(nrow, ncol)
    grid::pushViewport(grid::viewport(layout=gl))
    ic <- 1
    for(i in 1:nrow) {
        for(j in 1:ncol) {
            grid::pushViewport(grid::viewport(layout.pos.row=i, layout.pos.col=j))
            grid::grid.rect(gp= grid::gpar(fill=col[ic]))
            grid::grid.text(col[ic], gp=grid::gpar(col=txt.col))
            grid::upViewport()
            ic <- ic+1
        }
    }
    grid::upViewport()
    invisible(gl)
}

check.colour <- function(col) {
  if (!(col %in% grDevices::colours())) {
    plotCol(grep("([0-9]|medium|light|dark)",  grDevices::colors(), invert=TRUE, value=TRUE), nrow=20)
    return(FALSE)
  }
  return(TRUE)
}

##' Read data points from a file \code{dataponts.csv} in the directory
##' \code{dataset}. The CSV should contain two columns for every
##' dataset. Each pair of columns must contain a unique name in the
##' first cell of the first row  and a valid colour in the second
##' cell of the first row. In the remaining rows, the X coordinates of
##' data points should be in the first column and the Y coordinates
##' should be in the second column.
##'
##' @title Read data points in CSV format
##' @param dataset Path to directory containing \code{dataponts.csv}
##' @return List containing
##' \item{\code{Ds}}{List of sets of datapoints. Each set comprises a 2-column matrix and each set is named.}
##' \item{\code{cols}}{List of colours for each dataset. There is one element that corresponds to each element of \code{Ds} and which bears the same name.}
##' @author David Sterratt
read.datapoints <- function(dataset) {
  datfile <- file.path(dataset, "datapoints.csv")
  Ds <- list()
  cols <- c()
  if (file.exists(datfile)) {
    message("Reading datapoints")
    ## Read file. stringsAsFactors=FALSE prevents conversion to factors
    dat <- read.csv(file.path(datfile), stringsAsFactors=FALSE)

    ## Go through pairs of columns
    while(ncol(dat) >= 2) {
      ## Extract first two columns
      d <- dat[,1:2]
      dat <- dat[,-(1:2)]
      names <- colnames(d)

      ## Convert strings to numeric. Suppress warnings as sapply
      ## complains about coercion to NA
      suppressWarnings({d <- sapply(d, as.numeric, USE.NAMES=FALSE)})
      ## Force conversion to matrix, necessary when the data has only
      ## one row
      d <- matrix(d, ncol=2)
      
      ## Any strings (e.g. empty ones) that don't convert will be
      ## converted to NA. Get rid of these.
      d <- na.omit(d)
      attr(d, "na.action") <- NULL
      colnames(d) <- c("X", "Y")

      ## Add to lists with appropriate names
      
      D <- list(d)
      names(D) <- names[1]
      Ds <- c(Ds, D)

      col <- names[2]
      if (!(check.colour(col))) {
        stop("Invalid colour \"", col, "\" in datapoints.csv - see window for valid colour names")
      }
      names(col) <- names[1]
      cols <- c(cols, col)
    }
  }
  return(list(Ds=Ds, cols=cols))
}

##' Read data counts from a file \file{datacounts.csv} in the
##' directory \code{dataset}. The CSV file should contain two columns
##' for every dataset. Each pair of columns must contain a unique name
##' in the first cell of the first row and a valid colour in the
##' second cell of the first row. In the remaining rows, the X
##' coordinates of data counts should be in the first column and the Y
##' coordinates should be in the second column.
##'
##' @title Read data counts in CSV format
##' @param dataset Path to directory containing \code{dataponts.csv}
##' @return List containing
##' \item{\code{Ds}}{List of sets of data counts. Each set comprises a 2-column matrix and each set is named.}
##' \item{\code{cols}}{List of colours for each dataset. There is one element that corresponds to each element of \code{Ds} and which bears the same name.}
##' @author David Sterratt
read.datacounts <- function(dataset) {
  datfile <- file.path(dataset, "datacounts.csv")
  Gs <- list()
  cols <- c()
  if (file.exists(datfile)) {
    message("Reading datacounts")
    ## Read file. stringsAsFactors=FALSE prevents conversion to factors
    dat <- read.csv(file.path(datfile), stringsAsFactors=FALSE)

    ## Go through triples of columns
    while(ncol(dat) >= 3) {
      ## Extract first three columns
      d <- dat[,1:3]
      dat <- dat[,-(1:3)]
      names <- colnames(d)

      ## Convert strings to numeric. Suppress warnings as sapply
      ## complains about coercion to NA
      suppressWarnings({d <- sapply(d, as.numeric, USE.NAMES=FALSE)})
      ## Force conversion to matrix, necessary when the data has only
      ## one row
      d <- matrix(d, ncol=3)
      colnames(d) <- c("X", "Y", "C")
      
      ## Any strings (e.g. empty ones) that don't convert will be
      ## converted to NA. Get rid of these.
      d <- na.omit(d)
      attr(d, "na.action") <- NULL

      ## Add to lists with appropriate names
      
      G <- list(d)
      names(G) <- names[1]
      Gs <- c(Gs, G)

      col <- list(names[2])
      if (!(check.colour(col))) {
        stop("Invalid colour \"", col, "\" in datacounts.csv - see window for valid colour names")
      }
      names(col) <- names[1]
      cols <- c(cols, col)
    }
  }
  return(list(Gs=Gs, cols=cols))
}
