##' Check the whether  directory contains valid data 
##' @param dir Directory to check.
##' @return  \code{TRUE} if \code{dir} contains valid data;
##' \code{FALSE} otherwise.
##' @author David Sterratt
##' @export
checkDatadir <- function(dir=NULL) {
  if (idt.checkDatadir(dir))   { return("idt") }
  if (csv.checkDatadir(dir))   { return("csv") }
  if (ijroi.checkDatadir(dir)) { return("ijroi") }
  return(FALSE)
}

##' Read a retinal dataset in one of three formats; for information on
##' formats see \code{\link{idt.read.dataset}},
##' \code{\link{csv.read.dataset}} and
##' \code{\link{ijroi.read.dataset}}. The format is autodetected from
##' the files in the directory.
##' 
##' @title Read a retinal dataset
##' @param dataset Path to directory containing the files
##'   corresponding to each format.
##' @param report Function to report progress. Set to \code{FALSE} for
##'   no reporting.
##' @param ... Parameters passed to the format-specific functions.
##' @return A \code{\link{RetinalOutline}} object
##' @author David Sterratt
##' @export
retistruct.read.dataset <- function(dataset, report=message, ...) {
  ## Check to see if dataset is valid
  type <- checkDatadir(dataset)
  
  if (!is.function(report)) {
    if (report != FALSE) {
      stop("report must be FALSE or a function")
    }
    report <- function(...) {}
  }

  if (type=="idt")   { return(idt.read.dataset(dataset, report, ...))}
  if (type=="csv")   { return(csv.read.dataset(dataset, report, ...))}
  if (type=="ijroi") { return(ijroi.read.dataset(dataset, report, ...))}

  stop("No valid dataset format detected.")
}

##' Read the markup data contained in the files \file{markup.csv},
##' \file{P.csv} and \file{T.csv} in the directory \file{dataset},
##' which is specified in the reconstruction object \code{r}.
##'
##' The tear information is contained in the files \file{P.csv} and
##' \file{T.csv}. The first file contains the locations of outline
##' points that the tears were marked up on. The second file contains
##' the indices of the apex and backward and forward vertices of each
##' tear. It is necessary to have the file of points just in case the
##' algorithm that determines \code{P} in
##' \code{\link{retistruct.read.dataset}} has changed since the markup
##' of the tears.
##'
##' The remaining information is contained in the file
##' \file{markup.csv}.
##'
##' If \code{DVflip} is specified, the locations of points \code{P}
##' flipped in the \eqn{y}-direction. This operation also requires the
##' swapping of \code{gf}  and \code{gb} and \code{VF} and \code{VB}.
##' @title Read the markup data
##' @param a Dataset object, containing \code{dataset} path
##' @param error Function to run on error, by default \code{stop()}
##' @return o \code{RetinalDataset} object
##' \item{V0}{Indices in \code{P} of apices of tears}
##' \item{VB}{Indices in \code{P} of backward vertices of tears}
##' \item{VF}{Indices in \code{P} of backward vertices of tears}
##' \item{iN}{Index in \code{P} of nasal point, or \code{NA} if not marked}
##' \item{iD}{Index in \code{P} of dorsal point, or \code{NA} if not marked}
##' \item{iOD}{Index in \code{Ss} of optic disc }
##' \item{phi0}{Angle of rim in degrees}
##' \item{DVflip}{Boolean variable indicating if dorsoventral (DV) axis has been flipped}
##' @author David Sterratt
##' @importFrom utils read.csv 
##' @export
retistruct.read.markup <- function(a, error=stop) {
  ## Return index in P of closest point to x
  closest <- function(P, x) {
    if (any(is.na(x))) {
      return(NA)
    }
    return(which.min(vecnorm(t(t(P) - x))))
  }

  ## Function to map old tears (M.old) and old points (P.old) onto new
  ## points (P). It returns a new matrix of indices (M).
  convert.markup <- function(M.old, P.old, P) {
    M <- sapply(as.matrix(M.old), function(i) {
      ifelse(is.numeric(i), closest(P, P.old[i,]), i)
    })
    return(M)
  }

  ## Read in the old P data
  Pfile <- file.path(a$dataset, "P.csv")
  if (file.exists(Pfile)) {
    P.old <- as.matrix(read.csv(Pfile))
  } else {
    P.old <- a$getPoints()
  }
  
  ## Read in markup file
  markupfile <- file.path(a$dataset, "markup.csv")
  if (file.exists(markupfile)) {
    ## read.csv() returns a data frame
    M.df <- read.csv(markupfile)
    if (nrow(M.df) > 1) {
      stop("Markup file markup.csv contains more than one row")
    }
    M <- as.list(M.df)
    ## Get strings out
    if ("side" %in% names(M)) {
      a$side <- M[["side"]]
      ## Make sure that side is a valid value
      if (!(a$side %in% c("Left", "Right"))) {
        a$side <- "Right"
        warning("The side was neither Left nor Right; setting to Right.")
      }
    }
    if (!is.na(M[["iD"]])) {
      M[["iD"]] <- convert.markup(M[["iD"]], P.old, a$getPoints())
      a$setFixedPoint(M[["iD"]], "Dorsal")
    }
    if (!is.na(M[["iN"]])) {
      M[["iN"]] <- convert.markup(M[["iN"]], P.old, a$getPoints())
      a$setFixedPoint(M[["iN"]], "Nasal")
    }
    a$phi0 <- M[["phi0"]]*pi/180
    if ("iOD" %in% names(M)) {
      a$getFeatureSet("LandmarkSet")$setID(M[["iOD"]], "OD")
    }
    if ("DVflip" %in% names(M)) {
      a$DVflip <- M[["DVflip"]]
    }
  } else {
    error("Markup file markup.csv doesn't exist.")
  }
  
  ## Read in tearfile
  tearfile <- file.path(a$dataset, "T.csv")
  if (file.exists(tearfile)) {
    T.old <- read.csv(tearfile)
    cn <- colnames(T.old)
    T <- matrix(convert.markup(as.matrix(T.old), P.old, a$getPoints()), ncol=3)
    colnames(T) <- cn
    if (nrow(T) > 0) {
      for (i in 1:nrow(T)) {
        a$addTear(T[i,])
      }
    }
  } else {
    error("Tear file T.csv doesn't exist.")
  }
  return(a)
}

##' Check that markup such as tears and the nasal or dorsal points are present.
##' 
##' @title Retistruct check markup
##' @param o Outline object
##' @return If all markup is present, return \code{TRUE}. Otherwise
##' return \code{FALSE}.
##' @author David Sterratt
##' @export
retistruct.check.markup <- function(o) {
  return(!is.null(names(o$i0)))
}

##' Given an outline object with a \code{dataset} field, read the
##' reconstruction data from the file \file{\var{dataset}/r.Rdata}.
##'
##' @title Read the reconstruction data from file
##' @param o Outline object containing \code{dataset} field
##' @param check If \code{TRUE} check that the base information in the
##' reconstruction object is the same as the base data in source
##' files.
##' @return If the reconstruction data exists, return a reconstruction
##' object, else return the outline object \code{o}.
##' @author David Sterratt
##' @export
retistruct.read.recdata <- function(o, check=TRUE) {
  ## FIXME: Issue #27: Saving and reading recddata need to be reviewed
  r <- NULL
  recfile <- file.path(o$dataset, "r.Rdata")
  if (file.exists(recfile)) {
    r_list <- NULL
    load(recfile)                # This puts r_list in the environment
    ## If the algorithm in the codebase is newer than in the recdata
    ## file, reject the recfile data
    r <- list_to_R6(r_list)
    if (is.null(r$version) || (r$version != OutlineCommon$new()$version)) {
      unlink(recfile)
      warning("The algorithm has changed significantly since this retina was last reconstructed, so the cached reconstruction data has been deleted.")
      return(NULL)
    }
    ## If the base data doesn't match the recfile data, reject the
    ## recfile data
    if (check) {
      chk <- all.equal(o, r$ol0)
      if (!isTRUE(chk)) {
        print(chk)
        unlink(recfile)
        warning("The base data has changed since this retina was last reconstructed, so the cached reconstruction data has been deleted.")
        return(NULL)
      } 
    }
    
    ## Make sure the dataset information isn't overwritten
    ## r$dataset <- o$dataset
    return(r)
  }
  return(NULL)
}

##' Reconstruct a retina
##' @param a \code{\link{RetinalOutline}} object with tear and
##'   correspondence annotations
##' @param report Function to report progress. Set to \code{FALSE} for
##'   no reporting or to \code{NULL} to inherit from the argument given to \code{\link{retistruct.read.dataset}}
##' @param plot.3d If \code{TRUE} show progress in a 3D plot
##' @param dev.flat The ID of the device to which to plot the flat
##'   representation
##' @param dev.polar The ID of the device to which to plot the polar
##'   representation
##' @param debug If \code{TRUE} print extra debugging output
##' @param ... Parameters to be passed to
##'   \code{\link{RetinalReconstructedOutline}} constructor
##' @return A \code{\link{RetinalReconstructedOutline}} object
##' @author David Sterratt
##' @export
retistruct.reconstruct <- function(a, report=NULL,
                                   plot.3d=FALSE, dev.flat=NA, dev.polar=NA,
                                   debug=FALSE, ...) {
  o <- a$clone()
  
  ## Check that markup is there
  if (!retistruct.check.markup(o)) {
    stop("Neither dorsal nor nasal pole specified")
  }

  ## Check tears are valid
  ct <- o$checkTears()
  if (length(ct)) {
    stop(paste("Invalid tears", toString(ct), "marked up. Fix using \"Move Point\"."))
  }

  ## Set up fixed point 
  o$lambda0 <- 0
  ## Case of dorsal point specified...
  if (names(o$i0)[1]=="Dorsal") {
    o$lambda0 <- 90 * pi/180
  }
  if (names(o$i0)[1]=="Nasal") {
    if (o$side=="Right") {
      o$lambda0 <- 0
    } else {
      o$lambda0 <- pi
    }
  }

  ## Now do folding itself
  r <- NULL
  r <- RetinalReconstructedOutline$new()
  r$loadOutline(o, debug=debug)

  r$reconstruct(plot.3d=plot.3d, dev.flat=dev.flat,
                dev.polar=dev.polar, report=report,
                ...)
  if (!is.null(r)) {
    repstr <- paste("Mapping optimised. Deformation eL:", format(sqrt(r$E.l), 5),
                    ";", r$nflip, "flipped triangles.")
    if (is.null(r$EOD)) {
      repstr <- paste(repstr, "OD not marked up.")
    } else {
      repstr <- paste(repstr, "OD displacement:",
                      format(r$EOD, 2),
                      "degrees.")
    }
    report(repstr)      
  }
  return(r)
}

##' Save the markup in the \code{\link{RetinalOutline}} \code{a} to a
##' file called \code{markup.csv} in the directory \code{a$dataset}.
##'
##' @title Save markup
##' @param a \code{\link{RetinalOutline}} object
##' @author David Sterratt
##' @importFrom utils write.csv 
##' @export
retistruct.save.markup <- function(a) {
  ## Save the tear information and the outline
  write.csv(a$getTears(), file.path(a$dataset, "T.csv"),
            row.names=FALSE)
  write.csv(a$getPoints(), file.path(a$dataset, "P.csv"),
            row.names=FALSE)
      
  ## Save the dorsal and nasal locations and phi0 to markup.csv
  i0 <- a$getFixedPoint()
  iD <- ifelse(names(i0) == "Dorsal", i0, NA)
  iN <- ifelse(names(i0) == "Nasal" , i0, NA)

  iOD <- a$getFeatureSet("LandmarkSet")$getIndex("OD")
  if (length(iOD) == 0)
    iOD <- NA
  markup <- data.frame(iD=iD, iN=iN, phi0=a$phi0*180/pi, iOD=iOD, DVflip=a$DVflip, side=a$side)     
  write.csv(markup, file.path(a$dataset, "markup.csv"), row.names=FALSE)
}


##' Save the reconstruction data in an object \code{r} of class
##' \code{\link{RetinalReconstructedOutline}} to a file called
##' \code{r.Rdata} in the directory \code{r$dataset}.
##'
##' @title Save reconstruction data
##' @param r \code{\link{RetinalReconstructedOutline}} object
##' @author David Sterratt
##' @export
retistruct.save.recdata <- function(r) {
  ## FIXME: Issue #27: Saving and reading recddata need to be reviewed
  if ("RetinalReconstructedOutline" %in% class(r)) {
    ## Save the derived data
    ## r$version <- recfile.version        # Datafile version
    r_list <- R6_to_list(r)
    if (!is.null(r_list)) {
      save(r_list, file=file.path(r$ol$dataset, "r.Rdata"))
    }
  }
}

##' Save as a MATLAB object certain fields of an object \code{r} of
##' class\code{\link{RetinalReconstructedOutline}} to a file called
##' \code{r.mat} in the directory \code{r$dataset}.
##'
##' @title Save reconstruction data in MATLAB format
##' @param r \code{\link{RetinalReconstructedOutline}} object
##' @param filename Filename of output file. If not specified, is
##'   \code{r.mat} in the same directory as the input files
##' @author David Sterratt
##' @export
retistruct.export.matlab <- function(r, filename=NULL) {
  ## writeMat() doesn't seem to cope with hierarchical structures, so
  ## we have to unlist the KDE and KR objects using this function
  unlist.kernel.estimate <- function(KDE) {
    if (length(KDE) > 0) {
      for (i in 1:length(KDE)) {
        KDEi <- KDE[[i]]
        KDEi$flevels <- NULL
        KDEi$labels <- NULL
        KDEi$tot.contour.areas <- NULL
        KDEi <- unlist(KDEi, recursive=FALSE)
        KDE[[i]] <- list(flevels=KDE[[i]]$flevels, labels=KDE[[i]]$labels,
                         tot.contour.areas=KDE[[i]]$tot.contour.areas)
        KDE[[i]] <- c(KDE[[i]], KDEi)
        names(KDE[[i]]$contour.areas) <- c()
      }
      ## print(names(KDE))
      ## print(names(unlist(KDE, recursive=FALSE)))
      KDE <- unlist(KDE, recursive=FALSE)
      names(KDE) <- gsub('\\.', '_', names(KDE))
    }
    return(KDE)
  }

  if (!is.null(r)) {
    if (is.null(filename)) {
      filename <- file.path(r$ol$dataset, "r.mat")
    }
    message(paste("Saving", filename))
    KDE <- NULL
    if (!is.null(r$getFeatureSet("PointSet"))) {
      KDE <- unlist.kernel.estimate(r$getFeatureSet("PointSet")$getKDE())
    }
    KR <- NULL
    if (!is.null(r$getFeatureSet("CountSet"))) {
      KR <-  unlist.kernel.estimate(r$getFeatureSet("CountSet")$getKR())
    }
    R.matlab::writeMat(filename,
                       phi0=r$ol$phi0*180/pi,
                       Dss=r$getFeatureSet("PointSet")$Ps,
                       DssMean=r$getFeatureSet("PointSet")$getMean(),
                       DssHullarea=r$getFeatureSet("PointSet")$getHullarea(),
                       Sss=r$getFeatureSet("PointSet")$Ps,
                       Tss=r$getTearCoords(),
                       KDE=KDE,
                       KR=KR,
                       side=as.character(r$ol$side), DVflip=r$ol$DVflip
                       ## FIXME: Issue #25: Implement Wedge coords
                       ## Dsw=lapply(r$Dsc, function(x) {sphere.cart.to.sphere.wedge(x, r$phi0 + pi/2, r$R)}),
                       ## Dsdw=lapply(r$Dsc, function(x) {sphere.cart.to.sphere.dualwedge(x, r$phi0 + pi/2, r$R)})
                       )
  }
}

