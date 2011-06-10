## Global variables
recfile.version <- 3      # Version of reconstruction file data format

## Report function, with similar arguments to print
retistruct.report <- function(message, title="",...) {
  cat(paste(message, "\n", sep=""))
}

##' Read one of the Thompson lab's retinal datasets. Each dataset is a
##' folder containing a SYS file in SYSTAT format and a MAP file in
##' text format. The SYS file specifies the locations of the data
##' points and the MAP file specifies the outline. 
##'
##' The function returns the outline of the retina. In order to do so,
##' it has to join up the segments of the MAP file. The tracings are
##' not always precise; sometimes there are gaps between points that
##' are actually the same point. The parameter \code{d.close} specifies
##' how close points must be to count as the same point.
##' 
##' @title Read one of the Thompson lab's retinal datasets
##' @param dataset Path to directory containing as SYS and MAP file
##' @param d.close Maximum distance between points for them to count
##' as the same point. This is expressed as a fraction of the width of
##' the outline.
##' @return 
##' \item{dataset}{The path to the directory given as an argument}
##' \item{raw}{List containing\describe{
##'    \item{\code{map}}{The raw MAP data}
##'    \item{\code{sys}}{The raw SYS data}
##' }}
##' \item{P}{The points of the outline}
##' \item{gf}{Forward pointers along the outline}
##' \item{gb}{Backward pointers along the outline}
##' \item{Ds}{List of datapoints}
##' \item{Ss}{List of landmark lines}
##' }
##' @author David Sterratt
retistruct.read.dataset <- function(dataset, d.close=0.25) {
  ## Check to see if dataset is valid
  check.datadir(dataset)
  
  ## Read the raw data
  map <- read.map(dataset)
  sys <- read.sys(dataset)

  ## Extract datapoints
  ##
  ## At present, for the plotting functions to work, the name of each
  ## group has to be a valid colour.
  Ds <- list(green =cbind(na.omit(sys[,'XGREEN']), na.omit(sys[,'YGREEN'])),
             red   =cbind(
               na.omit(c(sys[,'XDOUBLE'], sys[,'XRED'])),
               na.omit(c(sys[,'YDOUBLE'], sys[,'YRED']))),
             double=cbind(na.omit(sys[,'XDOUBLE']),na.omit(sys[,'YDOUBLE'])))
  cols <- list(green="green",
               red="red",
               double="yellow",
               OD="blue",
               default="orange")
  
  ## Extract line data
  segs <- map.to.segments(map)

  ## Connect together segments that look to be joined
  Ss <- connect.segments(segs)

  ## Determine the lengths of the segments
  l <- c()
  for (i in 1:length(Ss)) {
    l[i] <- segment.length(Ss[[i]])
  }

  ## The outline (P) is the longest connected segment and the outline
  ## is removed from the list of segments
  P  <- Ss[[which.max(l)]]
  Ss <- Ss[-which.max(l)]
  ## if (length(Ss) > 0) {
  ##   names(Ss) <- 1:length(Ss)
  ## }
  
  ## Create forward and backward pointers
  o <- Outline(P)
  o <- simplify.outline(o)
  
  ## Check that P is more-or-less closed
  if (vecnorm(P[1,] - P[nrow(P),]) > (d.close * diff(range(P[,1])))) {
     plot.map(map, TRUE)
     points(P[c(1,nrow(P)),], col="black")
     stop("Unable to find a closed outline.")
  }

  d <- Dataset(o, dataset, Ds, Ss, cols=cols, raw=list(map=map, sys=sys))
  a <- AnnotatedOutline(d)
  a <- RetinalDataset(a)
  return(a)
}

##' Test the oputline object \code{o} for the prescense of
##' potential optice disc. This is done by checking that the list of
##' landmark lines \code{Ss} exists.
##'
##' @title Test for a potential optic disc
##' @param o Outline object
##' @return \code{TRUE} if an optic disc may be present; \code{FALSE} otherwise
##' @author David Sterratt
retistruct.potential.od <- function(o) {
  if (inherits(o, "dataset")) {
    return(with(o, exists("Ss")) && is.list(o$Ss) && (length(o$Ss) > 0))
  }
  return(FALSE)
}

##' Read the markup data contained in the files \file{markup.csv},
##' \file{P.csv} and \file{T.csv} in the directory \file{dataset},
##' which is specified in the reconstruction object \code{r}. 
##'
##' The tear information is contained in the files \file{P.csv} and
##' \file{T.csv}. The first file contains the locations of outline
##' points that the tears were marked up on. The second file contains
##' the indicies of the apicies and backward and forward verticies of
##' each tear. It is necessary to have the file of points just in case
##' the algorithm that determines \code{P} in
##' \code{\link{retistruct.read.dataset}} has changed since the markup
##' of the tears.
##'
##' The remaining information is contained  containted in the file
##' \file{markup.csv}.
##'
##' If \code{DVflip} is specified, the locations of points \code{P}
##' flipped in the \eqn{y}-direction. This operation also requires the
##' swapping of \code{gf}  and \code{gb} and \code{VF} and \code{VB}.
##' @title Read the markup data
##' @param o Dataset object, containing \code{dataset} path
##' @param error Function to run on error, by default \code{stop()}
##' @return o \code{RetinalDataset} object
##' \item{V0}{Indicies in \code{P} of apicies of tears}
##' \item{VB}{Indicies in \code{P} of backward verticies of tears}
##' \item{VF}{Indicies in \code{P} of backward verticies of tears}
##' \item{iN}{Index in \code{P} of nasal point, or \code{NA} if not marked}
##' \item{iD}{Index in \code{P} of dorsal point, or \code{NA} if not marked}
##' \item{iOD}{Index in \code{Ss} of optic disc }
##' \item{phi0}{Angle of rim in degrees}
##' \item{DVflip}{Boolean variable indicating if DV axis has been flipped}
##' @author David Sterratt
retistruct.read.markup <- function(a, error=stop) {
  ## Return index in P of closest point to x
  closest <- function(P, x) {
    if (any(is.na(x))) {
      return(NA)
    }
    return(which.min(vecnorm(t(t(P) - x))))
  }

  ## Function to map old tears (M.old) and old points (P.old) onto new
  ## points (P). It returns a new matrix of indicies (M).
  convert.markup <- function(M.old, P.old, P) {
    M <- sapply(M.old, function(i) {
      ifelse(is.numeric(i), closest(P, P.old[i,]), i)
    })
    return(M)
  }

  ## Read in the old P data
  Pfile <- file.path(a$dataset, "P.csv")
  if (file.exists(Pfile)) {
    P.old <- as.matrix(read.csv(Pfile))
  } else {
    P.old <- a$P
  }
  
  ## Read in markup file
  markupfile <- file.path(a$dataset, "markup.csv")
  if (file.exists(markupfile)) {
    M.old <- read.csv(markupfile)
    ## Convert to vector
    M <- sapply(M.old, function(x) x)
    if (!is.na(M["iD"])) {
      M["iD"] <- convert.markup(M["iD"], P.old, a$P)
      a <- setFixedPoint(a, M["iD"], "Dorsal")
    }
    if (!is.na(M["iN"])) {
      M["iN"] <- convert.markup(M["iN"], P.old, a$P)
      a <- setFixedPoint(a, M["iN"], "Nasal")
    }
    a$phi0 <- M["phi0"]*pi/180
    if ("iOD" %in% names(M)) {
      a <- nameLandmark(a, M["iOD"], "OD")
    }
    if ("DVflip" %in% names(M)) {
      a$DVflip <- M["DVflip"]
      if ("side" %in% names(M)) {
        a$side <- M["side"]
      }
    }
  } else {
    error("Markup file M.csv doesn't exist.")
  }
  
  ## Read in tearfile
  tearfile <- file.path(a$dataset, "T.csv")
  if (file.exists(tearfile)) {
    T.old <- read.csv(tearfile)
    cn <- colnames(T.old)
    T <- matrix(convert.markup(as.matrix(T.old), P.old, a$P), ncol=3)
    colnames(T) <- cn
    for (i in 1:nrow(T)) {
      a <- addTear(a, T[i,])
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
retistruct.check.markup <- function(o) {
  return(!is.null(names(o$i0)))
}

##' Given an outline object with a \code{dataset} field,  read the
##' reconstruction data from the file \file{\var{dataset}/r.Rdata}.
##'
##' @title Read the reconstruction data from file
##' @param o Outline object containing \code{dataset} field
##' @return If the reconstruction data exists, return a reconstruction
##' object, else return the outline object \code{o}.
##' @author David Sterratt
retistruct.read.recdata <- function(o) {
  recfile <- file.path(o$dataset, "r.Rdata")
  if (file.exists(recfile)) {
    load(recfile)                       # This puts r in the environment
    if (is.null(r$version) || (r$version != recfile.version)) {
      unlink(recfile)
      warning("The algorithm has changed significantly since this retina was last reconstructed, so the cached reconstruction data has been deleted.")
    } else {
      return(r)
    }
  }
  return(NULL)
}

##  Reconstructing the retina
retistruct.reconstruct <- function(o, report=retistruct.report,
                                   plot.3d=FALSE, dev.grid=NA, dev.polar=NA) {
  ## Check that markup is there
  if (!retistruct.check.markup(o)) {
    stop("Neither dorsal nor nasal pole specified")
  }

  ## Check tears are valid
  ct <- checkTears(o)
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
    o$lambda0 <- 0
  }
  if (!is.na(getLandmarkID(o, "OD"))) {
    o$Ds[["OD"]] <- with(o, matrix(colMeans(Ss[[getLandmarkID(o, "OD")]]), 1, 2))
  }

  ## Now do folding itself
  r <- NULL
  r <- ReconstructedOutline(o,
                    report=report,
                    plot.3d=plot.3d, dev.grid=dev.grid,
                    dev.polar=dev.polar)
  if (!is.null(r)) {
    r <- ReconstructedDataset(r, report=report)
    if (!is.na(getLandmarkID(r, "OD"))) {
      r$EOD <- 90 + r$Dss[["OD"]][1,"phi"] * 180/pi
    }
    report(paste("Mapping optimised. Error:", format(r$opt$value,5),
                 ";", r$nflip, "flipped triangles. OD displacement:",
                 format(r$EOD, 2),
                 "degrees."))
    r <- RetinalReconstructedDataset(r)
  }
  return(r)
}

## retistruct.save.markup() - save markup
retistruct.save.markup <- function(a) {
  if (inherits(a, "retinalDataset")) {
    with(a, {
      ## Save the tear information and the outline
      write.csv(cbind(V0, VB, VF), file.path(dataset, "T.csv"),
                row.names=FALSE)
      write.csv(P, file.path(dataset, "P.csv"), row.names=FALSE)
      
      ## Save the dorsal and nasal locations and phi0 to markup.csv
      iD <- ifelse(names(i0) == "Dorsal", i0, NA)
      iN <- ifelse(names(i0) == "Nasal" , i0, NA)
      iOD <- which(names(Ss)=="OD")
      if (length(iOD)==0)
        iOD <- NA
      markup <- data.frame(iD=iD, iN=iN, phi0=phi0*180/pi, iOD=iOD, DVflip=DVflip, side=side)     
      write.csv(markup, file.path(dataset, "markup.csv"), row.names=FALSE)
    })
  }
}

## retistruct.save.recdata() - save reconstruction data
retistruct.save.recdata <- function(r) {
  if (!is.null(r$dataset)) {
    ## Save the derived data
    r$version <- recfile.version        # Datafile version
    if (!is.null(r)) {
      save(r, file=file.path(r$dataset, "r.Rdata"))
    }
  }
}

## retistruct.export.matlab() - save selected reconstruction data to matlab
retistruct.export.matlab <- function(r) {
  if (!is.null(r$dataset)) {
    if (!is.null(r)) {
      f <- file.path(r$dataset, "r.mat")
      message(paste("Saving", f))
      writeMat(f,
               phi0=r$phi0*180/pi,
               Dss=getDss(r),
               DssMeans=getDss.mean(r),
               Sss=name.list(getSss(r)),
               Tss=name.list(getTss(r)),
               side=r$side, DVflip=r$DVflip)
    }
  }
}

