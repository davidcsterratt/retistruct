retistruct.global.revision <- 507

## Set initialised userdata in environment
retistruct.initialise.userdata <- function() {
  V0 <<- c()          # Indices of apices of tears
  VB <<- c()         # Indices of forward verticies of tears
  VF <<- c()         # Indices of backward verticies of tears
  phi0 <<- 50        # Height of rim of retina in degrees
  r <<- NULL         # Reconstruction object
  iN <<- NA          # Index of nasal point
  iD <<- NA          # Index of dorsal point
  iOD <<- NA         # Index of segment which is Optic Disc
}

## Set constants in environment
retistruct.initialise.consts <- function() {
  recfile.version <<- 2   # Version of reconstruction file data format
}

## Test message function, with similar arguments to gmessage
retistruct.mess <- function(message, title="",...) {
  cat(paste(title, ":", message, sep=""))
}

## Report function, with similar arguments to print
retistruct.report <- function(message, title="",...) {
  cat(paste(message, "\n", sep=""))
}

## Function to read a dataset
## 
## Changes the following global variables:
##   dataset - directory in which data is contained
##   map     - the map data
##   Ds      - list of datapoints
##   Ss      - list of landmark lines
##   P       - the outline
##   gf      - forward pointers
##   gb      - backward pointers
##   V0      - tear apices
##   VF      - tear forward verticies
##   VB      - tear backward verticies
##
## The function can throw various errors
##
retistruct.read.dataset <- function(mess=retistruct.mess, d.close=1000) {
  ## Initialise global variables
  retistruct.initialise.userdata()

  ## Read the raw data
  map <<- read.map(dataset)
  sys <<- read.sys(dataset)

  ## Extract datapoints
  ##
  ## At present, for the plotting functions to work, the name of each
  ## group has to be a valid colour.
  Ds <<- list(green =cbind(na.omit(sys[,'XGREEN']), na.omit(sys[,'YGREEN'])),
              red   =cbind(
                na.omit(c(sys[,'XDOUBLE'], sys[,'XRED'])),
                na.omit(c(sys[,'YDOUBLE'], sys[,'YRED']))),
              double=cbind(na.omit(sys[,'XDOUBLE']),na.omit(sys[,'YDOUBLE'])))
  D.cols <<- list(green="green",
                  red="red",
                  double="yellow",
                  od="blue")
    
  ## Extract line data
  segs <- map.to.segments(map)

  ## Connect together segments that look to be joined
  Ss <<- connect.segments(segs)

  ## Determine the lengths of the segments
  l <- c()
  for (i in 1:length(Ss)) {
    l[i] <- segment.length(Ss[[i]])
  }

  ## The outline (P) is the longest connected segment and the outline
  ## is removed from the list of segments
  P  <- Ss[[which.max(l)]]
  Ss <<- Ss[-which.max(l)]
  if (length(Ss) > 0) {
    names(Ss) <<- 1:length(Ss)
  }
  
  ## Create forward and backward pointers
  t <- triangulate.outline(P, n=NA)
  P <<- t$P
  gf <<- t$gf
  gb <<- t$gb

  
  ## Check that P is more-or-less closed
  if (vecnorm(P[1,] - P[nrow(P),]) > d.close) {
    plot.map(map, TRUE)
    points(P[c(1,nrow(P)),], col="black")
    stop("Unable to find a closed outline.")
  }
}

## Return logical value if optic disc may be present
retistruct.potential.od <- function() {
  return(exists("Ss") && is.list(Ss) && (length(Ss) > 0))
}

## Return index in P of closest point to x
retistruct.closest <- function(P, x) {
  if (any(is.na(x))) {
    return(NA)
  }
  return(which.min(vecnorm(t(t(P) - x))))
}

retistruct.convert.markup <- function(M.old, P.old, P) {
  M.old <- as.matrix(M.old)
  M <- matrix(sapply(M.old, function(i) {retistruct.closest(P, P.old[i,])}),
                ncol=ncol(M.old))
  colnames(M) <- colnames(M.old)
  return(M)
}

## retistruct.read.markup - read the markup data, if it exists
##
## Relies on global variable dataset
## 
retistruct.read.markup <- function(mess=retistruct.mess) {
  ## Read in the old P data
  Pfile <- file.path(dataset, "P.csv")
  if (file.exists(Pfile)) {
    P.old <- as.matrix(read.csv(Pfile))
  } else {
    P.old <- P
  }
  
  ## Read in tearfile
  tearfile <- file.path(dataset, "T.csv")
  if (file.exists(tearfile)) {
    T.old <- read.csv(tearfile)
    T <- retistruct.convert.markup(T.old, P.old, P)
    V0 <<- T[,1]                            # apicies of tears
    VB <<- T[,2]                           # forward verticies
    VF <<- T[,3]                           # backward verticies
  } else {
    stop("Tear file T.csv doesn't exist.")
  }
  markupfile <- file.path(dataset, "markup.csv")
  if (file.exists(markupfile)) {
    M.old <<- read.csv(markupfile)
    M <- retistruct.convert.markup(M.old, P.old, P)
    iD <<- M[1, "iD"]
    iN <<- M[1, "iN"]
    phi0 <<- M[1, "phi0"]
    if ("iOD" %in% colnames(M)) {
      iOD <<- M[1, "iOD"]
    }
  } else {
    stop("Markup file M.csv doesn't exist.")
  }
}

##' Check that all markup is present
##'
##' Check that markup such as tears and the nasal or dorsal points are present
##' @title Retistruct check markup
##' @param  
##' @return If all markup is present, return TRUE. Otherwise return an error string 
##' @author David Sterratt
retistruct.check.markup <- function() {
  if (is.na(iD) && is.na(iN)) {
    return(FALSE)
  }
  return(TRUE)
}

## retistruct.read.recdata - read the reconstruction data, if it exists
##
## Relies on global variable dataset
##
retistruct.read.recdata <- function(mess=retistruct.mess) {
  retistruct.initialise.consts()
  ## The format in f.Rdata is deprecated
  foldfile <- file.path(dataset, "f.Rdata")
  if (file.exists(foldfile)) {
    file.remove(foldfile)
    mess("The algorithm has changed significantly since this retina was last reconstructed, so the cached reconstruction data has been deleted.",
             title="Warning", icon="warning")
  }
  recfile <- file.path(dataset, "r.Rdata")
  if (file.exists(recfile)) {
    load(recfile, globalenv())
    if (r$version != recfile.version) {
      mess("The algorithm has changed significantly since this retina was last reconstructed, so the cached reconstruction data has been deleted.",
               title="Warning", icon="warning")
      r <<- NULL
    }
  } else {
    r <<- NULL
  }
}

##  Reconstructing the retina
retistruct.reconstruct <- function(mess=retistruct.mess,
                                   report=retistruct.report,
                                   plot.3d=FALSE, dev.grid=NA, dev.polar=NA) {
  if (!retistruct.check.markup()) {
    stop("Neither dorsal nor nasal pole specified")
  }
  ct <- check.tears(cbind(V0, VF, VB), gf, gb, P)
  if (length(ct)) {
    stop(paste("Invalid tears", toString(ct), "marked up. Fix using \"Move Point\"."))
  }
  i0 <- 0
  lambda0 <<- 0
  if (!is.na(iD)) {
    i0 <- iD
    lambda0 <<- 90
  }
  if (!is.na(iN)) {
    i0 <- iN
    lambda0 <<- 0
  }
  if (!is.na(iOD)) {
    Ds[["od"]] <<- matrix(colMeans(Ss[[iOD]]), 1, 2)
  }
  
  r <<- NULL
  r <<- fold.outline(P, V0, VB, VF, phi0, i0=i0, lambda0=lambda0,
                     report=report,
                     plot.3d=plot.3d, dev.grid=dev.grid, dev.polar=dev.polar)
  if (!is.null(r)) {
    r <<- infer.datapoint.landmark.coordinates(r, Ds=Ds, Ss=Ss,
                                               report=report)
    r <<- merge.lists(r, list(D.cols=D.cols))
    if (!is.na(iOD)) {
      r$EOD <<- 90 + r$Dss[["od"]][1,"phi"] * 180/pi
    }
    report(paste("Mapping optimised. Error:", format(r$opt$value,5),
                 ";", r$nflip, "flipped triangles. OD displacement:",
                 format(r$EOD, 2),
                 "degrees."))
  }
}

## retistruct.save.markup() - save markup
retistruct.save.markup <- function() {
  if (!is.null(dataset)) {
    ## Save the tear information and the outline
    write.csv(cbind(V0, VB, VF), file.path(dataset, "T.csv"),  row.names=FALSE)
    write.csv(P, file.path(dataset, "P.csv"), row.names=FALSE)

    ## Save the dorsal and nasal locations and phi0 to markup.csv
    markup <- data.frame(iD=iD, iN=iN, phi0=phi0, iOD=iOD)    
    write.csv(markup, file.path(dataset, "markup.csv"), row.names=FALSE)
  }
}

## retistruct.save.recdata() - save reconstruction data
retistruct.save.recdata <- function() {
  retistruct.initialise.consts()
  if (!is.null(dataset)) {
    ## Save the derived data
    r$version <- recfile.version        # Datafile version
    if (!is.null(r)) {
      save(r, file=file.path(dataset, "r.Rdata"))
    }
  }
}

## retistruct.export.matlab() - save selected reconstruction data to matlab
retistruct.export.matlab <- function() {
  if (!is.null(dataset)) {
    if (!is.null(r)) {
      r <<- infer.tear.coordinates(r)
      f <- file.path(dataset, "r.mat")
      print(paste("Saving", f))
      writeMat(f, phi0=r$phi0, Dss=r$Dss, Sss=name.list(r$Sss), Tss=name.list(r$Tss))

      save(r, file=file.path(dataset, "r.Rdata"))
    }
  }
}

