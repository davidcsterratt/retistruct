## Return initialised userdata list
retistruct.initialise.userdata <- function() {
  A <<- c()          # Indices of apices of tears
  VB <<- c()         # Indices of forward verticies of tears
  VF <<- c()         # Indices of backward verticies of tears
  phi0 <<- 50        # Height of rim of retina in degrees
  r <<- NULL         # Reconstruction object
  iN <<- NA          # Index of nasal point
  iD <<- NA          # Index of dorsal point
  recfile.version <<- 1   # Version of reconstruction file data format
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
##   A       - tear apices
##   VF      - tear forward verticies
##   VB      - tear backward verticies
retistruct.read.dataset <- function(mess=retistruct.mess) {
  retistruct.initialise.userdata()
  map <<- read.map(dataset)
  sys <<- read.sys(dataset)
  Ds <<- list(green=cbind(na.omit(sys[,'XGREEN']), na.omit(sys[,'YGREEN'])),
              red  =cbind(na.omit(sys[,'XRED'])  , na.omit(sys[,'YRED'])))
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
  P  <<- Ss[[which.max(l)]]
  Ss <<- Ss[-which.max(l)]
  names(Ss) <<- 1:length(Ss)

  ## Create forward and backward pointers
  t <- triangulate.outline(P, n=NA)
  gf <<- t$gf
  gb <<- t$gb

  ## Read in tearfile
  tearfile <- file.path(dataset, "T.csv")
  if (file.exists(tearfile)) {
    T <- read.csv(tearfile)
    A  <<- T[,1]                            # apicies of tears
    VB <<- T[,2]                           # forward verticies
    VF <<- T[,3]                           # backward verticies
  }
  markupfile <- file.path(dataset, "markup.csv")
  if (file.exists(markupfile)) {
    M <<- read.csv(markupfile)
    iD <<- M[1, "iD"]
    iN <<- M[1, "iN"]
    phi0 <<- M[1, "phi0"]
  }

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
  ct <- check.tears(cbind(A, VF, VB), gf, gb, P)
  if (length(ct)) {
    mess(paste("Invalid tears marked up. Somehow tears",
                   toString(ct),
                   "got messed up.  The red lines should always point in the clockwise direction and the orange ones in the anticlockwise direction. To fix, select \"Move Point\", and double click on tear",
                   toString(ct), "."), title="Error", icon="error")

    enable.widgets(TRUE)
    return()
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
  r <<- fold.outline(P, cbind(A, VB, VF), phi0, i0=i0, lambda0=lambda0,
                     Ds=Ds,
                     report=report,
                     plot.3d=plot.3d, dev.grid=dev.grid, dev.polar=dev.polar)
}

## retistruct.save() - save markup and state
retistruct.save <- function() {
  if (!is.null(dataset)) {
    ## Save the tear information and the outline
    write.csv(cbind(A, VB, VF), file.path(dataset, "T.csv"),  row.names=FALSE)
    write.csv(P, file.path(dataset, "P.csv"), row.names=FALSE)

    ## Save the dorsal and nasal locations and phi0 to markup.csv
    markup <- data.frame(iD=iD, iN=iN, phi0=phi0)    
    write.csv(markup, file.path(dataset, "markup.csv"))

    ## Save the derived data
    r$version <- recfile.version        # Datafile version
    if (!is.null(r)) {
      save(r, file=file.path(dataset, "r.Rdata"))
    }
  }
}
