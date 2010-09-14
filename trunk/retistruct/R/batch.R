batch.revision <- function() {
  return(as.integer(gsub("Rev: ", "" ,gsub("\\$", "", "$Rev: 259$"))))
}

##source("spheristruct.R")
##source("common.R")
##source("datafile-utils.R")

## This is the master function to process a retina
## Its only argument is a directory, which contains the following files
## * SYS.SYS
## * ALU.MAP
## * TEARMAT.csv - File containing the information on where cuts and tears are
## It returns various quantities which t
process.retina <- function(dir="data/Anatomy/ALU/M643-4/CONTRA") {
  ## Read in data
  sys <- read.sys(dir)
  map <- as.matrix(read.map(dir))

  ## Corner analysis
  segs <- map.to.segments(map)

  ## Some curation is required here. FIXME - this needs to go
  segs4 <- segs[[4]][23:1,]
  edge.path <- create.path(list(rbind(segs[[1]], segs[[3]], segs4)),
                           close=TRUE) 
  tearmat <- read.csv(paste(dir, "/TEARMAT.csv", sep=""), header=FALSE)
  colnames(tearmat) <- c("apex", "end1", "end2")
  print(tearmat)
  
  ## Actual munging routine
  P <- edge.path[-nrow(edge.path),1:2]

  s <- stitch.retina(P, tearmat)
  plot.stitch(s)

  t <- make.triangulation(s)
  with(t, trimesh(T, P, col="black"))

  m <- merge.points(t, s)
  ## Plot stiched retina in 2D (messy)
  ## trimesh(Tt, Pt, col="black")

  ## Plotting
  plot(P)
  with(s, plot.outline(P, gb))

  p <- project.to.sphere(m, t, phi0=50*pi/180)

  ## Initial plot in 3D space
  plot.retina(p$phi, p$lambda, p$R, m$Tt, m$Rsett)

  ## In order to plot the points, we need to know in which triangle they
  ## lie and where in that triangle they are. The first job is to create
  ## a triangulation, and then we can use tsearch to find the identity of the
  ## triangles and the location in the triangles.
  P.red   <- cbind(na.omit(sys[,"XRED"]), na.omit(sys[, "YRED"]))
  P.green <- cbind(na.omit(sys[,"XGREEN"]), na.omit(sys[, "YGREEN"]))
  P.gridcoo <- cbind(na.omit(sys[,"XGRIDCOO"]), na.omit(sys[, "YGRIDCOO"]))
  cb.red     <- with(t, tsearchn(P, T, P.red))
  cb.green   <- with(t, tsearchn(P, T, P.green))
  cb.gridcoo <- with(t, tsearchn(P, T, P.gridcoo))

  ## Optimisation procedure
  ## r <- optimise.mapping(p, m, t, s, E0.A=100)
  r <- optimise.mapping(p, m, t, s, E0.A=exp(3), k.A=1)
  p1 <- p
  p1$phi <- r$phi
  p1$lambda <- r$lambda
  r <- optimise.mapping(p1, m, t, s, E0.A=exp(10), k.A=20)

  return(list(s=t, t=t, m=m, p=p, r=r, dir=dir, sys=sys, map=map,
              cb.red=cb.red, cb.green=cb.green, cb.gridcoo=cb.gridcoo))
}

## Function to write all this data to file
## This produces three CSV files
## * SCRED.csv -     locations of red cell bodies in spherical
##                   phi-lambda coordinates
## * SCGREEN.csv -   locations of green cell bodies in spherical
##                   phi-lambda coordinates
## * SCGRIDCOO.csv - locations of centres of grid cells within
##                   retina in spherical phi-lambda coordinates, along with
##                   other information about the count of cells within
##                   those squares 
save.projection <- function(dat) {
  with(dat, {
    ## Save the positions of the cell bodies in spherical coordinates
    cs.red <- with(as.list(c(t, m, p)), cell.bodies.folded.sphere(r$phi, r$lambda, R, Tt, cb.red))
    write.csv(cs.red, paste(dir, "/SCRED.csv", sep=""))

    cs.green <- with(as.list(c(t, m, p)), cell.bodies.folded.sphere(r$phi, r$lambda, R, Tt, cb.green))
    write.csv(cs.green, paste(dir, "/SCGREEN.csv", sep=""))
    
    ## Save the positions of the grid coordinates and the counts
    ## in spherical coordinates
    cs.gridcoo <- with(as.list(c(t, m, p)), cell.bodies.folded.sphere(r$phi, r$lambda, R, Tt, cb.gridcoo))

    i <- which(!is.na(sys[, "YGRIDCOO"]))
    dat <- cbind(PHIGRID=cs.gridcoo$phi,
                 LAMBDAGRID=cs.gridcoo$lambda,
                 sys[i, c("XGRIDCOO", "YGRIDCOO", "COMPLETE",
                          "TOTALCEL", "TOTALRED", "TOTALGRE", "TOTALDOU")])
    
    write.csv(dat[!is.na(cs.gridcoo$phi),],
              paste(dir, "/SCGRIDCOO.csv", sep=""))
    
  })
}

