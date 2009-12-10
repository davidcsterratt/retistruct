require("foreign")

## Function to read the file containing the systat file with the
## locations of the cell bodies in it
read.sys <- function(dir=NULL) {
  read.systat(paste(dir, "SYS.SYS", sep="/"))
}

## Function to read the file containing the "map", i.e. the outline of
## the retina
read.map <- function(dir=NULL) {
  read.csv(paste(dir, "ALU.MAP", sep="/"), sep=" ", header=FALSE)
}

## Function to get corners of sample boxes in terms of lower and upper
## coordinates
## This is to be extended to draw boxes around the edge of the map
get.sys.boxes <- function(sys) {
  ci <- which(sys[,"COMPLETE"] == 1)
  lower <- rbind(sys[ci, "XGRIDCOO"]-sys[ci, "BOXSIZEX"],
                 sys[ci, "YGRIDCOO"]-sys[ci, "BOXSIZEY"])
  upper <- lower + rbind(2*sys[ci, "BOXSIZEX"],
                         2*sys[ci, "BOXSIZEY"])
  
  return(list(lower = lower, upper = upper))
}

## This function gets the lower and upper coordinates of boxes
## around the edge of the map, with width fw
get.map.boxes <- function(map, fw=1000) {
  ## Find the edges of the map
  mux <- max(map[,1])
  muy <- max(map[,2])
  mlx <- min(map[,1])
  mly <- min(map[,2])

  ## We need four boxes arround the edge
  lower <- rbind(c(mlx-fw, mlx   , mux   , mlx ),
                 c(mly-fw, muy   , mly-fw, mly-fw))
  upper <- rbind(c(mlx   , mux   , mux+fw, mux),
                 c(muy+fw, muy+fw, muy+fw, mly))
  return(list(lower = lower, upper = upper))
}

## This function gets the lower and upper coordinates of the "sys" and
## "map" boxes
get.boxes <- function(sys, map) {
  boxes.sys <- get.sys.boxes(sys)
  boxes.map <- get.map.boxes(map)
  return(list(lower=cbind(boxes.sys$lower, boxes.map$lower),
              upper=cbind(boxes.sys$upper, boxes.map$upper)))
}

## Convert map file to a list of segments
map.to.segments <- function(map) {
  colnames(map) <- c("X", "Y")
  segs <- list()
  i <- 1
  while(i < dim(map)[1]) {
    n <- map[i, 2]
    inds <- (i+1):(i+n)
    segs <- c(segs, list(map[inds,]))
    i <- i + n + 1
  }
  return(segs)
}

## Function to plot the "map", i.e. the outline of the retina
plot.map <- function(map, seginfo=FALSE) {
  par(mar=c(2,2,1.5,0.5))
  plot(NA, NA, xlim=lim(map[,1]), ylim=lim(map[,2]), xlab="", ylab="")
  ## Read first line of map; this tells us how many points
  ## are going to be plotted in this stroke
  segs <- map.to.segments(map)
  for (i in 1:length(segs)) {
    seg <- segs[[i]]
    col <- "black"
    if (seginfo) {
      col <- rainbow(10)[i]
      print(col)
      text(seg[1,1], seg[1,2], labels=i, col=col)
    }
    lines(seg[, 1], seg[, 2], lwd=2, col=col)    
  }
}

## Function to plot a box specified by a vector at the bottom left
## corner (lower) and a vector at the top right corner (upper)
plot.box <- function(lower, upper) {
  lines(c(lower[1], upper[1], upper[1], lower[1], lower[1]),
        c(lower[2], lower[2], upper[2], upper[2], lower[2]))
}

## Function to plot boxes specified by the corners specified by the
## columns of lower and upper (see plot.box())
plot.boxes <- function(lower, upper) {
  for(i in 1:(dim(lower)[2])) {
    plot.box(lower[,i], upper[,i])
  }
}

# Function to plot the unfolded retina, grid and labelled RG cells
plot.sys.map <- function(sys, map) {
  plot.map(map)
  points.sys(sys)
  
  boxes <- get.sys.boxes(sys) 
  plot.boxes(boxes$lower, boxes$upper)
}  

# Function to plot the labelled RG cells
points.sys <- function(sys) {
  points(sys[,'XGREEN'], sys[,'YGREEN'], col="green", pch=20,cex=0.5)
  points(sys[,'XRED'], sys[,'YRED'], col="red", pch=20,cex=0.5)
}

