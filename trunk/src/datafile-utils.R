require("foreign")
source("geometry.R")
source("misc.R")
library(geometry)

## Function to return minimum and maximum values of vector x
lim <- function(x) {
  return(c(min(x), max(x)))
}

## Function to read the file containing the systat file with the
## locations of the cell bodies in it
read.sys <- function(dir=NULL) {
  read.systat(paste(dir, "SYS.SYS", sep="/"))
}

## SYS.MAP might be better to use
## Function to read the file containing the "map", i.e. the outline of
## the retina
read.map <- function(dir=NULL) {
  map <- read.csv(paste(dir, "ALU.MAP", sep="/"), sep=" ", header=FALSE)
  return(as.matrix(map))
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

## Convert map matrix to a list of segments
map.to.segments <- function(map) {
  ## Give named columns of the map, and hence segment list
  colnames(map) <- c("X", "Y")
  segs <- list()                        # Initialise list of segments
  i <- 1                                # Current line of map
  while(i < dim(map)[1]) {
    ## Read current line of map; this tells us how many points
    ## belong to this stroke
    n <- map[i, 2]
    ## Add the points to the list of segments
    inds <- (i+1):(i+n)
    segs <- c(segs, list(map[inds,]))
    ## Update the current position
    i <- i + n + 1
  }
  return(segs)
}

## Connect segments whose ends lie close to each other
## Input arguments:
## segs - list of segments to connect
## Ouput:
## Ts   - list of connected segments
connect.segments <- function(segs) {
  N <- length(segs)                        # Number of segments
  ## First find the first and last points of each segment
  ## The start points are in columns 1:N of P, and the end points in
  ## (N+1):(2*N)
  P <- as.matrix((sapply(segs, function(S) {S[1,]} )))
  P <- cbind(P,
             as.matrix((sapply(segs, function(S) {S[nrow(S),]} ))))

  ## Now create vector of neighbours of points n
  ## Each element is the index of the closest point
  ## (apart from the point itself) 
  n <- c()                              # Initialisation
  for (i in 1:(2*N)) {
    d <- norm(t(P[,i] - P))             # Distances of P_i from all other points
    d[i] <- NA                          # Ignore distance to self
    n[i] <- which.min(d)                # Find index of closest point
  }

  ## Now create matrix T in which to store the sorted points,
  ## one per row
  i <- 1                                # Initial index
  Ts <- list()                          # New segment list
  T <- matrix(0, 0, 2)                  # Initialisation
  while(1) {
    j <- Mod(i + N, 2*N)                # Index at end of segment
    k <- n[j]                           # Closest index in another segment
    
    n[i] <- NA                   # Ignore these indicies in the future
    n[j] <- NA                   # Ignore these indicies in the future

    if (i<=N) {
      T <- rbind(T, segs[[i]])
      print(paste("Adding old segment", i))
    } else {
      T <- rbind(T, flipud(segs[[j]]))
      print(paste("Adding old segment", j))
    }
    ## If the segment connects back to a previously included point,
    ## store the segment and find a new starting index, if one exists
    if (is.na(k)) {
      Ts <- c(Ts, list(unique(T)))
      print(paste("Storing new segment", length(Ts)))
      T <- matrix(0, 0, 2)
      rem <- which(!is.na(n))
      ## Find a new starting index
      if (length(rem) > 0) {
        i <- n[rem[1]]
      } else {
        break
      }
    } else {
      ## Otherwise, continue
      i <- k
    }
  }
  return(Ts)
}

## Return the longest segment, which is hopefully the outline
## Input arguments:
## segs - list of segments to connect
## Ouput:
## P    - matrix containing the outline
segments.to.outline <- function(segs) {
  ## Connect together the segments
  Ts <- connect.segments(segs)

  ## Choose the longest segment
  l <- c()
  for (i in 1:length(Ts)) {
    print(i)
    seg <- Ts[[i]]
    v <- diff(rbind(seg, seg[1,]))
    print(dim(v))
    l[i] <- sum(sqrt(apply(v^2, 1, sum)))
  }
  print(l)
  return(Ts[[which.max(l)]])
}

## Function to plot the "map", i.e. the outline of the retina
plot.map <- function(map, seginfo=FALSE) {
  par(mar=c(2,2,1.5,0.5))
  plot(NA, NA, xlim=lim(map[,1]), ylim=lim(map[,2]), xlab="", ylab="")
  segs <- map.to.segments(map)
  for (i in 1:length(segs)) {
    seg <- segs[[i]]
    col <- "black"
    if (seginfo) {
      col <- rainbow(10)[i]
      print(col)
      text(seg[1,1], seg[1,2] + 100, labels=i, col=col)
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


