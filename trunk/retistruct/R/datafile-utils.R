## Function to check the wether the directory in question is a data
## directory or not. Returns TRUE if so, FALSE otherwise, and throws an
## error if the directory appears to be corrupt.
check.datadir <- function(dir=NULL) {
  sys <- file.exists(file.path(dir, "SYS.SYS"))
  map <- file.exists(file.path(dir, "ALU.MAP"))
  if (sys  & map ) return(TRUE)
  if (!sys & !map) return(FALSE)
  if (!sys) stop("SYS.SYS file doesn't exist")
  if (!map) stop("ALU.MAP file doesn't exist")
}

## Function to read the file containing the systat file with the
## locations of the cell bodies in it
read.sys <- function(dir=NULL) {
  read.systat(file.path(dir, "SYS.SYS"))
}

## SYS.MAP might be better to use
## Function to read the file containing the "map", i.e. the outline of
## the retina
read.map <- function(dir=NULL) {
  e <- function() {stop("Corrupt MAP file.")}
  map <- tryCatch(read.csv(file.path(dir, "ALU.MAP"), sep=" ", header=FALSE),
                  warning=e, error=e)
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
## segs -         list of segments to connect
## merge.rad -    Radius within which to merge start and end points of segments
##                NOT IMPLEMENTED YET
## Ouput:
## Ts   -         list of connected segments
connect.segments <- function(segs, merge.rad=10) {
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
  d <- as.matrix(dist(t(P), diag=TRUE))
  diag(d) <- NA
  n <- apply(d, 1, which.min)           # Find index of closest point

  ## Create a list Ts, to contain new segments
  Ts <- list()                          # New segment list
  added <- rep(FALSE, 2*N)              # Which points have been used
  while(any(!added)) {
    ## Find points remaning to be added, and select the first as the
    ## initial segment index
    rem <- which(!added)
    i <- rem[1]
    ## Now create matrix T in which to store the sorted points,
    ## one per row
    T <- matrix(0, 0, 2)

    ## If segment with i at one end hasn't been added, add it
    while(!added[i]) {
      j <- mod1(i + N, 2*N)       # Index into P of other end of segment i
      ## Add the segement to T
      if (i<=N) {
        T <- rbind(T, segs[[i]])
        print(paste("Adding old segment", i))
      } else {
        T <- rbind(T, flipud(segs[[j]]))
        print(paste("Adding old segment", j))
      }
      
      added[i] <- TRUE             # Ignore these indicies in the future
      added[j] <- TRUE             # Ignore these indicies in the future
      i <- n[j]                    # Closest index in another segment
    }
      
    ## If the segment connects back to a previously included point,
    ## store the segment and find a new starting index, if one exists
    T <- remove.identical.consecutive.rows(T)
    T <- remove.intersections(T)
    T <- remove.backtracks(T)
    Ts <- c(Ts, list(T))
    print(paste("Storing new segment", length(Ts)))
  }
  return(Ts)
}

## Return the length of a segment
## Input arguments:
## S    - segment to measure
## Ouput:
## l    - length of segment
segment.length <- function(S) {
  v <- diff(rbind(S, S[1,]))
  return(sum(sqrt(apply(v^2, 1, sum))))
}

## Function to plot the "map", i.e. the outline of the retina
plot.map <- function(map, seginfo=FALSE,
                     xlim=range(map[,1]), ylim=range(map[,2])) {
  par(mar=c(2,2,1.5,0.5))
  plot(NA, NA, xlim=xlim, ylim=ylim, xlab="", ylab="")
  segs <- map.to.segments(map)
  for (i in 1:length(segs)) {
    seg <- segs[[i]]
    col <- "black"
    if (seginfo) {
      col <- rainbow(10)[mod1(i,10)]
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

## Remove backtracks to the same point in a path
remove.backtracks <- function(P) {
  N <- nrow(P)
  ## Loop through P
  for (i in 2:N) {
    ## Find if we have already been here
    j <- which(apply(P[1:(i-1),,drop=FALSE], 1, function(x) {identical(x, P[i,])}))
    ## If so, remove all of P between there and here
    if (length(j)) {
      print(paste("Already been to ", j, "; Removing:"))
      ind <- j:(i-1)
      print(ind)
      return(remove.backtracks(P[-ind,]))
    }
  }
  return(P)
}


## Convert segment to pointers
## FIXME: At present this is not used. However, it might form part of a
## fix for Issue #178: error in segments2pointers()
segment.to.pointers <- function(P) {
  ## uniquify P
  N <- nrow(P)
  U <- unique(P)
  id <- 1:N
  for (i in 1:N) {
    id[i] <- which(apply(U, 1, function(x) {identical(x, P[i,])}))
  }
  M <- nrow(U)
  gf <- rep(NA, M)
  for (i in 1:(N-1)) {
    gf[id[i]] <- id[i+1]
  }
  gb <- rep(NA, M)
  for (i in N:2) {
    gb[id[i]] <- id[i-1]
  }
  return(list(id=id, P=U, gf=gf, gb=gb))
}


