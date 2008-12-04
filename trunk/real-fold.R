source("../../data/Anatomy/cluster-analysis.R")

map <- as.matrix(read.map("../../data/Anatomy/ALU/M643-4/CONTRA"))
plot.map(map)

map.to.segments <- function(map) {
  segs <- list()
  i <- 1
  while(i < dim(map)[1]) {
    n <- map[i, 2]
    inds <- (i+1):(i+n)
    segs <- c(segs, list(map[inds,]))
    lines(map[inds, 1], map[inds, 2], lwd=2)
    i <- i + n + 1
  }
  return(segs)
}

join.segments <- function(segs) {
  map <- matrix(0, 0, 2)
  for (i in 1:length(segs)) { 
    map <- rbind(map, segs[[i]])
  }
  return(map) 
}

norm2 <- function(v) {
  return(sqrt(as.vector(t(v) %*% v)))
}  

## Find angle between two vectors
find.angle <- function(v1, v2) {
  sp <- (t(v1) %*% v2) / (norm2(v1) * norm2(v2))
  print(sp)
  return(acos(sp))
}

find.corners <- function(map) {
  angle <- rep(NA, nrow(map))
  for (i in 2:(nrow(map)-1)) {
    v1 <- map[i,] - map[i-1,]
    print(v1)
    v2 <- map[i+1,] - map[i,]
    print(v2)
    angle[i] <-  find.angle(v1, v2) 
  }
  return(angle)
}

plot.triangle <- function(map, inds) {
  i <- inds[c(1:3,1)]
  lines(map[i,1], map[i,2])
}

## Check whether line from x1 to y1 and from x2 to y2 intersect
check.intersection <- function(x1, y1, x2, y2) {
  m <- cbind(y1-x1, -y2+x2)
  if (!(det(m) == 0)) {
    lambda <- solve(m) %*% (x2-x1)
    if (all((lambda<1) & (lambda>0))) {
      return(lambda)
    }
  }
  return(FALSE)
}

segs <- map.to.segments(map)
## map <- join.segments(list(segs[[1]], segs[[3]], segs[[4]]))
##map <- segs[[1]]
palette(rainbow(180))
angles <- find.corners(map) * 180/pi

##text(map[,1], map[,2], labels=format(angles,digits=2))
## points(map[,1], map[,2], col=angles, pch=20)


## Make a grid of indicies
## Specify maxium elevation and grid spacing in degrees
theta.max.deg <- 140 
dtheta.deg <- 10 
dphi.deg   <- 10 
r <- 10000                                # radius

## Make some indicies
is <- 1:(360/dphi.deg)
js <- 1:(theta.max.deg/dtheta.deg)

dphi <-   dphi.deg * pi/180
dtheta <- dtheta.deg * pi/180
theta.max <- theta.max.deg * pi/180

## The matricies X and Y specify the mapping from point i,j in polar space
## to cartesian coordinates
X <- matrix(NA, length(is), length(js))
Y <- matrix(NA, length(is), length(js))

gmap <- expand.grid(i=is, j=js)
gmap <- cbind(gmap,
              phi=  gmap[,"i"]*dphi,
              theta=gmap[,"j"]*dtheta)
gmap <- cbind(gmap, X=gmap[,"theta"]*cos(gmap[,"phi"]))
gmap <- cbind(gmap, Y=gmap[,"theta"]*sin(gmap[,"phi"]))

plot.mesh <- function(gmap) {
  ## Plot the tangential lines
  for (j in js) {
    inds <- (gmap[,"j"] == j)
    polygon(gmap[inds,"X"], gmap[inds,"Y"])
  }
  ## Plot the radial lines
  for (i in is) {
    inds <- (gmap[,"i"] == i)
    polygon(gmap[inds,"X"], gmap[inds,"Y"])
  }
}

## Plot the grid
## plot(NA,NA,xlim=c(-pi, pi), ylim=c(-pi, pi))
plot.mesh(gmap)

## Find the neigbours and put them in columns
## - nl - tangential left (in direction of increasing phi)
## - nr - tangential right (in direction of decreasing phi)
## - nu - radial up (in direction of increasing theta)
## - nd - radial down (in direction of decreasing theta)
M <- length(is)
N <- length(js)
gmap[,"nl"] <- (1:nrow(gmap))+ 1 - M * (((1:nrow(gmap)) %% M) == 0)
gmap[,"nr"] <- (1:nrow(gmap))- 1 + M * (((1:nrow(gmap)) %% M) == 1)
gmap[,"nu"] <- (1:nrow(gmap))+ M
gmap[gmap[,"nu"]>nrow(gmap),"nu"] <- NA
gmap[,"nd"] <- (1:nrow(gmap))- M
gmap[gmap[,"nd"]<1,"nd"] <- NA



## Computes the distance between grid points if they were
## spread over the surface of a sphere
## Find distance to left neigbours
gmap[,"Ll"] <- r * dphi * sin(gmap[,"theta"])
gmap[,"Lr"] <- r * dphi * sin(gmap[,"theta"])
gmap[,"Lu"] <- r * dtheta
gmap[is.na(gmap[,"nu"]),"Lu"] <- NA
gmap[,"Ld"] <- r * dtheta
gmap[is.na(gmap[,"nd"]),"Ld"] <- NA

## Convert to proximity matricies A and B
## Points 1:((N-1)*M) are variable points
## Points ((N-1)*M+1):(N*M) are fixed
C <- matrix(0, (N-1)*M, N*M)
mgmap <- as.matrix(gmap)
for (k in 1:(((N-1)*M))) {
  C[k,na.omit(mgmap[k,c("nl","nr","nu","nd")])] <- 1/na.omit(mgmap[k,c("Ll","Lr","Lu","Ld")])
}
A <- C[,1:((N-1)*M)]
B <- C[,((N-1)*M+1):(N*M)]
## Diagonal matrix of row sums of [ A B ]
D <- diag(apply(cbind(C), 1, sum))

## Initial test with P lying on circle
P <-  mgmap[((N-1)*M+1):(N*M),c("X","Y")]

## The matrix Q is an n by 2 matrix comprising the row vectors of the
## solution points
Q <- solve(D - A) %*% B %*% P

mgmap[,c("X","Y")] <- rbind(Q, P)

plot(NA,NA,xlim=c(-pi, pi), ylim=c(-pi, pi))
plot.mesh(mgmap)
x11()
## Compute force on each element
compute.force <- function(gmap) {
  ## Compute the force due to a component
  compute.force.component <- function(X1, Y1, n, l, L) {
  X2 <- gmap[gmap[,n],"X"]
  Y2 <- gmap[gmap[,n],"Y"]
  gmap[,l] <<- sqrt((X2 - X1)^2 + (Y2-Y1)^2)
  fac <- (gmap[,l] - gmap[,L]) / (gmap[,l] * gmap[,L])
  F.X <- fac * (X2 - X1) 
  F.Y <- fac * (Y2 - Y1)
  F.X[is.na(F.X)] <- 0
  F.Y[is.na(F.Y)] <- 0
  return(cbind(F.X, F.Y))
}
  X1 <- gmap[,"X"]
  Y1 <- gmap[,"Y"]

  Fl <- compute.force.component(X1, Y1, "nl", "ll", "Ll")
  Fr <- compute.force.component(X1, Y1, "nr", "lr", "Lr")
  Fu <- compute.force.component(X1, Y1, "nu", "lu", "Lu")
  Fd <- compute.force.component(X1, Y1, "nd", "ld", "Ld")

##  print(Fl)
  Ftot <- Fl + Fr + Fu + Fd
  
  gmap[,"F.X"] <- Ftot[,1]
  gmap[,"F.Y"] <- Ftot[,2]
  return(gmap)
}

## Compute energy
compute.energy <- function(X1, Y1, gmap) {
  ## Compute the energy due to a component
  compute.energy.component <- function(X1, Y1, n, l, L) {
    X2 <- X1[gmap[,n]]
    Y2 <- Y1[gmap[,n]]
    l <- sqrt((X2 - X1)^2 + (Y2 - Y1)^2)
    E <- (l - gmap[,L])^2 / (2 * gmap[,L])
    E[is.na(E)] <- 0
    return(E)
  }

  El <- compute.energy.component(X1, Y1, "nl", "ll", "Ll")
  Er <- compute.energy.component(X1, Y1, "nr", "lr", "Lr")
  Eu <- compute.energy.component(X1, Y1, "nu", "lu", "Lu")
  Ed <- compute.energy.component(X1, Y1, "nd", "ld", "Ld")

  Etot <- sum(El + Er + Eu + Ed)
  print(Etot)
  return(Etot)
}

compute.energy.gr <- function(X1, Y1, gmap) {
  compute.energy.gr.component <- function(X1, Y1, n, L) {
    X2 <- X1[gmap[,n]]
    Y2 <- Y1[gmap[,n]]
    l <- sqrt((X2 - X1)^2 + (Y2 - Y1)^2)
    L <- gmap[,L]
    dE <-   (l - L)/ (l * L) * c(X1 - X2, Y1 - Y2)
    dE[is.na(dE)] <- 0
    return(dE)
  }
  dEl <- compute.energy.gr.component(X1, Y1, "nl", "Ll")
  dEr <- compute.energy.gr.component(X1, Y1, "nr", "Lr")
  dEu <- compute.energy.gr.component(X1, Y1, "nu", "Lu")
  dEd <- compute.energy.gr.component(X1, Y1, "nd", "Ld")

  dE <- dEl + dEr + dEu + dEd
  return(dE)
}

compute.total.energy <- function(p, gmap) {
  X <- p[1:nrow(gmap)]
  Y <- p[(1:nrow(gmap))+nrow(gmap)]
  return(compute.energy(X, Y, gmap))
}

compute.total.energy.gr <- function(p, gmap) {
  X <- p[1:nrow(gmap)]
  Y <- p[(1:nrow(gmap))+nrow(gmap)]
  return(compute.energy.gr(X, Y, gmap))
}



dt <- 0.001
opt <- list(par = c(gmap[,"X"],gmap[,"Y"]) * r)
if (FALSE) {
for (time in 1:1000) {
  print(time)
##  plot.map(map)
  plot(NA,NA,xlim=r*c(-pi, pi), ylim=r*c(-pi, pi))
  plot.mesh(gmap)
  opt <- optim(opt$par,
               method="CG",
               fn=compute.total.energy,
               gr=compute.total.energy.gr,
               control=list(maxit=10), gmap=gmap)
  gmap[,"X"] <- opt$par[1:nrow(gmap)]
  gmap[,"Y"] <- opt$par[(1:nrow(gmap))+nrow(gmap)]
  ##  gmap <- compute.force(gmap)
##  gmap[,"X"] <- gmap[,"X"] + dt * gmap[,"F.X"]
##  gmap[,"Y"] <- gmap[,"Y"] + dt * gmap[,"F.Y"]
}
}
