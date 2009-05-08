require(rgl)
require(geometry)

fd <- function(p, ...) sqrt((p^2)%*%c(1,1,1)) - 1
                                        # also predefined as `mesh.dsphere'
fd <- function(p, phi, ...) {
  ## also predefined as `mesh.dsphere'
  d <- sqrt((p^2)%*%c(1,1,1)) - 1
  d <- d + (sqrt((p^2)%*%c(1,1,0)) < cos(phi)) * (p[,3] > 0) *
    (p[,3] - sin(phi) - (sqrt((p^2)%*%c(1,1,1)) - 1))
  return(d)
}

## fd <- function(p, ...) {
##   out <- sqrt((p^2)%*%c(1,1,1)) - 1 + (sqrt((p^2)%*%c(1,1,0)) < 0.5) * (p[,3] > 0)

  ## print((sqrt((p^2)%*%c(1,1,0)) < 0.5) * (p[,3] > 0))
  ## out <- sqrt((p^2)%*%c(1,1,1)) - 1
  ## print(dim(p))
  ## print(dim(out))
##  return(out)
##}
                                        # also predefined as `mesh.dsphere'
fh <- function(p, phi, ...) {
  h <- rep(1,nrow(p))
##   h <- h + 0.5 * (sqrt((p^2)%*%c(1,1,0)) < sin(phi)) * (p[,3] > 0)
  return(h)
}
                                        # also predefined as `mesh.hunif'

tri.area <- function(P, Pt) {
  A <- P[Pt[,1],]
  B <- P[Pt[,2],]
  C <- P[Pt[,3],]
  AB <- cbind(B-A, 0)
  AC <- cbind(C-A, 0)
  return(0.5 * sqrt(apply(extprod3d(AB, AC))^2, 1, sum))
}

## Find the normals of triangles
tri.normals <- function(P, Pt) {
  A <- P[Pt[,1],]
  B <- P[Pt[,2],]
  C <- P[Pt[,3],]

  AB <- cbind(B-A)
  AC <- cbind(C-A)
  print(dim(AB))
  print(dim(AC))  
  return(extprod3d(AB, AC))
}

## Create the hemisphere
phi <- 50 * pi/180
nfix <- round(cos(phi)* 2 * pi / 0.24 - 1)
pfix <- cbind(cos(phi) * sin(2*pi*(1:nfix)/nfix),
              cos(phi) * cos(2*pi*(1:nfix)/nfix),
              sin(phi))
bbox <- matrix(c(-1,1),2,3)
P <- distmeshnd(fd,fh,0.2,bbox, maxiter=100, pfix=pfix, phi=phi)
Pt <- delaunayn(P)

## This is the triagulation of the surface
Ptri = surf.tri(P, Pt)
normals <- tri.normals(P, Ptri)

## areas <- tri.area(P, Pt[,1:3])
ptri <- Ptri[atan2(normals[,3], sqrt(apply(normals[,1:2]^2, 1, sum))) < (pi/2 - 0.001),]
## ptri <- Ptri

##i.rm <- which(p[,3]>0.5)
##pt <- pt[!apply(matrix(pt %in% i.rm, nrow(pt), 4), 1, any),]
##p[i.rm,] <- c(NA, NA, NA)

## Create the asymmetric connectivity matrix
Cu <- rbind(ptri[,1:2], ptri[,2:3], ptri[,c(1,3)])
Cu <- Unique(Cu)

## Find lengths of connections
Ls <- sqrt(apply((P[Cu[,1],] - P[Cu[,2],])^2, 1, sum))

ptri <- t(ptri)
rgl.clear()
rgl.triangles(P[ptri,1], P[ptri,2], P[ptri,3], col="white", front="line")



