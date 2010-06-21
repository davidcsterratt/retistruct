require(rgl)
require(geometry)

## Find the normals of triangles
tri.normals <- function(P, Pt) {
  A <- P[Pt[,1],]
  B <- P[Pt[,2],]
  C <- P[Pt[,3],]
  AB <- cbind(B-A)
  AC <- cbind(C-A)
  return(extprod3d(AB, AC))
}

## Find the areas of triangles
tri.areas <- function(P, Pt) {
  return(0.5 * sqrt(apply(tri.normals(P, Pt))^2, 1, sum))
}

## Create a hemisphere cut off at longitude phi
create.hemisphere <- function(phi=40*pi/180) {
  ## Distance function for sphere which is flattened at a lattitude phi
  ## also predefined as `mesh.dsphere'
  fd <- function(p, phi, ...) {
    d <- sqrt((p^2)%*%c(1,1,1)) - 1
    d <- d + (sqrt((p^2)%*%c(1,1,0)) < cos(phi)) * (p[,3] > 0) *
      (p[,3] - sin(phi) - (sqrt((p^2)%*%c(1,1,1)) - 1))
    return(d)
  }

  ## Uniform length function
  ## also predefined as `mesh.hunif'
  fh <- function(p, phi, ...) {
    h <- rep(1,nrow(p))
    return(h)
  }

  ## A number of fixed points are chosen around the rim, 
  ## which is the line where the flattened part joins
  ## the sphere proper. The number of fixed points is chosen so
  ## that the intervals along the rim are approximately the same
  ## as the internode spacing
  ## nfix <- round(cos(phi)* 2 * pi / 0.24 - 1)
  nfix <- 100
  pfix <- cbind(cos(phi) * sin(2*pi*(1:nfix)/nfix),
                cos(phi) * cos(2*pi*(1:nfix)/nfix),
                sin(phi))
  bbox <- matrix(c(-1,1),2,3)

  ## This fuction from the geometry package creates the mesh, seeking
  ## to make the links as equal in length as possible
  P <- distmeshnd(fd,fh,0.2,bbox, maxiter=100, pfix=pfix, phi=phi)
  Pt <- delaunayn(P)

  ## This is the triagulation of the surface
  Ptri = surf.tri(P, Pt)
  normals <- tri.normals(P, Ptri)

  ## Remove upward-facing triangles as they are on the flattened part
  angles <- abs(atan2(normals[,3], sqrt(apply(normals[,1:2]^2, 1, sum))))
  ptri <- Ptri[angles < (pi/2 - 0.001),]

  ## Create the asymmetric connectivity matrix
  Cu <- rbind(ptri[,1:2], ptri[,2:3], ptri[,c(1,3)])
  Cu <- Unique(Cu, rows.are.sets=TRUE)

  ## Take out points that are no longer in the mesh
  k <- unique(sort(Cu))  # Indicies to keep
  old.to.new <- c()
  old.to.new[k] <- 1:length(k)
  P <- P[k,]
  Cu <- matrix(old.to.new[Cu], nrow(Cu), 2)
  ptri <- matrix(old.to.new[ptri], nrow(ptri), 3)
  
  ## Find lengths of connections
  Ls <- sqrt(apply((P[Cu[,1],] - P[Cu[,2],])^2, 1, sum))

  return(list(nfix=nfix, P=P, Ptri=ptri, Cu=Cu, Ls=Ls))
}

hs <- create.hemisphere()

ptri <- t(hs$Ptri)
P <- hs$P
rgl.clear()
rgl.triangles(P[ptri,1], P[ptri,2], P[ptri,3], col="white", front="line")



