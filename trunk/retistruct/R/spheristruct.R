## Return sequence of indicies in path between i and j, governed by
## pointer vector p
path <- function(i, j, g, h) {
  if (i == j) {
    return(i)
  } else {
    if (h[i] == i) {
      return(c(i, path(g[i], j, g, h)))
    } else {
      return(c(i, path(h[i], j, g, h)))
    }
  }
}

## Return sequence of indicies in path between i and j, governed by
## pointer vector p
path.length <- function(i, j, g, h, P) {
  if (any(is.na(c(i, j)))) {
    stop("i or j contains NA")
  }
  if (i == j) {
    return(0)
  } else {
    if (h[i] == i) {
      return(sqrt(sum((P[i,] - P[g[i],])^2)) + path.length(g[i], j, g, h, P))
    } else {
      return(path.length(h[i], j, g, h, P))
    }
  }
}

##' Label a set of three unlabelled points supposed to refer to the
##' apex and vertcies of a cut and tear with the V0 (Apex), VF
##' (forward vertex) and VB (backward vertex) labels.
##'
##' @title Label three outline point indicies as apicies and vertices of tear  
##' @param m the vector of three indicies
##' @param gf the forward pointer vector
##' @param gb the backward pointer vector
##' @param P the outline points arranged as a N-by-2 matrix
##' @return Vector of indicies labelled with V0, VF and VB
##' @author David Sterratt
label.tear.points <- function(m, gf, gb, P) {
  ## Each row of this matrix is a permutation of the markers
  p <- rbind(c(1, 2, 3),
             c(1, 3, 2),
             c(2, 1, 3),
             c(2, 3, 1),
             c(3, 1, 2),
             c(3, 2, 1))

  ## For each permuation of V0, VF, VB, measure the sum of length in
  ## the forwards direction from V0 to VF and in the backwards
  ## direction from V0 to VB. The permuation with the minimum distance
  ## is the correct one.
  tplmin <- Inf                      # The minimum path length
  h <- 1:nrow(P)                     # identity correspondence mapping
                                     # used for measuring distances
                                     # (this effectively ignores
                                     # sub-tears, but this doesn't
                                     # matter)
  for (i in 1:nrow(p)) {
    V0 <- m[p[i,1]]
    VF <- m[p[i,2]]
    VB <- m[p[i,3]]
    tpl <- path.length(V0, VF, gf, h, P) + path.length(V0, VB, gb, h, P)
    if (tpl < tplmin) {
      M <- m[p[i,]]
      tplmin <- tpl
    }
  }
  names(M) <- c("V0", "VF", "VB")
  return(M)
}

## Check that tears are all in the correct direction
##
## Given a tear matrix T with columns "V0", "VF", and "VB", check that
## all tears are correct.
##
## Output:
##   If all is OK, returns empty vector
##   If not, returns indicies of problematic tears
##
check.tears <- function(T, gf, gb, P) {
  out <- c()
  if (is.matrix(T)) {
    for (i in 1:nrow(T)) {
      ## Extract the markers for this row
       m <- T[i, c("V0", "VF", "VB")]
       M <- label.tear.points(m, gf, gb, P)
       if (!all(M == m)) {
         out <- c(out, i)
       }
     }
   }
   return(out)
 }

 ## order.Rset(Rset, gf, hf)
 ##
 ## It is nice to create Rset as an ordered set
 order.Rset <- function(Rset, gf, hf) {
   ## To to this, join the path from the first two members of the set.
   R12 <- path(Rset[1], Rset[2], gf, hf)
   R21 <- path(Rset[2], Rset[1], gf, hf)
   Rset <- c(R12[-1], R21[-1])
   return(Rset)
 }

 ## triangulate.outline(P, g=NULL, n=200, h=1:nrow(P))
 ##
 ## Create a triangulation of the outline defined by the points P,
 ## which are represented as N*2 matrix.  If the pointer list
 ## g is supplied, it is used to define the outline. There should be at least n
 ## triangles in the triangulation. Correspondences for any new points added
 ## are added to the h argument
 ## 
 ## Returns a list comprising:
 ## P   - The set of new points, with the existing points at the start
 ## T   - The triangulation
 ## a   - Array containing area of each triangle
 ## A   - Total area of outline
 ## Cu  - Unique set of M connections, as M*2 matrix
 ## L   - Length of each connection
 ## h   - Correspondances mapping
 ##
 triangulate.outline <- function(P, g=NULL, n=200, h=1:nrow(P),
                                suppress.external.steiner=FALSE) {
   ## By default, segments are outline of points in order
   S <- cbind(1:nrow(P), c(2:nrow(P), 1))
   if (!is.null(g)) {
     S <- pointers2segments(g)
   }
   ## Make initial triangulation
   out <- triangulate(pslg(V=P, S=S), Y=TRUE, j=TRUE, Q=TRUE)

   ## Sometimes a point exists which only belongs to one segment. The
   ## point to which it is connected, is itself connected by three
   ## segments. We want to get rid of these points, and the easiest way
   ## is to triangulate without the naughty points.
   i.bad <- which(table(out$S)==1)
   if (length(i.bad) > 0) {
     warning(paste("Bad points:", paste(i.bad, collapse=" ")))
     out <- triangulate(pslg(V=P[-i.bad,], S=S), Y=TRUE, j=TRUE, Q=TRUE)
     P <- out$V
   }

   ## Now determine the area
   A.tot <- sum(with(out, tri.area(P, T)))

   ## Produce refined triangulation
   if (!is.na(n)) {
     out <- triangulate(pslg(V=P, S=S), a=A.tot/n, q=20,q
                       Y=suppress.external.steiner, j=TRUE,
                       Q=TRUE)
  }
  if (any(P != out$P[1:nrow(P),])) {
    stop("Points changed in triangulation")
  }
  P <- out$V
  T <- out$T

  ## Create pointers from segments

  ## To ensure the correct orientaion, we use the fact that the
  ## triangles are all anticlockwise in orinentation, and that the
  ## orientation of the first row of the segment matrix determines the
  ## orientation of all the other rows.

  ## We therefore find the triangle which contains the first segment
  S <- out$S
  T1 <- which(apply(T, 1, function(x) {all(S[1,] %in% x)}))

  ## Then find out which of the vertices in the triangle is not the
  ## one we need
  i <- which((T[T1,] %in% S[1,]) == FALSE)
  if (i == 3) S[1,] <- T[T1,c(1,2)]
  if (i == 2) S[1,] <- T[T1,c(3,1)]
  if (i == 1) S[1,] <- T[T1,c(2,3)]

  ## Now create the pointers from the segments
  gf <- segments2pointers(S)
  gb <- gf
  gb[na.omit(gf)] <- which(!is.na(gf))
  Rset <- na.omit(gf)
  
  ## Derive edge matrix from triangulation
  Cu <- rbind(T[,1:2], T[,2:3], T[,c(3,1)])
  Cu <- Unique(Cu, TRUE)

  ## If we are in the business of refining triangles (i.e. specifying
  ## n), remove lines which join non-ajancent parts of the outline
  if (!is.na(n)) {
    for (i in 1:nrow(Cu)) {
      C1 <- Cu[i,1]
      C2 <- Cu[i,2]
      if (all(Cu[i,] %in% Rset)) {
        if (!((C1 == gf[C2]) ||
              (C2 == gf[C1]))) {
          ## Find triangles containing the line
          ## segments(P[C1,1], P[C1,2], P[C2,1], P[C2,2], col="yellow")
          Tind <- which(apply(T, 1 ,function(x) {(C1 %in% x) && (C2 %in% x)}))
          print(paste("Non-adjacent points in rim connected by line:", C1, C2))
          print(paste("In triangle:", Tind))
          ## Find points T1 & T2 in the two triangles which are not common
          ## with the edge
          T1 <- setdiff(T[Tind[1],], Cu[i,])
          T2 <- setdiff(T[Tind[2],], Cu[i,])
          print(paste("Other points in triangles:", T1, T2))
          ## Create a new point at the centroid of the four verticies
          ## C1, C2, T1, T2
          p <- apply(P[c(C1, C2, T1, T2),], 2, mean)
          P <- rbind(P, p)
          n <- nrow(P)
          ## Remove the two old triangles, and create the four new ones
          T[Tind[1],] <- c(n, C1, T1)
          T[Tind[2],] <- c(n, C1, T2)
          T <- rbind(T,
                     c(n, C2, T1),
                     c(n, C2, T2))
        }
      }
    }

    ## Add the new points to the correspondances vector
    h <- c(h, (length(h)+1):nrow(P))

    ## Create the edge matrix from the triangulation
    Cu <- rbind(T[,1:2], T[,2:3], T[,c(3,1)])
    Cu <- Unique(Cu, TRUE)
  }

  ## Swap orientation of triangles which have clockwise orientation
  A.signed <- tri.area.signed(P, T)
  T[A.signed<0,c(2,3)] <- T[A.signed<0,c(3,2)]
  A <- abs(A.signed)
    
  ## Find lengths of connections
  L <- vecnorm(P[Cu[,1],] - P[Cu[,2],])

  ## Check there are no zero-length lines
  if (any(L==0)) {
    print("WARNING: zero-length lines")
  }

  return(list(P=P, T=T, Cu=Cu, h=h,  A=A, L=L,
              A.signed=A.signed, A.tot=A.tot,
              gf=gf, gb=gb, S=out$S, E=out$E, EB=out$EB))
}

## Stitch together tears in an outline
##
## Input arguments:
## P     - the coordinates of points in a mesh, including those the outline
## gf    - the forward pointer list
## gb    - the backward pointer list
## V0    - indicies of the apex of each tear
## VF    - indicies of the forward vertex of each tear
## VB    - indicies of the backward vertex of each tear
## i0    - the index of the landmark; this needs to be in the rim
##
## The function returns a list contatining:
## Rset  - the set of points on the rim
## i0    - the index of the landmark
## P     - a new set of meshpoints
## V0    - indicies of the apex of each tear
## VF    - indicies of the forward vertex of each tear
## VB    - indicies of the backward vertex of each tear
## TFset - list containing indicies of points in each foward tear
## TBset - list containing indicies of points in each backward tear
## gf    - new forward pointer list
## gb    - new backward pointer list
## h     - correspondence mapping
## hf    - correspondence mapping in forward direction for
##         points on boundary
## hb    - correspondence mapping in backward direction for
##         points on boundary
##
stitch.outline <- function(P, gf, gb, V0, VB, VF, i0=NA) {
  ## Create initial sets of correspondances
  N <- nrow(P)                          # Number of points
  h <- 1:N                              # Initial correspondences
  hf <- h
  hb <- h
  M <- length(V0)                       # Number of tears
  i.parent <- rep(0, M)                 # Index of parent tear.
                                        # Is 0 if root otherwise
                                        # index of tear if in forward side
                                        # or negative index if in backward side 

  ## Initialise the set of points in the rim
  ## We don't assume that P is the entire set of points; instead
  ## get this information from the pointer list.
  Rset <- na.omit(gf)
  
  ## Create lists of forward and backward tears
  TFset <- list()
  TBset <- list()
  
  ## Iterate through the tears to create tear sets and rim set
  for (j in 1:M) {
    ## Create sets of points for each tear and remove these points from
    ## the rim set
    ##message(paste("Forward tear", j))
    TFset[[j]] <- mod1(path(V0[j], VF[j], gf, h), N)
    TBset[[j]] <- mod1(path(V0[j], VB[j], gb, h), N)
    Rset <- setdiff(Rset, setdiff(TFset[[j]], VF[j]))
    Rset <- setdiff(Rset, setdiff(TBset[[j]], VB[j]))
  }

  ## Search for parent tears
  ## Go through all tears
  for (j in 1:M) {
    for (k in setdiff(1:M, j)) {
      ## If this tear is contained in a forward tear
      if (all(c(V0[j], VF[j], VB[j]) %in% TFset[[k]])) {
        i.parent[j] <- k
        message(paste("Tear", j, "child of forward side of tear", k))
        ## Set the forward pointer
        hf[VB[j]] <- VF[j]
        ## Remove the child tear points from the parent
        TFset[[k]] <- setdiff(TFset[[k]],
                              setdiff(c(TBset[[j]], TFset[[j]]), c(VB[j], VF[j])))
        print(TFset[[k]])
      }
      ## If this tear is contained in a backward tear
      if (all(c(V0[j], VF[j], VB[j]) %in% TBset[[k]])) {
        i.parent[j] <- -k
        message(paste("Tear", j, "child of backward side of tear", k))
        ## Set the forward pointer
        hb[VF[j]] <- VB[j]
        ## Remove the child tear points from the parent
        TBset[[k]] <- setdiff(TBset[[k]],
                              setdiff(c(TBset[[j]], TFset[[j]]), c(VB[j], VF[j])))
      }
    }
    if (i.parent[j] == 0) {
      message(paste("Tear", j, "child of rim"))
      hf[VB[j]] <- VF[j]
      hb[VF[j]] <- VB[j]
    }
  }

  ## If not set, set the landmark marker index. Otherwise
  ## check it
  if (is.na(i0)) {
    i0 <- Rset[1]
  } else {
    if (!(i0 %in% Rset)) {
      return(NULL)
    }
  }

  ## Insert points on the backward tears corresponding to points on
  ## the forward tears
  sF <-      stitch.insert.points(P, V0, VF, VB, TFset, TBset,
                                  gf, gb, hf, hb, h,
                                  "Forwards")

  ## Insert points on the forward tears corresponding to points on
  ## the backward tears
  sB <- with(sF,
             stitch.insert.points(P, V0, VB, VF, TBset, TFset,
                                  gb, gf, hb, hf, h,
                                  "Backwards"))
  ## Extract data from object
  P <- sB$P
  gf <- sB$gb
  gb <- sB$gf
  hf <- sB$hb
  hb <- sB$hf
  h <- sB$h

  ## Link up points on rim
  h[Rset] <- hf[Rset]
  
  ## Make sure that there are no chains of correspondences
  while (!all(h==h[h])) {
   h <- h[h]
  }
  
  return(list(Rset=Rset, i0=i0,
              VF=VF, VB=VB, V0=V0,
              TFset=TFset, TBset=TBset,
              P=P, h=h, hf=hf, hb=hb,
              gf=gf, gb=gb))
}

## Inner function responsible for inserting the points
stitch.insert.points <- function(P, V0, VF, VB, TFset, TBset, gf, gb, hf, hb, h,
                                 dir) {
  M <- length(V0)                       # Number of tears
  ## Iterate through tears to insert new points
  for (j in 1:M) {
    ## Compute the total path length along each side of the tear
    Sf <- path.length(V0[j], VF[j], gf, hf, P)
    Sb <- path.length(V0[j], VB[j], gb, hb, P)
    message(paste("Tear", j, ": Sf =", Sf, "; Sb =", Sb))

    ## For each point in the forward path, create one in the backwards
    ## path at the same fractional location
    message(paste("   ", dir, " path", sep=""))
    for (i in setdiff(TFset[[j]], c(V0[j], VF[j]))) {
      sf <- path.length(V0[j], i, gf, hf, P)
      ## If the point isn't at the apex, insert a point
      if (sf > 0) {
        message(paste("    i =", i,
                                "; sf/Sf =", sf/Sf,
                                "; sf =", sf))
        for (k in TBset[[j]]) {
          sb <- path.length(V0[j], k, gb, hb, P)
          message(paste("      k =", format(k, width=4),
                                  "; sb/Sb =", sb/Sb,
                                  "; sb =", sb))
          if (sb/Sb > sf/Sf) {
            break;
          }
          k0 <- k
          sb0 <- sb
        }

        ## If this point does not point to another, create a new point
        if ((hf[i] == i)) {
          f <- (sf/Sf*Sb-sb0)/(sb-sb0)
          message(paste("      Creating new point: f =", f))
          ## browser(expr=is.infinite(f))
          p <- (1-f) * P[k0,] + f * P[k,]
          ## browser(expr=any(is.nan(p)))

          ## Find the index of any row of P that matches p
          n <- anyDuplicated(rbind(P, p), fromLast=TRUE) 
          if (n == 0) {
            ## If the point p doesn't exist
            P <- rbind(P, p)
            ## Update forward and backward pointers
            n <- nrow(P)                    # Index of new point
            gb[n]     <- k
            gf[n]     <- gf[k]
            gb[gf[k]] <- n
            gf[k]     <- n

            ## Update correspondences
            hf[n] <- n
            hb[n] <- n
            h[i] <- n
            h[n] <- n
          } else {
            warning(paste("Point", n, "already exists"))
          }
        } else {
          ## If not creating a point, set the point to point to the forward pointer 
          h[i] <- hf[i]
        }
      }
    } 
  }
  return(list(P=P, hf=hf, hb=hb, gf=gf, gb=gb, h=h))
}

## Convert a matrix containing on each line the indicies of the points
## forming a segment, and convert this to two sets of ordered pointers
segments2pointers <- function(S) {
  g <- c()
  j <- 1                                # Row of S
  k <- 1                                # Column of S
  while(nrow(S) > 0) {
    i <- S[j,3-k]                       # i is index of the next point
    g[S[j,k]] <- i                      # Set the pointer to i
    S <- S[-j,,drop=FALSE]              # We have used this row of S
    if (nrow(S) == 0) {
      return(g)
    }
    j <- which(S[,1] == i)            # Is i in the first column of S?
    if (length(j) > 1) {
      stop("The segment list is not valid as it contains an element more than twice.")
    }
    if (length(j)) {              # If so, set the current column to 1
      k <- 1
    } else {
      j <- which(S[,2] == i) # Otherwise, look for i in the second column
      k <- 2
      if (!length(j)) {
        stop(paste("No matching index for point", i, "in S."))
        return(NULL)
      }
    }
  }
  return(g)
}

## Convert a set of ordered pointers to a matrix containing on each
## line the indicies of the points forming a segment
pointers2segments <- function(g) {
  S1 <- which(!is.na(g))
  S2 <- g[S1]
  return(cbind(S1, S2))
}

## Merge points and edges before optimising the elastic energy
## The information to be merged is contained in the
## list t, and the index of the landmark i0
##
## The information includes 
## h    - the point correspondence mapping
## Rset - the set of points on the rim
## i0   - the index of the landmark
## T    - the triangulation
## Cu   - the edge matrix
## L    - the edge lengths
## P    - the point locations
##
## The function returns merged and transformed versions of all these
## objects (all suffixed with t), as well as a matrix Bt, which maps a
## binary vector representation of edge indicies onto a binary vector
## representation of the indicies of the points linked by the edge
merge.points.edges <- function(t) {
  h <- t$h
  T <- t$T
  Cu <- t$Cu
  L <- t$L
  P <- t$P
  gf <- t$gf
  
  ## Form the mapping from a new set of consecutive indicies
  ## the existing indicies onto the existing indicies
  u <- unique(h)

  ## Transform the point set into the new indicies
  Pt  <- P[u,]

  ## Transform the point correspondance mapping to the new index space  
  ht <- c()
  for (i in 1:length(h)) {
    ht[i] <- which(u == h[i])
  }

  ## DOESN'T WORK
  ## Form the inverse mapping from the existing indicies to the new
  ## set of consecutive indicies
  ## uinv <- c()
  ## uinv[u] <- 1:length(u)
  ## ht <- uinv[h[u]]

  ## Transform the triangulation to the new index space
  Tt  <- matrix(ht[T], ncol=3)

  ## Tansform the forward pointer into the new indicies
  gft <- ht[gf]

  ## Determine H, the mapping from edges onto corresponding edges
  Cut <- matrix(ht[Cu], ncol=2)
  Cut <- t(apply(Cut, 1, sort))
  M <- nrow(Cut)
  H <- rep(0, M)
  for (i in 1:M) {
    if (!H[i]) {
      H[i] <- i
      for (j in i:M) {
        if (identical(Cut[i,], Cut[j,])) {
          H[j] <- i
        }
      }
    }
  }

  ## Form the mapping from a new set of consecutive edge indicies
  ## onto the existing edge indicies
  U <- unique(H)

  ## Transform the edge set into the new indicies
  Cut <- Cut[U,]

  ## Transform the edge correspondance mapping to the new index space  
  Ht <- c()
  for (i in 1:length(H)) {
    Ht[i] <- which(U == H[i])
  }

  ## Create the lengths of the merge edges by averaging
  Lt <- c()
  for (k in 1:length(U)) {
    is <- which(Ht == k)
    ## if (length(is)>1) {
    ##   print(L[is])
    ## }
    Lt[k] <- mean(L[is])
  }

  ## Transform the rim set
  Rsett <- unique(ht[t$Rset])
  i0t <- ht[t$i0]

  ## Create the symmetric connection set
  Ct <- rbind(Cut, Cut[,2:1])

  ## Matrix to map line segments onto the points they link
  Bt <- matrix(0, nrow(Pt), nrow(Ct))
  for (i in 1:nrow(Ct)) {
    Bt[Ct[i,1],i] <- 1
  }
  
  return(list(Pt=Pt, Tt=Tt, Ct=Ct, Cut=Cut, Bt=Bt, Lt=Lt, ht=ht,
              Rsett=Rsett, i0t=i0t, P=P, H=H, Ht=Ht))
}

## Stretch the mesh in the flat retina to a circular outline
##
## Arguments:
## Cu - connection matrix
## L  - lengths
## i.fix - indicies of fixed points
## P.fix - coordinates of fixed points
##
## Returns:
## New matrix of point locations
stretch.mesh <- function(Cu, L, i.fix, P.fix) {
  N <- max(Cu)
  M <- length(L)
  C <- matrix(0, 2*N, 2*N)
  for (i in 1:M) {
    C[2*Cu[i,1]-1:0,2*Cu[i,2]-1:0] <- diag(2) / L[i]
    C[2*Cu[i,2]-1:0,2*Cu[i,1]-1:0] <- diag(2) / L[i]
  }

  ## FIXME: This debugging output should probably go once the "Lapack
  ## routine dgesv: system is exactly singular" bug has been
  ## vanquished
  print(det(C))
  dupC <- duplicated(C)
  if (any(dupC)) {
    i <- which(dupC)
    Ci <- C[i,]
    dups <- which(apply(C, 1, function(x) {identical(x, Ci)}))
    print(dups)
    print(Ci)
    for (d in dups) {
      print(d)
      print(which(C[d,]==1))
    }
  }
  
  ind <- as.vector(rbind(2*i.fix-1, 2*i.fix))
  A <- C[-ind, -ind]
  print(det(A))
  B <- C[-ind,  ind]
  P <- matrix(t(P.fix), ncol(B), 1)
  D <- diag(apply(cbind(A, 2*B), 1, sum))
  print(det(D))

  Q <- 2 * solve(D - A) %*% B %*% P
  Q <- matrix(Q, nrow(Q)/2, 2, byrow=TRUE)
  R <- matrix(0, nrow(Q) + length(i.fix), 2)
  R[i.fix,] <- P.fix
  R[-i.fix,] <- Q
  return(R)
}

## Project mesh points in the flat outline onto a sphere
##
## The information to be merged is contained in the list r and the
## lattitude at which the sphere is cut off is given in phi0
##
## The information includes:
## Pt     - the mesh point coordinates
## Rsett  - the set of points on the rim
## A.tot  - the area of the flat outline
##
## The function returns a list with the following members:
## phi    - lattitude of mesh points
## lambda - longitude of mesh points
## R      - radius of sphere
## phi0   - lattitude at which sphere is cut off (from input)
project.to.sphere <- function(r, phi0=50*pi/180, lambda0=lambda0) {
  Pt <- r$Pt
  Rsett <- r$Rsett
  i0t <- r$i0t
  A.tot <- r$A.tot
  Cut <- r$Cut
  Lt <- r$Lt
  
  Nt <- nrow(Pt)
  Nphi <- Nt - length(Rsett)

  ## From this we can infer what the radius should be from the formula
  ## for the area of a sphere which is cut off at a lattitude of phi0
  ## area = 2 * PI * R^2 * (sin(phi0)+1)
  R <- sqrt(A.tot/(2*pi*(sin(phi0)+1)))
  
  ## Stretch mesh points to circle
  Ps <- stretch.mesh(Cut, Lt, Rsett, circle(length(Rsett)))
  x <- Ps[,1]
  y <- Ps[,2]
  phi <- -pi/2 + sqrt(x^2 + y^2)*(phi0+pi/2)
  phi[Rsett] <- phi0
  lambda <- atan2(y, x)
  lambda <- lambda - lambda[i0t] + lambda0

  return(list(phi=phi, lambda=lambda, R=R, phi0=phi0, lambda0=lambda0, Ps=Ps))
}

##
## Energy/error functions
## 

## Formula for central angle
central.angle <- function(phi1, lambda1, phi2, lambda2) {
  return(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2)))
}

## Calculate lengths of connections on sphere
compute.lengths <- function(phi, lambda, Cu, R) {
  ## Use the upper triagular part of the connectivity matrix Cu
  phi1    <- phi[Cu[,1]]
  lambda1 <- lambda[Cu[,1]]
  phi2    <- phi[Cu[,2]]
  lambda2 <- lambda[Cu[,2]]
  l <- R*central.angle(phi1, lambda1, phi2, lambda2)

  return(l)
}

## Calculate lengths of connections on sphere
compute.areas <- function(phi, lambda, T, R) {
  P <- R * cbind(cos(phi)*cos(lambda),
                 cos(phi)*sin(lambda),
                 sin(phi))

  ## Find areas of all triangles
  areas <- -0.5/R * dot(P[T[,1],], extprod3d(P[T[,2],], P[T[,3],]))

  return(areas)
}

f <- function(x, x0=0.1, k=1) {
  y <- x

  c1 <- x <= 0
  c2 <- (0 < x) & (x < x0)
  c3 <- x0 <= x

  y[c1] <- -k*(x[c1] - x0/2)
  y[c2] <- k/2/x0*(x0 - x[c2])^2
  y[c3] <- 0

  return(y)
}

fp <- function(x, x0=0.1, k=1) {
  y <- x

  c1 <- x <= 0
  c2 <- (0 < x) & (x < x0)
  c3 <- x0 <= x

  y[c1] <- -k
  y[c2] <- -k/x0*(x0 - x[c2])
  y[c3] <- 0

  return(y)
}

## Now for the dreaded elastic error function....
E <- function(p, Cu, C, L, B, T, A, R, Rset, i0, phi0, lambda0, Nphi,
              E0.A=1, k.A=1, x0=0.1, N, verbose=FALSE) {
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- rep(lambda0, N)
  lambda[-i0] <- p[Nphi+1:(N-1)]

  ##
  ## Compute derivative of elastic energy
  ##
  ## Use the upper triagular part of the connectivity matrix Cu
  phi1    <- phi[Cu[,1]]
  lambda1 <- lambda[Cu[,1]]
  phi2    <- phi[Cu[,2]]
  lambda2 <- lambda[Cu[,2]]
  l <- R*central.angle(phi1, lambda1, phi2, lambda2)
  if (verbose==2) {
    print(l)
  }
  E.E <- 0.5 * sum((l - L)^2/L)
  ## E.E <- 0.5 * sum((l - L)^2)
  ## E.E <- 0.5*sum((l)^2/L)
  if (verbose>=1) {
    print(E.E)
  }

  E.A <- 0
  if (E0.A) {
    ##
    ## Compute areas
    ##
    P <- R * cbind(cos(phi)*cos(lambda),
                   cos(phi)*sin(lambda),
                   sin(phi))

    ## Find areas of all triangles
    a <- -0.5/R * dot(P[T[,1],], extprod3d(P[T[,2],], P[T[,3],]))
    ##E.A <- 0.5 * sum((areas - A)^2/A)
    ## E.A <- sum(exp(-k.A*areas/A))
    E.A <- sum(sqrt(A)*f(a/A, x0=x0))
  }
  
  return(E.E + E0.A*E.A)
}

## ... and the even more dreaded gradient of the elastic error
dE <- function(p, Cu, C, L, B, T, A, R, Rset, i0, phi0, lambda0, Nphi,
               E0.A=1, k.A=1, x0=0.1, N, verbose=FALSE) {
  phi <- rep(phi0, N)
  phi[-Rset] <- p[1:Nphi]
  lambda <- rep(lambda0, N)
  lambda[-i0] <- p[Nphi+1:(N-1)]

  ##
  ## Compute derivative of elastic energy
  ##
  phii   <- phi[C[,1]]
  lambdai <- lambda[C[,1]]
  phij    <- phi[C[,2]]
  lambdaj <- lambda[C[,2]]
  ## x is the argument of the acos in the central angle
  u <- sin(phii)*sin(phij) + cos(phii)*cos(phij)*cos(lambdai-lambdaj)
  ## the central angle
  l <- R * acos(u)
  if (verbose==2) {
    print(l)
  }
  fac <- R * (l - c(L, L))/(c(L, L) * sqrt(1-u^2))
  ## fac <- R * (l - c(L, L))/(sqrt(1-u^2))
  ## fac <- R * (l)/(c(L, L) * sqrt(1-u^2))
  dE.E.phii     <- B %*% (fac * (sin(phii)*cos(phij)*cos(lambdai-lambdaj)
                               - cos(phii)*sin(phij)))
  dE.E.dlambdai <- B %*% (fac * cos(phii)*cos(phij)*sin(lambdai-lambdaj))

  dE.A.dphi <- rep(0, N)
  dE.A.dlambda <- 0
  if (E0.A) {
    ##
    ## Compute derivative of areas
    ##
    P <- R * cbind(cos(phi)*cos(lambda),
                   cos(phi)*sin(lambda),
                   sin(phi))

    ## expand triangulation
    T <- rbind(T, T[,c(2,3,1)], T[,c(3,1,2)])
    A  <- c(A, A, A)

    ## Slow way of computing gradient
    dAdPt1 <- -0.5/R * extprod3d(P[T[,2],], P[T[,3],])

    ## Find areas of all triangles
    a <- dot(P[T[,1],], dAdPt1)

##     dEdPt1 <- (areas - A)/A * dAdPt1
##    dEdPt1 <- -k.A/A*exp(-k.A*areas/A) * dAdPt1
    dEdPt1 <- fp(a/A, x0=x0)/sqrt(A) * dAdPt1
    
    Pt1topi <- matrix(0, length(phi), nrow(T))
    for(m in 1:nrow(T)) {
      Pt1topi[T[m,1],m] <- 1
    }
    dEdpi <- Pt1topi %*% dEdPt1
    dpidphi <- R * cbind(-sin(phi) * cos(lambda),
                         -sin(phi) * sin(lambda),
                         cos(phi))
    dpidlambda <- R * cbind(-cos(phi) * sin(lambda),
                            cos(phi) * cos(lambda),
                            0)

    dE.A.dphi    <- rowSums(dEdpi * dpidphi)
    dE.A.dlambda <- rowSums(dEdpi * dpidlambda)
  }
  
  return(c(dE.E.phii[-Rset]      + E0.A * dE.A.dphi[-Rset],
           dE.E.dlambdai[-i0]    + E0.A * dE.A.dlambda[-i0]))
}

## flipped.triangles - Determine indicies of triangles that are flipped
flipped.triangles <- function(phi, lambda, Tt, R) {
  P <- sphere.spherical.to.sphere.cart(phi, lambda, R)
  ## Plot any flipped triangles
  ## First find verticies and find centres and normals of the triangles
  P1 <- P[Tt[,1],]
  P2 <- P[Tt[,2],]
  P3 <- P[Tt[,3],]
  cents <- (P1 + P2 + P3)/3
  normals <- 0.5 * extprod3d(P2 - P1, P3 - P2)

  areas <- -0.5/R * dot(P1, extprod3d(P2, P3))
  print(length(areas))
  
  flipped <- (-dot(cents, normals) < 0)
  return(list(flipped=flipped, cents=cents, areas=areas))
}
  
## Grand optimisation function
optimise.mapping <- function(r, E0.A=10, k.A=1, x0=0.5, method="BFGS",
                             plot.3d=FALSE, dev.grid=NA, dev.polar=NA) {
  phi <- r$phi
  lambda <- r$lambda
  R <- r$R
  phi0 <- r$phi0
  lambda0 <- r$lambda0
  Tt <- r$Tt
  A <- r$A
  Cut <- r$Cut
  Ct <- r$Ct
  Pt <- r$Pt
  Lt <- r$Lt
  Bt <- r$Bt
  Rsett <- r$Rsett
  i0t <- r$i0t
  Nt <- nrow(Pt)  
  Nphi <- Nt - length(Rsett)
  
  ## Optimisation and plotting 
  opt <- list()
  opt$p <- c(phi[-Rsett], lambda[-i0t])
  opt$conv <- 1

  ## iE <- E(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
  ##         E0.A=E0.A, k.A=k.A, N=Nt,
  ##         Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi)

  ## idE <- dE(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
  ##         E0.A=E0.A, k.A=k.A, N=Nt,
  ##         Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi)
  ## idE2 <- dE2(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
  ##         E0.A=E0.A, k.A=k.A, N=Nt,
  ##         Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi, verbose=TRUE)

  ## idE2c <- dE2c(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
  ##         E0.A=E0.A, k.A=k.A, N=Nt,
  ##         Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi)
  ## print(iE)
  
  ## print(idE[1:10])
  ## print("R dE2")
  ## print(idE2[1:10])
  ## print("C dE2")
  ## print(idE2c[1:10])
  ## print("Ratio")
  ## print(idE2c/idE2)
  
  while (opt$conv) {
    ## Optimise
    opt <- optim(opt$p, E, gr=dE,
                 method=method,
                 T=Tt, A=A, Cu=Cut, C=Ct, L=Lt, B=Bt, R=R,
                 E0.A=E0.A, k.A=k.A, N=Nt, x0=x0,
                 Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi,
                 verbose=FALSE)

    ## Report
    E.tot <- E(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
               E0.A=E0.A, k.A=k.A, N=Nt, x0=x0,
               Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi)
    E.l <- E(opt$p, Cu=Cut, C=Ct, L=Lt, B=Bt,  R=R, T=Tt, A=A,
               E0.A=0, k.A=k.A, N=Nt, x0=x0,
               Rset=Rsett, i0=i0t, phi0=phi0, lambda0=lambda0, Nphi=Nphi)
    print(paste("Total error:", E.tot, opt$value, "; Length error:", E.l))
    ft <- flipped.triangles(phi, lambda, Tt, R)
    nflip <- sum(ft$flipped)
    print(paste(nflip, "flipped triangles:"))
    print(which(ft$flipped))
    print("Areas")
    print(ft$areas[ft$flipped])
    print(A[ft$flipped])
    print(length(A))
    
    ## Decode p vector
    phi          <- rep(phi0, Nt)
    phi[-Rsett]  <- opt$p[1:Nphi]
    lambda       <- rep(lambda0, Nt)
    lambda[-i0t] <- opt$p[Nphi+1:(Nt-1)]

    ## Plot
    if (plot.3d) {
      plot.sphere.spherical(phi, lambda, R, Tt, Rsett)
      plot.outline.spherical(phi, lambda, R, r$gb, r$ht)
    }

    if (!is.na(dev.grid)) {
      dev.set(dev.grid)
      with(r, plot.outline.flat(P, gb))
      plot.gridlines.flat(r$P, r$T, phi, lambda, Tt, phi0*180/pi)
    }

    if (!is.na(dev.polar)) {
      dev.set(dev.polar)
      plot.polar(phi0 * 180/pi)
      r$phi <- phi
      r$lambda <- lambda
      plot.outline.polar(r)
    }
  }
  return(list(phi=phi, lambda=lambda, opt=opt, nflip=sum(ft$flipped),
              E.tot=E.tot, E.l=E.l))
}

## Function to plot the fractional change in length of connections 
compute.strain <- function(r) {
  ## Original lengths in flattened outline is a vector with
  ## M elements, the number of rows of Cu
  L <- r$L
  ## New lenghts in reconstructed object is a vector wtih Mt < M
  ## elements, the number of rows of Cut
  ls <- compute.lengths(r$phi, r$lambda, r$Cut, r$R)
  ## For each connection in the flattened object, we want the length of
  ## the corresponding connection in the reconstructed object
  ## The mapping Ht achieves this@
  l <- ls[r$Ht]
  strain <- l/L
  logstrainnorm <- -log(strain)/max(abs(log(strain)))
  return(data.frame(L=L, l=l, strain=strain, logstrainnorm=logstrainnorm))
}

##' Fold outline onto surface of sphere. Folding proceeds in a number
##' of stages:
##'
##' \enumerate{
##' 
##' \item The flat object is triangulated with at least \code{n}
##' triangles. This can introduce new vertices in the rim. 
##'
##' \item The triangulated object is stitched.
##'
##' \item The stitched object is triangulated again, but this time it
##' is not permitted to add extra vertices to the rim.
##'
##' \item The corresponding points determined by the stitching process
##' are merged to form a new set of merged points and a new
##' triangulation.
##'
##' \item The merged points are projected roughly to a sphere.
##'
##' \item The locations of the points on the sphere are moved so as to
##' minimise the energy function.
##' }
##'
##' @title Fold outline onto surface of sphere
##' @param o Outline list, containing the following information:\describe{
##' \item{\code{P}}{outline points as N-by-2 matrix}
##' \item{\code{V0}}{indicies of the apex of each tear}
##' \item{\code{VF}}{indicies of the forward vertex of each tear}
##' \item{\code{VB}}{indicies of the backward vertex of each tear}
##' \item{\code{i0}}{index of the landmark on the rim}
##' \item{\code{phi0}}{lattitude of rim of partial sphere}
##' \item{\code{lambda0}}{longitude of landmark on rim}
##' }
##' @param n Number of points in triangulation.
##' @param report Function used to report progress.
##' @param plot.3d Whether to show 3D picture during optimisation.
##' @param dev.grid Device to plot grid onto. Value of \code{NA} (default)
##' means no plotting.
##' @param dev.polar Device to plot polar plot onto. Value of NA
##' (default) means no plotting.
##' @return Reconstruction object containing the input information and
##' the following modified and extra information:
##' \item{\code{P}}{New set of points in flattened object}
##' \item{\code{gf}}{New set of forward pointers in flattened object}
##' \item{\code{gb}}{New set of backward pointers in flattened object}
##' \item{\code{phi}}{lattitude of new points on sphere}
##' \item{\code{lambda}}{longitude of new points on sphere}
##' \item{\code{Tt}}{New triangulation}
##' @author David Sterratt
fold.outline <- function(o, 
                         n=500,
                         report=print,
                         plot.3d=FALSE, dev.grid=NA, dev.polar=NA) {
  ## Clear polar plot, if it's required
  if (!is.na(dev.polar)) {
    dev.set(dev.polar)
    plot.polar(o$phi0)
  }
  
  report("Triangulating...")
  t <- with(o, triangulate.outline(P, h=1:nrow(P), n=n))
  t <- merge(t, o)
  if (!is.na(dev.grid)) {
    dev.set(dev.grid)
    with(t, trimesh(T, P, col="black"))
  }
    
  report("Stitching...")
  s <- with(t, stitch.outline(P, gf, gb, V0, VB, VF, i0))
  s <- merge(s, t)
  if (is.null(s)) {
    stop("Fixed point is not on the rim")
  }
  if (!is.na(dev.grid)) {
    dev.set(dev.grid)
    plot.stitch.flat(s)
  }

  report("Triangulating...")  
  t <- triangulate.outline(s$P, h=s$h, g=s$gf, n=n,
                           suppress.external.steiner=TRUE)
  r <- merge(t, s)
  if (!is.na(dev.grid)) {
    dev.set(dev.grid)
    plot.stitch.flat(r)
    with(r, trimesh(T, P, col="grey", add=TRUE))
  }
  r$Rset <- order.Rset(r$Rset, r$gf, r$hf)

  report("Merging points...")
  m <- merge.points.edges(r)
  r <- merge(m, r)

  report("Projecting to sphere...")
  p <- project.to.sphere(r, phi0=r$phi0*pi/180, lambda0=r$lambda0*pi/180)
  r <- merge(p, r)
  
  if (!is.na(dev.grid)) {
    ## Plot of initial gridlines
    dev.set(dev.grid)
    with(r, plot.outline.flat(P, gb))
    with(r, plot.gridlines.flat(P, T, phi, lambda, Tt, phi0*180/pi))
    
    ## Initial plot in 3D space
    plot.sphere.spherical(r$phi, r$lambda, r$R, r$Tt, r$Rsett)
    plot.outline.spherical(r$phi, r$lambda, r$R, r$gb, r$ht)
  }

  report("Optimising mapping...")
  ## This pass is experimental - ideally it wouldn't be here, but until
  ## we can fix problems with fixed triangles, it must stay.
  ##o <- optimise.mapping(r, E0.A=exp(1), k.A=2,
  ##                        plot.3d=plot.3d,
  ## dev.grid=dev.grid, dev.polar=dev.polar)
  ## r <- merge.lists(r, o)

  ## This pass is original
  ## o <- list(nflip=1)
  E0.A <- 32
  ##while(o$nflip>0) {
  o <- optimise.mapping(r, E0.A=E0.A,
                        plot.3d=plot.3d,
                        dev.grid=dev.grid, dev.polar=dev.polar)
  r <- merge(o, r)
  ## o <- fem.optimise.mapping(r, nu=0.45,
  ##                           plot.3d=plot.3d,
  ##                           dev.grid=dev.grid, dev.polar=dev.polar)
  ## r <- merge.lists(r, o)
  ##  E0.A <- E0.A * 2;
  ##  }
  
  report(paste("Mapping optimised. Error:", format(r$opt$value,5),
               ";", r$nflip, "flipped triangles."))
  return(r)
}

## Infer coordinates of datapoints
## Arguments:
## r    - object returned by fold.outline
## Ds   - list of sets of datapoints to plot
## Ss   - list of sets of landmarks to plot
##
## Returns:
## New object r, with new objects attached:
## 
## Dsb  - Datapoints in barycentric coordinates
## Dsc  - Datapoints on reconstructed sphere in cartesian coordinates
## Dss  - Datapoints on reconstructed sphere in spherical coordinates
## Ssb  - Landmarks in barycentric coordinates
## Ssc  - Landmarks on reconstructed sphere in cartesian coordinates
## Sss  - Landmarks on reconstructed sphere in spherical coordinates
##
infer.datapoint.landmark.coordinates <- function(r, report=print) {
  report("Inferring coordinates of datapoints")
  Dsb <- list() # Datapoints in barycentric coordinates
  Dsc <- list() # Datapoints on reconstructed sphere in cartesian coordinates
  Dss <- list() # Datapoints on reconstructed sphere in spherical coordinates
  if (!is.null(r$Ds)) {
    for (name in names(r$Ds)) {
      Dsb[[name]] <- tsearchn(r$P, r$T, r$Ds[[name]])
      Dsc[[name]] <- bary.to.sphere.cart(r$phi, r$lambda, r$R, r$Tt, Dsb[[name]])
      Dss[[name]] <- sphere.cart.to.sphere.spherical(Dsc[[name]], r$R)
    }
  }

  report("Inferring coordinates of landmarks")
  Ssb <- list() # Landmarks in barycentric coordinates
  Ssc <- list() # Landmarks on reconstructed sphere in cartesian coordinates
  Sss <- list() # Landmarks on reconstructed sphere in spherical coordinates
  if (!is.null(r$Ss)) {
    for (name in names(r$Ss)) {
      Ssb[[name]] <- with(r, tsearchn(P, T, r$Ss[[name]]))
      Ssc[[name]] <- bary.to.sphere.cart(r$phi, r$lambda, r$R, r$Tt, Ssb[[name]])
      Sss[[name]] <- sphere.cart.to.sphere.spherical(Ssc[[name]], r$R)
    }
  }

  r <- merge(list(Dsb=Dsb, Dsc=Dsc, Dss=Dss,
                  Ssb=Ssb, Ssc=Ssc, Sss=Sss), r)
}

## Infer coordinates of tears
## Arguments:
## r    - object returned by fold.outline
##
## Returns:
## New object r, with new object attached:
## 
## Tss  - Tear coordinates on reconstructed sphere in spherical coordinates
##
infer.tear.coordinates <- function(r,
                                   report=print) {
  report("Inferring coordinates of tears")
  Tss <- list() # Landmarks on reconstructed sphere in spherical coordinates
  if (!is.null(r$TFset)) {
    for (TF in r$TFset) {
      ## Convert indicies to the spherical frame of reference
      j <- r$ht[TF]
      Tss <- c(Tss, list(cbind(phi=r$phi[j], lambda=r$lambda[j])))
    }
  }
  r <- merge(list(Tss=Tss), r)
}
