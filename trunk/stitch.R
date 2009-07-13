source("common.R")

## Create points that lie on path p2 at the fraction distances of
## the points along p1
create.corresponding.points <- function(p1, p2, N) {
  ## Remove duplicated points
  ## p1 <- rbind(p1[!is.na(p1[,"l"]),], p1[nrow(p1),])
  ## p2 <- rbind(p2[!is.na(p2[,"l"]),], p2[nrow(p2),])

  s1 <- p1[,"s"]
  F1 <- s1/s1[length(s1)]                  # fractional distance along path 1
  F1 <- F1[(F1!=0) & (F1!=1)]
  s2 <- p2[,"s"]
  F2 <- s2/s2[length(s2)]                  # fractional distance along path 2
  
  P <- find.points.in.path(F1, p2)

  id <- N+1:nrow(P)
  breaks <- which(is.na(p1[,"l"]))
  breaks <- breaks[-length(breaks)]
  for (i in breaks) {
    print(paste("Link to break points detected. Putting in duplicate point", i))
    id[i:nrow(P)] <- id[i:nrow(P)]-1
  }
  P <- cbind(P, id=id)
  P <- cbind(P, Cf=NA)
  P <- cbind(P, Cb=NA)
  p2 <- cbind(p2[,1:ncol(P)], F=F2)
  P <-  cbind(P,  F=F1)
  P <- rbind(p2, P)

  s <- sort(c(F2, F1), index.return=TRUE)
  P <- P[s$ix,]
  return(P)
}

## Stich together a pair of paths tb (in the backwards direction) and
## tf (in the forwards direction) to create a new path and a set of
## correpondances
stitch <- function(tb, tf, P) {
  ## Create points in the backwards tear path that correspond to
  ## points in the forwards path
  Pb <- create.corresponding.points(tf, tb, nrow(P))

  ## Create the forward and backwards connections
  Pb[1:(nrow(Pb)-1),"Cb"] <- Pb[2:nrow(Pb),"id"]
  Pb[2:nrow(Pb),    "Cf"] <- Pb[1:(nrow(Pb)-1),"id"]
  print(Pb)

  ## Create the first column of the correspondences matrix
  corrs <- cbind(Pb[-1, "id"])

  ## Sort by id
  s <- sort(Pb[,"id"], index.return=TRUE)
  Pb <- Pb[s$ix,]
  
  Pb <- Pb[,1:(ncol(Pb)-1)]

  ## Put the new points in the master matrix and update the existing
  ## points, so they have the correct Cf. There should be more elegant
  ## way of computing N.new.
  N.new <- max(max(Pb[,"id"]) - nrow(P), 0)
  if (N.new) {
    P <- rbind(P, Pb[(nrow(Pb)-N.new)+1:N.new,])
  }
  P[Pb[,"id"],] <- Pb

  ## Create point in the forwards tear path that correspond to
  ## points in the backwards path
  Pf <- create.corresponding.points(tb, tf, nrow(P))

  ## Create the forward and backwards connections
  Pf[1:(nrow(Pf)-1),"Cf"] <- Pf[2:nrow(Pf),"id"]
  Pf[2:nrow(Pf),    "Cb"] <- Pf[1:(nrow(Pf)-1),"id"]
  print(Pf)

  if (nrow(Pb) != nrow(Pf)) {
    print("Correspondence error")
  }

  ## Add the second column of the correspondences matrix
  corrs <- cbind(corrs, Pf[-1, "id"])
  
  ## Sort by id
  s <- sort(Pf[,"id"], index.return=TRUE)
  Pf <- Pf[s$ix,]

  Pf <- Pf[,-ncol(Pf)]

  ## Number of new points - should be more elegant way of computing
  ## this
  N.new <- max(max(Pf[,"id"]) - nrow(P), 0)

  ## Put the new points in the master matrix
  ## and update the existing points, so they have the correct Cf
  if (N.new) {
    P <- rbind(P, Pf[(nrow(Pf)-N.new)+1:N.new,])
  }
  P[Pf[,"id"],] <- Pf

  return(list(P=P, corrs=corrs))
}

## Stitch up an entire retina given the tear matrix tm and the
## edge.path
## Each row of the tear matrix has the apex, the
## backwards-most point vertex of the tear and the forwards-most
## vertex
stitch.retina <- function(tm, edge.path) {
  N <- nrow(edge.path) - 1
  P <- edge.path[1:N, 1:2]
  colnames(P) <- c("X", "Y")
  ## Add an index column
  P <- cbind(P, id=1:N,  Cf=c(2:N, 1), Cb=(c(N, 1:(N-1))))
  corrs <- matrix(0, 0, 2)

  ifix <- 1:N
  
  for(i in 1:nrow(tm)) {
    ## Sort out edge1 and edge2 to account for circular nature of path
    ## Edge 1 goes from apex to end1, and is in the backwards (decreasing)
    ## direction
    if (tm[i,2] < tm[i,1]) {
      e1 <- tm[i,1]:tm[i,2]
    } else {
      e1 <- c(tm[i,1]:1, N:tm[i,2])
    }

    ## Edge 2 goes from apex to end1, and should be in the forwads (increasing)
    ## direction
    if (tm[i,3] > tm[i,1]) {
      e2 <- tm[i,1]:tm[i,3]               # Conventional
    } else {
      e2 <- c(tm[i,1]:N, 1:tm[i,3])       # Goes around end of path
    }
    ifix <- c(setdiff(ifix,
                      c(setdiff(e1, tm[i,2]),
                        setdiff(e2, tm[i,3])))) # Get rid of tears from fixed points
    
    ## Check to see if any other tears appear in this tear
    in.tears <- setdiff(which(tm[,1] %in% e1), i)
    for (j in in.tears) {
      Ne1 <- length(e1)
      e1 <- c(e1[1:which(e1==tm[j,3])], NA,
              e1[which(e1==tm[j,2]):Ne1]) 
      ##    print(e1)
    }
    in.tears <- setdiff(which(tm[,1] %in% e2), i)
    for (j in in.tears) {
      Ne2 <- length(e2)
      e2 <- c(e2[1:which(e2==tm[j,2])], NA,
              e2[which(e2==tm[j,3]):Ne2])
      ##    print(e2)
    }
    ##print(e1)  
    ##print(e2)

    p1 <- list()
    for (es in navec.to.list(e1)) {
      p1 <- c(p1, list(P[es,]))
      lines(P[es, 1:2], col="green")
    }
    ##  print(p1)

    p2 <- list()
    for (es in navec.to.list(e2)) {
      print(es)
      p2 <- c(p2, list(P[es,]))
      lines(P[es,], col="blue")    
    }
    ##  print(p2)

    s <- stitch(create.path(p1), create.path(p2), P)
    ## print(s$p)
    segments(s$P[s$corrs[,1],"X"], s$P[s$corrs[,1], "Y"],
             s$P[s$corrs[,2],"X"], s$P[s$corrs[,2], "Y"], col="red")
    P <- s$P
    corrs <- rbind(corrs, s$corrs)
    ## Insert the new points in the path
  }

  
  P <- reorder.path(P)
  trans <- 1:nrow(P)
  trans[P[,'id']] <- 1:nrow(P)
  corrs <- matrix(trans[corrs], ncol=2)
  ifix <- trans[ifix]
  
  segments(P[corrs[,1],"X"], P[corrs[,1], "Y"],
           P[corrs[,2],"X"], P[corrs[,2], "Y"], col="green")
  points(P[c(110, 135, 167), 1:2], col="blue")
  
  return(list(P=P, corrs=corrs, ifix=ifix))
}

## Put the path so that the points are in order around the path
reorder.path <- function(P) {
  p <- P
  id <- 1
  for (i in 1:nrow(P)) {
    p[i,] <- P[id,]
    id <- P[id, "Cf"]
  }
  return(p)
}
