M <- 20
N <- 20
G <- expand.grid(nx=1:M, ny=1:N)
G[,"x"] <- G[,"nx"]-G[,"ny"]/2
G[,"y"] <- G[,"ny"]*sqrt(3)/2

## Connect (x,y) to (x,y+1)
## Connect (x,y) to (x+1,y+1)
## Connect (x,y) to (x,y+1)
## Connectivity matrix: One row for each connection,
## which contains indicies of connected points
m <- M*N
## Start points & endpoints
sp <- (1:m)[(1:m %% M) != 0]
ep <- sp + 1
A <- cbind(sp, ep)
sp <- 1:(m-M)
ep <- sp + M
A <- rbind(A, cbind(sp, ep))
sp <- (1:(m-M))[(1:(m-M) %% M) != 0]
ep <- sp + M + 1
A <- rbind(A, cbind(sp, ep))

##           cbind(1:nrow(G), C[,2]),
##           cbind(1:nrow(G), C[,3]))

plot(G[,"x"], G[,"y"])
segments(G[A[,1],"x"], G[A[,1],"y"],
         G[A[,2],"x"], G[A[,2],"y"])

##x = rep(c(1:10, 0.5+1:10), 4)
##y = matrix(1:8, 8, 10)
##Q = cbind(x, as.vector(y))
##dtQ = delaunayn(Q)
##trimesh(dtQ, Q)

##plot(Q[,1], Q[,2])
##segments(Q[as.vector(t(tQ)),1], Q[as.vector(t(tQ)),2],
##         Q[as.vector(t(tQ[c(2,3,1),])),1], Q[as.vector(t(tQ[c(2,3,1),])),2])
