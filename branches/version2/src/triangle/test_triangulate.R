source("triangle.R")
if (TRUE) {
  P <- cbind(c(0, 0),
             c(0, 1.2),
             c(0.5, 1.2),
             c(0.7, 1.7),
             c(0.8, 1.1),
             c(1, 1),
             c(1, 0))
} else {
  P <- cbind(c(0, 0),
             c(0, 1),
             c(0.5, 1,5),
             c(1, 1),
             c(1, 0))
}
out <- triangulate(P, a=0.1)
print(out)
require(geometry)
trimesh(out$T, out$P)

with(out, segments(P[S[,1], 1], P[S[,1], 2],
                   P[S[,2], 1], P[S[,2], 2], col="red"))
