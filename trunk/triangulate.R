connections2triangulation <- function(C) {
  T <- matrix(NA, 0, 3)
  ## Convert a matrix containing links into a triangulation
  for (i in 1:(dim(C)[1])) {
    P1 <- C[i,1]
    P2 <- C[i,2]
    P1n <- setdiff(c(Cu[Cu[,2]==P1, 1], Cu[Cu[,1]==P1, 2]), P2)
    P2n <- setdiff(c(Cu[Cu[,2]==P2, 1], Cu[Cu[,1]==P2, 2]), P1)
##    print(P1n)
##    print(P2n)
    for (P3 in intersect(P1n, P2n)) {
      T.row <- sort(c(P1, P2, P3))
      if (length(unique(T.row)) != 3) {
        print("Error")
        print(T.row)
      } else {
        T <- rbind(T, T.row)
      }
    }
    if (dim(T)[1]==841) {
      print(paste(P1, P2, P3))
      print(paste(P1n))
      print(paste(P2n))
      print(i)
    }

  }
  return(T)
}
