AnnotatedDataset <- function(a, iN, iD, iOD, DVflip, side){
  class(a) <- c("annotatedDataset", class(a))
  a$iN <- iN
  a$iD <- iD
  a$iOD <- iOD
  a$DVflip <- DVflip
  a$side <- side
  return(a)
}

plot.flat.annotatedDataset <- function(a, axt="n", ylim=NULL, ...) {
  if (a$DVflip) {
    if (is.null(ylim)) {
      ylim <- c(max(a$P[,2]), min(a$P[,2]))
    } else {
      ylim <- sort(ylim, TRUE)
    }
  }
  plot.flat.dataset(a, axt=axt, ylim=ylim, ...)
}
