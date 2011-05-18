AnnotatedDataset <- function(a){
  class(a) <- c("annotatedDataset", class(a))
  a$iN <- NA
  a$iD <- NA
  a$iOD <- NA
  a$DVflip <- FALSE
  a$side <- "Right"
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
  NextMethod()
  ## plot.flat.dataset(a, axt=axt, ylim=ylim, ...)
}
