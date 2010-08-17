retistruct.cli <- function(dataset) {
  dataset <<- dataset
  retistruct.read.dataset()
  retistruct.reconstruct()
  retistruct.save()
}
