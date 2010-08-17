retistruct.cli <- function(dataset) {
  retistruct.read.dataset(dataset)
  retistruct.reconstruct()
  retistruct.save()
}
