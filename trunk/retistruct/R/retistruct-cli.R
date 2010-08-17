retistruct.cli <- function(dataset, cpu.time.limit=Inf) {
  setTimeLimit(cpu=cpu.time.limit)
  dataset <<- dataset
  out <- try(retistruct.cli.process())
  mess <- geterrmessage()
  if (inherits(out, "try-error")) {
    if (grepl("reached CPU time limit", mess)) {
      quit(status=1)
    } else {
      ## Unknown error
      quit(status=2)
    }
  }
  ## Success
  quit(status=0)
}

retistruct.cli.process <- function() {
  retistruct.read.dataset()
  retistruct.reconstruct()
  retistruct.save()
}
