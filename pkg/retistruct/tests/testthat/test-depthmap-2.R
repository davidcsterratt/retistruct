context("Real-world depthmap 2")

test_that("Real-world depthmap 2 works", {
  if (!is.null(getOption("retistruct.test.full"))) {
    if (length(system.file(package = "retistructdemos"))) {
      dataset <- file.path(tempdir(), "depthmap-PJN19-006")
      file.copy(file.path(system.file(package = "retistructdemos"), "extdata", "depthmap-PJN19-006"), tempdir(), recursive=TRUE, overwrite=TRUE)
      expect_warning(o <- retistruct.read.dataset(dataset), "^The background value")
      o <- retistruct.read.markup(o)
      r <- retistruct.reconstruct(o)
    }
  }
})
