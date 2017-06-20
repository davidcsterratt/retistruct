context("Checking reconstruction")
test_that("Retina reconstructs correctly", {
  dataset <- file.path(system.file(package = "retistruct"), "extdata", "GMB530/R-CONTRA")
  r <- retistruct.read.dataset(dataset)
  ## Load the human annotation of tears
  r <- retistruct.read.markup(r)
  ## Reconstruct
  r <- retistruct.reconstruct(r)

  ## Save as matlab
  r$dataset <- tempdir()
  retistruct.export.matlab(r)

  ## These values obtained from Retistruct 0.5.10
  expect_equal(r$E.l, 0.001450761, tolerance=0.00001)
  expect_equal(r$A.tot, 29191727)
  ## The Windows and Linux results for EOD are within 0.1 of each other. 
  ## It is not clear why, though it seems to occur during the minimisation
  ## using BFGS.
  expect_equal(r$EOD, 4.517927, tolerance=0.1) 
})
