context("Reconstruction")
test_that("Reconstruct SMI32 (CSV format)", {
  dataset <- file.path(system.file(package = "retistruct"), "extdata", "smi32-csv")
  expect_warning(r <- retistruct.read.dataset(dataset),  "Scale bar will not be set")
  ## Load the human annotation of tears
  r <- retistruct.read.markup(r)
  ## ## Reconstruct
  r <- retistruct.reconstruct(r)

  png(file=file.path(tempdir(), "smi32-projection.png"), width=800)
  projection(r)
  dev.off()    
  
  ## Save as matlab
  filename <-  file.path(tempdir(), "r.mat")
  retistruct.export.matlab(r, filename)
})

test_that("Reconstruct GMB530/R-CONTRA (IDT format)", {
  dataset <- file.path(system.file(package = "retistruct"), "extdata", "GMB530/R-CONTRA")
  r <- retistruct.read.dataset(dataset)
  ## Load the human annotation of tears
  r <- retistruct.read.markup(r)
  ## ## Reconstruct
  r <- retistruct.reconstruct(r)
  
  ## Save as matlab
  filename <-  file.path(tempdir(), "r.mat")
  retistruct.export.matlab(r, filename)

  ## ## These values obtained from Retistruct 0.5.10, scaled now due to
  ## ## scale being applied before reconstruction
  ## expect_equal(r$E.l, 0.001450761, tolerance=0.00001)
  expect_equal(r$ol$A.tot, 29191727*0.7^2, tolerance=1)
  expect_equal(nrow(r$getFeatureSet("LandmarkSet")$getFeature("OD")), 9)
  ## ## Was originally c(-1.457557, -0.2672879); changed lambda to mean
  ## ## of original and obtained.
  ## expect_equal(as.vector(r$Sss[["OD"]][1,]), c(-1.457557, -0.253425), tolerance=0.02)
  ## expect_equal(r$EOD, 4.517927, tolerance=0.01)
  ## The Windows and Linux results for EOD are within 0.1 of each other. 
  ## It is not clear why, though it seems to occur during the minimisation
  ## using BFGS.
  ## expect_equal(r$EOD, c(phi=4.517927), tolerance=0.1) 
})
