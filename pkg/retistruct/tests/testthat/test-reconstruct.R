context("Reconstruction")
test_that("Reconstruct SMI32 (CSV format)", {
  dataset <- file.path(tempdir(), "smi32-csv")
  file.copy(file.path(system.file(package = "retistruct"), "extdata", "smi32-csv"), tempdir(), recursive=TRUE, overwrite=TRUE)

  expect_warning(o <- retistruct.read.dataset(dataset),  "Scale bar will not be set")
  ## Load the human annotation of tears
  o <- retistruct.read.markup(o)
  ## Reconstruct
  r <- retistruct.reconstruct(o)

  ## Save as matlab
  filename <-  file.path(tempdir(), "r.mat")
  retistruct.export.matlab(r, filename)

  ## Nascent API
  expect_error(r$mapFlatToSpherical(0), "P must be matrix")
  expect_error(r$mapFlatToSpherical(c(0, 0)), "P must be matrix")
  expect_warning(r$mapFlatToSpherical(cbind(X=0, Y=0)), "1 points outwith the outline will be ignored")
  expect_equal(r$mapFlatToSpherical(cbind(X=200,Y=200)), cbind(phi=-0.1686439, lambda=-1.559814), tol=0.001)

  ## Serialisation
  r0 <- r$clone()
  retistruct.save.recdata(r)
  rm(r)
  r1 <- retistruct.read.recdata(o)
  ## The test below won't work if we plot a projection projection
  ## before, projection() results values being put into the private
  ## member variable r0$ims. However, this is variable is not stored.
  expect_equal(r0, r1)

  ## Test projection after the above test
  png(file=file.path(tempdir(), "smi32-projection.png"), width=800)
  projection(r0)
  dev.off()
})

test_that("Reconstruct GMB530/R-CONTRA (IDT format)", {
  dataset <- file.path(tempdir(), "GMB530", "R-CONTRA")
  if (!dir.exists(file.path(tempdir(), "GMB530"))) {
    dir.create(file.path(tempdir(), "GMB530"))
  }
  file.copy(file.path(system.file(package = "retistruct"), "extdata", "GMB530/R-CONTRA"), file.path(tempdir(), "GMB530"), recursive=TRUE, overwrite=TRUE)

  o <- retistruct.read.dataset(dataset)
  ## Load the human annotation of tears
  o <- retistruct.read.markup(o)
  ## Reconstruct
  r <- retistruct.reconstruct(o)

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


  ## Test serialisation - note we do this before saving to matlab,
  ## since KDE is cached and it messes up the comparsion
  r0 <- r$clone()
  retistruct.save.recdata(r)
  rm(r)
  r1 <- retistruct.read.recdata(o)
  expect_equal(r0$featureSets[[1]]$KDE, r1$featureSets[[1]]$KDE)
  expect_equal(r0$featureSets[[1]], r1$featureSets[[1]])
  expect_equal(r0$featureSets, r1$featureSets)
  expect_equal(r0, r1)

  ## Test getting features
  r0$getFeatureSet("PointSet")
  r0$getFeatureSet("LandmarkSet")
  cs <- r0$getFeatureSet("CountSet")
  expect_equal(colnames(cs$data[[1]]), c("phi", "lambda", "C"))

  r0$ol$DVflip  <- TRUE
  r0$getFeatureSet("PointSet")
  r0$getFeatureSet("LandmarkSet")
  cs <- r0$getFeatureSet("CountSet")
  expect_equal(colnames(cs$data[[1]]), c("phi", "lambda", "C"))

  ## Save as matlab
  filename <-  file.path(tempdir(), "r.mat")
  retistruct.export.matlab(r0, filename)
})

test_that("Serialisation works with a particular example", {
  dataset <- file.path(tempdir(), "GM509", "R-CONTRA")
  if (!dir.exists(file.path(tempdir(), "GM509"))) {
    dir.create(file.path(tempdir(), "GM509"))
  }
  file.copy(file.path(system.file(package = "retistruct"), "extdata", "GM509/R-CONTRA"), file.path(tempdir(), "GM509"), recursive=TRUE, overwrite=TRUE)

  a <- retistruct.read.dataset(dataset, report=FALSE)
  a <- retistruct.read.markup(a, error=message)
  r <- retistruct.reconstruct(a) ## plot.3d=getOption("show.sphere")
  retistruct.save.recdata(r)
  r1 <- retistruct.read.recdata(a)

  png(file=file.path(tempdir(), "flat.png"), width=800)
  flatplot(a)
  dev.off()

  png(file=file.path(tempdir(), "proj.png"), width=800)
  projection(r)
  dev.off()

  retistruct.save.recdata(r)
  r2 <- retistruct.read.recdata(a)

})
