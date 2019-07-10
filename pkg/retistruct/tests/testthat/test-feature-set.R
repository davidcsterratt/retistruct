context("FeatureSet")
test_that("FeatureSets work correctly", {
  ## This should work 
  Ds <- list(a = cbind(X=1:10, Y=1:10))
  cols <- "blue"
  fs <- PointSet$new(data=Ds, cols=cols)
  expect_equal(fs$getIDs(), "a")
  expect_equal(fs$getFeature("a"), Ds[[1]])
  
  ## If matrix columnames are not X and Y, an error should be thrown
  Ds <- list(a = cbind(x=1:10, y=1:10))
  cols <- "blue"
  expect_error(PointSet$new(data=Ds, cols=cols))


  ## Even if one matrix columnames in list are not X and Y, an error
  ## should be thrown
  Ds <- list(a = cbind(x=1:10, y=1:10),
             B = cbind(X=1:10, Y=1:10))
  cols <- "blue"
  expect_error(PointSet$new(data=Ds, cols=cols))
})

test_that("FeatureSets can be added to Outlines", {
  o <- Outline$new()
  
  Ds <- list(a = cbind(X=1:10, Y=1:10))
  cols <- "blue"
  fs <- PointSet$new(data=Ds, cols=cols)

  o$addFeatureSet(fs)
  
  expect_equal(o$getIDs(), "a")
  expect_equal(o$getFeatureSetTypes(), "PointSet")
  expect_error(o$addFeatureSet(fs), "There is already a PointSet attached to this outline")
})

test_that("FeatureSets work correctly when read from file", {
  dataset <- file.path(system.file(package = "retistruct"), "extdata", "GM114-4-RC")
  a <- retistruct.read.dataset(dataset)
  a <- retistruct.read.markup(a)
  ## A PointSet a LandmarkSet and a CountSet should be returned
  expect_equal(length(a$getFeatureSets()), 3)
  expect_equal(a$getFeatureSetTypes(), c("PointSet", "LandmarkSet", "CountSet"))
  expect_true(inherits(a$getFeatureSet("PointSet"), "PointSet"))
  expect_true(inherits(a$getFeatureSet("LandmarkSet"), "LandmarkSet"))
  expect_true(inherits(a$getFeatureSet("CountSet"), "CountSet"))

  r <- retistruct.reconstruct(a)
  ## FIXME - get retistruct.read.recdata() working
  ## r <- retistruct.read.recdata(a)

  ## A ReconstructedPointSet a ReconstructedLandmarkSet and a ReconstructedCountSet should be returned
  expect_equal(length(r$getFeatureSets()), 3)
  expect_equal(r$getFeatureSetTypes(), c("ReconstructedPointSet", "ReconstructedLandmarkSet", "ReconstructedCountSet"))
  expect_true(inherits(r$getFeatureSet("PointSet"), "ReconstructedPointSet"))
  expect_true(inherits(r$getFeatureSet("LandmarkSet"), "ReconstructedLandmarkSet"))
  expect_true(inherits(r$getFeatureSet("CountSet"), "ReconstructedCountSet"))

  ## Retistruct v0.5.x and earlier: Dss <- getDss(r)
  Dss <- r
  ## Retistruct v0.5.x and earlier: Sss <- getSss(r)
  Sss <- r$getFeatureSet("LandmarkSet")
  Dss.mean <- r$getFeatureSet("PointSet")$getMean()
  Dss.hullarea <- r$getFeatureSet("PointSet")$getHullarea()
  Dss.KDE <- r$getFeatureSet("PointSet")$getKDE()

  ## Images
  Ims <- r$getIms()
})
