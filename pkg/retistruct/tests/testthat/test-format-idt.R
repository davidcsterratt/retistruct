context("Format IDT")
test_that("IDT format is read correctly", {
  dataset <- file.path(system.file(package = "retistruct"), "extdata", "GMB530/R-CONTRA")
  r <- retistruct.read.dataset(dataset)
  ## Test that points are read in correctly
  P <- as.matrix(read.csv(file.path(dataset, "P.csv"), row.names=NULL))
  colnames(P) <- c("X", "Y")
  expect_that(r$getPointsXY(), equals(P))
  r <- retistruct.read.markup(r)
  ## Test that optic disc is read in correctly
  od <- rbind(c(X=1090, Y=1319),
              c(1230, 1359),
              c(1250, 1499),
              c(1160, 1559),
              c(1050, 1539),
              c( 980, 1519),
              c(1030, 1409),
              c(1070, 1339),
              c(1120, 1339))

  colnames(od) <- c("X", "Y")
  expect_that(r$getFeatureSet("LandmarkSet")$getFeature("OD"), equals(od))
  expect_that(r$scale, equals(0.7))

  expect_that(r$getFeatureSet("PointSet")$getIDs(), equals(c("green", "red", "double")))
  expect_true(is.matrix(r$getFeatureSet("PointSet")$getFeature("green")))
  expect_that(dim(r$getFeatureSet("PointSet")$getFeature("green")), equals(c(169, 2)))
  expect_that(colnames(r$getFeatureSet("PointSet")$getFeature("green")), equals(c("X", "Y")))

  expect_that(r$getFeatureSet("CountSet")$getIDs(), equals(c("green", "red", "double")))
  expect_true(is.matrix(r$getFeatureSet("CountSet")$getFeature("green")))
  expect_that(ncol(r$getFeatureSet("CountSet")$getFeature("green")), equals(3))
  expect_that(colnames(r$getFeatureSet("CountSet")$getFeature("green")), equals(c("X", "Y", "C")))
})
