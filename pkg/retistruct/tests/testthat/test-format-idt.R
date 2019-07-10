context("Format IDT")
test_that("IDT format is read correctly", {
  dataset <- file.path(system.file(package = "retistruct"), "extdata", "GMB530/R-CONTRA")
  r <- retistruct.read.dataset(dataset)
  ## Test that points are read in correctly 
  P <- as.matrix(read.csv(file.path(dataset, "P.csv"), row.names=NULL))
  colnames(P) <- c("X", "Y")
  expect_that(r$getPoints(), equals(P))
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
})
