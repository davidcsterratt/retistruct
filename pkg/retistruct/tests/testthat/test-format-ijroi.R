context("Format IJROI")
test_that("IJROI format with image is read correctly", {
  dataset <- file.path(system.file(package = "retistruct"), "extdata", "smi32")
  r <- expect_warning(retistruct.read.dataset(dataset), "Scale file \"scale.csv\" does not exist. Scale bar will not be set.")
  ## Test that points are read in correctly
  P <- as.matrix(read.csv(file.path(dataset, "P.csv")))
  expect_equal(r$getPointsXY(), P)
  ## Test that optic disc is read in correctly
  od <- rbind(c(X=354, Y=336),
              c(354, 347),
              c(359, 352),
              c(364, 347),
              c(364, 338),
              c(360, 333))
  expect_equal(r$getFeatureSet("LandmarkSet")$getFeature("OD"), od)
})

test_that("IJROI format with data points and data counts is read correctly", {
    dataset <- file.path(system.file(package = "retistruct"), "extdata", "ijroi1")
    r <- expect_warning(retistruct.read.dataset(dataset), "Scale file \"scale.csv\" does not exist. Scale bar will not be set.")
    ## Test that points are read in correctly
    P <- as.matrix(read.csv(file.path(dataset, "P.csv")))
    expect_equal(r$getPointsXY(), P)

    testcounts <- rbind(c(X=100,Y=200,C=4),
                        c(110,200,4),
                        c(120,200,4),
                        c(100,210,4),
                        c(110,210,7),
                        c(120,210,4),
                        c(100,220,4),
                        c(110,220,4),
                        c(120,220,4))

    ## Test that data counts are read correctly
    expect_equal(r$getFeatureSet("CountSet")$getFeature("testcounts"), testcounts)


})
