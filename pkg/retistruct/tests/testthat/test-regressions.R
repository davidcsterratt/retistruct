context("regressions")
test_that("No regression on Issue #21", {
  dataset <- file.path(system.file(package = "retistruct"), "extdata", "issue21")
  a <- expect_warning(retistruct.read.dataset(dataset), "Scale file \"scale.csv\" does not exist. Scale bar will not be set.")
  o <- retistruct.read.markup(a)
  expect_equal(nrow(o$P), length(a$gf))
  o$lambda0 <- 0
  n <- 500
  t <- TriangulatedOutline(o, n=n)
  s <- StitchedOutline(t)
  expect_equal(nrow(s$P), length(s$gf))
  r <- TriangulatedOutline(s, n=n,
                           suppress.external.steiner=TRUE)
  expect_equal(nrow(r$P), length(r$gf))
  r <- mergePointsEdges(r)
  expect_equal(nrow(r$P), length(r$gf))
  r <- projectToSphere(r)
})
