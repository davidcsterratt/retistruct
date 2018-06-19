context("TriangulatedOutline")
test_that("Triangulations have the right number of points", {
  dataset <- file.path(system.file(package = "retistruct"), "extdata", "issue21")
  a <- expect_warning(retistruct.read.dataset(dataset), "Scale file \"scale.csv\" does not exist. Scale bar will not be set.")
  o <- retistruct.read.markup(a)
  expect_equal(nrow(o$P), length(a$gf))
  o$lambda0 <- 0
  n <- 500
  t <- TriangulatedOutline(o, n=n)
  expect_equal(nrow(t$P), length(t$gf))
  expect_equal(nrow(t$P), length(t$gb))
  ## expect_equal(nrow(t$P), length(t$hf))
  ## expect_equal(nrow(t$P), length(t$hb))
  
})
