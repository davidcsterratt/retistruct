context("Real-world depthmap 1")
test_that("Real-world depthmap 1 works", {
  if (!is.null(getOption("retistruct.test.full"))) {
    dataset <- file.path(tempdir(), "depthmap-2020-08-21")
      file.copy(file.path(system.file(package = "retistructdemos"), "extdata", "depthmap-2020-08-21"), tempdir(), recursive=TRUE, overwrite=TRUE)

    expect_warning(o <- retistruct.read.dataset(dataset),  "Scale bar will not be set")
    ## Load the human annotation of tears
    o <- retistruct.read.markup(o)

    r <- ReconstructedOutline$new()
    r$loadOutline(o)

    expect_true(r$ol$isStitched())

    ## expect_equal(r$ol$hf[44], 409)
    ## rs  <- r$ol$getRimSet()
    ## i44  <- which(rs == 44)
    ## expect_equal(rs[i44 + 1], r$ol$hf[44])

    ## r$mergePointsEdges()
    ## expect_equal(length(r$Lt), nrow(r$Cut))
    ## expect_true(!any(is.na(r$Lt)))
    ## r$projectToSphere()
    ## expect_true(all(!is.na(r$phi)))
    ## expect_true(all(!is.na(r$lambda)))
    ## expect_true(is.numeric(r$R))
    ## r$getStrains()

    ## Reconstruct
    ## r <- retistruct.reconstruct(o)
  }
})
