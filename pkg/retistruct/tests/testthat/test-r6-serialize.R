context("R6 Serialize")
test_that("R6 serialization works", {
  ## Simple test class
  tc <- TestClass$new()
  tcl <- R6_to_list(tc)
  ## a is null, so shouldn't be created
  expect_equal(tc$b, tcl$b)
  expect_equal(tc$d, tcl$d)
  expect_equal(tc$e, tcl$e)
  expect_equal(tc$f, tcl$f)
  expect_equal(tc$g, tcl$g)
  expect_equal(tc$h, tcl$h)

  tclr6 <- list_to_R6(tcl)
  expect_equal(tc$b, tclr6$b)
  expect_equal(tc$d, tclr6$d)
  expect_equal(tc$e, tclr6$e)
  expect_equal(tc$f, tclr6$f)
  expect_equal(tc$g, tclr6$g)
  expect_equal(tc$h, tclr6$h)
  expect_equal(tc$h, tclr6$h)
  expect_equal(tc, tclr6)

  ## Test class with R6 object as member
  tc2 <- TestClass2$new()
  tc2l <- R6_to_list(tc2)
  tc2lr6 <- list_to_R6(tc2l)
  expect_equal(tc2, tc2lr6)

  ## Test class with list of R6 objects as member
  tc3 <- TestClass3$new()
  tc3l <- R6_to_list(tc3)
  tc3lr6 <- list_to_R6(tc3l)
  expect_equal(tc3, tc3lr6)

  ## Test serlization
  tc3l_file <- tempfile(fileext=".RData")
  save(tc3l, file=tc3l_file)
  rm(tc3l)
  rm(tc3lr6)
  load(tc3l_file)
  tc3lr6 <- list_to_R6(tc3l)
  expect_equal(tc3, tc3lr6)

  ## Test function
  tc4 <- TestClass4$new()
  tc4l <- R6_to_list(tc4)
  tc4lr6 <- list_to_R6(tc4l)

  ## Test function in list
  tc5 <- TestClass5$new()
  tc5l <- R6_to_list(tc5)
  tc5lr6 <- list_to_R6(tc5l)


  ## Outline object
  ol <- Outline$new()
  oll <- R6_to_list(ol)
  ollr6 <- list_to_R6(oll)
  expect_equal(ol, ollr6)


})
