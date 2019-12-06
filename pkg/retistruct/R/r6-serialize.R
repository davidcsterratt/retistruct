## Unfortunately we seem to have to create test classes here rather
## than in tests/testthat/test-r6-serialize.R

TestClass <-
  R6::R6Class("TestClass",
              public = list(
                a=NULL,
                b=NA,
                d="hello",
                e=45,
                f=1:10,
                g=list("a", 5, 1:9),
                h=list(a="a", b=5, c=1:9),
                i=list(a=NULL)
              ))

TestClass2 <-
  R6::R6Class("TestClass2",
              public = list(
                a=NULL,
                b=NA,
                d=NULL,
                e=45,
                f=1:10,
                g=list("a", 5, 1:9),
                initialize=function() {
                  self$d = TestClass$new()
                }
              ))

TestClass3 <-
  R6::R6Class("TestClass3",
              public = list(
                a=NULL,
                b=NA,
                d=NULL,
                e=45,
                f=1:10,
                g=list("a", 5, 1:9),
                initialize=function() {
                  self$d = list(TestClass$new(),
                           TestClass$new())
                }
              ))

TestClass4 <- 
  R6::R6Class("TestClass4",
              public = list(
                i = function() {}
              ))

TestClass5 <- 
  R6::R6Class("TestClass5",
              public = list(
                j = list(a=function() {}, b=7),
                k = list(function() {}, 7)
              ))

##' Convert an R6 object into a list, ignoring functions and
##' environments
##' @param r R6 object or list
##' @param path root of the path to the list - no need to supply. Not
##'   used but could be developed for pretty-printing
##' @param envs list of environments already encounted - do not set
##' @return List with structure mirroring the R6 object.
##' @author David Sterratt
R6_to_list <- function(r, path="", envs=list()) {
  ignore <- c("function", "environment", "NULL")
  if ("R6" %in% class(r)) {
    if (any(sapply(envs, function(x) identical(x, r$.__enclos_env__)))) {
      stop(path, " has already been serialised. Make sure there are no recursive links in the object hierarchy")
    }
    out <- list()
    attr(out, "class")  <- class(r)
    for (var in names(r)) {
      if (!any(class(r[[var]]) %in% ignore)) {
        out[[var]] <- R6_to_list(r[[var]],
                                 path=paste0(path, var, sep="$"),
                                 envs=c(envs, r$.__enclos_env__))
      }
    }
    if (is.null(out)) {
      return(list(NULL))
    }
    return(out)
  }
  if ("list" %in% class(r)) {
    out <- list()
    if (length(r) > 0) {
      length(out) <- length(r)
      for (i in 1:length(r)) {
        if (!any(class((r[[i]])) %in% ignore)) {
          out[[i]] <- R6_to_list(r[[i]],
                           path=paste0(path, i, sep="$"),
                           envs=envs)
        }
      }
      if (length(out) != length(r)) {
        stop("Length mismatch: Input r has length ", length(r), " ; Output has length ", length(out))
      }
      if (!is.null(names(r))) {
        names(out) <- names(r)
      }
    }
    return(out)
  }
  return(r)
}

##' Convert an list created by R6_to_list() into an R6 object. 
##' @param l list created by R6_to_list()
##' @return R6 object or list list
##' @author David Sterratt
list_to_R6 <- function(l) {
  cl <- attr(l, "class")
  if ("R6" %in% cl) {
    str <-  paste0("r <- ", cl[1], "$new()")
    eval(parse(text=str))
    for (var in names(l)) {
      str <- paste0("r$", var, " <- list_to_R6(l$", var, ")")
      eval(parse(text=str))
    }
    return(r)
  }
  if (is.list(l)) {
    r <- list()
    if (length(l) > 0) {
      length(r) <- length(l)
      for (i in 1:length(l)) {
        if (is.null(l[[i]])) {
          r[i] <- list(NULL)
        } else {
          r[[i]] <- list_to_R6(l[[i]])
        }
      }
      if (!is.null(names(l))) {
        names(r) <- names(l)
      }
    }
    class(r) <- cl
    return(r)
  }
  return(l)
}
