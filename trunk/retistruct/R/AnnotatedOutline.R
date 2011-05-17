AnnotatedOutline <- function(o){
  a <- o
  class(a) <- c("annotatedOutline", class(o))
  a$V0 <- c()
  a$VB <- c()
  a$VF <- c()
  a$phi0 <- 0
  return(a)
}

##' Label a set of three unlabelled points supposed to refer to the
##' apex and vertcies of a cut and tear with the V0 (Apex), VF
##' (forward vertex) and VB (backward vertex) labels.
##'
##' @title Label three outline point indicies as apicies and vertices of tear
##' @param m the vector of three indicies
##' @param o Outline object
##' @return Vector of indicies labelled with V0, VF and VB
##' @author David Sterratt
labelTearPoints <- function(o, m) {
  with(o, {
    ## Each row of this matrix is a permutation of the markers
    p <- rbind(c(1, 2, 3),
               c(1, 3, 2),
               c(2, 1, 3),
               c(2, 3, 1),
               c(3, 1, 2),
               c(3, 2, 1))

    ## For each permuation of V0, VF, VB, measure the sum of length in
    ## the forwards direction from V0 to VF and in the backwards
    ## direction from V0 to VB. The permuation with the minimum distance
    ## is the correct one.
    tplmin <- Inf                      # The minimum path length
    h <- 1:nrow(P)                     # identity correspondence mapping
                                        # used for measuring distances
                                        # (this effectively ignores
                                        # sub-tears, but this doesn't
                                        # matter)
    for (i in 1:nrow(p)) {
      V0 <- m[p[i,1]]
      VF <- m[p[i,2]]
      VB <- m[p[i,3]]
      tpl <- path.length(V0, VF, gf, h, P) + path.length(V0, VB, gb, h, P)
      if (tpl < tplmin) {
        M <- m[p[i,]]
        tplmin <- tpl
      }
    }
    names(M) <- c("V0", "VF", "VB")
    return(M)})
}

addTear <- function(o, id) {
  M <- labelTearPoints(r, id)
  o$V0 <- c(o$V0, M["V0"])
  o$VF <- c(o$VF, M["VF"])
  o$VB <- c(o$VB, M["VB"])
  return(o)
}
