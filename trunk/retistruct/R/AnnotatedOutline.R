AnnotatedOutline <- function(o, V0, VB, VF, phi0){
  a <- o
  class(a) <- c("annotatedOutline", class(o))
  a$V0 <- V0
  a$VB <- VB
  a$VF <- VF
  a$phi0 <- phi0
  return(a)
}
