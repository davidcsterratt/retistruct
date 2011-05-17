AnnotatedOutline <- function(o){
  a <- o
  class(a) <- c("annotatedOutline", class(o))
  a$V0 <- c()
  a$VB <- c()
  a$VF <- c()
  a$phi0 <- 0
  return(a)
}
