## Modulus function, that returns
## i , if i <= N
## i - N, if i > N
Mod <- function(i, N) {
  return((i - 1) %% N + 1)
}

## Flip matrix between up and down
flipud <- function(M) {
  return(M[nrow(M):1,])
}
