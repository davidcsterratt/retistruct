## Test of kernel density clustering

## Formula to comput the distance between two points on the sphere
## with coordinates (phi1, lambda1), (phi2, lambda2)
central.angle <- function(phi1, lambda1, phi2, lambda2) {
  return(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2)))
}

## Metric function - compute the distance from one point (r) to a
## group of points in the rows of the matrix s
metric <- function(r, s) {
  return(central.angle(r[1], r[2], s[,1], s[,2]))
}

## Generate some test data

N <- 200                                # Number of points

lambda <- 2*pi*runif(N)                 # Latitudes
phi <- -pi/2 + 0.2*(rnorm(N))^2         # Longitudes

## Bind the lattitudes and longitudes together
mu <- cbind(phi, lambda)                

## Project the test data to a polar representation, and plot it
r <- cbind(cos(phi)*cos(lambda), cos(phi)*sin(lambda))
plot(r)

d <- matrix(NA, N, N)
## Compute distance matrix
for (i in 1:N) {
  for (j in 1:N) {
    d[i, j] <- ifelse(i !=j,
                      metric(mu[i], mu[j]),
                      0)

  }
}

kappa <- function(x, mu, sigma) {
  return(1/sqrt(2*pi)/sigma*exp(-1/2*(metric(x, mu))^2/sigma^2))
}

K <- function(x, mu, sigma) {
  return(1/nrow(mu)*sum(kappa(x, mu, sigma)))
}

Ki <- function(mu, sigma, i) {
  return(K(mu[i,], mu[-i,], sigma))
}

logP <- function(mu, sigma) {
  return(sum(sapply(1:nrow(mu), function(i) {log(Ki(mu, sigma, i))})))
}


opt <- optimise(function(sigma) {logP(mu, sigma)}, interval=c(0.1, 5), maximum=TRUE)

sigma <- opt$maximum

