## Test of kernel density clustering

## Metric function - compute the distance from one point (r) to a
## group of points in the rows of the matrix s
metric <- function(r, s) {
  return(central.angle(r[1], r[2], s[,1], s[,2]))
}

## Bind the lattitudes and longitudes together to create the matrix of
## data points, with one data point on each row
## mu <- cbind(phi=phi, lambda=lambda)                

## Project the test data to a polar representation, and plot it
## r <- cbind(cos(phi)*cos(lambda), cos(phi)*sin(lambda))
## plot(r)

## Kernel function - just a Gaussian
kappa <- function(x, mu, sigma) {
  return(1/sqrt(2*pi)/sigma*exp(-1/2*(metric(x, mu))^2/sigma^2))
}

## Estimate of density at x, given points at mu and a bandwidth of
## sigma
K <- function(x, mu, sigma) {
  return(1/nrow(mu)*sum(kappa(x, mu, sigma)))
}

## Estimate of the log probability of the points mu given a particular
## value of sigma
logP <- function(mu, sigma) {
  ## We get clever here, and define a funtion within a function. This
  ## is the kernel density for data point i if the density is
  ## determined by all the other points. Note that mu[-i,] means all
  ## the rows of mu apart from row i
  logKi <- function(i) {
    return(log(K(mu[i,], mu[-i,], sigma)))
  }

  ## We now pass this function to sapply, which creates a vector of
  ## the result of logKi for every value of i
  return(sum(sapply(1:nrow(mu), logKi)))
}

## Find the value of sigma that optimises logP
compute.bandwidth <- function(mu, K) {
  opt <- optimise(function(sigma) {logP(mu, sigma)}, interval=c(0.1, 5),
                  maximum=TRUE)
  return(opt$maximum)
}

## gcb <- polar.to.cart(gp)
## image(xs, ys, k)
## points(r)
## cs <- contourLines(xs, ys, k, levels=klevels)


