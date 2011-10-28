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

## Project the test data to a polar representation, and plot it
r <- cbind(cos(phi)*cos(lambda), cos(phi)*sin(lambda))
plot(r)

## Bind the lattitudes and longitudes together to create the matrix of
## data points, with one data point on each row
mu <- cbind(phi, lambda)                

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

## Now find the value of sigma that optimises logP
opt <- optimise(function(sigma) {logP(mu, sigma)}, interval=c(0.1, 5),
                maximum=TRUE)
sigma <- opt$maximum

## Now we've found sigma, let's try to estimate and display the
## density over our polar representation of the data points

## First create a grid in Cartesian coordinates
xlim <- range(r[,1])
ylim <- range(r[,2])
xs <- seq(xlim[1], xlim[2], len=100)
ys <- seq(ylim[1], ylim[2], len=101)

## Create grid
gxs <- outer(ys*0, xs, "+")
gys <- outer(ys, xs*0, "+")
             
## Now convert it to polar coordinates
gphi <- sqrt(gxs^2 + gys^2)
glambda <- atan2(gys, gxs)




## Make space for output
kg <- gxs * 0
for (i in 1:nrow(gphi)) {
  for (j in 1:ncol(gphi)) {
    gmu <- c(gphi[i, j], glambda[i, j])
    kg[i, j] <- K(gmu, mu, sigma)
  }
}

image(xs, ys, t(kg))
points(r)

rg <- cbind(as.vector(cos(gphi)*cos(glambda)),
            as.vector(cos(gphi)*sin(glambda)))
plot(rg)
