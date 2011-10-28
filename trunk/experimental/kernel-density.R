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

cart.to.polar <- function(r) {
  rho <- sqrt(r[,"x"]^2 + r[,"y"]^2)
  phi <- asin(rho - 1)
  lambda <- atan2(r[,"y"], r[,"x"])
  return(cbind(phi=phi, lambda=lambda))
}

polar.to.cart <- function(r) {
  x <- (1 + sin(r[,"phi"]))*cos(r[,"lambda"])
  y <- (1 + sin(r[,"phi"]))*sin(r[,"lambda"])
  return(cbind(x=x, y=y))
}

## Generate some test data
N <- 100                                # Number of points
lambda <- 2*pi*runif(N)                 # Latitudes
phi <- -pi/2 + 0.2*(rnorm(N))^2         # Longitudes

## Bind the lattitudes and longitudes together to create the matrix of
## data points, with one data point on each row
mu <- cbind(phi=phi, lambda=lambda)                

## Project the test data to a polar representation, and plot it
r <- cbind(cos(phi)*cos(lambda), cos(phi)*sin(lambda))
plot(r)


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
gxs <- outer(xs*0, ys, "+")
gys <- outer(xs, ys*0, "+")

## gxs and gys are both 101 by 100 matrixes We now combine both
## matrices as a 101*100 by 2 matrix. The conversion as.vector() goes
## down the columns of the matrices gxs and gys
gc <- cbind(x=as.vector(gxs), y=as.vector(gys))

## Now convert the cartesian coordinates to polar coordinates
gp <- cart.to.polar(gc)

## gcb <- polar.to.cart(gp)

## Make space for the kernel density estimates
kg <- rep(0, nrow(gp))
for (i in 1:nrow(gp)) {
  kg[i] <- K(gp[i,], mu, sigma)
}

## Put the estimates back into a matrix. The matrix is filled up
## column-wise, so the matrix elements should match the elements of
## gxs and gys
kg <- matrix(kg, 100, 101)
image(xs, ys, kg)
points(r)


