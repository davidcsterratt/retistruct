h <- 0.1                               # Time step
epsilon <- 1E-5                         # Tolerance

M <- diag(c(1, 1))                       # Mass
gamma <- 1.5                              # Friction coeff
q0 <- c(-1, 0)                           # Initial position
v0 <- c( 0, 0)                           # Initial velocity

t <- 0                                  # Initial time
T <- 20                                  # Maximum time

w0 <- c(0, 0, 0)
w <- w0

q <- matrix(NA, 2, T/h)
v <- matrix(NA, 2, T/h)
q[,1:2] <- q0
v[,1:2] <- v0

n <- 2

## Q <- function(q, v) { 
Q.v <- function(q, v) { diag(gamma*c(-1, -1)) }
Phi <- function(q) { matrix(sum(q^2) - 1, 1, 1) }
Phi.q <- function(q) { matrix(2*q, 1, 2) }
Phi.q.q <- function(q) { diag(2, 2) }

while (n*h < T) {
  iter <- 20
  while (iter > 0) {
    q[,n+1] <- 4/3*q[,n] - 1/3*q[,n-1] + h*(8/9*v[,n] - 2/9*v[,n-1]) + 4/9*h^2*w[1:2]
    v[,n+1] <- 4/3*v[,n] - 1/3*v[,n-1] + 2/3*h*w[1:2]

    Q.v <- diag(gamma*c(-1, -1))
    Phi <- matrix(sum(q[,n+1]^2) - 1, 1, 1)
    Phi.q <- matrix(2*q[,n+1], 1, 2)
    Phi.q.q <- diag(2, 2)
    P <- 4/9*h^2*(Phi.q.q * w[3]) - 2/3*h*Q.v
    J <- rbind(cbind(M + P, t(Phi.q)),
               cbind(Phi.q, 0))
    Y <- c(M %*% w[1:2] + t(Phi.q) %*% w[3] - (c(0, -1) - gamma*v[,n+1]),
           9/4/h^2*Phi)
    Deltaw <- solve(J, Y)
    ## print(Deltaw)
    if (all(abs(Deltaw) < epsilon)) {
      break
    }
    w <- w - Deltaw
    iter <- iter - 1
  }
  if (iter == 0) {
    print("Warning: maxiumum number of iterations reached")
  }
  n <- n + 1
}

par(mfrow=c(3,1))
plot(q[1,], type="l", col="red")
lines(q[2,], col="blue")
plot(v[1,], type="l", col="red")
lines(v[2,], col="blue")

plot(atan2(q[1,], -q[2,]), type="l")
