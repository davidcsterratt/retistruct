## NSTIFF - numerical integration of system with constraints
##
## Based on NSTIFF algorithm presented in Negrut, Jay & Khude
## (2009). Journal of ##Computational and Nonlinear Dynamics 4
## 
nstiff <- function(M, Q, Q.v, Phi, Phi.q, Phi.q.q,
                   q0, v0,
                   h=0.02,              # Time step
                   epslion=1E-5,        # Tolerance
                   T=20,                # Maximum time
                   kmax=10 # Maximum number of iterations per time step
                   ) {  
  
  n <- length(q0)                  # Number of generalised coordinates
  m <- nrow(Phi(q0))                  # Number of Lagrange multipliers

  I <- T/h + 1                          # Number of steps
  t <- seq(from=0, by=h, len=I)         # Time store
  q <- matrix(NA, N, I)                 # Position store
  v <- matrix(NA, N, I)                 # Velocity store
  q[,1:2] <- q0                         # Initial values
  v[,1:2] <- v0                         # Initial values
  w <- rep(0, n + m)                    # Vector of dvdt and lambda
  i <- 2                                # Time index

  while (i < I) {
    k <- 0
    while (k < kmax) {
      q[,i+1] <- 4/3*q[,i] - 1/3*q[,i-1] + h*(8/9*v[,i] - 2/9*v[,i-1]) + 4/9*h^2*w[1:n]
      v[,i+1] <- 4/3*v[,i] - 1/3*v[,i-1] + 2/3*h*w[1:n]

      P <- 4/9*h^2*(Phi.q.q(q[,i+1]) * w[n+1:m]) - 2/3*h*Q.v(q[,i+1], v[,i+1])
      J <- rbind(cbind(M + P,          t(Phi.q(q[,i+1]))),
                 cbind(Phi.q(q[,i+1]), 0))

      Y <- c(M %*% w[1:n] + t(Phi.q(q[,i+1])) %*% w[n+1:m] - Q(q[,i+1], v[,i+1]),
             9/4/h^2*Phi(q[,i+1]))

      Deltaw <- solve(J, Y)
      if (all(abs(Deltaw) < epsilon)) {
        break
      }
      w <- w - Deltaw
      k <- k + 1
    }
    if (k == kmax) {
      print("Warning: maxiumum number of iterations reached")
    }
    i <- i + 1
  }
  return(list(t=t, q=q, v=v))
}

## Demo of nstiff()
nstiff.demo <- function() {
  ## The system
  gamma <- 1.5                              # Friction coeff

  out <- nstiff(M = diag(c(1, 1)),                       # Mass
                Q = function(q, v)    { c(0, -1) - gamma*v },
                Q.v = function(q, v)  { diag(gamma*c(-1, -1)) },
                Phi = function(q)     { matrix(sum(q^2) - 1, 1, 1) },
                Phi.q = function(q)   { matrix(2*q, 1, 2) },
                Phi.q.q = function(q) { diag(2, 2) },
                q0 <- c(-1, 0),                  # Initial position
                v0 <- c( 0, 0)                   # Initial velocity
                )

  with(out, {
    par(mfrow=c(3,1))
    plot(t, q[1,], type="l", col="red")
    lines(t, q[2,], col="blue")
    plot(t, v[1,], type="l", col="red")
    lines(t, v[2,], col="blue")
    plot(t, atan2(q[1,], -q[2,]), type="l")
  })
}
