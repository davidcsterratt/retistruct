
E <- 1                                  # Modulus
nu <- 0                 # Poission's ratio - can be in range (0, 0.5)

## Start of with one triangle
x <- c(0, 1, 0.5)
y <- c(0, 0, 1) 

## Area of triangle
Delta <- 0.5 * det(cbind(1, x, y))

## Construct B matrix
a <- c()
b <- c()
c <- c()

a[1] <- x[2]*y[3] - y[2]*x[3]
a[2] <- x[3]*y[1] - y[3]*x[1]
a[3] <- x[1]*y[2] - y[1]*x[2]
b[1] <- y[2] - y[3]
b[2] <- y[3] - y[1]
b[3] <- y[1] - y[2]
c[1] <- x[3] - x[2]
c[2] <- x[1] - x[3]
c[3] <- x[2] - x[1]

B1 <- 1/(2*Delta) * rbind(c(b[1], 0),
                          c(0,    c[1]),
                          c(c[1], b[1]))
B2 <- 1/(2*Delta) * rbind(c(b[2], 0),
                          c(0,    c[2]),
                          c(c[2], b[2]))
B3 <- 1/(2*Delta) * rbind(c(b[3], 0),
                          c(0,    c[3]),
                          c(c[3], b[3]))

B <- cbind(B1, B2, B3)

## Elasticity matrix for plane stress in an isotropic material
D <- E * (1-nu^2) * rbind(c(1 , nu, 0       ),
                          c(nu,  1, 0       ),
                          c(0 ,  0, (1-nu)/2))


print(t(B1) %*% D %*% B1)

K <- t(B) %*% D %*% B

rotate <- function(theta) {
  return(rbind(c(cos(theta), -sin(theta)),
               c(sin(theta),  cos(theta))))
}

rotate3 <- function(theta) {
  M <- matrix(0, 6, 6)
  M[1:2, 1:2] <- rotate(theta)
  M[3:4, 3:4] <- rotate(theta)
  M[5:6, 5:6] <- rotate(theta)
  return(M)
}

shear <- function(k) {
  return(rbind(c(1, k),
               c(0, 1)))
}

shear3 <- function(k) {
  M <- matrix(0, 6, 6)
  M[1:2, 1:2] <- shear(k)
  M[3:4, 3:4] <- shear(k)
  M[5:6, 5:6] <- shear(k)
  return(M)
}


## Original position
v <- matrix((rbind(x, y)), 6, 1)

## Translation
u <- v + 1
print("Translation")
a <- u - v                              # Displacement
print("Forces")
print(K %*% a)
print("Strain")
print(B %*% a)

## Expansion
print("Expansion")
u <- 2 * v
a <- u - v
plot( u[c(1,3,5,1)], u[c(2,4,6,2)], type='l', col="red")
lines(v[c(1,3,5,1)], v[c(2,4,6,2)])
print("Forces")
print(K %*% a)
print("Strain")
print(B %*% a)

## Rotation
print("Rotation")
u <- rotate3(pi/2) %*% v
a <- u - v
plot( u[c(1,3,5,1)], u[c(2,4,6,2)], type='l', col="red", xlim=c(-1, 1), ylim=c(0,2))
lines(v[c(1,3,5,1)], v[c(2,4,6,2)])
print("Forces")
print(K %*% a)
print("Strain")
print(B %*% a)

## Shear
print("Shear")
u <- shear3(0.5) %*% v
a <- u - v
plot( u[c(1,3,5,1)], u[c(2,4,6,2)], type='l', col="red", xlim=c(-1, 1), ylim=c(0,2))
lines(v[c(1,3,5,1)], v[c(2,4,6,2)])
print("Forces")
print(K %*% a)
print("Strain")
print(B %*% a)
