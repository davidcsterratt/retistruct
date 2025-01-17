R.phi0 <- function(phi0, A.tot=1) {
  sqrt(A.tot/(2*pi*(sin(phi0)+1)))
}

dR.dphi0 <- function(phi0, A.tot=1) {
  (pi*cos(phi0*sqrt(A.tot)))/((2*pi*(sin(phi0)+1))^1.5)
}

phi0 <- (0:90)/180*pi

par(mfrow=c(2,1))
par(mar=c(3, 3, 0.5, 0.5))
par(mgp=c(1.5, 0.5, 0))

plot(phi0, R.phi0(phi0),type='l', ylim=c(0, 0.5),
     xlab=expression(phi[0]),
     ylab="R")
plot(phi0, dR.dphi0(phi0), type='l',
     xlab=expression(phi[0]),
     ylab=expression(paste(dR/d, phi[0])))

dev.copy2pdf(file="R-phi.pdf", width=3, height=3)




