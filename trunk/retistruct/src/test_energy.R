if (is.loaded("energy")) dyn.unload("energy.so")
dyn.load("energy.so")

p <- 1:17
Cu <- matrix(0, 10, 2)
Cu <- matrix(1:3, 20, 2)
L <- 1:20
B <- NULL
T <- NULL
A <- NULL
R <- 1
Rset <- c(2,3)
i0 <- 10
phi0 <- 30
Nphi <- 8
E0A <- 0
kA <- 1
N <- 10
verbose <- FALSE

for (i in 1:10000) {
  .Call("spheristruct_E",
      p, Cu, C, L, B, T,
      as.double(A), as.double(R),
      as.integer(Rset), as.integer(i0), phi0, Nphi, 
      E0A, kA, as.integer(N), verbose)
}
print("C finished")

for (i in 1:10000) {
E(p, Cu, C, L, B, T,
      A, R,
      Rset, i0, phi0, Nphi, 
      E0A, kA, N, verbose)
}
print("R finished")

for (i in 1:1000) 
E(p, f$m$Cut, f$m$Ct, f$m$Lt, f$m$Bt, f$m$Tt, f$t$a, f$p$R, f$m$Rsett, f$m$i0, f$p$phi0, nrow(f$m$Pt)-length(f$m$Rsett), 0, 1, nrow(f$m$Pt), FALSE)
