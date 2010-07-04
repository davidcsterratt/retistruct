#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

SEXP spheristruct_E(SEXP p, SEXP Cu, SEXP C, SEXP L, SEXP B, SEXP T,
                    SEXP A, SEXP R,
                    SEXP Rset, SEXP i0, SEXP phi0, SEXP Nphi, 
                    SEXP E0A, SEXP kA, SEXP N, SEXP verbose)
{
  SEXP phi, lambda, ans;
  double *xp, *xphi, *xlambda;
  int *xRset;
  /* phi <- rep(phi0, N) */

  PROTECT(phi = allocVector(REALSXP, *INTEGER(N)));
  xp = REAL(p);
  xphi = REAL(phi);
  xRset = INTEGER(Rset);
  int j=0;
  int k=0;
  for (int i=0; i<*INTEGER(N); i++) {
    if (i==xRset[j]) {
      xphi[i] = *REAL(phi0);
      j++;
    } else {
      xphi[i] = xp[k];
      k++;
    }
  }

  PROTECT(lambda = allocVector(REALSXP, *INTEGER(N)));
  xlambda = REAL(lambda);
  j=0;
  k=0;
  for (int i=0; i<*INTEGER(N); i++) {
    if (i==*INTEGER(i0)) {
      xlambda[i] = 0;
      j++;
    } else {
      xlambda[i] = xp[k];
      k++;
    }
  }
  

  /*
  phi[-Rset] <- p[1:Nphi]
  lambda <- rep(0, N)
  lambda[-i0] <- p[Nphi+1:(N-1)]

  ##
  ## Compute derivative of elastic energy
  ##
  ## Use the upper triagular part of the connectivity matrix Cu
  phi1    <- phi[Cu[,1]]
  lambda1 <- lambda[Cu[,1]]
  phi2    <- phi[Cu[,2]]
  lambda2 <- lambda[Cu[,2]]
  l <- R*central.angle(phi1, lambda1, phi2, lambda2)
  if (verbose==2) {
    print(l)
  }
  E.E <- 0.5 * sum((l - L)^2/L)
  ## E.E <- 0.5 * sum((l - L)^2)
  ## E.E <- 0.5*sum((l)^2/L)
  if (verbose>=1) {
    print(E.E)
  }

  E.A <- 0
  if (E0.A) {
    ##
    ## Compute areas
    ##
    P <- R * cbind(cos(phi)*cos(lambda),
                   cos(phi)*sin(lambda),
                   sin(phi))

    ## Find areas of all triangles
    areas <- -0.5/R * dot(P[T[,1],], extprod3d(P[T[,2],], P[T[,3],]))
    ##E.A <- 0.5 * sum((areas - A)^2/A)
    E.A <- sum(exp(-k.A*areas/A))
  }
  
  return(E.E + E0.A*E.A) */

  PROTECT(ans = allocVector(REALSXP, 1));
  REAL(ans)[0] = 45;
  UNPROTECT(2);
  
  return(lambda);
}
