#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

double central_angle(double phi1, double lambda1, double phi2, double lambda2) {
  return(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2)));
}

/* Given a M by 3 matrix dF of force components that act on the first
   vertex of M elements, add these to the N by 3 force vector F using the
   vertex indices supplied in C */
SEXP sum_force_components(SEXP sdF, SEXP sC, SEXP sF) {
  double *xdF = REAL(sdF);
  int    *xC  = INTEGER(sC);
  int    N    = nrows(sF);
  int    M    = nrows(sC);

  /* PROTECT(sF = allocMatrix(REALSXP, N, 3)); */
  double *xF = REAL(sF);
  int j;

  for (int i=0; i<M; i++) {
    j = xC[i] - 1;
    xF[j]       += xdF[i];
    xF[j + N]   += xdF[i + M];
    xF[j + 2*N] += xdF[i + 2*M];
  }
  /* UNPROTECT(1); */
  return(sF);
}
