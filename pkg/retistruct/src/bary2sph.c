#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>


/* Given a M by 3 matrix dF of force components that act on the first
   vertex of M elements, add these to the N by 3 force vector F using the
   vertex indices supplied in SEXP */
SEXP bary2sph(SEXP sidx, SEXP sb, SEXP sT, SEXP sP) {
  /* Read in points */
  int    *xidx = INTEGER(sidx);
  double *xb   = REAL(sb);
  int    *xT   = INTEGER(sT);
  double *xP   = REAL(sP);

  
	if (!isMatrix(sb) || !isReal(sb)) {
		error("sb should be a real matrix.");
	}
  if (ncols(sb) != 3) {
    error("b should have 3 columns.");
  }
  if (LENGTH(sidx) != nrows(sb)) {
    error("idx should have the same length as the number of rows in sb");
  }
	if (!isMatrix(sT) || !isInteger(sT)) {
		error("T should be a integer matrix.");
	}
  if (ncols(sT) != 3) {
    error("T should have 3 columns.");
  }
  
  /* Number of points */
  int    Nb   = nrows(sb);
  int    NT   = nrows(sT);
  int    NP   = nrows(sP);
  
  /* Make space for output */
  SEXP sPs;
  PROTECT(sPs = allocMatrix(REALSXP, Nb, 2));
  double *xPs = REAL(sPs);

  int i, j;
  int T;
  double X[3], Y[3], Z[3];
  double beta[3];
  double Xc, Yc, Zc;

  /* For each barycentric coordinate*/
  for (i=0; i<Nb; i++) {
    if (ISNA(xidx[i]) | ISNA(xb[i]) | ISNA(xb[i+Nb]) | ISNA(xb[i+2*Nb])) {
      xPs[i] = NA_REAL;
      xPs[i + Nb] = NA_REAL; 
    } else {
      /* Find the coordinates of these points */
      for (j=0; j<3; j++) {
        /* Find vertices of triangle idx[i] */
        T = xT[xidx[i] - 1 + j*NT] - 1;
        X[j] = xP[T];
        Y[j] = xP[T + NP];
        Z[j] = xP[T + 2*NP];
      }
      /* Get barycentric coordinates within triangle */
      for (j=0; j<3; j++) {
        beta[j] = xb[i + Nb*j];
      }
      /* Convert to Cartesian 3D coordinate */
      /* Pc = xbeta %*% xX */
      Xc = Yc = Zc = 0;
      for (j=0; j<3; j++) {
        Xc += beta[j]*X[j];
        Yc += beta[j]*Y[j];
        Zc += beta[j]*Z[j];
      }
      xPs[i] = atan2(Zc, sqrt(Xc*Xc + Yc*Yc));
      xPs[i + Nb] = atan2(Yc, Xc);
    }
  }
  SEXP colnames, dimnames;
  PROTECT(colnames = allocVector(STRSXP, 2));
  SET_STRING_ELT(colnames, 0, mkChar("phi"));
  SET_STRING_ELT(colnames, 1, mkChar("lambda"));
  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 1, colnames);
  dimnamesgets(sPs, dimnames);
  UNPROTECT(3);
  return(sPs);
}
