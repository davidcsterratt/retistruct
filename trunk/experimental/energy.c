#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

void spheristruct_p_to_lambda_phi(double *xp, int *Rset, int N, int Nset) {

}

double central_angle(double phi1, double lambda1, double phi2, double lambda2) {
  return(acos(sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lambda1-lambda2)));
}

SEXP spheristruct_E(SEXP p, SEXP Cu, SEXP C, SEXP L, SEXP B, SEXP T,
                    SEXP A, SEXP R,
                    SEXP Rset, SEXP i0, SEXP phi0, SEXP Nphi, 
                    SEXP E0A, SEXP kA, SEXP N, SEXP verbose)
{
  SEXP phi, lambda, ans;
  double *xp, *xphi, *xlambda, *xL;
  int *xRset;
  int NRset;
  int *xCu;
  /* phi <- rep(phi0, N) */

  PROTECT(phi = allocVector(REALSXP, *INTEGER(N)));
  PROTECT(p = AS_NUMERIC(p));
  xp = REAL(p);
  xphi = REAL(phi);
  xRset = INTEGER(Rset);
  NRset = LENGTH(Rset);
  int j=0;
  int k=0;
  for (int i=0; i<*INTEGER(N); i++) {
    if (i==(xRset[j]-1)) {
      xphi[i] = *REAL(phi0);
      if (j + 1 < NRset) j++;
    } else {
      xphi[i] = xp[k];
      k++;
    }
  }

  PROTECT(lambda = allocVector(REALSXP, *INTEGER(N)));
  xlambda = REAL(lambda);
  j=0;
  for (int i=0; i<*INTEGER(N); i++) {
    if (i==(*INTEGER(i0)-1)) {
      xlambda[i] = 0;
      j++;
    } else {
      xlambda[i] = xp[k];
      k++;
    }
  }
  
  PROTECT(Cu = AS_INTEGER(Cu));
  xCu = INTEGER(Cu);
  PROTECT(L = AS_NUMERIC(L));
  xL = REAL(L);
  int M = LENGTH(Cu)/2;

  /* 
   *  Compute derivative of elastic energy
   *  Use the upper triagular part of the connectivity matrix Cu 
   */
  double E=0;
  double l;
  for (int i=0; i<M; i++) {
    /*    printf("%i %i %f %i %i %f %f\n", i, xCu[i], xphi[xCu[i]], 
           i+M, xCu[i+M], xphi[xCu[i+M]],
           central_angle(xphi[xCu[i]],   xlambda[xCu[i]], 
           xphi[xCu[i+M]], xlambda[xCu[i+M]])); */
    l = *REAL(R)*central_angle(xphi[xCu[i]],   xlambda[xCu[i]], 
                        xphi[xCu[i+M]], xlambda[xCu[i+M]]);
    E += (l - xL[i])*(l - xL[i])/xL[i];
  }
  E *= 0.5;
  /*
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
  REAL(ans)[0] = E;
  UNPROTECT(6);
  
  return(ans);
}

SEXP spheristruct_E2(SEXP sphi, SEXP slambda, SEXP sCu, SEXP sL, SEXP sR) {
  int i, j, k;
  double dx, dy, dz, l, dll;
  double dlidpj, dlidpk, dlidlj, dlidlk;
  static double    x[10000],    y[10000],    z[10000];
  static double dxdp[10000], dydp[10000], dzdp[10000];
  static double dxdl[10000], dydl[10000];
  static int frun = 0;
  SEXP dE;
  double *xdE;

  double *xp =  REAL(sphi);
  double *xl =  REAL(slambda);
  double  R  = *REAL(sR);
  double *L  =  REAL(sL);
  int    *C  =  INTEGER(sCu);
  int     N  =  LENGTH(slambda);
  int     M  =  LENGTH(sCu)/2;

  PROTECT(dE   = allocVector(REALSXP, 2*N));
  xdE   = REAL(dE);

  for (i=0; i<N; i++) { 
    x[i]    =  R*cos(xp[i])*cos(xl[i]);
    y[i]    =  R*cos(xp[i])*sin(xl[i]);
    z[i]    =  R*sin(xp[i]);
    if((i == 9) && !frun) {
      printf("x=%f y=%f z=%f\n", x[i], y[i], z[i]);
    }

    dxdp[i] = -R*sin(xp[i])*cos(xl[i]);
    dydp[i] = -R*sin(xp[i])*sin(xl[i]);
    dzdp[i] =  R*cos(xp[i]);
    
    dxdl[i] = -R*cos(xp[i])*sin(xl[i]);
    dydl[i] =  R*cos(xp[i])*cos(xl[i]);
    /* dzdl =  rep(0, N) */
    
    xdE[i]   = 0;
    xdE[i+N] = 0;
  }
  
  for (i=0; i<M; i++) {
    j = C[i]-1;
    k = C[i+M]-1;

    dx = x[j]-x[k];
    dy = y[j]-y[k];
    dz = z[j]-z[k];
    
    l = sqrt(dx*dx + dy*dy + dz*dz);
    
    dlidpj =  1.0/l*(dx*dxdp[j] + dy*dydp[j] + dz*dzdp[j]);
    dlidpk = -1.0/l*(dx*dxdp[k] + dy*dydp[k] + dz*dzdp[k]);
    dlidlj =  1.0/l*(dx*dxdl[j] + dy*dydl[j]);
    dlidlk = -1.0/l*(dx*dxdl[k] + dy*dydl[k]);

    dll = (l - L[i])/L[i];
      
    xdE[j]   = xdE[j] +   dlidpj*dll;
    xdE[k]   = xdE[k] +   dlidpk*dll;
    xdE[j+N] = xdE[j+N] + dlidlj*dll;
    xdE[k+N] = xdE[k+N] + dlidlk*dll;
  }
  UNPROTECT(1);
  frun = 1;
  return(dE);
}
