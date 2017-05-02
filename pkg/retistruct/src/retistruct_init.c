#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP sum_force_components(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"sum_force_components", (DL_FUNC) &sum_force_components, 3},
    {NULL, NULL, 0}
};

void R_init_retistruct(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
