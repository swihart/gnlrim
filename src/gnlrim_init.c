#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

/* .C calls */
/* None yet */


/* .Call()  calls*/
SEXP romberg_sexp(SEXP fcn, SEXP a, SEXP b, SEXP len, SEXP eps,
                   SEXP pts, SEXP max, SEXP err, SEXP envir);

static const R_CallMethodDef callMethods[]  = {
  {"romberg_sexp", (DL_FUNC) &romberg_sexp, 9},
  {NULL, NULL, 0}
};


/* .Fortran calls */
/* None yet */

void R_init_gnlrim(DllInfo *dll)
{
//    R_registerRoutines(dll, CEntries, callMethods, FortranEntries, NULL);
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
