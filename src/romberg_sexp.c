#include <R.h>
#include <Rinternals.h>



SEXP romberg_sexp(SEXP fcn, SEXP state, SEXP envir)
{


  SEXP call, result;
  PROTECT(call   = lang2(fcn,state));
  PROTECT(result = eval(call,envir));
  SEXP foo;
  PROTECT(foo = coerceVector(result, REALSXP));
  int len=LENGTH(foo);
  for (int i = 0; i < len; i++)
    if (! R_finite(REAL(foo)[i]))
      error("function returned vector with non-finite element");
    UNPROTECT(3);
    return foo;
  }
