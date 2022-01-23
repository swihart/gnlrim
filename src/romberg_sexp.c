#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stddef.h>


SEXP romberg_sexp(SEXP fcn, SEXP state, SEXP len,
		  SEXP pts, SEXP max, int *err, SEXP envir)
{

  int i, j, j1, finish;
  double errsum, *tab1, *tab2, *x, *fx, *sum, *tmpsum, *zz, *pnt1, *pnt2;

  int *PTS = INTEGER(pts);
  int *MAX = INTEGER(max);
  int *LEN = INTEGER(len);
  
  /*allocate vectors according to member of points and length of vector*/
    tab1 = (double*)R_alloc((size_t)(*PTS)         ,sizeof(double));
    tab2 = (double*)R_alloc((size_t)(*PTS)         ,sizeof(double));
       x = (double*)R_alloc((size_t)(*MAX*(*LEN+1)),sizeof(double));
      fx = (double*)R_alloc((size_t)(*MAX*(*LEN+1)),sizeof(double));
     sum = (double*)R_alloc((size_t)(*LEN)         ,sizeof(double));
  tmpsum = (double*)R_alloc((size_t)(*LEN)         ,sizeof(double));
      zz = (double*)R_alloc((size_t)(*LEN)         ,sizeof(double));
    pnt1 = (double*)R_alloc((size_t)(*LEN)         ,sizeof(double));
    pnt2 = (double*)R_alloc((size_t)(*LEN)         ,sizeof(double));

    if(!tab1||!tab2||!x||!fx||!sum||!tmpsum||!zz||!pnt1||!pnt2){
      *err=1;
      Rf_error("*err is now 1");}
    *err=0;
    for(i=0;i<*LEN;i++)x[i**MAX]=1.0;

  
  SEXP call, result;
  PROTECT(call   = Rf_lang2(fcn,state));
  PROTECT(result = Rf_eval(call,envir));
  SEXP foo;
  PROTECT(foo = Rf_coerceVector(result, REALSXP));
  int len2 = LENGTH(foo);
  for (int i = 0; i < len2; i++)
    if (! R_finite(REAL(foo)[i]))
      Rf_error("function returned vector with non-finite element");
    UNPROTECT(3);
    return foo;
  }
