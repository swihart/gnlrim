#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stddef.h>

SEXP evalRfn_sexp(SEXP fcn, double a[], double b[], int n, int len,
		double sum[], double tmpsum[], double zz[],
		  double pnt1[], double pnt2[], SEXP envir)
{
  double *x, nn;
  int i,j,k;


  SEXP call, result;

  /*Notes to self:  in CJGEYERs `QUX.c`*/
  /* there are the 3 lines:*/
  /* SEXP call, result;
     Rf_protect(call=lang2(func,state));
     Rf_protect(result=eval(call,envir)); */
  /* I think this for if I were to push
     `state` from R (say, the x's I wanted
     ff(x) to evaluate.  But that's not the
     case for my romberg integration.  I just
     want the R funcition ff(x) here in C-land
     with the ability to through in values for
     x that supplied within C-land*/


  /*Rf_protect(call=Rf_lang1(fcn));*/
  /*Rf_protect(result=Rf_eval(fcn, envir));*/




  if(n==1){
    /*one trapezoid*/
    for(k=0;k<len;k++) zz[k]=0.5*(a[k]+b[k]);

    /* Rf_protect(call=Rf_lang2(fcn,zz)); */
    /* Rf_protect(result=Rf_eval(call, envir)); */

    /*Rf_protect(result=Rf_eval(call, envir));*/
    /*Rf_protect(x = Rf_coerceVector(result,REALSXP));*/
    for(k=0;k<len;k++) sum[k]=(b[k]-a[k])*x[k];

    /*return 1;*/
  }
  /*Rf_unprotect(3);   */
  return 1;
}

SEXP romberg_sexp(SEXP fcn, SEXP a, SEXP b, SEXP state, SEXP len, SEXP eps,
		  SEXP pts, SEXP max, int *err, SEXP envir)
{

  int i, j, j1, finish;
  double errsum, *tab1, *tab2, *x, *fx, *sum, *tmpsum, *zz, *pnt1, *pnt2;

  int *PTS = INTEGER(pts);
  int *MAX = INTEGER(max);
  int *LEN = INTEGER(len);
  double *EPS = REAL(eps);
  double *A = REAL(a);
  double *B = REAL(b);

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

  /* iterate, decreasing step size, until convergence or max number of steps */
  for(j=0;j<*MAX;j++){
    j1=j+1;
    /* evaluate R function and calculate sum of trapezoids */
    evalRfn_sexp(fcn,A,B,j1,*LEN,sum,tmpsum,zz,pnt1,pnt2,envir);
     finish=(j1>=*PTS?1:0);
    /*repeatedly call polynomial interpolation routine*/
     for(i=0;i<*LEN;i++){
       fx[j+i**MAX]=sum[i];
       if(j1>=*PTS){
    /* 	interp(&x[j1-*PTS+i**MAX],&fx[j1-*PTS+i**MAX],*PTS,tab1,tab2,&sumlen[i],&errsum,err); */
    if(*err)Rf_error("*err is now 1");
    /*  check convergence  */
	 /* 	if(fabs(errsum)>*EPS*fabs(sumlen[i]))finish=0;*/}
       /* decrease step size */
       x[j1+i**MAX]=x[j+i**MAX]/9.0;
       fx[j1+i**MAX]=fx[j+i**MAX];}
                     }

  SEXP call, result;
  /* try a state2 def and put in c(1,2,3,9) and see if you can get it*/
  SEXP abcd;
  abcd = PROTECT(Rf_allocVector(REALSXP, 4));
  REAL(abcd)[0] = 1;
  REAL(abcd)[1] = 2;
  REAL(abcd)[2] = 2.5;
  REAL(abcd)[3] = 3.11;
  /*...*/

  PROTECT(call   = Rf_lang2(fcn,abcd));
  PROTECT(result = Rf_eval(call,envir));
  SEXP foo;
  PROTECT(foo = Rf_coerceVector(result, REALSXP));
  int len2 = LENGTH(foo);
  for (int i = 0; i < len2; i++)
    if (! R_finite(REAL(foo)[i]))
      Rf_error("function returned vector with non-finite element");
    UNPROTECT(4);
    return foo;
  }
