#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C() */
extern void romberg(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"romberg",    (DL_FUNC) &romberg,     9},
    {NULL, NULL, 0}
};



/* .Call()  */
extern void romberg_sexp(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef callMethods[]  = {
    {"romberg_sexp", (DL_FUNC) &romberg_sexp, 10},
    {NULL, NULL, 0}
};

void R_init_gnlrim(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}