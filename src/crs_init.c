#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
  The following symbols/expressions for .NAME have been omitted

    snomadRInfo
    snomadRSolve
    smultinomadRSolve

  Most likely possible values need to be added below.
*/

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void gsl_bspline(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gsl_bspline_deriv(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RuniqueCombs(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"gsl_bspline",       (DL_FUNC) &gsl_bspline,        9},
    {"gsl_bspline_deriv", (DL_FUNC) &gsl_bspline_deriv, 11},
    {"RuniqueCombs",      (DL_FUNC) &RuniqueCombs,       4},
    {NULL, NULL, 0}
};

void R_init_crs(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
