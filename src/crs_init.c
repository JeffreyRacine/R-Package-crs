#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "mgcv.h"

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void RuniqueCombs(void *, void *, void *, void *);
extern void gsl_bspline(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gsl_bspline_deriv(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP smultinomadRSolve(SEXP);
extern SEXP snomadRInfo(SEXP);
extern SEXP snomadRSolve(SEXP);

static const R_CMethodDef CEntries[] = {
    {"RuniqueCombs",      (DL_FUNC) &RuniqueCombs,       4},
    {"gsl_bspline",       (DL_FUNC) &gsl_bspline,        9},
    {"gsl_bspline_deriv", (DL_FUNC) &gsl_bspline_deriv, 11},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"glp_model_tmm",     (DL_FUNC) &glp_model_tmm,     7}, 
    {"mgcv_tmm",          (DL_FUNC) &mgcv_tmm,          5}, 
    {"smultinomadRSolve", (DL_FUNC) &smultinomadRSolve, 1},
    {"snomadRInfo",       (DL_FUNC) &snomadRInfo,       1},
    {"snomadRSolve",      (DL_FUNC) &snomadRSolve,      1},
    {NULL, NULL, 0}
};

void R_init_crs(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
