#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP crs_prod_kernel_matrix(SEXP Z_s, SEXP z_s, SEXP lambda_s, SEXP ordered_s)
{
    SEXP Z = PROTECT(Rf_coerceVector(Z_s, REALSXP));
    SEXP z = PROTECT(Rf_coerceVector(z_s, REALSXP));
    SEXP lambda = PROTECT(Rf_coerceVector(lambda_s, REALSXP));
    SEXP ordered = PROTECT(Rf_coerceVector(ordered_s, LGLSXP));
    SEXP dim = PROTECT(Rf_getAttrib(Z, R_DimSymbol));

    if (dim == R_NilValue || XLENGTH(dim) != 2) {
        UNPROTECT(5);
        Rf_error("Z must be a matrix");
    }

    const int *dims = INTEGER(dim);
    const R_xlen_t n = (R_xlen_t) dims[0];
    const R_xlen_t p = (R_xlen_t) dims[1];

    if (XLENGTH(z) < p || XLENGTH(lambda) < p || XLENGTH(ordered) < p) {
        UNPROTECT(5);
        Rf_error("incompatible dimensions for Z, z, lambda, and is.ordered.z");
    }

    SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
    double *ans = REAL(out);
    const double *Zp = REAL(Z);
    const double *zp = REAL(z);
    const double *lp = REAL(lambda);
    const int *op = LOGICAL(ordered);

    for (R_xlen_t i = 0; i < n; ++i) ans[i] = 1.0;

    for (R_xlen_t j = 0; j < p; ++j) {
        const double zj = zp[j];
        const double lj = lp[j];
        const int is_ordered = (op[j] == TRUE);
        const R_xlen_t offset = j * n;

        if (!is_ordered) {
            for (R_xlen_t i = 0; i < n; ++i) {
                const double Zij = Zp[offset + i];
                double kij = lj;
                if (!ISNAN(Zij) && !ISNAN(zj) && Zij == zj) kij = 1.0;
                ans[i] *= kij;
            }
        } else {
            for (R_xlen_t i = 0; i < n; ++i) {
                const double Zij = Zp[offset + i];
                double kij = R_pow(lj, fabs(Zij - zj));
                if (!ISNAN(Zij) && !ISNAN(zj) && Zij == zj) kij = 1.0;
                ans[i] *= kij;
            }
        }
    }

    UNPROTECT(6);
    return out;
}
