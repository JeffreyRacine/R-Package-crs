#include <R.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>
#include <string.h>

SEXP crs_hat_diag(SEXP qr, SEXP qraux, SEXP rank) {
  if (!isReal(qr)) error("qr must be a real matrix");
  if (!isReal(qraux)) error("qraux must be a real vector");
  if (!isInteger(rank)) error("rank must be an integer");

  SEXP dim = getAttrib(qr, R_DimSymbol);
  if (isNull(dim) || LENGTH(dim) != 2) error("qr must be a matrix");

  int n = INTEGER(dim)[0];
  int k = INTEGER(rank)[0];
  if (k < 0) k = 0;
  if (k > n) k = n;

  SEXP res = PROTECT(allocVector(REALSXP, n));
  double *out = REAL(res);

  if (k == 0) {
    for (int i = 0; i < n; ++i) out[i] = 0.0;
    UNPROTECT(1);
    return res;
  }

  double *y = (double *) R_alloc((size_t)n * (size_t)k, sizeof(double));
  double *qy = (double *) R_alloc((size_t)n * (size_t)k, sizeof(double));
  memset(y, 0, (size_t)n * (size_t)k * sizeof(double));

  for (int i = 0; i < k; ++i) {
    y[i + (size_t)i * n] = 1.0;
  }

  int ny = k;
  F77_CALL(dqrqy)(REAL(qr), &n, &k, REAL(qraux), y, &ny, qy);

  for (int i = 0; i < n; ++i) {
    double sum = 0.0;
    for (int j = 0; j < k; ++j) {
      double v = qy[i + (size_t)j * n];
      sum += v * v;
    }
    out[i] = sum;
  }

  UNPROTECT(1);
  return res;
}
