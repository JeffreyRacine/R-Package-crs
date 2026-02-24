/* We would like to acknowledge the contributions of the GNU GSL
   authors. In particular, we adapt the GNU GSL B-spline routine
   gsl_bspline.c adding automated support for quantile knots (in
   addition to uniform knots), providing missing functionality for
   derivatives, and for extending the splines beyond their
   endpoints. The source files were downloaded from
   http://www.gnu.org/software/gsl/ version 1.14.*, distributed under
   the terms of the GPL, version 2 or later. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gsl_bspline.h"
#include <R.h>
#include <Rinternals.h>

/* Code to replicate bs() in splines package. Note that feeding in
   x_min and x_max is necessary if you want to replicate predict.bs()
   - use x_min/x_max for your training data and let x be the
   evaluation data. 

   knots_int is an integer (0=uniform knots, 1=quantile knots) and
   quantile_vector a vector of knots.*/

int gsl_bspline(double *x,
                int *n,
                int *degree,
                int *nbreak,
                double *x_min,
                double *x_max,
                double *quantile_vector,
                int *knots_int,
                double *Bx)
{

  int k = *degree + 1; /* k in gsl */
  int ncoeffs;
  int i, j;

  gsl_bspline_workspace *bw = gsl_bspline_alloc(k, *nbreak);
  ncoeffs = (int)gsl_bspline_ncoeffs(bw);  /* *nbreak+k-2 */
  gsl_vector *B = gsl_vector_alloc(ncoeffs);
  gsl_vector *quantile_vec = gsl_vector_alloc(*nbreak);

  /* 7/12/10 added support for quantile knots */

  if(*knots_int == 0) {
    gsl_bspline_knots_uniform(*x_min, *x_max, bw);
  } else {
    for(i = 0; i < *nbreak; i++) gsl_vector_set(quantile_vec, i, quantile_vector[i]);
    gsl_bspline_knots(quantile_vec, bw);
  }

  for (i = 0; i < *n; ++i)
    {

      /* compute B_j(xi) for all j */
      gsl_bspline_eval(x[i], B, bw);
      
      /* fill in row i of Bx */
      for (j = 0; j < ncoeffs; ++j)
        {
          double Bj = gsl_vector_get(B, j);
          Bx[i*ncoeffs+j] = Bj; /* Bx:*n-by-(*nbreak+*degree-1) */
        }
    }
  
  gsl_bspline_free(bw);
  gsl_vector_free(B);
  gsl_vector_free(quantile_vec);

  return(0);

} /* main() */

/* Provide missing functionality derivative bs() function in package
   splines */

int gsl_bspline_deriv(double *x,
                      int *n,
                      int *degree,
                      int *nbreak,
                      int *order,
											int *order_max, 
                      double *x_min,
                      double *x_max,
                      double *quantile_vector,
                      int *knots_int,
                      double *Bx)
{

  int k = *degree + 1; /* k in gsl */
  int ncoeffs;
  size_t i, j;

  gsl_bspline_workspace *bw = gsl_bspline_alloc(k, *nbreak);
  ncoeffs = (int)gsl_bspline_ncoeffs(bw);
  gsl_vector *dBorder = gsl_vector_alloc(ncoeffs);
  gsl_bspline_deriv_workspace *derivWS = gsl_bspline_deriv_alloc(k);
	gsl_matrix *dB = gsl_matrix_alloc(ncoeffs, *order_max+1);

	gsl_vector *quantile_vec = gsl_vector_alloc(*nbreak);

	/* 7/12/10 added support for quantile knots */

	if(*knots_int == 0) {
			gsl_bspline_knots_uniform(*x_min, *x_max, bw);
	} else {
			for(i = 0; i < *nbreak; i++) gsl_vector_set(quantile_vec, i, quantile_vector[i]);
			gsl_bspline_knots(quantile_vec, bw);
	}

	for (i = 0; i < *n; ++i)
	{

			/* compute B_j(xi) for all j */
			gsl_bspline_deriv_eval(x[i], order[i], dB, bw, derivWS);

			/* fill in row i of Bx */
			gsl_matrix_get_col(dBorder, dB, order[i]);

			for (j = 0; j < ncoeffs; ++j)
			{
					double Bj = gsl_vector_get(dBorder, j);
					Bx[i*ncoeffs+j] = Bj;
			}
	}

	gsl_bspline_free(bw);
	gsl_vector_free(dBorder);
	gsl_matrix_free(dB);
	gsl_vector_free(quantile_vec);
	gsl_bspline_deriv_free(derivWS);

	return(0);

} /* main() */

SEXP crs_gsl_bspline_call(SEXP x,
                          SEXP degree,
                          SEXP nbreak,
                          SEXP x_min,
                          SEXP x_max,
                          SEXP knots)
{
  SEXP x_r = PROTECT(coerceVector(x, REALSXP));
  SEXP degree_i = PROTECT(coerceVector(degree, INTSXP));
  SEXP nbreak_i = PROTECT(coerceVector(nbreak, INTSXP));
  SEXP x_min_r = PROTECT(coerceVector(x_min, REALSXP));
  SEXP x_max_r = PROTECT(coerceVector(x_max, REALSXP));
  SEXP out = R_NilValue;
  SEXP knots_r = R_NilValue;

  int n = LENGTH(x_r);
  int deg = INTEGER(degree_i)[0];
  int nb = INTEGER(nbreak_i)[0];
  int knots_int = isNull(knots) ? 0 : 1;
  int ncol = nb + deg - 1;
  int i, j;
  double dummy = 0.0;
  double *qvec = &dummy;
  double *bx_rowmajor;

  if (n < 1) {
    UNPROTECT(5);
    Rf_error("x must have positive length");
  }
  if (deg < 1) {
    UNPROTECT(5);
    Rf_error("degree must be positive");
  }
  if (nb < 2) {
    UNPROTECT(5);
    Rf_error("nbreak must be at least 2");
  }
  if (knots_int) {
    knots_r = PROTECT(coerceVector(knots, REALSXP));
    if (LENGTH(knots_r) < nb) {
      UNPROTECT(6);
      Rf_error("knots length is smaller than nbreak");
    }
    qvec = REAL(knots_r);
  }

  bx_rowmajor = (double *) R_alloc((size_t)n * (size_t)ncol, sizeof(double));
  gsl_bspline(REAL(x_r), &n, &deg, &nb, REAL(x_min_r), REAL(x_max_r), qvec, &knots_int, bx_rowmajor);

  PROTECT(out = allocMatrix(REALSXP, n, ncol));
  for (i = 0; i < n; i++) {
    for (j = 0; j < ncol; j++) {
      REAL(out)[i + n * j] = bx_rowmajor[i * ncol + j];
    }
  }

  if (knots_int) {
    UNPROTECT(7);
  } else {
    UNPROTECT(6);
  }
  return out;
}

SEXP crs_gsl_bspline_deriv_call(SEXP x,
                                SEXP degree,
                                SEXP nbreak,
                                SEXP deriv,
                                SEXP x_min,
                                SEXP x_max,
                                SEXP knots)
{
  SEXP x_r = PROTECT(coerceVector(x, REALSXP));
  SEXP degree_i = PROTECT(coerceVector(degree, INTSXP));
  SEXP nbreak_i = PROTECT(coerceVector(nbreak, INTSXP));
  SEXP deriv_i = PROTECT(coerceVector(deriv, INTSXP));
  SEXP x_min_r = PROTECT(coerceVector(x_min, REALSXP));
  SEXP x_max_r = PROTECT(coerceVector(x_max, REALSXP));
  SEXP out = R_NilValue;
  SEXP knots_r = R_NilValue;

  int n = LENGTH(x_r);
  int deg = INTEGER(degree_i)[0];
  int nb = INTEGER(nbreak_i)[0];
  int knots_int = isNull(knots) ? 0 : 1;
  int ncol = nb + deg - 1;
  int i, j;
  int order_max = 0;
  double dummy = 0.0;
  double *qvec = &dummy;
  double *bx_rowmajor;

  if (n < 1) {
    UNPROTECT(6);
    Rf_error("x must have positive length");
  }
  if (LENGTH(deriv_i) != n) {
    UNPROTECT(6);
    Rf_error("deriv must have same length as x");
  }
  if (deg < 1) {
    UNPROTECT(6);
    Rf_error("degree must be positive");
  }
  if (nb < 2) {
    UNPROTECT(6);
    Rf_error("nbreak must be at least 2");
  }

  for (i = 0; i < n; i++) {
    if (INTEGER(deriv_i)[i] < 0) {
      UNPROTECT(6);
      Rf_error("deriv must be non-negative");
    }
    if (INTEGER(deriv_i)[i] > order_max) {
      order_max = INTEGER(deriv_i)[i];
    }
  }

  if (knots_int) {
    knots_r = PROTECT(coerceVector(knots, REALSXP));
    if (LENGTH(knots_r) < nb) {
      UNPROTECT(7);
      Rf_error("knots length is smaller than nbreak");
    }
    qvec = REAL(knots_r);
  }

  bx_rowmajor = (double *) R_alloc((size_t)n * (size_t)ncol, sizeof(double));
  gsl_bspline_deriv(REAL(x_r), &n, &deg, &nb, INTEGER(deriv_i), &order_max,
                    REAL(x_min_r), REAL(x_max_r), qvec, &knots_int, bx_rowmajor);

  PROTECT(out = allocMatrix(REALSXP, n, ncol));
  for (i = 0; i < n; i++) {
    for (j = 0; j < ncol; j++) {
      REAL(out)[i + n * j] = bx_rowmajor[i * ncol + j];
    }
  }

  if (knots_int) {
    UNPROTECT(8);
  } else {
    UNPROTECT(7);
  }
  return out;
}
