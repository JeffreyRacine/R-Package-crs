/* Copyright (C) 2000-2012 Simon N. Wood  simon.wood@r-project.org

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
(www.gnu.org/copyleft/gpl.html)

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
USA. */

#include "mgcv.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>

void mgcv_tensor_mm(double *X,double *T,int *d,int *m,int *n) {
/* Code for efficient production of row tensor product model matrices.
   X contains rows of matrices to be producted. Contains m matrices,
   d[i] is number of columns in ith matrix (which are stored in ascending 
   order). Each column has n rows. T is the target matrix, with n rows
   and \prod_i d[i] columns.   
*/
  ptrdiff_t start, i,j,k, tp=1, xp=0,pd;
  double *Xj,*Xj1,*Xk, *Tk,*p,*p1,*p2;
  /*Rprintf("m = %d  n = %d  d = ",*m,*n);
    for (i=0;i<*m;i++) Rprintf(" %d,",d[i]);*/
  /* compute total columns in X, xp, and total columns in T, tp ... */
  for (i=0;i<*m;i++) { xp += d[i];tp *= d[i];} 
  Xk = X + (xp-d[*m-1]) * *n; /* start of last matrix in X */
  Tk = T + (tp-d[*m-1]) * *n; /* start of last (filled) block in T */
  /* initialize by putting final matrix in X into end block of T... */
  p = Xk;p1 = Xk + *n * (ptrdiff_t) d[*m-1];p2 = Tk;
  for (;p<p1;p++,p2++) *p2 = *p;
  pd = d[*m-1]; /* cols in filled (end) block of T */
  for (i = *m - 2;i>=0;i--) { /* work down through matrices stored in X */
    Xk -= *n * (ptrdiff_t) d[i]; /* start of ith X matrix */
    Xj = Xk;
    start = tp - pd * d[i]; /* start column of target block in T */
    p = T + start * *n; /* start location in T */
    for (j=0;j<d[i];j++) { /* work through columns of ith X matrix */
      p1 = Tk; /* start of T block to be multiplied by jth col of ith X matrix*/
      Xj1 = Xj + *n; /* end of col j of ith X matrix */
      /* following multiplies end T block by jth col of ith X matrix, storing 
         result in appropriate T block... */
      for (k=0;k<pd;k++) for (p2 = Xj;p2<Xj1;p2++,p1++,p++) *p = *p1 * *p2;
      Xj += *n; /* start of col j+1 of ith X matrix */
    }
    pd *= d[i]; /* number of cols in filled in T block */
    Tk = T + (tp - pd) * *n; /* start of filled in T block */
  }  
} /* mgcv_tensor_mm */

SEXP mgcv_tmm(SEXP x,SEXP t,SEXP D,SEXP M, SEXP N) {
  /* wrapper for calling mgcv_tensor_mm using .Call */
  double *X,*T;
  int *d,*m,*n;
  X = REAL(x);T=REAL(t);
  d = INTEGER(D);
  m = INTEGER(M);
  n = INTEGER(N);
  mgcv_tensor_mm(X,T,d,m,n);
  return(R_NilValue);
}

/* 
   x contains rows of matrices to be producted. Contains m matrices,
   d[i] is number of columns in ith matrix. Each column has n rows. 
	 t is the target matrix, with n rows and zn columns.
	 z contains a matrix of integers with m columns and zn rows.
*/
void glp_model_mm(double *x, int *z, double *t, int *d, int *m, int *n, int *zn)
{
		ptrdiff_t i,j,k, tsize;
		double *xi;  // the ith column of a given matrix.

		tsize = *zn * *n;  // the size of the matrix t.
		// the first column of z. fill with the first x matrix.
		for(i = 0; i < *zn; ++i)
		{
				if (*z == 0)
				{
						for(j = 0; j < *n; ++j)
								*t++ = 1.0;
				}
				else
				{
						xi = x + (*z - 1) * *n;
						for(j = 0; j < *n; ++j)
								*t++ = *xi++;
				}
				++z;
		}
		t -= tsize; // go back to the start pointer.

		if(*m > 1)
		{
				// the kth matrix
				for(k = 1; k < *m; ++k)
				{
						x += d[k-1] * *n; // go the next matrix.

						for(i = 0; i < *zn; ++i)
						{
								if (*z == 0)
								{
												t += *n;
								}
								else
								{
										xi = x + (*z - 1) * *n;
										for(j = 0; j < *n; ++j, ++t, ++xi)
												*t *= *xi;
								}
								++z;
						}
						t -= tsize; // go back to the start pointer.
				}
		}
		return;
}

SEXP glp_model_tmm(SEXP X,SEXP Z, SEXP T, SEXP D,SEXP M, SEXP N, SEXP ZN)
{
		/* wrapper for calling glp_tmodel_mm using .Call */
		double *x,*t;
		int *d,*m,*n, *zn, *z;

		x = REAL(X);
		z = INTEGER(Z);
		t=REAL(T);
		d = INTEGER(D);
		m = INTEGER(M);
		n = INTEGER(N);
		zn = INTEGER(ZN);

		glp_model_mm(x, z, t, d, m, n, zn);
                return(R_NilValue);

}
