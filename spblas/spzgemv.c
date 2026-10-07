/* spzgemv.c
 * 
 * Copyright (C) 2026 The GSL Team
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <config.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_blas.h>

/*
gsl_spblas_zgemv()
  Multiply a complex sparse matrix and a complex vector

Inputs: alpha - complex scalar factor
        A     - complex sparse matrix
        x     - dense complex vector
        beta  - complex scalar factor
        y     - (input/output) dense complex vector

Return: y = alpha*op(A)*x + beta*y
*/

int
gsl_spblas_zgemv(const CBLAS_TRANSPOSE_t TransA, const gsl_complex alpha,
                 const gsl_spmatrix_complex *A, const gsl_vector_complex *x,
                 const gsl_complex beta, gsl_vector_complex *y)
{
  const size_t M = A->size1;
  const size_t N = A->size2;
  const int conj = (TransA == CblasConjTrans);

  if ((TransA == CblasNoTrans && N != x->size) ||
      (TransA != CblasNoTrans && M != x->size))
    {
      GSL_ERROR("invalid length of x vector", GSL_EBADLEN);
    }
  else if ((TransA == CblasNoTrans && M != y->size) ||
           (TransA != CblasNoTrans && N != y->size))
    {
      GSL_ERROR("invalid length of y vector", GSL_EBADLEN);
    }
  else
    {
      size_t j;
      size_t incX, incY;
      size_t lenX, lenY;
      double *X, *Y;
      double *Ad;
      int *Ap, *Ai, *Aj;
      int p;
      const double ar = GSL_REAL(alpha);
      const double ai = GSL_IMAG(alpha);
      const double br = GSL_REAL(beta);
      const double bi = GSL_IMAG(beta);

      if (TransA == CblasNoTrans)
        {
          lenX = N;
          lenY = M;
        }
      else
        {
          lenX = M;
          lenY = N;
        }

      /* form y := beta*y */

      Y = y->data;
      incY = y->stride;

      if (br == 0.0 && bi == 0.0)
        {
          size_t jy = 0;
          for (j = 0; j < lenY; ++j)
            {
              Y[2 * jy] = 0.0;
              Y[2 * jy + 1] = 0.0;
              jy += incY;
            }
        }
      else if (!(br == 1.0 && bi == 0.0))
        {
          size_t jy = 0;
          for (j = 0; j < lenY; ++j)
            {
              double y0 = Y[2 * jy];
              double y1 = Y[2 * jy + 1];

              Y[2 * jy] = br * y0 - bi * y1;
              Y[2 * jy + 1] = bi * y0 + br * y1;
              jy += incY;
            }
        }

      if (ar == 0.0 && ai == 0.0)
        return GSL_SUCCESS;

      /* form y := alpha*op(A)*x + y */
      Ap = A->p;
      Ad = A->data;
      X = x->data;
      incX = x->stride;

      if ((GSL_SPMATRIX_ISCCS(A) && TransA == CblasNoTrans) ||
          (GSL_SPMATRIX_ISCRS(A) && TransA != CblasNoTrans))
        {
          Ai = A->i;

          for (j = 0; j < lenX; ++j)
            {
              const double xr = X[2 * j * incX];
              const double xi = X[2 * j * incX + 1];

              for (p = Ap[j]; p < Ap[j + 1]; ++p)
                {
                  const double mr = Ad[2 * p];
                  const double mi = conj ? -Ad[2 * p + 1] : Ad[2 * p + 1];
                  const double tr = mr * xr - mi * xi;
                  const double ti = mr * xi + mi * xr;
                  const size_t idx = 2 * Ai[p] * incY;

                  Y[idx] += ar * tr - ai * ti;
                  Y[idx + 1] += ar * ti + ai * tr;
                }
            }
        }
      else if ((GSL_SPMATRIX_ISCCS(A) && TransA != CblasNoTrans) ||
               (GSL_SPMATRIX_ISCRS(A) && TransA == CblasNoTrans))
        {
          Ai = A->i;

          for (j = 0; j < lenY; ++j)
            {
              const size_t idx = 2 * j * incY;

              for (p = Ap[j]; p < Ap[j + 1]; ++p)
                {
                  const double xr = X[2 * Ai[p] * incX];
                  const double xi = X[2 * Ai[p] * incX + 1];
                  const double mr = Ad[2 * p];
                  const double mi = conj ? -Ad[2 * p + 1] : Ad[2 * p + 1];
                  const double tr = mr * xr - mi * xi;
                  const double ti = mr * xi + mi * xr;

                  Y[idx] += ar * tr - ai * ti;
                  Y[idx + 1] += ar * ti + ai * tr;
                }
            }
        }
      else if (GSL_SPMATRIX_ISTRIPLET(A))
        {
          if (TransA == CblasNoTrans)
            {
              Ai = A->i;
              Aj = A->p;
            }
          else
            {
              Ai = A->p;
              Aj = A->i;
            }

          for (p = 0; p < (int) A->nz; ++p)
            {
              const double xr = X[2 * Aj[p] * incX];
              const double xi = X[2 * Aj[p] * incX + 1];
              const double mr = Ad[2 * p];
              const double mi = conj ? -Ad[2 * p + 1] : Ad[2 * p + 1];
              const double tr = mr * xr - mi * xi;
              const double ti = mr * xi + mi * xr;
              const size_t idx = 2 * Ai[p] * incY;

              Y[idx] += ar * tr - ai * ti;
              Y[idx + 1] += ar * ti + ai * tr;
            }
        }
      else
        {
          GSL_ERROR("unsupported matrix type", GSL_EINVAL);
        }

      return GSL_SUCCESS;
    }
} /* gsl_spblas_zgemv() */
