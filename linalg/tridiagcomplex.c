/* linalg/tridiagcomplex.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2002, 2004, 2007 Gerard Jungman, Brian Gough, David Necas
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

/* Author: G. Jungman */

#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include "tridiag.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex_math.h>

#define COMPLEX_EQ(z1, z2) (GSL_REAL(z1) == GSL_REAL(z2) && GSL_IMAG(z1) == GSL_IMAG(z2))

/* for description of method see [Engeln-Mullges + Uhlig, p. 92]
 *
 *     diag[0]  offdiag[0]             0   .....
 *  offdiag[0]     diag[1]    offdiag[1]   .....
 *           0  offdiag[1]       diag[2]
 *           0           0    offdiag[2]   .....
 */
static
int
solve_tridiag(
  const double diagdata[], size_t d_stride,
  const double offdiagdata[], size_t o_stride,
  const double bdata[], size_t b_stride,
  double xdata[], size_t x_stride,
  size_t N)
{
  int status = GSL_SUCCESS;
  gsl_complex *gamma = (gsl_complex *) malloc (N * sizeof (gsl_complex));
  gsl_complex *alpha = (gsl_complex *) malloc (N * sizeof (gsl_complex));
  gsl_complex *c = (gsl_complex *) malloc (N * sizeof (gsl_complex));
  gsl_complex *z = (gsl_complex *) malloc (N * sizeof (gsl_complex));
  const gsl_complex *diag = (const gsl_complex *)diagdata;
  const gsl_complex *offdiag = (const gsl_complex *)offdiagdata;
  const gsl_complex *b = (const gsl_complex *)bdata;
  gsl_complex *x = (gsl_complex *)xdata;

  if (gamma == 0 || alpha == 0 || c == 0 || z == 0)
    {
      GSL_ERROR("failed to allocate working space", GSL_ENOMEM);
    }
  else
    {
      size_t i, j;

      /* Cholesky decomposition
         A = L.D.L^t
         lower_diag(L) = gamma
         diag(D) = alpha
       */
      alpha[0] = diag[0];
      gamma[0] = gsl_complex_div(offdiag[0], alpha[0]);

      if (COMPLEX_EQ(alpha[0], GSL_COMPLEX_ZERO)) {
        status = GSL_EZERODIV;
      }

      for (i = 1; i < N - 1; i++)
        {
          alpha[i] = gsl_complex_sub(diag[d_stride * i],
            gsl_complex_mul(offdiag[o_stride*(i - 1)], gamma[i - 1]));
          gamma[i] = gsl_complex_div(offdiag[o_stride * i], alpha[i]);
          if (COMPLEX_EQ(alpha[i], GSL_COMPLEX_ZERO)) {
            status = GSL_EZERODIV;
          }
        }

      if (N > 1)
        {
          alpha[N - 1] = gsl_complex_sub(diag[d_stride * (N - 1)],
            gsl_complex_mul(offdiag[o_stride*(N - 2)], gamma[N - 2]));
        }

      /* update RHS */
      z[0] = b[0];
      for (i = 1; i < N; i++)
        {
          z[i] = gsl_complex_sub(b[b_stride * i], gsl_complex_mul(gamma[i - 1], z[i - 1]));
        }
      for (i = 0; i < N; i++)
        {
          c[i] = gsl_complex_div(z[i], alpha[i]);
        }

      /* backsubstitution */
      x[x_stride * (N - 1)] = c[N - 1];
      if (N >= 2)
        {
          for (i = N - 2, j = 0; j <= N - 2; j++, i--)
            {
              x[x_stride * i] = gsl_complex_sub(c[i], gsl_complex_mul(gamma[i], x[x_stride * (i + 1)]));
            }
        }
    }

  if (z != 0)
    free (z);
  if (c != 0)
    free (c);
  if (alpha != 0)
    free (alpha);
  if (gamma != 0)
    free (gamma);

  if (status == GSL_EZERODIV) {
    GSL_ERROR ("matrix must be positive definite", status);
  }

  return status;
}

/* plain gauss elimination, only not bothering with the zeroes
 *
 *       diag[0]  abovediag[0]             0   .....
 *  belowdiag[0]       diag[1]  abovediag[1]   .....
 *             0  belowdiag[1]       diag[2]
 *             0             0  belowdiag[2]   .....
 */
static
int
solve_tridiag_nonsym(
  const double diagdata[], size_t d_stride,
  const double abovediagdata[], size_t a_stride,
  const double belowdiagdata[], size_t b_stride,
  const double rhsdata[], size_t r_stride,
  double xdata[], size_t x_stride,
  size_t N)
{
  int status = GSL_SUCCESS;
  gsl_complex *alpha = (gsl_complex *) malloc (N * sizeof (gsl_complex));
  gsl_complex *z = (gsl_complex *) malloc (N * sizeof (gsl_complex));
  const gsl_complex *diag = (const gsl_complex *)diagdata;
  const gsl_complex *abovediag = (const gsl_complex *)abovediagdata;
  const gsl_complex *belowdiag = (const gsl_complex *)belowdiagdata;
  const gsl_complex *rhs = (const gsl_complex *)rhsdata;
  gsl_complex *x = (gsl_complex *)xdata;

  if (alpha == 0 || z == 0)
    {
      GSL_ERROR("failed to allocate working space", GSL_ENOMEM);
    }
  else
    {
      size_t i, j;

      /* Bidiagonalization (eliminating belowdiag)
         & rhs update
         diag' = alpha
         rhs' = z
       */
      alpha[0] = diag[0];
      z[0] = rhs[0];

      if (COMPLEX_EQ(alpha[0], GSL_COMPLEX_ZERO)) {
        status = GSL_EZERODIV;
      }

      for (i = 1; i < N; i++)
        {
          const gsl_complex t = gsl_complex_div(belowdiag[b_stride*(i - 1)], alpha[i-1]);
          alpha[i] = gsl_complex_sub(diag[d_stride*i], gsl_complex_mul(t, abovediag[a_stride*(i - 1)]));
          z[i] = gsl_complex_sub(rhs[r_stride*i], gsl_complex_mul(t, z[i-1]));
          if (COMPLEX_EQ(alpha[i], GSL_COMPLEX_ZERO)) {
            status = GSL_EZERODIV;
          }
        }

      /* backsubstitution */
      x[x_stride * (N - 1)] = gsl_complex_div(z[N - 1], alpha[N - 1]);
      if (N >= 2)
        {
          for (i = N - 2, j = 0; j <= N - 2; j++, i--)
            {
              x[x_stride * i] = gsl_complex_div(gsl_complex_sub(z[i], gsl_complex_mul(abovediag[a_stride*i], x[x_stride * (i + 1)])), alpha[i]);
            }
        }
    }

  if (z != 0)
    free (z);
  if (alpha != 0)
    free (alpha);

  if (status == GSL_EZERODIV) {
    GSL_ERROR ("matrix must be positive definite", status);
  }

  return status;
}

/* for description of method see [Engeln-Mullges + Uhlig, p. 96]
 *
 *      diag[0]  offdiag[0]             0   .....  offdiag[N-1]
 *   offdiag[0]     diag[1]    offdiag[1]   .....
 *            0  offdiag[1]       diag[2]
 *            0           0    offdiag[2]   .....
 *          ...         ...
 * offdiag[N-1]         ...
 *
 */
static
int
solve_cyc_tridiag(
  const double diagdata[], size_t d_stride,
  const double offdiagdata[], size_t o_stride,
  const double bdata[], size_t b_stride,
  double xdata[], size_t x_stride,
  size_t N)
{
  int status = GSL_SUCCESS;
  gsl_complex * delta = (gsl_complex *) malloc (N * sizeof (gsl_complex));
  gsl_complex * gamma = (gsl_complex *) malloc (N * sizeof (gsl_complex));
  gsl_complex * alpha = (gsl_complex *) malloc (N * sizeof (gsl_complex));
  gsl_complex * c = (gsl_complex *) malloc (N * sizeof (gsl_complex));
  gsl_complex * z = (gsl_complex *) malloc (N * sizeof (gsl_complex));
  const gsl_complex *diag = (const gsl_complex *)diagdata;
  const gsl_complex *offdiag = (const gsl_complex *)offdiagdata;
  const gsl_complex *b = (const gsl_complex *)bdata;
  gsl_complex *x = (gsl_complex *)xdata;

  if (delta == 0 || gamma == 0 || alpha == 0 || c == 0 || z == 0)
    {
      GSL_ERROR("failed to allocate working space", GSL_ENOMEM);
    }
  else
    {
      size_t i, j;
      gsl_complex sum = GSL_COMPLEX_ZERO;

      /* factor */

      if (N == 1)
        {
          x[0] = gsl_complex_div(b[0], diag[0]);
          free(delta);
          free(gamma);
          free(alpha);
          free(c);
          free(z);
          return GSL_SUCCESS;
        }

      alpha[0] = diag[0];
      gamma[0] = gsl_complex_div(offdiag[0], alpha[0]);
      delta[0] = gsl_complex_div(offdiag[o_stride * (N-1)], alpha[0]);

      if (COMPLEX_EQ(alpha[0], GSL_COMPLEX_ZERO)) {
        status = GSL_EZERODIV;
      }

      for (i = 1; i < N - 2; i++)
        {
          alpha[i] = gsl_complex_sub(diag[d_stride * i], gsl_complex_mul(offdiag[o_stride * (i-1)], gamma[i - 1]));
          gamma[i] = gsl_complex_div(offdiag[o_stride * i], alpha[i]);
          delta[i] = gsl_complex_mul_real(gsl_complex_mul(delta[i - 1], gsl_complex_div(offdiag[o_stride * (i-1)], alpha[i])), -1.0);
          if (COMPLEX_EQ(alpha[i], GSL_COMPLEX_ZERO)) {
            status = GSL_EZERODIV;
          }
        }

      for (i = 0; i < N - 2; i++)
        {
          sum = gsl_complex_add(sum, gsl_complex_mul(gsl_complex_mul(alpha[i], delta[i]), delta[i]));
        }

      alpha[N - 2] = gsl_complex_sub(diag[d_stride * (N - 2)], gsl_complex_mul(offdiag[o_stride * (N - 3)], gamma[N - 3]));

      gamma[N - 2] = gsl_complex_div(gsl_complex_sub(offdiag[o_stride * (N - 2)], gsl_complex_mul(offdiag[o_stride * (N - 3)], delta[N - 3])), alpha[N - 2]);

      alpha[N - 1] = gsl_complex_sub(gsl_complex_sub(diag[d_stride * (N - 1)],
        sum),
        gsl_complex_mul(gsl_complex_mul(alpha[(N - 2)], gamma[N - 2]), gamma[N - 2]));

      /* update */
      z[0] = b[0];
      for (i = 1; i < N - 1; i++)
        {
          z[i] = gsl_complex_sub(b[b_stride * i], gsl_complex_mul(z[i - 1], gamma[i - 1]));
        }
      sum = GSL_COMPLEX_ZERO;
      for (i = 0; i < N - 2; i++)
        {
          sum = gsl_complex_add(sum, gsl_complex_mul(delta[i], z[i]));
        }
      z[N - 1] = gsl_complex_sub(gsl_complex_sub(b[b_stride * (N - 1)],
        sum),
        gsl_complex_mul(gamma[N - 2], z[N - 2]));
      for (i = 0; i < N; i++)
        {
          c[i] = gsl_complex_div(z[i], alpha[i]);
        }

      /* backsubstitution */
      x[x_stride * (N - 1)] = c[N - 1];
      x[x_stride * (N - 2)] = gsl_complex_sub(c[N - 2],
        gsl_complex_mul(gamma[N - 2], x[x_stride * (N - 1)]));
      if (N >= 3)
        {
          for (i = N - 3, j = 0; j <= N - 3; j++, i--)
            {
              x[x_stride * i] = gsl_complex_sub(gsl_complex_sub(c[i],
                gsl_complex_mul(gamma[i], x[x_stride * (i + 1)])),
                gsl_complex_mul(delta[i], x[x_stride * (N - 1)]));
            }
        }
    }

  if (z != 0)
    free (z);
  if (c != 0)
    free (c);
  if (alpha != 0)
    free (alpha);
  if (gamma != 0)
    free (gamma);
  if (delta != 0)
    free (delta);

  if (status == GSL_EZERODIV) {
    GSL_ERROR ("matrix must be positive definite", status);
  }

  return status;
}

/* solve following system w/o the corner elements and then use
 * Sherman-Morrison formula to compensate for them
 *
 *        diag[0]  abovediag[0]             0   .....  belowdiag[N-1]
 *   belowdiag[0]       diag[1]  abovediag[1]   .....
 *              0  belowdiag[1]       diag[2]
 *              0             0  belowdiag[2]   .....
 *            ...           ...
 * abovediag[N-1]           ...
 */
static
int solve_cyc_tridiag_nonsym(
  const double diagdata[], size_t d_stride,
  const double abovediagdata[], size_t a_stride,
  const double belowdiagdata[], size_t b_stride,
  const double rhsdata[], size_t r_stride,
  double xdata[], size_t x_stride,
  size_t N)
{
  int status = GSL_SUCCESS;
  gsl_complex *alpha = (gsl_complex *) malloc (N * sizeof (gsl_complex));
  gsl_complex *zb = (gsl_complex *) malloc (N * sizeof (gsl_complex));
  gsl_complex *zu = (gsl_complex *) malloc (N * sizeof (gsl_complex));
  gsl_complex *w = (gsl_complex *) malloc (N * sizeof (gsl_complex));
  const gsl_complex *diag = (const gsl_complex *)diagdata;
  const gsl_complex *abovediag = (const gsl_complex *)abovediagdata;
  const gsl_complex *belowdiag = (const gsl_complex *)belowdiagdata;
  const gsl_complex *rhs = (const gsl_complex *)rhsdata;
  gsl_complex *x = (gsl_complex *)xdata;

  if (alpha == 0 || zb == 0 || zu == 0 || w == 0)
    {
      GSL_ERROR("failed to allocate working space", GSL_ENOMEM);
    }
  else
    {
      gsl_complex beta;

      /* Bidiagonalization (eliminating belowdiag)
         & rhs update
         diag' = alpha
         rhs' = zb
         rhs' for Aq=u is zu
       */
      zb[0] = rhs[0];
      if (!COMPLEX_EQ(diag[0], GSL_COMPLEX_ZERO)) beta = gsl_complex_mul_real(diag[0], -1.0); else beta = GSL_COMPLEX_ONE;
      {
        const gsl_complex q = gsl_complex_sub(GSL_COMPLEX_ONE,
          gsl_complex_mul(abovediag[0],
            gsl_complex_div(belowdiag[0], gsl_complex_mul(diag[0], diag[d_stride])))
        );
        if (gsl_complex_abs2(gsl_complex_div(q, beta)) > 0.25 && gsl_complex_abs2(gsl_complex_div(q, beta)) < 4) {
          beta = gsl_complex_mul(beta, (gsl_complex_abs2(gsl_complex_div(q, beta)) < 1) ? gsl_complex_rect(0.5,0) : gsl_complex_rect(2,0));
        }
      }
      zu[0] = beta;
      alpha[0] = gsl_complex_sub(diag[0], beta);

      if (COMPLEX_EQ(alpha[0], GSL_COMPLEX_ZERO)) {
        status = GSL_EZERODIV;
      }

      {
        size_t i;
        for (i = 1; i+1 < N; i++)
        {
          const gsl_complex t = gsl_complex_div(belowdiag[b_stride*(i - 1)], alpha[i-1]);
          alpha[i] = gsl_complex_sub(diag[d_stride*i], gsl_complex_mul(t, abovediag[a_stride*(i - 1)]));
          zb[i] = gsl_complex_sub(rhs[r_stride*i], gsl_complex_mul(t, zb[i-1]));
          zu[i] = gsl_complex_mul(gsl_complex_mul_real(t, -1.0), zu[i-1]);
          /* FIXME!!! */
          if (COMPLEX_EQ(alpha[i], GSL_COMPLEX_ZERO)) {
            status = GSL_EZERODIV;
          }
        }
      }

      {
        const size_t i = N-1;
        const gsl_complex t = gsl_complex_div(belowdiag[b_stride*(i - 1)], alpha[i-1]);
        alpha[i] = gsl_complex_sub(gsl_complex_sub(diag[d_stride*i],
                      gsl_complex_div(
                        gsl_complex_mul(abovediag[a_stride*i], belowdiag[b_stride*i]),
                        beta)),
                      gsl_complex_mul(t, abovediag[a_stride*(i - 1)])
                   );
        zb[i] = gsl_complex_sub(rhs[r_stride*i], gsl_complex_mul(t, zb[i-1]));
        zu[i] = gsl_complex_sub(abovediag[a_stride*i], gsl_complex_mul(t, zu[i-1]));
        /* FIXME!!! */
        if (COMPLEX_EQ(alpha[i], GSL_COMPLEX_ZERO)) {
          status = GSL_EZERODIV;
        }
      }

      /* backsubstitution */
      {
        size_t i, j;
        w[N-1] = gsl_complex_div(zu[N-1], alpha[N-1]);
        x[x_stride*(N-1)] = gsl_complex_div(zb[N-1], alpha[N-1]);
        for (i = N - 2, j = 0; j <= N - 2; j++, i--)
          {
            w[i] = gsl_complex_div(gsl_complex_sub(zu[i], gsl_complex_mul(abovediag[a_stride*i], w[i+1])), alpha[i]);
            x[i*x_stride] = gsl_complex_div(gsl_complex_sub(zb[i], gsl_complex_mul(abovediag[a_stride*i], x[x_stride*(i + 1)])), alpha[i]);
          }
      }

      /* Sherman-Morrison */
      {
        const gsl_complex vw = gsl_complex_add(w[0], gsl_complex_mul(gsl_complex_div(belowdiag[b_stride*(N - 1)], beta), w[N-1]));
        const gsl_complex vx = gsl_complex_add(x[0], gsl_complex_mul(gsl_complex_div(belowdiag[b_stride*(N - 1)], beta), x[x_stride*(N - 1)]));
        /* FIXME!!! */
        if (COMPLEX_EQ(gsl_complex_add_real(vw, 1.0), GSL_COMPLEX_ZERO)) {
          status = GSL_EZERODIV;
        }

        {
          size_t i;
          for (i = 0; i < N; i++)
            x[i*x_stride] = gsl_complex_sub(
              x[i*x_stride],
              gsl_complex_mul(
                gsl_complex_div(vx,
                  gsl_complex_add_real(vw, 1.0)),
                w[i])
            );
        }
      }
    }

  if (zb != 0)
    free (zb);
  if (zu != 0)
    free (zu);
  if (w != 0)
    free (w);
  if (alpha != 0)
    free (alpha);

  if (status == GSL_EZERODIV) {
    GSL_ERROR ("matrix must be positive definite", status);
  }

  return status;
}

int
gsl_linalg_complex_solve_symm_tridiag(
  const gsl_vector_complex * diag,
  const gsl_vector_complex * offdiag,
  const gsl_vector_complex * rhs,
  gsl_vector_complex * solution)
{
  if(diag->size != rhs->size)
    {
      GSL_ERROR ("size of diag must match rhs", GSL_EBADLEN);
    }
  else if (offdiag->size != rhs->size-1)
    {
      GSL_ERROR ("size of offdiag must match rhs-1", GSL_EBADLEN);
    }
  else if (solution->size != rhs->size)
    {
      GSL_ERROR ("size of solution must match rhs", GSL_EBADLEN);
    }
  else
    {
      return solve_tridiag(diag->data, diag->stride,
                           offdiag->data, offdiag->stride,
                           rhs->data, rhs->stride,
                           solution->data, solution->stride,
                           diag->size);

    }
}

int
gsl_linalg_complex_solve_tridiag(
  const gsl_vector_complex * diag,
  const gsl_vector_complex * abovediag,
  const gsl_vector_complex * belowdiag,
  const gsl_vector_complex * rhs,
  gsl_vector_complex * solution)
{
  if(diag->size != rhs->size)
    {
      GSL_ERROR ("size of diag must match rhs", GSL_EBADLEN);
    }
  else if (abovediag->size != rhs->size-1)
    {
      GSL_ERROR ("size of abovediag must match rhs-1", GSL_EBADLEN);
    }
  else if (belowdiag->size != rhs->size-1)
    {
      GSL_ERROR ("size of belowdiag must match rhs-1", GSL_EBADLEN);
    }
  else if (solution->size != rhs->size)
    {
      GSL_ERROR ("size of solution must match rhs", GSL_EBADLEN);
    }
  else
    {
      return solve_tridiag_nonsym(diag->data, diag->stride,
                                  abovediag->data, abovediag->stride,
                                  belowdiag->data, belowdiag->stride,
                                  rhs->data, rhs->stride,
                                  solution->data, solution->stride,
                                  diag->size);
    }
}


int
gsl_linalg_complex_solve_symm_cyc_tridiag(
  const gsl_vector_complex * diag,
  const gsl_vector_complex * offdiag,
  const gsl_vector_complex * rhs,
  gsl_vector_complex * solution)
{
  if(diag->size != rhs->size)
    {
      GSL_ERROR ("size of diag must match rhs", GSL_EBADLEN);
    }
  else if (offdiag->size != rhs->size)
    {
      GSL_ERROR ("size of offdiag must match rhs", GSL_EBADLEN);
    }
  else if (solution->size != rhs->size)
    {
      GSL_ERROR ("size of solution must match rhs", GSL_EBADLEN);
    }
  else if (diag->size < 3)
    {
      GSL_ERROR ("size of cyclic system must be 3 or more", GSL_EBADLEN);
    }
  else
    {
      return solve_cyc_tridiag(diag->data, diag->stride,
                               offdiag->data, offdiag->stride,
                               rhs->data, rhs->stride,
                               solution->data, solution->stride,
                               diag->size);
    }
}

int
gsl_linalg_complex_solve_cyc_tridiag(
  const gsl_vector_complex * diag,
  const gsl_vector_complex * abovediag,
  const gsl_vector_complex * belowdiag,
  const gsl_vector_complex * rhs,
  gsl_vector_complex * solution)
{
  if(diag->size != rhs->size)
    {
      GSL_ERROR ("size of diag must match rhs", GSL_EBADLEN);
    }
  else if (abovediag->size != rhs->size)
    {
      GSL_ERROR ("size of abovediag must match rhs", GSL_EBADLEN);
    }
  else if (belowdiag->size != rhs->size)
    {
      GSL_ERROR ("size of belowdiag must match rhs", GSL_EBADLEN);
    }
  else if (solution->size != rhs->size)
    {
      GSL_ERROR ("size of solution must match rhs", GSL_EBADLEN);
    }
  else if (diag->size < 3)
    {
      GSL_ERROR ("size of cyclic system must be 3 or more", GSL_EBADLEN);
    }
  else
    {
      return solve_cyc_tridiag_nonsym(diag->data, diag->stride,
                                      abovediag->data, abovediag->stride,
                                      belowdiag->data, belowdiag->stride,
                                      rhs->data, rhs->stride,
                                      solution->data, solution->stride,
                                      diag->size);
    }
}
