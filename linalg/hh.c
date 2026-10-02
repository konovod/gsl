/* linalg/hh.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
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

/* Author:  G. Jungman */

/* Originally [Engeln-Mullges + Uhlig, Alg. 4.42], which only handled square
 * systems.  The dimension tests were transposed, so overdetermined systems
 * (M > N) were rejected as "underdetermined" and the in-place solve mixed up
 * the row and column counts; see Savannah bug #42472.  The routine below is a
 * standard tall Householder QR reduction and solves the square and
 * least-squares cases uniformly. */

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#define REAL double

/* Solve min ||A x - b|| for an M x N matrix A with M >= N.  b has length M
 * and holds the right-hand side; x has length N.  A is overwritten with the
 * R factor and the reflectors, and b with Q^T b on exit (so for the in-place
 * case b and x may be the same vector only when M == N).
 */

static int
HH_solve_QR (gsl_matrix * A, gsl_vector * b, gsl_vector * x)
{
  const size_t M = A->size1;    /* number of rows    */
  const size_t N = A->size2;    /* number of columns */
  const size_t rank = GSL_MIN (M, N);
  const size_t vsize = (M > 0) ? M : 1;
  size_t i, j, k;
  REAL *v = (REAL *) malloc (vsize * sizeof (REAL));
  REAL *Rdiag = (REAL *) malloc (((N > 0) ? N : 1) * sizeof (REAL));

  if (v == 0 || Rdiag == 0)
    {
      free (v);
      free (Rdiag);
      GSL_ERROR ("could not allocate memory for workspace", GSL_ENOMEM);
    }

  /* Reduce the columns of A to upper triangular form with Householder
   * reflectors, applying each reflector to the right-hand side as it is
   * formed.  For M > N only the first N columns are reduced; the trailing
   * M - N rows of the transformed A are the residual and are not used. */

  for (i = 0; i < rank; i++)
    {
      const REAL aii = gsl_matrix_get (A, i, i);
      REAL alpha, beta, vnorm2;

      /* Householder vector for column i, rows i..M-1 */

      vnorm2 = 0.0;
      for (k = i; k < M; k++)
        {
          const REAL aki = gsl_matrix_get (A, k, i);
          vnorm2 += aki * aki;
        }
      vnorm2 = sqrt (vnorm2);

      if (vnorm2 == 0.0)
        {
          free (v);
          free (Rdiag);
          GSL_ERROR ("matrix is rank deficient", GSL_ESING);
        }

      /* Choose the sign that avoids cancellation in v[i]. */

      alpha = (aii >= 0.0) ? -vnorm2 : vnorm2;

      v[i] = aii - alpha;
      for (k = i + 1; k < M; k++)
        v[k] = gsl_matrix_get (A, k, i);

      vnorm2 = 0.0;
      for (k = i; k < M; k++)
        vnorm2 += v[k] * v[k];

      if (vnorm2 == 0.0)
        {
          free (v);
          free (Rdiag);
          GSL_ERROR ("apparent singularity detected", GSL_ESING);
        }

      beta = 2.0 / vnorm2;

      /* Apply H = I - beta v v^T to the trailing columns of A. */

      for (k = i; k < N; k++)
        {
          REAL dot = 0.0;
          REAL f;
          for (j = i; j < M; j++)
            dot += v[j] * gsl_matrix_get (A, j, k);
          f = beta * dot;
          for (j = i; j < M; j++)
            {
              REAL ajk = gsl_matrix_get (A, j, k);
              gsl_matrix_set (A, j, k, ajk - f * v[j]);
            }
        }

      /* Apply H to the right-hand side. */

      {
        REAL dot = 0.0;
        REAL f;
        for (j = i; j < M; j++)
          dot += v[j] * gsl_vector_get (b, j);
        f = beta * dot;
        for (j = i; j < M; j++)
          gsl_vector_set (b, j, gsl_vector_get (b, j) - f * v[j]);
      }

      Rdiag[i] = alpha;

      /* Store the R diagonal back in A. The subdiagonal entries of column i
       * still hold the reflector components; they are overwritten by the
       * next reflectors and are not read again. */

      gsl_matrix_set (A, i, i, alpha);
    }

  /* Solve R x = (Q^T b)[0..N-1].  When b and x alias (in-place, M == N) we
   * must go from the bottom up: the update reads b[i] before x[i] is written,
   * and only reads x[k] for k > i, which are already the solution values. */

  for (i = N; i-- > 0;)
    {
      REAL sum = gsl_vector_get (b, i);
      for (k = i + 1; k < N; k++)
        sum -= gsl_matrix_get (A, i, k) * gsl_vector_get (x, k);

      if (Rdiag[i] == 0.0)
        {
          free (v);
          free (Rdiag);
          GSL_ERROR ("matrix is rank deficient", GSL_ESING);
        }

      gsl_vector_set (x, i, sum / Rdiag[i]);
    }

  free (v);
  free (Rdiag);
  return GSL_SUCCESS;
}

int
gsl_linalg_HH_solve (gsl_matrix * A, const gsl_vector * b, gsl_vector * x)
{
  const size_t M = A->size1;    /* number of rows    */
  const size_t N = A->size2;    /* number of columns */

  if (N > M)
    {
      /* Fewer equations than unknowns: there is no unique solution. */

      GSL_ERROR ("System is underdetermined", GSL_EINVAL);
    }
  else if (b->size != M || x->size != N)
    {
      GSL_ERROR ("matrix and vector sizes must be equal", GSL_EBADLEN);
    }
  else if (N == 0)
    {
      return GSL_SUCCESS;
    }
  else
    {
      int status;
      gsl_vector * bwork = gsl_vector_alloc (M);

      if (bwork == 0)
        {
          GSL_ERROR ("could not allocate memory for workspace", GSL_ENOMEM);
        }

      /* b may be longer than x (and must not be modified), so work on a copy. */

      gsl_vector_memcpy (bwork, b);

      status = HH_solve_QR (A, bwork, x);

      gsl_vector_free (bwork);

      return status;
    }
}

int
gsl_linalg_HH_svx (gsl_matrix * A, gsl_vector * x)
{
  const size_t M = A->size1;    /* number of rows    */
  const size_t N = A->size2;    /* number of columns */

  if (N > M)
    {
      /* Fewer equations than unknowns: there is no unique solution. */

      GSL_ERROR ("System is underdetermined", GSL_EINVAL);
    }
  else if (N != x->size)
    {
      GSL_ERROR ("matrix and vector sizes must be equal", GSL_EBADLEN);
    }
  else if (M != N)
    {
      /* In-place x has room only for the N solution values, so it cannot
       * hold the M-element right-hand side of an overdetermined system;
       * use gsl_linalg_HH_solve for M > N. */

      GSL_ERROR ("in-place solve requires a square matrix", GSL_EBADLEN);
    }
  else if (N == 0)
    {
      return GSL_SUCCESS;
    }
  else
    {
      return HH_solve_QR (A, x, x);
    }
}
