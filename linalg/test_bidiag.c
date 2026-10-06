/* linalg/test_bidiag.c
 *
 * Copyright (C) 2024 Christian Krueger
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
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

int
test_bidiag_unpack_eps (gsl_matrix * A, double eps, const char *desc)
{
  int s = 0;
  unsigned long i, j, M = A->size1, N = A->size2;

  gsl_matrix *u = gsl_matrix_alloc (M, N);
  gsl_matrix *v = gsl_matrix_alloc (N, N);
  gsl_matrix *v2 = gsl_matrix_alloc (N, N);

  gsl_vector *tau1 = gsl_vector_alloc (N);
  gsl_vector *tau2 = gsl_vector_alloc (N - 1);
  gsl_vector *d = gsl_vector_alloc (N);
  gsl_vector *sd = gsl_vector_alloc (N - 1);


  s += gsl_linalg_bidiag_decomp (A, tau1, tau2);
  s += gsl_linalg_bidiag_unpack (A, tau1, u, tau2, v, d, sd);
  s += gsl_linalg_bidiag_unpack2 (A, tau1, tau2, v2);

  /* compare u and A */
  for (i = 0; i < M; i++)
    {
      for (j = 0; j < N; j++)
        {
          double uij = gsl_matrix_get (u, i, j);
          double Aij = gsl_matrix_get (A, i, j);
          gsl_test_rel (uij, Aij, eps,
                        "%s uA (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                        desc, M, N, i, j, uij, Aij);
        }
    }

  /* compare v and v2 */
  for (i = 0; i < N; i++)
    {
      for (j = 0; j < N; j++)
        {
          double vij = gsl_matrix_get (v, i, j);
          double v2ij = gsl_matrix_get (v2, i, j);
          gsl_test_rel (vij, v2ij, eps,
                        "%s vv2 (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g\n",
                        desc, N, N, i, j, vij, v2ij);
        }
    }

  /* compare d and tau1 */
  for (i = 0; i < N; i++)
    {
      double di = gsl_vector_get (d, i);
      double t1i = gsl_vector_get (tau1, i);
      gsl_test_rel (di, t1i, eps, "%s dtau1 (%3lu)[%lu]: %22.18g   %22.18g\n",
                    desc, N, i, di, t1i);
    }

  /* compare sd and tau2 */
  for (i = 0; i < N - 1; i++)
    {
      double sdi = gsl_vector_get (sd, i);
      double t2i = gsl_vector_get (tau2, i);
      gsl_test_rel (sdi, t2i, eps,
                    "%s sdtau2 (%3lu)[%lu]: %22.18g   %22.18g\n", desc, N - 1,
                    i, sdi, t2i);
    }

  gsl_matrix_free (v2);
  gsl_matrix_free (v);
  gsl_matrix_free (u);
  gsl_vector_free (tau1);
  gsl_vector_free (tau2);
  gsl_vector_free (d);
  gsl_vector_free (sd);

  return s;
}

int
test_bidiag_unpack (gsl_rng * r)
{
  int f;
  int s = 0;
  size_t M, N;
  double acc;

  for (M = 1; M <= 50; M++)
    {
      for (N = 1; N <= M; N++)
        {
          gsl_matrix *A = gsl_matrix_alloc (M, N);

          create_random_matrix (A, r);
          s += test_bidiag_unpack_eps (A, GSL_DBL_EPSILON, "  bidiag_unpack");

          gsl_matrix_free (A);
        }
    }

  return s;
}
