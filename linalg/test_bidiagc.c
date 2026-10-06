/* linalg/test_bidiagc.c
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
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>

int
test_bidiag_complex_decomp_dim (const gsl_matrix_complex * m, double eps)
{
  int s = 0;
  unsigned long i, j, k, r, M = m->size1, N = m->size2;

  gsl_matrix_complex *A = gsl_matrix_complex_alloc (M, N);
  gsl_matrix_complex *a = gsl_matrix_complex_alloc (M, N);
  gsl_matrix *b = gsl_matrix_alloc (N, N);

  gsl_matrix_complex *u = gsl_matrix_complex_alloc (M, N);
  gsl_matrix_complex *v = gsl_matrix_complex_alloc (N, N);

  gsl_vector_complex *tau1 = gsl_vector_complex_alloc (N);
  gsl_vector_complex *tau2 = gsl_vector_complex_alloc (N - 1);
  gsl_vector *d = gsl_vector_alloc (N);
  gsl_vector *sd = gsl_vector_alloc (N - 1);

  gsl_vector_complex * w = gsl_vector_complex_alloc(N);

  gsl_matrix_complex_memcpy (A, m);


  s += gsl_linalg_complex_bidiag_decomp (A, tau1, tau2);
  s += gsl_linalg_complex_bidiag_unpack (A, tau1, u, tau2, v, d, sd, w);

  gsl_matrix_set_zero (b);
  for (i = 0; i < N; i++)
    gsl_matrix_set (b, i, i, gsl_vector_get (d, i));
  for (i = 0; i < N - 1; i++)
    gsl_matrix_set (b, i, i + 1, gsl_vector_get (sd, i));


  /* Compute A = U B V^H */

  for (i = 0; i < M; i++)
    {
      for (j = 0; j < N; j++)
        {
          gsl_complex sum = GSL_COMPLEX_ZERO;
          gsl_complex h1;

          for (k = 0; k < N; k++)
            {
              for (r = 0; r < N; r++)
                {
                  h1 = gsl_matrix_complex_get (u, i, k);
                  h1 = gsl_complex_mul_real (h1, gsl_matrix_get (b, k, r));
                  h1 =
                    gsl_complex_mul (h1,
                                     gsl_complex_conjugate
                                     (gsl_matrix_complex_get (v, j, r)));
                  sum = gsl_complex_add (sum, h1);
                }
            }
          gsl_matrix_complex_set (a, i, j, sum);
        }
    }


  /* Compare original matrix with the decomposed and reconstructed one */

  {
    double scale = 0.0;

    for (i = 0; i < M; i++)
      {
        for (j = 0; j < N; j++)
          {
            double t = gsl_complex_abs (gsl_matrix_complex_get (m, i, j));
            if (t > scale) scale = t;
          }
      }

    for (i = 0; i < M; i++)
      {
        for (j = 0; j < N; j++)
          {
            gsl_complex aij = gsl_matrix_complex_get (a, i, j);
            gsl_complex mij = gsl_matrix_complex_get (m, i, j);
            double d = gsl_complex_abs (gsl_complex_sub (aij, mij));
            int foo = (d > eps * (gsl_complex_abs (mij) + scale));
            if (foo)
              {
                printf ("(%3lu,%3lu)[%lu,%lu]: %22.18g%+22.18gi   %22.18g%+22.18gi\n",
                        M, N, i, j, GSL_REAL (aij), GSL_IMAG (aij),
                        GSL_REAL (mij), GSL_IMAG (mij));
              }
            s += foo;
          }
      }
  }

  gsl_matrix_complex_free (A);
  gsl_matrix_complex_free (a);
  gsl_matrix_complex_free (u);
  gsl_matrix_complex_free (v);
  gsl_vector_complex_free (w);
  gsl_matrix_free (b);
  gsl_vector_complex_free (tau1);
  gsl_vector_complex_free (tau2);
  gsl_vector_free (d);
  gsl_vector_free (sd);

  return s;
}

int
test_bidiag_complex_decomp (gsl_rng * r)
{
  int f;
  int s = 0;
  size_t M, N;
  double acc;

  for (M = 1; M <= 50; M++)
    {
      for (N = 1; N <= M; N++)
        {
          gsl_matrix_complex *A = gsl_matrix_complex_alloc (M, N);

          acc = (M < 10) ? 1e4 * GSL_DBL_EPSILON
              : (M < 48) ? 1e6 * GSL_DBL_EPSILON
                         : 1e8 * GSL_DBL_EPSILON;

          create_random_complex_matrix (A, r);
          f = test_bidiag_complex_decomp_dim (A, acc);
          gsl_test (f, "  bidiag_complex_decomp");
          s += f;
        }
    }

  return s;
}

int
test_bidiag_complex_unpack_eps (gsl_matrix_complex * A, double eps, const char *desc)
{
  int s = 0;
  unsigned long i, j, M = A->size1, N = A->size2;

  gsl_matrix_complex *u = gsl_matrix_complex_alloc (M, N);
  gsl_matrix_complex *v = gsl_matrix_complex_alloc (N, N);
  gsl_matrix_complex *v2 = gsl_matrix_complex_alloc (N, N);

  gsl_vector_complex *tau1 = gsl_vector_complex_alloc (N);
  gsl_vector_complex *tau2 = gsl_vector_complex_alloc (N - 1);
  gsl_vector *d = gsl_vector_alloc (N);
  gsl_vector *sd = gsl_vector_alloc (N - 1);
  gsl_vector *d2 = gsl_vector_alloc (N);
  gsl_vector *sd2 = gsl_vector_alloc (N - 1);

  gsl_vector_complex *work = gsl_vector_complex_alloc (N);


  s += gsl_linalg_complex_bidiag_decomp (A, tau1, tau2);
  s += gsl_linalg_complex_bidiag_unpack (A, tau1, u, tau2, v, d, sd, work);
  s += gsl_linalg_complex_bidiag_unpack2 (A, tau1, tau2, v2, d2, sd2);

  /* compare u and A */
  {
    double scale = 0.0;

    for (i = 0; i < M; i++)
      {
        for (j = 0; j < N; j++)
          {
            double t = gsl_complex_abs (gsl_matrix_complex_get (A, i, j));
            if (t > scale) scale = t;
          }
      }

    for (i = 0; i < M; i++)
      {
        for (j = 0; j < N; j++)
          {
            gsl_complex uij = gsl_matrix_complex_get (u, i, j);
            gsl_complex Aij = gsl_matrix_complex_get (A, i, j);
            double d = gsl_complex_abs (gsl_complex_sub (uij, Aij));
            int foo = (d > eps * (gsl_complex_abs (Aij) + scale));
            if (foo)
              {
                printf ("%s uA (%3lu,%3lu)[%lu,%lu]: %22.18g%+22.18gi   %22.18g%+22.18gi\n",
                        desc, M, N, i, j, GSL_REAL (uij), GSL_IMAG (uij),
                        GSL_REAL (Aij), GSL_IMAG (Aij));
              }
            s += foo;
          }
      }
  }

  /* compare v and v2 */
  {
    double scale = 0.0;

    for (i = 0; i < N; i++)
      {
        for (j = 0; j < N; j++)
          {
            double t = gsl_complex_abs (gsl_matrix_complex_get (v, i, j));
            if (t > scale) scale = t;
          }
      }

    for (i = 0; i < N; i++)
      {
        for (j = 0; j < N; j++)
          {
            gsl_complex vij = gsl_matrix_complex_get (v, i, j);
            gsl_complex v2ij = gsl_matrix_complex_get (v2, i, j);
            double d = gsl_complex_abs (gsl_complex_sub (vij, v2ij));
            int foo = (d > eps * (gsl_complex_abs (v2ij) + scale));
            if (foo)
              {
                printf ("%s vv2 (%3lu,%3lu)[%lu,%lu]: %22.18g%+22.18gi   %22.18g%+22.18gi\n",
                        desc, N, N, i, j, GSL_REAL (vij), GSL_IMAG (vij),
                        GSL_REAL (v2ij), GSL_IMAG (v2ij));
              }
            s += foo;
          }
      }
  }

  /* compare d and d2, and sd and sd2, with an absolute scale */
  {
    double scale = 0.0;

    for (i = 0; i < N; i++)
      {
        double t = fabs (gsl_vector_get (d, i));
        if (t > scale) scale = t;
      }
    for (i = 0; i < N - 1; i++)
      {
        double t = fabs (gsl_vector_get (sd, i));
        if (t > scale) scale = t;
      }

    for (i = 0; i < N; i++)
      {
        double di = gsl_vector_get (d, i);
        double d2i = gsl_vector_get (d2, i);
        gsl_test (fabs (di - d2i) > eps * (fabs (d2i) + scale),
                  "%s dtau1 (%3lu)[%lu]: %22.18g   %22.18g", desc, N, i, di, d2i);
      }

    for (i = 0; i < N - 1; i++)
      {
        double sdi = gsl_vector_get (sd, i);
        double sd2i = gsl_vector_get (sd2, i);
        gsl_test (fabs (sdi - sd2i) > eps * (fabs (sd2i) + scale),
                  "%s sdtau2 (%3lu)[%lu]: %22.18g   %22.18g", desc, N - 1, i, sdi, sd2i);
      }
  }

  gsl_matrix_complex_free (v2);
  gsl_matrix_complex_free (v);
  gsl_matrix_complex_free (u);
  gsl_vector_complex_free (tau1);
  gsl_vector_complex_free (tau2);
  gsl_vector_free (d);
  gsl_vector_free (sd);
  gsl_vector_free (d2);
  gsl_vector_free (sd2);

  return s;
}

int
test_bidiag_complex_unpack (gsl_rng * r)
{
  int f;
  int s = 0;
  size_t M, N;
  double acc;

  for (M = 1; M <= 50; M++)
    {
      for (N = 1; N <= M; N++)
        {
          gsl_matrix_complex *A = gsl_matrix_complex_alloc (M, N);

          acc = (M <  7) ? 1e2 * GSL_DBL_EPSILON
              : (M < 11) ? 1e4 * GSL_DBL_EPSILON
              : (M < 34) ? 1e5 * GSL_DBL_EPSILON
              :            1e6 * GSL_DBL_EPSILON;

          create_random_complex_matrix (A, r);
          s += test_bidiag_complex_unpack_eps (A, acc, "  bidiag_complex_unpack");

          gsl_matrix_complex_free (A);
        }
    }

  return s;
}
