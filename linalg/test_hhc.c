/* linalg/test_hhc.c
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
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>

static int
test_HH_complex_left_eps(const gsl_matrix_complex * m, const gsl_vector_complex * w, const double eps, const char * desc)
{
  int s = 0;
  const size_t M = m->size1;
  const size_t N = m->size2;
  size_t i, j;

  gsl_complex tau;

  gsl_matrix_complex * A = gsl_matrix_complex_alloc(M, N);
  gsl_vector_complex * w1 = gsl_vector_complex_alloc(w->size);
  gsl_vector_complex * work = gsl_vector_complex_alloc(N);

  gsl_matrix_complex_memcpy(A, m);
  gsl_vector_complex_memcpy(w1, w);

  tau = gsl_linalg_complex_householder_transform(w1);

  s += gsl_linalg_complex_householder_left(tau, w1, A, work);
  tau = gsl_complex_conjugate(tau);
  s += gsl_linalg_complex_householder_left(tau, w1, A, work);

  {
    double scale = 0.0;

    for (i = 0; i < M; i++)
      {
        for (j = 0; j < N; j++)
          {
            double a = gsl_complex_abs(gsl_matrix_complex_get(m, i, j));
            if (a > scale) scale = a;
          }
      }

    for (i = 0; i < M; i++)
      {
        for (j = 0; j < N; j++)
          {
            gsl_complex aij = gsl_matrix_complex_get(A, i, j);
            gsl_complex mij = gsl_matrix_complex_get(m, i, j);
            double d = gsl_complex_abs(gsl_complex_sub(aij, mij));
            int foo = (d > eps * (gsl_complex_abs(mij) + scale));
            if (foo)
              {
                printf("%s (%3lu,%3lu)[%lu,%lu]: %22.18g%+22.18gi   %22.18g%+22.18gi\n",
                       desc, M, N, i, j, GSL_REAL(aij), GSL_IMAG(aij),
                       GSL_REAL(mij), GSL_IMAG(mij));
              }
            s += foo;
          }
      }
  }

  gsl_vector_complex_free(work);
  gsl_matrix_complex_free(A);

  return s;
}

static int
test_HH_complex_left(gsl_rng * r)
{
  int s = 0;
  size_t M, N;

  for (M = 1; M <= 50; ++M)
    {
      for (N = 1; N <= M; ++N)
        {
          gsl_matrix_complex * A = gsl_matrix_complex_alloc(M, N);
          gsl_vector_complex * w = gsl_vector_complex_alloc(M);

          create_random_complex_matrix(A, r);
          create_random_complex_vector(w, r);
          s += test_HH_complex_left_eps(A, w, 1.0e5 * M * GSL_DBL_EPSILON, "HH_complex_left random");

          gsl_vector_complex_free(w);
          gsl_matrix_complex_free(A);
        }
    }

  return s;
}


static int
test_HH_complex_right_eps(const gsl_matrix_complex * m, const gsl_vector_complex * w, const double eps, const char * desc)
{
  int s = 0;
  const size_t M = m->size1;
  const size_t N = m->size2;
  size_t i, j;

  gsl_complex tau;

  gsl_matrix_complex * A = gsl_matrix_complex_alloc(M, N);
  gsl_vector_complex * w1 = gsl_vector_complex_alloc(w->size);
  gsl_vector_complex * work = gsl_vector_complex_alloc(M);

  gsl_matrix_complex_memcpy(A, m);
  gsl_vector_complex_memcpy(w1, w);

  tau = gsl_linalg_complex_householder_transform(w1);

  s += gsl_linalg_complex_householder_right(tau, w1, A, work);
  tau = gsl_complex_conjugate(tau);
  s += gsl_linalg_complex_householder_right(tau, w1, A, work);

  {
    double scale = 0.0;

    for (i = 0; i < M; i++)
      {
        for (j = 0; j < N; j++)
          {
            double a = gsl_complex_abs(gsl_matrix_complex_get(m, i, j));
            if (a > scale) scale = a;
          }
      }

    for (i = 0; i < M; i++)
      {
        for (j = 0; j < N; j++)
          {
            gsl_complex aij = gsl_matrix_complex_get(A, i, j);
            gsl_complex mij = gsl_matrix_complex_get(m, i, j);
            double d = gsl_complex_abs(gsl_complex_sub(aij, mij));
            int foo = (d > eps * (gsl_complex_abs(mij) + scale));
            if (foo)
              {
                printf("%s (%3lu,%3lu)[%lu,%lu]: %22.18g%+22.18gi   %22.18g%+22.18gi\n",
                       desc, M, N, i, j, GSL_REAL(aij), GSL_IMAG(aij),
                       GSL_REAL(mij), GSL_IMAG(mij));
              }
            s += foo;
          }
      }
  }

  gsl_vector_complex_free(work);
  gsl_matrix_complex_free(A);

  return s;
}

static int
test_HH_complex_right(gsl_rng * r)
{
  int s = 0;
  size_t M, N;

  for (M = 1; M <= 50; ++M)
    {
      for (N = 1; N <= M; ++N)
        {
          gsl_matrix_complex * A = gsl_matrix_complex_alloc(M, N);
          gsl_vector_complex * w = gsl_vector_complex_alloc(N);

          create_random_complex_matrix(A, r);
          create_random_complex_vector(w, r);
          s += test_HH_complex_right_eps(A, w, 1.0e5 * M * GSL_DBL_EPSILON, "HH_complex_right random");

          gsl_vector_complex_free(w);
          gsl_matrix_complex_free(A);
        }
    }

  return s;
}

/* Complex Givens rotations: the rotation must zero the second component
 * and preserve the norm of (a,b), with c real and |c|^2 + |s|^2 = 1. */

static int
test_givens_complex_eps(const gsl_complex a, const gsl_complex b, const char * desc)
{
  int s = 0;
  double c;
  double norm2 = gsl_complex_abs2(a) + gsl_complex_abs2(b);
  gsl_complex sn, r, z;
  gsl_vector_complex * v = gsl_vector_complex_alloc(2);

  gsl_linalg_complex_givens(a, b, &c, &sn);

  gsl_vector_complex_set(v, 0, a);
  gsl_vector_complex_set(v, 1, b);
  gsl_linalg_complex_givens_gv(v, 0, 1, c, sn);

  r = gsl_vector_complex_get(v, 0);
  z = gsl_vector_complex_get(v, 1);

  gsl_test(gsl_complex_abs(z) > 1.0e-14 * (sqrt(norm2) + 1.0),
           "%s: second component %g not zero", desc, gsl_complex_abs(z));
  gsl_test(fabs(gsl_complex_abs2(r) - norm2) > 1.0e-13 * (norm2 + 1.0),
           "%s: norm %g, expected %g", desc, gsl_complex_abs2(r), norm2);
  gsl_test(fabs(c * c + gsl_complex_abs2(sn) - 1.0) > 1.0e-13,
           "%s: c^2 + |s|^2 = %g", desc, c * c + gsl_complex_abs2(sn));

  gsl_vector_complex_free(v);

  return s;
}

static int
test_givens_complex(gsl_rng * r)
{
  int s = 0;
  size_t i;

  for (i = 0; i < 1000; ++i)
    {
      gsl_complex a = gsl_complex_rect(2.0 * (gsl_rng_uniform(r) - 0.5),
                                       2.0 * (gsl_rng_uniform(r) - 0.5));
      gsl_complex b = gsl_complex_rect(2.0 * (gsl_rng_uniform(r) - 0.5),
                                       2.0 * (gsl_rng_uniform(r) - 0.5));
      s += test_givens_complex_eps(a, b, "givens complex random");
    }

  s += test_givens_complex_eps(gsl_complex_rect(0.0, 0.0), gsl_complex_rect(1.0, 2.0),
                               "givens complex a = 0");
  s += test_givens_complex_eps(gsl_complex_rect(3.0, -1.0), gsl_complex_rect(0.0, 0.0),
                               "givens complex b = 0");

  return s;
}

