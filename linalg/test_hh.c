/* linalg/test_hh.c
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
#include <math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>

static int
test_HH_left_eps(const gsl_matrix * m, const gsl_vector * w, const double eps, const char * desc)
{
  int s = 0;
  const size_t M = m->size1;
  const size_t N = m->size2;
  size_t i, j;

  double tau;

  gsl_matrix * A = gsl_matrix_alloc(M, N);
  gsl_vector * w1 = gsl_vector_alloc(w->size);
  gsl_vector * work = gsl_vector_alloc(N);

  gsl_matrix_memcpy(A, m);
  gsl_vector_memcpy(w1, w);

  tau = gsl_linalg_householder_transform(w1);
  gsl_vector_set(w1, 0, 1.0);

  s += gsl_linalg_householder_left(tau, w1, A, work);
  s += gsl_linalg_householder_left(tau, w1, A, work);

  {
    double scale = 0.0;

    for (i = 0; i < M; i++)
      {
        for (j = 0; j < N; j++)
          {
            double a = fabs(gsl_matrix_get(m, i, j));
            if (a > scale) scale = a;
          }
      }

    for (i = 0; i < M; i++)
      {
        for (j = 0; j < N; j++)
          {
            double hij = gsl_matrix_get(A, i, j);
            double mij = gsl_matrix_get(m, i, j);

            gsl_test(fabs(hij - mij) > eps * (fabs(mij) + scale),
                     "%s (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g",
                     desc, M, N, i, j, hij, mij);
          }
      }
  }

  gsl_vector_free(work);
  gsl_matrix_free(A);

  return s;
}

static int
test_HH_left(gsl_rng * r)
{
  int s = 0;
  size_t M, N;

  for (M = 1; M <= 50; ++M)
    {
      for (N = 1; N <= M; ++N)
        {
          gsl_matrix * A = gsl_matrix_alloc(M, N);
          gsl_vector * w = gsl_vector_alloc(M);

          create_random_matrix(A, r);
          create_random_vector(w, r);
          s += test_HH_left_eps(A, w, 1.0e4 * M * GSL_DBL_EPSILON, "HH_left random");

          gsl_vector_free(w);
          gsl_matrix_free(A);
        }
    }

  return s;
}


static int
test_HH_right_eps(const gsl_matrix * m, const gsl_vector * w, const double eps, const char * desc)
{
  int s = 0;
  const size_t M = m->size1;
  const size_t N = m->size2;
  size_t i, j;

  double tau;

  gsl_matrix * A = gsl_matrix_alloc(M, N);
  gsl_vector * w1 = gsl_vector_alloc(w->size);
  gsl_vector * work = gsl_vector_alloc(M);

  gsl_matrix_memcpy(A, m);
  gsl_vector_memcpy(w1, w);

  tau = gsl_linalg_householder_transform(w1);

  s += gsl_linalg_householder_right(tau, w1, A, work);
  s += gsl_linalg_householder_right(tau, w1, A, work);

  {
    double scale = 0.0;

    for (i = 0; i < M; i++)
      {
        for (j = 0; j < N; j++)
          {
            double a = fabs(gsl_matrix_get(m, i, j));
            if (a > scale) scale = a;
          }
      }

    for (i = 0; i < M; i++)
      {
        for (j = 0; j < N; j++)
          {
            double hij = gsl_matrix_get(A, i, j);
            double mij = gsl_matrix_get(m, i, j);

            gsl_test(fabs(hij - mij) > eps * (fabs(mij) + scale),
                     "%s (%3lu,%3lu)[%lu,%lu]: %22.18g   %22.18g",
                     desc, M, N, i, j, hij, mij);
          }
      }
  }

  gsl_vector_free(work);
  gsl_matrix_free(A);

  return s;
}

static int
test_HH_right(gsl_rng * r)
{
  int s = 0;
  size_t M, N;

  for (M = 1; M <= 50; ++M)
    {
      for (N = 1; N <= M; ++N)
        {
          gsl_matrix * A = gsl_matrix_alloc(M, N);
          gsl_vector * w = gsl_vector_alloc(N);

          create_random_matrix(A, r);
          create_random_vector(w, r);
          s += test_HH_right_eps(A, w, 1.0e4 * M * GSL_DBL_EPSILON, "HH_right random");

          gsl_vector_free(w);
          gsl_matrix_free(A);
        }
    }

  return s;
}

