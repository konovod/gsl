/* integration/lobatto.c
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

/* Gauss-Lobatto quadrature with n nodes, including both endpoints.
 * The interior nodes are the roots of P'_{n-1}, which are the Gauss
 * nodes for the Jacobi weight (1-x^2); the weights follow from
 *   w_i = 2 / (n (n-1) P_{n-1}(x_i)^2)
 * on [-1,1].  The rule is exact for polynomials of degree <= 2n-3. */

#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

static int
lobatto_check(const size_t n, const gsl_integration_fixed_params * params)
{
  if (n < 2)
    {
      GSL_ERROR("lobatto quadrature requires n >= 2", GSL_EDOM);
    }
  else if (fabs(params->b - params->a) <= GSL_DBL_EPSILON)
    {
      GSL_ERROR("|b - a| too small", GSL_EDOM);
    }
  else
    {
      return GSL_SUCCESS;
    }
}

/* Legendre polynomial P_n(x) by the three-term recurrence */
static double
legendre_p(int n, double x)
{
  double p0 = 1.0;
  double p1 = x;
  int k;

  if (n == 0)
    return 1.0;

  for (k = 2; k <= n; k++)
    {
      double p2 = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / k;
      p0 = p1;
      p1 = p2;
    }

  return p1;
}

static int
lobatto_nodes(const size_t n, double * x, double * w, gsl_integration_fixed_params * params)
{
  const size_t m = n - 2;
  const double shft = 0.5 * (params->b + params->a);
  const double slp = 0.5 * (params->b - params->a);
  size_t i, j;

  x[0] = -1.0;
  x[n - 1] = 1.0;

  if (m > 0)
    {
      gsl_matrix *J = gsl_matrix_calloc(m, m);
      gsl_vector *eval = gsl_vector_alloc(m);
      gsl_eigen_symm_workspace *eworkspace = gsl_eigen_symm_alloc(m);

      for (i = 0; i + 1 < m; i++)
        {
          const double k = (double) (i + 1);
          const double b = sqrt(k * (k + 2.0) / ((2.0 * k + 1.0) * (2.0 * k + 3.0)));
          gsl_matrix_set(J, i, i + 1, b);
          gsl_matrix_set(J, i + 1, i, b);
        }

      gsl_eigen_symm(J, eval, eworkspace);

      for (i = 0; i < m; i++)
        x[i + 1] = gsl_vector_get(eval, i);

      /* sort the interior nodes */
      for (i = 1; i < n - 1; i++)
        {
          for (j = i + 1; j < n - 1; j++)
            {
              if (x[j] < x[i])
                {
                  double t = x[i];
                  x[i] = x[j];
                  x[j] = t;
                }
            }
        }

      gsl_eigen_symm_free(eworkspace);
      gsl_vector_free(eval);
      gsl_matrix_free(J);
    }

  /* weights on [-1,1] */
  for (i = 0; i < n; i++)
    {
      double p = legendre_p((int) (n - 1), x[i]);
      w[i] = 2.0 / (n * (n - 1.0) * p * p);
    }

  /* scale nodes and weights to [a,b] */
  for (i = 0; i < n; i++)
    {
      x[i] = shft + slp * x[i];
      w[i] *= slp;
    }

  return GSL_SUCCESS;
}

static const gsl_integration_fixed_type lobatto_type =
{
  lobatto_check,
  NULL,
  lobatto_nodes
};

const gsl_integration_fixed_type *gsl_integration_fixed_lobatto = &lobatto_type;
