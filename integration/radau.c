/* integration/radau.c
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

/* Gauss-Radau quadrature with n nodes, including the right endpoint b.
 * The nodes are the eigenvalues of the Legendre Jacobi matrix with the
 * last diagonal element replaced by n/(2n-1), and the weights follow
 * from the Golub-Welsch formula.  The rule is exact for polynomials of
 * degree <= 2n-2. */

#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

static int
radau_check(const size_t n, const gsl_integration_fixed_params * params)
{
  if (n < 1)
    {
      GSL_ERROR("radau quadrature requires n >= 1", GSL_EDOM);
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

static int
radau_init(const size_t n, double * diag, double * subdiag, gsl_integration_fixed_params * params)
{
  size_t i;

  /* Legendre Jacobi matrix, with the last diagonal modified so that the
     right endpoint is a node */
  for (i = 1; i <= n; i++)
    {
      diag[i - 1] = 0.0;
      subdiag[i - 1] = (double) i / sqrt(4.0 * i * i - 1.0);
    }

  diag[n - 1] = (double) n / (2.0 * n - 1.0);

  params->zemu = 2.0;
  params->shft = 0.5 * (params->b + params->a);
  params->slp = 0.5 * (params->b - params->a);
  params->al = 0.0;
  params->be = 0.0;

  return GSL_SUCCESS;
}

static const gsl_integration_fixed_type radau_type =
{
  radau_check,
  radau_init,
  NULL
};

const gsl_integration_fixed_type *gsl_integration_fixed_radau = &radau_type;
