/* integration/clenshaw_curtis.c
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

/* Clenshaw-Curtis quadrature with n nodes, including both endpoints.
 * The nodes are x_k = cos(k pi/(n-1)) and the weights follow the
 * closed-form cosine series (Trefethen, "Is Gauss quadrature better
 * than Clenshaw-Curtis?", SIAM Rev. 50 (2008) 67-87). */

#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

static int
clenshaw_curtis_check(const size_t n, const gsl_integration_fixed_params * params)
{
  if (n < 2)
    {
      GSL_ERROR("clenshaw-curtis quadrature requires n >= 2", GSL_EDOM);
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
clenshaw_curtis_nodes(const size_t n, double * x, double * w,
                      gsl_integration_fixed_params * params)
{
  const double shft = 0.5 * (params->b + params->a);
  const double slp = 0.5 * (params->b - params->a);
  const size_t nm1 = n - 1;
  size_t k, j;

  for (k = 0; k < n; k++)
    {
      const double theta = M_PI * (double) k / (double) nm1;
      const double ck = (k == 0 || k == nm1) ? 1.0 : 2.0;
      double sum = 0.0;

      x[k] = cos(theta);

      for (j = 1; j <= nm1 / 2; j++)
        {
          const double bj = (2 * j == (int) nm1) ? 1.0 : 2.0;
          sum += bj / (4.0 * j * j - 1.0) * cos(2.0 * j * theta);
        }

      w[k] = ck / (double) nm1 * (1.0 - sum);
    }

  /* scale nodes and weights to [a,b] */
  for (k = 0; k < n; k++)
    {
      x[k] = shft + slp * x[k];
      w[k] *= slp;
    }

  return GSL_SUCCESS;
}

static const gsl_integration_fixed_type clenshaw_curtis_type =
{
  clenshaw_curtis_check,
  NULL,
  clenshaw_curtis_nodes
};

const gsl_integration_fixed_type *gsl_integration_fixed_clenshaw_curtis = &clenshaw_curtis_type;
