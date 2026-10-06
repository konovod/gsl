/* integration/fejer.c
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

/* Fejer's first rule with n nodes, none of which are the endpoints.
 * The nodes are x_k = cos((2k+1) pi/(2n)) and the weights follow the
 * closed-form cosine series. */

#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

static int
fejer_check(const size_t n, const gsl_integration_fixed_params * params)
{
  if (n < 1)
    {
      GSL_ERROR("fejer quadrature requires n >= 1", GSL_EDOM);
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
fejer_nodes(const size_t n, double * x, double * w,
            gsl_integration_fixed_params * params)
{
  const double shft = 0.5 * (params->b + params->a);
  const double slp = 0.5 * (params->b - params->a);
  size_t k, j;

  for (k = 0; k < n; k++)
    {
      const double theta = (2.0 * k + 1.0) * M_PI / (2.0 * n);
      double sum = 0.0;

      x[k] = cos(theta);

      for (j = 1; j <= n / 2; j++)
        {
          sum += 1.0 / (4.0 * j * j - 1.0) * cos(2.0 * j * theta);
        }

      w[k] = 2.0 / (double) n * (1.0 - 2.0 * sum);
    }

  /* scale nodes and weights to [a,b] */
  for (k = 0; k < n; k++)
    {
      x[k] = shft + slp * x[k];
      w[k] *= slp;
    }

  return GSL_SUCCESS;
}

static const gsl_integration_fixed_type fejer_type =
{
  fejer_check,
  NULL,
  fejer_nodes
};

const gsl_integration_fixed_type *gsl_integration_fixed_fejer = &fejer_type;
