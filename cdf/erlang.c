/* cdf/erlang.c
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

#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cdf.h>

#include "error.h"

/* The Erlang distribution is the gamma distribution with an integer
   shape parameter k and rate lambda,

     p(x) dx = lambda^k x^(k-1) exp(-lambda x) / (k-1)! dx,   x > 0,

   i.e. a gamma variate with shape a = k and scale b = 1/lambda.  The
   cumulative distribution functions are therefore those of the gamma
   distribution with a = k and b = 1/lambda.

   See Savannah bug #66949. */

double
gsl_cdf_erlang_P (const double x, const unsigned int k, const double lambda)
{
  if (k < 1 || lambda <= 0.0)
    {
      CDF_ERROR ("k < 1 or lambda <= 0", GSL_EDOM);
    }

  return gsl_cdf_gamma_P (x, (double) k, 1.0 / lambda);
}

double
gsl_cdf_erlang_Q (const double x, const unsigned int k, const double lambda)
{
  if (k < 1 || lambda <= 0.0)
    {
      CDF_ERROR ("k < 1 or lambda <= 0", GSL_EDOM);
    }

  return gsl_cdf_gamma_Q (x, (double) k, 1.0 / lambda);
}
