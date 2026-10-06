/* cdf/binomialinv.c
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

/* The cumulative distribution function P(k) is non-decreasing in k and
   Q(k) is non-increasing, so the quantiles are found by bisection over
   the integers 0..n.  The result is the smallest k with P(k) >= P (or
   Q(k) <= Q).  See Savannah bug #66767. */

double
gsl_cdf_binomial_Pinv (const double P, const double p, const unsigned int n)
{
  unsigned int lo, hi;

  if (p < 0.0 || p > 1.0)
    {
      CDF_ERROR ("p < 0 or p > 1", GSL_EDOM);
    }

  if (P <= 0.0)
    {
      return 0.0;
    }
  if (P >= 1.0)
    {
      return (double) n;
    }

  lo = 0;
  hi = n;

  while (lo < hi)
    {
      unsigned int mid = lo + (hi - lo) / 2;

      if (gsl_cdf_binomial_P (mid, p, n) >= P)
        {
          hi = mid;
        }
      else
        {
          lo = mid + 1;
        }
    }

  return (double) lo;
}

double
gsl_cdf_binomial_Qinv (const double Q, const double p, const unsigned int n)
{
  unsigned int lo, hi;

  if (p < 0.0 || p > 1.0)
    {
      CDF_ERROR ("p < 0 or p > 1", GSL_EDOM);
    }

  if (Q <= 0.0)
    {
      return (double) n;
    }
  if (Q >= 1.0)
    {
      return 0.0;
    }

  lo = 0;
  hi = n;

  while (lo < hi)
    {
      unsigned int mid = lo + (hi - lo) / 2;

      if (gsl_cdf_binomial_Q (mid, p, n) <= Q)
        {
          hi = mid;
        }
      else
        {
          lo = mid + 1;
        }
    }

  return (double) lo;
}
