/* cdf/poissoninv.c
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
#include <limits.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cdf.h>

#include "error.h"

/* The Poisson support is unbounded, so first bracket the quantile by
   doubling an upper limit until the tail condition holds, then bisect.
   The result is the smallest k with P(k) >= P (or Q(k) <= Q).  Because
   the cumulative distribution function takes an unsigned int argument,
   the quantile saturates at UINT_MAX.

   See Savannah bug #66775. */

double
gsl_cdf_poisson_Pinv (const double P, const double mu)
{
  unsigned int lo, hi;

  if (mu <= 0.0)
    {
      CDF_ERROR ("mu <= 0", GSL_EDOM);
    }

  if (P <= 0.0)
    {
      return 0.0;
    }
  if (P >= 1.0)
    {
      return GSL_POSINF;
    }

  hi = 1;
  while (hi < UINT_MAX && gsl_cdf_poisson_P (hi, mu) < P)
    {
      hi = (hi > UINT_MAX / 2) ? UINT_MAX : 2 * hi;
    }

  lo = 0;

  while (lo < hi)
    {
      unsigned int mid = lo + (hi - lo) / 2;

      if (gsl_cdf_poisson_P (mid, mu) >= P)
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
gsl_cdf_poisson_Qinv (const double Q, const double mu)
{
  unsigned int lo, hi;

  if (mu <= 0.0)
    {
      CDF_ERROR ("mu <= 0", GSL_EDOM);
    }

  if (Q <= 0.0)
    {
      return GSL_POSINF;
    }
  if (Q >= 1.0)
    {
      return 0.0;
    }

  hi = 1;
  while (hi < UINT_MAX && gsl_cdf_poisson_Q (hi, mu) > Q)
    {
      hi = (hi > UINT_MAX / 2) ? UINT_MAX : 2 * hi;
    }

  lo = 0;

  while (lo < hi)
    {
      unsigned int mid = lo + (hi - lo) / 2;

      if (gsl_cdf_poisson_Q (mid, mu) <= Q)
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
