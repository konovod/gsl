/* randist/nakagami.c
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
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* The Nakagami distribution with shape mu > 0 and scale omega > 0 is
   defined by

     p(x) dx = 2 mu^mu / (Gamma(mu) omega^mu) x^(2mu-1)
               exp(-mu x^2 / omega) dx,   x >= 0.

   If Y is gamma distributed with shape mu and scale omega/mu then
   X = sqrt(Y) is Nakagami(mu, omega).  The distribution is normally
   restricted to mu >= 1/2.

   See Savannah bug #66816. */

double
gsl_ran_nakagami (const gsl_rng * r, const double mu, const double omega)
{
  double y = gsl_ran_gamma (r, mu, omega / mu);
  return sqrt (y);
}

double
gsl_ran_nakagami_pdf (const double x, const double mu, const double omega)
{
  if (x < 0.0)
    {
      return 0.0;
    }
  else if (x == 0.0 && mu == 0.5)
    {
      /* 0 * log(0) in the exponent below; the limit is finite. */
      return sqrt (2.0 / (M_PI * omega));
    }
  else
    {
      double lngamma = gsl_sf_lngamma (mu);
      return 2.0 * exp (mu * log (mu / omega) - lngamma
                        + (2.0 * mu - 1.0) * log (x)
                        - mu * x * x / omega);
    }
}
