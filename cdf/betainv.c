/* cdf/betainv.c
 *
 * Copyright (C) 2004 Free Software Foundation, Inc.
 * Copyright (C) 2006, 2007 Brian Gough
 * Written by Jason H. Stover.
 * Modified for GSL by Brian Gough
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
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * Invert the Beta distribution. 
 * 
 * References:
 *
 * Roger W. Abernathy and Robert P. Smith. "Applying Series Expansion
 * to the Inverse Beta Distribution to Find Percentiles of the
 * F-Distribution," ACM Transactions on Mathematical Software, volume
 * 19, number 4, December 1993, pages 474-480.
 *
 * G.W. Hill and A.W. Davis. "Generalized asymptotic expansions of a
 * Cornish-Fisher type," Annals of Mathematical Statistics, volume 39,
 * number 8, August 1968, pages 1264-1273.
 */

#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#include "error.h"

/* First approximation to the lower quantile, P(x;a,b) = P, for P <= 1/2.
   For small P the leading term P(x) ~ x^a / (a B(a,b)) is inverted and a
   correction factor is applied; otherwise the mean is used.  This is a
   starting point only: the solver below no longer relies on its accuracy. */

static double
beta_initial (const double P, const double a, const double b)
{
  double mean = a / (a + b);
  double x;

  if (P < 0.1)
    {
      /* small x */

      double lg_ab = gsl_sf_lngamma (a + b);
      double lg_a = gsl_sf_lngamma (a);
      double lg_b = gsl_sf_lngamma (b);

      double lx = (log (a) + lg_a + lg_b - lg_ab + log (P)) / a;
      if (lx <= 0)
        {
          x = exp (lx);                     /* first approximation */
          x *= pow (1 - x, -(b - 1) / a);   /* second approximation */
        }
      else
        {
          x = mean;
        }

      if (x > mean)
        x = mean;
    }
  else
    {
      /* Use expected value as first guess */
      x = mean;
    }

  if (!(x > 0.0) || !(x < 1.0))
    x = 0.5;

  return x;
}

/* Invert the regularized incomplete beta function.

   beta_inverse solves P(x;a,b) = target when upper is zero and
   Q(x;a,b) = target when upper is non-zero, for 0 < target <= 1/2.  The
   iteration is carried out on t = logit(x), so that the many orders of
   magnitude a quantile can span when a or b is small are all reached.

   The residual is increasing in t in both cases:

     upper == 0 :  F(t) = P(x(t);a,b) - target,   F -> -target,     1-target
     upper == 1 :  F(t) = target - Q(x(t);a,b),   F -> target - 1, target

   so [log(DBL_MIN), -log(DBL_MIN)] is always a bracket.  A Newton step
   from the density is taken when it stays inside the bracket, otherwise
   the bracket is halved; convergence is guaranteed for every a,b > 0.

   The earlier implementation bisected to an absolute tolerance in P and
   then ran an unsafeguarded Newton iteration in x.  For large a and b
   the initial approximation can already satisfy that tolerance while
   still being far from the root, and for a << 1 the quantile can lie
   hundreds of decades below the mean, so the routine reported failure
   (NaN) for much of the parameter space. */

static double
beta_inverse (const double target, const double a, const double b,
              const int upper, double x)
{
  const double tmin = log (GSL_DBL_MIN);
  const double tmax = -tmin;
  const double ltarget = upper ? 1.0 - target : target;
  double t, tlo, thi;
  int converged = 0;
  unsigned int n;

  /* The quantile is below the smallest normal number: zero is the only
     representable answer. */
  if (gsl_cdf_beta_P (GSL_DBL_MIN, a, b) > ltarget)
    return 0.0;

  if (!(x > 0.0))
    x = GSL_DBL_MIN;
  else if (!(x < 1.0))
    x = 1.0;

  t = log (x / (1.0 - x));
  if (!(t > tmin))
    t = tmin;
  else if (t > tmax)
    t = tmax;
  x = 1.0 / (1.0 + exp (-t));

  tlo = tmin;
  thi = tmax;

  for (n = 0; n < 200; n++)
    {
      const double xlo = 1.0 / (1.0 + exp (-tlo));
      const double xhi = 1.0 / (1.0 + exp (-thi));
      const double F = upper ? target - gsl_cdf_beta_Q (x, a, b)
                             : gsl_cdf_beta_P (x, a, b) - target;

      /* Stop when the bracket has closed in x.  Near x = 1 (or 0) many
         values of t map to the same double x, so the t-bracket alone can
         stall while x is already the representable root. */
      if (F == 0.0 || xhi <= xlo
          || thi - tlo <= 4.0 * GSL_DBL_EPSILON * GSL_MAX (1.0, fabs (t)))
        {
          converged = 1;
          break;
        }

      if (F < 0.0)
        tlo = t;
      else
        thi = t;

      {
        const double dpdt = gsl_ran_beta_pdf (x, a, b) * x * (1.0 - x);
        const double tn = t - F / dpdt;
        const double dtol = 4.0 * GSL_DBL_EPSILON * GSL_MAX (1.0, fabs (t));
        int bisect = 1;

        if (gsl_finite (tn) && tn > tlo && tn < thi)
          {
            /* A Newton step smaller than the resolution of t means the
               root has been reached, even if the opposite end of the
               bracket was set many iterations ago. */
            if (fabs (tn - t) <= dtol)
              {
                t = tn;
                converged = 1;
                break;
              }

            /* Otherwise take the step only if it moves x: once x is at
               the resolution of a double a further step in t cannot
               change it, so fall back to bisection. */
            if (1.0 / (1.0 + exp (-tn)) != x)
              {
                t = tn;
                x = 1.0 / (1.0 + exp (-t));
                bisect = 0;
              }
          }

        if (bisect)
          {
            const double tm = 0.5 * (tlo + thi);

            /* The bracket cannot be narrowed further. */
            if (tm <= tlo || tm >= thi)
              {
                converged = 1;
                break;
              }

            t = tm;
            x = 1.0 / (1.0 + exp (-t));
          }
      }
    }

  if (!converged)
    {
      GSL_ERROR_VAL ("inverse failed to converge", GSL_EFAILED, GSL_NAN);
    }

  x = 1.0 / (1.0 + exp (-t));
  if (x < GSL_DBL_MIN)
    return 0.0;

  return x;
}

double
gsl_cdf_beta_Pinv (const double P, const double a, const double b)
{
  if (P < 0.0 || P > 1.0)
    {
      CDF_ERROR ("P must be in range 0 < P < 1", GSL_EDOM);
    }

  if (a < 0.0)
    {
      CDF_ERROR ("a < 0", GSL_EDOM);
    }

  if (b < 0.0)
    {
      CDF_ERROR ("b < 0", GSL_EDOM);
    }

  if (P == 0.0)
    {
      return 0.0;
    }

  if (P == 1.0)
    {
      return 1.0;
    }

  /* Work with whichever tail is no larger than one half, so that the
     inverted function is evaluated where it is relatively accurate. */
  if (P > 0.5)
    {
      return beta_inverse (1.0 - P, a, b, 1,
                           1.0 - beta_initial (1.0 - P, b, a));
    }

  return beta_inverse (P, a, b, 0, beta_initial (P, a, b));
}

double
gsl_cdf_beta_Qinv (const double Q, const double a, const double b)
{
  if (Q < 0.0 || Q > 1.0)
    {
      CDF_ERROR ("Q must be inside range 0 < Q < 1", GSL_EDOM);
    }

  if (a < 0.0)
    {
      CDF_ERROR ("a < 0", GSL_EDOM);
    }

  if (b < 0.0)
    {
      CDF_ERROR ("b < 0", GSL_EDOM);
    }

  if (Q == 0.0)
    {
      return 1.0;
    }

  if (Q == 1.0)
    {
      return 0.0;
    }

  /* Work with whichever tail is no larger than one half. */
  if (Q > 0.5)
    {
      return beta_inverse (1.0 - Q, a, b, 0,
                           beta_initial (1.0 - Q, a, b));
    }

  return beta_inverse (Q, a, b, 1, 1.0 - beta_initial (Q, b, a));
}
