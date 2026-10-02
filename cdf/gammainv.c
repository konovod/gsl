/* cdf/gammainv.c
 * 
 * Copyright (C) 2003, 2007 Brian Gough
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
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>

/* The forward regularized incomplete gamma functions, evaluated through
   the routines that report a status rather than raising a GSL error.

   The inverse solver probes x well beyond the root while bracketing.
   Some of those probes land in the region x > 1e6, a <= x, a > x/5, where
   the large-x expansion of Q does not converge and GSL raises an error
   (even from the status-returning routine).  There the complementary
   continued fraction for P is used instead, which is stable. */

static int
gamma_inc_bad_large_x (const double a, const double x)
{
  return x > 1.0e6 && a > 0.2 * x;
}

static double
gamma_inc_P (const double a, const double x)
{
  gsl_sf_result r;

  if (x > a && !gamma_inc_bad_large_x (a, x)
      && gsl_sf_gamma_inc_Q_e (a, x, &r) == GSL_SUCCESS)
    return 1.0 - r.val;

  gsl_sf_gamma_inc_P_e (a, x, &r);
  return r.val;
}

static double
gamma_inc_Q (const double a, const double x)
{
  gsl_sf_result r;

  if (x < a)
    {
      if (gsl_sf_gamma_inc_P_e (a, x, &r) == GSL_SUCCESS)
        return 1.0 - r.val;
    }
  else if (gamma_inc_bad_large_x (a, x))
    {
      gsl_sf_gamma_inc_P_e (a, x, &r);
      return 1.0 - r.val;
    }

  if (gsl_sf_gamma_inc_Q_e (a, x, &r) == GSL_SUCCESS)
    return r.val;

  gsl_sf_gamma_inc_P_e (a, x, &r);
  return 1.0 - r.val;
}

/* Signed residual of the inverse problem, as a function of t = log(x):

     upper == 0 : H(t) = P(a,x) - target
     upper == 1 : H(t) = target - Q(a,x)

   P is increasing and Q is decreasing in x, so H is increasing in t in
   both cases, with H -> -target as t -> -infinity and H -> 1 - target as
   t -> +infinity.  The root therefore always admits a bracket in the
   representable range. */

static double
gamma_inverse_residual (const double t, const double a, const double target,
                        const int upper)
{
  double x = exp (t);

  if (!gsl_finite (x))
    x = GSL_DBL_MAX;

  return upper ? target - gamma_inc_Q (a, x) : gamma_inc_P (a, x) - target;
}

/* Solve P(a,x) = target (upper == 0) or Q(a,x) = target (upper == 1)
   for a scale of one and 0 < target <= 1/2.

   The root is bracketed by expanding around an initial estimate and the
   bracket is then closed with the Pegasus variant of regula falsi, which
   is guaranteed to converge for a monotone residual.  The iteration is
   carried out on t = log(x) so that the search spans the many orders of
   magnitude a gamma quantile can have when a is small. */

static double
gamma_invert (const double target, const double a, const int upper)
{
  const double tmin = log (GSL_DBL_MIN);
  const double tmax = log (GSL_DBL_MAX);
  double lg1 = gsl_sf_lngamma (a + 1.0);
  double ltarget = upper ? log1p (-target) : log (target);
  double t, tlo, thi, Hlo, Hhi, step;
  unsigned int n;

  if (lg1 + ltarget < 0.0)
    {
      /* Small root: P(a,x) ~ x^a / Gamma(a+1). */
      t = (lg1 + ltarget) / a;
    }
  else
    {
      double xg = upper ? gsl_cdf_ugaussian_Qinv (target)
                        : gsl_cdf_ugaussian_Pinv (target);
      double x0 = (xg < -0.5 * sqrt (a)) ? a : sqrt (a) * xg + a;
      t = (x0 > 0.0) ? log (x0) : log (a);
    }

  if (!(t > tmin))
    t = tmin;
  else if (t > tmax)
    t = tmax;

  /* Bracket the root. */
  step = 1.0;
  if (gamma_inverse_residual (t, a, target, upper) < 0.0)
    {
      tlo = t;
      thi = t + step;
      if (thi > tmax)
        thi = tmax;

      while (gamma_inverse_residual (thi, a, target, upper) < 0.0)
        {
          if (thi >= tmax)
            break;
          step *= 2.0;
          thi = t + step;
          if (thi > tmax)
            thi = tmax;
        }
    }
  else
    {
      thi = t;
      tlo = t - step;
      if (tlo < tmin)
        tlo = tmin;

      while (gamma_inverse_residual (tlo, a, target, upper) > 0.0)
        {
          if (tlo <= tmin)
            break;
          step *= 2.0;
          tlo = t - step;
          if (tlo < tmin)
            tlo = tmin;
        }
    }

  Hlo = gamma_inverse_residual (tlo, a, target, upper);
  Hhi = gamma_inverse_residual (thi, a, target, upper);

  /* No sign change down to the smallest normal number: the quantile
     underflows and zero is the only representable answer. */
  if (Hlo > 0.0)
    return 0.0;

  /* Pegasus iteration. */
  t = 0.5 * (tlo + thi);
  for (n = 0; n < 200; n++)
    {
      double H;

      if (thi - tlo <= GSL_DBL_EPSILON * GSL_MAX (1.0, fabs (t)))
        break;

      t = tlo + (thi - tlo) * (Hlo / (Hlo - Hhi));
      if (!gsl_finite (t) || t <= tlo || t >= thi)
        t = 0.5 * (tlo + thi);

      H = gamma_inverse_residual (t, a, target, upper);

      if (H < 0.0)
        {
          if (Hlo + H != 0.0)
            Hhi *= Hlo / (Hlo + H);
          tlo = t;
          Hlo = H;
        }
      else if (H > 0.0)
        {
          if (Hhi + H != 0.0)
            Hlo *= Hhi / (Hhi + H);
          thi = t;
          Hhi = H;
        }
      else
        {
          break;
        }
    }

  {
    double x = exp (t);

    /* The root is below the smallest normal number: the quantile
       underflows and the only representable answer is zero. */
    if (x < GSL_DBL_MIN)
      return 0.0;

    return x;
  }
}

double
gsl_cdf_gamma_Pinv (const double P, const double a, const double b)
{
  double x;

  if (P == 1.0)
    {
      return GSL_POSINF;
    }
  else if (P == 0.0)
    {
      return 0.0;
    }

  /* Work with whichever tail is no larger than one half, so that the
     inverted function is evaluated where it is relatively accurate. */
  x = (P <= 0.5) ? gamma_invert (P, a, 0) : gamma_invert (1.0 - P, a, 1);

  return b * x;
}

double
gsl_cdf_gamma_Qinv (const double Q, const double a, const double b)
{
  double x;

  if (Q == 1.0)
    {
      return 0.0;
    }
  else if (Q == 0.0)
    {
      return GSL_POSINF;
    }

  x = (Q <= 0.5) ? gamma_invert (Q, a, 1) : gamma_invert (1.0 - Q, a, 0);

  return b * x;
}
