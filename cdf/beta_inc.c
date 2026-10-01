/* specfunc/beta_inc.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

/* Author:  G. Jungman */
/* Modified for cdfs by Brian Gough, June 2003 */

#include <gsl/gsl_sf_gamma.h>

static double
beta_cont_frac (const double a, const double b, const double x,
                const double epsabs)
{
  /* 512 is not enough once the compiler is allowed to contract
   * multiplications and additions into fused operations: the terms round
   * slightly differently, the convergence test stalls, and the iteration
   * count for gsl_cdf_fdist_Q(3.786820954867802, 1, 100000) rises from 35
   * to 532, at which point this returned GSL_NAN and the caller returned
   * NaN instead of 0.0516604706391814.  See Savannah bug #64613.  A normal
   * evaluation converges in a few tens of steps, so the extra headroom
   * costs nothing. */
  const unsigned int max_iter = 4096;   /* control iterations      */
  const double cutoff = 2.0 * GSL_DBL_MIN;      /* control the zero cutoff */
  unsigned int iter_count = 0;
  double cf;

  /* standard initialization for continued fraction */
  double num_term = 1.0;
  double den_term = 1.0 - (a + b) * x / (a + 1.0);

  /* Clamp to the cutoff rather than poisoning the term: den_term is
   * divided by on the next line, so a NaN here propagates into cf and out
   * to the caller.  specfunc/beta_inc.c already does this. */
  if (fabs (den_term) < cutoff)
    den_term = cutoff;

  den_term = 1.0 / den_term;
  cf = den_term;

  while (iter_count < max_iter)
    {
      const int k = iter_count + 1;
      double coeff = k * (b - k) * x / (((a - 1.0) + 2 * k) * (a + 2 * k));
      double delta_frac;

      /* first step */
      den_term = 1.0 + coeff * den_term;
      num_term = 1.0 + coeff / num_term;

      if (fabs (den_term) < cutoff)
        den_term = cutoff;

      if (fabs (num_term) < cutoff)
        num_term = cutoff;

      den_term = 1.0 / den_term;

      delta_frac = den_term * num_term;
      cf *= delta_frac;

      coeff = -(a + k) * (a + b + k) * x / ((a + 2 * k) * (a + 2 * k + 1.0));

      /* second step */
      den_term = 1.0 + coeff * den_term;
      num_term = 1.0 + coeff / num_term;

      if (fabs (den_term) < cutoff)
        den_term = cutoff;

      if (fabs (num_term) < cutoff)
        num_term = cutoff;

      den_term = 1.0 / den_term;

      delta_frac = den_term * num_term;
      cf *= delta_frac;

      if (fabs (delta_frac - 1.0) < 2.0 * GSL_DBL_EPSILON)
        break;

      if (cf * fabs (delta_frac - 1.0) < epsabs)
        break;

      ++iter_count;
    }

  if (iter_count >= max_iter)
    return GSL_NAN;

  return cf;
}

/* The function beta_inc_AXPY(A,Y,a,b,x) computes A * beta_inc(a,b,x)
   + Y taking account of possible cancellations when using the
   hypergeometric transformation beta_inc(a,b,x)=1-beta_inc(b,a,1-x).

   It also adjusts the accuracy of beta_inc() to fit the overall
   absolute error when A*beta_inc is added to Y. (e.g. if Y >>
   A*beta_inc then the accuracy of beta_inc can be reduced) */



static double
beta_inc_AXPY (const double A, const double Y,
               const double a, const double b, const double x)
{
  if (x == 0.0)
    {
      return A * 0 + Y;
    }
  else if (x == 1.0)
    {
      return A * 1 + Y;
    }
  else if (a > 1e5 && b < 10 && x > a / (a + b))
    {
      /* Handle asymptotic regime, large a, small b, x > peak [AS 26.5.17] */
      double N = a + (b - 1.0) / 2.0;
      return A * gsl_sf_gamma_inc_Q (b, -N * log (x)) + Y;
    }
  else if (b > 1e5 && a < 10 && x < b / (a + b))
    {
      /* Handle asymptotic regime, small a, large b, x < peak [AS 26.5.17] */
      double N = b + (a - 1.0) / 2.0;
      return A * gsl_sf_gamma_inc_P (a, -N * log1p (-x)) + Y;
    }
  else
    {
      double ln_beta = gsl_sf_lnbeta (a, b);
      double ln_pre = -ln_beta + a * log (x) + b * log1p (-x);

      double prefactor = exp (ln_pre);

      if (x < (a + 1.0) / (a + b + 2.0))
        {
          /* Apply continued fraction directly. */
          double epsabs = fabs (Y / (A * prefactor / a)) * GSL_DBL_EPSILON;

          double cf = beta_cont_frac (a, b, x, epsabs);

          return A * (prefactor * cf / a) + Y;
        }
      else
        {
          /* Apply continued fraction after hypergeometric transformation. */
          double epsabs =
            fabs ((A + Y) / (A * prefactor / b)) * GSL_DBL_EPSILON;
          double cf;
          double term;

          /* With A == -Y the A + Y above cancels exactly and epsabs
           * becomes 0, which disables the second convergence test in
           * beta_cont_frac() and leaves only the relative one.  That is
           * harmless on its own - the relative test is the stronger of the
           * two - but cf can be far larger than 1 here, and then the
           * relative test needs many more iterations to fire.  Floor epsabs
           * at one epsilon, which keeps the absolute test meaningful for
           * cf of order 1 and leaves large cf to be settled by the relative
           * test.  See Savannah bug #64613. */
          if (epsabs < GSL_DBL_EPSILON)
            epsabs = GSL_DBL_EPSILON;

          cf = beta_cont_frac (b, a, 1.0 - x, epsabs);
          term = prefactor * cf / b;

          if (A == -Y)
            {
              return -A * term;
            }
          else
            {
              return A * (1 - term) + Y;
            }
        }
    }
}

/* Direct series evaluation for testing purposes only */

#if 0
static double
beta_series (const double a, const double b, const double x,
             const double epsabs)
{
  double f = x / (1 - x);
  double c = (b - 1) / (a + 1) * f;
  double s = 1;
  double n = 0;

  s += c;

  do
    {
      n++;
      c *= -f * (2 + n - b) / (2 + n + a);
      s += c;
    }
  while (n < 512 && fabs (c) > GSL_DBL_EPSILON * fabs (s) + epsabs);

  s /= (1 - x);

  return s;
}
#endif
