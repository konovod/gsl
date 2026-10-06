/* specfunc/zeta_complex.c
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

/* Complex Riemann and Hurwitz zeta functions, and the complex eta
 * function.
 *
 * For Re(s) >= 0 the Hurwitz zeta function is evaluated with the
 * Euler-Maclaurin summation formula [Moshier, p. 400], generalised to
 * complex s.  For Re(s) < 0 the Riemann zeta function uses the
 * reflection formula
 *
 *   zeta(s) = 2^s pi^(s-1) sin(pi s/2) Gamma(1-s) zeta(1-s),
 *
 * where zeta(1-s) again has Re(1-s) > 1 and is computed with the
 * Euler-Maclaurin formula.  The trivial zeros at negative even integers
 * are returned exactly.
 */

#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_trig.h>

/* B_{2j}/(2j)!, j = 1..14 */
static const double zeta_em_c[15] = {
  0.0,
  0.083333333333333333333333333333,
  -0.00138888888888888888888888888889,
  0.000033068783068783068783068783069,
  -8.2671957671957671957671957672e-07,
  2.0876756987868098979210090321e-08,
  -5.2841901386874931848476822022e-10,
  1.3382536530684678832826980975e-11,
  -3.3896802963225828668301953912e-13,
  8.5860620562778445641359054504e-15,
  -2.1748686985580618730415164239e-16,
  5.5090028283602295152026526089e-18,
  -1.3954464685812523340707686264e-19,
  3.5347070396294674716932299778e-21,
  -8.9535174270375468504026113181e-23
};

/* Euler-Maclaurin evaluation of the complex Hurwitz zeta function
 * zeta(s, q) for q > 0, Re(s) >= 0 and s != 1.  The number of leading
 * terms is scaled with |s|, which is what the critical strip needs. */

static int
complex_hzeta_em(double x, double y, double q, gsl_complex * z, double * err)
{
  const int J = 14;
  const int K = 12 + (int) ceil(fabs(x)) + (int) ceil(fabs(y));
  const gsl_complex s = gsl_complex_rect(x, y);
  const gsl_complex ms = gsl_complex_negative(s);
  const gsl_complex Kq = gsl_complex_rect(K + q, 0.0);
  const gsl_complex invKq = gsl_complex_inverse(Kq);
  const gsl_complex invKq2 = gsl_complex_mul(invKq, invKq);
  gsl_complex pmax, pcp, scp, ans;
  double amag, last = 0.0;
  int k, j;

  pmax = gsl_complex_pow(Kq, ms);
  pcp = gsl_complex_mul(pmax, invKq);
  ans = gsl_complex_mul(pmax,
                        gsl_complex_add(gsl_complex_div(Kq, gsl_complex_sub(s, GSL_COMPLEX_ONE)),
                                        gsl_complex_rect(0.5, 0.0)));

  for (k = 0; k < K; k++)
    {
      gsl_complex t = gsl_complex_pow(gsl_complex_rect(k + q, 0.0), ms);
      ans = gsl_complex_add(ans, t);
    }

  scp = s;
  for (j = 0; j < J; j++)
    {
      gsl_complex delta = gsl_complex_mul_real(gsl_complex_mul(scp, pcp), zeta_em_c[j + 1]);

      ans = gsl_complex_add(ans, delta);
      last = gsl_complex_abs(delta);

      {
        gsl_complex a = gsl_complex_add_real(s, 2 * j + 1);
        gsl_complex b = gsl_complex_add_real(s, 2 * j + 2);
        scp = gsl_complex_mul(scp, gsl_complex_mul(a, b));
      }
      pcp = gsl_complex_mul(pcp, invKq2);
    }

  amag = gsl_complex_abs(ans);
  *err = 8.0 * (K + J + 4) * GSL_DBL_EPSILON * amag + last;
  *z = ans;

  return GSL_SUCCESS;
}

/* Reflection formula for Re(s) < 0. */

static int
complex_zeta_reflect(double x, double y, gsl_complex * z, double * err)
{
  const gsl_complex s = gsl_complex_rect(x, y);
  gsl_complex zw, gamma1ms, sinhalf, pref;
  double zwerr, gmag, amag;
  gsl_sf_result lnr, arg, sr, si;
  int status;

  status = complex_hzeta_em(1.0 - x, -y, 1.0, &zw, &zwerr);
  if (status)
    return status;

  status = gsl_sf_lngamma_complex_e(1.0 - x, -y, &lnr, &arg);
  if (status)
    return status;

  gmag = exp(lnr.val);
  gamma1ms = gsl_complex_rect(gmag * cos(arg.val), gmag * sin(arg.val));

  status = gsl_sf_complex_sin_e(0.5 * M_PI * x, 0.5 * M_PI * y, &sr, &si);
  if (status)
    return status;

  sinhalf = gsl_complex_rect(sr.val, si.val);

  pref = gsl_complex_mul(gsl_complex_pow(gsl_complex_rect(2.0, 0.0), s),
                         gsl_complex_pow(gsl_complex_rect(M_PI, 0.0),
                                         gsl_complex_sub(s, GSL_COMPLEX_ONE)));
  pref = gsl_complex_mul(pref, sinhalf);
  pref = gsl_complex_mul(pref, gamma1ms);

  *z = gsl_complex_mul(pref, zw);

  amag = gsl_complex_abs(*z);
  *err = 16.0 * (fabs(x) + fabs(y) + 4.0) * GSL_DBL_EPSILON * amag
    + gsl_complex_abs(pref) * zwerr
    + fabs(lnr.err) * amag + fabs(arg.err) * amag;

  return GSL_SUCCESS;
}

static void
zeta_domain_error(gsl_sf_result * result_re, gsl_sf_result * result_im)
{
  result_re->val = GSL_NAN;
  result_re->err = GSL_NAN;
  result_im->val = GSL_NAN;
  result_im->err = GSL_NAN;
}

int
gsl_sf_complex_zeta_e(double x, double y,
                      gsl_sf_result * result_re, gsl_sf_result * result_im)
{
  gsl_complex z;
  double err;
  int status;

  if (x == 1.0 && y == 0.0)
    {
      zeta_domain_error(result_re, result_im);
      GSL_ERROR ("zeta has a pole at s = 1", GSL_EDOM);
    }

  /* trivial zeros at the negative even integers */
  if (y == 0.0 && x < 0.0 && x == floor(x) && fmod(x, 2.0) == 0.0)
    {
      result_re->val = 0.0;
      result_re->err = 0.0;
      result_im->val = 0.0;
      result_im->err = 0.0;
      return GSL_SUCCESS;
    }

  if (x >= 0.0)
    status = complex_hzeta_em(x, y, 1.0, &z, &err);
  else
    status = complex_zeta_reflect(x, y, &z, &err);

  if (status)
    {
      zeta_domain_error(result_re, result_im);
      return status;
    }

  result_re->val = GSL_REAL(z);
  result_im->val = GSL_IMAG(z);
  result_re->err = err + GSL_DBL_EPSILON * fabs(result_re->val);
  result_im->err = err + GSL_DBL_EPSILON * fabs(result_im->val);

  return GSL_SUCCESS;
}

int
gsl_sf_complex_hzeta_e(double x, double y, double q,
                       gsl_sf_result * result_re, gsl_sf_result * result_im)
{
  gsl_complex z;
  double err;

  if (q <= 0.0)
    {
      zeta_domain_error(result_re, result_im);
      GSL_ERROR ("q must be positive", GSL_EDOM);
    }
  else if (x <= 1.0)
    {
      /* the real gsl_sf_hzeta_e is defined for s > 1 only */
      zeta_domain_error(result_re, result_im);
      GSL_ERROR ("Re(s) must be greater than 1", GSL_EDOM);
    }

  complex_hzeta_em(x, y, q, &z, &err);

  result_re->val = GSL_REAL(z);
  result_im->val = GSL_IMAG(z);
  result_re->err = err + GSL_DBL_EPSILON * fabs(result_re->val);
  result_im->err = err + GSL_DBL_EPSILON * fabs(result_im->val);

  return GSL_SUCCESS;
}

int
gsl_sf_complex_eta_e(double x, double y,
                     gsl_sf_result * result_re, gsl_sf_result * result_im)
{
  if (x == 1.0 && y == 0.0)
    {
      /* eta(1) = log 2 */
      result_re->val = M_LN2;
      result_re->err = GSL_DBL_EPSILON * M_LN2;
      result_im->val = 0.0;
      result_im->err = 0.0;
      return GSL_SUCCESS;
    }
  else
    {
      gsl_sf_result zr, zi;
      gsl_complex s, factor, z, e;
      double fmag;
      int status = gsl_sf_complex_zeta_e(x, y, &zr, &zi);

      if (status)
        {
          zeta_domain_error(result_re, result_im);
          return status;
        }

      /* eta(s) = (1 - 2^(1-s)) zeta(s) */
      s = gsl_complex_rect(x, y);
      factor = gsl_complex_sub(GSL_COMPLEX_ONE,
                               gsl_complex_pow(gsl_complex_rect(2.0, 0.0),
                                               gsl_complex_sub(GSL_COMPLEX_ONE, s)));
      z = gsl_complex_rect(zr.val, zi.val);
      e = gsl_complex_mul(factor, z);

      fmag = gsl_complex_abs(factor);
      result_re->val = GSL_REAL(e);
      result_im->val = GSL_IMAG(e);
      result_re->err = fmag * zr.err + GSL_DBL_EPSILON * fabs(result_re->val);
      result_im->err = fmag * zi.err + GSL_DBL_EPSILON * fabs(result_im->val);

      return GSL_SUCCESS;
    }
}
