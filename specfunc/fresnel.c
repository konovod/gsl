/* specfunc/fresnel.c
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

/* The Fresnel integrals
 *
 *   S(x) = Integral_0^x sin(pi t^2 / 2) dt
 *   C(x) = Integral_0^x cos(pi t^2 / 2) dt
 *
 * both odd in x.  For 0 <= x <= 3/2 the Maclaurin series is summed
 * directly.  For x > 3/2 the auxiliary functions f and g are used,
 *
 *   C(x) = 1/2 + f(x) sin(theta) - g(x) cos(theta)
 *   S(x) = 1/2 - f(x) cos(theta) - g(x) sin(theta),   theta = pi x^2 / 2,
 *
 * where P(w) = pi x f(x) and Q(w) = pi x g(x), w = 1/(pi x^2), are
 * smooth.  For 3/2 < x <= 6 (w in [1/(36 pi), 1/(2.25 pi)]) P and Q are
 * evaluated from Chebyshev series; for x > 6 the asymptotic series
 *
 *   P(w) = Sum_{m>=0} (-1)^m (4m-1)!! w^(2m),
 *   Q(w) = Sum_{m>=0} (-1)^m (4m+1)!! w^(2m+1)
 *
 * is summed to its smallest term (about 1e-16 for x >= 6).
 *
 * The Chebyshev coefficients below were generated with a 60-digit
 * reference and give an absolute error below 3e-16 in P and Q on the
 * stated interval, so C and S are accurate to a few 1e-15 for all
 * arguments up to the point where the phase theta = pi x^2 / 2 can no
 * longer be represented; beyond that the reported error grows with
 * theta, as it must.
 */

#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_fresnel.h>

#include "error.h"
#include "chebyshev.h"
#include "cheb_eval.c"

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* Chebyshev series for P(w) and Q(w), w in [1/(36 pi), 1/(2.25 pi)],
   i.e. 3/2 <= x <= 6.  The argument is mapped to [-1,1] by cheb_eval_e
   through the interval limits stored in the cheb_series. */

static double fp_data[25] = {
   1.96464209597920991e+00, -2.11175635204337679e-02, -3.10882145498337118e-03,
   5.32888810920020394e-04, -3.81714347493584049e-05, -4.48478338332769354e-06,
   2.32403553165567628e-06, -5.07865779292447650e-07,  6.16275470185251256e-08,
   3.84813511352086020e-09, -4.93714978020152674e-09,  1.76294653523080445e-09,
  -4.23353973488526074e-10,  6.31799954242335477e-11,  2.80427259123239281e-12,
  -6.44810536471669038e-12,  2.98774669778080142e-12, -9.68341593995321255e-13,
   2.36878229133717811e-13, -3.57267676620086869e-14, -3.68048334749179347e-15,
   5.75324097636166090e-15, -2.94148974487640684e-15,  1.10942178939193286e-15,
  -3.26018268554680196e-16
};

static cheb_series fp_cs = {
  fp_data,
  24,
  8.84194128288307533e-03, 1.41471060526129205e-01,
  12
};

static double gp_data[25] = {
   1.33390971622754179e-01,  5.46099334797397820e-02, -3.42190578169660961e-03,
  -8.52840032281207519e-05,  6.91906742851440705e-05, -1.27058894503927135e-05,
   1.09249672782785422e-06,  1.62921153891352867e-07, -1.01432166919102399e-07,
   2.76934766701423698e-08, -4.83639204805630952e-09,  2.49713910202835754e-10,
   2.16630704534444391e-10, -1.14895674887978797e-10,  3.70858955959258709e-11,
  -8.69036532480454585e-12,  1.20598681965615150e-12,  1.48513233802384235e-13,
  -1.89990422054053576e-13,  8.88911336788242673e-14, -3.05025550651316223e-14,
   8.18605563295600020e-15, -1.51516341698596002e-15,  1.54746852240262444e-17,
   1.24697882945424735e-16
};

static cheb_series gp_cs = {
  gp_data,
  24,
  8.84194128288307533e-03, 1.41471060526129205e-01,
  12
};

/* Maclaurin series for C and S, for 0 <= x <= 3/2.  Both series are
   alternating with monotonically decreasing terms over this range, so
   there is no cancellation. */
static int
fresnel_series(const double x,
               gsl_sf_result * result_c, gsl_sf_result * result_s)
{
  const double a2 = M_PI_2;
  const double x2 = x * x;
  const double x4 = x2 * x2;
  double c = 0.0, s = 0.0;
  double tc = x;                /* n = 0 term of C */
  double ts = a2 * x2 * x / 3.0; /* n = 0 term of S */
  double c_last = tc, s_last = ts;
  int n;

  for (n = 0; n < 60; n++)
    {
      c += tc;
      s += ts;

      /* tc_{n+1}/tc_n = -a2^2 x^4 (4n+1) / ((4n+5)(2n+1)(2n+2))
         ts_{n+1}/ts_n = -a2^2 x^4 (4n+3) / ((4n+7)(2n+2)(2n+3)) */
      tc *= -a2 * a2 * x4 * (4.0 * n + 1.0)
            / ((4.0 * n + 5.0) * (2.0 * n + 1.0) * (2.0 * n + 2.0));
      ts *= -a2 * a2 * x4 * (4.0 * n + 3.0)
            / ((4.0 * n + 7.0) * (2.0 * n + 2.0) * (2.0 * n + 3.0));

      c_last = fabs(tc);
      s_last = fabs(ts);

      if (tc == 0.0 && ts == 0.0)
        break;
      if (fabs(tc) <= GSL_DBL_EPSILON * fabs(c)
          && fabs(ts) <= GSL_DBL_EPSILON * fabs(s))
        break;
    }

  result_c->val = c;
  result_c->err = c_last + 4.0 * GSL_DBL_EPSILON * fabs(c);
  result_s->val = s;
  result_s->err = s_last + 4.0 * GSL_DBL_EPSILON * fabs(s);

  return GSL_SUCCESS;
}

/* P(w) and Q(w) for x > 6 by the asymptotic series, truncated at the
   smallest term. */
static void
fresnel_asymptotic(const double x, gsl_sf_result * result_p,
                   gsl_sf_result * result_q)
{
  const double w = 1.0 / (M_PI * x * x);
  const double w2 = w * w;
  double p = 1.0, q = w;
  double tp = 1.0, tq = w;
  double perr = 0.0, qerr = 0.0;
  int m;

  for (m = 0; m < 200; m++)
    {
      const double np = -tp * (4.0 * m + 1.0) * (4.0 * m + 3.0) * w2;
      const double nq = -tq * (4.0 * m + 3.0) * (4.0 * m + 5.0) * w2;

      if (fabs(np) >= fabs(tp) || fabs(nq) >= fabs(tq))
        {
          perr = fabs(np);
          qerr = fabs(nq);
          break;
        }

      p += np;
      q += nq;
      tp = np;
      tq = nq;
    }

  result_p->val = p;
  result_p->err = perr + 2.0 * GSL_DBL_EPSILON * fabs(p);
  result_q->val = q;
  result_q->err = qerr + 2.0 * GSL_DBL_EPSILON * fabs(q);
}

/* C and S for x > 3/2 via the auxiliary functions. */
static int
fresnel_aux(const double x,
            gsl_sf_result * result_c, gsl_sf_result * result_s)
{
  const double w = 1.0 / (M_PI * x * x);
  const double theta = M_PI * x * x / 2.0;
  const double invpix = 1.0 / (M_PI * x);
  gsl_sf_result rp, rq;
  double f, g, ferr, gerr, st, ct;
  double phase_err;

  if (w >= 8.84194128288307533e-03)
    {
      /* 3/2 < x <= 6: Chebyshev series */
      cheb_eval_e(&fp_cs, w, &rp);
      cheb_eval_e(&gp_cs, w, &rq);
    }
  else
    {
      fresnel_asymptotic(x, &rp, &rq);
    }

  f = rp.val * invpix;
  g = rq.val * invpix;
  ferr = rp.err * invpix + 2.0 * GSL_DBL_EPSILON * fabs(f);
  gerr = rq.err * invpix + 2.0 * GSL_DBL_EPSILON * fabs(g);

  st = sin(theta);
  ct = cos(theta);

  /* theta is formed in double precision, so its rounding is an
     irreducible part of the error for large x */
  phase_err = (fabs(f) + fabs(g)) * fabs(theta) * GSL_DBL_EPSILON;

  result_c->val = 0.5 + f * st - g * ct;
  result_c->err = ferr + gerr + phase_err
                  + 2.0 * GSL_DBL_EPSILON * fabs(result_c->val);

  result_s->val = 0.5 - f * ct - g * st;
  result_s->err = ferr + gerr + phase_err
                  + 2.0 * GSL_DBL_EPSILON * fabs(result_s->val);

  return GSL_SUCCESS;
}

/* C(|x|) and S(|x|) */
static int
fresnel_cs(const double x,
           gsl_sf_result * result_c, gsl_sf_result * result_s)
{
  if (x <= 1.5)
    return fresnel_series(x, result_c, result_s);
  else
    return fresnel_aux(x, result_c, result_s);
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_fresnel_c_e(const double x, gsl_sf_result * result)
{
  gsl_sf_result rc, rs;
  int status;

  if (gsl_isnan(x))
    {
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      return GSL_SUCCESS;
    }
  else if (gsl_isinf(x))
    {
      DOMAIN_ERROR(result);
    }

  status = fresnel_cs(fabs(x), &rc, &rs);
  if (status != GSL_SUCCESS)
    return status;

  result->val = (x < 0.0) ? -rc.val : rc.val;
  result->err = rc.err;
  return GSL_SUCCESS;
}

int
gsl_sf_fresnel_s_e(const double x, gsl_sf_result * result)
{
  gsl_sf_result rc, rs;
  int status;

  if (gsl_isnan(x))
    {
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      return GSL_SUCCESS;
    }
  else if (gsl_isinf(x))
    {
      DOMAIN_ERROR(result);
    }

  status = fresnel_cs(fabs(x), &rc, &rs);
  if (status != GSL_SUCCESS)
    return status;

  result->val = (x < 0.0) ? -rs.val : rs.val;
  result->err = rs.err;
  return GSL_SUCCESS;
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double gsl_sf_fresnel_c(const double x)
{
  EVAL_RESULT(gsl_sf_fresnel_c_e(x, &result));
}

double gsl_sf_fresnel_s(const double x)
{
  EVAL_RESULT(gsl_sf_fresnel_s_e(x, &result));
}
