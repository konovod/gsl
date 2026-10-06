/* specfunc/elljac.c
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

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_elljac.h>

#include "eval.h"


/* GJ: See [Thompson, Atlas for Computing Mathematical Functions] */

/* BJG 2005-07: New algorithm based on Algorithm 5 from Numerische
   Mathematik 7, 78-90 (1965) "Numerical Calculation of Elliptic
   Integrals and Elliptic Functions" R. Bulirsch.

   Minor tweak is to avoid division by zero when sin(x u_l) = 0 by
   computing reflected values sn(K-u) cn(K-u) dn(K-u) and using
   transformation from Abramowitz & Stegun table 16.8 column "K-u"*/

int
gsl_sf_elljac_e(double u, double m, double * sn, double * cn, double * dn)
{
  if(fabs(m) > 1.0) {
    *sn = 0.0;
    *cn = 0.0;
    *dn = 0.0;
    GSL_ERROR ("|m| > 1.0", GSL_EDOM);
  }
  else if(fabs(m) < 2.0*GSL_DBL_EPSILON) {
    *sn = sin(u);
    *cn = cos(u);
    *dn = 1.0;
    return GSL_SUCCESS;
  }
  else if(fabs(m - 1.0) < 2.0*GSL_DBL_EPSILON) {
    *sn = tanh(u);
    *cn = 1.0/cosh(u);
    *dn = *cn;
    return GSL_SUCCESS;
  }
  else {
    int status = GSL_SUCCESS;
    const int N = 16;
    double mu[16];
    double nu[16];
    double c[16];
    double d[16];
    double sin_umu, cos_umu, t, r;
    int n = 0;

    mu[0] = 1.0;
    nu[0] = sqrt(1.0 - m);

    while( fabs(mu[n] - nu[n]) > 4.0 * GSL_DBL_EPSILON * fabs(mu[n]+nu[n])) {
      mu[n+1] = 0.5 * (mu[n] + nu[n]);
      nu[n+1] = sqrt(mu[n] * nu[n]);
      ++n;
      if(n >= N - 1) {
        status = GSL_EMAXITER;
        break;
      }
    }

    sin_umu = sin(u * mu[n]);
    cos_umu = cos(u * mu[n]);

    /* Since sin(u*mu(n)) can be zero we switch to computing sn(K-u),
       cn(K-u), dn(K-u) when |sin| < |cos| */

    if (fabs(sin_umu) < fabs(cos_umu))
      {
        t = sin_umu / cos_umu;
        
        c[n] = mu[n] * t;
        d[n] = 1.0;
        
        while(n > 0) {
          n--;
          c[n] = d[n+1] * c[n+1];
          r = (c[n+1] * c[n+1]) / mu[n+1];
          d[n] = (r + nu[n]) / (r + mu[n]);
          }
        
        *dn = sqrt(1.0-m) / d[n];
        *cn = (*dn) * GSL_SIGN(cos_umu) / gsl_hypot(1.0, c[n]);
        *sn = (*cn) * c[n] /sqrt(1.0-m);
      }
    else
      {
        t = cos_umu / sin_umu;
        
        c[n] = mu[n] * t;
        d[n] = 1.0;
        
        while(n > 0) {
          --n;
          c[n] = d[n+1] * c[n+1];
          r = (c[n+1] * c[n+1]) / mu[n+1];
          d[n] = (r + nu[n]) / (r + mu[n]);
        }
        
        *dn = d[n];
        *sn = GSL_SIGN(sin_umu) / gsl_hypot(1.0, c[n]);
        *cn = c[n] * (*sn);
      }
    
    return status;
  }
}


/* Inverse Jacobi elliptic functions.
 *
 * Each function returns the principal value u with 0 <= u <= 2 K(m)
 * such that the corresponding direct function equals the argument.
 * The direct functions are evaluated with the parameter m (not the
 * modulus k), matching gsl_sf_elljac_e().
 */

/* F(phi|m) = sin(phi) * RF(cos^2 phi, 1 - m sin^2 phi, 1), for a
   reduced amplitude |phi| <= pi/2.  Unlike gsl_sf_ellint_F_e() this
   takes the parameter m rather than the modulus k, so it also covers
   m < 0. */

static int
elljac_F_e(double sin_phi, double m, gsl_sf_result * result)
{
  const double sin2_phi = sin_phi * sin_phi;
  gsl_sf_result rf;
  const int status = gsl_sf_ellint_RF_e(1.0 - sin2_phi, 1.0 - m * sin2_phi,
                                        1.0, GSL_PREC_DOUBLE, &rf);

  result->val = sin_phi * rf.val;
  result->err = GSL_DBL_EPSILON * fabs(result->val) + fabs(sin_phi * rf.err);

  return status;
}

/* K(m) = RF(0, 1 - m, 1), real for every m <= 1 including m < 0. */

static int
elljac_K_e(double m, gsl_sf_result * result)
{
  return gsl_sf_ellint_RF_e(0.0, 1.0 - m, 1.0, GSL_PREC_DOUBLE, result);
}

int
gsl_sf_elljac_arcsn_e(double sn, double m, gsl_sf_result * result)
{
  if (fabs(m) > 1.0)
    {
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      GSL_ERROR ("|m| > 1.0", GSL_EDOM);
    }
  else if (fabs(sn) > 1.0)
    {
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      GSL_ERROR ("|sn| > 1.0", GSL_EDOM);
    }
  else
    {
      /* u = F(asin(sn)|m) */
      return elljac_F_e(sn, m, result);
    }
}

int
gsl_sf_elljac_arccn_e(double cn, double m, gsl_sf_result * result)
{
  if (fabs(m) > 1.0)
    {
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      GSL_ERROR ("|m| > 1.0", GSL_EDOM);
    }
  else if (fabs(cn) > 1.0)
    {
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      GSL_ERROR ("|cn| > 1.0", GSL_EDOM);
    }
  else
    {
      /* sin(phi) = sqrt(1 - cn^2), with phi = acos(cn) */
      const double sin_phi = sqrt((1.0 - cn) * (1.0 + cn));

      if (cn >= 0.0)
        {
          /* u = F(phi|m), phi <= pi/2 */
          return elljac_F_e(sin_phi, m, result);
        }
      else
        {
          /* u = 2 K(m) - F(pi - phi|m) */
          gsl_sf_result K, F;
          const int status_K = elljac_K_e(m, &K);
          const int status_F = elljac_F_e(sin_phi, m, &F);

          result->val = 2.0 * K.val - F.val;
          result->err = 2.0 * K.err + F.err;

          return GSL_ERROR_SELECT_2(status_K, status_F);
        }
    }
}

int
gsl_sf_elljac_arcdn_e(double dn, double m, gsl_sf_result * result)
{
  if (fabs(m) > 1.0)
    {
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      GSL_ERROR ("|m| > 1.0", GSL_EDOM);
    }
  else if (m == 0.0)
    {
      /* dn(u|0) = 1 for every u, so the inverse is not unique */
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      GSL_ERROR ("m = 0, inverse of dn is not unique", GSL_EDOM);
    }
  else
    {
      /* sin^2(phi) = (1 - dn^2)/m */
      const double sin2_phi = (1.0 - dn) * (1.0 + dn) / m;

      if (sin2_phi < 0.0 || sin2_phi > 1.0)
        {
          result->val = GSL_NAN;
          result->err = GSL_NAN;
          GSL_ERROR ("dn outside the range of dn(u|m)", GSL_EDOM);
        }
      else
        {
          return elljac_F_e(sqrt(sin2_phi), m, result);
        }
    }
}

double
gsl_sf_elljac_arcsn(double sn, double m)
{
  EVAL_RESULT(gsl_sf_elljac_arcsn_e(sn, m, &result));
}

double
gsl_sf_elljac_arccn(double cn, double m)
{
  EVAL_RESULT(gsl_sf_elljac_arccn_e(cn, m, &result));
}

double
gsl_sf_elljac_arcdn(double dn, double m)
{
  EVAL_RESULT(gsl_sf_elljac_arcdn_e(dn, m, &result));
}
