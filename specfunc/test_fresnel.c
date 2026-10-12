/* specfunc/test_fresnel.c
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_sf.h>
#include "test_sf.h"

/* Reference values from mpmath (40 decimal digits). */
static const double fresnel_x[] = {
  0.0, 0.05, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5,
  4.0, 5.0, 6.0, 8.0, 10.0, 20.0, 50.0, 100.0
};
static const double fresnel_c_ref[] = {
   0.0000000000000000e+00,  4.9999922893770664e-02,  2.4975915035654320e-01,
   4.9234422587144638e-01,  7.7989340037682287e-01,  4.4526117603982152e-01,
   4.8825340607534073e-01,  4.5741300964177706e-01,  6.0572078929768558e-01,
   5.3257243502800089e-01,  4.9842603303817762e-01,  5.6363118870401219e-01,
   4.9953146785550112e-01,  4.9980218037719715e-01,  4.9989869420551575e-01,
   4.9998733497234438e-01,  4.9999918943072796e-01,  4.9999989867881789e-01
};
static const double fresnel_s_ref[] = {
   0.0000000000000000e+00,  6.5449774855615424e-05,  8.1756002357777566e-03,
   6.4732432859999273e-02,  4.3825914739035476e-01,  6.9750496008209306e-01,
   3.4341567836369824e-01,  6.1918175581959289e-01,  4.9631299896737502e-01,
   4.1524801197243755e-01,  4.2051575424692844e-01,  4.9919138191711687e-01,
   4.4696076123693029e-01,  4.6021421439301446e-01,  4.6816997858488224e-01,
   4.8408453592595391e-01,  4.9363380258593875e-01,  4.9681690114783755e-01
};
/* negative-argument references, from the odd symmetry */
static const double fresnel_xn[] = { -0.5, -1.5, -3.0, -7.0, -25.0 };
static const double fresnel_cn[] = {
  -4.9234422587144638e-01, -4.4526117603982152e-01, -6.0572078929768558e-01,
  -5.4546709254696979e-01, -5.1273238553977019e-01
};
static const double fresnel_sn[] = {
  -6.4732432859999273e-02, -6.9750496008209306e-01, -4.9631299896737502e-01,
  -4.9970478945344676e-01, -4.9999351546947618e-01
};

int
test_fresnel(void)
{
  gsl_sf_result r, rc, rs;
  int s = 0;
  size_t i;

  for (i = 0; i < sizeof(fresnel_x) / sizeof(fresnel_x[0]); i++)
    {
      double x = fresnel_x[i];

      TEST_SF(s, gsl_sf_fresnel_c_e, (x, &r), fresnel_c_ref[i], TEST_TOL2,
              GSL_SUCCESS);
      TEST_SF(s, gsl_sf_fresnel_s_e, (x, &r), fresnel_s_ref[i], TEST_TOL2,
              GSL_SUCCESS);
    }

  for (i = 0; i < sizeof(fresnel_xn) / sizeof(fresnel_xn[0]); i++)
    {
      double x = fresnel_xn[i];

      TEST_SF(s, gsl_sf_fresnel_c_e, (x, &r), fresnel_cn[i], TEST_TOL2,
              GSL_SUCCESS);
      TEST_SF(s, gsl_sf_fresnel_s_e, (x, &r), fresnel_sn[i], TEST_TOL2,
              GSL_SUCCESS);
    }

  /* the natural prototypes agree with the error-handling forms */
  for (i = 0; i < sizeof(fresnel_x) / sizeof(fresnel_x[0]); i++)
    {
      double x = fresnel_x[i];
      TEST_SF_VAL(s, gsl_sf_fresnel_c, (x), fresnel_c_ref[i], TEST_TOL2);
      TEST_SF_VAL(s, gsl_sf_fresnel_s, (x), fresnel_s_ref[i], TEST_TOL2);
    }

  /* symmetry: C and S are odd */
  {
    static const double sym_x[] = { 0.5, 1.5, 3.7, 9.0 };
    for (i = 0; i < sizeof(sym_x) / sizeof(sym_x[0]); i++)
      {
        double x = sym_x[i];
        gsl_test_abs(gsl_sf_fresnel_c(-x), -gsl_sf_fresnel_c(x), 1e-15,
                     "gsl_sf_fresnel_c(%g) is odd", x);
        gsl_test_abs(gsl_sf_fresnel_s(-x), -gsl_sf_fresnel_s(x), 1e-15,
                     "gsl_sf_fresnel_s(%g) is odd", x);
      }
  }

  /* exceptional returns: no finite limit at infinity, NaN propagates */
  TEST_SF_RETURN(s, gsl_sf_fresnel_c_e, (GSL_POSINF, &r), GSL_EDOM);
  TEST_SF_RETURN(s, gsl_sf_fresnel_s_e, (GSL_NEGINF, &r), GSL_EDOM);
  TEST_SF_RETURN(s, gsl_sf_fresnel_c_e, (GSL_NAN, &r), GSL_SUCCESS);
  gsl_test(!gsl_isnan(r.val), "gsl_sf_fresnel_c_e(NaN) returns NaN");

  return s;
}
