/* complex/test_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * Copyright (C) 2021 Patrick Alken
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

struct f
{
  char *name;
  double (*f) (gsl_complex z);
  double x;
  double y;
  double fx;
  double fy;
};

struct fz
{
  char *name;
  gsl_complex (*f) (gsl_complex z);
  double x;
  double y;
  double fx;
  double fy;
};

struct fzz
{
  char *name;
  gsl_complex (*f) (gsl_complex z1, gsl_complex z2);
  double x1;
  double y1;
  double x2;
  double y2;
  double fx;
  double fy;
};


struct freal
{
  char *name;
  gsl_complex (*f) (double x);
  double x;
  double fx;
  double fy;
};

struct fzreal
{
  char *name;
  gsl_complex (*f) (gsl_complex z, double a);
  double x;
  double y;
  double a;
  double fx;
  double fy;
};

#define FN(x) "gsl_complex_" #x, gsl_complex_ ## x
#define ARG(x,y) x, y
#define RES(x,y) x, y

struct f list[] =
{
#include "results1.h"
  {"", 0, 0, 0, 0, 0}
};


struct fz listz[] =
{
#include "results.h"
  {"", 0, 0, 0, 0, 0}
};

struct fzz listzz[] =
{
  {FN (pow), ARG(0.0,0.0), ARG(0.0,0.0), RES(1.0, 0.0)},
#include "results2.h"
  {"", 0, 0, 0, 0, 0, 0, 0}
};

struct freal listreal[] =
{
#include "results_real.h"
  {"", 0, 0, 0, 0}
};


struct fzreal listzreal[] =
{
#include "results_zreal.h"
  {"", 0, 0, 0, 0, 0, 0}
};


#ifndef TEST_FACTOR
#ifdef RELEASED
#define TEST_FACTOR 100.0  
#else
#define TEST_FACTOR 1.0
#endif
#endif

static int
FUNCTION (test, all) ()
{
  size_t i = 0;
#if defined(BASE_DOUBLE)
  const double tol = GSL_DBL_EPSILON;
#endif

  gsl_ieee_env_setup();

  for (i = 0 ; i < 10; i++) 
    {
      ATOMIC x = (i - 5.0) * 0.3;
      ATOMIC y = (i + 2.1) * 0.5;
      TYPE(gsl_complex) z = gsl_complex_rect(x, y);

      gsl_test_rel (GSL_REAL(z), x, tol, "gsl_complex_rect real part at (x=%g,y=%g)", x, y);
      gsl_test_rel (GSL_IMAG(z), y, tol, "gsl_complex_rect imag part at (x=%g,y=%g)", x, y);

      GSL_REAL(z) = y;
      gsl_test_rel (GSL_REAL(z), y, tol, "assignment real part (%g)", y);

      GSL_IMAG(z) = x;
      gsl_test_rel (GSL_IMAG(z), x, tol, "assignment imag part (%g)", x);
    }

  for (i = 0 ; i < 10; i++) 
    {
      ATOMIC r = (i - 5.0) * 0.3 ;
      ATOMIC t = 2.0 * M_PI * i / 5 ;
      ATOMIC x = r * cos(t), y = r * sin(t) ;
      TYPE(gsl_complex) z = gsl_complex_polar (r, t) ;
      gsl_test_rel (GSL_REAL(z), x, tol, "gsl_complex_polar real part at (r=%g,t=%g)", r, t);
      gsl_test_rel (GSL_IMAG(z), y, tol, "gsl_complex_polar imag part at (r=%g,t=%g)", r, t);
    }

    i = 0;

  while (list[i].f)
    {
      struct f t = list[i];
      TYPE(gsl_complex) z = gsl_complex_rect (t.x, t.y);
      ATOMIC f = (t.f) (z);
      gsl_test_rel (f, t.fx, tol, "%s at (%g,%g)", t.name, t.x, t.y);
      i++;
    }

  i = 0;

  while (listz[i].f)
    {
      struct fz t = listz[i];
      TYPE(gsl_complex) z = gsl_complex_rect (t.x, t.y);
      TYPE(gsl_complex) fz = (t.f) (z);
      ATOMIC fx = GSL_REAL (fz), fy = GSL_IMAG (fz);

#ifdef DEBUG
      printf("x = "); gsl_ieee_fprintf_double (stdout, &t.x); printf("\n");
      printf("y = "); gsl_ieee_fprintf_double (stdout, &t.y); printf("\n");
      printf("fx = "); gsl_ieee_fprintf_double (stdout, &fx); printf("\n");
      printf("ex = "); gsl_ieee_fprintf_double (stdout, &t.fx); printf("\n");
      printf("fy = "); gsl_ieee_fprintf_double (stdout, &fy); printf("\n");
      printf("ey = "); gsl_ieee_fprintf_double (stdout, &t.fy); printf("\n");
#endif

      gsl_test_rel (fx, t.fx, 10.0 * tol, "%s real part at (%g,%g)", t.name, t.x, t.y);
      gsl_test_rel (fy, t.fy, 10.0 * tol, "%s imag part at (%g,%g)", t.name, t.x, t.y);
      i++;
    }

  i = 0;

  while (listzz[i].f)
    {
      struct fzz t = listzz[i];
      TYPE(gsl_complex) z1 = gsl_complex_rect (t.x1, t.y1);
      TYPE(gsl_complex) z2 = gsl_complex_rect (t.x2, t.y2);
      TYPE(gsl_complex) fz = (t.f) (z1, z2);
      ATOMIC fx = GSL_REAL (fz), fy = GSL_IMAG (fz);

#ifdef DEBUG
      printf("x1 = "); gsl_ieee_fprintf_double (stdout, &t.x1); printf("\n");
      printf("y1 = "); gsl_ieee_fprintf_double (stdout, &t.y1); printf("\n");
      printf("x2 = "); gsl_ieee_fprintf_double (stdout, &t.x2); printf("\n");
      printf("y2 = "); gsl_ieee_fprintf_double (stdout, &t.y2); printf("\n");
      printf("fx = "); gsl_ieee_fprintf_double (stdout, &fx); printf("\n");
      printf("ex = "); gsl_ieee_fprintf_double (stdout, &t.fx); printf("\n");
      printf("fy = "); gsl_ieee_fprintf_double (stdout, &fy); printf("\n");
      printf("ey = "); gsl_ieee_fprintf_double (stdout, &t.fy); printf("\n");
#endif

      gsl_test_rel (fx, t.fx, 1.0e3 * tol, "%s real part at (%g,%g;%g,%g)", t.name, t.x1, t.y1, t.x2, t.y2);
      gsl_test_rel (fy, t.fy, 1.0e3 * tol, "%s imag part at (%g,%g;%g,%g)", t.name, t.x1, t.y1, t.x2, t.y2);
      i++;
    }

  i = 0;

  while (listreal[i].f)
    {
      struct freal t = listreal[i];
      TYPE(gsl_complex) fz = (t.f) (t.x);
      ATOMIC fx = GSL_REAL (fz), fy = GSL_IMAG (fz);

#ifdef DEBUG
      printf("x = "); gsl_ieee_fprintf_double (stdout, &t.x); printf("\n");
      printf("fx = "); gsl_ieee_fprintf_double (stdout, &fx); printf("\n");
      printf("ex = "); gsl_ieee_fprintf_double (stdout, &t.fx); printf("\n");
      printf("fy = "); gsl_ieee_fprintf_double (stdout, &fy); printf("\n");
      printf("ey = "); gsl_ieee_fprintf_double (stdout, &t.fy); printf("\n");
#endif

      gsl_test_rel (fx, t.fx, tol, "%s real part at (%g,0)", t.name, t.x);
      gsl_test_rel (fy, t.fy, tol, "%s imag part at (%g,0)", t.name, t.x);
      i++;
    }

  i = 0;

  while (listzreal[i].f)
    {
      struct fzreal t = listzreal[i];
      TYPE(gsl_complex) z = gsl_complex_rect (t.x, t.y);
      TYPE(gsl_complex) fz = (t.f) (z, t.a);
      ATOMIC fx = GSL_REAL (fz), fy = GSL_IMAG (fz);

#ifdef DEBUG
      printf("x = "); gsl_ieee_fprintf_double (stdout, &t.x); printf("\n");
      printf("y = "); gsl_ieee_fprintf_double (stdout, &t.y); printf("\n");
      printf("a = "); gsl_ieee_fprintf_double (stdout, &t.a); printf("\n");
      printf("fx = "); gsl_ieee_fprintf_double (stdout, &fx); printf("\n");
      printf("ex = "); gsl_ieee_fprintf_double (stdout, &t.fx); printf("\n");
      printf("fy = "); gsl_ieee_fprintf_double (stdout, &fy); printf("\n");
      printf("ey = "); gsl_ieee_fprintf_double (stdout, &t.fy); printf("\n");
#endif

      gsl_test_rel (fx, t.fx, 10.0 * tol, "%s real part at (%g,0)", t.name, t.x);
      gsl_test_rel (fy, t.fy, 10.0 * tol, "%s imag part at (%g,0)", t.name, t.x);
      i++;
    }

  /* gsl_complex_arccosh() for a real argument must agree with
   * gsl_complex_arccosh_real(), in particular it must return +0.0 (not
   * -0.0) for the real part on the positive real axis (bug #47345).  The
   * tables above cannot see this because they compare magnitudes. */

  {
    const double rx[] = { -10.0, -2.0, -1.0, 0.5, 1.0, 2.0 };
    size_t k;

    for (k = 0; k < sizeof (rx) / sizeof (rx[0]); ++k)
      {
        TYPE (gsl_complex) fz = gsl_complex_arccosh (gsl_complex_rect (rx[k], 0.0));
        TYPE (gsl_complex) rz = gsl_complex_arccosh_real (rx[k]);
        ATOMIC fx = GSL_REAL (fz), fy = GSL_IMAG (fz);

        gsl_test_rel (fx, GSL_REAL (rz), tol,
                      "gsl_complex_arccosh real part at (%g,0)", rx[k]);
        gsl_test_rel (fy, GSL_IMAG (rz), tol,
                      "gsl_complex_arccosh imag part at (%g,0)", rx[k]);
        gsl_test (signbit (fx) ? 1 : 0,
                  "gsl_complex_arccosh real part sign at (%g,0)", rx[k]);
      }
  }

  /* GSL_SET_COMPLEX() must assign the two components directly rather
   * than evaluate x + I*y: the complex arithmetic in that expression
   * rounds away the sign of a zero component and computes 0*inf for a
   * non-finite part, turning it into NaN (bug #65912).  The tables above
   * compare magnitudes and cannot see this. */

  {
    gsl_complex z;
    gsl_complex_float zf;
    gsl_complex_long_double zl;
    const long double eps = 1.0L + LDBL_EPSILON;

    GSL_SET_COMPLEX (&z, -0.0, 0.0);
    gsl_test (!signbit (GSL_REAL (z)),
              "GSL_SET_COMPLEX real sign at (-0.0, +0.0)");
    gsl_test (GSL_IMAG (z) != 0.0 || signbit (GSL_IMAG (z)),
              "GSL_SET_COMPLEX imag sign at (-0.0, +0.0)");

    GSL_SET_COMPLEX (&z, 0.0, -0.0);
    gsl_test (!signbit (GSL_IMAG (z)),
              "GSL_SET_COMPLEX imag sign at (+0.0, -0.0)");

    GSL_SET_COMPLEX (&z, 0.0, GSL_POSINF);
    gsl_test (GSL_REAL (z) != 0.0 || signbit (GSL_REAL (z)),
              "GSL_SET_COMPLEX real part at (+0.0, +inf)");
    gsl_test (gsl_isinf (GSL_IMAG (z)) <= 0,
              "GSL_SET_COMPLEX imag part at (+0.0, +inf)");

    GSL_SET_COMPLEX (&z, GSL_NEGINF, 1.0);
    gsl_test (gsl_isinf (GSL_REAL (z)) >= 0 || GSL_IMAG (z) != 1.0,
              "GSL_SET_COMPLEX at (-inf, 1.0)");

    GSL_SET_COMPLEX (&z, GSL_NAN, 2.0);
    gsl_test (!gsl_isnan (GSL_REAL (z)) || GSL_IMAG (z) != 2.0,
              "GSL_SET_COMPLEX at (nan, 2.0)");

    GSL_SET_COMPLEX (&zf, -0.0f, 0.0f);
    gsl_test (!signbit (GSL_REAL (zf)),
              "GSL_SET_COMPLEX float real sign at (-0.0f, +0.0f)");

    GSL_SET_COMPLEX (&zl, -0.0L, 0.0L);
    gsl_test (!signbit (GSL_REAL (zl)),
              "GSL_SET_COMPLEX long double real sign at (-0.0L, +0.0L)");
    GSL_SET_COMPLEX (&zl, eps, 0.0L);
    gsl_test (GSL_REAL (zl) != eps,
              "GSL_SET_COMPLEX long double precision at (1+eps, 0.0)");
  }

  return 0;
}
