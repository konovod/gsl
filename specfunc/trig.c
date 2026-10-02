/* specfunc/trig.c
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
#include <stdint.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_trig.h>

#include "error.h"

#include "chebyshev.h"
#include "cheb_eval.c"

/* sinh(x) series
 * double-precision for |x| < 1.0
 */
inline
static
int
sinh_series(const double x, double * result)
{
  const double y = x*x;
  const double c0 = 1.0/6.0;
  const double c1 = 1.0/120.0;
  const double c2 = 1.0/5040.0;
  const double c3 = 1.0/362880.0;
  const double c4 = 1.0/39916800.0;
  const double c5 = 1.0/6227020800.0;
  const double c6 = 1.0/1307674368000.0;
  const double c7 = 1.0/355687428096000.0;
  *result = x*(1.0 + y*(c0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*c7))))))));
  return GSL_SUCCESS;
}


/* cosh(x)-1 series
 * double-precision for |x| < 1.0
 */
inline
static
int
cosh_m1_series(const double x, double * result)
{
  const double y = x*x;
  const double c0 = 0.5;
  const double c1 = 1.0/24.0;
  const double c2 = 1.0/720.0;
  const double c3 = 1.0/40320.0;
  const double c4 = 1.0/3628800.0;
  const double c5 = 1.0/479001600.0;
  const double c6 = 1.0/87178291200.0;
  const double c7 = 1.0/20922789888000.0;
  const double c8 = 1.0/6402373705728000.0;
  *result = y*(c0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*(c7+y*c8))))))));
  return GSL_SUCCESS;
}


/* Chebyshev expansion for f(t) = sinc((t+1)/2), -1 < t < 1
 */
static double sinc_data[17] = {
  1.133648177811747875422,
 -0.532677564732557348781,
 -0.068293048346633177859,
  0.033403684226353715020,
  0.001485679893925747818,
 -0.000734421305768455295,
 -0.000016837282388837229,
  0.000008359950146618018,
  0.000000117382095601192,
 -0.000000058413665922724,
 -0.000000000554763755743,
  0.000000000276434190426,
  0.000000000001895374892,
 -0.000000000000945237101,
 -0.000000000000004900690,
  0.000000000000002445383,
  0.000000000000000009925
};
static cheb_series sinc_cs = {
  sinc_data,
  16,
  -1, 1,
  10
};


/* Chebyshev expansion for f(t) = g((t+1)Pi/8), -1<t<1
 * g(x) = (sin(x)/x - 1)/(x*x)
 */
static double sin_data[12] = {
  -0.3295190160663511504173,
   0.0025374284671667991990,
   0.0006261928782647355874,
  -4.6495547521854042157541e-06,
  -5.6917531549379706526677e-07,
   3.7283335140973803627866e-09,
   3.0267376484747473727186e-10,
  -1.7400875016436622322022e-12,
  -1.0554678305790849834462e-13,
   5.3701981409132410797062e-16,
   2.5984137983099020336115e-17,
  -1.1821555255364833468288e-19
};
static cheb_series sin_cs = {
  sin_data,
  11,
  -1, 1,
  11
};

/* Chebyshev expansion for f(t) = g((t+1)Pi/8), -1<t<1
 * g(x) = (2(cos(x) - 1)/(x^2) + 1) / x^2
 */
static double cos_data[11] = {
  0.165391825637921473505668118136,
 -0.00084852883845000173671196530195,
 -0.000210086507222940730213625768083,
  1.16582269619760204299639757584e-6,
  1.43319375856259870334412701165e-7,
 -7.4770883429007141617951330184e-10,
 -6.0969994944584252706997438007e-11,
  2.90748249201909353949854872638e-13,
  1.77126739876261435667156490461e-14,
 -7.6896421502815579078577263149e-17,
 -3.7363121133079412079201377318e-18
};
static cheb_series cos_cs = {
  cos_data,
  10,
  -1, 1,
  10
};


/*-*-*-*-*-*-*-*-*-*-*-* Argument reduction *-*-*-*-*-*-*-*-*-*-*/

/* The reduction below is adapted from fdlibm's __rem_pio2_large
 * (k_rem_pio2.c).  fdlibm carries the following notice:
 *
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunSoft, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */

/* Above this |x| the three-term reduction below loses the reduced angle
 * (its error in the quotient `y = x/(pi/4)` reaches a full ulp of the
 * result around here), so the exact reduction is used instead.  The
 * value matches the threshold at which fdlibm switches to its
 * Payne-Hanek path; there the fast reduction is still exact. */
#define REM_PIO2_THRESHOLD (2097152.0)   /* 2^21 */

/* 2/pi after the binary point, in 24-bit chunks.  A double has at most
 * 1024 exponent bits, so at most (1024 - 3)/24 + 4 = 47 chunks are
 * needed; the table (from fdlibm) carries 66, which covers the whole
 * double range. */
static const int32_t two_over_pi[] = {
  0xA2F983, 0x6E4E44, 0x1529FC, 0x2757D1, 0xF534DD, 0xC0DB62,
  0x95993C, 0x439041, 0xFE5163, 0xABDEBB, 0xC561B7, 0x246E3A,
  0x424DD2, 0xE00649, 0x2EEA09, 0xD1921C, 0xFE1DEB, 0x1CB129,
  0xA73EE8, 0x8235F5, 0x2EBB44, 0x84E99C, 0x7026B4, 0x5F7E41,
  0x3991D6, 0x398353, 0x39F49C, 0x845F8B, 0xBDF928, 0x3B1FF8,
  0x97FFDE, 0x05980F, 0xEF2F11, 0x8B5A0A, 0x6D1F6D, 0x367ECF,
  0x27CB09, 0xB74F46, 0x3F669E, 0x5FEA2D, 0x7527BA, 0xC7EBE5,
  0xF17B3D, 0x0739F7, 0x8A5292, 0xEA6BFB, 0x5FB11F, 0x8D5D08,
  0x560330, 0x46FC7B, 0x6BABF0, 0xCFBC20, 0x9AF436, 0x1DA9E3,
  0x91615E, 0xE61B08, 0x659985, 0x5F14A0, 0x68408D, 0xFFD880,
  0x4D7327, 0x310606, 0x1556CA, 0x73A8C9, 0x60E27B, 0xC08C6B
};

/* pi/2 cut into 24-bit chunks (fdlibm). */
static const double PIo2[] = {
  1.57079625129699707031e+00,
  7.54978941586159635335e-08,
  5.39030252995776476554e-15,
  3.28200341580791294123e-22,
  1.27065575308067607349e-29,
  1.22933308981111328932e-36,
  2.73370053816464559624e-44,
  2.16741683877804819444e-51
};

/* fdlibm __rem_pio2_large: the integer (mod 8) and fraction parts of
 * x*(2/pi), without doing the full multiplication.  The input x is
 * positive, given as nx 24-bit chunks with the scaled exponent e0, and
 * on return y[0]+y[1] is x - n*(pi/2).  The return value is n mod 8. */
static int
rem_pio2_large(const double *x, double *y, int e0, int nx, int prec)
{
  static const int init_jk[] = { 3, 4, 4, 6 };
  int32_t jz, jx, jv, jp, jk, carry, n, iq[20], i, j, k, m, q0, ih;
  double z, fw, f[20], fq[20], q[20];

  jk = init_jk[prec];
  jp = jk;

  jx = nx - 1;
  jv = (e0 - 3) / 24;  if(jv < 0) jv = 0;
  q0 = e0 - 24*(jv + 1);

  j = jv - jx;  m = jx + jk;
  for(i = 0; i <= m; i++, j++)
    f[i] = j < 0 ? 0.0 : (double)two_over_pi[j];

  for(i = 0; i <= jk; i++) {
    for(j = 0, fw = 0.0; j <= jx; j++)
      fw += x[j]*f[jx + i - j];
    q[i] = fw;
  }

  jz = jk;
recompute:
  for(i = 0, j = jz, z = q[jz]; j > 0; i++, j--) {
    fw    = (double)(int32_t)(z*(1.0/16777216.0));
    iq[i] = (int32_t)(z - 16777216.0*fw);
    z     = q[j-1] + fw;
  }

  z  = ldexp(z, (int)q0);
  z -= 8.0*floor(z*0.125);
  n  = (int32_t)z;
  z -= (double)n;
  ih = 0;
  if(q0 > 0) {
    i  = iq[jz-1] >> (24 - q0);
    n += i;
    iq[jz-1] -= i << (24 - q0);
    ih = iq[jz-1] >> (23 - q0);
  }
  else if(q0 == 0) {
    ih = iq[jz-1] >> 23;
  }
  else if(z >= 0.5) {
    ih = 2;
  }

  if(ih > 0) {
    n += 1;  carry = 0;
    for(i = 0; i < jz; i++) {
      j = iq[i];
      if(carry == 0) {
        if(j != 0) {
          carry = 1;
          iq[i] = 0x1000000 - j;
        }
      }
      else {
        iq[i] = 0xffffff - j;
      }
    }
    if(q0 > 0) {
      switch(q0) {
      case 1:
        iq[jz-1] &= 0x7fffff;  break;
      case 2:
        iq[jz-1] &= 0x3fffff;  break;
      }
    }
    if(ih == 2) {
      z = 1.0 - z;
      if(carry != 0)
        z -= ldexp(1.0, (int)q0);
    }
  }

  if(z == 0.0) {
    j = 0;
    for(i = jz-1; i >= jk; i--) j |= iq[i];
    if(j == 0) {
      for(k = 1; iq[jk-k] == 0; k++);
      for(i = jz+1; i <= jz+k; i++) {
        f[jx+i] = (double)two_over_pi[jv+i];
        for(j = 0, fw = 0.0; j <= jx; j++)
          fw += x[j]*f[jx+i-j];
        q[i] = fw;
      }
      jz += k;
      goto recompute;
    }
  }

  if(z == 0.0) {
    jz -= 1;
    q0 -= 24;
    while(iq[jz] == 0) { jz--; q0 -= 24; }
  }
  else {
    z = ldexp(z, (int)-q0);
    if(z >= 16777216.0) {
      fw = (double)(int32_t)(z*(1.0/16777216.0));
      iq[jz] = (int32_t)(z - 16777216.0*fw);
      jz += 1;
      q0 += 24;
      iq[jz] = (int32_t)fw;
    }
    else {
      iq[jz] = (int32_t)z;
    }
  }

  fw = ldexp(1.0, (int)q0);
  for(i = jz; i >= 0; i--) {
    q[i] = fw*(double)iq[i];
    fw *= 1.0/16777216.0;
  }

  for(i = jz; i >= 0; i--) {
    for(fw = 0.0, k = 0; k <= jp && k <= jz-i; k++)
      fw += PIo2[k]*q[i+k];
    fq[jz-i] = fw;
  }

  switch(prec) {
  case 0:
    fw = 0.0;
    for(i = jz; i >= 0; i--) fw += fq[i];
    y[0] = ih == 0 ? fw : -fw;
    break;
  case 1:
  case 2:
    fw = 0.0;
    for(i = jz; i >= 0; i--) fw += fq[i];
    y[0] = ih == 0 ? fw : -fw;
    fw = fq[0] - fw;
    for(i = 1; i <= jz; i++) fw += fq[i];
    y[1] = ih == 0 ? fw : -fw;
    break;
  case 3:
    for(i = jz; i > 0; i--) {
      fw     = fq[i-1] + fq[i];
      fq[i] += fq[i-1] - fw;
      fq[i-1] = fw;
    }
    for(i = jz; i > 1; i--) {
      fw     = fq[i-1] + fq[i];
      fq[i] += fq[i-1] - fw;
      fq[i-1] = fw;
    }
    for(fw = 0.0, i = jz; i >= 2; i--) fw += fq[i];
    if(ih == 0) {
      y[0] =  fq[0];  y[1] =  fq[1];  y[2] =  fw;
    }
    else {
      y[0] = -fq[0];  y[1] = -fq[1];  y[2] = -fw;
    }
  }
  return n & 7;
}


/* Reduce a finite |x| > REM_PIO2_THRESHOLD modulo pi/2.  Returns n with
 * |x| = n*(pi/2) + r, |r| <= pi/4, the remainder r as the double-double
 * y0 + y1. */
static int
rem_pio2(double x, int *n, double *y0, double *y1)
{
  double tx[3], ty[2];
  double ax = fabs(x);
  double m, z;
  int ex, i;

  m = frexp(ax, &ex);           /* ax = m * 2^ex, 0.5 <= m < 1 */
  z = ldexp(m, 24);             /* 2^23 <= z < 2^24 */
  for(i = 0; i < 2; i++) {
    tx[i] = (double)(int32_t)z;
    z     = (z - tx[i])*16777216.0;
  }
  tx[2] = z;
  i = 2;
  while(tx[i] == 0.0) i--;

  *n = rem_pio2_large(tx, ty, ex - 24, i+1, 1);
  if(x < 0.0) {
    *y0 = -ty[0];
    *y1 = -ty[1];
    *n  = -*n;
  }
  else {
    *y0 = ty[0];
    *y1 = ty[1];
  }
  return GSL_SUCCESS;
}


/* sin(z) and cos(z) for |z| <= pi/4, from the Chebyshev series */
static int
sin_cos_pio4(double z, double *s, double *c)
{
  gsl_sf_result sin_result;
  gsl_sf_result cos_result;
  const double t  = 8.0*fabs(z)/M_PI - 1.0;
  const double z2 = z*z;
  const int stat_s = cheb_eval_e(&sin_cs, t, &sin_result);
  const int stat_c = cheb_eval_e(&cos_cs, t, &cos_result);

  *s = z * (1.0 + z2*sin_result.val);
  *c = 1.0 - 0.5*z2*(1.0 - z2*cos_result.val);
  return GSL_ERROR_SELECT_2(stat_s, stat_c);
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

/* I would have prefered just using the library sin() function.
 * But after some experimentation I decided that there was
 * no good way to understand the error; library sin() is just a black box.
 * So we have to roll our own.
 */
int
gsl_sf_sin_e(double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    const double P1 = 7.85398125648498535156e-1;
    const double P2 = 3.77489470793079817668e-8;
    const double P3 = 2.69515142907905952645e-15;

    const double sgn_x = GSL_SIGN(x);
    const double abs_x = fabs(x);

    if(abs_x < GSL_ROOT4_DBL_EPSILON) {
      const double x2 = x*x;
      result->val = x * (1.0 - x2/6.0);
      result->err = fabs(x*x2*x2 / 100.0);
      return GSL_SUCCESS;
    }
    else if(abs_x > REM_PIO2_THRESHOLD && gsl_finite(x)) {
      int n;
      double y0, y1;
      double s0, c0, sz, cz;
      const int stat_red = rem_pio2(x, &n, &y0, &y1);
      const int stat_cs  = sin_cos_pio4(y0, &s0, &c0);

      /* fold in the low word of the reduced argument to first order */
      sz = s0 + y1*c0;
      cz = c0 - y1*s0;

      switch(n & 3) {
      case 0:  result->val =  sz;  break;
      case 1:  result->val =  cz;  break;
      case 2:  result->val = -sz;  break;
      default: result->val = -cz;  break;
      }

      result->err = 2.0 * GSL_DBL_EPSILON * (fabs(result->val) + fabs(y1));
      return GSL_ERROR_SELECT_2(stat_red, stat_cs);
    }
    else {
      double sgn_result = sgn_x;
      double y = floor(abs_x/(0.25*M_PI));
      int octant = y - ldexp(floor(ldexp(y,-3)),3);
      int stat_cs;
      double z;

      if(GSL_IS_ODD(octant)) {
        octant += 1;
        octant &= 07;
        y += 1.0;
      }

      if(octant > 3) {
        octant -= 4;
        sgn_result = -sgn_result;
      }
      
      z = ((abs_x - y * P1) - y * P2) - y * P3;

      if(octant == 0) {
        gsl_sf_result sin_cs_result;
        const double t = 8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&sin_cs, t, &sin_cs_result);
        result->val = z * (1.0 + z*z * sin_cs_result.val);
      }
      else { /* octant == 2 */
        gsl_sf_result cos_cs_result;
        const double t = 8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&cos_cs, t, &cos_cs_result);
        result->val = 1.0 - 0.5*z*z * (1.0 - z*z * cos_cs_result.val);
      }

      result->val *= sgn_result;

      if(abs_x > 1.0/GSL_DBL_EPSILON) {
        result->err = fabs(result->val);
      }
      else if(abs_x > 100.0/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * abs_x * GSL_DBL_EPSILON * fabs(result->val);
      }
      else if(abs_x > 0.1/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * GSL_SQRT_DBL_EPSILON * fabs(result->val);
      }
      else {
        result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      }

      return stat_cs;
    }
  }
}


int
gsl_sf_cos_e(double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    const double P1 = 7.85398125648498535156e-1;
    const double P2 = 3.77489470793079817668e-8;
    const double P3 = 2.69515142907905952645e-15;

    const double abs_x = fabs(x);

    if(abs_x < GSL_ROOT4_DBL_EPSILON) {
      const double x2 = x*x;
      result->val = 1.0 - 0.5*x2;
      result->err = fabs(x2*x2/12.0);
      return GSL_SUCCESS;
    }
    else if(abs_x > REM_PIO2_THRESHOLD && gsl_finite(x)) {
      int n;
      double y0, y1;
      double s0, c0, sz, cz;
      const int stat_red = rem_pio2(x, &n, &y0, &y1);
      const int stat_cs  = sin_cos_pio4(y0, &s0, &c0);

      /* fold in the low word of the reduced argument to first order */
      sz = s0 + y1*c0;
      cz = c0 - y1*s0;

      switch(n & 3) {
      case 0:  result->val =  cz;  break;
      case 1:  result->val = -sz;  break;
      case 2:  result->val = -cz;  break;
      default: result->val =  sz;  break;
      }

      result->err = 2.0 * GSL_DBL_EPSILON * (fabs(result->val) + fabs(y1));
      return GSL_ERROR_SELECT_2(stat_red, stat_cs);
    }
    else {
      double sgn_result = 1.0;
      double y = floor(abs_x/(0.25*M_PI));
      int octant = y - ldexp(floor(ldexp(y,-3)),3);
      int stat_cs;
      double z;

      if(GSL_IS_ODD(octant)) {
        octant += 1;
        octant &= 07;
        y += 1.0;
      }

      if(octant > 3) {
        octant -= 4;
        sgn_result = -sgn_result;
      }

      if(octant > 1) {
        sgn_result = -sgn_result;
      }

      z = ((abs_x - y * P1) - y * P2) - y * P3;

      if(octant == 0) {
        gsl_sf_result cos_cs_result;
        const double t = 8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&cos_cs, t, &cos_cs_result);
        result->val = 1.0 - 0.5*z*z * (1.0 - z*z * cos_cs_result.val);
      }
      else { /* octant == 2 */
        gsl_sf_result sin_cs_result;
        const double t = 8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&sin_cs, t, &sin_cs_result);
        result->val = z * (1.0 + z*z * sin_cs_result.val);
      }

      result->val *= sgn_result;

      if(abs_x > 1.0/GSL_DBL_EPSILON) {
        result->err = fabs(result->val);
      }
      else if(abs_x > 100.0/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * abs_x * GSL_DBL_EPSILON * fabs(result->val);
      }
      else if(abs_x > 0.1/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * GSL_SQRT_DBL_EPSILON * fabs(result->val);
      }
      else {
        result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      }

      return stat_cs;
    }
  }
}


int
gsl_sf_hypot_e(const double x, const double y, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  /* As for gsl_hypot, an infinite argument gives +Inf even if the other
     argument is NaN; a NaN argument with no infinite argument propagates.
     The magnitude branch below would otherwise divide by, or evaluate
     sqrt of, the NaN and report an overflow instead. */
  if(gsl_isinf(x) || gsl_isinf(y)) {
    result->val = GSL_POSINF;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(gsl_isnan(x) || gsl_isnan(y)) {
    result->val = GSL_NAN;
    result->err = GSL_NAN;
    return GSL_SUCCESS;
  }
  else if(x == 0.0 && y == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    const double a = fabs(x);
    const double b = fabs(y);
    const double min = GSL_MIN_DBL(a,b);
    const double max = GSL_MAX_DBL(a,b);
    const double rat = min/max;
    const double root_term = sqrt(1.0 + rat*rat);

    if(max < GSL_DBL_MAX/root_term) {
      result->val = max * root_term;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      OVERFLOW_ERROR(result);
    }
  }
}


int
gsl_sf_complex_sin_e(const double zr, const double zi,
                        gsl_sf_result * szr, gsl_sf_result * szi)
{
  /* CHECK_POINTER(szr) */
  /* CHECK_POINTER(szi) */

  if(fabs(zi) < 1.0) {
    double ch_m1, sh;
    sinh_series(zi, &sh);
    cosh_m1_series(zi, &ch_m1);
    szr->val = sin(zr)*(ch_m1 + 1.0);
    szi->val = cos(zr)*sh;
    szr->err = 2.0 * GSL_DBL_EPSILON * fabs(szr->val);
    szi->err = 2.0 * GSL_DBL_EPSILON * fabs(szi->val);
    return GSL_SUCCESS;
  }
  else if(fabs(zi) < GSL_LOG_DBL_MAX) {
    double ex = exp(zi);
    double ch = 0.5*(ex+1.0/ex);
    double sh = 0.5*(ex-1.0/ex);
    szr->val = sin(zr)*ch;
    szi->val = cos(zr)*sh;
    szr->err = 2.0 * GSL_DBL_EPSILON * fabs(szr->val);
    szi->err = 2.0 * GSL_DBL_EPSILON * fabs(szi->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR_2(szr, szi);
  }
}


int
gsl_sf_complex_cos_e(const double zr, const double zi,
                        gsl_sf_result * czr, gsl_sf_result * czi)
{
  /* CHECK_POINTER(czr) */
  /* CHECK_POINTER(czi) */

  if(fabs(zi) < 1.0) {
    double ch_m1, sh;
    sinh_series(zi, &sh);
    cosh_m1_series(zi, &ch_m1);
    czr->val =  cos(zr)*(ch_m1 + 1.0);
    czi->val = -sin(zr)*sh;
    czr->err = 2.0 * GSL_DBL_EPSILON * fabs(czr->val);
    czi->err = 2.0 * GSL_DBL_EPSILON * fabs(czi->val);
    return GSL_SUCCESS;
  }
  else if(fabs(zi) < GSL_LOG_DBL_MAX) {
    double ex = exp(zi);
    double ch = 0.5*(ex+1.0/ex);
    double sh = 0.5*(ex-1.0/ex);
    czr->val =  cos(zr)*ch;
    czi->val = -sin(zr)*sh;
    czr->err = 2.0 * GSL_DBL_EPSILON * fabs(czr->val);
    czi->err = 2.0 * GSL_DBL_EPSILON * fabs(czi->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR_2(czr,czi);
  }
}


int
gsl_sf_complex_logsin_e(const double zr, const double zi,
                           gsl_sf_result * lszr, gsl_sf_result * lszi)
{
  /* CHECK_POINTER(lszr) */
  /* CHECK_POINTER(lszi) */

  if(zi > 60.0) {
    lszr->val = -M_LN2 + zi;
    lszi->val =  0.5*M_PI - zr;
    lszr->err = 2.0 * GSL_DBL_EPSILON * fabs(lszr->val);
    lszi->err = 2.0 * GSL_DBL_EPSILON * fabs(lszi->val);
  }
  else if(zi < -60.0) {
    lszr->val = -M_LN2 - zi;
    lszi->val = -0.5*M_PI + zr;
    lszr->err = 2.0 * GSL_DBL_EPSILON * fabs(lszr->val);
    lszi->err = 2.0 * GSL_DBL_EPSILON * fabs(lszi->val);
  }
  else {
    gsl_sf_result sin_r, sin_i;
    int status;
    gsl_sf_complex_sin_e(zr, zi, &sin_r, &sin_i); /* ok by construction */
    status = gsl_sf_complex_log_e(sin_r.val, sin_i.val, lszr, lszi);
    if(status == GSL_EDOM) {
      DOMAIN_ERROR_2(lszr, lszi);
    }
  }
  return gsl_sf_angle_restrict_symm_e(&(lszi->val));
}


int
gsl_sf_lnsinh_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(fabs(x) < 1.0) {
    double eps;
    sinh_series(x, &eps);
    result->val = log(eps);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < -0.5*GSL_LOG_DBL_EPSILON) {
    result->val = x + log(0.5*(1.0 - exp(-2.0*x)));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = -M_LN2 + x;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int gsl_sf_lncosh_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(fabs(x) < 1.0) {
    double eps;
    cosh_m1_series(x, &eps);
    return gsl_sf_log_1plusx_e(eps, result);
  }
  else if(fabs(x) < -0.5*GSL_LOG_DBL_EPSILON) {
    result->val = fabs(x) + log(0.5*(1.0 + exp(-2.0*fabs(x))));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = -M_LN2 + fabs(x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


/*
inline int gsl_sf_sincos_e(const double theta, double * s, double * c)
{
  double tan_half = tan(0.5 * theta);
  double den = 1. + tan_half*tan_half;
  double cos_theta = (1.0 - tan_half*tan_half) / den;
  double sin_theta = 2.0 * tan_half / den;
}
*/

int
gsl_sf_polar_to_rect(const double r, const double theta,
                          gsl_sf_result * x, gsl_sf_result * y)
{
  double t   = theta;
  int status = gsl_sf_angle_restrict_symm_e(&t);
  double c = cos(t);
  double s = sin(t);
  x->val = r * cos(t);
  y->val = r * sin(t);
  x->err  = r * fabs(s * GSL_DBL_EPSILON * t);
  x->err += 2.0 * GSL_DBL_EPSILON * fabs(x->val);
  y->err  = r * fabs(c * GSL_DBL_EPSILON * t);
  y->err += 2.0 * GSL_DBL_EPSILON * fabs(y->val);
  return status;
}


int
gsl_sf_rect_to_polar(const double x, const double y,
                          gsl_sf_result * r, gsl_sf_result * theta)
{
  int stat_h = gsl_sf_hypot_e(x, y, r);
  if(r->val > 0.0) {
    theta->val = atan2(y, x);
    theta->err = 2.0 * GSL_DBL_EPSILON * fabs(theta->val);
    return stat_h;
  }
  else {
    DOMAIN_ERROR(theta);
  }
}


/* Reduce theta modulo 2 pi into [0, 2 pi) (`pos') or (-pi, pi] (`symm').
 *
 * Above |theta| = 2^21 the three-term representation of 2 pi used by the
 * small-argument path loses the reduced angle -- the same defect that was
 * fixed in gsl_sf_sin_e/gsl_sf_cos_e (Savannah bug #45746, see rem_pio2
 * above).  The large branch reuses the exact pi/2 reduction and maps the
 * n mod 8 octant plus the remainder onto the wanted interval. */
static int
angle_restrict_pio2(double theta, int pos, double *val, double *err)
{
  int n;
  double y0, y1;
  double r, s, a;
  const double ax = fabs(theta);

  rem_pio2(ax, &n, &y0, &y1);

  /* r = |theta| - n*(pi/2) in [-pi/4, pi/4] */
  r = y0 + y1;

  /* Build the symmetric value s in (-pi, pi] from the positive argument,
   * then restore the sign of theta.  With |theta| = n*(pi/2) + r:
   *   n mod 4 = 0: |theta| = r
   *             1: |theta| = pi/2 + r
   *             2: |theta| = pi + r
   *             3: |theta| = 3*pi/2 + r
   */
  switch(n & 3) {
  case 0:  s = r;              break;
  case 1:  s = M_PI/2 + r;     break;
  case 2:  s = M_PI + r;       break;
  default: s = 3.0*M_PI/2 + r; break;
  }
  if(s > M_PI) s -= 2.0*M_PI;

  if(theta < 0.0) s = -s;

  if(pos) {
    *val = (s < 0.0) ? s + 2.0*M_PI : s;
  }
  else {
    *val = s;
  }

  a = fabs(*val);
  *err = 2.0 * GSL_DBL_EPSILON * (a > M_PI ? M_PI : a) + GSL_DBL_EPSILON;
  return GSL_SUCCESS;
}


int gsl_sf_angle_restrict_symm_err_e(const double theta, gsl_sf_result * result)
{
  if(fabs(theta) > REM_PIO2_THRESHOLD && gsl_finite(theta)) {
    int stat = angle_restrict_pio2(theta, 0, &result->val, &result->err);
    if(fabs(theta) > 0.0625/GSL_DBL_EPSILON) {
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      GSL_ERROR ("error", GSL_ELOSS);
    }
    return stat;
  }
  else {
    /* synthetic extended precision constants */
    const double P1 = 4 * 7.8539812564849853515625e-01;
    const double P2 = 4 * 3.7748947079307981766760e-08;
    const double P3 = 4 * 2.6951514290790594840552e-15;
    const double TwoPi = 2*(P1 + P2 + P3);

    const double y = GSL_SIGN(theta) * 2 * floor(fabs(theta)/TwoPi);
    double r = ((theta - y*P1) - y*P2) - y*P3;

    if(r >  M_PI) { r = (((r-2*P1)-2*P2)-2*P3); }  /* r-TwoPi */
    else if (r < -M_PI) r = (((r+2*P1)+2*P2)+2*P3); /* r+TwoPi */

    result->val = r;

    if(fabs(theta) > 0.0625/GSL_DBL_EPSILON) {
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      GSL_ERROR ("error", GSL_ELOSS);
    }
    else if(fabs(theta) > 0.0625/GSL_SQRT_DBL_EPSILON) {
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val - theta);
      return GSL_SUCCESS;
    }
    else {
      double delta = fabs(result->val - theta);
      result->err = 2.0 * GSL_DBL_EPSILON * ((delta < M_PI) ? delta : M_PI);
      return GSL_SUCCESS;
    }
  }
}


int gsl_sf_angle_restrict_pos_err_e(const double theta, gsl_sf_result * result)
{
  if(fabs(theta) > REM_PIO2_THRESHOLD && gsl_finite(theta)) {
    int stat = angle_restrict_pio2(theta, 1, &result->val, &result->err);
    if(fabs(theta) > 0.0625/GSL_DBL_EPSILON) {
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      GSL_ERROR ("error", GSL_ELOSS);
    }
    return stat;
  }
  else {
    /* synthetic extended precision constants */
    const double P1 = 4 * 7.85398125648498535156e-01;
    const double P2 = 4 * 3.77489470793079817668e-08;
    const double P3 = 4 * 2.69515142907905952645e-15;
    const double TwoPi = 2*(P1 + P2 + P3);

    const double y = 2*floor(theta/TwoPi);

    double r = ((theta - y*P1) - y*P2) - y*P3;

    if(r > TwoPi) {r = (((r-2*P1)-2*P2)-2*P3); }  /* r-TwoPi */
    else if (r < 0) { /* may happen due to FP rounding */
      r = (((r+2*P1)+2*P2)+2*P3); /* r+TwoPi */
    }

    result->val = r;

    if(fabs(theta) > 0.0625/GSL_DBL_EPSILON) {
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      GSL_ERROR ("error", GSL_ELOSS);
    }
    else if(fabs(theta) > 0.0625/GSL_SQRT_DBL_EPSILON) {
      result->err = GSL_DBL_EPSILON * fabs(result->val - theta);
      return GSL_SUCCESS;
    }
    else {
      double delta = fabs(result->val - theta);
      result->err = 2.0 * GSL_DBL_EPSILON * ((delta < M_PI) ? delta : M_PI);
      return GSL_SUCCESS;
    }
  }
}


int gsl_sf_angle_restrict_symm_e(double * theta)
{
  gsl_sf_result r;
  int stat = gsl_sf_angle_restrict_symm_err_e(*theta, &r);
  *theta = r.val;
  return stat;
}


int gsl_sf_angle_restrict_pos_e(double * theta)
{
  gsl_sf_result r;
  int stat = gsl_sf_angle_restrict_pos_err_e(*theta, &r);
  *theta = r.val;
  return stat;
}


int gsl_sf_sin_err_e(const double x, const double dx, gsl_sf_result * result)
{
  int stat_s = gsl_sf_sin_e(x, result);
  result->err += fabs(cos(x) * dx);
  result->err += GSL_DBL_EPSILON * fabs(result->val);
  return stat_s;
}


int gsl_sf_cos_err_e(const double x, const double dx, gsl_sf_result * result)
{
  int stat_c = gsl_sf_cos_e(x, result);
  result->err += fabs(sin(x) * dx);
  result->err += GSL_DBL_EPSILON * fabs(result->val);
  return stat_c;
}


#if 0
int
gsl_sf_sin_pi_x_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(-100.0 < x && x < 100.0) {
    result->val = sin(M_PI * x) / (M_PI * x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double N = floor(x + 0.5);
    const double f = x - N;

    if(N < INT_MAX && N > INT_MIN) {
      /* Make it an integer if we can. Saves another
       * call to floor().
       */
      const int intN    = (int)N;
      const double sign = ( GSL_IS_ODD(intN) ? -1.0 : 1.0 );
      result->val = sign * sin(M_PI * f);
      result->err = GSL_DBL_EPSILON * fabs(result->val);
    }
    else if(N > 2.0/GSL_DBL_EPSILON || N < -2.0/GSL_DBL_EPSILON) {
      /* All integer-valued floating point numbers
       * bigger than 2/eps=2^53 are actually even.
       */
      result->val = 0.0;
      result->err = 0.0;
    }
    else {
      const double resN = N - 2.0*floor(0.5*N); /* 0 for even N, 1 for odd N */
      const double sign = ( fabs(resN) > 0.5 ? -1.0 : 1.0 );
      result->val = sign * sin(M_PI*f);
      result->err = GSL_DBL_EPSILON * fabs(result->val);
    }

    return GSL_SUCCESS;
  }
}
#endif


int gsl_sf_sinc_e(double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    const double ax = fabs(x);

    if(ax < 0.8) {
      /* Do not go to the limit of the fit since
       * there is a zero there and the Chebyshev
       * accuracy will go to zero.
       */
      return cheb_eval_e(&sinc_cs, 2.0*ax-1.0, result);
    }
    else if(ax < 100.0) {
      /* Small arguments are no problem.
       * We trust the library sin() to
       * roughly machine precision.
       */
      result->val = sin(M_PI * ax)/(M_PI * ax);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      /* Large arguments must be handled separately.
       */
      const double r = M_PI*ax;
      gsl_sf_result s;
      int stat_s = gsl_sf_sin_e(r, &s);
      result->val = s.val/r;
      result->err = s.err/r + 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return stat_s;
    }
  }
}



/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double gsl_sf_sin(const double x)
{
  EVAL_RESULT(gsl_sf_sin_e(x, &result));
}

double gsl_sf_cos(const double x)
{
  EVAL_RESULT(gsl_sf_cos_e(x, &result));
}

double gsl_sf_hypot(const double x, const double y)
{
  EVAL_RESULT(gsl_sf_hypot_e(x, y, &result));
}

double gsl_sf_lnsinh(const double x)
{
  EVAL_RESULT(gsl_sf_lnsinh_e(x, &result));
}

double gsl_sf_lncosh(const double x)
{
  EVAL_RESULT(gsl_sf_lncosh_e(x, &result));
}

double gsl_sf_angle_restrict_symm(const double theta)
{
  double result = theta;
  EVAL_DOUBLE(gsl_sf_angle_restrict_symm_e(&result));
}

double gsl_sf_angle_restrict_pos(const double theta)
{
  double result = theta;
  EVAL_DOUBLE(gsl_sf_angle_restrict_pos_e(&result));
}

#if 0
double gsl_sf_sin_pi_x(const double x)
{
  EVAL_RESULT(gsl_sf_sin_pi_x_e(x, &result));
}
#endif

double gsl_sf_sinc(const double x)
{
  EVAL_RESULT(gsl_sf_sinc_e(x, &result));
}
