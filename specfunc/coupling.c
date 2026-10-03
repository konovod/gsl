/* specfunc/coupling.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002 Gerard Jungman
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

/* Author:  G. Jungman 
 *          A. A. Illarionov */

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>

#include "error.h"

inline
static
int locMax3(const int a, const int b, const int c)
{
  int d = GSL_MAX(a, b);
  return GSL_MAX(d, c);
}

inline
static
int locMin3(const int a, const int b, const int c)
{
  int d = GSL_MIN(a, b);
  return GSL_MIN(d, c);
}

/* Calculation of the 6j symbol by the Schulten-Gordon three-term
 * recurrence in the first angular momentum, with a forward and a
 * backward recursion matched in the interior.  This avoids both the
 * factorial overflow of the Racah sum and the catastrophic cancellation
 * it suffers for large angular momenta.
 *
 * Adapted from the public-domain SLATEC routine DRC6J, which implements
 * Schulten & Gordon, J. Math. Phys. 16 (1975) 1961 and
 * Comput. Phys. Commun. 11 (1976) 269.  The recurrence variable is
 * two_ja; the other five arguments are held fixed.
 */

static
int
coupling_6j_recurrence(int two_ja, int two_jb, int two_jc,
                       int two_jd, int two_je, int two_jf,
                       gsl_sf_result * result)
{
  const double L2 = 0.5 * two_jb;
  const double L3 = 0.5 * two_jc;
  const double L4 = 0.5 * two_jd;
  const double L5 = 0.5 * two_je;
  const double L6 = 0.5 * two_jf;
  const double L1 = 0.5 * two_ja;
  const double huge = sqrt(GSL_DBL_MAX / 20.0);
  const double srhuge = sqrt(huge);
  const double tiny = 1.0 / huge;
  const double srtiny = 1.0 / srhuge;
  const double l1min = GSL_MAX(fabs(L2 - L3), fabs(L5 - L6));
  const double l1max = GSL_MIN(L2 + L3, L5 + L6);
  const int nfin = (int)(l1max - l1min + 1.01);
  double stack[64];
  double * six;
  double l1, newfac, oldfac, c1, c1old, c2 = 0.0, dv, denom, x = 0.0, y = 0.0;
  double x1 = 0.0, x2 = 0.0, x3 = 0.0, y1 = 0.0, y2 = 0.0, y3 = 0.0;
  double sum1 = 0.0, sum2 = 0.0, sumfor = 0.0, sumbac = 0.0, sumuni = 0.0;
  double ratio, cnorm, sign1, sign2;
  int lstep, nstep2, nfinp2, nfinp3, nlim, n, have_sumuni;
  const int lsum = (two_jb + two_jc + two_je + two_jf) / 2;

  if(nfin <= 0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }

  if(nfin <= 1) {
    /* The triangle conditions leave a single value for two_ja. */
    result->val = (GSL_IS_ODD(lsum) ? -1.0 : 1.0)
                  / sqrt((l1min + l1min + 1.0) * (L4 + L4 + 1.0));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }

  if(nfin <= (int)(sizeof(stack) / sizeof(stack[0]))) {
    six = stack;
  }
  else {
    six = (double *) malloc(nfin * sizeof(double));
    if(six == NULL) {
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      GSL_ERROR("allocation failed", GSL_ENOMEM);
    }
  }

  /* Forward recursion from l1min. */
  l1 = l1min;
  newfac = 0.0;
  c1 = 0.0;
  c1old = 0.0;
  denom = 0.0;
  x = 0.0;
  have_sumuni = 0;
  six[0] = srtiny;
  sum1 = (l1 + l1 + 1.0) * tiny;
  lstep = 1;

  for(;;) {
    lstep++;
    l1 += 1.0;
    oldfac = newfac;
    {
      double a1 = (l1 + L2 + L3 + 1.0) * (l1 - L2 + L3)
                * (l1 + L2 - L3) * (-l1 + L2 + L3 + 1.0);
      double a2 = (l1 + L5 + L6 + 1.0) * (l1 - L5 + L6)
                * (l1 + L5 - L6) * (-l1 + L5 + L6 + 1.0);
      newfac = sqrt(a1 * a2);
    }

    if(l1 < 1.01) {
      c1 = -2.0 * (L2 * (L2 + 1.0) + L5 * (L5 + 1.0) - L4 * (L4 + 1.0))
           / newfac;
    }
    else {
      dv = 2.0 * (L2 * (L2 + 1.0) * L5 * (L5 + 1.0)
                  + L3 * (L3 + 1.0) * L6 * (L6 + 1.0)
                  - l1 * (l1 - 1.0) * L4 * (L4 + 1.0))
         - (L2 * (L2 + 1.0) + L3 * (L3 + 1.0) - l1 * (l1 - 1.0))
           * (L5 * (L5 + 1.0) + L6 * (L6 + 1.0) - l1 * (l1 - 1.0));
      denom = (l1 - 1.0) * newfac;
      if(lstep > 2)
        c1old = fabs(c1);
      c1 = -(l1 + l1 - 1.0) * dv / denom;
    }

    if(lstep <= 2) {
      x = srtiny * c1;
      six[1] = x;
      sum1 += tiny * (l1 + l1 + 1.0) * c1 * c1;
      if(lstep == nfin) {
        sumuni = sum1;
        have_sumuni = 1;
        break;
      }
      continue;
    }

    c2 = -l1 * oldfac / denom;
    x = c1 * six[lstep - 2] + c2 * six[lstep - 3];
    six[lstep - 1] = x;
    sumfor = sum1;
    sum1 += (l1 + l1 + 1.0) * x * x;
    if(lstep == nfin)
      break;

    if(fabs(x) >= srhuge) {
      for(n = 0; n < lstep; n++) {
        if(fabs(six[n]) < srtiny)
          six[n] = 0.0;
        six[n] /= srhuge;
      }
      sum1 /= huge;
      sumfor /= huge;
      x /= srhuge;
    }

    /* Continue only while |c1| decreases; the recurrence is stable in
     * the downhill direction. */
    if(c1old - fabs(c1) <= 0.0)
      break;
  }

  x1 = x;
  x2 = six[lstep - 2];
  x3 = six[lstep - 3];

  if(!have_sumuni) {
    /* Backward recursion from l1max, matched to the forward one. */
    nfinp2 = nfin + 2;
    nfinp3 = nfin + 3;
    nstep2 = nfin - lstep + 3;
    l1 = l1max;
    six[nfin - 1] = srtiny;
    sum2 = (l1 + l1 + 1.0) * tiny;
    l1 += 2.0;
    lstep = 1;
    y = 0.0;

    for(;;) {
      lstep++;
      l1 -= 1.0;
      oldfac = newfac;
      {
        double a1s = (l1 + L2 + L3) * (l1 - L2 + L3 - 1.0)
                   * (l1 + L2 - L3 - 1.0) * (-l1 + L2 + L3 + 2.0);
        double a2s = (l1 + L5 + L6) * (l1 - L5 + L6 - 1.0)
                   * (l1 + L5 - L6 - 1.0) * (-l1 + L5 + L6 + 2.0);
        newfac = sqrt(a1s * a2s);
      }
      dv = 2.0 * (L2 * (L2 + 1.0) * L5 * (L5 + 1.0)
                  + L3 * (L3 + 1.0) * L6 * (L6 + 1.0)
                  - l1 * (l1 - 1.0) * L4 * (L4 + 1.0))
         - (L2 * (L2 + 1.0) + L3 * (L3 + 1.0) - l1 * (l1 - 1.0))
           * (L5 * (L5 + 1.0) + L6 * (L6 + 1.0) - l1 * (l1 - 1.0));
      denom = l1 * newfac;
      c1 = -(l1 + l1 - 1.0) * dv / denom;

      if(lstep <= 2) {
        y = srtiny * c1;
        six[nfin - 2] = y;
        if(lstep == nstep2)
          break;
        sumbac = sum2;
        sum2 += tiny * (l1 + l1 - 3.0) * c1 * c1;
        continue;
      }

      c2 = -(l1 - 1.0) * oldfac / denom;
      y = c1 * six[nfinp2 - lstep - 1] + c2 * six[nfinp3 - lstep - 1];
      if(lstep == nstep2)
        break;
      six[nfin - lstep] = y;
      sumbac = sum2;
      sum2 += (l1 + l1 - 3.0) * y * y;

      if(fabs(y) >= srhuge) {
        for(n = 0; n < lstep; n++) {
          int index = nfin - n;
          if(fabs(six[index - 1]) < srtiny)
            six[index - 1] = 0.0;
          six[index - 1] /= srhuge;
        }
        sumbac /= huge;
        sum2 /= huge;
      }
    }

    y3 = y;
    y2 = six[nfinp2 - lstep - 1];
    y1 = six[nfinp3 - lstep - 1];

    ratio = (x1 * y1 + x2 * y2 + x3 * y3)
            / (x1 * x1 + x2 * x2 + x3 * x3);
    nlim = nfin - nstep2 + 1;

    if(fabs(ratio) >= 1.0) {
      for(n = 0; n < nlim; n++)
        six[n] *= ratio;
      sumuni = ratio * ratio * sumfor + sumbac;
    }
    else {
      nlim += 1;
      ratio = 1.0 / ratio;
      for(n = nlim - 1; n < nfin; n++)
        six[n] *= ratio;
      sumuni = sumfor + ratio * ratio * sumbac;
    }
  }

  /* Normalise to (2 L4 + 1) sum (2 L1 + 1) {6j}^2 = 1 and fix the
   * overall sign by the phase of the last coefficient. */
  cnorm = 1.0 / sqrt((L4 + L4 + 1.0) * sumuni);
  sign1 = (six[nfin - 1] >= 0.0) ? 1.0 : -1.0;
  sign2 = GSL_IS_ODD(lsum) ? -1.0 : 1.0;
  if(sign1 * sign2 <= 0.0)
    cnorm = -cnorm;

  {
    int index = (int)(L1 - l1min + 0.5);
    const double val = cnorm * six[index];
    result->val = val;
    result->err = (8.0 + 2.0 * nfin) * GSL_DBL_EPSILON * fabs(val);
  }

  if(six != stack)
    free(six);

  return GSL_SUCCESS;
}


static
int
triangle_selection_fails(int two_ja, int two_jb, int two_jc)
{
  /*
   * enough to check the triangle condition for one spin vs. the other two
   */
  return ( (two_jb < abs(two_ja - two_jc)) || (two_jb > two_ja + two_jc) ||
           GSL_IS_ODD(two_ja + two_jb + two_jc) );
}


static
int
m_selection_fails(int two_ja, int two_jb, int two_jc,
                  int two_ma, int two_mb, int two_mc)
{
  return (
         abs(two_ma) > two_ja 
      || abs(two_mb) > two_jb
      || abs(two_mc) > two_jc
      || GSL_IS_ODD(two_ja + two_ma)
      || GSL_IS_ODD(two_jb + two_mb)
      || GSL_IS_ODD(two_jc + two_mc)
      || (two_ma + two_mb + two_mc) != 0
          );
}

/* Calculation of 3j-symbol is based on edge-recursion algorithm 
 * described by Tuzun, Comp Phys Comm 112 (1998) */

extern int regge3j_index_permutations_table[720];

static inline
int
regge3j_find_ordered_params(
    int twoj1, int twoj2, int twoj3, int twom1, int twom2, int twom3, 
    int regge[9], int *  sign_param)
{
    int tmp[9] = {
        (-twoj1 + twoj2 + twoj3)/2, (twoj1 - twoj2 + twoj3)/ 2, (twoj1 + twoj2 - twoj3)/2, 
        (twoj1 - twom1)/2, (twoj2 - twom2)/2, (twoj3 - twom3)/2, 
        (twoj1 + twom1)/2, (twoj2 + twom2)/2 , (twoj3 + twom3)/2
    };
    int j;
    int * i = regge3j_index_permutations_table, 
        * ie = regge3j_index_permutations_table + 720;

    for (; i < ie ; i += 10) {
        if (tmp[i[0]] > tmp[i[4]]) continue;
        if (tmp[i[4]] > tmp[i[8]]) continue;
        if (tmp[i[8]] > tmp[i[3]]) continue;
        if (tmp[i[3]] > tmp[i[1]]) continue;
        if (tmp[i[4]] > tmp[i[7]]) continue;
        if ((tmp[i[4]] == tmp[i[7]]) && (tmp[i[5]] > tmp[i[8]])) continue;
        
        for(j = 0; j < 9 ; ++j)
            regge[j] = tmp[i[j]];
        *sign_param = i[9];
        return 0;
    }

    return 1;
}

static inline
int
m1_pow(int n)
{
  return 1 - 2*(n&1);
}

static inline
void
coupling_3j_sum_e(int r[9], gsl_sf_result * result)
{
  gsl_sf_result pp, p;
  double tmp;
  int k;

  gsl_sf_choose_e(r[2], r[7], &pp);
  pp.val *= m1_pow(r[7]);

  if (r[0] == 0) {
    *result = pp;
    return ;
  }

  tmp = ((double)(r[7] * r[5])) / (r[2] - r[7] + 1) ; 
  p.val = pp.val * (r[8] - tmp);
  p.err = pp.err * fabs(r[8] - tmp) + fabs(pp.val) * 
        (abs(r[8]) + fabs(tmp)) * 2.0 * GSL_DBL_EPSILON;
  if (r[0] == 1) {
    *result = p;
    return;
  }

  for (k = r[0] - 2; k >= 0; --k) {
    tmp = (r[0] - k)*(r[4] - k);
    double tmp1 = (-r[2]*(r[5] + 1) + (r[4] - k - 1)*(r[1] + k + 2) + 
        (r[3] + k + 1)*(r[0] - k)) / tmp;
    double tmp2 = (r[3] + k + 2)*(r[1] + k + 2) / tmp;

    result->val = p.val * tmp1 - pp.val * tmp2;
    result->err = p.err * tmp1 + pp.err * tmp2 + 
        2.0 * GSL_DBL_EPSILON * (fabs(p.val * tmp1) + fabs(pp.val*tmp2));

    pp = p;
    p = *result;
  }
    
  return;
}

static inline
int
plus_log_pochhammer_ratio_e(gsl_sf_result * result, const int a, const int b, int n)
{
  double r = 1.;
  int status;
  gsl_sf_result f;

  int i;
  for(i = n ; i > 0 ; --i) 
    r *= ((double)(a + i - 1))/((double)(b + i - 1));   

  status = gsl_sf_log_e(r, &f);
  f.err += (n * 2) * GSL_DBL_EPSILON;

  result->val += f.val;
  result->err += f.err; // here we use the fact that all values have the same sign

  return status;
}


static inline
int
coupling_3j_log_factor_e(int r[9], gsl_sf_result * result)
{
  int status;
  
  status = gsl_sf_log_e(r[0] + r[1] + r[2] + 1., result);
  result->val *= -1;

  status += plus_log_pochhammer_ratio_e(result, r[8] + 1, r[6] + r[8] + 1, r[2] - r[4]);
  status += plus_log_pochhammer_ratio_e(result, r[5] + 1, r[7] + 1, r[6] - r[5]);
  status += plus_log_pochhammer_ratio_e(result, r[0] + 1, r[8] + r[5] + 1, r[6] - r[5]);
  status += plus_log_pochhammer_ratio_e(result, 1, r[1] + 1, r[0]);
  status += plus_log_pochhammer_ratio_e(result, 1, r[3] + r[6] + 1, r[0]);

  return status;
}

static inline
void
coupling_3j_0_0_0_abs_e(int a, int b, int c, gsl_sf_result * result)
/* 
 * a = (twoj1 + twoj3 - twoj2) / 2 and even
 * b = (twoj1 + twoj2 - twoj3) / 2 and even
 * c = (twoj2 + twoj3 - twoj1) / 2 and even
 * a >= b >= c
 */
{
  int a2 = a / 2, b2 = b / 2, c2 = c / 2;
  int i;
  double val = 1.;

  for (i = 1 ; i <= b2 ; ++i) {
    val *= (b2 + i) * (a2 + i) * (a2 + i);
    val /= (a  + i) * (a + b2 + i) * i;
  }
  for (i = 1 ; i <= c2 ; ++i) {
    val *= (c2 + i) * (a2 + b2 + i) * (a2 + b2 + i);
    val /= (a + b + i) * (a + b + c2 + i) * i;
  }

  result->val = sqrt(val / (a + b + c + 1.));
  result->err = (b2 + c2 + 1) * GSL_DBL_EPSILON * fabs(result->val);

  return;
}

static inline
void
sort_3(int * a0, int * a1 , int * a2)
{
#ifndef locswap
#define locswap(a, b) {int t = (b); (b) = (a); (a) = t;}
  if ( *a1 > *a2 ) locswap(*a1, *a2);
  if ( *a0 > *a1 ) locswap(*a0, *a1);
  if ( *a1 > *a2 ) locswap(*a1, *a2);
#undef locswap
#endif
}



/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_coupling_3j_e (int two_ja, int two_jb, int two_jc,
                      int two_ma, int two_mb, int two_mc,
                      gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(two_ja < 0 || two_jb < 0 || two_jc < 0) {
    DOMAIN_ERROR(result);
  }
  else if (   triangle_selection_fails(two_ja, two_jb, two_jc)
           || m_selection_fails(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc)
     ) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    int j_sum = (two_ja + two_jb + two_jc) / 2;
    if ( two_ma == 0 && two_mb == 0 && two_mc == 0 ) {
      if (j_sum % 2 != 0) { 
        result->val = 0.0;
        result->err = 0.0;
        return GSL_SUCCESS;
      } else {
        int a = (two_ja + two_jb - two_jc) / 2, 
            b = (two_jb + two_jc - two_ja) / 2, 
            c = (two_jc + two_ja - two_jb) / 2;
        sort_3(&c, &b, &a);
        
        coupling_3j_0_0_0_abs_e(a, b, c, result);
        result->val *= m1_pow(j_sum/2);
        return GSL_SUCCESS;
      }
    } else {
      int regge[9] = {0}, sign_param;
        
      regge3j_find_ordered_params(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc, 
            regge, &sign_param);

      if (regge[3] == regge[6] && regge[4] == regge[7] && regge[5] == regge[8]) {
        if (j_sum % 2 != 0) {
          result->val = 0.;
          result->err = 0.;
          return GSL_SUCCESS;
        } else {
          coupling_3j_0_0_0_abs_e(regge[1], regge[2], regge[0], result);
          result->val *= m1_pow(j_sum/2 + sign_param * 
              (regge[0] + regge[1] + regge[2]));
          return GSL_SUCCESS;
        }
      } else {
        int status;
        gsl_sf_result factor;

        /* we have to compute factor anyway for proper error estimation */
        status = coupling_3j_log_factor_e(regge, &factor);
        factor.val *= 0.5;
        factor.err *= 0.5;

        coupling_3j_sum_e(regge, result);

        if (result->val == 0.) {
          result->err *= exp(factor.val);
          return status;
        } else {
          double rel_err = result->err / fabs(result->val);

          sign_param = GSL_SIGN(result->val) * m1_pow((regge[1] - regge[8]) + 
              (regge[0] + regge[1] + regge[2]) * sign_param); 
          /* now contains the sign */

          status += gsl_sf_log_abs_e(result->val, result);
          result->err += rel_err + factor.err + 2.0 * GSL_DBL_EPSILON * 
              (fabs(result->val) + fabs(factor.val));
          result->val += factor.val;

          result->val = exp(result->val);
          result->err *= result->val;

          result->val *= sign_param;

          return status;
        }
      }
    }
  }
}

#ifndef GSL_DISABLE_DEPRECATED

int
gsl_sf_coupling_6j_INCORRECT_e(int two_ja, int two_jb, int two_jc,
                               int two_jd, int two_je, int two_jf,
                               gsl_sf_result * result)
{
  return gsl_sf_coupling_6j_e(two_ja, two_jb, two_je, two_jd, two_jc, two_jf, result);
}

#endif


int
gsl_sf_coupling_6j_e(int two_ja, int two_jb, int two_jc,
                     int two_jd, int two_je, int two_jf,
                     gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(   two_ja < 0 || two_jb < 0 || two_jc < 0
     || two_jd < 0 || two_je < 0 || two_jf < 0
     ) {
    DOMAIN_ERROR(result);
  }
  else if(   triangle_selection_fails(two_ja, two_jb, two_jc)
          || triangle_selection_fails(two_ja, two_je, two_jf)
          || triangle_selection_fails(two_jb, two_jd, two_jf)
          || triangle_selection_fails(two_je, two_jd, two_jc)
     ) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    return coupling_6j_recurrence(two_ja, two_jb, two_jc,
                                  two_jd, two_je, two_jf, result);
  }
}


int
gsl_sf_coupling_RacahW_e(int two_ja, int two_jb, int two_jc,
                         int two_jd, int two_je, int two_jf,
                         gsl_sf_result * result)
{
  int status = gsl_sf_coupling_6j_e(two_ja, two_jb, two_je, two_jd, two_jc, two_jf, result);
  int phase_sum = (two_ja + two_jb + two_jc + two_jd)/2;
  result->val *= ( GSL_IS_ODD(phase_sum) ? -1.0 : 1.0 );
  return status;
}


int
gsl_sf_coupling_9j_e(int two_ja, int two_jb, int two_jc,
                     int two_jd, int two_je, int two_jf,
                     int two_jg, int two_jh, int two_ji,
                     gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(   two_ja < 0 || two_jb < 0 || two_jc < 0
     || two_jd < 0 || two_je < 0 || two_jf < 0
     || two_jg < 0 || two_jh < 0 || two_ji < 0
     ) {
    DOMAIN_ERROR(result);
  }
  else if(   triangle_selection_fails(two_ja, two_jb, two_jc)
          || triangle_selection_fails(two_jd, two_je, two_jf)
          || triangle_selection_fails(two_jg, two_jh, two_ji)
          || triangle_selection_fails(two_ja, two_jd, two_jg)
          || triangle_selection_fails(two_jb, two_je, two_jh)
          || triangle_selection_fails(two_jc, two_jf, two_ji)
     ) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    int tk;
    int tkmin = locMax3(abs(two_ja-two_ji), abs(two_jh-two_jd), abs(two_jb-two_jf));
    int tkmax = locMin3(two_ja + two_ji, two_jh + two_jd, two_jb + two_jf);
    double sum_pos = 0.0;
    double sum_neg = 0.0;
    double sumsq_err = 0.0;
    double phase;
    for(tk=tkmin; tk<=tkmax; tk += 2) {
      gsl_sf_result s1, s2, s3;
      double term;
      double term_err;
      int status = 0;

      status += gsl_sf_coupling_6j_e(two_ja, two_ji, tk,  two_jh, two_jd, two_jg,  &s1);
      status += gsl_sf_coupling_6j_e(two_jb, two_jf, tk,  two_jd, two_jh, two_je,  &s2);
      status += gsl_sf_coupling_6j_e(two_ja, two_ji, tk,  two_jf, two_jb, two_jc,  &s3);

      if(status != GSL_SUCCESS) {
        OVERFLOW_ERROR(result);
      }
      term = s1.val * s2.val * s3.val;
      term_err  = s1.err * fabs(s2.val*s3.val);
      term_err += s2.err * fabs(s1.val*s3.val);
      term_err += s3.err * fabs(s1.val*s2.val);

      if(term >= 0.0) {
        sum_pos += (tk + 1) * term;
      }
      else {
        sum_neg -= (tk + 1) * term;
      }

      sumsq_err += ((tk+1) * term_err) * ((tk+1) * term_err);
    }

    phase = GSL_IS_ODD(tkmin) ? -1.0 : 1.0;

    result->val  = phase * (sum_pos - sum_neg);
    result->err  = 2.0 * GSL_DBL_EPSILON * (sum_pos + sum_neg);
    result->err += sqrt(sumsq_err / (0.5*(tkmax-tkmin)+1.0));
    result->err += 2.0 * GSL_DBL_EPSILON * (tkmax-tkmin + 2.0) * fabs(result->val);

    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double gsl_sf_coupling_3j(int two_ja, int two_jb, int two_jc,
                          int two_ma, int two_mb, int two_mc)
{
  EVAL_RESULT(gsl_sf_coupling_3j_e(two_ja, two_jb, two_jc,
                                   two_ma, two_mb, two_mc,
                                   &result));
}

#ifndef GSL_DISABLE_DEPRECATED

double gsl_sf_coupling_6j_INCORRECT(int two_ja, int two_jb, int two_jc,
                                    int two_jd, int two_je, int two_jf)
{
  EVAL_RESULT(gsl_sf_coupling_6j_INCORRECT_e(two_ja, two_jb, two_jc,
                                             two_jd, two_je, two_jf,
                                             &result));
}

#endif


double gsl_sf_coupling_6j(int two_ja, int two_jb, int two_jc,
                          int two_jd, int two_je, int two_jf)
{
  EVAL_RESULT(gsl_sf_coupling_6j_e(two_ja, two_jb, two_jc,
                                   two_jd, two_je, two_jf,
                                   &result));
}


double gsl_sf_coupling_RacahW(int two_ja, int two_jb, int two_jc,
                          int two_jd, int two_je, int two_jf)
{
  EVAL_RESULT(gsl_sf_coupling_RacahW_e(two_ja, two_jb, two_jc,
                                      two_jd, two_je, two_jf,
                                      &result));
}


double gsl_sf_coupling_9j(int two_ja, int two_jb, int two_jc,
                          int two_jd, int two_je, int two_jf,
                          int two_jg, int two_jh, int two_ji)
{
  EVAL_RESULT(gsl_sf_coupling_9j_e(two_ja, two_jb, two_jc,
                                   two_jd, two_je, two_jf,
                                   two_jg, two_jh, two_ji,
                                   &result));
}
