/* specfunc/hyperg_2F1.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 * Copyright (C) 2009 Brian Gough
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
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_hyperg.h>

#include "error.h"

#define locEPS (1000.0*GSL_DBL_EPSILON)


/* Assumes c != negative integer.
 */
static int
hyperg_2F1_series(const double a, const double b, const double c,
                  const double x, 
                  gsl_sf_result * result
                  )
{
  double sum_pos = 1.0;
  double sum_neg = 0.0;
  double del_pos = 1.0;
  double del_neg = 0.0;
  double del = 1.0;
  double absdel = 1.0;
  double absdel_prev = 0.0;
  double absdel_prev_prev = 0.0;
  double absval = 1.0;
  double err_est = 0.0;
  double k = 0.0;
  int i = 0;

  if(fabs(c) < GSL_DBL_EPSILON) {
    result->val = 0.0; /* FIXME: ?? */
    result->err = 1.0;
    GSL_ERROR ("error", GSL_EDOM);
  }

  do {
    if(++i > 30000) {
      result->val  = sum_pos - sum_neg;
      result->err  = del_pos + del_neg;
      result->err += 2.0 * GSL_DBL_EPSILON * (sum_pos + sum_neg);
      result->err += 2.0 * GSL_DBL_EPSILON * (2.0*sqrt(k)+1.0) * fabs(result->val);
      GSL_ERROR ("error", GSL_EMAXITER);
    }
    del *= (a+k)*(b+k) * x / ((c+k) * (k+1.0));  /* Gauss series */
    absdel_prev_prev = absdel_prev;
    absdel_prev = absdel;
    absdel = fabs(del);

    if(del > 0.0) {
      del_pos  =  del;
      sum_pos +=  del;
    }
    else if(del == 0.0) {
      /* Exact termination (a or b was a negative integer).
       */
      err_est = 0.0;
      break;
    }
    else {
      del_neg  = -del;
      sum_neg -=  del;
    }
    absval = fabs(sum_pos - sum_neg);

    /*
     * This stopping criteria is taken from the thesis
     * "Computation of Hypergeometic Functions" by J. Pearson, pg. 31
     * (http://people.maths.ox.ac.uk/porterm/research/pearson_final.pdf)
     * and fixes bug #45926
     */
    /* Pearson suggests a three-term stopping criteria for 2F1. */
    if (absdel_prev_prev < GSL_DBL_EPSILON*absval &&
        absdel_prev < GSL_DBL_EPSILON*absval &&
        absdel < GSL_DBL_EPSILON*absval )
      {
        /* Pick largest of three last terms as error estimate */
        err_est = GSL_MAX(absdel, GSL_MAX(absdel_prev, absdel_prev_prev));
        break;
      }

    k += 1.0;
    err_est = del_pos + del_neg;
  } while(err_est > GSL_DBL_EPSILON*absval);

  result->val = sum_pos - sum_neg;
  result->err = err_est;
  result->err += 2.0 * GSL_DBL_EPSILON * (sum_pos + sum_neg);
  result->err += 2.0 * GSL_DBL_EPSILON * (2.0*sqrt(k) + 1.0) * absval;

  return GSL_SUCCESS;
}


/* a = aR + i aI, b = aR - i aI */
static
int
hyperg_2F1_conj_series(const double aR, const double aI, const double c,
                       double x,
                       gsl_sf_result * result)
{
  if(c == 0.0) {
    result->val = 0.0; /* FIXME: should be Inf */
    result->err = 0.0;
    GSL_ERROR ("error", GSL_EDOM);
  }
  else {
    double sum_pos = 1.0;
    double sum_neg = 0.0;
    double del_pos = 1.0;
    double del_neg = 0.0;
    double del = 1.0;
    double k = 0.0;
    do {
      del *= ((aR+k)*(aR+k) + aI*aI)/((k+1.0)*(c+k)) * x;

      if(del >= 0.0) {
        del_pos  =  del;
        sum_pos +=  del;
      }
      else {
        del_neg  = -del;
        sum_neg -=  del;
      }

      if(k > 30000) {
        result->val  = sum_pos - sum_neg;
        result->err  = del_pos + del_neg;
        result->err += 2.0 * GSL_DBL_EPSILON * (sum_pos + sum_neg);
        result->err += 2.0 * GSL_DBL_EPSILON * (2.0*sqrt(k)+1.0) * fabs(result->val);
        GSL_ERROR ("error", GSL_EMAXITER);
      }

      k += 1.0;
    } while(fabs((del_pos + del_neg)/(sum_pos - sum_neg)) > GSL_DBL_EPSILON);

    result->val  = sum_pos - sum_neg;
    result->err  = del_pos + del_neg;
    result->err += 2.0 * GSL_DBL_EPSILON * (sum_pos + sum_neg);
    result->err += 2.0 * GSL_DBL_EPSILON * (2.0*sqrt(k) + 1.0) * fabs(result->val);

    return GSL_SUCCESS;
  }
}


/* Luke's rational approximation. The most accesible
 * discussion is in [Kolbig, CPC 23, 51 (1981)].
 * The convergence is supposedly guaranteed for x < 0.
 * You have to read Luke's books to see this and other
 * results. Unfortunately, the stability is not so
 * clear to me, although it seems very efficient when
 * it works.
 */
static
int
hyperg_2F1_luke(const double a, const double b, const double c,
                const double xin, 
                gsl_sf_result * result)
{
  int stat_iter;
  const double RECUR_BIG = 1.0e+50;
  const int nmax = 20000;
  int n = 3;
  const double x  = -xin;
  const double x3 = x*x*x;
  const double t0 = a*b/c;
  const double t1 = (a+1.0)*(b+1.0)/(2.0*c);
  const double t2 = (a+2.0)*(b+2.0)/(2.0*(c+1.0));
  double F = 1.0;
  double prec;

  double Bnm3 = 1.0;                                  /* B0 */
  double Bnm2 = 1.0 + t1 * x;                         /* B1 */
  double Bnm1 = 1.0 + t2 * x * (1.0 + t1/3.0 * x);    /* B2 */
 
  double Anm3 = 1.0;                                                      /* A0 */
  double Anm2 = Bnm2 - t0 * x;                                            /* A1 */
  double Anm1 = Bnm1 - t0*(1.0 + t2*x)*x + t0 * t1 * (c/(c+1.0)) * x*x;   /* A2 */

  while(1) {
    double npam1 = n + a - 1;
    double npbm1 = n + b - 1;
    double npcm1 = n + c - 1;
    double npam2 = n + a - 2;
    double npbm2 = n + b - 2;
    double npcm2 = n + c - 2;
    double tnm1  = 2*n - 1;
    double tnm3  = 2*n - 3;
    double tnm5  = 2*n - 5;
    double n2 = n*n;
    double F1 =  (3.0*n2 + (a+b-6)*n + 2 - a*b - 2*(a+b)) / (2*tnm3*npcm1);
    double F2 = -(3.0*n2 - (a+b+6)*n + 2 - a*b)*npam1*npbm1/(4*tnm1*tnm3*npcm2*npcm1);
    double F3 = (npam2*npam1*npbm2*npbm1*(n-a-2)*(n-b-2)) / (8*tnm3*tnm3*tnm5*(n+c-3)*npcm2*npcm1);
    double E  = -npam1*npbm1*(n-c-1) / (2*tnm3*npcm2*npcm1);

    double An = (1.0+F1*x)*Anm1 + (E + F2*x)*x*Anm2 + F3*x3*Anm3;
    double Bn = (1.0+F1*x)*Bnm1 + (E + F2*x)*x*Bnm2 + F3*x3*Bnm3;
    double r = An/Bn;

    prec = fabs((F - r)/F);
    F = r;

    if(prec < GSL_DBL_EPSILON || n > nmax) break;

    if(fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG) {
      An   /= RECUR_BIG;
      Bn   /= RECUR_BIG;
      Anm1 /= RECUR_BIG;
      Bnm1 /= RECUR_BIG;
      Anm2 /= RECUR_BIG;
      Bnm2 /= RECUR_BIG;
      Anm3 /= RECUR_BIG;
      Bnm3 /= RECUR_BIG;
    }
    else if(fabs(An) < 1.0/RECUR_BIG || fabs(Bn) < 1.0/RECUR_BIG) {
      An   *= RECUR_BIG;
      Bn   *= RECUR_BIG;
      Anm1 *= RECUR_BIG;
      Bnm1 *= RECUR_BIG;
      Anm2 *= RECUR_BIG;
      Bnm2 *= RECUR_BIG;
      Anm3 *= RECUR_BIG;
      Bnm3 *= RECUR_BIG;
    }

    n++;
    Bnm3 = Bnm2;
    Bnm2 = Bnm1;
    Bnm1 = Bn;
    Anm3 = Anm2;
    Anm2 = Anm1;
    Anm1 = An;
  }

  result->val  = F;
  result->err  = 2.0 * fabs(prec * F);
  result->err += 2.0 * GSL_DBL_EPSILON * (n+1.0) * fabs(F);

  /* FIXME: just a hack: there's a lot of shit going on here */
  result->err *= 8.0 * (fabs(a) + fabs(b) + 1.0);

  stat_iter = (n >= nmax ? GSL_EMAXITER : GSL_SUCCESS );

  return stat_iter;
}


/* Luke's rational approximation for the
 * case a = aR + i aI, b = aR - i aI.
 */
static
int
hyperg_2F1_conj_luke(const double aR, const double aI, const double c,
                     const double xin, 
                     gsl_sf_result * result)
{
  int stat_iter;
  const double RECUR_BIG = 1.0e+50;
  const int nmax = 10000;
  int n = 3;
  const double x = -xin;
  const double x3 = x*x*x;
  const double atimesb = aR*aR + aI*aI;
  const double apb     = 2.0*aR;
  const double t0 = atimesb/c;
  const double t1 = (atimesb +     apb + 1.0)/(2.0*c);
  const double t2 = (atimesb + 2.0*apb + 4.0)/(2.0*(c+1.0));
  double F = 1.0;
  double prec;

  double Bnm3 = 1.0;                                  /* B0 */
  double Bnm2 = 1.0 + t1 * x;                         /* B1 */
  double Bnm1 = 1.0 + t2 * x * (1.0 + t1/3.0 * x);    /* B2 */
 
  double Anm3 = 1.0;                                                      /* A0 */
  double Anm2 = Bnm2 - t0 * x;                                            /* A1 */
  double Anm1 = Bnm1 - t0*(1.0 + t2*x)*x + t0 * t1 * (c/(c+1.0)) * x*x;   /* A2 */

  while(1) {
    double nm1 = n - 1;
    double nm2 = n - 2;
    double npam1_npbm1 = atimesb + nm1*apb + nm1*nm1;
    double npam2_npbm2 = atimesb + nm2*apb + nm2*nm2;
    double npcm1 = nm1 + c;
    double npcm2 = nm2 + c;
    double tnm1  = 2*n - 1;
    double tnm3  = 2*n - 3;
    double tnm5  = 2*n - 5;
    double n2 = n*n;
    double F1 =  (3.0*n2 + (apb-6)*n + 2 - atimesb - 2*apb) / (2*tnm3*npcm1);
    double F2 = -(3.0*n2 - (apb+6)*n + 2 - atimesb)*npam1_npbm1/(4*tnm1*tnm3*npcm2*npcm1);
    double F3 = (npam2_npbm2*npam1_npbm1*(nm2*nm2 - nm2*apb + atimesb)) / (8*tnm3*tnm3*tnm5*(n+c-3)*npcm2*npcm1);
    double E  = -npam1_npbm1*(n-c-1) / (2*tnm3*npcm2*npcm1);

    double An = (1.0+F1*x)*Anm1 + (E + F2*x)*x*Anm2 + F3*x3*Anm3;
    double Bn = (1.0+F1*x)*Bnm1 + (E + F2*x)*x*Bnm2 + F3*x3*Bnm3;
    double r = An/Bn;

    prec = fabs(F - r)/fabs(F);
    F = r;

    if(prec < GSL_DBL_EPSILON || n > nmax) break;

    if(fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG) {
      An   /= RECUR_BIG;
      Bn   /= RECUR_BIG;
      Anm1 /= RECUR_BIG;
      Bnm1 /= RECUR_BIG;
      Anm2 /= RECUR_BIG;
      Bnm2 /= RECUR_BIG;
      Anm3 /= RECUR_BIG;
      Bnm3 /= RECUR_BIG;
    }
    else if(fabs(An) < 1.0/RECUR_BIG || fabs(Bn) < 1.0/RECUR_BIG) {
      An   *= RECUR_BIG;
      Bn   *= RECUR_BIG;
      Anm1 *= RECUR_BIG;
      Bnm1 *= RECUR_BIG;
      Anm2 *= RECUR_BIG;
      Bnm2 *= RECUR_BIG;
      Anm3 *= RECUR_BIG;
      Bnm3 *= RECUR_BIG;
    }

    n++;
    Bnm3 = Bnm2;
    Bnm2 = Bnm1;
    Bnm1 = Bn;
    Anm3 = Anm2;
    Anm2 = Anm1;
    Anm1 = An;
  }
  
  result->val  = F;
  result->err  = 2.0 * fabs(prec * F);
  result->err += 2.0 * GSL_DBL_EPSILON * (n+1.0) * fabs(F);

  /* FIXME: see above */
  result->err *= 8.0 * (fabs(aR) + fabs(aI) + 1.0);

  stat_iter = (n >= nmax ? GSL_EMAXITER : GSL_SUCCESS );

  return stat_iter;
}


/* Integer d = c-a-b.  The connection formula then carries logarithmic
 * terms; it is the limit of [A&S 15.3.6] as c-a-b -> m, given by
 * [A&S 15.3.10] for m = 0 and [A&S 15.3.11] for m >= 1:
 *
 *   F(a,b;a+b+m;z)
 *     = Gamma(m)Gamma(a+b+m)/(Gamma(a+m)Gamma(b+m))
 *         sum_{n=0}^{m-1} (a)_n(b)_n (m-n-1)!/n! (z-1)^n
 *     - Gamma(a+b+m)/(Gamma(a)Gamma(b)m!) (z-1)^m
 *         sum_{n=0}^inf (a+m)_n(b+m)_n/(n!(n+m)!)
 *           (1-z)^n [ln(1-z) - psi(n+1) - psi(n+m+1)
 *                    + psi(a+n+m) + psi(b+n+m)].
 *
 * The factor 1/m! on the second term comes from the residue of
 * Gamma(-d) at d = -m; the formula as usually printed omits it and is
 * then only correct for m = 0 and m = 1.  Each gamma factor carries its
 * sign explicitly, which the Moshier construction did not do.
 *
 * Assumes m = rint(c-a-b) >= 0, x >= 0.995, and that the caller's
 * dispatch has excluded the pole cases (a, b, c and a+m, b+m not
 * non-positive integers).
 */
static int
hyperg_2F1_reflect_dint(const double a, const double b, const double c,
                        const int m, const double x, gsl_sf_result * result)
{
  const int maxiter = 2000;
  const double omx = 1.0 - x;
  const double ln_omx = log(omx);
  const double zm1 = x - 1.0;

  gsl_sf_result lng_c, lng_a, lng_b, lng_am, lng_bm;
  double sg_c, sg_a, sg_b, sg_am, sg_bm;
  int stat = 0;

  stat |= gsl_sf_lngamma_sgn_e(c,     &lng_c,  &sg_c);
  stat |= gsl_sf_lngamma_sgn_e(a,     &lng_a,  &sg_a);
  stat |= gsl_sf_lngamma_sgn_e(b,     &lng_b,  &sg_b);
  stat |= gsl_sf_lngamma_sgn_e(a + m, &lng_am, &sg_am);
  stat |= gsl_sf_lngamma_sgn_e(b + m, &lng_bm, &sg_bm);

  if(m == 0) {
    /* [A&S 15.3.10] */
    const double ln_pre = lng_c.val - lng_a.val - lng_b.val;
    const double ln_pre_err = lng_c.err + lng_a.err + lng_b.err
                            + GSL_DBL_EPSILON * fabs(ln_pre);
    const double sgn_pre = sg_c * sg_a * sg_b;
    double term = 1.0;
    double sum = 0.0;
    double sum_err = 0.0;
    double abs_sum = 0.0;
    gsl_sf_result pre;
    int stat_e;
    int n;

    for(n = 0; n < maxiter; n++) {
      gsl_sf_result psi_n1, psi_an, psi_bn;
      double bracket, term_val;

      stat |= gsl_sf_psi_e((double)(n+1), &psi_n1);
      stat |= gsl_sf_psi_e(a + n,        &psi_an);
      stat |= gsl_sf_psi_e(b + n,        &psi_bn);

      bracket = 2.0*psi_n1.val - psi_an.val - psi_bn.val - ln_omx;
      term_val = term * bracket;
      sum += term_val;
      abs_sum += fabs(term_val);
      sum_err += fabs(term) * (psi_n1.err + psi_an.err + psi_bn.err);

      if(n > 20 && fabs(term) < GSL_DBL_EPSILON * fabs(sum)) break;
      term *= (a + n) * (b + n) / ((double)(n+1) * (double)(n+1)) * omx;
    }
    sum_err += 2.0 * GSL_DBL_EPSILON * abs_sum;

    stat_e = gsl_sf_exp_err_e(ln_pre, ln_pre_err, &pre);
    pre.val *= sgn_pre;
    result->val = pre.val * sum;
    result->err = fabs(pre.val) * sum_err + fabs(sum) * pre.err
                + 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat, stat_e);
  }
  else {
    /* [A&S 15.3.11] */
    double lgamma_m;
    double fact_m_inv = 1.0;
    double ln_pre1, ln_pre1_err, ln_pre2, ln_pre2_err;
    double sgn_pre1, sgn_pre2;
    double term, sum1 = 0.0, sum1_err = 0.0, abs_sum1 = 0.0;
    double sum2 = 0.0, sum2_err = 0.0, abs_sum2 = 0.0;
    gsl_sf_result pre1, pre2;
    int stat_e1, stat_e2;
    int n;

    lgamma_m = gsl_sf_lngamma((double) m);

    /* P1 = Gamma(c)/(Gamma(a+m)Gamma(b+m)); the Gamma(m) of the
     * printed formula is already carried by the finite sum, whose
     * first term is (m-1)!. */
    ln_pre1 = lng_c.val - lng_am.val - lng_bm.val;
    ln_pre1_err = lng_c.err + lng_am.err + lng_bm.err
                + 2.0 * GSL_DBL_EPSILON * fabs(ln_pre1);
    sgn_pre1 = sg_c * sg_am * sg_bm;

    /* P2 = -Gamma(c)/(Gamma(a)Gamma(b)) (z-1)^m.  The 1/m! of the
     * printed second series is carried by its first term. */
    ln_pre2 = lng_c.val - lng_a.val - lng_b.val
            + m * log(fabs(zm1));
    ln_pre2_err = lng_c.err + lng_a.err + lng_b.err
                + 2.0 * GSL_DBL_EPSILON * (fabs(ln_pre2) + m);
    sgn_pre2 = -sg_c * sg_a * sg_b;
    if(GSL_IS_ODD(m)) sgn_pre2 = -sgn_pre2;

    /* Finite sum: term_n = (a)_n(b)_n (m-n-1)!/n! (z-1)^n. */
    term = exp(lgamma_m);  /* (m-1)! */
    for(n = 0; n < m; n++) {
      sum1 += term;
      abs_sum1 += fabs(term);
      if(n + 1 < m)
        term *= (a + n) * (b + n) / ((double)(m - n - 1) * (double)(n + 1)) * zm1;
    }
    sum1_err = 2.0 * GSL_DBL_EPSILON * abs_sum1;

    /* Infinite sum: term_n = (a+m)_n(b+m)_n/(n!(n+m)!) (1-z)^n. */
    for(n = 0; n <= m; n++)  /* 1/m! */
      fact_m_inv /= (double)(n == 0 ? 1 : n);
    term = fact_m_inv;
    for(n = 0; n < maxiter; n++) {
      gsl_sf_result psi_n1, psi_nm1, psi_anm, psi_bnm;
      double bracket, delta;

      stat |= gsl_sf_psi_e((double)(n+1),     &psi_n1);
      stat |= gsl_sf_psi_e((double)(n+m+1),   &psi_nm1);
      stat |= gsl_sf_psi_e(a + n + m,         &psi_anm);
      stat |= gsl_sf_psi_e(b + n + m,         &psi_bnm);

      bracket = ln_omx - psi_n1.val - psi_nm1.val + psi_anm.val + psi_bnm.val;
      delta = term * bracket;
      sum2 += delta;
      abs_sum2 += fabs(delta);
      sum2_err += fabs(term) * (psi_n1.err + psi_nm1.err + psi_anm.err + psi_bnm.err);

      if(n > 20 && fabs(term) < GSL_DBL_EPSILON * fabs(sum2)) break;
      term *= (a + m + n) * (b + m + n)
            / ((double)(n+1) * (double)(n + m + 1)) * omx;
    }
    sum2_err += 2.0 * GSL_DBL_EPSILON * abs_sum2;

    stat_e1 = gsl_sf_exp_err_e(ln_pre1, ln_pre1_err, &pre1);
    stat_e2 = gsl_sf_exp_err_e(ln_pre2, ln_pre2_err, &pre2);
    pre1.val *= sgn_pre1;
    pre2.val *= sgn_pre2;

    result->val = pre1.val * sum1 + pre2.val * sum2;
    result->err = fabs(pre1.val) * sum1_err + fabs(pre2.val) * sum2_err
                + fabs(sum1) * pre1.err + fabs(sum2) * pre2.err
                + 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_4(stat, stat_e1, stat_e2, GSL_SUCCESS);
  }
}


static int pow_omx(const double x, const double p, gsl_sf_result * result);


/* Do the reflection described in [Moshier, p. 334].
 * Assumes a,b,c != neg integer.
 */
static
int
hyperg_2F1_reflect(const double a, const double b, const double c,
                   const double x, gsl_sf_result * result)
{
  const double d = c - a - b;
  const int intd  = rint(d);
  const int d_integer = ( fabs(d - intd) < locEPS );

  if(d_integer) {
    if(intd == 0) {
      return hyperg_2F1_reflect_dint(a, b, c, 0, x, result);
    }
    else if(intd > 0) {
      return hyperg_2F1_reflect_dint(a, b, c, intd, x, result);
    }
    else {
      /* c-a-b = -m < 0: apply the [A&S 15.3.12] shift (1-z)^m first,
       * which turns it into the m > 0 case. */
      const int m = -intd;
      const double am = a - m;
      const double bm = b - m;
      const double rintam = rint(am);
      const double rintbm = rint(bm);
      gsl_sf_result F;
      gsl_sf_result p;
      int stat_F;
      int stat_p;
      if(   (am <= 0.0 && fabs(am - rintam) < locEPS)
         || (bm <= 0.0 && fabs(bm - rintbm) < locEPS)) {
        /* The shifted function terminates: a-m or b-m is zero or a
         * negative integer.  hyperg_2F1_reflect_dint() forms its
         * prefactors from lngamma at those arguments, which are poles
         * there, so evaluate the terminating series directly instead. */
        stat_F = hyperg_2F1_series(am, bm, c, x, &F);
      }
      else {
        stat_F = hyperg_2F1_reflect_dint(am, bm, c, m, x, &F);
      }
      stat_p = pow_omx(x, -(double)m, &p);
      result->val = p.val * F.val;
      result->err = fabs(p.val) * F.err + fabs(F.val) * p.err
                  + 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_ERROR_SELECT_2(stat_F, stat_p);
    }
  }
  else {
    /* d not an integer */

    gsl_sf_result pre1, pre2;
    double sgn1, sgn2;
    gsl_sf_result F1, F2;
    int status_F1, status_F2;

    /* These gamma functions appear in the denominator, so we
     * catch their harmless domain errors and set the terms to zero.
     */
    gsl_sf_result ln_g1ca,  ln_g1cb,  ln_g2a,  ln_g2b;
    double sgn_g1ca, sgn_g1cb, sgn_g2a, sgn_g2b;
    int stat_1ca = gsl_sf_lngamma_sgn_e(c-a, &ln_g1ca, &sgn_g1ca);
    int stat_1cb = gsl_sf_lngamma_sgn_e(c-b, &ln_g1cb, &sgn_g1cb);
    int stat_2a  = gsl_sf_lngamma_sgn_e(a, &ln_g2a, &sgn_g2a);
    int stat_2b  = gsl_sf_lngamma_sgn_e(b, &ln_g2b, &sgn_g2b);
    int ok1 = (stat_1ca == GSL_SUCCESS && stat_1cb == GSL_SUCCESS);
    int ok2 = (stat_2a  == GSL_SUCCESS && stat_2b  == GSL_SUCCESS);
    
    gsl_sf_result ln_gc,  ln_gd,  ln_gmd;
    double sgn_gc, sgn_gd, sgn_gmd;
    gsl_sf_lngamma_sgn_e( c, &ln_gc,  &sgn_gc);
    gsl_sf_lngamma_sgn_e( d, &ln_gd,  &sgn_gd);
    gsl_sf_lngamma_sgn_e(-d, &ln_gmd, &sgn_gmd);
    
    sgn1 = sgn_gc * sgn_gd  * sgn_g1ca * sgn_g1cb;
    sgn2 = sgn_gc * sgn_gmd * sgn_g2a  * sgn_g2b;

    if(ok1 && ok2) {
      double ln_pre1_val = ln_gc.val + ln_gd.val  - ln_g1ca.val - ln_g1cb.val;
      double ln_pre2_val = ln_gc.val + ln_gmd.val - ln_g2a.val  - ln_g2b.val + d*log(1.0-x);
      double ln_pre1_err = ln_gc.err + ln_gd.err + ln_g1ca.err + ln_g1cb.err;
      double ln_pre2_err = ln_gc.err + ln_gmd.err + ln_g2a.err  + ln_g2b.err;
      if(ln_pre1_val < GSL_LOG_DBL_MAX && ln_pre2_val < GSL_LOG_DBL_MAX) {
        gsl_sf_exp_err_e(ln_pre1_val, ln_pre1_err, &pre1);
        gsl_sf_exp_err_e(ln_pre2_val, ln_pre2_err, &pre2);
        pre1.val *= sgn1;
        pre2.val *= sgn2;
      }
      else {
        OVERFLOW_ERROR(result);
      }
    }
    else if(ok1 && !ok2) {
      double ln_pre1_val = ln_gc.val + ln_gd.val - ln_g1ca.val - ln_g1cb.val;
      double ln_pre1_err = ln_gc.err + ln_gd.err + ln_g1ca.err + ln_g1cb.err;
      if(ln_pre1_val < GSL_LOG_DBL_MAX) {
        gsl_sf_exp_err_e(ln_pre1_val, ln_pre1_err, &pre1);
        pre1.val *= sgn1;
        pre2.val = 0.0;
        pre2.err = 0.0;
      }
      else {
        OVERFLOW_ERROR(result);
      }
    }
    else if(!ok1 && ok2) {
      double ln_pre2_val = ln_gc.val + ln_gmd.val - ln_g2a.val - ln_g2b.val + d*log(1.0-x);
      double ln_pre2_err = ln_gc.err + ln_gmd.err + ln_g2a.err + ln_g2b.err;
      if(ln_pre2_val < GSL_LOG_DBL_MAX) {
        pre1.val = 0.0;
        pre1.err = 0.0;
        gsl_sf_exp_err_e(ln_pre2_val, ln_pre2_err, &pre2);
        pre2.val *= sgn2;
      }
      else {
        OVERFLOW_ERROR(result);
      }
    }
    else {
      pre1.val = 0.0;
      pre2.val = 0.0;
      UNDERFLOW_ERROR(result);
    }

    status_F1 = hyperg_2F1_series(  a,   b, 1.0-d, 1.0-x, &F1);
    status_F2 = hyperg_2F1_series(c-a, c-b, 1.0+d, 1.0-x, &F2);

    result->val  = pre1.val*F1.val + pre2.val*F2.val;
    result->err  = fabs(pre1.val*F1.err) + fabs(pre2.val*F2.err);
    result->err += fabs(pre1.err*F1.val) + fabs(pre2.err*F2.val);
    result->err += 2.0 * GSL_DBL_EPSILON * (fabs(pre1.val*F1.val) + fabs(pre2.val*F2.val));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

    if (status_F1)
      return status_F1;

    if (status_F2)
      return status_F2;

    return GSL_SUCCESS;
  }
}


static int pow_omx(const double x, const double p, gsl_sf_result * result)
{
  double ln_omx;
  double ln_result;
  if(fabs(x) < GSL_ROOT5_DBL_EPSILON) {
    ln_omx = -x*(1.0 + x*(1.0/2.0 + x*(1.0/3.0 + x/4.0 + x*x/5.0)));
  }
  else {
    ln_omx = log(1.0-x);
  }
  ln_result = p * ln_omx;
  return gsl_sf_exp_err_e(ln_result, GSL_DBL_EPSILON * fabs(ln_result), result);
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_hyperg_2F1_e(double a, double b, const double c,
                       const double x,
                       gsl_sf_result * result)
{
  const double d = c - a - b;
  const double rinta = rint(a);
  const double rintb = rint(b);
  const double rintc = rint(c);
  const double rintd = rint(d);
  const int a_neg_integer = ( a < 0.0  &&  fabs(a - rinta) < locEPS );
  const int b_neg_integer = ( b < 0.0  &&  fabs(b - rintb) < locEPS );
  const int c_neg_integer = ( c < 0.0  &&  fabs(c - rintc) < locEPS );
  const int d_integer     = ( fabs(d - rintd) < locEPS );

  result->val = 0.0;
  result->err = 0.0;

   /* Handle x == 1.0 RJM */

  if (fabs (x - 1.0) < locEPS && (c - a - b) > 0 && c != 0 && !c_neg_integer) {
    gsl_sf_result lngamc, lngamcab, lngamca, lngamcb;
    double lngamc_sgn, lngamca_sgn, lngamcb_sgn;
    int status;
    int stat1 = gsl_sf_lngamma_sgn_e (c, &lngamc, &lngamc_sgn);
    int stat2 = gsl_sf_lngamma_e (c - a - b, &lngamcab);
    int stat3 = gsl_sf_lngamma_sgn_e (c - a, &lngamca, &lngamca_sgn);
    int stat4 = gsl_sf_lngamma_sgn_e (c - b, &lngamcb, &lngamcb_sgn);
    
    if (stat1 != GSL_SUCCESS || stat2 != GSL_SUCCESS
        || stat3 != GSL_SUCCESS || stat4 != GSL_SUCCESS)
      {
        DOMAIN_ERROR (result);
      }
    
    status =
      gsl_sf_exp_err_e (lngamc.val + lngamcab.val - lngamca.val - lngamcb.val,
                        lngamc.err + lngamcab.err + lngamca.err + lngamcb.err,
                        result);
    
    result->val *= lngamc_sgn / (lngamca_sgn * lngamcb_sgn);
      return status;
  }
  
  /* 2F1 reduces to a polynomial when a or b is zero or a negative
   * integer.  That polynomial is finite for every x, including |x| >= 1
   * where the defining Gauss series does not converge, so evaluate it
   * with the terminating series before the convergence-domain check.
   * The existing dispatch only reached the terminating series through
   * the |a|,|b| < 10 branch below, so e.g. 2F1(-1,-10,1,1/2) was
   * rejected with GSL_EUNIMPL (Savannah bugs #50711, #53905, #21835).
   */
  if(a == 0.0 || b == 0.0 || a_neg_integer || b_neg_integer) {
    const double ap = (a_neg_integer ? rinta : a);
    const double bp = (b_neg_integer ? rintb : b);

    if(c_neg_integer) {
      /* c is a negative integer; a or b must terminate before the
       * factor (c + k) vanishes, exactly as in the general path. */
      if(! (a_neg_integer && a > c + 0.1) && ! (b_neg_integer && b > c + 0.1)) {
        DOMAIN_ERROR(result);
      }
    }

    return hyperg_2F1_series(ap, bp, c, x, result);
  }

  /* The Gauss series only converges for |x| < 1, but the Pfaff
   * transformations [DLMF 15.8.1 and 15.8.2]
   *
   *   2F1(a,b;c;x) = (1-x)^-a 2F1(a,   c-b; c, x/(x-1))
   *                = (1-x)^-b 2F1(c-a, b;   c, x/(x-1))
   *
   * continue the function to every x < 1.  For x < -1 the transformed
   * argument z = x/(x-1) lies in (1/2,1), where the dispatch below is
   * at its best.  Savannah bug #30324.
   */
  if(x < -1.0) {
    const double z = x / (x - 1.0);
    const double ln_omx = log(1.0 - x);
    double a1, b1, p;
    gsl_sf_result F;
    int stat_F;
    int use_a;

    /* Prefer a form whose transformed parameter terminates the series
     * exactly; otherwise use the smaller prefactor exponent.  The
     * sign-exchange 2F1(a,b;c;x) = 2F1(b,a;c;x) has already been
     * folded into the a/b naming by the caller, so either is valid.
     */
    if(fabs(c - b - rint(c - b)) < locEPS && c - b <= 0.0) {
      use_a = 1;
    }
    else if(fabs(c - a - rint(c - a)) < locEPS && c - a <= 0.0) {
      use_a = 0;
    }
    else {
      use_a = (fabs(a) <= fabs(b));
    }

    if(use_a) {
      a1 = a;
      b1 = c - b;
      p  = -a;
    }
    else {
      a1 = c - a;
      b1 = b;
      p  = -b;
    }

    stat_F = gsl_sf_hyperg_2F1_e(a1, b1, c, z, &F);
    {
      const double ln_pre = p * ln_omx;
      const double ln_pre_err = GSL_DBL_EPSILON * (fabs(ln_pre) + fabs(ln_omx));
      const int stat_e = gsl_sf_exp_mult_err_e(ln_pre, ln_pre_err,
                                               F.val, F.err, result);
      return GSL_ERROR_SELECT_2(stat_e, stat_F);
    }
  }

  if(1.0 <= x) {
    DOMAIN_ERROR(result);
  }

  if(c_neg_integer) {
    /* If c is a negative integer, then either a or b must be a
       negative integer of smaller magnitude than c to ensure
       cancellation of the series. */
    if(! (a_neg_integer && a > c + 0.1) && ! (b_neg_integer && b > c + 0.1)) {
      DOMAIN_ERROR(result);
    }
  }

  /* Contiguous relations for the shifts c = a-1 and c = b-1.  Solving
     DLMF 15.5.17
       (a-1+(b+1-c)z) F(a,b;c;z) + (c-a) F(a-1,b;c;z)
         - (c-1)(1-z) F(a,b;c-1;z) = 0
     for F(a,b;c-1), then setting c = a and using F(a,b;a;z) = (1-z)^-b,
     gives the closed forms

       2F1(a, b, a-1, x) = (a-1 + (b+1-a) x) / ((a-1) (1-x)^(b+1))
       2F1(a, b, b-1, x) = (b-1 + (a+1-b) x) / ((b-1) (1-x)^(a+1))

     valid for all |x| < 1.  The previous code reached these through the
     general path, where d = c-a-b+1 = 0 made the reflection formula fail
     and returned a spurious 0 (or GSL_EUNIMPL) for c < 0.  The a = 0 and
     b = 0 cases were handled above.  See Savannah bug #66850. */
  if(fabs(c-(a-1.0)) < locEPS) {
    gsl_sf_result p;
    gsl_sf_result val;
    const double num = (a-1.0) + (b+1.0-a)*x;
    int stat_p = pow_omx(x, b+1.0, &p);
    val.val = num / ((a-1.0) * p.val);
    val.err = fabs(val.val) * (p.err / fabs(p.val));
    val.err += 2.0 * GSL_DBL_EPSILON * fabs(val.val);
    result->val = val.val;
    result->err = val.err;
    return stat_p;
  }
  if(fabs(c-(b-1.0)) < locEPS) {
    gsl_sf_result p;
    gsl_sf_result val;
    const double num = (b-1.0) + (a+1.0-b)*x;
    int stat_p = pow_omx(x, a+1.0, &p);
    val.val = num / ((b-1.0) * p.val);
    val.err = fabs(val.val) * (p.err / fabs(p.val));
    val.err += 2.0 * GSL_DBL_EPSILON * fabs(val.val);
    result->val = val.val;
    result->err = val.err;
    return stat_p;
  }

  if(fabs(c-a) < locEPS) {
    return pow_omx(x, d, result);  /* 2F1(a,b,a,x) = (1-x)^(c-a-b) */
  }

  /* When c-a-b is an integer the connection formula carries logarithmic
   * terms.  For x close to 1 the Gauss series converges too slowly to
   * meet its iteration cap, but the [A&S 15.3.10] (d = 0), [A&S 15.3.11]
   * (d = m > 0) and [A&S 15.3.12] (d = -m < 0) limits implemented in
   * hyperg_2F1_reflect() expand in (1-x) and converge rapidly there.
   * Savannah bug #21835.
   *
   * Two guards keep the reflection within the region where it is
   * accurate.  |d| bounds the (m-1)! finite sum of the m > 0 form; the
   * product (|a|+|d|)(|b|+|d|)(1-x) bounds the growth of the (1-x)
   * terms relative to their sum, so that the cancellation stays below
   * about 1e-12.  Outside it the existing dispatch (Gauss series, Luke
   * or the large-parameter branches) is retained unchanged. */
  if(d_integer && x >= 0.995
     && fabs(rintd) < 100.0
     && GSL_MAX_DBL(1.0, (fabs(a)+fabs(rintd))*(fabs(b)+fabs(rintd)))
          * (1.0 - x) < 20.0) {
    return hyperg_2F1_reflect(a, b, c, x, result);
  }

  if(a >= 0.0 && b >= 0.0 && c >=0.0 && x >= 0.0 && x < 0.995) {
    /* Series has all positive definite
     * terms and x is not close to 1.
     */
    return hyperg_2F1_series(a, b, c, x, result);
  }

  if(fabs(a) < 10.0 && fabs(b) < 10.0) {
    /* a and b are not too large, so we attempt
     * variations on the series summation.
     */
    if(a_neg_integer) {
      return hyperg_2F1_series(rinta, b, c, x, result);
    }
    if(b_neg_integer) {
      return hyperg_2F1_series(a, rintb, c, x, result);
    }

    if(x < -0.25) {
      return hyperg_2F1_luke(a, b, c, x, result);
    }
    else if(x < 0.5) {
      return hyperg_2F1_series(a, b, c, x, result);
    }
    else {
      if(fabs(c) > 10.0 || (d_integer && x < 0.995)) {
        /* The reflection formula is not reliable when c-a-b is an
         * integer; use the Gauss series instead while it still
         * converges quickly. */
        return hyperg_2F1_series(a, b, c, x, result);
      }
      else {
        return hyperg_2F1_reflect(a, b, c, x, result);
      }
    }
  }
  else {
    /* Either a or b or both large.
     * Introduce some new variables ap,bp so that bp is
     * the larger in magnitude.
     */
    double ap, bp; 
    if(fabs(a) > fabs(b)) {
      bp = a;
      ap = b;
    }
    else {
      bp = b;
      ap = a;
    }

    if(x < 0.0) {
      /* What the hell, maybe Luke will converge.
       */
      return hyperg_2F1_luke(a, b, c, x, result);
    }

    if(GSL_MAX_DBL(fabs(ap),1.0)*fabs(bp)*fabs(x) < 2.0*fabs(c)) {
      /* If c is large enough or x is small enough,
       * we can attempt the series anyway.
       */
      return hyperg_2F1_series(a, b, c, x, result);
    }

    if(fabs(bp*bp*x*x) < 0.001*fabs(bp) && fabs(ap) < 10.0) {
      /* The famous but nearly worthless "large b" asymptotic.
       */
      int stat = gsl_sf_hyperg_1F1_e(ap, c, bp*x, result);
      result->err = 0.001 * fabs(result->val);
      return stat;
    }

    /* We give up. */
    result->val = 0.0;
    result->err = 0.0;
    GSL_ERROR ("error", GSL_EUNIMPL);
  }
}


/* Pfaff continuation of the conjugate 2F1 to x < -1:
 *
 *   2F1(a,a*;c;x) = (1-x)^-a 2F1(a, c-a*; c, x/(x-1)),  a = aR + i aI.
 *
 * The transformed pair (aR + i aI, (c-aR) + i aI) is not a conjugate
 * pair, so the Gauss series is summed directly in complex arithmetic and
 * the real part of the product with the complex prefactor is returned
 * (the imaginary part vanishes identically).  Savannah bug #30324.
 */
static int
hyperg_2F1_conj_xlt_m1(const double aR, const double aI, const double c,
                       const double x, gsl_sf_result * result)
{
  const int max_iter = 30000;
  const double z = x / (x - 1.0);
  const double bR = c - aR;
  double sum_re = 1.0;
  double sum_im = 0.0;
  double term_re = 1.0;
  double term_im = 0.0;
  double abs_sum = 1.0;
  double last_mag = 0.0;
  int k;

  if(aI == 0.0) {
    /* The pair is real; the real routine continues it. */
    return gsl_sf_hyperg_2F1_e(aR, aR, c, x, result);
  }

  for(k = 0; k < max_iter; k++) {
    const double ar = aR + k;
    const double br = bR + k;
    const double den = (c + k) * (k + 1.0);
    const double f_re = (ar * br - aI * aI) / den * z;
    const double f_im = aI * (ar + br) / den * z;
    const double t_re = term_re * f_re - term_im * f_im;
    const double t_im = term_re * f_im + term_im * f_re;

    term_re = t_re;
    term_im = t_im;
    sum_re += term_re;
    sum_im += term_im;

    last_mag = fabs(term_re) + fabs(term_im);
    abs_sum += last_mag;

    if(!gsl_finite(last_mag)) {
      result->val = 0.0;
      result->err = 0.0;
      GSL_ERROR ("error", GSL_EMAXITER);
    }

    if(last_mag <= GSL_DBL_EPSILON * (fabs(sum_re) + fabs(sum_im)))
      break;
  }

  if(k >= max_iter) {
    result->val = 0.0;
    result->err = 0.0;
    GSL_ERROR ("error", GSL_EMAXITER);
  }

  {
    const double L = log(1.0 - x);
    const double mag_p = exp(-aR * L);
    const double phase_re = cos(aI * L);
    const double phase_im = -sin(aI * L);
    const double H_abs = fabs(sum_re) + fabs(sum_im);
    const double err_H = 2.0 * GSL_DBL_EPSILON * abs_sum + last_mag;
    const double err_p = GSL_DBL_EPSILON
      * (fabs(aR * L) + fabs(aI * L) + 1.0) * mag_p;
    const double val = mag_p * (phase_re * sum_re - phase_im * sum_im);
    double err = mag_p * err_H + err_p * H_abs;

    err += 2.0 * GSL_DBL_EPSILON * fabs(val);
    result->val = val;
    result->err = err;
  }

  return GSL_SUCCESS;
}


int
gsl_sf_hyperg_2F1_conj_e(const double aR, const double aI, const double c,
                            const double x,
                            gsl_sf_result * result)
{
  const double ax = fabs(x);
  const double rintc = rint(c);
  const int c_neg_integer = ( c < 0.0  &&  fabs(c - rintc) < locEPS );

  result->val = 0.0;
  result->err = 0.0;

  if(c_neg_integer || c == 0.0) {
    DOMAIN_ERROR(result);
  }

  /* Continue to x < -1 with the Pfaff transformation, as for the real
   * 2F1 (Savannah bug #30324).  The transformed pair is not conjugate,
   * so this goes through a complex series; see below. */
  if(x < -1.0) {
    return hyperg_2F1_conj_xlt_m1(aR, aI, c, x, result);
  }

  if(ax >= 1.0) {
    DOMAIN_ERROR(result);
  }

  if(   (ax < 0.25 && fabs(aR) < 20.0 && fabs(aI) < 20.0)
     || (c > 0.0 && x > 0.0)
    ) {
    return hyperg_2F1_conj_series(aR, aI, c, x, result);
  }
  else if(fabs(aR) < 10.0 && fabs(aI) < 10.0) {
    if(x < -0.25) {
      return hyperg_2F1_conj_luke(aR, aI, c, x, result);
    }
    else {
      return hyperg_2F1_conj_series(aR, aI, c, x, result);
    }
  }
  else {
    if(x < 0.0) {
      /* What the hell, maybe Luke will converge.
       */
      return hyperg_2F1_conj_luke(aR, aI, c, x, result); 
    }

    /* Give up. */
    result->val = 0.0;
    result->err = 0.0;
    GSL_ERROR ("error", GSL_EUNIMPL);
  }
}


int
gsl_sf_hyperg_2F1_renorm_e(const double a, const double b, const double c,
                              const double x,
                              gsl_sf_result * result
                              )
{
  const double rinta = rint(a);
  const double rintb = rint(b);
  const double rintc = rint(c);
  const int a_neg_integer = ( a < 0.0  &&  fabs(a - rinta) < locEPS );
  const int b_neg_integer = ( b < 0.0  &&  fabs(b - rintb) < locEPS );
  const int c_neg_integer = ( c < 0.0  &&  fabs(c - rintc) < locEPS );
  
  if(c_neg_integer) {
    if((a_neg_integer && a > c+0.1) || (b_neg_integer && b > c+0.1)) {
      /* 2F1 terminates early */
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else {
      /* 2F1 does not terminate early enough, so something survives */
      /* [Abramowitz+Stegun, 15.1.2] */
      gsl_sf_result g1, g2, g3, g4, g5;
      double s1, s2, s3, s4, s5;
      int stat = 0;
      stat += gsl_sf_lngamma_sgn_e(a-c+1, &g1, &s1);
      stat += gsl_sf_lngamma_sgn_e(b-c+1, &g2, &s2);
      stat += gsl_sf_lngamma_sgn_e(a, &g3, &s3);
      stat += gsl_sf_lngamma_sgn_e(b, &g4, &s4);
      stat += gsl_sf_lngamma_sgn_e(-c+2, &g5, &s5);
      if(stat != 0) {
        DOMAIN_ERROR(result);
      }
      else {
        /* [Abramowitz+Stegun, 15.1.2] carries a factor x^(1-c), which the
         * original code dropped; c is a negative integer here, so 1-c is
         * a positive integer and gsl_sf_pow_int_e applies
         * (Savannah bug #53876). */
        gsl_sf_result F;
        gsl_sf_result powx;
        int stat_F = gsl_sf_hyperg_2F1_e(a-c+1, b-c+1, -c+2, x, &F);
        int stat_p = gsl_sf_pow_int_e(x, 1 - (int)rintc, &powx);
        double ln_pre_val = g1.val + g2.val - g3.val - g4.val - g5.val;
        double ln_pre_err = g1.err + g2.err + g3.err + g4.err + g5.err;
        double sg  = s1 * s2 * s3 * s4 * s5;
        double pre_val = sg * powx.val * F.val;
        double pre_err = fabs(sg * powx.val) * F.err
                       + fabs(sg * F.val) * powx.err
                       + 2.0 * GSL_DBL_EPSILON * fabs(pre_val);
        int stat_e = gsl_sf_exp_mult_err_e(ln_pre_val, ln_pre_err,
                                              pre_val, pre_err,
                                              result);
        return GSL_ERROR_SELECT_3(stat_e, stat_F, stat_p);
      }
    }
  }
  else {
    /* generic c */
    gsl_sf_result F;
    gsl_sf_result lng;
    double sgn;
    int stat_g = gsl_sf_lngamma_sgn_e(c, &lng, &sgn);
    int stat_F = gsl_sf_hyperg_2F1_e(a, b, c, x, &F);
    int stat_e = gsl_sf_exp_mult_err_e(-lng.val, lng.err,
                                          sgn*F.val, F.err,
                                          result);
    return GSL_ERROR_SELECT_3(stat_e, stat_F, stat_g);
  }
}


int
gsl_sf_hyperg_2F1_conj_renorm_e(const double aR, const double aI, const double c,
                                   const double x,
                                   gsl_sf_result * result
                                   )
{
  const double rintc = rint(c);
  const double rinta = rint(aR);
  const int a_neg_integer = ( aR < 0.0 && fabs(aR-rinta) < locEPS && aI == 0.0);
  const int c_neg_integer = (  c < 0.0 && fabs(c - rintc) < locEPS );

  if(c_neg_integer) {
    if(a_neg_integer && aR > c+0.1) {
      /* 2F1 terminates early */
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else {
      /* 2F1 does not terminate early enough, so something survives */
      /* [Abramowitz+Stegun, 15.1.2] */
      gsl_sf_result g1, g2;
      gsl_sf_result g3;
      gsl_sf_result a1, a2;
      int stat = 0;
      stat += gsl_sf_lngamma_complex_e(aR-c+1, aI, &g1, &a1);
      stat += gsl_sf_lngamma_complex_e(aR, aI, &g2, &a2);
      stat += gsl_sf_lngamma_e(-c+2.0, &g3);
      if(stat != 0) {
        DOMAIN_ERROR(result);
      }
      else {
        /* As above, the x^(1-c) factor of [Abramowitz+Stegun, 15.1.2]
         * was missing (Savannah bug #53876). */
        gsl_sf_result F;
        gsl_sf_result powx;
        int stat_F = gsl_sf_hyperg_2F1_conj_e(aR-c+1, aI, -c+2, x, &F);
        int stat_p = gsl_sf_pow_int_e(x, 1 - (int)rintc, &powx);
        double ln_pre_val = 2.0*(g1.val - g2.val) - g3.val;
        double ln_pre_err = 2.0 * (g1.err + g2.err) + g3.err;
        double pre_val = powx.val * F.val;
        double pre_err = fabs(powx.val) * F.err
                       + fabs(F.val) * powx.err
                       + 2.0 * GSL_DBL_EPSILON * fabs(pre_val);
        int stat_e = gsl_sf_exp_mult_err_e(ln_pre_val, ln_pre_err,
                                              pre_val, pre_err,
                                              result);
        return GSL_ERROR_SELECT_3(stat_e, stat_F, stat_p);
      }
    }
  }
  else {
    /* generic c */
    gsl_sf_result F;
    gsl_sf_result lng;
    double sgn;
    int stat_g = gsl_sf_lngamma_sgn_e(c, &lng, &sgn);
    int stat_F = gsl_sf_hyperg_2F1_conj_e(aR, aI, c, x, &F);
    int stat_e = gsl_sf_exp_mult_err_e(-lng.val, lng.err,
                                          sgn*F.val, F.err,
                                          result);
    return GSL_ERROR_SELECT_3(stat_e, stat_F, stat_g);
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double gsl_sf_hyperg_2F1(double a, double b, double c, double x)
{
  EVAL_RESULT(gsl_sf_hyperg_2F1_e(a, b, c, x, &result));
}

double gsl_sf_hyperg_2F1_conj(double aR, double aI, double c, double x)
{
  EVAL_RESULT(gsl_sf_hyperg_2F1_conj_e(aR, aI, c, x, &result));
}

double gsl_sf_hyperg_2F1_renorm(double a, double b, double c, double x)
{
  EVAL_RESULT(gsl_sf_hyperg_2F1_renorm_e(a, b, c, x, &result));
}

double gsl_sf_hyperg_2F1_conj_renorm(double aR, double aI, double c, double x)
{
  EVAL_RESULT(gsl_sf_hyperg_2F1_conj_renorm_e(aR, aI, c, x, &result));
}
