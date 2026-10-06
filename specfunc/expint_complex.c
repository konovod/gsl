/* specfunc/expint_complex.c
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

/* Exponential integrals
 *
 *              E_n(z) = Integral_1^Infinity exp(-z t) / t^n dt
 *
 * for a complex argument z in the cut plane -Pi < arg(z) <= Pi, and the
 * associated integral
 *
 *              Ei(z) = - PV Integral_{-z}^Infinity exp(-t) / t dt.
 *
 * The implementation is a C translation of the double-precision routines
 * of
 *
 *   D. E. Amos, "Computation of Exponential Integrals of a Complex
 *   Argument", ACM Trans. Math. Software 16(2), 169-177 (1990),
 *   Algorithm 683, "A Portable FORTRAN Subroutine for Exponential
 *   Integrals of a Complex Argument", ibid. 178-182.
 *
 * It combines a power series for |z| <= 2, a backward (Miller)
 * recurrence for the confluent hypergeometric representation away from
 * the negative real axis, and Taylor continuation in a strip about the
 * negative real axis where the backward recurrence converges only
 * slowly.  The sequence length of the original routine is fixed at one
 * here: the public interface returns a single order.
 *
 * Conventions:
 *   - the branch cut of E_n runs along the negative real axis and is
 *     approached from above, i.e. -Pi < arg(z) <= Pi;
 *   - E_n_scaled(z) = exp(z) E_n(z);
 *   - a NaN argument propagates; a finite complex limit for zr = +Inf
 *     is returned for E_n (zero), and any other non-finite argument is
 *     a domain error.
 *
 * The error estimates are not produced by the algorithm itself; they
 * are estimated from the size of the computed result and are estimates
 * only.
 */

#include <config.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_psi.h>

#include "error.h"

/* storage for the coefficients of the backward recurrence; the original
 * routine allows for the recurrence to run past this without storing,
 * up to a hard limit for the termination test. */
#define EXPINT_ICDIM 250
#define EXPINT_ICMAX 2000

/*--------------------------------------------------------------------*/
/* complex helpers, matching ZMLT / ZDIV / ZSQRT / ZEXP / ZLOG / ZABS  */
/*--------------------------------------------------------------------*/

static double
zc_abs(double zr, double zi)
{
  return gsl_complex_abs(gsl_complex_rect(zr, zi));
}

static void
zc_div(double ar, double ai, double br, double bi, double *cr, double *ci)
{
  gsl_complex c = gsl_complex_div(gsl_complex_rect(ar, ai),
                                  gsl_complex_rect(br, bi));
  *cr = GSL_REAL(c);
  *ci = GSL_IMAG(c);
}

static void
zc_sqrt(double ar, double ai, double *br, double *bi)
{
  gsl_complex c = gsl_complex_sqrt(gsl_complex_rect(ar, ai));
  *br = GSL_REAL(c);
  *bi = GSL_IMAG(c);
}

static void
zc_exp(double ar, double ai, double *br, double *bi)
{
  gsl_complex c = gsl_complex_exp(gsl_complex_rect(ar, ai));
  *br = GSL_REAL(c);
  *bi = GSL_IMAG(c);
}

static void
zc_log(double ar, double ai, double *br, double *bi)
{
  /* principal branch, arg in (-Pi,Pi]; ar == ai == 0 is never reached
     by the callers below. */
  gsl_complex c = gsl_complex_log(gsl_complex_rect(ar, ai));
  *br = GSL_REAL(c);
  *bi = GSL_IMAG(c);
}

/*--------------------------------------------------------------------*/
/* machine parameters, matching the I1MACH / D1MACH block of ZEXINT   */
/*--------------------------------------------------------------------*/

static void
expint_machine_params(double *elim, double *alim, double *rbr_y, double *urnd)
{
  const double r1m5 = 0.301029995663981195;      /* log10(2) = D1MACH(5) */
  const int k = GSL_MIN(-DBL_MIN_EXP, DBL_MAX_EXP);

  *urnd = GSL_MAX(GSL_DBL_EPSILON, 1.0e-18);

  *elim = 2.303 * ((double)k * r1m5 - 3.0);
  {
    const double aa = r1m5 * (double)(DBL_MANT_DIG - 1) * 2.303;
    *alim = *elim + GSL_MAX(-aa, -41.45);
  }

  *rbr_y = (*urnd > 1.0e-8) ? 1.0 : 2.0;
}

/*--------------------------------------------------------------------*/
/* ZEXENZ: power series and backward recurrence                        */
/*--------------------------------------------------------------------*/

static void
expint_zexenz(int n, double zr, double zi, int scale,
              double rbr_y, double tol, double elim, double alim,
              int *ierr, double *cyr, double *cyi)
{
  const double az = zc_abs(zr, zi);

  if (az <= rbr_y) {
    /*------------------------------------------------------------------*/
    /* series for E(n,z) for |z| <= rbr_y                                */
    /*------------------------------------------------------------------*/
    const int iz  = (int)(az + 0.5);
    int icase = 2;
    int nm, nd, i, ik;
    double fnm, csr = 0.0, csi = 0.0, xtol, aam, caar, caai, aa, ak;

    if (iz > n)
      icase = 1;
    nm  = n - icase + 1;
    nd  = nm + 1;
    fnm = (double)nm;

    xtol = 0.3333 * tol;
    aam  = 1.0;
    if (nd != 1) {
      aam = 1.0 / fnm;
      csr = aam;
      csi = 0.0;
    }

    caar = 1.0;
    caai = 0.0;
    aa   = 1.0;
    ak   = 1.0;

    ik = (az < xtol * aam) ? 1 : 35;

    for (i = 1; i <= ik; i++) {
      double at = 1.0 / ak;
      double ap1;

      ap1  = -(caar * zr - caai * zi) * at;
      caai = -(caar * zi + caai * zr) * at;
      caar = ap1;
      aa   = aa * az * at;

      if (i == nm) {
        double catr, cati;
        zc_log(zr, zi, &catr, &cati);
        catr = -catr + gsl_sf_psi_int(nd);
        cati = -cati;
        csr += caar * catr - caai * cati;
        csi += caar * cati + caai * catr;
      }
      else {
        at   = 1.0 / (ak - fnm);
        csr -= caar * at;
        csi -= caai * at;
      }

      if (aa <= xtol * zc_abs(csr, csi))
        break;
      ak += 1.0;
    }

    if (nd == 1) {
      double catr, cati;
      zc_log(zr, zi, &catr, &cati);
      catr = -catr - M_EULER;
      cati = -cati;
      csr += catr;
      csi += cati;
    }

    if (scale) {
      double catr, cati, t;
      zc_exp(zr, zi, &catr, &cati);
      t   = csr * catr - csi * cati;
      csi = csr * cati + csi * catr;
      csr = t;
    }

    if (icase == 1) {
      double emzr = 1.0, emzi = 0.0;
      if (!scale)
        zc_exp(-zr, -zi, &emzr, &emzi);
      {
        double catr = emzr - csr;
        double cati = emzi - csi;
        zc_div(catr, cati, zr, zi, &csr, &csi);
      }
    }

    *cyr = csr;
    *cyi = csi;
    return;
  }

  /*--------------------------------------------------------------------*/
  /* backward recursive Miller algorithm for                            */
  /*         E(n,z) = exp(-z) z^(n-1) U(n,n,z)                          */
  /*--------------------------------------------------------------------*/
  {
    double scler = 1.0, r_scler = 1.0;
    double emzr = 1.0, emzi = 0.0;
    int iz, kn, icase, ks, jset, ah;
    double xaa, caar, caai, aam, aams, tzr, tzi, fzr, fzi, bk, ck, xtol;
    double cp1r, cp1i, cp2r, cp2i, ctr, cti;
    double cakr, caki, aem, at, bt, err, ap1, fc, ak;
    double cy1r, cy1i, cy2r, cy2i;
    double catr, cati, cbrr, cbri, cptr, cpti, cnrmr, cnrmi;
    double cyyr[3], cyyi[3];
    double ca[EXPINT_ICDIM + 1];
    double cbr[EXPINT_ICDIM + 1];
    double cbi[EXPINT_ICDIM + 1];
    int ic, ict, k;

    /* scale near exponent extremes on the unscaled variant */
    if (!scale) {
      int kflag;
      if (zr >= 0.0) {
        const double ct = zr + (double)n;
        const double aa2 = zc_abs(ct, zi);
        at = zr + log(aa2);
        if (at > elim) {
          *ierr = 2;
          *cyr = 0.0;
          *cyi = 0.0;
          return;
        }
        kflag = 1;
      }
      else {
        at = zr;
        if (at < -elim) {
          *ierr = 3;
          return;
        }
        kflag = 2;
      }

      if (fabs(at) >= alim) {
        const double tola = exp(alim - elim);
        const double r_tola = 1.0 / tola;
        if (kflag == 2) {
          scler   = tola;
          r_scler = r_tola;
        }
        else {
          scler   = r_tola;
          r_scler = tola;
        }
        emzr = scler;
        emzi = 0.0;
      }

      {
        double t;
        zc_exp(-zr, -zi, &catr, &cati);
        t    = emzr * catr - emzi * cati;
        emzi = emzr * cati + emzi * catr;
        emzr = t;
      }
    }

    iz = (int)(az + 0.5);
    kn = n;

    /* select the starting order and direction */
    if (kn <= iz) {
      icase = 1;
      ks = kn;
      if (kn > 1)
        goto start;
      ks = 2;
      icase = 3;
      goto start;
    }
    else if (n < iz && iz < kn) {
      *ierr = 7;
      return;
    }
    else if (n >= iz) {
      icase = 2;
      ks = n;
      if (n > 1)
        goto start;
      if (kn == 1) {
        ks = 2;
        icase = 3;
        goto start;
      }
      iz = 2;
    }
    else {
      *ierr = 7;
      return;
    }

  start:
    ah   = ks / 2;
    jset = 1 + ks - 2 * ah;
    ic   = 0;

    xaa  = (double)ah + (double)ah;
    caar = xaa;
    caai = 0.0;
    aam  = xaa - 1.0;
    aams = aam * aam;
    tzr  = zr + zr;
    tzi  = zi + zi;
    fzr  = tzr + tzr;
    fzi  = tzi + tzi;
    ak   = (double)ah;
    xtol = tol;
    ctr  = aams + fzr * (double)ah;
    cti  = fzi * (double)ah;
    cakr = zr + caar;
    caki = zi + caai;
    aem  = ((double)ah + 1.0) / xtol;
    aem /= zc_abs(cakr, caki);
    bk   = xaa;
    ck   = (double)ah * (double)ah;

    cp1r = 0.0; cp1i = 0.0;
    cp2r = 1.0; cp2i = 0.0;

    /* forward recursion for P(ic), P(ic+1) and the index ic for the
       backward recursion; the coefficients are stored while they fit */
    for (;;) {
      ic++;
      if (ic > EXPINT_ICDIM)
        break;

      ak += 1.0;
      ck += 1.0;
      at = bk / (bk + ak + ck);
      bk = bk + ak + ak;
      ca[ic]  = at;
      bt      = 1.0 / (ak + 1.0);
      cbr[ic] = (ak + ak + zr) * bt;
      cbi[ic] = zi * bt;

      cptr = cp2r;
      cpti = cp2i;
      ap1  = cbr[ic] * cp2r - cbi[ic] * cp2i - cp1r * at;
      cp2i = cbr[ic] * cp2i + cbi[ic] * cp2r - cp1i * at;
      cp2r = ap1;
      cp1r = cptr;
      cp1i = cpti;

      ctr += fzr;
      cti += fzi;
      aem *= at;

      bt  = zc_abs(ctr, cti);
      err = aem / sqrt(bt);
      ap1 = zc_abs(cp1r, cp1i);
      if (!(err * (ak + 1.0) / ap1 > ap1))
        break;
    }

    if (ic > EXPINT_ICDIM) {
      ic--;
      /* continue the forward recurrence unindexed */
      for (;;) {
        ic++;
        if (ic > EXPINT_ICMAX) {
          *ierr = 6;
          return;
        }

        ak += 1.0;
        ck += 1.0;
        at = bk / (bk + ak + ck);
        bk = bk + ak + ak;
        bt = 1.0 / (ak + 1.0);

        cptr = cp2r;
        cpti = cp2i;
        ap1  = (ak + ak + zr) * bt * cp2r - (zi * bt) * cp2i - cp1r * at;
        cp2i = (ak + ak + zr) * bt * cp2i + (zi * bt) * cp2r - cp1i * at;
        cp2r = ap1;
        cp1r = cptr;
        cp1i = cpti;

        ctr += fzr;
        cti += fzi;
        aem *= at;

        bt  = zc_abs(ctr, cti);
        err = aem / sqrt(bt);
        ap1 = zc_abs(cp1r, cp1i);
        if (!(err * (ak + 1.0) / ap1 > ap1))
          break;
      }
    }

    fc  = (double)ic;
    at  = ((fc + 1.0) / (ak + 1.0)) * ((ak + (double)ah) / (ak + 1.0));
    catr = ctr + fzr;
    cati = cti + fzi;
    zc_div(catr, cati, catr, cati, &catr, &cati);
    zc_sqrt(catr, cati, &catr, &cati);
    catr *= at;
    cati *= at;
    zc_div(cp1r, cp1i, cp2r, cp2i, &cbrr, &cbri);
    cy2r = catr * cbrr - cati * cbri;
    cy2i = catr * cbri + cati * cbrr;
    cy1r = 1.0;
    cy1i = 0.0;

    /* backward recurrence for
     *      cy1 = C U(a,a,z),  cy2 = C (a/(1+a/2)) U(a+1,a,z)
     */
    if (ic > EXPINT_ICDIM) {
      bt = xaa + zr;
      for (;;) {
        double bk2, ck2, dk2;
        ak  = (double)ah + fc;
        bk2 = ak + 1.0;
        ck2 = aam + fc;
        dk2 = bt + fc + fc;
        at  = (ak / fc) * (bk2 / ck2);
        cbrr = dk2 / bk2;
        cbri = zi / bk2;

        cptr = cy1r;
        cpti = cy1i;
        ap1  = (cbrr * cy1r - cbri * cy1i - cy2r) * at;
        cy1i = (cbrr * cy1i + cbri * cy1r - cy2i) * at;
        cy1r = ap1;
        cy2r = cptr;
        cy2i = cpti;

        fc -= 1.0;
        ic -= 1;
        if (ic <= EXPINT_ICDIM)
          break;
      }
    }

    ict = ic;
    for (k = 1; k <= ict; k++) {
      at   = 1.0 / ca[ic];
      cptr = cy1r;
      cpti = cy1i;
      ap1  = (cbr[ic] * cy1r - cbi[ic] * cy1i - cy2r) * at;
      cy1i = (cbr[ic] * cy1i + cbi[ic] * cy1r - cy2i) * at;
      cy1r = ap1;
      cy2r = cptr;
      cy2i = cpti;
      ic--;
    }

    /* contiguous relation
     *      z U(b,c+1,z) = (c-b) U(b,c,z) + U(b-1,c,z)
     * with b = a+1, c = a, giving C U(a+1,a+1,z)
     */
    zc_div(cy2r, cy2i, cy1r, cy1i, &cptr, &cpti);
    bt    = ((double)ah + 1.0) / xaa;
    cnrmr = 1.0 - cptr * bt;
    cnrmi = -cpti * bt;
    catr  = cnrmr * caar - cnrmi * caai + zr;
    cati  = cnrmr * caai + cnrmi * caar + zi;
    zc_div(1.0, 0.0, catr, cati, &cyyr[1], &cyyi[1]);
    cyyr[2] = cnrmr * cyyr[1] - cnrmi * cyyi[1];
    cyyi[2] = cnrmr * cyyi[1] + cnrmi * cyyr[1];

    if (icase == 3) {
      catr = 1.0 - cyyr[1];
      cati = -cyyi[1];
      {
        double cbtr = emzr * catr - emzi * cati;
        double cbti = emzr * cati + emzi * catr;
        zc_div(cbtr, cbti, zr, zi, &ctr, &cti);
      }
      *cyr = ctr * r_scler;
      *cyi = cti * r_scler;
      return;
    }

    ctr = emzr * cyyr[jset] - emzi * cyyi[jset];
    cti = emzr * cyyi[jset] + emzi * cyyr[jset];
    *cyr = ctr * r_scler;
    *cyi = cti * r_scler;
  }
}

/*--------------------------------------------------------------------*/
/* ZACEXI: analytic continuation in a strip about the negative axis    */
/*--------------------------------------------------------------------*/

static void
expint_zacexi(int nu, double zr, double zi, int scale, double yb,
              double rbr_y, double tol, double elim, double alim,
              int *ierr, double *yr, double *yi)
{
  double scler = 1.0, r_scler = 1.0;
  double zid = zi;
  int kyb = 0;
  const double az = zc_abs(zr, zi);
  const int iaz = (int)(az + 0.5);
  int nub, nb, nflg;

  if (zi < 0.0)
    zid = -zid;

  /*   nflg = 1: continue up from nub < nu
   *   nflg = 2: nub = nu, nothing to continue
   * (with a sequence length of one the branch nub = iaz < nu cannot
   * occur, so no downward continuation is needed) */
  if (nu >= iaz) {
    nub  = GSL_MAX(iaz, 1);
    nb   = nu - nub;
    nflg = 1;
  }
  else {
    nub  = nu;
    nb   = 0;
    nflg = 2;
  }

  for (;;) {
    double del, htol, xtol, ztr = 0.0, zti = 0.0;
    double sumr, sumi, trmr, trmi, ceztr, cezti, zwr, zwi;
    double zpr[65], zpi[65], rq[65];
    double fk, fj;
    double cexr, cexi, yyr = 0.0, yyi = 0.0;
    int iy, k, kmax, i;

    del = yb - zid;

    /* make del large enough to avoid underflow in the powers */
    if (fabs(del) <= 1.0e-4) {
      yb += 1.0e-4;
      del = yb - zid;
    }

    htol = 0.125 * tol;

    expint_zexenz(nub, zr, yb, 1, rbr_y, htol, elim, alim,
                  ierr, &yyr, &yyi);
    if (*ierr == 6) {
      yb += 0.5;
      kyb++;
      if (kyb > 10)
        return;
      *ierr = 0;
      continue;
    }

    /* analytic continuation by Taylor series for order nub */
    iy   = (int)(del + del);
    {
      const double yt = del / (double)(iy + 1);
      sumr = yyr;
      sumi = yyi;
      htol = 0.25 * tol;
      trmr = sumr;
      trmi = sumi;
      ceztr = cos(yt);
      cezti = -sin(yt);
      zwr = 1.0;
      zwi = 0.0;
      zpr[1] = 1.0;
      zpi[1] = 0.0;
      fk = 1.0;
      fj = (double)(nub - 1);
      zc_div(1.0, 0.0, zr, yb, &ztr, &zti);

      for (k = 2; k <= 64; k++) {
        double rw, atrm, asum, ctr, cti;

        rw   = yt / fk;
        atrm = -zwi * rw;
        zwi  = zwr * rw;
        zwr  = atrm;
        rw   = fj * rw;

        ctr  = zwr + trmi * rw;
        cti  = zwi - trmr * rw;
        trmr = ctr * ztr - cti * zti;
        trmi = ctr * zti + cti * ztr;
        sumr += trmr;
        sumi += trmi;

        zpr[k] = zwr;
        zpi[k] = zwi;
        rq[k]  = rw;

        asum = zc_abs(sumr, sumi);
        atrm = zc_abs(trmr, trmi);
        if (atrm < htol * asum)
          break;

        fk += 1.0;
        fj -= 1.0;
      }
      if (k > 64)
        k = 64;
      kmax = k;

      {
        double atrm = sumr * ceztr - sumi * cezti;
        sumi = sumr * cezti + sumi * ceztr;
        sumr = atrm;
      }

      for (i = 1; i <= iy; i++) {
        const double rzi = (double)(iy - i + 1) * yt + zid;
        double atrm, asum;

        zc_div(1.0, 0.0, zr, rzi, &ztr, &zti);
        trmr = sumr;
        trmi = sumi;

        for (k = 2; k <= kmax; k++) {
          double ctr = zpr[k] + trmi * rq[k];
          double cti = zpi[k] - trmr * rq[k];
          trmr = ctr * ztr - cti * zti;
          trmi = ctr * zti + cti * ztr;
          sumr += trmr;
          sumi += trmi;
        }

        atrm = zc_abs(trmr, trmi);
        asum = zc_abs(sumr, sumi);
        xtol = htol * asum;
        if (atrm < xtol)
          goto next_i;
        if (kmax >= 64)
          goto next_i;

        kmax++;
        for (k = kmax; k <= 64; k++) {
          double rw, ctr, cti;
          rw   = yt / fk;
          atrm = -zwi * rw;
          zwi  = zwr * rw;
          zwr  = atrm;
          rw   = fj * rw;

          ctr  = zwr + trmi * rw;
          cti  = zwi - trmr * rw;
          trmr = ctr * ztr - cti * zti;
          trmi = ctr * zti + cti * ztr;
          sumr += trmr;
          sumi += trmi;

          zpr[k] = zwr;
          zpi[k] = zwi;
          rq[k]  = rw;

          atrm = zc_abs(trmr, trmi);
          if (atrm < xtol)
            break;

          fk += 1.0;
          fj -= 1.0;
        }
        if (k > 64)
          k = 64;
        kmax = k;

      next_i:
        {
          double atrm = sumr * ceztr - sumi * cezti;
          sumi = sumr * cezti + sumi * ceztr;
          sumr = atrm;
        }
      }

      if (zi < 0.0)
        sumi = -sumi;

      cexr = 1.0;
      cexi = 0.0;

      /* scale near the overflow limit on the unscaled variant */
      if (!scale) {
        if (fabs(zr) >= alim) {
          if (fabs(zr) > elim) {
            *ierr = 3;
            return;
          }
          {
            const double tola = exp(alim - elim);
            const double r_tola = 1.0 / tola;
            scler   = tola;
            r_scler = r_tola;
          }
        }
        {
          double ctr, cti;
          zc_exp(-zr, -zi, &ctr, &cti);
          cexr = scler * ctr;
          cexi = scler * cti;
        }
      }

      {
        const double trm_r = sumr * cexr - sumi * cexi;
        const double trm_i = sumr * cexi + sumi * cexr;

        if (nflg != 1) {
          *yr = trm_r * r_scler;
          *yi = trm_i * r_scler;
          return;
        }

        /* continue upward from nub to nu */
        {
          double trsr = trm_r, trsi = trm_i;
          double fk2 = (double)nub;

          for (k = 1; k <= nb; k++) {
            const double at = 1.0 / fk2;
            const double atrm = (cexr - (zr * trsr - zi * trsi)) * at;
            trsi = (cexi - (zr * trsi + zi * trsr)) * at;
            trsr = atrm;
            fk2 += 1.0;
          }

          *yr = trsr * r_scler;
          *yi = trsi * r_scler;
          return;
        }
      }
    }
  }
}

/*--------------------------------------------------------------------*/
/* ZEXINT: dispatch                                                    */
/*--------------------------------------------------------------------*/

static void
expint_zexint(int n, double zr, double zi, int scale,
              int *ierr, double *cyr, double *cyi)
{
  const double tol = GSL_DBL_EPSILON;
  double elim, alim, rbr_y, urnd;
  double az, fn, xaa;

  *ierr = 0;
  if (zr == 0.0 && zi == 0.0 && n == 1)
    *ierr = 1;
  if (n < 1)
    *ierr = 1;
  if (*ierr != 0)
    return;

  expint_machine_params(&elim, &alim, &rbr_y, &urnd);

  if (zr == 0.0 && zi == 0.0 && n > 1) {
    *cyr = 1.0 / (double)(n - 1);
    *cyi = 0.0;
    return;
  }

  az = zc_abs(zr, zi);
  fn = (double)n;

  xaa = 0.5 / urnd;
  xaa = GSL_MIN(xaa, 0.5 * (double)INT_MAX);
  if (az > xaa) {
    *ierr = 5;
    return;
  }
  if (fn > xaa) {
    *ierr = 5;
    return;
  }
  xaa = sqrt(xaa);
  if (az > xaa)
    *ierr = 4;
  if (fn > xaa)
    *ierr = 4;

  if (zr >= 0.0) {
    expint_zexenz(n, zr, zi, scale, rbr_y, tol, elim, alim, ierr, cyr, cyi);
    return;
  }

  if (az <= rbr_y) {
    expint_zexenz(n, zr, zi, scale, rbr_y, tol, elim, alim, ierr, cyr, cyi);
    return;
  }

  {
    const double d  = -0.4342945 * log10(tol);
    const double yb = 10.5 - 0.538460 * (18.0 - d);

    if (fabs(zi) < yb) {
      expint_zacexi(n, zr, zi, scale, yb, rbr_y, tol, elim, alim,
                    ierr, cyr, cyi);
    }
    else {
      const double htol = 0.125 * tol;
      expint_zexenz(n, zr, zi, scale, rbr_y, htol, elim, alim,
                    ierr, cyr, cyi);
    }
  }
}

/*--------------------------------------------------------------------*/
/* error estimates                                                     */
/*--------------------------------------------------------------------*/

/* The algorithm carries no error estimate.  The dominant error is the
 * truncation controlled by its convergence test, which is relative to
 * the size of the result; a handful of ulps of the combined magnitude
 * is a fair, if conservative, bar.  Underflow/overflow are reported
 * separately. */
#define EXPINT_ERR_FACTOR 64.0

static void
expint_set_err(gsl_sf_result *re, gsl_sf_result *im)
{
  double err = EXPINT_ERR_FACTOR * GSL_DBL_EPSILON
             * (fabs(re->val) + fabs(im->val));
  if (err == 0.0)
    err = GSL_DBL_MIN;
  re->err = err;
  im->err = err;
}

/*--------------------------------------------------------------------*/
/* entry point for E_n                                                 */
/*--------------------------------------------------------------------*/

static int
expint_complex_En_impl(int n, double zr, double zi, int scale,
                       gsl_sf_result *re, gsl_sf_result *im)
{
  double cyr = 0.0, cyi = 0.0;
  int ierr = 0;

  if (gsl_isnan(zr) || gsl_isnan(zi)) {
    re->val = GSL_NAN; re->err = GSL_NAN;
    im->val = GSL_NAN; im->err = GSL_NAN;
    return GSL_SUCCESS;
  }

  if (zi == 0.0)
    zi = 0.0;                    /* normalise -0 to the upper side */

  if (gsl_isinf(zr) || gsl_isinf(zi)) {
    if (gsl_isinf(zr) && zr > 0.0 && gsl_finite(zi)) {
      re->val = 0.0; re->err = 0.0;
      im->val = 0.0; im->err = 0.0;
      return GSL_SUCCESS;
    }
    DOMAIN_ERROR_2(re, im);
  }

  if (zr == 0.0 && zi == 0.0) {
    if (n == 1)
      DOMAIN_ERROR_2(re, im);
    re->val = 1.0 / (double)(n - 1);
    re->err = 0.0;
    im->val = 0.0;
    im->err = 0.0;
    return GSL_SUCCESS;
  }

  expint_zexint(n, zr, zi, scale, &ierr, &cyr, &cyi);

  switch (ierr) {
  case 0:
    re->val = cyr;
    im->val = cyi;
    expint_set_err(re, im);
    return GSL_SUCCESS;

  case 4:
    re->val = cyr;
    im->val = cyi;
    expint_set_err(re, im);
    return GSL_ELOSS;

  case 2:
    re->val = 0.0;  re->err = GSL_DBL_MIN;
    im->val = 0.0;  im->err = GSL_DBL_MIN;
    return GSL_EUNDRFLW;

  case 3:
    OVERFLOW_ERROR_2(re, im);

  case 6:
    re->val = GSL_NAN; re->err = GSL_NAN;
    im->val = GSL_NAN; im->err = GSL_NAN;
    GSL_ERROR("iteration limit exceeded", GSL_EMAXITER);

  case 5:
    re->val = GSL_NAN; re->err = GSL_NAN;
    im->val = GSL_NAN; im->err = GSL_NAN;
    return GSL_ELOSS;

  case 7:
  default:
    re->val = GSL_NAN; re->err = GSL_NAN;
    im->val = GSL_NAN; im->err = GSL_NAN;
    GSL_ERROR("algorithm failure", GSL_EFAILED);
  }
}

/*--------------------------------------------------------------------*/
/* Functions with error codes                                          */
/*--------------------------------------------------------------------*/

int
gsl_sf_complex_expint_En_e(int n, double zr, double zi,
                           gsl_sf_result * re, gsl_sf_result * im)
{
  if (n < 1)
    DOMAIN_ERROR_2(re, im);
  return expint_complex_En_impl(n, zr, zi, 0, re, im);
}

int
gsl_sf_complex_expint_En_scaled_e(int n, double zr, double zi,
                                  gsl_sf_result * re, gsl_sf_result * im)
{
  if (n < 1)
    DOMAIN_ERROR_2(re, im);
  return expint_complex_En_impl(n, zr, zi, 1, re, im);
}

int
gsl_sf_complex_expint_E1_e(double zr, double zi,
                           gsl_sf_result * re, gsl_sf_result * im)
{
  return expint_complex_En_impl(1, zr, zi, 0, re, im);
}

int
gsl_sf_complex_expint_E1_scaled_e(double zr, double zi,
                                  gsl_sf_result * re, gsl_sf_result * im)
{
  return expint_complex_En_impl(1, zr, zi, 1, re, im);
}

/* Ei(z) = -E1(-z) + (log z - log(-z)), and log z - log(-z) is +i Pi on
 * the upper side of the cut and -i Pi on the lower side (including the
 * positive real axis).  This keeps Ei real on the positive real axis
 * and gives the upper-side limit on the negative real axis. */
int
gsl_sf_complex_expint_Ei_e(double zr, double zi,
                           gsl_sf_result * re, gsl_sf_result * im)
{
  gsl_sf_result e1r, e1i;
  double theta, d;
  int status;

  if (gsl_isnan(zr) || gsl_isnan(zi)) {
    re->val = GSL_NAN; re->err = GSL_NAN;
    im->val = GSL_NAN; im->err = GSL_NAN;
    return GSL_SUCCESS;
  }

  if (zi == 0.0)
    zi = 0.0;

  if (zr == 0.0 && zi == 0.0)
    DOMAIN_ERROR_2(re, im);

  /* Ei has no finite limit as z -> infinity in any direction */
  if (gsl_isinf(zr) || gsl_isinf(zi))
    DOMAIN_ERROR_2(re, im);

  theta = atan2(zi, zr);
  d     = (theta > 0.0) ? M_PI : -M_PI;

  status = expint_complex_En_impl(1, -zr, -zi, 0, &e1r, &e1i);

  re->val = -e1r.val;
  im->val = -e1i.val + d;
  re->err = e1r.err + 2.0 * GSL_DBL_EPSILON * fabs(re->val);
  im->err = e1i.err + 2.0 * GSL_DBL_EPSILON * (fabs(im->val) + M_PI);

  return status;
}
