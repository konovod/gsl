/* specfunc/gsl_sf_expint.h
 * 
 * Copyright (C) 2007 Brian Gough
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

/* Author: G. Jungman */

#ifndef GSL_SF_EXPINT_H__
#define GSL_SF_EXPINT_H__

#include <gsl/gsl_sf_result.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS


/* E_1(x) := Re[ Integrate[ Exp[-xt]/t, {t,1,Infinity}] ]
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_E1_e(const double x, gsl_sf_result * result);
double  gsl_sf_expint_E1(const double x);


/* E_2(x) := Re[ Integrate[ Exp[-xt]/t^2, {t,1,Infinity}] ]
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_E2_e(const double x, gsl_sf_result * result);
double  gsl_sf_expint_E2(const double x);


/* E_n(x) := Re[ Integrate[ Exp[-xt]/t^n, {t,1,Infinity}] ]
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_En_e(const int n, const double x, gsl_sf_result * result);
double  gsl_sf_expint_En(const int n, const double x);


/* E_1_scaled(x) := exp(x) E_1(x)
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_E1_scaled_e(const double x, gsl_sf_result * result);
double  gsl_sf_expint_E1_scaled(const double x);


/* E_2_scaled(x) := exp(x) E_2(x)
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_E2_scaled_e(const double x, gsl_sf_result * result);
double  gsl_sf_expint_E2_scaled(const double x);

/* E_n_scaled(x) := exp(x) E_n(x)
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_En_scaled_e(const int n, const double x, gsl_sf_result * result);
double  gsl_sf_expint_En_scaled(const int n, const double x);


/* Ei(x) := - PV Integrate[ Exp[-t]/t, {t,-x,Infinity}]
 *       :=   PV Integrate[ Exp[t]/t, {t,-Infinity,x}]
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_Ei_e(const double x, gsl_sf_result * result);
double  gsl_sf_expint_Ei(const double x);


/* Ei_scaled(x) := exp(-x) Ei(x)
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_expint_Ei_scaled_e(const double x, gsl_sf_result * result);
double  gsl_sf_expint_Ei_scaled(const double x);


/* Shi(x) := Integrate[ Sinh[t]/t, {t,0,x}]
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_Shi_e(const double x, gsl_sf_result * result);
double  gsl_sf_Shi(const double x);


/* Chi(x) := Re[ M_EULER + log(x) + Integrate[(Cosh[t]-1)/t, {t,0,x}] ]
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_Chi_e(const double x, gsl_sf_result * result);
double  gsl_sf_Chi(const double x);


/* Ei_3(x) := Integral[ Exp[-t^3], {t,0,x}]
 *
 * x >= 0.0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_expint_3_e(const double x, gsl_sf_result * result);
double  gsl_sf_expint_3(double x);


/* Si(x) := Integrate[ Sin[t]/t, {t,0,x}]
 *
 * exceptions: none
 */
int     gsl_sf_Si_e(const double x, gsl_sf_result * result);
double  gsl_sf_Si(const double x);


/* Ci(x) := -Integrate[ Cos[t]/t, {t,x,Infinity}]
 *
 * x > 0.0
 * exceptions: GSL_EDOM 
 */
int     gsl_sf_Ci_e(const double x, gsl_sf_result * result);
double  gsl_sf_Ci(const double x);


/* AtanInt(x) := Integral[ Arctan[t]/t, {t,0,x}]
 *
 *
 * exceptions:
 */
int     gsl_sf_atanint_e(const double x, gsl_sf_result * result);
double  gsl_sf_atanint(const double x);


/* Complex exponential integrals
 *
 * E_n(z) := Integral_1^Infinity exp(-z t) / t^n dt
 *
 * analytically continued to the cut plane -Pi < arg(z) <= Pi; the
 * negative real axis is approached from above.  z != 0 when n == 1.
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW, GSL_ELOSS,
 *             GSL_EMAXITER, GSL_EFAILED
 */
int     gsl_sf_complex_expint_E1_e(double zr, double zi,
                                   gsl_sf_result * re, gsl_sf_result * im);
int     gsl_sf_complex_expint_En_e(int n, double zr, double zi,
                                   gsl_sf_result * re, gsl_sf_result * im);


/* E_n_scaled(z) := exp(z) E_n(z)
 *
 * exceptions: as for gsl_sf_complex_expint_E1_e / _En_e
 */
int     gsl_sf_complex_expint_E1_scaled_e(double zr, double zi,
                                          gsl_sf_result * re, gsl_sf_result * im);
int     gsl_sf_complex_expint_En_scaled_e(int n, double zr, double zi,
                                          gsl_sf_result * re, gsl_sf_result * im);


/* Ei(z) := - PV Integral_{-z}^Infinity exp(-t) / t dt
 *
 * for z in the cut plane -Pi < arg(z) <= Pi.  z != 0.
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW, GSL_ELOSS,
 *             GSL_EMAXITER, GSL_EFAILED
 */
int     gsl_sf_complex_expint_Ei_e(double zr, double zi,
                                   gsl_sf_result * re, gsl_sf_result * im);


__END_DECLS

#endif /* GSL_SF_EXPINT_H__ */
