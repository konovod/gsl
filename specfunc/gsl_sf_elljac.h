/* specfunc/gsl_sf_elljac.h
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

#ifndef GSL_SF_ELLJAC_H__
#define GSL_SF_ELLJAC_H__

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


/* Jacobian elliptic functions sn, dn, cn,
 * by descending Landen transformations
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_elljac_e(double u, double m, double * sn, double * cn, double * dn);

/* Inverse Jacobi elliptic functions.
 *
 * For the parameter m with |m| <= 1 these return the principal value
 * u such that the corresponding direct function equals the argument,
 * i.e. gsl_sf_elljac_e(u, m, &sn, &cn, &dn) reproduces the argument.
 * The principal value satisfies 0 <= u <= 2 K(m).  The inverse of dn
 * is not defined for m = 0, where dn(u|0) = 1 for every u.
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_elljac_arcsn_e(double sn, double m, gsl_sf_result * result);
double gsl_sf_elljac_arcsn(double sn, double m);

int gsl_sf_elljac_arccn_e(double cn, double m, gsl_sf_result * result);
double gsl_sf_elljac_arccn(double cn, double m);

int gsl_sf_elljac_arcdn_e(double dn, double m, gsl_sf_result * result);
double gsl_sf_elljac_arcdn(double dn, double m);


__END_DECLS

#endif /* GSL_SF_ELLJAC_H__ */
