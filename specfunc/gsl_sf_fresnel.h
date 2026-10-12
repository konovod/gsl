/* specfunc/gsl_sf_fresnel.h
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

#ifndef GSL_SF_FRESNEL_H__
#define GSL_SF_FRESNEL_H__

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

/* Fresnel cosine integral
 *
 *   C(x) = Integral_0^x cos(pi t^2 / 2) dt
 *
 * C(-x) = -C(x).
 * exceptions: GSL_EDOM (x is infinite)
 */
int gsl_sf_fresnel_c_e(double x, gsl_sf_result * result);
double gsl_sf_fresnel_c(double x);

/* Fresnel sine integral
 *
 *   S(x) = Integral_0^x sin(pi t^2 / 2) dt
 *
 * S(-x) = -S(x).
 * exceptions: GSL_EDOM (x is infinite)
 */
int gsl_sf_fresnel_s_e(double x, gsl_sf_result * result);
double gsl_sf_fresnel_s(double x);

__END_DECLS

#endif /* GSL_SF_FRESNEL_H__ */
