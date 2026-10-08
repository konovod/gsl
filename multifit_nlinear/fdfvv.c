/* multifit_nlinear/fdfvv.c
 * 
 * Copyright (C) 2015 Patrick Alken
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

/*
fdfvv()
  Compute approximate second directional derivative using
finite differences.

See Eq. 19 of:

M. K. Transtrum, J. P. Sethna, Improvements to the Levenberg
Marquardt algorithm for nonlinear least-squares minimization,
arXiv:1201.5885, 2012.

Inputs: h     - step size for finite difference
        x     - parameter vector, size p
        v     - geodesic velocity, size p
        f     - vector of function values f_i(x), size n
        J     - Jacobian matrix J(x), n-by-p
        swts  - data weights
        fdf   - fdf struct
        fvv   - (output) approximate second directional derivative
                vector D_v^2 f(x)
        work  - workspace, size p

Return: success or error
*/

static int
fdfvv(const double h, const gsl_vector *x, const gsl_vector *v,
      const gsl_vector *f, const gsl_matrix *J, const gsl_vector *swts,
      gsl_multifit_nlinear_fdf *fdf, gsl_vector *fvv, gsl_vector *work)
{
  int status;
  const size_t n = fdf->n;
  const size_t p = fdf->p;
  const double vnorm = gsl_blas_dnrm2(v);
  double hstep = h;
  double hinv;
  size_t i;

  /*
   * The step h scales the velocity v, so the perturbation is x + h v.
   * When ||v|| is very small (near convergence, or when the LM parameter
   * is large) that perturbation can fall below the rounding unit of x;
   * f(x + h v) is then evaluated at a point indistinguishable from x and
   * the resulting fvv is dominated by roundoff rather than by the second
   * directional derivative.  Ensure the perturbation retains a meaningful
   * size relative to the parameters, as the finite difference Jacobian
   * does for each parameter in fdjac.c.
   *
   * If v is exactly zero the acceleration contributes nothing, so fvv = 0.
   */
  if (vnorm == 0.0)
    {
      gsl_vector_set_zero(fvv);
      return GSL_SUCCESS;
    }
  else
    {
      const double xnorm = gsl_blas_dnrm2(x);
      const double hmin = GSL_SQRT_DBL_EPSILON * (1.0 + xnorm) / vnorm;

      /* hmin can overflow for a denormal ||v||; leave h unchanged then,
       * since such a velocity is numerically zero anyway */
      if (hstep < hmin && gsl_finite(hmin))
        hstep = hmin;
    }

  hinv = 1.0 / hstep;

  /* compute work = x + hstep*v */
  for (i = 0; i < p; ++i)
    {
      double xi = gsl_vector_get(x, i);
      double vi = gsl_vector_get(v, i);

      gsl_vector_set(work, i, xi + hstep * vi);
    }

  /* compute f(x + hstep*v) */
  status = gsl_multifit_nlinear_eval_f (fdf, work, swts, fvv);
  if (status)
    return status;

  for (i = 0; i < n; ++i)
    {
      double fi = gsl_vector_get(f, i);    /* f_i(x) */
      double fip = gsl_vector_get(fvv, i); /* f_i(x + hstep*v) */
      gsl_vector_const_view row = gsl_matrix_const_row(J, i);
      double u, fvvi;

      /* compute u = sum_j J_{ij} v_j */
      gsl_blas_ddot(&row.vector, v, &u);

      fvvi = (2.0 * hinv) * ((fip - fi) * hinv - u);

      gsl_vector_set(fvv, i, fvvi);
    }

  return status;
}

/*
gsl_multifit_nlinear_fdfvv()
  Compute approximate second directional derivative
using finite differences

Inputs: h    - step size for finite difference
        x    - parameter vector, size p
        v    - geodesic velocity, size p
        f    - function values f_i(x), size n
        J    - Jacobian matrix J(x), n-by-p
        swts - sqrt data weights (set to NULL if not needed)
        fdf  - fdf
        fvv  - (output) approximate (weighted) second directional derivative
               vector, size n, sqrt(W) fvv
        work - workspace, size p

Return: success or error
*/

int
gsl_multifit_nlinear_fdfvv(const double h, const gsl_vector *x, const gsl_vector *v,
                           const gsl_vector *f, const gsl_matrix *J,
                           const gsl_vector *swts, gsl_multifit_nlinear_fdf *fdf,
                           gsl_vector *fvv, gsl_vector *work)
{
  return fdfvv(h, x, v, f, J, swts, fdf, fvv, work);
}
