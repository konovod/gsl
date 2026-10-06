/* linalg/givensc.c
 * 
 * Copyright (C) 2024 Christian Krueger
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
#include <gsl/gsl_complex_math.h>

#include <gsl/gsl_linalg.h>

/* Generate a Givens rotation which takes v=(a,b) to (|v|,0)
   From Bindel et al. On computing Givens rotations reliably
   and efficiently.
   https://www.netlib.org/lapack/lawnspdf/lawn148.pdf
   Not all special cases that improve accuracy are implemented. */
void
gsl_linalg_complex_givens (const gsl_complex a, const gsl_complex b,
                           double *c, gsl_complex *s)
{
  if (GSL_REAL(b) == 0 && GSL_IMAG(b) == 0)
    {
      *c = 1;
      *s = GSL_COMPLEX_ZERO;
    }
  else if (GSL_REAL(a) == 0 && GSL_IMAG(a) == 0)
    {
      *c = 0;
      double d1 = gsl_complex_abs(b);
      *s = gsl_complex_conjugate(b);
      *s = gsl_complex_div_real(*s, d1);
    }
  else
    {
      double a2 = gsl_complex_abs2(a);
      double b2 = gsl_complex_abs2(b);
      double ab2 = a2 + b2;
      double d1 = 1.0 / sqrt(a2*ab2);
      *c = a2 * d1;
      ab2 *= d1;
      *s = gsl_complex_mul_real(a, d1);
      *s = gsl_complex_mul(gsl_complex_conjugate(b), *s);
    }
} /* gsl_linalg_complex_givens() */

void gsl_linalg_complex_givens_gv (gsl_vector_complex * v,
                                   const size_t i, const size_t j,
                                   const double c, const gsl_complex s)
{
  gsl_complex h1, h2;

  gsl_complex vi = gsl_vector_complex_get (v, i);
  gsl_complex vj = gsl_vector_complex_get (v, j);

  h1 = gsl_complex_mul_real(vi, c);
  h2 = gsl_complex_mul(vj, s);
  h1 = gsl_complex_add(h1, h2);
  gsl_vector_complex_set (v, i, h1);

  h1 = gsl_complex_mul_real(vj, c);
  h2 = gsl_complex_mul(vi, gsl_complex_conjugate(s));
  h1 = gsl_complex_sub(h1, h2);
  gsl_vector_complex_set (v, j, h1);
  // TODO Why not set this component to zero as mathematically
  // that is the case?
}
