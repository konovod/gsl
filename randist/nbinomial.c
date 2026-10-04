/* randist/nbinomial.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 James Theiler, Brian Gough
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
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sys.h>

/* The negative binomial distribution has the form,

   prob(k) =  Gamma(n + k)/(Gamma(n) Gamma(k + 1))  p^n (1-p)^k 

   for k = 0, 1, ... . Note that n does not have to be an integer.

   This is the Leger's algorithm (given in the answers in Knuth) */

unsigned int
gsl_ran_negative_binomial (const gsl_rng * r, double p, double n)
{
  double X = gsl_ran_gamma (r, n, 1.0) ;
  unsigned int k = gsl_ran_poisson (r, X*(1-p)/p) ;
  return k ;
}

double
gsl_ran_negative_binomial_pdf (const unsigned int k, const double p, double n)
{
  if (p == 1.0)
    {
      /* The distribution degenerates to k == 0: p^n (1-p)^k -> 1 for
         k = 0 and 0 for k > 0.  Taken before the general formula, where
         n * log(p) = 0 and k * log1p(-p) = -inf give NaN at k = 0.
         Savannah bug #66800. */
      return (k == 0) ? 1.0 : 0.0;
    }
  else if (p == 0.0)
    {
      /* Here (1-p)^k = 1 and p^n = 0, so everything collapses to k == 0
         as well (the 0^0 weight at k = 0 is taken as 1). */
      return (k == 0) ? 1.0 : 0.0;
    }
  else
    {
      double f = gsl_sf_lngamma (k + n) ;
      double a = gsl_sf_lngamma (n) ;
      double b = gsl_sf_lngamma (k + 1.0) ;

      return exp(f - a - b + n * log(p) + k * log1p(-p));
    }
}
