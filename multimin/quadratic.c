/* multimin/quadratic.c
 * 
 * Copyright (C) 2011 Jussi Lehtola
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

/* A "brute force" multidimensional minimiser which needs no derivatives.
   It minimises along each coordinate direction in turn: for the current
   coordinate it fits a parabola through the current point and the two
   points one step away in the positive and negative directions, and
   moves to the predicted minimum if that improves the value; otherwise
   the step is halved and the fit retried.  On a quadratic the
   parabolic step is exact, so the method converges in a few sweeps.

   Suggested by Jussi Lehtola, see Savannah bug #32776. */

#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>

typedef struct
{
  gsl_vector *step_size;        /* current step size in each coordinate */
  gsl_vector *xl;               /* trial point x - step */
  gsl_vector *xr;               /* trial point x + step */
  gsl_vector *xtrial;           /* trial point from the parabolic fit */
}
quadratic_state_t;


static int
quadratic_alloc (void *vstate, size_t n)
{
  quadratic_state_t *state = (quadratic_state_t *) vstate;

  if (n == 0)
    {
      GSL_ERROR ("invalid number of parameters specified", GSL_EINVAL);
    }

  state->step_size = gsl_vector_alloc (n);

  if (state->step_size == NULL)
    {
      GSL_ERROR ("failed to allocate space for step_size", GSL_ENOMEM);
    }

  state->xl = gsl_vector_alloc (n);

  if (state->xl == NULL)
    {
      gsl_vector_free (state->step_size);
      GSL_ERROR ("failed to allocate space for xl", GSL_ENOMEM);
    }

  state->xr = gsl_vector_alloc (n);

  if (state->xr == NULL)
    {
      gsl_vector_free (state->xl);
      gsl_vector_free (state->step_size);
      GSL_ERROR ("failed to allocate space for xr", GSL_ENOMEM);
    }

  state->xtrial = gsl_vector_alloc (n);

  if (state->xtrial == NULL)
    {
      gsl_vector_free (state->xr);
      gsl_vector_free (state->xl);
      gsl_vector_free (state->step_size);
      GSL_ERROR ("failed to allocate space for xtrial", GSL_ENOMEM);
    }

  return GSL_SUCCESS;
}


static double
quadratic_size (const quadratic_state_t * state)
{
  /* the size of the minimiser is the norm of the step sizes divided by
     the number of parameters */
  return gsl_blas_dnrm2 (state->step_size) / (double) state->step_size->size;
}


static int
quadratic_set (void *vstate, gsl_multimin_function * f,
               const gsl_vector * x, double * size,
               const gsl_vector * step_size)
{
  quadratic_state_t *state = (quadratic_state_t *) vstate;
  double val;

  if (state->step_size->size != x->size)
    {
      GSL_ERROR ("incompatible size of x", GSL_EINVAL);
    }

  if (state->step_size->size != step_size->size)
    {
      GSL_ERROR ("incompatible size of step_size", GSL_EINVAL);
    }

  val = GSL_MULTIMIN_FN_EVAL (f, x);

  if (!gsl_finite (val))
    {
      GSL_ERROR ("non-finite function value encountered", GSL_EBADFUNC);
    }

  gsl_vector_memcpy (state->step_size, step_size);

  *size = quadratic_size (state);

  return GSL_SUCCESS;
}


static int
quadratic_iterate (void *vstate, gsl_multimin_function * f,
                   gsl_vector * x, double * size, double * fval)
{
  quadratic_state_t *state = (quadratic_state_t *) vstate;
  gsl_vector *xl = state->xl;
  gsl_vector *xr = state->xr;
  gsl_vector *xtrial = state->xtrial;
  const size_t n = x->size;
  size_t i;
  double fm;

  fm = GSL_MULTIMIN_FN_EVAL (f, x);

  if (!gsl_finite (fm))
    {
      GSL_ERROR ("non-finite function value encountered", GSL_EBADFUNC);
    }

  /* Increase the step sizes: after a successful parabolic step the
     step can be reduced again, and this allows the minimiser to follow
     the minimum as it moves. */
  gsl_vector_scale (state->step_size, 4.0);

  for (i = 0; i < n; i++)
    {
      while (1)
        {
          const double pval = gsl_vector_get (x, i);
          const double ss = gsl_vector_get (state->step_size, i);
          double fl, fr, ftrial, delta, curv;

          gsl_vector_memcpy (xl, x);
          gsl_vector_set (xl, i, pval - ss);

          gsl_vector_memcpy (xr, x);
          gsl_vector_set (xr, i, pval + ss);

          fl = GSL_MULTIMIN_FN_EVAL (f, xl);
          fr = GSL_MULTIMIN_FN_EVAL (f, xr);

          if (!gsl_finite (fl) || !gsl_finite (fr))
            {
              GSL_ERROR ("non-finite function value encountered", GSL_EBADFUNC);
            }

          /* curvature of the parabola fitted through (pval-ss, fl),
             (pval, fm) and (pval+ss, fr) */
          curv = fr - 2.0 * fm + fl;

          if (fabs (curv) > DBL_EPSILON)
            {
              /* position of the fitted minimum relative to pval */
              delta = -0.5 * ss * (fr - fl) / curv;

              if (delta > ss)
                {
                  gsl_vector_memcpy (xtrial, xr);
                  ftrial = fr;
                }
              else if (delta < -ss)
                {
                  gsl_vector_memcpy (xtrial, xl);
                  ftrial = fl;
                }
              else
                {
                  gsl_vector_memcpy (xtrial, x);
                  gsl_vector_set (xtrial, i, pval + delta);
                  ftrial = GSL_MULTIMIN_FN_EVAL (f, xtrial);

                  if (!gsl_finite (ftrial))
                    {
                      GSL_ERROR ("non-finite function value encountered", GSL_EBADFUNC);
                    }
                }

              if (ftrial < fm)
                {
                  gsl_vector_memcpy (x, xtrial);
                  fm = ftrial;
                  break;
                }
              else
                {
                  /* no improvement: retry with a smaller step */
                  gsl_vector_set (state->step_size, i, 0.5 * ss);
                }
            }
          else
            {
              /* the function is flat along this coordinate */
              break;
            }
        }
    }

  *fval = fm;
  *size = quadratic_size (state);

  return GSL_SUCCESS;
}


static void
quadratic_free (void *vstate)
{
  quadratic_state_t *state = (quadratic_state_t *) vstate;

  gsl_vector_free (state->step_size);
  gsl_vector_free (state->xl);
  gsl_vector_free (state->xr);
  gsl_vector_free (state->xtrial);
}


static const gsl_multimin_fminimizer_type quadratic_type =
{ "quadratic",  /* name */
  sizeof (quadratic_state_t),
  &quadratic_alloc,
  &quadratic_set,
  &quadratic_iterate,
  &quadratic_free
};

const gsl_multimin_fminimizer_type
  * gsl_multimin_fminimizer_quadratic = &quadratic_type;
