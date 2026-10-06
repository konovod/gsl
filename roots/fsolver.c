/* roots/fsolver.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Reid Priedhorsky, Brian Gough
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
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

/* Closure used by gsl_root_fsolver_set_with_values().  The solver's
   init routine asks for f(x_lower) and f(x_upper); the wrapper returns
   the caller-supplied values for those two points instead of evaluating
   the (possibly expensive) function again.  Any other point is
   forwarded to the real function. */

typedef struct
{
  gsl_function * f;
  double x_lower, f_lower;
  double x_upper, f_upper;
}
fsolver_wv_closure;

static double
fsolver_wv_eval (double x, void * params)
{
  fsolver_wv_closure * c = (fsolver_wv_closure *) params;

  if (x == c->x_lower)
    return c->f_lower;
  else if (x == c->x_upper)
    return c->f_upper;
  else
    return GSL_FN_EVAL (c->f, x);
}

gsl_root_fsolver *
gsl_root_fsolver_alloc (const gsl_root_fsolver_type * T)
{
  gsl_root_fsolver * s = (gsl_root_fsolver *) malloc (sizeof (gsl_root_fsolver));

  if (s == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for root solver struct",
                        GSL_ENOMEM, 0);
    };

  s->state = malloc (T->size);

  if (s->state == 0)
    {
      free (s);         /* exception in constructor, avoid memory leak */

      GSL_ERROR_VAL ("failed to allocate space for root solver state",
                        GSL_ENOMEM, 0);
    };

  s->type = T ;
  s->function = NULL ;

  return s;
}

int
gsl_root_fsolver_set (gsl_root_fsolver * s, gsl_function * f, double x_lower, double x_upper)
{
  if (x_lower > x_upper)
    {
      GSL_ERROR ("invalid interval (lower > upper)", GSL_EINVAL);
    }

  s->function = f;
  s->root = 0.5 * (x_lower + x_upper);  /* initial estimate */
  s->x_lower = x_lower;
  s->x_upper = x_upper;

  return (s->type->set) (s->state, s->function, &(s->root), x_lower, x_upper);
}

int
gsl_root_fsolver_set_with_values (gsl_root_fsolver * s, gsl_function * f,
                                  double x_lower, double f_lower,
                                  double x_upper, double f_upper)
{
  fsolver_wv_closure c;
  gsl_function wf;

  if (x_lower > x_upper)
    {
      GSL_ERROR ("invalid interval (lower > upper)", GSL_EINVAL);
    }

  s->function = f;
  s->root = 0.5 * (x_lower + x_upper);  /* initial estimate */
  s->x_lower = x_lower;
  s->x_upper = x_upper;

  c.f = f;
  c.x_lower = x_lower;
  c.f_lower = f_lower;
  c.x_upper = x_upper;
  c.f_upper = f_upper;

  wf.function = &fsolver_wv_eval;
  wf.params = &c;

  return (s->type->set) (s->state, &wf, &(s->root), x_lower, x_upper);
}

int
gsl_root_fsolver_iterate (gsl_root_fsolver * s)
{
  return (s->type->iterate) (s->state, 
                             s->function, &(s->root), 
                             &(s->x_lower), &(s->x_upper));
}

void
gsl_root_fsolver_free (gsl_root_fsolver * s)
{
  RETURN_IF_NULL (s);
  free (s->state);
  free (s);
}

const char *
gsl_root_fsolver_name (const gsl_root_fsolver * s)
{
  return s->type->name;
}

double
gsl_root_fsolver_root (const gsl_root_fsolver * s)
{
  return s->root;
}

double
gsl_root_fsolver_x_lower (const gsl_root_fsolver * s)
{
  return s->x_lower;
}

double
gsl_root_fsolver_x_upper (const gsl_root_fsolver * s)
{
  return s->x_upper;
}

