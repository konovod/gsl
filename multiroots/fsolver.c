/* multiroots/fsolver.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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
#include <gsl/gsl_multiroots.h>

gsl_multiroot_fsolver * 
gsl_multiroot_fsolver_alloc (const gsl_multiroot_fsolver_type * T,
                                     size_t n)
{
  int status;

  gsl_multiroot_fsolver * s;

  s = (gsl_multiroot_fsolver *) malloc (sizeof (gsl_multiroot_fsolver));

  if (s == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for multiroot solver struct",
                        GSL_ENOMEM, 0);
    }

  s->x = gsl_vector_calloc (n);

  if (s->x == 0) 
    {
      free (s);
      GSL_ERROR_VAL ("failed to allocate space for x", GSL_ENOMEM, 0);
    }

  s->f = gsl_vector_calloc (n);

  if (s->f == 0) 
    {
      gsl_vector_free (s->x);
      free (s);
      GSL_ERROR_VAL ("failed to allocate space for f", GSL_ENOMEM, 0);
    }

  s->dx = gsl_vector_calloc (n);

  if (s->dx == 0) 
    {
      gsl_vector_free (s->x);
      gsl_vector_free (s->f);
      free (s);
      GSL_ERROR_VAL ("failed to allocate space for dx", GSL_ENOMEM, 0);
    }

  s->state = malloc (T->size);

  if (s->state == 0)
    {
      gsl_vector_free (s->dx);
      gsl_vector_free (s->x);
      gsl_vector_free (s->f);
      free (s);         /* exception in constructor, avoid memory leak */
      
      GSL_ERROR_VAL ("failed to allocate space for multiroot solver state",
                        GSL_ENOMEM, 0);
    }

  s->type = T ;

  s->fdjac_epsrel = GSL_SQRT_DBL_EPSILON;

  status = (s->type->alloc)(s->state, n);

  if (status != GSL_SUCCESS)
    {
      (s->type->free)(s->state);
      free (s->state);
      gsl_vector_free (s->dx);
      gsl_vector_free (s->x);
      gsl_vector_free (s->f);
      free (s);         /* exception in constructor, avoid memory leak */
      
      GSL_ERROR_VAL ("failed to set solver", status, 0);
    }

  s->function = NULL;
  
  return s;
}

int
gsl_multiroot_fsolver_set (gsl_multiroot_fsolver * s, 
                           gsl_multiroot_function * f, 
                           const gsl_vector * x)
{
  if (s->x->size != f->n)
    {
      GSL_ERROR ("function incompatible with solver size", GSL_EBADLEN);
    }
  
  if (x->size != f->n)
    {
      GSL_ERROR ("vector length not compatible with function", GSL_EBADLEN);
    }  
  
  s->function = f;
  gsl_vector_memcpy(s->x,x);
  
  return (s->type->set) (s->state, s->function, s->x, s->f, s->dx,
                         s->fdjac_epsrel);
}

int
gsl_multiroot_fsolver_iterate (gsl_multiroot_fsolver * s)
{
  return (s->type->iterate) (s->state, s->function, s->x, s->f, s->dx,
                             s->fdjac_epsrel);
}

int
gsl_multiroot_fsolver_set_fdjac_epsrel (gsl_multiroot_fsolver * s,
                                        double epsrel)
{
  if (epsrel <= 0.0)
    {
      GSL_ERROR ("fdjac epsrel must be positive", GSL_EDOM);
    }

  s->fdjac_epsrel = epsrel;

  return GSL_SUCCESS;
}

double
gsl_multiroot_fsolver_fdjac_epsrel (const gsl_multiroot_fsolver * s)
{
  return s->fdjac_epsrel;
}

void
gsl_multiroot_fsolver_free (gsl_multiroot_fsolver * s)
{
  RETURN_IF_NULL (s);
  (s->type->free) (s->state);
  free (s->state);
  gsl_vector_free (s->dx);
  gsl_vector_free (s->x);
  gsl_vector_free (s->f);
  free (s);
}

const char *
gsl_multiroot_fsolver_name (const gsl_multiroot_fsolver * s)
{
  return s->type->name;
}

gsl_vector *
gsl_multiroot_fsolver_root (const gsl_multiroot_fsolver * s)
{
  return s->x;
}

gsl_vector *
gsl_multiroot_fsolver_dx (const gsl_multiroot_fsolver * s)
{
  return s->dx;
}

gsl_vector *
gsl_multiroot_fsolver_f (const gsl_multiroot_fsolver * s)
{
  return s->f;
}
