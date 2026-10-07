/* ode-initval/cfxd.c
 * 
 * Copyright (C) 2010 Sergey B Kirpichev
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

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>

/* The fixed-step controller never changes the step size, so it carries
   no state; the alloc routine reserves a single byte. */

static void *
fc_control_alloc (void)
{
  char * s = (char *) malloc (1);

  if (s == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for fixed control state",
                      GSL_ENOMEM);
    }

  return s;
}

static int
fc_control_init (void * vstate, double eps_abs, double eps_rel,
                 double a_y, double a_dydt)
{
  (void) vstate;
  (void) eps_abs;
  (void) eps_rel;
  (void) a_y;
  (void) a_dydt;

  return GSL_SUCCESS;
}

static int
fc_control_hadjust (void * vstate, size_t dim, unsigned int ord,
                    const double y[], const double yerr[],
                    const double yp[], double * h)
{
  (void) vstate;
  (void) dim;
  (void) ord;
  (void) y;
  (void) yerr;
  (void) yp;
  (void) h;

  return GSL_ODEIV_HADJ_NIL;    /* never change the step size */
}

static void
fc_control_free (void * vstate)
{
  free (vstate);
}

static const gsl_odeiv_control_type fc_control_type = {
  "fixed",                      /* name */
  &fc_control_alloc,
  &fc_control_init,
  &fc_control_hadjust,
  &fc_control_free
};

const gsl_odeiv_control_type *gsl_odeiv_control_fixed = &fc_control_type;

gsl_odeiv_control *
gsl_odeiv_control_fixed_new (void)
{
  return gsl_odeiv_control_alloc (gsl_odeiv_control_fixed);
}
