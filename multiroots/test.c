/* multiroots/test.c
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
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_multiroots.h>

#include <gsl/gsl_ieee_utils.h>

#include "test_funcs.h"
int test_fdf (const char * desc, gsl_multiroot_function_fdf * function, initpt_function initpt, double factor, const gsl_multiroot_fdfsolver_type * T);
int test_f (const char * desc, gsl_multiroot_function_fdf * fdf, initpt_function initpt, double factor, const gsl_multiroot_fsolver_type * T);
int test_fdjac_epsrel (void);

#if defined(_MSC_VER)
#include <float.h>
#elif defined(__GLIBC__)
#include <fenv.h>
extern int feenableexcept (int excepts);
#endif

/* Regression tests for the division by zero when the initial point is
   already a root (Savannah bugs #42219 and #42220).

   #42219: gnewton takes phi0 from the fdf callback.  When that returns
   an exact zero residual while the plain f callback returns a nonzero
   one, the relative step reduction divides by phi0.

   #42220: the hybrid methods take a zero dogleg step at a root, and the
   rank-1 update then divides by the zero step norm.  This need not
   change the returned root, so the invalid-operation trap is enabled
   while these cases run to make a regression observable. */

static int
root_point_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  gsl_vector_set (f, 0, gsl_vector_get (x, 0) - 1.0);
  (void) params;
  return GSL_SUCCESS;
}

static int
root_point_df (const gsl_vector * x, void *params, gsl_matrix * df)
{
  gsl_matrix_set (df, 0, 0, 1.0);
  (void) x;
  (void) params;
  return GSL_SUCCESS;
}

static int
root_point_fdf (const gsl_vector * x, void *params,
                gsl_vector * f, gsl_matrix * df)
{
  root_point_f (x, params, f);
  root_point_df (x, params, df);
  return GSL_SUCCESS;
}

static void
root_point_initpt (gsl_vector * x)
{
  gsl_vector_set (x, 0, 1.0);
}

static gsl_multiroot_function_fdf root_point =
  { &root_point_f, &root_point_df, &root_point_fdf, 1, 0 };

static int
gnewton_mismatch_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  gsl_vector_set (f, 0, gsl_vector_get (x, 0) - 1.0 + GSL_DBL_EPSILON);
  (void) params;
  return GSL_SUCCESS;
}

static int
gnewton_mismatch_df (const gsl_vector * x, void *params, gsl_matrix * df)
{
  gsl_matrix_set (df, 0, 0, 1.0);
  (void) x;
  (void) params;
  return GSL_SUCCESS;
}

static int
gnewton_mismatch_fdf (const gsl_vector * x, void *params,
                      gsl_vector * f, gsl_matrix * df)
{
  gsl_vector_set (f, 0, gsl_vector_get (x, 0) - 1.0);
  gsl_matrix_set (df, 0, 0, 1.0);
  (void) params;
  return GSL_SUCCESS;
}

static void
gnewton_mismatch_initpt (gsl_vector * x)
{
  gsl_vector_set (x, 0, 1.0);
}

static gsl_multiroot_function_fdf gnewton_mismatch =
  { &gnewton_mismatch_f, &gnewton_mismatch_df, &gnewton_mismatch_fdf, 1, 0 };

/* A function whose evaluation carries a smooth but high-frequency
   perturbation of order 1e-5, modelling an imprecise function
   (Savannah bug #45782).  The default finite-difference step sees the
   perturbation's large derivative; a larger step averages it out. */

static int
noisy_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double xv = gsl_vector_get (x, 0);
  (void) params;
  gsl_vector_set (f, 0, (xv - 1.0) + 1.0e-5 * sin (1.0e7 * xv));
  return GSL_SUCCESS;
}

static void
noisy_initpt (gsl_vector * x)
{
  gsl_vector_set (x, 0, M_PI / 1.0e7);
}

static void
test_root_point (void)
{
  test_fdf ("mismatched f and fdf", &gnewton_mismatch,
            gnewton_mismatch_initpt, 1.0, gsl_multiroot_fdfsolver_gnewton);

  test_f ("initial point is root", &root_point, root_point_initpt, 1.0,
          gsl_multiroot_fsolver_hybrid);
  test_f ("initial point is root", &root_point, root_point_initpt, 1.0,
          gsl_multiroot_fsolver_hybrids);
  test_fdf ("initial point is root", &root_point, root_point_initpt, 1.0,
            gsl_multiroot_fdfsolver_hybridj);
  test_fdf ("initial point is root", &root_point, root_point_initpt, 1.0,
            gsl_multiroot_fdfsolver_hybridsj);
}


int 
main (void)
{
  const gsl_multiroot_fsolver_type * fsolvers[5] ;
  const gsl_multiroot_fsolver_type ** T1 ;

  const gsl_multiroot_fdfsolver_type * fdfsolvers[5] ;
  const gsl_multiroot_fdfsolver_type ** T2 ;

  double f;

  fsolvers[0] = gsl_multiroot_fsolver_dnewton;
  fsolvers[1] = gsl_multiroot_fsolver_broyden;
  fsolvers[2] = gsl_multiroot_fsolver_hybrid;
  fsolvers[3] = gsl_multiroot_fsolver_hybrids;
  fsolvers[4] = 0;

  fdfsolvers[0] = gsl_multiroot_fdfsolver_newton;
  fdfsolvers[1] = gsl_multiroot_fdfsolver_gnewton;
  fdfsolvers[2] = gsl_multiroot_fdfsolver_hybridj;
  fdfsolvers[3] = gsl_multiroot_fdfsolver_hybridsj;
  fdfsolvers[4] = 0;

  gsl_ieee_env_setup();


  f = 1.0 ;
  
  T1 = fsolvers ;
  
  while (*T1 != 0) 
    {
      test_f ("Rosenbrock", &rosenbrock, rosenbrock_initpt, f, *T1);
      test_f ("Roth", &roth, roth_initpt, f, *T1);
      test_f ("Powell badly scaled", &powellscal, powellscal_initpt, f, *T1);
      test_f ("Brown badly scaled", &brownscal, brownscal_initpt, f, *T1);
      test_f ("Powell singular", &powellsing, powellsing_initpt, f, *T1);
      test_f ("Wood", &wood, wood_initpt, f, *T1);
      test_f ("Helical", &helical, helical_initpt, f, *T1);
      test_f ("Discrete BVP", &dbv, dbv_initpt, f, *T1);
      test_f ("Trig", &trig, trig_initpt, f, *T1);
      T1++;
    }
  
  T2 = fdfsolvers ;
  
  while (*T2 != 0) 
    {
      test_fdf ("Rosenbrock", &rosenbrock, rosenbrock_initpt, f, *T2);
      test_fdf ("Roth", &roth, roth_initpt, f, *T2);
      test_fdf ("Powell badly scaled", &powellscal, powellscal_initpt, f, *T2);
      test_fdf ("Brown badly scaled", &brownscal, brownscal_initpt, f, *T2);
      test_fdf ("Powell singular", &powellsing, powellsing_initpt, f, *T2);
      test_fdf ("Wood", &wood, wood_initpt, f, *T2);
      test_fdf ("Helical", &helical, helical_initpt, f, *T2);
      test_fdf ("Discrete BVP", &dbv, dbv_initpt, f, *T2);
      test_fdf ("Trig", &trig, trig_initpt, f, *T2);
      T2++;
    }

  /* Savannah #42219 and #42220: the solvers must not raise an invalid
     operation when the initial point is already a root.  Enable the
     trap so that the division by zero is observable. */

#if defined(_MSC_VER)
  {
    unsigned int old_cw;
    _controlfp_s (&old_cw, 0, _EM_INVALID | _EM_ZERODIVIDE);
    test_root_point ();
    _controlfp_s (&old_cw, old_cw, _MCW_EM);
  }
#elif defined(__GLIBC__)
  {
    fenv_t env;
    fegetenv (&env);
    feclearexcept (FE_ALL_EXCEPT);
    feenableexcept (FE_INVALID | FE_DIVBYZERO);
    test_root_point ();
    fesetenv (&env);
  }
#else
  test_root_point ();
#endif

  test_fdjac_epsrel ();

  exit (gsl_test_summary ());
}

void scale (gsl_vector * x, double factor);

void
scale (gsl_vector * x, double factor)
{
  size_t i, n = x->size;

  if (gsl_vector_isnull(x))
    {
      for (i = 0; i < n; i++)
        {
          gsl_vector_set (x, i, factor);
        }
    }
  else
    {
      for (i = 0; i < n; i++)
        {
          double xi = gsl_vector_get(x, i);
          gsl_vector_set(x, i, factor * xi);
        }
    } 
}

int
test_fdf (const char * desc, gsl_multiroot_function_fdf * function, 
          initpt_function initpt, double factor,
          const gsl_multiroot_fdfsolver_type * T)
{
  int status;
  double residual = 0;
  size_t i, n = function->n, iter = 0;
  
  gsl_vector *x = gsl_vector_alloc (n);
  gsl_matrix *J = gsl_matrix_alloc (n, n);

  gsl_multiroot_fdfsolver *s;

  (*initpt) (x);

  if (factor != 1.0) scale(x, factor);

  s = gsl_multiroot_fdfsolver_alloc (T, n);
  gsl_multiroot_fdfsolver_set (s, function, x);
 
  do
    {
      iter++;
      status = gsl_multiroot_fdfsolver_iterate (s);
      
      if (status)
        break ;

      status = gsl_multiroot_test_residual (s->f, 0.0000001);
    }
  while (status == GSL_CONTINUE && iter < 1000);

#ifdef DEBUG
  printf("x "); gsl_vector_fprintf (stdout, s->x, "%g"); printf("\n");
  printf("f "); gsl_vector_fprintf (stdout, s->f, "%g"); printf("\n");
#endif


#ifdef TEST_JACOBIAN
 {
    double r,sum; size_t j;

    gsl_multiroot_function f1 ;
    f1.f = function->f ;
    f1.n = function->n ;
    f1.params = function->params ;
    
    gsl_multiroot_fdjacobian (&f1, s->x, s->f, GSL_SQRT_DBL_EPSILON, J);
  
    /* compare J and s->J */
    
    r=0;sum=0;
    for (i = 0; i < n; i++)
      for (j = 0; j< n ; j++)
        {
          double u = gsl_matrix_get(J,i,j);
          double su = gsl_matrix_get(s->J, i, j);
          r = fabs(u - su)/(1e-6 + 1e-6 * fabs(u)); sum+=r;
          if (fabs(u - su) > 1e-6 + 1e-6 * fabs(u))
            printf("broken jacobian %g\n", r);
        }
    printf("avg r = %g\n", sum/(n*n));
  }
#endif

  for (i = 0; i < n ; i++)
    {
      residual += fabs(gsl_vector_get(s->f, i));
    }

  gsl_multiroot_fdfsolver_free (s);
  gsl_matrix_free(J);
  gsl_vector_free(x);

  gsl_test(status, "%s on %s (%g), %u iterations, residual = %.2g", T->name, desc, factor, iter, residual);

  return status;
}


int
test_f (const char * desc, gsl_multiroot_function_fdf * fdf, 
        initpt_function initpt, double factor,
        const gsl_multiroot_fsolver_type * T)
{
  int status;
  size_t i, n = fdf->n, iter = 0;
  double residual = 0;

  gsl_vector *x;

  gsl_multiroot_fsolver *s;
  gsl_multiroot_function function;

  function.f = fdf->f;
  function.params = fdf->params;
  function.n = n ;

  x = gsl_vector_alloc (n);

  (*initpt) (x);

  if (factor != 1.0) scale(x, factor);

  s = gsl_multiroot_fsolver_alloc (T, n);
  gsl_multiroot_fsolver_set (s, &function, x);

/*   printf("x "); gsl_vector_fprintf (stdout, s->x, "%g"); printf("\n"); */
/*   printf("f "); gsl_vector_fprintf (stdout, s->f, "%g"); printf("\n"); */

  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);
      
      if (status)
        break ;

      status = gsl_multiroot_test_residual (s->f, 0.0000001);
    }
  while (status == GSL_CONTINUE && iter < 1000);

#ifdef DEBUG
  printf("x "); gsl_vector_fprintf (stdout, s->x, "%g"); printf("\n");
  printf("f "); gsl_vector_fprintf (stdout, s->f, "%g"); printf("\n");
#endif

  for (i = 0; i < n ; i++)
    {
      residual += fabs(gsl_vector_get(s->f, i));
    }

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free(x);

  gsl_test(status, "%s on %s (%g), %u iterations, residual = %.2g", T->name, desc, factor, iter, residual);

  return status;
}

static int
noisy_solve (double epsrel, int use_epsrel, double * x_final)
{
  int status = GSL_SUCCESS;
  size_t i;
  gsl_multiroot_function function;
  gsl_multiroot_fsolver * s;
  gsl_vector * x = gsl_vector_alloc (1);

  function.f = &noisy_f;
  function.params = 0;
  function.n = 1;

  noisy_initpt (x);

  s = gsl_multiroot_fsolver_alloc (gsl_multiroot_fsolver_dnewton, 1);

  if (use_epsrel)
    {
      status = gsl_multiroot_fsolver_set_fdjac_epsrel (s, epsrel);
      if (status != GSL_SUCCESS)
        {
          gsl_multiroot_fsolver_free (s);
          gsl_vector_free (x);
          return status;
        }
    }

  status = gsl_multiroot_fsolver_set (s, &function, x);

  if (status == GSL_SUCCESS)
    {
      for (i = 0; i < 200; i++)
        {
          status = gsl_multiroot_fsolver_iterate (s);
          if (status != GSL_SUCCESS)
            break;
          if (fabs (gsl_vector_get (s->x, 0) - 1.0) < 1.0e-6)
            break;
        }
    }

  *x_final = gsl_vector_get (s->x, 0);

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);

  return status;
}

int
test_fdjac_epsrel (void)
{
  double x_default = 0.0, x_tuned = 0.0;
  gsl_multiroot_fsolver * solver;

  /* default value, validation and storage */
  solver = gsl_multiroot_fsolver_alloc (gsl_multiroot_fsolver_dnewton, 1);
  gsl_test (gsl_multiroot_fsolver_fdjac_epsrel (solver) != GSL_SQRT_DBL_EPSILON,
            "fsolver default fdjac epsrel");
  {
    gsl_error_handler_t * old = gsl_set_error_handler_off ();
    gsl_test (gsl_multiroot_fsolver_set_fdjac_epsrel (solver, -1.0) != GSL_EDOM,
              "fsolver rejects a non-positive fdjac epsrel");
    gsl_set_error_handler (old);
  }
  gsl_test (gsl_multiroot_fsolver_set_fdjac_epsrel (solver, 1.0e-3) != GSL_SUCCESS,
            "fsolver accepts a positive fdjac epsrel");
  gsl_test (gsl_multiroot_fsolver_fdjac_epsrel (solver) != 1.0e-3,
            "fsolver stores the fdjac epsrel");
  gsl_multiroot_fsolver_free (solver);

  /* a noisy function defeats the default step but not a larger one */
  noisy_solve (0.0, 0, &x_default);
  noisy_solve (1.0e-3, 1, &x_tuned);

  gsl_test (fabs (x_tuned - 1.0) > 1.0e-3,
            "tuned fdjac epsrel reaches the root: x = %.17g", x_tuned);
  gsl_test (fabs (x_default - 1.0) < 1.0e-3,
            "default fdjac epsrel is deflected by the noise: x = %.17g", x_default);

  return 0;
}
