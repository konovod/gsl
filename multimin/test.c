/* multimin/test.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Fabrice Rossi
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

/* Modified by Tuomo Keskitalo to add Nelder Mead Simplex test suite */

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_ieee_utils.h>

#include "test_funcs.h"

unsigned int fcount, gcount;

int
test_fdf(const char * desc, gsl_multimin_function_fdf *f, 
         initpt_function initpt, const gsl_multimin_fdfminimizer_type *T);

int
test_f(const char * desc, gsl_multimin_function *f, initpt_function initpt,
       const gsl_multimin_fminimizer_type *T);

int
test_quadratic(void);

int
test_error(const gsl_multimin_fminimizer_type * T);

int
main (void)
{
  gsl_ieee_env_setup ();

  {
    const gsl_multimin_fdfminimizer_type *fdfminimizers[6];
    const gsl_multimin_fdfminimizer_type ** T;

    fdfminimizers[0] = gsl_multimin_fdfminimizer_steepest_descent;
    fdfminimizers[1] = gsl_multimin_fdfminimizer_conjugate_pr;
    fdfminimizers[2] = gsl_multimin_fdfminimizer_conjugate_fr;
    fdfminimizers[3] = gsl_multimin_fdfminimizer_vector_bfgs;
    fdfminimizers[4] = gsl_multimin_fdfminimizer_vector_bfgs2;
    fdfminimizers[5] = 0;

    T = fdfminimizers;
    
    while (*T != 0) 
      {
        test_fdf("Roth", &roth, roth_initpt,*T);
        test_fdf("Wood", &wood, wood_initpt,*T);
        test_fdf("Rosenbrock", &rosenbrock, rosenbrock_initpt,*T);
        test_fdf("Rosenbrock1", &rosenbrock, rosenbrock_initpt1,*T);
        test_fdf("SimpleAbs", &simpleabs, simpleabs_initpt,*T);
        T++;
      }
    
    T = fdfminimizers;
    
    while (*T != 0) 
      {
        test_fdf("NRoth", &Nroth, roth_initpt,*T);
        test_fdf("NWood", &Nwood, wood_initpt,*T);
        test_fdf("NRosenbrock", &Nrosenbrock, rosenbrock_initpt,*T);
        T++;
      }
  }


  {
    const gsl_multimin_fminimizer_type *fminimizers[4];
    const gsl_multimin_fminimizer_type ** T;

    fminimizers[0] = gsl_multimin_fminimizer_nmsimplex;
    fminimizers[1] = gsl_multimin_fminimizer_nmsimplex2;
    fminimizers[2] = gsl_multimin_fminimizer_nmsimplex2rand;
    fminimizers[3] = 0;
    
    T = fminimizers;
    
    while (*T != 0) 
      {
        test_f("Roth", &roth_fmin, roth_initpt,*T);
        test_f("Wood", &wood_fmin, wood_initpt,*T);
        test_f("Rosenbrock", &rosenbrock_fmin, rosenbrock_initpt,*T);
        test_f("Spring", &spring_fmin, spring_initpt,*T);
        T++;
      }
  }

  /* The quadratic minimiser is exact on a quadratic function (the
     parabolic step along each coordinate lands on the minimum), so it
     gets its own test rather than being run on the curved-valley test
     functions above. */
  test_quadratic ();

  /* Savannah bug #41527: an error in the objective function must be
     reported with its real error code, not masked as GSL_EFAILED. */
  test_error (gsl_multimin_fminimizer_nmsimplex);
  test_error (gsl_multimin_fminimizer_nmsimplex2);

  exit (gsl_test_summary());
}


/* f(x) = sum_i (i+1) (x_i - 1)^2, minimum 0 at x_i = 1 */
static double
quadratic_fn (const gsl_vector * v, void * params)
{
  size_t i;
  double result = 0.0;

  (void) params;

  for (i = 0; i < v->size; i++)
    {
      const double x = gsl_vector_get (v, i);
      result += (i + 1.0) * (x - 1.0) * (x - 1.0);
    }

  return result;
}

int
test_quadratic (void)
{
  const size_t n = 4;
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_quadratic;
  gsl_multimin_function f;
  gsl_multimin_fminimizer *s;
  gsl_vector *x = gsl_vector_alloc (n);
  gsl_vector *step_size = gsl_vector_alloc (n);
  size_t i, iter = 0;
  int status = GSL_CONTINUE;
  double maxerr = 0.0;

  f.f = &quadratic_fn;
  f.n = n;
  f.params = 0;

  for (i = 0; i < n; i++)
    {
      gsl_vector_set (x, i, 5.0 + i);
      gsl_vector_set (step_size, i, 1.0);
    }

  s = gsl_multimin_fminimizer_alloc (T, n);
  gsl_multimin_fminimizer_set (s, &f, x, step_size);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate (s);
      if (status)
        break;
      status = gsl_multimin_test_size (gsl_multimin_fminimizer_size (s), 1e-8);
    }
  while (status == GSL_CONTINUE && iter < 1000);

  for (i = 0; i < n; i++)
    {
      const double xi = gsl_vector_get (gsl_multimin_fminimizer_x (s), i);
      maxerr = GSL_MAX (maxerr, fabs (xi - 1.0));
    }

  gsl_test (status || maxerr > 1e-6,
            "quadratic, n=%d: %d iter, max |x_i - 1| = %g, f(x)=%g",
            (int) n, (int) iter, maxerr, gsl_multimin_fminimizer_minimum (s));

  gsl_multimin_fminimizer_free (s);
  gsl_vector_free (x);
  gsl_vector_free (step_size);

  return status || maxerr > 1e-6;
}

int
test_fdf(const char * desc, 
         gsl_multimin_function_fdf *f,
         initpt_function initpt,
         const gsl_multimin_fdfminimizer_type *T)
{
  int status;
  size_t iter = 0;
  double step_size;
  
  gsl_vector *x = gsl_vector_alloc (f->n);

  gsl_multimin_fdfminimizer *s;
  fcount = 0; gcount = 0;

  (*initpt) (x);

  step_size = 0.1 * gsl_blas_dnrm2 (x);

  s = gsl_multimin_fdfminimizer_alloc(T, f->n);

  gsl_multimin_fdfminimizer_set (s, f, x, step_size, 0.1);

#ifdef DEBUG
  printf("x "); gsl_vector_fprintf (stdout, s->x, "%g"); 
  printf("g "); gsl_vector_fprintf (stdout, s->gradient, "%g"); 
#endif

  do 
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate(s);

#ifdef DEBUG
      printf("%i: \n",iter);
      printf("x "); gsl_vector_fprintf (stdout, s->x, "%g"); 
      printf("g "); gsl_vector_fprintf (stdout, s->gradient, "%g"); 
      printf("f(x) %g\n",s->f);
      printf("dx %g\n",gsl_blas_dnrm2(s->dx));
      printf("status=%d\n", status);
      printf("\n");
#endif
      if (status == GSL_ENOPROG)
        break;

      status = gsl_multimin_test_gradient(s->gradient,1e-3);
    }
  while (iter < 5000 && status == GSL_CONTINUE);

  /* If no error in iteration, test for numerical convergence */
  if (status == GSL_CONTINUE || status == GSL_ENOPROG) 
    status = (fabs(s->f) > 1e-5);

  gsl_test(status, "%s, on %s: %i iters (fn+g=%d+%d), f(x)=%g",
           gsl_multimin_fdfminimizer_name(s),desc, iter, fcount, gcount, s->f);

  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);

  return status;
}

int
test_f(const char * desc, gsl_multimin_function *f, initpt_function initpt,
       const gsl_multimin_fminimizer_type *T)
{
  int status;
  size_t i, iter = 0;

  gsl_vector *x = gsl_vector_alloc (f->n);

  gsl_vector *step_size = gsl_vector_alloc (f->n);

  gsl_multimin_fminimizer *s;

  fcount = 0; gcount = 0;
  (*initpt) (x);

  for (i = 0; i < f->n; i++) 
    gsl_vector_set (step_size, i, 1);

  s = gsl_multimin_fminimizer_alloc(T, f->n);

  gsl_multimin_fminimizer_set (s, f, x, step_size);

#ifdef DEBUG
  printf("x "); gsl_vector_fprintf (stdout, s->x, "%g"); 
#endif

  do 
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);

#ifdef DEBUG
      printf("%i: \n",iter);
      printf("x "); gsl_vector_fprintf (stdout, s->x, "%g"); 
      printf("f(x) %g\n", gsl_multimin_fminimizer_minimum (s));
      printf("size: %g\n", gsl_multimin_fminimizer_size (s));
      printf("\n");
#endif

      status = gsl_multimin_test_size (gsl_multimin_fminimizer_size (s),
                                       1e-3);
    }
  while (iter < 5000 && status == GSL_CONTINUE);

  status |= (fabs(s->fval) > 1e-5);

  gsl_test(status, "%s, on %s: %d iter (fn=%d), f(x)=%g",
           gsl_multimin_fminimizer_name(s),desc, iter, fcount, s->fval);

  gsl_multimin_fminimizer_free(s);
  gsl_vector_free(x);
  gsl_vector_free(step_size);

  return status;
}

/* Objective that is finite only outside a central square.  The simplex
   reflection and the one-dimensional contraction both land inside the
   square, so the whole-simplex contraction is reached and evaluates the
   midpoint of two corners there, yielding a non-finite value. */
static double
error_fn (const gsl_vector * v, void * params)
{
  const double x0 = gsl_vector_get (v, 0);
  const double x1 = gsl_vector_get (v, 1);

  (void) params;

  if (fabs (x0) < 1.4 && fabs (x1) < 1.4)
    return GSL_NAN;

  return x0 + 2.0 * x1;
}

/* Savannah bug #41527: a non-finite objective encountered during the
   contraction must be reported as GSL_EBADFUNC, not swallowed and
   replaced by GSL_EFAILED. */
int
test_error (const gsl_multimin_fminimizer_type * T)
{
  gsl_multimin_function f;
  gsl_multimin_fminimizer *s;
  gsl_vector *x = gsl_vector_alloc (2);
  gsl_vector *step_size = gsl_vector_alloc (2);
  gsl_error_handler_t *old_handler;
  int status;

  f.f = &error_fn;
  f.n = 2;
  f.params = 0;

  gsl_vector_set (x, 0, 1.5);
  gsl_vector_set (x, 1, 1.5);
  gsl_vector_set (step_size, 0, -1.0);
  gsl_vector_set (step_size, 1, -1.0);

  /* the reported use case runs with the error handler disabled so that
     the returned status can be inspected */
  old_handler = gsl_set_error_handler_off ();

  s = gsl_multimin_fminimizer_alloc (T, 2);
  gsl_multimin_fminimizer_set (s, &f, x, step_size);

  status = gsl_multimin_fminimizer_iterate (s);

  gsl_test (status != GSL_EBADFUNC,
            "%s: non-finite objective gives GSL_EBADFUNC (got %d)",
            gsl_multimin_fminimizer_name (s), status);

  gsl_multimin_fminimizer_free (s);
  gsl_vector_free (x);
  gsl_vector_free (step_size);

  gsl_set_error_handler (old_handler);

  return status != GSL_EBADFUNC;
}
