/* linalg/lsi.c
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

/* Linear least squares with linear inequality constraints, after
 * C. L. Lawson and R. J. Hanson, "Solving Least Squares Problems"
 * (SIAM, 1995), chapters 23 and 24:
 *
 *   gsl_linalg_nnls  min ||A x - b||   subject to x >= 0
 *   gsl_linalg_ldp   min ||x||         subject to G x >= h
 *   gsl_linalg_lsi   min ||A x - b||   subject to G x >= h
 *
 * NNLS is the classical active-set algorithm.  LDP is reduced to NNLS
 * as in Algorithm 23.17: with E = [G^T; h^T] ((n+1) by m) and
 * f = [0; ...; 0; 1] (length n+1), the NNLS solution y of
 * min ||E y - f||, y >= 0, gives x = (G^T y) / (1 - h^T y).
 * LSI runs a primal active-set method for the convex quadratic program
 * min ||A x - b||^2 s.t. G x >= h, starting from the LDP feasible
 * point.
 */

#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

/* Solve min ||A_P s - b|| over the columns of A in the index list
 * `set` (of length p), placing the result in s (length p).  `sub` and
 * `tau` and `res` are caller-supplied working storage of the sizes
 * m x p, p and m respectively. */
static int
solve_ls_subset(const gsl_matrix * A, const gsl_vector * b,
                const size_t * set, size_t p,
                gsl_matrix * sub, gsl_vector * tau,
                gsl_vector * s, gsl_vector * res)
{
  size_t i, j;

  for (j = 0; j < p; j++)
    {
      gsl_vector_const_view col = gsl_matrix_const_column (A, set[j]);
      gsl_matrix_set_col (sub, j, &col.vector);
    }

  {
    gsl_matrix_view sv = gsl_matrix_submatrix (sub, 0, 0, A->size1, p);
    gsl_vector_view tv = gsl_vector_subvector (tau, 0, p);
    gsl_vector_view xv = gsl_vector_subvector (s, 0, p);
    int status;

    status = gsl_linalg_QR_decomp (&sv.matrix, &tv.vector);
    if (status != GSL_SUCCESS)
      return status;

    status = gsl_linalg_QR_lssolve (&sv.matrix, &tv.vector, b,
                                    &xv.vector, res);
    if (status != GSL_SUCCESS)
      return status;

    /* the remaining entries of s are meaningless; the caller zeroes
       them for the variables outside the set */
    (void) i;
  }

  return GSL_SUCCESS;
}

/* Lawson-Hanson nonnegative least squares: min ||A x - b||, x >= 0. */
static int
nnls_core(const gsl_matrix * A, const gsl_vector * b, gsl_vector * x)
{
  const size_t m = A->size1;
  const size_t n = A->size2;
  int status = GSL_SUCCESS;
  size_t * set = (size_t *) malloc (n * sizeof (size_t));
  unsigned char * in_set = (unsigned char *) calloc (n, 1);
  gsl_matrix * sub = gsl_matrix_alloc (m, n > 0 ? n : 1);
  gsl_vector * tau = gsl_vector_alloc (n > 0 ? n : 1);
  gsl_vector * s = gsl_vector_alloc (n > 0 ? n : 1);
  gsl_vector * res = gsl_vector_alloc (m);
  gsl_vector * r = gsl_vector_alloc (m);
  gsl_vector * w = gsl_vector_alloc (n);
  gsl_vector * Atb = gsl_vector_alloc (n);
  double scale, tol;
  size_t p = 0, outer, inner;
  size_t j;

  if (set == 0 || in_set == 0 || sub == 0 || tau == 0 || s == 0 || res == 0
      || r == 0 || w == 0 || Atb == 0)
    {
      GSL_ERROR ("failed to allocate working space", GSL_ENOMEM);
    }

  gsl_vector_set_zero (x);

  /* scale for the termination test: ||A^T b||_inf */
  {
    gsl_blas_dgemv (CblasTrans, 1.0, A, b, 0.0, Atb);
    scale = 0.0;
    for (j = 0; j < n; j++)
      {
        double a = fabs (gsl_vector_get (Atb, j));
        if (a > scale)
          scale = a;
      }
  }
  tol = 1.0e-13 * (1.0 + scale);

  for (outer = 0; outer < 3 * n + 10 && status == GSL_SUCCESS; outer++)
    {
      size_t t = 0;
      double wmax = 0.0;

      /* dual vector w = A^T (b - A x) */
      gsl_blas_dgemv (CblasNoTrans, 1.0, A, x, 0.0, r);
      gsl_vector_sub (r, b);            /* r = A x - b */
      gsl_blas_dgemv (CblasTrans, -1.0, A, r, 0.0, w);  /* w = -A^T r */

      /* pick the index outside the set with the largest positive w */
      for (j = 0; j < n; j++)
        {
          if (!in_set[j] && gsl_vector_get (w, j) > wmax)
            {
              wmax = gsl_vector_get (w, j);
              t = j;
            }
        }

      if (wmax <= tol)
        break;                          /* KKT conditions satisfied */

      /* move t into the passive set */
      set[p] = t;
      in_set[t] = 1;
      p++;

      for (inner = 0; inner < n + 2; inner++)
        {
          size_t k;
          double alpha = 2.0;

          status = solve_ls_subset (A, b, set, p, sub, tau, s, res);
          if (status != GSL_SUCCESS)
            break;

          /* check feasibility of the unconstrained solution */
          for (k = 0; k < p; k++)
            {
              double sk = gsl_vector_get (s, k);
              double xk = gsl_vector_get (x, set[k]);
              if (sk <= 0.0)
                {
                  double a = xk / (xk - sk);
                  if (a < alpha)
                    {
                      alpha = a;
                    }
                }
            }

          if (alpha == 2.0)
            {
              /* all positive: accept the solution */
              for (k = 0; k < p; k++)
                gsl_vector_set (x, set[k], gsl_vector_get (s, k));
              break;
            }

          /* interpolate towards the feasible boundary */
          for (k = 0; k < p; k++)
            {
              double xk = gsl_vector_get (x, set[k]);
              double sk = gsl_vector_get (s, k);
              gsl_vector_set (x, set[k], xk + alpha * (sk - xk));
            }

          /* move the variables that reached zero back to the active set */
          for (k = 0; k < p; )
            {
              size_t v = set[k];
              if (gsl_vector_get (x, v) <= tol)
                {
                  gsl_vector_set (x, v, 0.0);
                  in_set[v] = 0;
                  set[k] = set[p - 1];
                  p--;
                }
              else
                k++;
            }

          if (p == 0)
            break;                      /* nothing left: x stays feasible */
        }

      if (inner == n + 2)
        status = GSL_EMAXITER;
    }

  if (outer == 3 * n + 10)
    status = GSL_EMAXITER;

  free (set);
  free (in_set);
  gsl_matrix_free (sub);
  gsl_vector_free (tau);
  gsl_vector_free (s);
  gsl_vector_free (res);
  gsl_vector_free (r);
  gsl_vector_free (w);
  gsl_vector_free (Atb);

  return status;
}

/* least distance programming: min ||x|| subject to G x >= h */
int
gsl_linalg_ldp(const gsl_matrix * G, const gsl_vector * h, gsl_vector * x,
               gsl_vector * residual)
{
  const size_t m = G->size1;    /* number of constraints */
  const size_t n = G->size2;    /* number of variables */
  gsl_matrix * E;
  gsl_vector * f;
  gsl_vector * y;
  int status;
  size_t i, j;
  double fac;

  if (G->size1 != h->size)
    GSL_ERROR ("constraint matrix and vector sizes must agree", GSL_EBADLEN);
  else if (x->size != n)
    GSL_ERROR ("solution vector has wrong size", GSL_EBADLEN);
  else if (n == 0)
    GSL_ERROR ("number of variables must be positive", GSL_EBADLEN);

  E = gsl_matrix_alloc (n + 1, m > 0 ? m : 1);
  f = gsl_vector_alloc (n + 1);
  y = gsl_vector_alloc (m > 0 ? m : 1);

  if (E == 0 || f == 0 || y == 0)
    {
      gsl_matrix_free (E);
      gsl_vector_free (f);
      gsl_vector_free (y);
      GSL_ERROR ("failed to allocate working space", GSL_ENOMEM);
    }

  /* E = [G^T; h^T], f = [0; ...; 0; 1] */
  for (j = 0; j < n; j++)
    for (i = 0; i < m; i++)
      gsl_matrix_set (E, j, i, gsl_matrix_get (G, i, j));
  for (i = 0; i < m; i++)
    gsl_matrix_set (E, n, i, gsl_vector_get (h, i));

  gsl_vector_set_zero (f);
  gsl_vector_set (f, n, 1.0);

  status = nnls_core (E, f, y);
  if (status != GSL_SUCCESS)
    {
      gsl_matrix_free (E);
      gsl_vector_free (f);
      gsl_vector_free (y);
      return status;
    }

  /* fac = 1 - h^T y */
  fac = 1.0;
  for (i = 0; i < m; i++)
    fac -= gsl_vector_get (h, i) * gsl_vector_get (y, i);

  if (fabs (fac) <= GSL_DBL_EPSILON * (1.0 + fabs (fac)))
    {
      /* 1 - h^T y vanishes: the constraints are not compatible */
      gsl_matrix_free (E);
      gsl_vector_free (f);
      gsl_vector_free (y);
      GSL_ERROR ("constraints are not compatible", GSL_EDOM);
    }

  /* x = (G^T y) / fac */
  for (j = 0; j < n; j++)
    {
      double s = 0.0;
      for (i = 0; i < m; i++)
        s += gsl_matrix_get (G, i, j) * gsl_vector_get (y, i);
      gsl_vector_set (x, j, s / fac);
    }

  if (residual != 0)
    {
      if (residual->size != m)
        {
          gsl_matrix_free (E);
          gsl_vector_free (f);
          gsl_vector_free (y);
          GSL_ERROR ("residual vector has wrong size", GSL_EBADLEN);
        }
      for (i = 0; i < m; i++)
        {
          double s = 0.0;
          for (j = 0; j < n; j++)
            s += gsl_matrix_get (G, i, j) * gsl_vector_get (x, j);
          gsl_vector_set (residual, i, s - gsl_vector_get (h, i));
        }
    }

  gsl_matrix_free (E);
  gsl_vector_free (f);
  gsl_vector_free (y);
  return GSL_SUCCESS;
}

/* linear least squares with inequality constraints:
   min ||A x - b|| subject to G x >= h */
int
gsl_linalg_lsi(const gsl_matrix * A, const gsl_vector * b,
               const gsl_matrix * G, const gsl_vector * h,
               gsl_vector * x, gsl_vector * residual)
{
  const size_t m = A->size1;
  const size_t n = A->size2;
  const size_t r = G->size1;
  int status = GSL_SUCCESS;
  unsigned char * active = (unsigned char *) calloc (r > 0 ? r : 1, 1);
  gsl_matrix * H = gsl_matrix_alloc (n, n);
  gsl_matrix * K = gsl_matrix_alloc (n + r, n + r);
  gsl_matrix * V = gsl_matrix_alloc (n + r, n + r);
  gsl_vector * S = gsl_vector_alloc (n + r);
  gsl_vector * work = gsl_vector_alloc (n + r);
  gsl_vector * rhs = gsl_vector_alloc (n + r);
  gsl_vector * sol = gsl_vector_alloc (n + r);
  gsl_vector * p = gsl_vector_alloc (n);
  gsl_vector * rr = gsl_vector_alloc (m);
  gsl_vector * grad = gsl_vector_alloc (n);
  gsl_vector * Gx = gsl_vector_alloc (r > 0 ? r : 1);
  gsl_vector * Gp = gsl_vector_alloc (r > 0 ? r : 1);
  double gscale, tol_step, tol_lam;
  size_t iter;

  if (active == 0 || H == 0 || K == 0 || V == 0 || S == 0 || work == 0
      || rhs == 0 || sol == 0 || p == 0 || rr == 0 || grad == 0 || Gx == 0
      || Gp == 0)
    {
      GSL_ERROR ("failed to allocate working space", GSL_ENOMEM);
    }

  /* H = A^T A */
  gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, A, A, 0.0, H);

  /* phase I: the LDP minimiser is a feasible point */
  status = gsl_linalg_ldp (G, h, x, 0);
  if (status != GSL_SUCCESS)
    {
      if (status == GSL_EDOM)
        GSL_ERROR ("constraints are not compatible", GSL_EDOM);
      return status;
    }

  /* the active set is the constraints that the LDP point meets */
  {
    size_t i;
    gsl_blas_dgemv (CblasNoTrans, 1.0, G, x, 0.0, Gx);
    for (i = 0; i < r; i++)
      {
        double slack = gsl_vector_get (Gx, i) - gsl_vector_get (h, i);
        if (slack <= 1.0e-10 * (1.0 + fabs (gsl_vector_get (h, i))))
          active[i] = 1;
      }
  }

  gscale = 1.0;
  for (iter = 0; iter < 3 * (n + r) + 30 && status == GSL_SUCCESS; iter++)
    {
      size_t i, nw = 0;
      double pnorm;
      double mu_max = -1.0e300;
      size_t mu_idx = 0;
      int have_mu = 0;

      gsl_blas_dgemv (CblasNoTrans, 1.0, A, x, 0.0, rr);
      gsl_vector_sub (rr, b);                 /* rr = A x - b */
      gsl_blas_dgemv (CblasTrans, 1.0, A, rr, 0.0, grad);  /* A^T rr */

      /* assemble the KKT system */
      gsl_matrix_set_zero (K);
      {
        gsl_matrix_view Hv = gsl_matrix_submatrix (K, 0, 0, n, n);
        gsl_matrix_memcpy (&Hv.matrix, H);
      }
      for (i = 0; i < r; i++)
        {
          if (active[i])
            {
              size_t j;
              for (j = 0; j < n; j++)
                {
                  double g = gsl_matrix_get (G, i, j);
                  gsl_matrix_set (K, j, n + nw, g);
                  gsl_matrix_set (K, n + nw, j, g);
                }
              nw++;
            }
        }
      gsl_vector_set_zero (rhs);
      {
        gsl_vector_view rv = gsl_vector_subvector (rhs, 0, n);
        gsl_vector_memcpy (&rv.vector, grad);
        gsl_vector_scale (&rv.vector, -1.0);
      }

      status = gsl_linalg_SV_decomp (K, V, S, work);
      if (status != GSL_SUCCESS)
        break;
      status = gsl_linalg_SV_solve (K, V, S, rhs, sol);
      if (status != GSL_SUCCESS)
        break;

      {
        gsl_vector_view pv = gsl_vector_subvector (sol, 0, n);
        gsl_vector_memcpy (p, &pv.vector);
      }

      pnorm = gsl_blas_dnrm2 (p);

      /* multipliers for the active constraints: the inequality
         multiplier is lambda = -mu, so a negative lambda (to be
         dropped) is a positive mu */
      for (i = 0; i < nw; i++)
        {
          double mu = gsl_vector_get (sol, n + i);
          if (!have_mu || mu > mu_max)
            {
              mu_max = mu;
              mu_idx = i;
              have_mu = 1;
            }
        }

      tol_step = 1.0e-12 * (1.0 + gsl_blas_dnrm2 (x));
      tol_lam = 1.0e-10 * (1.0 + gscale);

      if (pnorm <= tol_step)
        {
          if (have_mu && mu_max > tol_lam)
            {
              /* drop the constraint with the most negative multiplier */
              size_t seen = 0;
              for (i = 0; i < r; i++)
                {
                  if (active[i])
                    {
                      if (seen == mu_idx)
                        {
                          active[i] = 0;
                          break;
                        }
                      seen++;
                    }
                }
              continue;
            }
          break;                              /* optimal */
        }

      /* step to the nearest blocking constraint */
      {
        double alpha = 1.0;
        double gpmax = 0.0;
        size_t block = (size_t) -1;

        gsl_blas_dgemv (CblasNoTrans, 1.0, G, p, 0.0, Gp);
        gsl_blas_dgemv (CblasNoTrans, 1.0, G, x, 0.0, Gx);

        for (i = 0; i < r; i++)
          {
            double a = fabs (gsl_vector_get (Gp, i));
            if (a > gpmax)
              gpmax = a;
          }

        for (i = 0; i < r; i++)
          {
            double gpi = gsl_vector_get (Gp, i);
            if (!active[i] && gpi < -1.0e-13 * (1.0 + gpmax))
              {
                double slack = gsl_vector_get (Gx, i) - gsl_vector_get (h, i);
                double t = -slack / gpi;
                if (t < alpha)
                  {
                    alpha = t;
                    block = i;
                  }
              }
          }

        if (alpha < 0.0)
          alpha = 0.0;
        gsl_vector_scale (p, alpha);
        gsl_vector_add (x, p);

        if (alpha < 1.0 && block != (size_t) -1)
          active[block] = 1;
      }
    }

  if (iter == 3 * (n + r) + 30)
    status = GSL_EMAXITER;

  if (status == GSL_SUCCESS && residual != 0)
    {
      if (residual->size != m)
        {
          free (active);
          gsl_matrix_free (H);
          gsl_matrix_free (K);
          gsl_matrix_free (V);
          gsl_vector_free (S);
          gsl_vector_free (work);
          gsl_vector_free (rhs);
          gsl_vector_free (sol);
          gsl_vector_free (p);
          gsl_vector_free (rr);
          gsl_vector_free (grad);
          gsl_vector_free (Gx);
          gsl_vector_free (Gp);
          GSL_ERROR ("residual vector has wrong size", GSL_EBADLEN);
        }
      gsl_blas_dgemv (CblasNoTrans, 1.0, A, x, 0.0, residual);
      gsl_vector_sub (residual, b);
    }

  free (active);
  gsl_matrix_free (H);
  gsl_matrix_free (K);
  gsl_matrix_free (V);
  gsl_vector_free (S);
  gsl_vector_free (work);
  gsl_vector_free (rhs);
  gsl_vector_free (sol);
  gsl_vector_free (p);
  gsl_vector_free (rr);
  gsl_vector_free (grad);
  gsl_vector_free (Gx);
  gsl_vector_free (Gp);

  return status;
}

int
gsl_linalg_nnls(const gsl_matrix * A, const gsl_vector * b, gsl_vector * x,
                gsl_vector * residual)
{
  int status;

  if (A->size1 != b->size)
    GSL_ERROR ("matrix and vector sizes must agree", GSL_EBADLEN);
  else if (A->size2 != x->size)
    GSL_ERROR ("solution vector has wrong size", GSL_EBADLEN);

  status = nnls_core (A, b, x);
  if (status != GSL_SUCCESS)
    return status;

  if (residual != 0)
    {
      if (residual->size != A->size1)
        GSL_ERROR ("residual vector has wrong size", GSL_EBADLEN);
      gsl_blas_dgemv (CblasNoTrans, 1.0, A, x, 0.0, residual);
      gsl_vector_sub (residual, b);
    }

  return GSL_SUCCESS;
}
