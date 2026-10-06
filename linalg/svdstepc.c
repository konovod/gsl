/* linalg/svdstepc.c 
 *
 * Copyright (C) 2007, 2010 Brian Gough
 * Copyright (C) 2024 Christian Krueger (Adjusted)
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

static void
svd2_complex (gsl_vector *d, gsl_vector *f, gsl_matrix_complex *U,
              gsl_matrix_complex *V)
{
  size_t i;
  double c, s, a11, a12, a21, a22;

  const size_t M = U->size1;
  const size_t N = V->size1;

  double d0 = gsl_vector_get (d, 0);
  double f0 = gsl_vector_get (f, 0);

  double d1 = gsl_vector_get (d, 1);

  if (d0 == 0.0)
    {
      /* Eliminate off-diagonal element in [0,f0;0,d1] to make [d,0;0,0] */

      gsl_linalg_givens (f0, d1, &c, &s);

      /* compute B <= G^T B X,  where X = [0,1;1,0] */

      gsl_vector_set (d, 0, c * f0 - s * d1);
      gsl_vector_set (f, 0, s * f0 + c * d1);
      gsl_vector_set (d, 1, 0.0);

      /* Compute U <= U G */

      for (i = 0; i < M; i++)
        {
          gsl_complex Uip = gsl_matrix_complex_get (U, i, 0);
          gsl_complex Uiq = gsl_matrix_complex_get (U, i, 1);

          gsl_complex h1 = gsl_complex_mul_real (Uip, c);
          gsl_complex h2 = gsl_complex_mul_real (Uiq, s);
          gsl_matrix_complex_set (U, i, 0, gsl_complex_sub (h1, h2));

          h1 = gsl_complex_mul_real (Uip, s);
          h2 = gsl_complex_mul_real (Uiq, c);
          gsl_matrix_complex_set (U, i, 1, gsl_complex_add (h1, h2));
        }

      /* Compute V <= V X */

      gsl_matrix_complex_swap_columns (V, 0, 1);

      return;
    }
  else if (d1 == 0.0)
    {
      /* Eliminate off-diagonal element in [d0,f0;0,0] */

      gsl_linalg_givens (d0, f0, &c, &s);

      /* compute B <= B G */

      gsl_vector_set (d, 0, d0 * c - f0 * s);
      gsl_vector_set (f, 0, 0.0);

      /* Compute V <= V G */

      for (i = 0; i < N; i++)
        {
          gsl_complex Vip = gsl_matrix_complex_get (V, i, 0);
          gsl_complex Viq = gsl_matrix_complex_get (V, i, 1);

          gsl_complex h1 = gsl_complex_mul_real (Vip, c);
          gsl_complex h2 = gsl_complex_mul_real (Viq, s);
          gsl_matrix_complex_set (V, i, 0, gsl_complex_sub (h1, h2));

          h1 = gsl_complex_mul_real (Vip, s);
          h2 = gsl_complex_mul_real (Viq, c);
          gsl_matrix_complex_set (V, i, 1, gsl_complex_add (h1, h2));
        }

      return;
    }
  else
    {
      /* Make columns orthogonal, A = [d0, f0; 0, d1] * G */

      create_schur (d0, f0, d1, &c, &s);

      /* compute B <= B G */

      a11 = c * d0 - s * f0;
      a21 = -s * d1;

      a12 = s * d0 + c * f0;
      a22 = c * d1;

      /* Compute V <= V G */

      for (i = 0; i < N; i++)
        {
          gsl_complex Vip = gsl_matrix_complex_get (V, i, 0);
          gsl_complex Viq = gsl_matrix_complex_get (V, i, 1);

          gsl_complex h1 = gsl_complex_mul_real (Vip, c);
          gsl_complex h2 = gsl_complex_mul_real (Viq, s);
          gsl_matrix_complex_set (V, i, 0, gsl_complex_sub (h1, h2));

          h1 = gsl_complex_mul_real (Vip, s);
          h2 = gsl_complex_mul_real (Viq, c);
          gsl_matrix_complex_set (V, i, 1, gsl_complex_add (h1, h2));
        }

      /* Eliminate off-diagonal elements, bring column with largest
         norm to first column */

      if (hypot (a11, a21) < hypot (a12, a22))
        {
          double t1, t2;

          /* B <= B X */

          t1 = a11;
          a11 = a12;
          a12 = t1;
          t2 = a21;
          a21 = a22;
          a22 = t2;

          /* V <= V X */

          gsl_matrix_complex_swap_columns (V, 0, 1);
        }

      gsl_linalg_givens (a11, a21, &c, &s);

      /* compute B <= G^T B */

      gsl_vector_set (d, 0, c * a11 - s * a21);
      gsl_vector_set (f, 0, c * a12 - s * a22);
      gsl_vector_set (d, 1, s * a12 + c * a22);

      /* Compute U <= U G */

      for (i = 0; i < M; i++)
        {
          gsl_complex Uip = gsl_matrix_complex_get (U, i, 0);
          gsl_complex Uiq = gsl_matrix_complex_get (U, i, 1);

          gsl_complex h1 = gsl_complex_mul_real (Uip, c);
          gsl_complex h2 = gsl_complex_mul_real (Uiq, s);
          gsl_matrix_complex_set (U, i, 0, gsl_complex_sub (h1, h2));

          h1 = gsl_complex_mul_real (Uip, s);
          h2 = gsl_complex_mul_real (Uiq, c);
          gsl_matrix_complex_set (U, i, 1, gsl_complex_add (h1, h2));
        }

      return;
    }
}


static void
chase_out_intermediate_zero_complex (gsl_vector *d, gsl_vector *f,
                             gsl_matrix_complex *U, size_t k0)
{
#if !USE_BLAS
  const size_t M = U->size1;
#endif
  const size_t n = d->size;
  double c, s;
  double x, y;
  size_t k;

  x = gsl_vector_get (f, k0);
  y = gsl_vector_get (d, k0 + 1);

  for (k = k0; k < n - 1; k++)
    {
      gsl_linalg_givens (y, -x, &c, &s);

      /* Compute U <= U G */

// #ifdef USE_BLAS
//       {
//         gsl_vector_complex_view Uk0 = gsl_matrix_complex_column (U, k0);
//         gsl_vector_complex_view Ukp1 = gsl_matrix_complex_column (U, k + 1);
//         gsl_blas_zdrot (&Uk0.vector, &Ukp1.vector, c, -s);
//         // TODO blas_zdrot does not exist yet
//       }
// #else
      {
        size_t i;

        for (i = 0; i < M; i++)
          {
            gsl_complex Uip = gsl_matrix_complex_get (U, i, k0);
            gsl_complex Uiq = gsl_matrix_complex_get (U, i, k + 1);

            gsl_complex h1 = gsl_complex_mul_real (Uip, c);
            gsl_complex h2 = gsl_complex_mul_real (Uiq, s);
            gsl_matrix_complex_set (U, i, k0, gsl_complex_sub (h1, h2));

            h1 = gsl_complex_mul_real (Uip, s);
            h2 = gsl_complex_mul_real (Uiq, c);
            gsl_matrix_complex_set (U, i, k + 1, gsl_complex_add (h1, h2));
          }
      }
// #endif

      /* compute B <= G^T B */

      gsl_vector_set (d, k + 1, s * x + c * y);

      if (k == k0)
        gsl_vector_set (f, k, c * x - s * y);

      if (k < n - 2)
        {
          double z = gsl_vector_get (f, k + 1);
          gsl_vector_set (f, k + 1, c * z);

          x = -s * z;
          y = gsl_vector_get (d, k + 2);
        }
    }
}

static void
chase_out_trailing_zero_complex (gsl_vector *d, gsl_vector *f, gsl_matrix_complex *V)
{
#if !USE_BLAS
  const size_t N = V->size1;
#endif
  const size_t n = d->size;
  double c, s;
  double x, y;
  size_t k;

  x = gsl_vector_get (d, n - 2);
  y = gsl_vector_get (f, n - 2);

  for (k = n - 1; k-- > 0;)
    {
      gsl_linalg_givens (x, y, &c, &s);

      /* Compute V <= V G where G = [c, s ; -s, c] */

// #ifdef USE_BLAS
//       {
//         gsl_vector_complex_view Vp = gsl_matrix_complex_column (V, k);
//         gsl_vector_complex_view Vq = gsl_matrix_complex_column (V, n - 1);
//         gsl_blas_zdrot (&Vp.vector, &Vq.vector, c, -s);
//         // TODO blas_zdrot does not exist yet
//       }
// #else
      {
        size_t i;

        for (i = 0; i < N; i++)
          {
            gsl_complex Vip = gsl_matrix_complex_get (V, i, k);
            gsl_complex Viq = gsl_matrix_complex_get (V, i, n - 1);

            gsl_complex h1 = gsl_complex_mul_real (Vip, c);
            gsl_complex h2 = gsl_complex_mul_real (Viq, s);
            gsl_matrix_complex_set (V, i, k, gsl_complex_sub (h1, h2));

            h1 = gsl_complex_mul_real (Vip, s);
            h2 = gsl_complex_mul_real (Viq, c);
            gsl_matrix_complex_set (V, i, n - 1, gsl_complex_add (h1, h2));
          }
      }
// #endif

      /* compute B <= B G */

      gsl_vector_set (d, k, c * x - s * y);

      if (k == n - 2)
        gsl_vector_set (f, k, s * x + c * y);

      if (k > 0)
        {
          double z = gsl_vector_get (f, k - 1);
          gsl_vector_set (f, k - 1, c * z);

          x = gsl_vector_get (d, k - 1);
          y = s * z;
        }
    }
}

static void
qrstep_complex (gsl_vector *d, gsl_vector *f, gsl_matrix_complex *U,
                gsl_matrix_complex *V)
{
#if !USE_BLAS
  const size_t M = U->size1;
  const size_t N = V->size1;
#endif
  const size_t n = d->size;
  double y, z;
  double ak, bk, zk, ap, bp, aq;
  size_t i, k;

  if (n == 1)
    return;                     /* shouldn't happen */

  /* Compute 2x2 svd directly */

  if (n == 2)
    {
      svd2_complex (d, f, U, V);
      return;
    }

  /* Chase out any zeroes on the diagonal */

  for (i = 0; i < n - 1; i++)
    {
      double d_i = gsl_vector_get (d, i);

      if (d_i == 0.0)
        {
          chase_out_intermediate_zero_complex (d, f, U, i);
          return;
        }
    }

  /* Chase out any zero at the end of the diagonal */

  {
    double d_nm1 = gsl_vector_get (d, n - 1);

    if (d_nm1 == 0.0)
      {
        chase_out_trailing_zero_complex (d, f, V);
        return;
      }
  }


  /* Apply QR reduction steps to the diagonal and offdiagonal */

  {
    double d0 = gsl_vector_get (d, 0);
    double f0 = gsl_vector_get (f, 0);

    double d1 = gsl_vector_get (d, 1);

    {
      double mu = trailing_eigenvalue (d, f);

      y = d0 * d0 - mu;
      z = d0 * f0;
    }

    /* Set up the recurrence for Givens rotations on a bidiagonal matrix */

    ak = 0;
    bk = 0;

    ap = d0;
    bp = f0;

    aq = d1;
  }

  for (k = 0; k < n - 1; k++)
    {
      double c, s;
      gsl_linalg_givens (y, z, &c, &s);

      /* Compute V <= V G */

// #ifdef USE_BLAS
//       {
//         gsl_vector_complex_view Vk = gsl_matrix_complex_column (V, k);
//         gsl_vector_complex_view Vkp1 = gsl_matrix_complex_column (V, k + 1);
//         gsl_blas_zdrot (&Vk.vector, &Vkp1.vector, c, -s);
//         // TODO blas_zdrot does not exist yet
//       }
// #else
      for (i = 0; i < N; i++)
        {
          gsl_complex Vip = gsl_matrix_complex_get (V, i, k);
          gsl_complex Viq = gsl_matrix_complex_get (V, i, k + 1);

          gsl_complex h1 = gsl_complex_mul_real (Vip, c);
          gsl_complex h2 = gsl_complex_mul_real (Viq, s);
          gsl_matrix_complex_set (V, i, k, gsl_complex_sub (h1, h2));

          h1 = gsl_complex_mul_real (Vip, s);
          h2 = gsl_complex_mul_real (Viq, c);
          gsl_matrix_complex_set (V, i, k + 1, gsl_complex_add (h1, h2));
        }
// #endif

      /* compute B <= B G */

      {
        double bk1 = c * bk - s * z;

        double ap1 = c * ap - s * bp;
        double bp1 = s * ap + c * bp;
        double zp1 = -s * aq;

        double aq1 = c * aq;

        if (k > 0)
          {
            gsl_vector_set (f, k - 1, bk1);
          }

        ak = ap1;
        bk = bp1;
        zk = zp1;

        ap = aq1;

        if (k < n - 2)
          {
            bp = gsl_vector_get (f, k + 1);
          }
        else
          {
            bp = 0.0;
          }

        y = ak;
        z = zk;
      }

      gsl_linalg_givens (y, z, &c, &s);

      /* Compute U <= U G */

// #ifdef USE_BLAS
//       {
//         gsl_vector_complex_view Uk = gsl_matrix_complex_column (U, k);
//         gsl_vector_complex_view Ukp1 = gsl_matrix_complex_column (U, k + 1);
//         gsl_blas_zdrot (&Uk.vector, &Ukp1.vector, c, -s);
//         // TODO blas_zdrot does not exist yet
//       }
// #else
      for (i = 0; i < M; i++)
        {
          gsl_complex Uip = gsl_matrix_complex_get (U, i, k);
          gsl_complex Uiq = gsl_matrix_complex_get (U, i, k + 1);

          gsl_complex h1 = gsl_complex_mul_real (Uip, c);
          gsl_complex h2 = gsl_complex_mul_real (Uiq, s);
          gsl_matrix_complex_set (U, i, k, gsl_complex_sub (h1, h2));

          h1 = gsl_complex_mul_real (Uip, s);
          h2 = gsl_complex_mul_real (Uiq, c);
          gsl_matrix_complex_set (U, i, k + 1, gsl_complex_add (h1, h2));
        }
// #endif

      /* compute B <= G^T B */

      {
        double ak1 = c * ak - s * zk;
        double bk1 = c * bk - s * ap;
        double zk1 = -s * bp;

        double ap1 = s * bk + c * ap;
        double bp1 = c * bp;

        gsl_vector_set (d, k, ak1);

        ak = ak1;
        bk = bk1;
        zk = zk1;

        ap = ap1;
        bp = bp1;

        if (k < n - 2)
          {
            aq = gsl_vector_get (d, k + 2);
          }
        else
          {
            aq = 0.0;
          }

        y = bk;
        z = zk;
      }
    }

  gsl_vector_set (f, n - 2, bk);
  gsl_vector_set (d, n - 1, ap);
}
