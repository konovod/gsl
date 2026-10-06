#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>

int
test_SV_complex_decomp_dim(const gsl_matrix_complex * m, double eps)
{
  int s = 0;
  double di1;
  unsigned long i,j, M = m->size1, N = m->size2;
  unsigned long input_nans = 0;

  gsl_matrix_complex * v  = gsl_matrix_complex_alloc(M,N);
  gsl_matrix_complex * a  = gsl_matrix_complex_alloc(M,N);
  gsl_matrix_complex * q  = gsl_matrix_complex_alloc(N,N);
  gsl_matrix_complex * dqt  = gsl_matrix_complex_alloc(N,N);
  gsl_vector * d  = gsl_vector_alloc(N);
  gsl_vector_complex * w1  = gsl_vector_complex_alloc(N);
  gsl_vector_complex * w2  = gsl_vector_complex_alloc(N-1);

  gsl_matrix_complex_memcpy(v,m);

  /* Check for nans in the input */
  for (i = 0; i<M; i++) {
    for (j = 0; j<N; j++) {
      gsl_complex m_ij = gsl_matrix_complex_get (m, i, j);
      if (gsl_isnan (GSL_REAL(m_ij)) || gsl_isnan (GSL_IMAG(m_ij))) input_nans++;
    }
  }

  s = gsl_linalg_complex_SV_decomp(v, q, d, w1, w2);

  if (s) printf("returned error code %d = %s\n", s, gsl_strerror(s));

  /* Check that singular values are non-negative and in non-decreasing
     order */
  
  di1 = 0.0;

  for (i = 0; i < N; i++)
    {
      double di = gsl_vector_get (d, i);

      if (gsl_isnan (di))
        {
          if (input_nans > 0) 
            continue;  /* skip NaNs if present in input */
          else
            {
              s++;
              printf("bad singular value %lu = %22.18g\n", i, di);
            }
        }

      if (di < 0) {
        s++;
        printf("singular value %lu = %22.18g < 0\n", i, di);
      }

      if(i > 0 && di > di1) {
        s++;
        printf("singular value %lu = %22.18g vs previous %22.18g\n", i, di, di1);
      }

      di1 = di;
    }      
  
  /* Scale dqt = D Q^T */
  
  for (i = 0; i < N ; i++)
    {
      double di = gsl_vector_get (d, i);

      for (j = 0; j < N; j++)
        {
          gsl_complex qji = gsl_matrix_complex_get(q, j, i);
          qji = gsl_complex_mul_real(gsl_complex_conjugate(qji), di);
          gsl_matrix_complex_set (dqt, i, j, qji);
        }
    }
            
  /* compute a = v dqt */
  gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, v, dqt, GSL_COMPLEX_ZERO, a);

  {
    double scale = 0.0;

    for(i=0; i<M; i++)
      for(j=0; j<N; j++)
        {
          double t = gsl_complex_abs(gsl_matrix_complex_get(m, i, j));
          if (t > scale) scale = t;
        }

    for(i=0; i<M; i++) {
      for(j=0; j<N; j++) {
        gsl_complex aij = gsl_matrix_complex_get(a, i, j);
        gsl_complex mij = gsl_matrix_complex_get(m, i, j);
        double dd = gsl_complex_abs(gsl_complex_sub(aij, mij));
        int foo = (dd > eps * (gsl_complex_abs(mij) + scale));
        if(foo) {
          printf("(%3lu,%3lu)[%lu,%lu]: %22.18g%+22.18gi   %22.18g%+22.18gi\n", M, N, i,j,
                 GSL_REAL(aij), GSL_IMAG(aij), GSL_REAL(mij), GSL_IMAG(mij));
        }
        s += foo;
      }
    }
  }gsl_vector_complex_free(w1);
  gsl_vector_complex_free(w2);
  gsl_vector_free(d);
  gsl_matrix_complex_free(v);
  gsl_matrix_complex_free(a);
  gsl_matrix_complex_free(q);
  gsl_matrix_complex_free(dqt);

  return s;
}

int
test_SV_complex_decomp (gsl_rng * r)
{
  int f;
  int s = 0;
  size_t M, N;
  double acc;

  for (M = 1; M <= 50; M++)
    {
      for (N = 1; N <= M; N++)
        {
          gsl_matrix_complex *A = gsl_matrix_complex_alloc (M, N);

          acc = (M <  8) ? 1e4 * GSL_DBL_EPSILON
              : (M < 20) ? 1e5 * GSL_DBL_EPSILON
              : (M < 30) ? 1e7 * GSL_DBL_EPSILON
                         : 1e8 * GSL_DBL_EPSILON;

          create_random_complex_matrix (A, r);
          f = test_SV_complex_decomp_dim (A, acc);
          gsl_test (f, "  svd_complex_decomp");
          s += f;
        }
    }

  return s;
}

/* A fixed 4x3 complex matrix with singular values from numpy.linalg.svd. */

int
test_SV_complex_reference(void)
{
  int s = 0;
  unsigned long i, j;
  const double expected[3] = {
    4.4207273904778104, 3.6334962276217313, 1.7218810943087512
  };
  const double re[4][3] = {
    { 1.0, 3.0, 0.5 },
    { 2.0, 0.0, -1.0 },
    { 0.3, 1.0, 2.0 },
    { 1.0, -0.5, 0.7 }
  };
  const double im[4][3] = {
    { 2.0, -1.0, 0.5 },
    { -1.0, 1.0, 2.0 },
    { 0.2, -1.0, 0.0 },
    { 1.0, 0.5, -0.3 }
  };
  gsl_matrix_complex *m = gsl_matrix_complex_alloc(4, 3);
  gsl_matrix_complex *v = gsl_matrix_complex_alloc(4, 3);
  gsl_matrix_complex *q = gsl_matrix_complex_alloc(3, 3);
  gsl_vector *d = gsl_vector_alloc(3);
  gsl_vector_complex *w1 = gsl_vector_complex_alloc(3);
  gsl_vector_complex *w2 = gsl_vector_complex_alloc(2);

  for (i = 0; i < 4; i++)
    for (j = 0; j < 3; j++)
      gsl_matrix_complex_set(m, i, j, gsl_complex_rect(re[i][j], im[i][j]));

  gsl_matrix_complex_memcpy(v, m);
  s += gsl_linalg_complex_SV_decomp(v, q, d, w1, w2);

  for (i = 0; i < 3; i++)
    {
      double di = gsl_vector_get(d, i);
      gsl_test(fabs(di - expected[i]) > 1.0e-13 * expected[i],
               "  svd_complex reference singular value %lu: %.17g vs %.17g",
               i, di, expected[i]);
    }

  gsl_vector_complex_free(w2);
  gsl_vector_complex_free(w1);
  gsl_vector_free(d);
  gsl_matrix_complex_free(q);
  gsl_matrix_complex_free(v);
  gsl_matrix_complex_free(m);

  return s;
}
