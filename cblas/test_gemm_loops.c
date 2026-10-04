/* cblas/test_gemm_loops.c
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

/* Regression test for Savannah bug #54925.
 *
 * cblas/source_gemm_r.h was reordered to traverse the output row by row
 * so that a row of C stays resident across the whole k loop.  That is a
 * pure change of loop nesting: for each C(i,j) the terms are still added
 * in the same k order, so the result is unchanged.  This test checks the
 * remaining property directly, by comparing cblas_dgemm against a
 * straightforward triple-loop evaluation of the same mapping, over all
 * four TransA/TransB combinations and both storage orders.
 */

#include <config.h>
#include <stddef.h>

#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>

#include "tests.h"

static void
gemm_reference (int order, int transA, int transB, int M, int N, int K,
                double alpha, const double *A, int lda,
                const double *B, int ldb,
                double beta, double *C, int ldc)
{
  int i, j, k;
  int n1, n2, ldf, ldg, TransF, TransG;
  const double *F, *G;

  /* This mirrors the mapping used by source_gemm_r.h. */

  if (order == CblasRowMajor)
    {
      n1 = M;
      n2 = N;
      F = A;
      ldf = lda;
      TransF = (transA == CblasConjTrans) ? CblasTrans : transA;
      G = B;
      ldg = ldb;
      TransG = (transB == CblasConjTrans) ? CblasTrans : transB;
    }
  else
    {
      n1 = N;
      n2 = M;
      F = B;
      ldf = ldb;
      TransF = (transB == CblasConjTrans) ? CblasTrans : transB;
      G = A;
      ldg = lda;
      TransG = (transA == CblasConjTrans) ? CblasTrans : transA;
    }

  for (i = 0; i < n1; i++)
    {
      for (j = 0; j < n2; j++)
        {
          double s = 0.0;

          for (k = 0; k < K; k++)
            {
              double f = (TransF == CblasNoTrans) ? F[ldf * i + k]
                                                  : F[ldf * k + i];
              double g = (TransG == CblasNoTrans) ? G[ldg * k + j]
                                                  : G[ldg * j + k];
              s += f * g;
            }

          C[ldc * i + j] = beta * C[ldc * i + j] + alpha * s;
        }
    }
}

static void
gemm_loops_case (int order, int transA, int transB, int n,
                 double alpha, double beta)
{
  double A[64], B[64], C[64], Cref[64];
  int i;

  for (i = 0; i < n * n; i++)
    {
      A[i] = 0.17 * (i + 1) - 0.53;
      B[i] = 0.11 * (i + 1) - 0.29;
      C[i] = 0.05 * (i + 1) - 0.37;
      Cref[i] = C[i];
    }

  cblas_dgemm (order, transA, transB, n, n, n, alpha, A, n, B, n,
               beta, C, n);
  gemm_reference (order, transA, transB, n, n, n, alpha, A, n, B, n,
                  beta, Cref, n);

  for (i = 0; i < n * n; i++)
    {
      gsl_test_rel (C[i], Cref[i], 1e-12,
                    "dgemm loops(order=%d, transA=%d, transB=%d, n=%d) [%d]",
                    order, transA, transB, n, i);
    }
}

void
test_gemm_loops (void)
{
  const int orders[2] = { CblasRowMajor, CblasColMajor };
  const int trans[3] = { CblasNoTrans, CblasTrans, CblasConjTrans };
  int o, a, b;

  for (o = 0; o < 2; o++)
    {
      for (a = 0; a < 3; a++)
        {
          for (b = 0; b < 3; b++)
            {
              gemm_loops_case (orders[o], trans[a], trans[b], 4, 0.75, -0.25);
              gemm_loops_case (orders[o], trans[a], trans[b], 5, -1.5, 0.0);
              gemm_loops_case (orders[o], trans[a], trans[b], 3, 1.0, 1.0);
            }
        }
    }
}
