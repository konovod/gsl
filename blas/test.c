/* blas/test.c
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

/* The gsl_blas_* wrappers reject an enumerated argument that is not one
 * of the CBLAS values with GSL_EINVAL, using the same predicates as the
 * CBLAS kernels.  Before the wrappers checked, such a value reached
 * cblas_xerbla() and aborted the process unconditionally.  The error
 * handler is switched off so the return code can be observed.  See
 * Savannah bug #29834. */

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>

#define BAD (99)

static gsl_matrix *ds3;
static gsl_vector *dv3;
static gsl_matrix_float *fs3;
static gsl_vector_float *fv3;
static gsl_matrix_complex *czs3;
static gsl_vector_complex *czv3;
static gsl_matrix_complex_float *cfs3;
static gsl_vector_complex_float *cfv3;

static void
setup (void)
{
  gsl_complex cz = gsl_complex_rect (1.0, 0.0);
  gsl_complex_float cf;
  GSL_SET_COMPLEX (&cf, 1.0f, 0.0f);

  ds3 = gsl_matrix_alloc (3, 3);
  dv3 = gsl_vector_alloc (3);
  fs3 = gsl_matrix_float_alloc (3, 3);
  fv3 = gsl_vector_float_alloc (3);
  czs3 = gsl_matrix_complex_alloc (3, 3);
  czv3 = gsl_vector_complex_alloc (3);
  cfs3 = gsl_matrix_complex_float_alloc (3, 3);
  cfv3 = gsl_vector_complex_float_alloc (3);

  gsl_matrix_set_identity (ds3);
  gsl_vector_set_all (dv3, 1.0);
  gsl_matrix_float_set_identity (fs3);
  gsl_vector_float_set_all (fv3, 1.0f);
  gsl_matrix_complex_set_identity (czs3);
  gsl_vector_complex_set_all (czv3, cz);
  gsl_matrix_complex_float_set_identity (cfs3);
  gsl_vector_complex_float_set_all (cfv3, cf);
}

int
main (void)
{
  gsl_complex cz_one = gsl_complex_rect (1.0, 0.0);
  gsl_complex_float cf_one;
  GSL_SET_COMPLEX (&cf_one, 1.0f, 0.0f);

  setup ();
  gsl_set_error_handler_off ();

  gsl_test (gsl_blas_sgemv (CblasNoTrans, 1.0f, fs3, fv3, 1.0f, fv3) != GSL_SUCCESS,
            "gsl_blas_sgemv valid arguments");
  gsl_test (gsl_blas_sgemv ((CBLAS_TRANSPOSE_t) BAD, 1.0f, fs3, fv3, 1.0f, fv3) != GSL_EINVAL,
            "gsl_blas_sgemv invalid TransA");

  gsl_test (gsl_blas_dgemv (CblasNoTrans, 1.0, ds3, dv3, 1.0, dv3) != GSL_SUCCESS,
            "gsl_blas_dgemv valid arguments");
  gsl_test (gsl_blas_dgemv ((CBLAS_TRANSPOSE_t) BAD, 1.0, ds3, dv3, 1.0, dv3) != GSL_EINVAL,
            "gsl_blas_dgemv invalid TransA");

  gsl_test (gsl_blas_cgemv (CblasNoTrans, cf_one, cfs3, cfv3, cf_one, cfv3) != GSL_SUCCESS,
            "gsl_blas_cgemv valid arguments");
  gsl_test (gsl_blas_cgemv ((CBLAS_TRANSPOSE_t) BAD, cf_one, cfs3, cfv3, cf_one, cfv3) != GSL_EINVAL,
            "gsl_blas_cgemv invalid TransA");

  gsl_test (gsl_blas_zgemv (CblasNoTrans, cz_one, czs3, czv3, cz_one, czv3) != GSL_SUCCESS,
            "gsl_blas_zgemv valid arguments");
  gsl_test (gsl_blas_zgemv ((CBLAS_TRANSPOSE_t) BAD, cz_one, czs3, czv3, cz_one, czv3) != GSL_EINVAL,
            "gsl_blas_zgemv invalid TransA");

  gsl_test (gsl_blas_chemv (CblasUpper, cf_one, cfs3, cfv3, cf_one, cfv3) != GSL_SUCCESS,
            "gsl_blas_chemv valid arguments");
  gsl_test (gsl_blas_chemv ((CBLAS_UPLO_t) BAD, cf_one, cfs3, cfv3, cf_one, cfv3) != GSL_EINVAL,
            "gsl_blas_chemv invalid Uplo");

  gsl_test (gsl_blas_zhemv (CblasUpper, cz_one, czs3, czv3, cz_one, czv3) != GSL_SUCCESS,
            "gsl_blas_zhemv valid arguments");
  gsl_test (gsl_blas_zhemv ((CBLAS_UPLO_t) BAD, cz_one, czs3, czv3, cz_one, czv3) != GSL_EINVAL,
            "gsl_blas_zhemv invalid Uplo");

  gsl_test (gsl_blas_ssymv (CblasUpper, 1.0f, fs3, fv3, 1.0f, fv3) != GSL_SUCCESS,
            "gsl_blas_ssymv valid arguments");
  gsl_test (gsl_blas_ssymv ((CBLAS_UPLO_t) BAD, 1.0f, fs3, fv3, 1.0f, fv3) != GSL_EINVAL,
            "gsl_blas_ssymv invalid Uplo");

  gsl_test (gsl_blas_dsymv (CblasUpper, 1.0, ds3, dv3, 1.0, dv3) != GSL_SUCCESS,
            "gsl_blas_dsymv valid arguments");
  gsl_test (gsl_blas_dsymv ((CBLAS_UPLO_t) BAD, 1.0, ds3, dv3, 1.0, dv3) != GSL_EINVAL,
            "gsl_blas_dsymv invalid Uplo");

  gsl_test (gsl_blas_strmv (CblasUpper, CblasNoTrans, CblasNonUnit, fs3, fv3) != GSL_SUCCESS,
            "gsl_blas_strmv valid arguments");
  gsl_test (gsl_blas_strmv ((CBLAS_UPLO_t) BAD, CblasNoTrans, CblasNonUnit, fs3, fv3) != GSL_EINVAL,
            "gsl_blas_strmv invalid Uplo");
  gsl_test (gsl_blas_strmv (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, CblasNonUnit, fs3, fv3) != GSL_EINVAL,
            "gsl_blas_strmv invalid TransA");
  gsl_test (gsl_blas_strmv (CblasUpper, CblasNoTrans, (CBLAS_DIAG_t) BAD, fs3, fv3) != GSL_EINVAL,
            "gsl_blas_strmv invalid Diag");

  gsl_test (gsl_blas_dtrmv (CblasUpper, CblasNoTrans, CblasNonUnit, ds3, dv3) != GSL_SUCCESS,
            "gsl_blas_dtrmv valid arguments");
  gsl_test (gsl_blas_dtrmv ((CBLAS_UPLO_t) BAD, CblasNoTrans, CblasNonUnit, ds3, dv3) != GSL_EINVAL,
            "gsl_blas_dtrmv invalid Uplo");
  gsl_test (gsl_blas_dtrmv (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, CblasNonUnit, ds3, dv3) != GSL_EINVAL,
            "gsl_blas_dtrmv invalid TransA");
  gsl_test (gsl_blas_dtrmv (CblasUpper, CblasNoTrans, (CBLAS_DIAG_t) BAD, ds3, dv3) != GSL_EINVAL,
            "gsl_blas_dtrmv invalid Diag");

  gsl_test (gsl_blas_ctrmv (CblasUpper, CblasNoTrans, CblasNonUnit, cfs3, cfv3) != GSL_SUCCESS,
            "gsl_blas_ctrmv valid arguments");
  gsl_test (gsl_blas_ctrmv ((CBLAS_UPLO_t) BAD, CblasNoTrans, CblasNonUnit, cfs3, cfv3) != GSL_EINVAL,
            "gsl_blas_ctrmv invalid Uplo");
  gsl_test (gsl_blas_ctrmv (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, CblasNonUnit, cfs3, cfv3) != GSL_EINVAL,
            "gsl_blas_ctrmv invalid TransA");
  gsl_test (gsl_blas_ctrmv (CblasUpper, CblasNoTrans, (CBLAS_DIAG_t) BAD, cfs3, cfv3) != GSL_EINVAL,
            "gsl_blas_ctrmv invalid Diag");

  gsl_test (gsl_blas_ztrmv (CblasUpper, CblasNoTrans, CblasNonUnit, czs3, czv3) != GSL_SUCCESS,
            "gsl_blas_ztrmv valid arguments");
  gsl_test (gsl_blas_ztrmv ((CBLAS_UPLO_t) BAD, CblasNoTrans, CblasNonUnit, czs3, czv3) != GSL_EINVAL,
            "gsl_blas_ztrmv invalid Uplo");
  gsl_test (gsl_blas_ztrmv (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, CblasNonUnit, czs3, czv3) != GSL_EINVAL,
            "gsl_blas_ztrmv invalid TransA");
  gsl_test (gsl_blas_ztrmv (CblasUpper, CblasNoTrans, (CBLAS_DIAG_t) BAD, czs3, czv3) != GSL_EINVAL,
            "gsl_blas_ztrmv invalid Diag");

  gsl_test (gsl_blas_strsv (CblasUpper, CblasNoTrans, CblasNonUnit, fs3, fv3) != GSL_SUCCESS,
            "gsl_blas_strsv valid arguments");
  gsl_test (gsl_blas_strsv ((CBLAS_UPLO_t) BAD, CblasNoTrans, CblasNonUnit, fs3, fv3) != GSL_EINVAL,
            "gsl_blas_strsv invalid Uplo");
  gsl_test (gsl_blas_strsv (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, CblasNonUnit, fs3, fv3) != GSL_EINVAL,
            "gsl_blas_strsv invalid TransA");
  gsl_test (gsl_blas_strsv (CblasUpper, CblasNoTrans, (CBLAS_DIAG_t) BAD, fs3, fv3) != GSL_EINVAL,
            "gsl_blas_strsv invalid Diag");

  gsl_test (gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, ds3, dv3) != GSL_SUCCESS,
            "gsl_blas_dtrsv valid arguments");
  gsl_test (gsl_blas_dtrsv ((CBLAS_UPLO_t) BAD, CblasNoTrans, CblasNonUnit, ds3, dv3) != GSL_EINVAL,
            "gsl_blas_dtrsv invalid Uplo");
  gsl_test (gsl_blas_dtrsv (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, CblasNonUnit, ds3, dv3) != GSL_EINVAL,
            "gsl_blas_dtrsv invalid TransA");
  gsl_test (gsl_blas_dtrsv (CblasUpper, CblasNoTrans, (CBLAS_DIAG_t) BAD, ds3, dv3) != GSL_EINVAL,
            "gsl_blas_dtrsv invalid Diag");

  gsl_test (gsl_blas_ctrsv (CblasUpper, CblasNoTrans, CblasNonUnit, cfs3, cfv3) != GSL_SUCCESS,
            "gsl_blas_ctrsv valid arguments");
  gsl_test (gsl_blas_ctrsv ((CBLAS_UPLO_t) BAD, CblasNoTrans, CblasNonUnit, cfs3, cfv3) != GSL_EINVAL,
            "gsl_blas_ctrsv invalid Uplo");
  gsl_test (gsl_blas_ctrsv (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, CblasNonUnit, cfs3, cfv3) != GSL_EINVAL,
            "gsl_blas_ctrsv invalid TransA");
  gsl_test (gsl_blas_ctrsv (CblasUpper, CblasNoTrans, (CBLAS_DIAG_t) BAD, cfs3, cfv3) != GSL_EINVAL,
            "gsl_blas_ctrsv invalid Diag");

  gsl_test (gsl_blas_ztrsv (CblasUpper, CblasNoTrans, CblasNonUnit, czs3, czv3) != GSL_SUCCESS,
            "gsl_blas_ztrsv valid arguments");
  gsl_test (gsl_blas_ztrsv ((CBLAS_UPLO_t) BAD, CblasNoTrans, CblasNonUnit, czs3, czv3) != GSL_EINVAL,
            "gsl_blas_ztrsv invalid Uplo");
  gsl_test (gsl_blas_ztrsv (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, CblasNonUnit, czs3, czv3) != GSL_EINVAL,
            "gsl_blas_ztrsv invalid TransA");
  gsl_test (gsl_blas_ztrsv (CblasUpper, CblasNoTrans, (CBLAS_DIAG_t) BAD, czs3, czv3) != GSL_EINVAL,
            "gsl_blas_ztrsv invalid Diag");

  gsl_test (gsl_blas_cher (CblasUpper, 1.0f, cfv3, cfs3) != GSL_SUCCESS,
            "gsl_blas_cher valid arguments");
  gsl_test (gsl_blas_cher ((CBLAS_UPLO_t) BAD, 1.0f, cfv3, cfs3) != GSL_EINVAL,
            "gsl_blas_cher invalid Uplo");

  gsl_test (gsl_blas_zher (CblasUpper, 1.0, czv3, czs3) != GSL_SUCCESS,
            "gsl_blas_zher valid arguments");
  gsl_test (gsl_blas_zher ((CBLAS_UPLO_t) BAD, 1.0, czv3, czs3) != GSL_EINVAL,
            "gsl_blas_zher invalid Uplo");

  gsl_test (gsl_blas_cher2 (CblasUpper, cf_one, cfv3, cfv3, cfs3) != GSL_SUCCESS,
            "gsl_blas_cher2 valid arguments");
  gsl_test (gsl_blas_cher2 ((CBLAS_UPLO_t) BAD, cf_one, cfv3, cfv3, cfs3) != GSL_EINVAL,
            "gsl_blas_cher2 invalid Uplo");

  gsl_test (gsl_blas_zher2 (CblasUpper, cz_one, czv3, czv3, czs3) != GSL_SUCCESS,
            "gsl_blas_zher2 valid arguments");
  gsl_test (gsl_blas_zher2 ((CBLAS_UPLO_t) BAD, cz_one, czv3, czv3, czs3) != GSL_EINVAL,
            "gsl_blas_zher2 invalid Uplo");

  gsl_test (gsl_blas_ssyr (CblasUpper, 1.0f, fv3, fs3) != GSL_SUCCESS,
            "gsl_blas_ssyr valid arguments");
  gsl_test (gsl_blas_ssyr ((CBLAS_UPLO_t) BAD, 1.0f, fv3, fs3) != GSL_EINVAL,
            "gsl_blas_ssyr invalid Uplo");

  gsl_test (gsl_blas_dsyr (CblasUpper, 1.0, dv3, ds3) != GSL_SUCCESS,
            "gsl_blas_dsyr valid arguments");
  gsl_test (gsl_blas_dsyr ((CBLAS_UPLO_t) BAD, 1.0, dv3, ds3) != GSL_EINVAL,
            "gsl_blas_dsyr invalid Uplo");

  gsl_test (gsl_blas_ssyr2 (CblasUpper, 1.0f, fv3, fv3, fs3) != GSL_SUCCESS,
            "gsl_blas_ssyr2 valid arguments");
  gsl_test (gsl_blas_ssyr2 ((CBLAS_UPLO_t) BAD, 1.0f, fv3, fv3, fs3) != GSL_EINVAL,
            "gsl_blas_ssyr2 invalid Uplo");

  gsl_test (gsl_blas_dsyr2 (CblasUpper, 1.0, dv3, dv3, ds3) != GSL_SUCCESS,
            "gsl_blas_dsyr2 valid arguments");
  gsl_test (gsl_blas_dsyr2 ((CBLAS_UPLO_t) BAD, 1.0, dv3, dv3, ds3) != GSL_EINVAL,
            "gsl_blas_dsyr2 invalid Uplo");

  gsl_test (gsl_blas_sgemm (CblasNoTrans, CblasNoTrans, 1.0f, fs3, fs3, 1.0f, fs3) != GSL_SUCCESS,
            "gsl_blas_sgemm valid arguments");
  gsl_test (gsl_blas_sgemm ((CBLAS_TRANSPOSE_t) BAD, CblasNoTrans, 1.0f, fs3, fs3, 1.0f, fs3) != GSL_EINVAL,
            "gsl_blas_sgemm invalid TransA");
  gsl_test (gsl_blas_sgemm (CblasNoTrans, (CBLAS_TRANSPOSE_t) BAD, 1.0f, fs3, fs3, 1.0f, fs3) != GSL_EINVAL,
            "gsl_blas_sgemm invalid TransB");

  gsl_test (gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, ds3, ds3, 1.0, ds3) != GSL_SUCCESS,
            "gsl_blas_dgemm valid arguments");
  gsl_test (gsl_blas_dgemm ((CBLAS_TRANSPOSE_t) BAD, CblasNoTrans, 1.0, ds3, ds3, 1.0, ds3) != GSL_EINVAL,
            "gsl_blas_dgemm invalid TransA");
  gsl_test (gsl_blas_dgemm (CblasNoTrans, (CBLAS_TRANSPOSE_t) BAD, 1.0, ds3, ds3, 1.0, ds3) != GSL_EINVAL,
            "gsl_blas_dgemm invalid TransB");

  gsl_test (gsl_blas_cgemm (CblasNoTrans, CblasNoTrans, cf_one, cfs3, cfs3, cf_one, cfs3) != GSL_SUCCESS,
            "gsl_blas_cgemm valid arguments");
  gsl_test (gsl_blas_cgemm ((CBLAS_TRANSPOSE_t) BAD, CblasNoTrans, cf_one, cfs3, cfs3, cf_one, cfs3) != GSL_EINVAL,
            "gsl_blas_cgemm invalid TransA");
  gsl_test (gsl_blas_cgemm (CblasNoTrans, (CBLAS_TRANSPOSE_t) BAD, cf_one, cfs3, cfs3, cf_one, cfs3) != GSL_EINVAL,
            "gsl_blas_cgemm invalid TransB");

  gsl_test (gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, cz_one, czs3, czs3, cz_one, czs3) != GSL_SUCCESS,
            "gsl_blas_zgemm valid arguments");
  gsl_test (gsl_blas_zgemm ((CBLAS_TRANSPOSE_t) BAD, CblasNoTrans, cz_one, czs3, czs3, cz_one, czs3) != GSL_EINVAL,
            "gsl_blas_zgemm invalid TransA");
  gsl_test (gsl_blas_zgemm (CblasNoTrans, (CBLAS_TRANSPOSE_t) BAD, cz_one, czs3, czs3, cz_one, czs3) != GSL_EINVAL,
            "gsl_blas_zgemm invalid TransB");

  gsl_test (gsl_blas_ssymm (CblasLeft, CblasUpper, 1.0f, fs3, fs3, 1.0f, fs3) != GSL_SUCCESS,
            "gsl_blas_ssymm valid arguments");
  gsl_test (gsl_blas_ssymm ((CBLAS_SIDE_t) BAD, CblasUpper, 1.0f, fs3, fs3, 1.0f, fs3) != GSL_EINVAL,
            "gsl_blas_ssymm invalid Side");
  gsl_test (gsl_blas_ssymm (CblasLeft, (CBLAS_UPLO_t) BAD, 1.0f, fs3, fs3, 1.0f, fs3) != GSL_EINVAL,
            "gsl_blas_ssymm invalid Uplo");

  gsl_test (gsl_blas_dsymm (CblasLeft, CblasUpper, 1.0, ds3, ds3, 1.0, ds3) != GSL_SUCCESS,
            "gsl_blas_dsymm valid arguments");
  gsl_test (gsl_blas_dsymm ((CBLAS_SIDE_t) BAD, CblasUpper, 1.0, ds3, ds3, 1.0, ds3) != GSL_EINVAL,
            "gsl_blas_dsymm invalid Side");
  gsl_test (gsl_blas_dsymm (CblasLeft, (CBLAS_UPLO_t) BAD, 1.0, ds3, ds3, 1.0, ds3) != GSL_EINVAL,
            "gsl_blas_dsymm invalid Uplo");

  gsl_test (gsl_blas_csymm (CblasLeft, CblasUpper, cf_one, cfs3, cfs3, cf_one, cfs3) != GSL_SUCCESS,
            "gsl_blas_csymm valid arguments");
  gsl_test (gsl_blas_csymm ((CBLAS_SIDE_t) BAD, CblasUpper, cf_one, cfs3, cfs3, cf_one, cfs3) != GSL_EINVAL,
            "gsl_blas_csymm invalid Side");
  gsl_test (gsl_blas_csymm (CblasLeft, (CBLAS_UPLO_t) BAD, cf_one, cfs3, cfs3, cf_one, cfs3) != GSL_EINVAL,
            "gsl_blas_csymm invalid Uplo");

  gsl_test (gsl_blas_zsymm (CblasLeft, CblasUpper, cz_one, czs3, czs3, cz_one, czs3) != GSL_SUCCESS,
            "gsl_blas_zsymm valid arguments");
  gsl_test (gsl_blas_zsymm ((CBLAS_SIDE_t) BAD, CblasUpper, cz_one, czs3, czs3, cz_one, czs3) != GSL_EINVAL,
            "gsl_blas_zsymm invalid Side");
  gsl_test (gsl_blas_zsymm (CblasLeft, (CBLAS_UPLO_t) BAD, cz_one, czs3, czs3, cz_one, czs3) != GSL_EINVAL,
            "gsl_blas_zsymm invalid Uplo");

  gsl_test (gsl_blas_chemm (CblasLeft, CblasUpper, cf_one, cfs3, cfs3, cf_one, cfs3) != GSL_SUCCESS,
            "gsl_blas_chemm valid arguments");
  gsl_test (gsl_blas_chemm ((CBLAS_SIDE_t) BAD, CblasUpper, cf_one, cfs3, cfs3, cf_one, cfs3) != GSL_EINVAL,
            "gsl_blas_chemm invalid Side");
  gsl_test (gsl_blas_chemm (CblasLeft, (CBLAS_UPLO_t) BAD, cf_one, cfs3, cfs3, cf_one, cfs3) != GSL_EINVAL,
            "gsl_blas_chemm invalid Uplo");

  gsl_test (gsl_blas_zhemm (CblasLeft, CblasUpper, cz_one, czs3, czs3, cz_one, czs3) != GSL_SUCCESS,
            "gsl_blas_zhemm valid arguments");
  gsl_test (gsl_blas_zhemm ((CBLAS_SIDE_t) BAD, CblasUpper, cz_one, czs3, czs3, cz_one, czs3) != GSL_EINVAL,
            "gsl_blas_zhemm invalid Side");
  gsl_test (gsl_blas_zhemm (CblasLeft, (CBLAS_UPLO_t) BAD, cz_one, czs3, czs3, cz_one, czs3) != GSL_EINVAL,
            "gsl_blas_zhemm invalid Uplo");

  gsl_test (gsl_blas_ssyrk (CblasUpper, CblasNoTrans, 1.0f, fs3, 1.0f, fs3) != GSL_SUCCESS,
            "gsl_blas_ssyrk valid arguments");
  gsl_test (gsl_blas_ssyrk ((CBLAS_UPLO_t) BAD, CblasNoTrans, 1.0f, fs3, 1.0f, fs3) != GSL_EINVAL,
            "gsl_blas_ssyrk invalid Uplo");
  gsl_test (gsl_blas_ssyrk (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, 1.0f, fs3, 1.0f, fs3) != GSL_EINVAL,
            "gsl_blas_ssyrk invalid Trans");

  gsl_test (gsl_blas_dsyrk (CblasUpper, CblasNoTrans, 1.0, ds3, 1.0, ds3) != GSL_SUCCESS,
            "gsl_blas_dsyrk valid arguments");
  gsl_test (gsl_blas_dsyrk ((CBLAS_UPLO_t) BAD, CblasNoTrans, 1.0, ds3, 1.0, ds3) != GSL_EINVAL,
            "gsl_blas_dsyrk invalid Uplo");
  gsl_test (gsl_blas_dsyrk (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, 1.0, ds3, 1.0, ds3) != GSL_EINVAL,
            "gsl_blas_dsyrk invalid Trans");

  gsl_test (gsl_blas_csyrk (CblasUpper, CblasNoTrans, cf_one, cfs3, cf_one, cfs3) != GSL_SUCCESS,
            "gsl_blas_csyrk valid arguments");
  gsl_test (gsl_blas_csyrk ((CBLAS_UPLO_t) BAD, CblasNoTrans, cf_one, cfs3, cf_one, cfs3) != GSL_EINVAL,
            "gsl_blas_csyrk invalid Uplo");
  gsl_test (gsl_blas_csyrk (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, cf_one, cfs3, cf_one, cfs3) != GSL_EINVAL,
            "gsl_blas_csyrk invalid Trans");

  gsl_test (gsl_blas_zsyrk (CblasUpper, CblasNoTrans, cz_one, czs3, cz_one, czs3) != GSL_SUCCESS,
            "gsl_blas_zsyrk valid arguments");
  gsl_test (gsl_blas_zsyrk ((CBLAS_UPLO_t) BAD, CblasNoTrans, cz_one, czs3, cz_one, czs3) != GSL_EINVAL,
            "gsl_blas_zsyrk invalid Uplo");
  gsl_test (gsl_blas_zsyrk (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, cz_one, czs3, cz_one, czs3) != GSL_EINVAL,
            "gsl_blas_zsyrk invalid Trans");

  gsl_test (gsl_blas_cherk (CblasUpper, CblasNoTrans, 1.0f, cfs3, 1.0f, cfs3) != GSL_SUCCESS,
            "gsl_blas_cherk valid arguments");
  gsl_test (gsl_blas_cherk ((CBLAS_UPLO_t) BAD, CblasNoTrans, 1.0f, cfs3, 1.0f, cfs3) != GSL_EINVAL,
            "gsl_blas_cherk invalid Uplo");
  gsl_test (gsl_blas_cherk (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, 1.0f, cfs3, 1.0f, cfs3) != GSL_EINVAL,
            "gsl_blas_cherk invalid Trans");

  gsl_test (gsl_blas_zherk (CblasUpper, CblasNoTrans, 1.0, czs3, 1.0, czs3) != GSL_SUCCESS,
            "gsl_blas_zherk valid arguments");
  gsl_test (gsl_blas_zherk ((CBLAS_UPLO_t) BAD, CblasNoTrans, 1.0, czs3, 1.0, czs3) != GSL_EINVAL,
            "gsl_blas_zherk invalid Uplo");
  gsl_test (gsl_blas_zherk (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, 1.0, czs3, 1.0, czs3) != GSL_EINVAL,
            "gsl_blas_zherk invalid Trans");

  gsl_test (gsl_blas_ssyr2k (CblasUpper, CblasNoTrans, 1.0f, fs3, fs3, 1.0f, fs3) != GSL_SUCCESS,
            "gsl_blas_ssyr2k valid arguments");
  gsl_test (gsl_blas_ssyr2k ((CBLAS_UPLO_t) BAD, CblasNoTrans, 1.0f, fs3, fs3, 1.0f, fs3) != GSL_EINVAL,
            "gsl_blas_ssyr2k invalid Uplo");
  gsl_test (gsl_blas_ssyr2k (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, 1.0f, fs3, fs3, 1.0f, fs3) != GSL_EINVAL,
            "gsl_blas_ssyr2k invalid Trans");

  gsl_test (gsl_blas_dsyr2k (CblasUpper, CblasNoTrans, 1.0, ds3, ds3, 1.0, ds3) != GSL_SUCCESS,
            "gsl_blas_dsyr2k valid arguments");
  gsl_test (gsl_blas_dsyr2k ((CBLAS_UPLO_t) BAD, CblasNoTrans, 1.0, ds3, ds3, 1.0, ds3) != GSL_EINVAL,
            "gsl_blas_dsyr2k invalid Uplo");
  gsl_test (gsl_blas_dsyr2k (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, 1.0, ds3, ds3, 1.0, ds3) != GSL_EINVAL,
            "gsl_blas_dsyr2k invalid Trans");

  gsl_test (gsl_blas_csyr2k (CblasUpper, CblasNoTrans, cf_one, cfs3, cfs3, cf_one, cfs3) != GSL_SUCCESS,
            "gsl_blas_csyr2k valid arguments");
  gsl_test (gsl_blas_csyr2k ((CBLAS_UPLO_t) BAD, CblasNoTrans, cf_one, cfs3, cfs3, cf_one, cfs3) != GSL_EINVAL,
            "gsl_blas_csyr2k invalid Uplo");
  gsl_test (gsl_blas_csyr2k (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, cf_one, cfs3, cfs3, cf_one, cfs3) != GSL_EINVAL,
            "gsl_blas_csyr2k invalid Trans");

  gsl_test (gsl_blas_zsyr2k (CblasUpper, CblasNoTrans, cz_one, czs3, czs3, cz_one, czs3) != GSL_SUCCESS,
            "gsl_blas_zsyr2k valid arguments");
  gsl_test (gsl_blas_zsyr2k ((CBLAS_UPLO_t) BAD, CblasNoTrans, cz_one, czs3, czs3, cz_one, czs3) != GSL_EINVAL,
            "gsl_blas_zsyr2k invalid Uplo");
  gsl_test (gsl_blas_zsyr2k (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, cz_one, czs3, czs3, cz_one, czs3) != GSL_EINVAL,
            "gsl_blas_zsyr2k invalid Trans");

  gsl_test (gsl_blas_cher2k (CblasUpper, CblasNoTrans, cf_one, cfs3, cfs3, 1.0f, cfs3) != GSL_SUCCESS,
            "gsl_blas_cher2k valid arguments");
  gsl_test (gsl_blas_cher2k ((CBLAS_UPLO_t) BAD, CblasNoTrans, cf_one, cfs3, cfs3, 1.0f, cfs3) != GSL_EINVAL,
            "gsl_blas_cher2k invalid Uplo");
  gsl_test (gsl_blas_cher2k (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, cf_one, cfs3, cfs3, 1.0f, cfs3) != GSL_EINVAL,
            "gsl_blas_cher2k invalid Trans");

  gsl_test (gsl_blas_zher2k (CblasUpper, CblasNoTrans, cz_one, czs3, czs3, 1.0, czs3) != GSL_SUCCESS,
            "gsl_blas_zher2k valid arguments");
  gsl_test (gsl_blas_zher2k ((CBLAS_UPLO_t) BAD, CblasNoTrans, cz_one, czs3, czs3, 1.0, czs3) != GSL_EINVAL,
            "gsl_blas_zher2k invalid Uplo");
  gsl_test (gsl_blas_zher2k (CblasUpper, (CBLAS_TRANSPOSE_t) BAD, cz_one, czs3, czs3, 1.0, czs3) != GSL_EINVAL,
            "gsl_blas_zher2k invalid Trans");

  gsl_test (gsl_blas_strmm (CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0f, fs3, fs3) != GSL_SUCCESS,
            "gsl_blas_strmm valid arguments");
  gsl_test (gsl_blas_strmm ((CBLAS_SIDE_t) BAD, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0f, fs3, fs3) != GSL_EINVAL,
            "gsl_blas_strmm invalid Side");
  gsl_test (gsl_blas_strmm (CblasLeft, (CBLAS_UPLO_t) BAD, CblasNoTrans, CblasNonUnit, 1.0f, fs3, fs3) != GSL_EINVAL,
            "gsl_blas_strmm invalid Uplo");
  gsl_test (gsl_blas_strmm (CblasLeft, CblasUpper, (CBLAS_TRANSPOSE_t) BAD, CblasNonUnit, 1.0f, fs3, fs3) != GSL_EINVAL,
            "gsl_blas_strmm invalid TransA");
  gsl_test (gsl_blas_strmm (CblasLeft, CblasUpper, CblasNoTrans, (CBLAS_DIAG_t) BAD, 1.0f, fs3, fs3) != GSL_EINVAL,
            "gsl_blas_strmm invalid Diag");

  gsl_test (gsl_blas_dtrmm (CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, ds3, ds3) != GSL_SUCCESS,
            "gsl_blas_dtrmm valid arguments");
  gsl_test (gsl_blas_dtrmm ((CBLAS_SIDE_t) BAD, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, ds3, ds3) != GSL_EINVAL,
            "gsl_blas_dtrmm invalid Side");
  gsl_test (gsl_blas_dtrmm (CblasLeft, (CBLAS_UPLO_t) BAD, CblasNoTrans, CblasNonUnit, 1.0, ds3, ds3) != GSL_EINVAL,
            "gsl_blas_dtrmm invalid Uplo");
  gsl_test (gsl_blas_dtrmm (CblasLeft, CblasUpper, (CBLAS_TRANSPOSE_t) BAD, CblasNonUnit, 1.0, ds3, ds3) != GSL_EINVAL,
            "gsl_blas_dtrmm invalid TransA");
  gsl_test (gsl_blas_dtrmm (CblasLeft, CblasUpper, CblasNoTrans, (CBLAS_DIAG_t) BAD, 1.0, ds3, ds3) != GSL_EINVAL,
            "gsl_blas_dtrmm invalid Diag");

  gsl_test (gsl_blas_ctrmm (CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, cf_one, cfs3, cfs3) != GSL_SUCCESS,
            "gsl_blas_ctrmm valid arguments");
  gsl_test (gsl_blas_ctrmm ((CBLAS_SIDE_t) BAD, CblasUpper, CblasNoTrans, CblasNonUnit, cf_one, cfs3, cfs3) != GSL_EINVAL,
            "gsl_blas_ctrmm invalid Side");
  gsl_test (gsl_blas_ctrmm (CblasLeft, (CBLAS_UPLO_t) BAD, CblasNoTrans, CblasNonUnit, cf_one, cfs3, cfs3) != GSL_EINVAL,
            "gsl_blas_ctrmm invalid Uplo");
  gsl_test (gsl_blas_ctrmm (CblasLeft, CblasUpper, (CBLAS_TRANSPOSE_t) BAD, CblasNonUnit, cf_one, cfs3, cfs3) != GSL_EINVAL,
            "gsl_blas_ctrmm invalid TransA");
  gsl_test (gsl_blas_ctrmm (CblasLeft, CblasUpper, CblasNoTrans, (CBLAS_DIAG_t) BAD, cf_one, cfs3, cfs3) != GSL_EINVAL,
            "gsl_blas_ctrmm invalid Diag");

  gsl_test (gsl_blas_ztrmm (CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, cz_one, czs3, czs3) != GSL_SUCCESS,
            "gsl_blas_ztrmm valid arguments");
  gsl_test (gsl_blas_ztrmm ((CBLAS_SIDE_t) BAD, CblasUpper, CblasNoTrans, CblasNonUnit, cz_one, czs3, czs3) != GSL_EINVAL,
            "gsl_blas_ztrmm invalid Side");
  gsl_test (gsl_blas_ztrmm (CblasLeft, (CBLAS_UPLO_t) BAD, CblasNoTrans, CblasNonUnit, cz_one, czs3, czs3) != GSL_EINVAL,
            "gsl_blas_ztrmm invalid Uplo");
  gsl_test (gsl_blas_ztrmm (CblasLeft, CblasUpper, (CBLAS_TRANSPOSE_t) BAD, CblasNonUnit, cz_one, czs3, czs3) != GSL_EINVAL,
            "gsl_blas_ztrmm invalid TransA");
  gsl_test (gsl_blas_ztrmm (CblasLeft, CblasUpper, CblasNoTrans, (CBLAS_DIAG_t) BAD, cz_one, czs3, czs3) != GSL_EINVAL,
            "gsl_blas_ztrmm invalid Diag");

  gsl_test (gsl_blas_strsm (CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0f, fs3, fs3) != GSL_SUCCESS,
            "gsl_blas_strsm valid arguments");
  gsl_test (gsl_blas_strsm ((CBLAS_SIDE_t) BAD, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0f, fs3, fs3) != GSL_EINVAL,
            "gsl_blas_strsm invalid Side");
  gsl_test (gsl_blas_strsm (CblasLeft, (CBLAS_UPLO_t) BAD, CblasNoTrans, CblasNonUnit, 1.0f, fs3, fs3) != GSL_EINVAL,
            "gsl_blas_strsm invalid Uplo");
  gsl_test (gsl_blas_strsm (CblasLeft, CblasUpper, (CBLAS_TRANSPOSE_t) BAD, CblasNonUnit, 1.0f, fs3, fs3) != GSL_EINVAL,
            "gsl_blas_strsm invalid TransA");
  gsl_test (gsl_blas_strsm (CblasLeft, CblasUpper, CblasNoTrans, (CBLAS_DIAG_t) BAD, 1.0f, fs3, fs3) != GSL_EINVAL,
            "gsl_blas_strsm invalid Diag");

  gsl_test (gsl_blas_dtrsm (CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, ds3, ds3) != GSL_SUCCESS,
            "gsl_blas_dtrsm valid arguments");
  gsl_test (gsl_blas_dtrsm ((CBLAS_SIDE_t) BAD, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, ds3, ds3) != GSL_EINVAL,
            "gsl_blas_dtrsm invalid Side");
  gsl_test (gsl_blas_dtrsm (CblasLeft, (CBLAS_UPLO_t) BAD, CblasNoTrans, CblasNonUnit, 1.0, ds3, ds3) != GSL_EINVAL,
            "gsl_blas_dtrsm invalid Uplo");
  gsl_test (gsl_blas_dtrsm (CblasLeft, CblasUpper, (CBLAS_TRANSPOSE_t) BAD, CblasNonUnit, 1.0, ds3, ds3) != GSL_EINVAL,
            "gsl_blas_dtrsm invalid TransA");
  gsl_test (gsl_blas_dtrsm (CblasLeft, CblasUpper, CblasNoTrans, (CBLAS_DIAG_t) BAD, 1.0, ds3, ds3) != GSL_EINVAL,
            "gsl_blas_dtrsm invalid Diag");

  gsl_test (gsl_blas_ctrsm (CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, cf_one, cfs3, cfs3) != GSL_SUCCESS,
            "gsl_blas_ctrsm valid arguments");
  gsl_test (gsl_blas_ctrsm ((CBLAS_SIDE_t) BAD, CblasUpper, CblasNoTrans, CblasNonUnit, cf_one, cfs3, cfs3) != GSL_EINVAL,
            "gsl_blas_ctrsm invalid Side");
  gsl_test (gsl_blas_ctrsm (CblasLeft, (CBLAS_UPLO_t) BAD, CblasNoTrans, CblasNonUnit, cf_one, cfs3, cfs3) != GSL_EINVAL,
            "gsl_blas_ctrsm invalid Uplo");
  gsl_test (gsl_blas_ctrsm (CblasLeft, CblasUpper, (CBLAS_TRANSPOSE_t) BAD, CblasNonUnit, cf_one, cfs3, cfs3) != GSL_EINVAL,
            "gsl_blas_ctrsm invalid TransA");
  gsl_test (gsl_blas_ctrsm (CblasLeft, CblasUpper, CblasNoTrans, (CBLAS_DIAG_t) BAD, cf_one, cfs3, cfs3) != GSL_EINVAL,
            "gsl_blas_ctrsm invalid Diag");

  gsl_test (gsl_blas_ztrsm (CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, cz_one, czs3, czs3) != GSL_SUCCESS,
            "gsl_blas_ztrsm valid arguments");
  gsl_test (gsl_blas_ztrsm ((CBLAS_SIDE_t) BAD, CblasUpper, CblasNoTrans, CblasNonUnit, cz_one, czs3, czs3) != GSL_EINVAL,
            "gsl_blas_ztrsm invalid Side");
  gsl_test (gsl_blas_ztrsm (CblasLeft, (CBLAS_UPLO_t) BAD, CblasNoTrans, CblasNonUnit, cz_one, czs3, czs3) != GSL_EINVAL,
            "gsl_blas_ztrsm invalid Uplo");
  gsl_test (gsl_blas_ztrsm (CblasLeft, CblasUpper, (CBLAS_TRANSPOSE_t) BAD, CblasNonUnit, cz_one, czs3, czs3) != GSL_EINVAL,
            "gsl_blas_ztrsm invalid TransA");
  gsl_test (gsl_blas_ztrsm (CblasLeft, CblasUpper, CblasNoTrans, (CBLAS_DIAG_t) BAD, cz_one, czs3, czs3) != GSL_EINVAL,
            "gsl_blas_ztrsm invalid Diag");

  exit (gsl_test_summary ());
}
