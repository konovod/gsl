/* cblas/error_cblas.h
 *
 * Copyright (C) 2010 José Luis García Pallero
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

#ifndef ERROR_CBLAS_H__
#define ERROR_CBLAS_H__


#define CHECK_ARGS_X(FUNCTION,VAR,ARGS) do { int VAR = 0 ;      \
    CBLAS_ERROR_##FUNCTION ARGS ; \
    if (VAR) cblas_xerbla(pos,__FILE__,""); } while (0)

#define CHECK_ARGS7(FUNCTION,A1,A2,A3,A4,A5,A6,A7) \
  CHECK_ARGS_X(FUNCTION,pos,(pos,A1,A2,A3,A4,A5,A6,A7))

#define CHECK_ARGS8(FUNCTION,A1,A2,A3,A4,A5,A6,A7,A8) \
  CHECK_ARGS_X(FUNCTION,pos,(pos,A1,A2,A3,A4,A5,A6,A7,A8))

#define CHECK_ARGS9(FUNCTION,A1,A2,A3,A4,A5,A6,A7,A8,A9) \
  CHECK_ARGS_X(FUNCTION,pos,(pos,A1,A2,A3,A4,A5,A6,A7,A8,A9))

#define CHECK_ARGS10(FUNCTION,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10) \
  CHECK_ARGS_X(FUNCTION,pos,(pos,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10))

#define CHECK_ARGS11(FUNCTION,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11) \
  CHECK_ARGS_X(FUNCTION,pos,(pos,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11))

#define CHECK_ARGS12(FUNCTION,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12) \
  CHECK_ARGS_X(FUNCTION,pos,(pos,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12))

#define CHECK_ARGS13(FUNCTION,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13) \
  CHECK_ARGS_X(FUNCTION,pos,(pos,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13))

#define CHECK_ARGS14(FUNCTION,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14) \
  CHECK_ARGS_X(FUNCTION,pos,(pos,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14))

/* Predicates for the CBLAS argument types.  These are the single
 * definition of "what counts as a valid argument" and are shared with the
 * GSL BLAS wrappers, which reuse them to validate their enumerated
 * arguments before calling into the CBLAS kernels.  The CHECK_* macros
 * below and the gsl_blas_* wrappers must not disagree. */

/* check if CBLAS_ORDER is correct */
#define CBLAS_ORDER_IS_VALID(order) \
  (((order)==CblasRowMajor)||((order)==CblasColMajor))

/* check if CBLAS_TRANSPOSE is correct */
#define CBLAS_TRANSPOSE_IS_VALID(Trans) \
  (((Trans)==CblasNoTrans)||((Trans)==CblasTrans)||((Trans)==CblasConjTrans))

/* check if CBLAS_UPLO is correct */
#define CBLAS_UPLO_IS_VALID(Uplo) \
  (((Uplo)==CblasUpper)||((Uplo)==CblasLower))

/* check if CBLAS_DIAG is correct */
#define CBLAS_DIAG_IS_VALID(Diag) \
  (((Diag)==CblasNonUnit)||((Diag)==CblasUnit))

/* check if CBLAS_SIDE is correct */
#define CBLAS_SIDE_IS_VALID(Side) \
  (((Side)==CblasLeft)||((Side)==CblasRight))

/* check if a dimension argument is correct */
#define CBLAS_DIM_IS_VALID(dim) ((dim) >= 0)

/* check if a stride argument is correct */
#define CBLAS_STRIDE_IS_VALID(stride) ((stride) != 0)

#define CHECK_ORDER(pos,posIfError,order) \
if(!CBLAS_ORDER_IS_VALID(order)) \
     pos = posIfError;

#define CHECK_TRANSPOSE(pos,posIfError,Trans) \
if(!CBLAS_TRANSPOSE_IS_VALID(Trans)) \
    pos = posIfError;

#define CHECK_UPLO(pos,posIfError,Uplo) \
if(!CBLAS_UPLO_IS_VALID(Uplo)) \
    pos = posIfError;

#define CHECK_DIAG(pos,posIfError,Diag) \
if(!CBLAS_DIAG_IS_VALID(Diag)) \
    pos = posIfError;

#define CHECK_SIDE(pos,posIfError,Side) \
if(!CBLAS_SIDE_IS_VALID(Side)) \
    pos = posIfError;

#define CHECK_DIM(pos,posIfError,dim) \
if(!CBLAS_DIM_IS_VALID(dim)) \
    pos = posIfError;

#define CHECK_STRIDE(pos,posIfError,stride) \
if(!CBLAS_STRIDE_IS_VALID(stride)) \
    pos = posIfError;

#endif /* ERROR_CBLAS_H__ */
