/* -----------------------------------------------------------------
 * Complex Dense SUNMatrix implementation for SUNDIALS.
 *
 * Uses C99 double _Complex as the base element type.
 * Column-major storage, following SUNDIALS SUNMatrix interface.
 * -----------------------------------------------------------------*/

#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include "sunmatrix_complex_dense.h"

/* ================================================================
 * Constructors / Destructors
 * ================================================================*/

SUNMatrix SUNMatNew_ComplexDense(sunindextype M, sunindextype N,
                                 SUNContext sunctx)
{
  SUNMatrix A;
  SUNMatrixContent_ComplexDense content;
  sunindextype j;

  if (M <= 0 || N <= 0) { return NULL; }

  /* Create generic matrix shell */
  A = SUNMatNewEmpty(sunctx);
  if (A == NULL) { return NULL; }

  /* Attach ops */
  A->ops->getid      = SUNMatGetID_ComplexDense;
  A->ops->clone      = SUNMatClone_ComplexDense;
  A->ops->destroy    = SUNMatDestroy_ComplexDense;
  A->ops->zero       = SUNMatZero_ComplexDense;
  A->ops->copy       = SUNMatCopy_ComplexDense;
  A->ops->scaleadd   = SUNMatScaleAdd_ComplexDense;
  A->ops->scaleaddi  = SUNMatScaleAddI_ComplexDense;
  A->ops->matvec     = SUNMatMatvec_ComplexDense;
  A->ops->space      = SUNMatSpace_ComplexDense;
  A->ops->mathermitiantransposevec = SUNMatHermitianTransposeVec_ComplexDense;

  /* Allocate content */
  content = (SUNMatrixContent_ComplexDense)malloc(sizeof *content);
  if (content == NULL) { SUNMatDestroy(A); return NULL; }

  content->M     = M;
  content->N     = N;
  content->ldata = M * N;

  /* Allocate data array */
  content->data = (suncomplextype *)calloc((size_t)(M * N),
                                           sizeof(suncomplextype));
  if (content->data == NULL) { free(content); SUNMatDestroy(A); return NULL; }

  /* Allocate column pointer array */
  content->cols = (suncomplextype **)malloc((size_t)N * sizeof(suncomplextype *));
  if (content->cols == NULL)
  {
    free(content->data);
    free(content);
    SUNMatDestroy(A);
    return NULL;
  }

  /* Set column pointers */
  for (j = 0; j < N; j++) { content->cols[j] = content->data + j * M; }

  A->content = content;

  return A;
}

SUNMatrix SUNMatClone_ComplexDense(SUNMatrix A)
{
  return SUNMatNew_ComplexDense(SM_ROWS_ZD(A), SM_COLUMNS_ZD(A), A->sunctx);
}

void SUNMatDestroy_ComplexDense(SUNMatrix A)
{
  if (A == NULL) { return; }

  if (A->content != NULL)
  {
    SUNMatrixContent_ComplexDense content = SM_CONTENT_ZD(A);
    if (content->data != NULL) { free(content->data); }
    if (content->cols != NULL) { free(content->cols); }
    free(content);
    A->content = NULL;
  }

  if (A->ops != NULL) { free(A->ops); A->ops = NULL; }
  free(A);
}

/* ================================================================
 * SUNMatrix ops implementations
 * ================================================================*/

SUNMatrix_ID SUNMatGetID_ComplexDense(SUNMatrix A)
{
  (void)A;
  return SUNMATRIX_CUSTOM;
}

SUNErrCode SUNMatZero_ComplexDense(SUNMatrix A)
{
  memset(SM_DATA_ZD(A), 0,
         (size_t)SM_LDATA_ZD(A) * sizeof(suncomplextype));
  return SUN_SUCCESS;
}

SUNErrCode SUNMatCopy_ComplexDense(SUNMatrix A, SUNMatrix B)
{
  memcpy(SM_DATA_ZD(B), SM_DATA_ZD(A),
         (size_t)SM_LDATA_ZD(A) * sizeof(suncomplextype));
  return SUN_SUCCESS;
}

/* A = c*A + B (c is real) */
SUNErrCode SUNMatScaleAdd_ComplexDense(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
  sunindextype i, ldata;
  suncomplextype *Ad, *Bd;

  ldata = SM_LDATA_ZD(A);
  Ad    = SM_DATA_ZD(A);
  Bd    = SM_DATA_ZD(B);

  for (i = 0; i < ldata; i++) { Ad[i] = c * Ad[i] + Bd[i]; }

  return SUN_SUCCESS;
}

/* A = c*A + I (c is real) */
SUNErrCode SUNMatScaleAddI_ComplexDense(sunrealtype c, SUNMatrix A)
{
  sunindextype i, j, M, N;
  suncomplextype *Ad;

  M  = SM_ROWS_ZD(A);
  N  = SM_COLUMNS_ZD(A);
  Ad = SM_DATA_ZD(A);

  /* Scale all elements */
  for (i = 0; i < M * N; i++) { Ad[i] *= c; }

  /* Add identity to diagonal */
  for (j = 0; j < (M < N ? M : N); j++)
  {
    SM_ELEMENT_ZD(A, j, j) += 1.0;
  }

  return SUN_SUCCESS;
}

/* y = A*x */
SUNErrCode SUNMatMatvec_ComplexDense(SUNMatrix A, N_Vector x, N_Vector y)
{
  sunindextype i, j, M, N;
  suncomplextype *xd, *yd, *col_j;

  M  = SM_ROWS_ZD(A);
  N  = SM_COLUMNS_ZD(A);
  xd = NV_DATA_ZS(x);
  yd = NV_DATA_ZS(y);

  /* Zero y */
  memset(yd, 0, (size_t)M * sizeof(suncomplextype));

  /* y += A[:,j] * x[j] */
  for (j = 0; j < N; j++)
  {
    col_j = SM_COLUMN_ZD(A, j);
    suncomplextype xj = xd[j];
    for (i = 0; i < M; i++) { yd[i] += col_j[i] * xj; }
  }

  return SUN_SUCCESS;
}

/* y = A^H * x (conjugate transpose) */
SUNErrCode SUNMatHermitianTransposeVec_ComplexDense(SUNMatrix A, N_Vector x,
                                                    N_Vector y)
{
  sunindextype i, j, M, N;
  suncomplextype *xd, *yd, *col_j;

  M  = SM_ROWS_ZD(A);
  N  = SM_COLUMNS_ZD(A);
  xd = NV_DATA_ZS(x);
  yd = NV_DATA_ZS(y);

  /* y[j] = sum_i conj(A[i,j]) * x[i] */
  for (j = 0; j < N; j++)
  {
    col_j = SM_COLUMN_ZD(A, j);
    suncomplextype sum = 0.0;
    for (i = 0; i < M; i++) { sum += conj(col_j[i]) * xd[i]; }
    yd[j] = sum;
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNMatSpace_ComplexDense(SUNMatrix A, long int *lenrw,
                                    long int *leniw)
{
  *lenrw = 2 * (long int)SM_LDATA_ZD(A);  /* each complex = 2 reals */
  *leniw = 3;  /* M, N, ldata */
  return SUN_SUCCESS;
}

/* ================================================================
 * Complex-specific operations
 * ================================================================*/

/* B = c*A + B where c is complex */
SUNErrCode SUNMatComplexScaleAdd_ComplexDense(suncomplextype c, SUNMatrix A,
                                              SUNMatrix B)
{
  sunindextype i, ldata;
  suncomplextype *Ad, *Bd;

  ldata = SM_LDATA_ZD(A);
  Ad    = SM_DATA_ZD(A);
  Bd    = SM_DATA_ZD(B);

  for (i = 0; i < ldata; i++) { Bd[i] = c * Ad[i] + Bd[i]; }

  return SUN_SUCCESS;
}

/* A = c*A + I where c is complex */
SUNErrCode SUNMatComplexScaleAddI_ComplexDense(suncomplextype c, SUNMatrix A)
{
  sunindextype i, j, M, N;
  suncomplextype *Ad;

  M  = SM_ROWS_ZD(A);
  N  = SM_COLUMNS_ZD(A);
  Ad = SM_DATA_ZD(A);

  /* Scale all elements */
  for (i = 0; i < M * N; i++) { Ad[i] *= c; }

  /* Add identity to diagonal */
  for (j = 0; j < (M < N ? M : N); j++)
  {
    SM_ELEMENT_ZD(A, j, j) += 1.0;
  }

  return SUN_SUCCESS;
}
