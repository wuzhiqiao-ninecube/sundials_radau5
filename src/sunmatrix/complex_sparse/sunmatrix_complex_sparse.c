/* -----------------------------------------------------------------
 * Complex Sparse (CSC) SUNMatrix implementation for SUNDIALS.
 *
 * Uses C99 double _Complex as the base element type.
 * CSC (Compressed Sparse Column) format, compatible with KLU.
 * -----------------------------------------------------------------*/

#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include "sunmatrix_complex_sparse.h"

/* ================================================================
 * Constructors / Destructors
 * ================================================================*/

SUNMatrix SUNMatNew_ComplexSparse(sunindextype M, sunindextype N,
                                   sunindextype NNZ, SUNContext sunctx)
{
  SUNMatrix A;
  SUNMatrixContent_ComplexSparse content;

  if (M <= 0 || N <= 0 || NNZ < 0) { return NULL; }

  /* Create generic matrix shell */
  A = SUNMatNewEmpty(sunctx);
  if (A == NULL) { return NULL; }

  /* Attach ops */
  A->ops->getid     = SUNMatGetID_ComplexSparse;
  A->ops->clone     = SUNMatClone_ComplexSparse;
  A->ops->destroy   = SUNMatDestroy_ComplexSparse;
  A->ops->zero      = SUNMatZero_ComplexSparse;
  A->ops->copy      = SUNMatCopy_ComplexSparse;
  A->ops->scaleadd  = SUNMatScaleAdd_ComplexSparse;
  A->ops->scaleaddi = SUNMatScaleAddI_ComplexSparse;
  A->ops->matvec    = SUNMatMatvec_ComplexSparse;
  A->ops->space     = SUNMatSpace_ComplexSparse;

  /* Allocate content */
  content = (SUNMatrixContent_ComplexSparse)malloc(sizeof *content);
  if (content == NULL) { SUNMatDestroy(A); return NULL; }

  content->M   = M;
  content->N   = N;
  content->NNZ = NNZ;
  content->NP  = N;

  /* Allocate data array */
  content->data = (NNZ > 0)
      ? (suncomplextype *)calloc((size_t)NNZ, sizeof(suncomplextype))
      : NULL;
  if (NNZ > 0 && content->data == NULL)
  {
    free(content); SUNMatDestroy(A); return NULL;
  }

  /* Allocate row indices */
  content->indexvals = (NNZ > 0)
      ? (sunindextype *)calloc((size_t)NNZ, sizeof(sunindextype))
      : NULL;
  if (NNZ > 0 && content->indexvals == NULL)
  {
    free(content->data); free(content); SUNMatDestroy(A); return NULL;
  }

  /* Allocate column pointers */
  content->indexptrs = (sunindextype *)calloc((size_t)(N + 1),
                                              sizeof(sunindextype));
  if (content->indexptrs == NULL)
  {
    if (content->indexvals) { free(content->indexvals); }
    if (content->data) { free(content->data); }
    free(content); SUNMatDestroy(A); return NULL;
  }

  A->content = content;

  return A;
}

SUNMatrix SUNMatNewFromPattern_ComplexSparse(sunindextype M, sunindextype N,
                                             sunindextype *colptrs,
                                             sunindextype *rowinds,
                                             SUNContext sunctx)
{
  sunindextype NNZ = colptrs[N];
  SUNMatrix A = SUNMatNew_ComplexSparse(M, N, NNZ, sunctx);
  if (A == NULL) { return NULL; }

  /* Copy sparsity pattern */
  memcpy(SM_INDEXPTRS_ZS(A), colptrs, (size_t)(N + 1) * sizeof(sunindextype));
  memcpy(SM_INDEXVALS_ZS(A), rowinds, (size_t)NNZ * sizeof(sunindextype));

  return A;
}

SUNMatrix SUNMatClone_ComplexSparse(SUNMatrix A)
{
  SUNMatrix B;
  sunindextype M, N, NNZ;

  M   = SM_ROWS_ZS(A);
  N   = SM_COLUMNS_ZS(A);
  NNZ = SM_NNZ_ZS(A);

  B = SUNMatNew_ComplexSparse(M, N, NNZ, A->sunctx);
  if (B == NULL) { return NULL; }

  /* Copy sparsity pattern */
  memcpy(SM_INDEXPTRS_ZS(B), SM_INDEXPTRS_ZS(A),
         (size_t)(N + 1) * sizeof(sunindextype));
  memcpy(SM_INDEXVALS_ZS(B), SM_INDEXVALS_ZS(A),
         (size_t)NNZ * sizeof(sunindextype));

  return B;
}

void SUNMatDestroy_ComplexSparse(SUNMatrix A)
{
  if (A == NULL) { return; }

  if (A->content != NULL)
  {
    SUNMatrixContent_ComplexSparse content = SM_CONTENT_ZS(A);
    if (content->data != NULL) { free(content->data); }
    if (content->indexvals != NULL) { free(content->indexvals); }
    if (content->indexptrs != NULL) { free(content->indexptrs); }
    free(content);
    A->content = NULL;
  }

  if (A->ops != NULL) { free(A->ops); A->ops = NULL; }
  free(A);
}

/* ================================================================
 * SUNMatrix ops implementations
 * ================================================================*/

SUNMatrix_ID SUNMatGetID_ComplexSparse(SUNMatrix A)
{
  (void)A;
  return SUNMATRIX_CUSTOM;
}

SUNErrCode SUNMatZero_ComplexSparse(SUNMatrix A)
{
  sunindextype NNZ = SM_NNZ_ZS(A);
  if (NNZ > 0)
  {
    memset(SM_DATA_ZS(A), 0, (size_t)NNZ * sizeof(suncomplextype));
  }
  return SUN_SUCCESS;
}

SUNErrCode SUNMatCopy_ComplexSparse(SUNMatrix A, SUNMatrix B)
{
  sunindextype NNZ = SM_INDEXPTRS_ZS(A)[SM_COLUMNS_ZS(A)];

  /* Assumes B has same sparsity pattern as A */
  memcpy(SM_DATA_ZS(B), SM_DATA_ZS(A), (size_t)NNZ * sizeof(suncomplextype));

  return SUN_SUCCESS;
}

/* A = c*A + B (c is real, assumes same sparsity pattern) */
SUNErrCode SUNMatScaleAdd_ComplexSparse(sunrealtype c, SUNMatrix A,
                                         SUNMatrix B)
{
  sunindextype k, NNZ;
  suncomplextype *Ad, *Bd;

  NNZ = SM_INDEXPTRS_ZS(A)[SM_COLUMNS_ZS(A)];
  Ad  = SM_DATA_ZS(A);
  Bd  = SM_DATA_ZS(B);

  for (k = 0; k < NNZ; k++) { Ad[k] = c * Ad[k] + Bd[k]; }

  return SUN_SUCCESS;
}

/* A = c*A + I (c is real) */
SUNErrCode SUNMatScaleAddI_ComplexSparse(sunrealtype c, SUNMatrix A)
{
  sunindextype j, k, N;
  suncomplextype *Ad;
  sunindextype *colptrs, *rowinds;

  N       = SM_COLUMNS_ZS(A);
  Ad      = SM_DATA_ZS(A);
  colptrs = SM_INDEXPTRS_ZS(A);
  rowinds = SM_INDEXVALS_ZS(A);

  /* Scale all elements */
  sunindextype NNZ = colptrs[N];
  for (k = 0; k < NNZ; k++) { Ad[k] *= c; }

  /* Add 1 to diagonal entries */
  for (j = 0; j < N; j++)
  {
    for (k = colptrs[j]; k < colptrs[j + 1]; k++)
    {
      if (rowinds[k] == j)
      {
        Ad[k] += 1.0;
        break;
      }
    }
  }

  return SUN_SUCCESS;
}

/* y = A*x (CSC format) */
SUNErrCode SUNMatMatvec_ComplexSparse(SUNMatrix A, N_Vector x, N_Vector y)
{
  sunindextype j, k, N, M;
  suncomplextype *xd, *yd, *Ad;
  sunindextype *colptrs, *rowinds;

  M       = SM_ROWS_ZS(A);
  N       = SM_COLUMNS_ZS(A);
  xd      = NV_DATA_ZS(x);
  yd      = NV_DATA_ZS(y);
  Ad      = SM_DATA_ZS(A);
  colptrs = SM_INDEXPTRS_ZS(A);
  rowinds = SM_INDEXVALS_ZS(A);

  /* Zero y */
  memset(yd, 0, (size_t)M * sizeof(suncomplextype));

  /* y[row] += A[row,col] * x[col] */
  for (j = 0; j < N; j++)
  {
    suncomplextype xj = xd[j];
    for (k = colptrs[j]; k < colptrs[j + 1]; k++)
    {
      yd[rowinds[k]] += Ad[k] * xj;
    }
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNMatSpace_ComplexSparse(SUNMatrix A, long int *lenrw,
                                      long int *leniw)
{
  sunindextype NNZ = SM_NNZ_ZS(A);
  sunindextype N   = SM_COLUMNS_ZS(A);
  *lenrw = 2 * (long int)NNZ;       /* complex values = 2*NNZ reals */
  *leniw = (long int)(NNZ + N + 1); /* rowinds + colptrs */
  return SUN_SUCCESS;
}

/* ================================================================
 * Complex-specific operations
 * ================================================================*/

/* A = c*A + I where c is complex */
SUNErrCode SUNMatComplexScaleAddI_ComplexSparse(suncomplextype c, SUNMatrix A)
{
  sunindextype j, k, N;
  suncomplextype *Ad;
  sunindextype *colptrs, *rowinds;

  N       = SM_COLUMNS_ZS(A);
  Ad      = SM_DATA_ZS(A);
  colptrs = SM_INDEXPTRS_ZS(A);
  rowinds = SM_INDEXVALS_ZS(A);

  /* Scale all elements by complex c */
  sunindextype NNZ = colptrs[N];
  for (k = 0; k < NNZ; k++) { Ad[k] *= c; }

  /* Add 1 to diagonal entries */
  for (j = 0; j < N; j++)
  {
    for (k = colptrs[j]; k < colptrs[j + 1]; k++)
    {
      if (rowinds[k] == j)
      {
        Ad[k] += 1.0;
        break;
      }
    }
  }

  return SUN_SUCCESS;
}
