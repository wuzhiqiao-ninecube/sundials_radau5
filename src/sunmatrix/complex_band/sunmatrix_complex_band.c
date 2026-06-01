/* -----------------------------------------------------------------
 * Complex Band SUNMatrix implementation for SUNDIALS.
 *
 * Uses C99 double _Complex as the base element type.
 * LAPACK-compatible band storage with fill-in rows for LU.
 * -----------------------------------------------------------------*/

#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include "sunmatrix_complex_band.h"

/* ================================================================
 * Constructors / Destructors
 * ================================================================*/

SUNMatrix SUNMatNew_ComplexBand(sunindextype N, sunindextype mu,
                                sunindextype ml, SUNContext sunctx)
{
  /* Default: s_mu = min(N-1, mu+ml) for LU fill-in */
  sunindextype s_mu = (mu + ml < N - 1) ? mu + ml : N - 1;
  return SUNMatNew_ComplexBandStorage(N, mu, ml, s_mu, sunctx);
}

SUNMatrix SUNMatNew_ComplexBandStorage(sunindextype N, sunindextype mu,
                                       sunindextype ml, sunindextype s_mu,
                                       SUNContext sunctx)
{
  SUNMatrix A;
  SUNMatrixContent_ComplexBand content;
  sunindextype j, ldim;

  if (N <= 0 || mu < 0 || ml < 0 || s_mu < mu) { return NULL; }

  ldim = s_mu + ml + 1;

  /* Create generic matrix shell */
  A = SUNMatNewEmpty(sunctx);
  if (A == NULL) { return NULL; }

  /* Attach ops */
  A->ops->getid     = SUNMatGetID_ComplexBand;
  A->ops->clone     = SUNMatClone_ComplexBand;
  A->ops->destroy   = SUNMatDestroy_ComplexBand;
  A->ops->zero      = SUNMatZero_ComplexBand;
  A->ops->copy      = SUNMatCopy_ComplexBand;
  A->ops->scaleadd  = SUNMatScaleAdd_ComplexBand;
  A->ops->scaleaddi = SUNMatScaleAddI_ComplexBand;
  A->ops->matvec    = SUNMatMatvec_ComplexBand;
  A->ops->space     = SUNMatSpace_ComplexBand;

  /* Allocate content */
  content = (SUNMatrixContent_ComplexBand)malloc(sizeof *content);
  if (content == NULL) { SUNMatDestroy(A); return NULL; }

  content->M     = N;
  content->N     = N;
  content->mu    = mu;
  content->ml    = ml;
  content->s_mu  = s_mu;
  content->ldim  = ldim;
  content->ldata = N * ldim;

  /* Allocate data array */
  content->data = (suncomplextype *)calloc((size_t)(N * ldim),
                                           sizeof(suncomplextype));
  if (content->data == NULL)
  {
    free(content);
    SUNMatDestroy(A);
    return NULL;
  }

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
  for (j = 0; j < N; j++) { content->cols[j] = content->data + j * ldim; }

  A->content = content;

  return A;
}

SUNMatrix SUNMatClone_ComplexBand(SUNMatrix A)
{
  return SUNMatNew_ComplexBandStorage(SM_COLUMNS_ZB(A), SM_UBAND_ZB(A),
                                      SM_LBAND_ZB(A), SM_SUBAND_ZB(A),
                                      A->sunctx);
}

void SUNMatDestroy_ComplexBand(SUNMatrix A)
{
  if (A == NULL) { return; }

  if (A->content != NULL)
  {
    SUNMatrixContent_ComplexBand content = SM_CONTENT_ZB(A);
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

SUNMatrix_ID SUNMatGetID_ComplexBand(SUNMatrix A)
{
  (void)A;
  return SUNMATRIX_CUSTOM;
}

SUNErrCode SUNMatZero_ComplexBand(SUNMatrix A)
{
  memset(SM_DATA_ZB(A), 0,
         (size_t)SM_LDATA_ZB(A) * sizeof(suncomplextype));
  return SUN_SUCCESS;
}

SUNErrCode SUNMatCopy_ComplexBand(SUNMatrix A, SUNMatrix B)
{
  sunindextype j, N, s_mu_A, ml_A, s_mu_B;
  suncomplextype **A_cols, **B_cols;

  N      = SM_COLUMNS_ZB(A);
  s_mu_A = SM_SUBAND_ZB(A);
  ml_A   = SM_LBAND_ZB(A);
  s_mu_B = SM_SUBAND_ZB(B);
  A_cols = SM_COLS_ZB(A);
  B_cols = SM_COLS_ZB(B);

  /* Zero B first */
  SUNMatZero_ComplexBand(B);

  /* Copy band elements */
  for (j = 0; j < N; j++)
  {
    sunindextype istart = (-s_mu_A > -s_mu_B) ? -s_mu_A : -s_mu_B;
    sunindextype iend   = ml_A;
    sunindextype i;
    suncomplextype *colA = A_cols[j] + s_mu_A;
    suncomplextype *colB = B_cols[j] + s_mu_B;
    for (i = istart; i <= iend; i++) { colB[i] = colA[i]; }
  }

  return SUN_SUCCESS;
}

/* A = c*A + B (c is real) */
SUNErrCode SUNMatScaleAdd_ComplexBand(sunrealtype c, SUNMatrix A, SUNMatrix B)
{
  sunindextype j, i, N, s_mu, ml;
  suncomplextype **A_cols, **B_cols;

  N      = SM_COLUMNS_ZB(A);
  s_mu   = SM_SUBAND_ZB(A);
  ml     = SM_LBAND_ZB(A);
  A_cols = SM_COLS_ZB(A);
  B_cols = SM_COLS_ZB(B);

  for (j = 0; j < N; j++)
  {
    suncomplextype *colA = A_cols[j] + s_mu;
    suncomplextype *colB = B_cols[j] + SM_SUBAND_ZB(B);
    sunindextype lo = (j - s_mu > 0) ? j - s_mu : 0;
    sunindextype hi = (j + ml < N - 1) ? j + ml : N - 1;
    for (i = lo; i <= hi; i++)
    {
      colA[i - j] = c * colA[i - j] + colB[i - j];
    }
  }

  return SUN_SUCCESS;
}

/* A = c*A + I (c is real) */
SUNErrCode SUNMatScaleAddI_ComplexBand(sunrealtype c, SUNMatrix A)
{
  sunindextype j, i, N;
  suncomplextype *Ad;

  N  = SM_COLUMNS_ZB(A);
  Ad = SM_DATA_ZB(A);

  /* Scale all stored elements */
  for (i = 0; i < SM_LDATA_ZB(A); i++) { Ad[i] *= c; }

  /* Add identity to diagonal */
  for (j = 0; j < N; j++)
  {
    SM_ELEMENT_ZB(A, j, j) += 1.0;
  }

  return SUN_SUCCESS;
}

/* y = A*x */
SUNErrCode SUNMatMatvec_ComplexBand(SUNMatrix A, N_Vector x, N_Vector y)
{
  sunindextype i, j, N, mu, ml;
  suncomplextype *xd, *yd, *col_j;

  N  = SM_COLUMNS_ZB(A);
  mu = SM_UBAND_ZB(A);
  ml = SM_LBAND_ZB(A);
  xd = NV_DATA_ZS(x);
  yd = NV_DATA_ZS(y);

  /* Zero y */
  memset(yd, 0, (size_t)N * sizeof(suncomplextype));

  /* Accumulate: y[i] += A[i,j] * x[j] for i in band of column j */
  for (j = 0; j < N; j++)
  {
    col_j = SM_COLUMN_ZB(A, j);  /* points to diagonal of column j */
    sunindextype lo = (j - (sunindextype)mu > 0) ? j - mu : 0;
    sunindextype hi = (j + ml < N - 1) ? j + ml : N - 1;
    suncomplextype xj = xd[j];
    for (i = lo; i <= hi; i++)
    {
      yd[i] += col_j[i - j] * xj;
    }
  }

  return SUN_SUCCESS;
}

SUNErrCode SUNMatSpace_ComplexBand(SUNMatrix A, long int *lenrw,
                                   long int *leniw)
{
  *lenrw = 2 * (long int)SM_LDATA_ZB(A);
  *leniw = 7;  /* M, N, mu, ml, s_mu, ldim, ldata */
  return SUN_SUCCESS;
}

/* ================================================================
 * Complex-specific operations
 * ================================================================*/

/* A = c*A + I where c is complex */
SUNErrCode SUNMatComplexScaleAddI_ComplexBand(suncomplextype c, SUNMatrix A)
{
  sunindextype j, i, N;
  suncomplextype *Ad;

  N  = SM_COLUMNS_ZB(A);
  Ad = SM_DATA_ZB(A);

  /* Scale all stored elements by complex c */
  for (i = 0; i < SM_LDATA_ZB(A); i++) { Ad[i] *= c; }

  /* Add identity to diagonal */
  for (j = 0; j < N; j++)
  {
    SM_ELEMENT_ZB(A, j, j) += 1.0;
  }

  return SUN_SUCCESS;
}
