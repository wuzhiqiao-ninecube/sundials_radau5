/* -----------------------------------------------------------------
 * Complex Dense SUNMatrix implementation for SUNDIALS.
 *
 * Uses C99 double _Complex as the base element type.
 * Column-major storage, following SUNDIALS SUNMatrix interface.
 * -----------------------------------------------------------------*/

#ifndef _SUNMATRIX_COMPLEX_DENSE_H
#define _SUNMATRIX_COMPLEX_DENSE_H

#include <complex.h>
#include <sundials/sundials_matrix.h>
#include "nvector_complex_serial.h"

#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------------------------------------------
 * Content structure
 * -----------------------------------------------------------------*/

struct _SUNMatrixContent_ComplexDense
{
  sunindextype M;            /* number of rows */
  sunindextype N;            /* number of columns */
  suncomplextype *data;      /* column-major data array [M*N] */
  sunindextype ldata;        /* M * N */
  suncomplextype **cols;     /* column pointer array [N] */
};

typedef struct _SUNMatrixContent_ComplexDense *SUNMatrixContent_ComplexDense;

/* -----------------------------------------------------------------
 * Access macros
 * -----------------------------------------------------------------*/

#define SM_CONTENT_ZD(A)      ((SUNMatrixContent_ComplexDense)(A->content))
#define SM_ROWS_ZD(A)         (SM_CONTENT_ZD(A)->M)
#define SM_COLUMNS_ZD(A)      (SM_CONTENT_ZD(A)->N)
#define SM_DATA_ZD(A)         (SM_CONTENT_ZD(A)->data)
#define SM_LDATA_ZD(A)        (SM_CONTENT_ZD(A)->ldata)
#define SM_COLS_ZD(A)         (SM_CONTENT_ZD(A)->cols)
#define SM_COLUMN_ZD(A, j)    (SM_CONTENT_ZD(A)->cols[j])
#define SM_ELEMENT_ZD(A, i, j) (SM_CONTENT_ZD(A)->cols[j][i])

/* -----------------------------------------------------------------
 * Constructors / destructors
 * -----------------------------------------------------------------*/

SUNMatrix SUNMatNew_ComplexDense(sunindextype M, sunindextype N,
                                 SUNContext sunctx);
SUNMatrix SUNMatClone_ComplexDense(SUNMatrix A);
void SUNMatDestroy_ComplexDense(SUNMatrix A);

/* -----------------------------------------------------------------
 * SUNMatrix ops implementations
 * -----------------------------------------------------------------*/

SUNMatrix_ID SUNMatGetID_ComplexDense(SUNMatrix A);
SUNErrCode SUNMatZero_ComplexDense(SUNMatrix A);
SUNErrCode SUNMatCopy_ComplexDense(SUNMatrix A, SUNMatrix B);
SUNErrCode SUNMatScaleAdd_ComplexDense(sunrealtype c, SUNMatrix A, SUNMatrix B);
SUNErrCode SUNMatScaleAddI_ComplexDense(sunrealtype c, SUNMatrix A);
SUNErrCode SUNMatMatvec_ComplexDense(SUNMatrix A, N_Vector x, N_Vector y);
SUNErrCode SUNMatHermitianTransposeVec_ComplexDense(SUNMatrix A, N_Vector x,
                                                    N_Vector y);
SUNErrCode SUNMatSpace_ComplexDense(SUNMatrix A, long int *lenrw,
                                    long int *leniw);

/* -----------------------------------------------------------------
 * Complex-specific operations
 * -----------------------------------------------------------------*/

/* B = c*A + B where c is complex */
SUNErrCode SUNMatComplexScaleAdd_ComplexDense(suncomplextype c, SUNMatrix A,
                                              SUNMatrix B);

/* A = c*A + I where c is complex */
SUNErrCode SUNMatComplexScaleAddI_ComplexDense(suncomplextype c, SUNMatrix A);

#ifdef __cplusplus
}
#endif

#endif /* _SUNMATRIX_COMPLEX_DENSE_H */
