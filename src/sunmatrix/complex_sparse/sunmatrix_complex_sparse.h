/* -----------------------------------------------------------------
 * Complex Sparse (CSC) SUNMatrix implementation for SUNDIALS.
 *
 * Uses C99 double _Complex as the base element type.
 * CSC (Compressed Sparse Column) format, compatible with KLU.
 * -----------------------------------------------------------------*/

#ifndef _SUNMATRIX_COMPLEX_SPARSE_H
#define _SUNMATRIX_COMPLEX_SPARSE_H

#include <complex.h>
#include <sundials/sundials_matrix.h>
#include "nvector_complex_serial.h"

#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------------------------------------------
 * Content structure
 * -----------------------------------------------------------------*/

struct _SUNMatrixContent_ComplexSparse
{
  sunindextype M;            /* number of rows */
  sunindextype N;            /* number of columns */
  sunindextype NNZ;          /* number of nonzeros allocated */
  sunindextype NP;           /* N (number of column pointers - 1) */
  suncomplextype *data;      /* nonzero values [NNZ] */
  sunindextype *indexvals;   /* row indices [NNZ] */
  sunindextype *indexptrs;   /* column pointers [N+1] */
};

typedef struct _SUNMatrixContent_ComplexSparse *SUNMatrixContent_ComplexSparse;

/* -----------------------------------------------------------------
 * Access macros
 * -----------------------------------------------------------------*/

#define SM_CONTENT_ZS(A)      ((SUNMatrixContent_ComplexSparse)(A->content))
#define SM_ROWS_ZS(A)         (SM_CONTENT_ZS(A)->M)
#define SM_COLUMNS_ZS(A)      (SM_CONTENT_ZS(A)->N)
#define SM_NNZ_ZS(A)          (SM_CONTENT_ZS(A)->NNZ)
#define SM_NP_ZS(A)           (SM_CONTENT_ZS(A)->NP)
#define SM_DATA_ZS(A)         (SM_CONTENT_ZS(A)->data)
#define SM_INDEXVALS_ZS(A)    (SM_CONTENT_ZS(A)->indexvals)
#define SM_INDEXPTRS_ZS(A)    (SM_CONTENT_ZS(A)->indexptrs)

/* -----------------------------------------------------------------
 * Constructors / destructors
 * -----------------------------------------------------------------*/

/* Create sparse complex matrix with given dimensions and NNZ capacity */
SUNMatrix SUNMatNew_ComplexSparse(sunindextype M, sunindextype N,
                                   sunindextype NNZ, SUNContext sunctx);

/* Create from existing CSC structure (copies pattern, zeros data) */
SUNMatrix SUNMatNewFromPattern_ComplexSparse(sunindextype M, sunindextype N,
                                             sunindextype *colptrs,
                                             sunindextype *rowinds,
                                             SUNContext sunctx);

SUNMatrix SUNMatClone_ComplexSparse(SUNMatrix A);
void SUNMatDestroy_ComplexSparse(SUNMatrix A);

/* -----------------------------------------------------------------
 * SUNMatrix ops implementations
 * -----------------------------------------------------------------*/

SUNMatrix_ID SUNMatGetID_ComplexSparse(SUNMatrix A);
SUNErrCode SUNMatZero_ComplexSparse(SUNMatrix A);
SUNErrCode SUNMatCopy_ComplexSparse(SUNMatrix A, SUNMatrix B);
SUNErrCode SUNMatScaleAdd_ComplexSparse(sunrealtype c, SUNMatrix A,
                                         SUNMatrix B);
SUNErrCode SUNMatScaleAddI_ComplexSparse(sunrealtype c, SUNMatrix A);
SUNErrCode SUNMatMatvec_ComplexSparse(SUNMatrix A, N_Vector x, N_Vector y);
SUNErrCode SUNMatSpace_ComplexSparse(SUNMatrix A, long int *lenrw,
                                      long int *leniw);

/* -----------------------------------------------------------------
 * Complex-specific operations
 * -----------------------------------------------------------------*/

/* A = c*A + I where c is complex */
SUNErrCode SUNMatComplexScaleAddI_ComplexSparse(suncomplextype c, SUNMatrix A);

#ifdef __cplusplus
}
#endif

#endif /* _SUNMATRIX_COMPLEX_SPARSE_H */
