/* -----------------------------------------------------------------
 * Complex Band SUNMatrix implementation for SUNDIALS.
 *
 * Uses C99 double _Complex as the base element type.
 * LAPACK-compatible band storage with fill-in rows for LU.
 * -----------------------------------------------------------------*/

#ifndef _SUNMATRIX_COMPLEX_BAND_H
#define _SUNMATRIX_COMPLEX_BAND_H

#include <complex.h>
#include <sundials/sundials_matrix.h>
#include "nvector_complex_serial.h"

#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------------------------------------------
 * Content structure
 * -----------------------------------------------------------------*/

struct _SUNMatrixContent_ComplexBand
{
  sunindextype M;            /* number of rows */
  sunindextype N;            /* number of columns */
  sunindextype mu;           /* upper bandwidth (actual) */
  sunindextype ml;           /* lower bandwidth (actual) */
  sunindextype s_mu;         /* stored upper bandwidth (>= mu, for LU fill-in) */
  sunindextype ldim;         /* leading dimension = s_mu + ml + 1 */
  suncomplextype *data;      /* column-major data [N * ldim] */
  sunindextype ldata;        /* N * ldim */
  suncomplextype **cols;     /* column pointer array [N] */
};

typedef struct _SUNMatrixContent_ComplexBand *SUNMatrixContent_ComplexBand;

/* -----------------------------------------------------------------
 * Access macros
 * -----------------------------------------------------------------*/

#define SM_CONTENT_ZB(A)      ((SUNMatrixContent_ComplexBand)(A->content))
#define SM_ROWS_ZB(A)         (SM_CONTENT_ZB(A)->M)
#define SM_COLUMNS_ZB(A)      (SM_CONTENT_ZB(A)->N)
#define SM_UBAND_ZB(A)        (SM_CONTENT_ZB(A)->mu)
#define SM_LBAND_ZB(A)        (SM_CONTENT_ZB(A)->ml)
#define SM_SUBAND_ZB(A)       (SM_CONTENT_ZB(A)->s_mu)
#define SM_LDIM_ZB(A)         (SM_CONTENT_ZB(A)->ldim)
#define SM_DATA_ZB(A)         (SM_CONTENT_ZB(A)->data)
#define SM_LDATA_ZB(A)        (SM_CONTENT_ZB(A)->ldata)
#define SM_COLS_ZB(A)         (SM_CONTENT_ZB(A)->cols)
#define SM_COLUMN_ZB(A, j)    (SM_CONTENT_ZB(A)->cols[j] + SM_SUBAND_ZB(A))
#define SM_ELEMENT_ZB(A, i, j) (SM_CONTENT_ZB(A)->cols[j][(i)-(j)+SM_SUBAND_ZB(A)])

/* -----------------------------------------------------------------
 * Constructors / destructors
 * -----------------------------------------------------------------*/

/* Create band matrix: s_mu = mu + ml (default storage for LU fill-in) */
SUNMatrix SUNMatNew_ComplexBand(sunindextype N, sunindextype mu,
                                sunindextype ml, SUNContext sunctx);

/* Create with explicit stored upper bandwidth */
SUNMatrix SUNMatNew_ComplexBandStorage(sunindextype N, sunindextype mu,
                                       sunindextype ml, sunindextype s_mu,
                                       SUNContext sunctx);

SUNMatrix SUNMatClone_ComplexBand(SUNMatrix A);
void SUNMatDestroy_ComplexBand(SUNMatrix A);

/* -----------------------------------------------------------------
 * SUNMatrix ops implementations
 * -----------------------------------------------------------------*/

SUNMatrix_ID SUNMatGetID_ComplexBand(SUNMatrix A);
SUNErrCode SUNMatZero_ComplexBand(SUNMatrix A);
SUNErrCode SUNMatCopy_ComplexBand(SUNMatrix A, SUNMatrix B);
SUNErrCode SUNMatScaleAdd_ComplexBand(sunrealtype c, SUNMatrix A, SUNMatrix B);
SUNErrCode SUNMatScaleAddI_ComplexBand(sunrealtype c, SUNMatrix A);
SUNErrCode SUNMatMatvec_ComplexBand(SUNMatrix A, N_Vector x, N_Vector y);
SUNErrCode SUNMatSpace_ComplexBand(SUNMatrix A, long int *lenrw,
                                   long int *leniw);

/* -----------------------------------------------------------------
 * Complex-specific operations
 * -----------------------------------------------------------------*/

SUNErrCode SUNMatComplexScaleAddI_ComplexBand(suncomplextype c, SUNMatrix A);

#ifdef __cplusplus
}
#endif

#endif /* _SUNMATRIX_COMPLEX_BAND_H */
