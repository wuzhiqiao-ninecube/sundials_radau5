/* -----------------------------------------------------------------
 * Complex Band SUNLinearSolver implementation for SUNDIALS.
 *
 * Uses LAPACK zgbtrf/zgbtrs for complex band LU factorization and solve.
 * -----------------------------------------------------------------*/

#ifndef _SUNLINSOL_COMPLEX_BAND_H
#define _SUNLINSOL_COMPLEX_BAND_H

#include <sundials/sundials_linearsolver.h>
#include "sunmatrix_complex_band.h"

#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------------------------------------------
 * Content structure
 * -----------------------------------------------------------------*/

struct _SUNLinearSolverContent_ComplexBand
{
  sunindextype N;            /* matrix dimension */
  sunindextype *pivots;      /* pivot array for LU (LAPACK ipiv) */
  sunindextype last_flag;    /* last error flag */
};

typedef struct _SUNLinearSolverContent_ComplexBand *SUNLinearSolverContent_ComplexBand;

/* -----------------------------------------------------------------
 * Constructor
 * -----------------------------------------------------------------*/

SUNLinearSolver SUNLinSol_ComplexBand(N_Vector y, SUNMatrix A,
                                      SUNContext sunctx);

/* -----------------------------------------------------------------
 * SUNLinearSolver ops implementations
 * -----------------------------------------------------------------*/

SUNLinearSolver_Type SUNLinSolGetType_ComplexBand(SUNLinearSolver S);
SUNLinearSolver_ID SUNLinSolGetID_ComplexBand(SUNLinearSolver S);
SUNErrCode SUNLinSolInitialize_ComplexBand(SUNLinearSolver S);
int SUNLinSolSetup_ComplexBand(SUNLinearSolver S, SUNMatrix A);
int SUNLinSolSolve_ComplexBand(SUNLinearSolver S, SUNMatrix A,
                               N_Vector x, N_Vector b, sunrealtype tol);
SUNErrCode SUNLinSolSpace_ComplexBand(SUNLinearSolver S, long int *lenrw,
                                      long int *leniw);
SUNErrCode SUNLinSolFree_ComplexBand(SUNLinearSolver S);
sunindextype SUNLinSolLastFlag_ComplexBand(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

#endif /* _SUNLINSOL_COMPLEX_BAND_H */
