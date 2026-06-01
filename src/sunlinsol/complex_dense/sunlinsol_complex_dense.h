/* -----------------------------------------------------------------
 * Complex Dense SUNLinearSolver implementation for SUNDIALS.
 *
 * Uses LAPACK zgetrf/zgetrs for complex LU factorization and solve.
 * -----------------------------------------------------------------*/

#ifndef _SUNLINSOL_COMPLEX_DENSE_H
#define _SUNLINSOL_COMPLEX_DENSE_H

#include <sundials/sundials_linearsolver.h>
#include "sunmatrix_complex_dense.h"

#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------------------------------------------
 * Content structure
 * -----------------------------------------------------------------*/

struct _SUNLinearSolverContent_ComplexDense
{
  sunindextype N;            /* matrix dimension */
  sunindextype *pivots;      /* pivot array for LU (LAPACK ipiv) */
  sunindextype last_flag;    /* last error flag */
};

typedef struct _SUNLinearSolverContent_ComplexDense *SUNLinearSolverContent_ComplexDense;

/* -----------------------------------------------------------------
 * Constructor
 * -----------------------------------------------------------------*/

SUNLinearSolver SUNLinSol_ComplexDense(N_Vector y, SUNMatrix A,
                                       SUNContext sunctx);

/* -----------------------------------------------------------------
 * SUNLinearSolver ops implementations
 * -----------------------------------------------------------------*/

SUNLinearSolver_Type SUNLinSolGetType_ComplexDense(SUNLinearSolver S);
SUNLinearSolver_ID SUNLinSolGetID_ComplexDense(SUNLinearSolver S);
SUNErrCode SUNLinSolInitialize_ComplexDense(SUNLinearSolver S);
int SUNLinSolSetup_ComplexDense(SUNLinearSolver S, SUNMatrix A);
int SUNLinSolSolve_ComplexDense(SUNLinearSolver S, SUNMatrix A,
                                N_Vector x, N_Vector b, sunrealtype tol);
SUNErrCode SUNLinSolSpace_ComplexDense(SUNLinearSolver S, long int *lenrw,
                                       long int *leniw);
SUNErrCode SUNLinSolFree_ComplexDense(SUNLinearSolver S);
sunindextype SUNLinSolLastFlag_ComplexDense(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

#endif /* _SUNLINSOL_COMPLEX_DENSE_H */
