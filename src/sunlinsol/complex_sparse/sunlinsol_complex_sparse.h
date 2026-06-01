/* -----------------------------------------------------------------
 * Complex Sparse SUNLinearSolver implementation for SUNDIALS.
 *
 * Uses KLU klu_z_factor/klu_z_solve for complex sparse factorization.
 * -----------------------------------------------------------------*/

#ifndef _SUNLINSOL_COMPLEX_SPARSE_H
#define _SUNLINSOL_COMPLEX_SPARSE_H

#include <sundials/sundials_linearsolver.h>
#include <klu.h>
#include "sunmatrix_complex_sparse.h"

#ifdef __cplusplus
extern "C" {
#endif

/* -----------------------------------------------------------------
 * Content structure
 * -----------------------------------------------------------------*/

struct _SUNLinearSolverContent_ComplexSparse
{
  sunindextype N;              /* matrix dimension */
  klu_symbolic *Symbolic;      /* KLU symbolic factorization */
  klu_numeric *Numeric;        /* KLU numeric factorization */
  klu_common Common;           /* KLU common parameters */
  sunbooleantype first_factor; /* flag: first factorization needed */
  sunindextype last_flag;      /* last error flag */
};

typedef struct _SUNLinearSolverContent_ComplexSparse *SUNLinearSolverContent_ComplexSparse;

/* -----------------------------------------------------------------
 * Constructor
 * -----------------------------------------------------------------*/

SUNLinearSolver SUNLinSol_ComplexSparse(N_Vector y, SUNMatrix A,
                                        SUNContext sunctx);

/* -----------------------------------------------------------------
 * SUNLinearSolver ops implementations
 * -----------------------------------------------------------------*/

SUNLinearSolver_Type SUNLinSolGetType_ComplexSparse(SUNLinearSolver S);
SUNLinearSolver_ID SUNLinSolGetID_ComplexSparse(SUNLinearSolver S);
SUNErrCode SUNLinSolInitialize_ComplexSparse(SUNLinearSolver S);
int SUNLinSolSetup_ComplexSparse(SUNLinearSolver S, SUNMatrix A);
int SUNLinSolSolve_ComplexSparse(SUNLinearSolver S, SUNMatrix A,
                                  N_Vector x, N_Vector b, sunrealtype tol);
SUNErrCode SUNLinSolSpace_ComplexSparse(SUNLinearSolver S, long int *lenrw,
                                         long int *leniw);
SUNErrCode SUNLinSolFree_ComplexSparse(SUNLinearSolver S);
sunindextype SUNLinSolLastFlag_ComplexSparse(SUNLinearSolver S);

#ifdef __cplusplus
}
#endif

#endif /* _SUNLINSOL_COMPLEX_SPARSE_H */
