/* -----------------------------------------------------------------
 * Complex Dense SUNLinearSolver implementation for SUNDIALS.
 *
 * Uses LAPACK zgetrf/zgetrs for complex LU factorization and solve.
 * -----------------------------------------------------------------*/

#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include "sunlinsol_complex_dense.h"

/* LAPACK prototypes */
extern void zgetrf_(sunindextype *m, sunindextype *n, double _Complex *a,
                    sunindextype *lda, sunindextype *ipiv, sunindextype *info);
extern void zgetrs_(const char *trans, sunindextype *n, sunindextype *nrhs,
                    double _Complex *a, sunindextype *lda, sunindextype *ipiv,
                    double _Complex *b, sunindextype *ldb, sunindextype *info);

/* ================================================================
 * Constructor
 * ================================================================*/

SUNLinearSolver SUNLinSol_ComplexDense(N_Vector y, SUNMatrix A,
                                       SUNContext sunctx)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_ComplexDense content;
  sunindextype N;

  (void)y; /* used only for type/size validation if needed */

  N = SM_ROWS_ZD(A);

  /* Create generic linear solver shell */
  S = SUNLinSolNewEmpty(sunctx);
  if (S == NULL) { return NULL; }

  /* Attach ops */
  S->ops->gettype    = SUNLinSolGetType_ComplexDense;
  S->ops->getid      = SUNLinSolGetID_ComplexDense;
  S->ops->initialize = SUNLinSolInitialize_ComplexDense;
  S->ops->setup      = SUNLinSolSetup_ComplexDense;
  S->ops->solve      = SUNLinSolSolve_ComplexDense;
  S->ops->space      = SUNLinSolSpace_ComplexDense;
  S->ops->free       = SUNLinSolFree_ComplexDense;
  S->ops->lastflag   = SUNLinSolLastFlag_ComplexDense;

  /* Allocate content */
  content = (SUNLinearSolverContent_ComplexDense)malloc(sizeof *content);
  if (content == NULL) { SUNLinSolFree(S); return NULL; }

  content->N         = N;
  content->last_flag = 0;

  /* Allocate pivot array */
  content->pivots = (sunindextype *)calloc((size_t)N, sizeof(sunindextype));
  if (content->pivots == NULL)
  {
    free(content);
    SUNLinSolFree(S);
    return NULL;
  }

  S->content = content;

  return S;
}

/* ================================================================
 * SUNLinearSolver ops implementations
 * ================================================================*/

SUNLinearSolver_Type SUNLinSolGetType_ComplexDense(SUNLinearSolver S)
{
  (void)S;
  return SUNLINEARSOLVER_DIRECT;
}

SUNLinearSolver_ID SUNLinSolGetID_ComplexDense(SUNLinearSolver S)
{
  (void)S;
  return SUNLINEARSOLVER_CUSTOM;
}

SUNErrCode SUNLinSolInitialize_ComplexDense(SUNLinearSolver S)
{
  SUNLinearSolverContent_ComplexDense content =
      (SUNLinearSolverContent_ComplexDense)(S->content);
  content->last_flag = 0;
  return SUN_SUCCESS;
}

/* LU factorization via zgetrf */
int SUNLinSolSetup_ComplexDense(SUNLinearSolver S, SUNMatrix A)
{
  SUNLinearSolverContent_ComplexDense content =
      (SUNLinearSolverContent_ComplexDense)(S->content);
  sunindextype N = content->N;
  sunindextype info;

  zgetrf_(&N, &N, SM_DATA_ZD(A), &N, content->pivots, &info);

  content->last_flag = info;

  if (info > 0) { return SUNLS_LUFACT_FAIL; }
  if (info < 0) { return SUN_ERR_EXT_FAIL; }

  return SUN_SUCCESS;
}

/* Solve A*x = b using factored A (from Setup) */
int SUNLinSolSolve_ComplexDense(SUNLinearSolver S, SUNMatrix A,
                                N_Vector x, N_Vector b, sunrealtype tol)
{
  SUNLinearSolverContent_ComplexDense content =
      (SUNLinearSolverContent_ComplexDense)(S->content);
  sunindextype N = content->N;
  sunindextype nrhs = 1;
  sunindextype info;
  suncomplextype *xd, *bd;

  (void)tol;

  /* Copy b into x (zgetrs solves in-place) */
  xd = NV_DATA_ZS(x);
  bd = NV_DATA_ZS(b);
  memcpy(xd, bd, (size_t)N * sizeof(suncomplextype));

  zgetrs_("N", &N, &nrhs, SM_DATA_ZD(A), &N, content->pivots, xd, &N, &info);

  content->last_flag = info;

  if (info != 0) { return SUN_ERR_EXT_FAIL; }

  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSpace_ComplexDense(SUNLinearSolver S, long int *lenrw,
                                       long int *leniw)
{
  SUNLinearSolverContent_ComplexDense content =
      (SUNLinearSolverContent_ComplexDense)(S->content);
  *lenrw = 0;
  *leniw = (long int)content->N;  /* pivot array */
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolFree_ComplexDense(SUNLinearSolver S)
{
  if (S == NULL) { return SUN_SUCCESS; }

  if (S->content != NULL)
  {
    SUNLinearSolverContent_ComplexDense content =
        (SUNLinearSolverContent_ComplexDense)(S->content);
    if (content->pivots != NULL) { free(content->pivots); }
    free(content);
    S->content = NULL;
  }

  if (S->ops != NULL) { free(S->ops); S->ops = NULL; }
  free(S);

  return SUN_SUCCESS;
}

sunindextype SUNLinSolLastFlag_ComplexDense(SUNLinearSolver S)
{
  SUNLinearSolverContent_ComplexDense content =
      (SUNLinearSolverContent_ComplexDense)(S->content);
  return content->last_flag;
}
