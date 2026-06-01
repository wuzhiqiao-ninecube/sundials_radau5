/* -----------------------------------------------------------------
 * Complex Sparse SUNLinearSolver implementation for SUNDIALS.
 *
 * Uses KLU klu_z_factor/klu_z_solve for complex sparse factorization.
 * C99 double _Complex is layout-compatible with double[2] (C11 §6.2.5¶13),
 * so we cast directly to (double*) for KLU's complex interface.
 * -----------------------------------------------------------------*/

#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include "sunlinsol_complex_sparse.h"

/* ================================================================
 * Constructor
 * ================================================================*/

SUNLinearSolver SUNLinSol_ComplexSparse(N_Vector y, SUNMatrix A,
                                        SUNContext sunctx)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_ComplexSparse content;
  sunindextype N;

  (void)y;

  N = SM_COLUMNS_ZS(A);

  /* Create generic linear solver shell */
  S = SUNLinSolNewEmpty(sunctx);
  if (S == NULL) { return NULL; }

  /* Attach ops */
  S->ops->gettype    = SUNLinSolGetType_ComplexSparse;
  S->ops->getid      = SUNLinSolGetID_ComplexSparse;
  S->ops->initialize = SUNLinSolInitialize_ComplexSparse;
  S->ops->setup      = SUNLinSolSetup_ComplexSparse;
  S->ops->solve      = SUNLinSolSolve_ComplexSparse;
  S->ops->space      = SUNLinSolSpace_ComplexSparse;
  S->ops->free       = SUNLinSolFree_ComplexSparse;
  S->ops->lastflag   = SUNLinSolLastFlag_ComplexSparse;

  /* Allocate content */
  content = (SUNLinearSolverContent_ComplexSparse)malloc(sizeof *content);
  if (content == NULL) { SUNLinSolFree(S); return NULL; }

  content->N            = N;
  content->Symbolic     = NULL;
  content->Numeric      = NULL;
  content->first_factor = SUNTRUE;
  content->last_flag    = 0;

  /* Initialize KLU common */
  klu_defaults(&content->Common);

  /* Perform symbolic factorization (ordering) */
  content->Symbolic = klu_analyze((int)N, (int *)SM_INDEXPTRS_ZS(A),
                                  (int *)SM_INDEXVALS_ZS(A),
                                  &content->Common);
  if (content->Symbolic == NULL)
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

SUNLinearSolver_Type SUNLinSolGetType_ComplexSparse(SUNLinearSolver S)
{
  (void)S;
  return SUNLINEARSOLVER_DIRECT;
}

SUNLinearSolver_ID SUNLinSolGetID_ComplexSparse(SUNLinearSolver S)
{
  (void)S;
  return SUNLINEARSOLVER_CUSTOM;
}

SUNErrCode SUNLinSolInitialize_ComplexSparse(SUNLinearSolver S)
{
  SUNLinearSolverContent_ComplexSparse content =
      (SUNLinearSolverContent_ComplexSparse)(S->content);

  /* Free any existing numeric factorization */
  if (content->Numeric != NULL)
  {
    klu_z_free_numeric(&content->Numeric, &content->Common);
    content->Numeric = NULL;
  }

  content->first_factor = SUNTRUE;
  content->last_flag    = 0;

  return SUN_SUCCESS;
}

/* Numeric factorization via klu_z_factor / klu_z_refactor */
int SUNLinSolSetup_ComplexSparse(SUNLinearSolver S, SUNMatrix A)
{
  SUNLinearSolverContent_ComplexSparse content =
      (SUNLinearSolverContent_ComplexSparse)(S->content);
  int *colptrs = (int *)SM_INDEXPTRS_ZS(A);
  int *rowinds = (int *)SM_INDEXVALS_ZS(A);
  double *vals = (double *)SM_DATA_ZS(A);  /* interleaved re/im */

  if (content->first_factor || content->Numeric == NULL)
  {
    /* First factorization */
    if (content->Numeric != NULL)
    {
      klu_z_free_numeric(&content->Numeric, &content->Common);
    }
    content->Numeric = klu_z_factor(colptrs, rowinds, vals,
                                    content->Symbolic, &content->Common);
    content->first_factor = SUNFALSE;
  }
  else
  {
    /* Refactorization (same sparsity pattern, new values) */
    int ok = klu_z_refactor(colptrs, rowinds, vals,
                            content->Symbolic, content->Numeric,
                            &content->Common);
    if (!ok)
    {
      /* Refactor failed — try full factor */
      klu_z_free_numeric(&content->Numeric, &content->Common);
      content->Numeric = klu_z_factor(colptrs, rowinds, vals,
                                      content->Symbolic, &content->Common);
    }
  }

  if (content->Numeric == NULL)
  {
    content->last_flag = 1;
    return SUNLS_LUFACT_FAIL;
  }

  content->last_flag = 0;
  return SUN_SUCCESS;
}

/* Solve A*x = b using factored A */
int SUNLinSolSolve_ComplexSparse(SUNLinearSolver S, SUNMatrix A,
                                  N_Vector x, N_Vector b, sunrealtype tol)
{
  SUNLinearSolverContent_ComplexSparse content =
      (SUNLinearSolverContent_ComplexSparse)(S->content);
  sunindextype N = content->N;
  suncomplextype *xd, *bd;
  int ok;

  (void)A;
  (void)tol;

  /* Copy b into x (KLU solves in-place) */
  xd = NV_DATA_ZS(x);
  bd = NV_DATA_ZS(b);
  memcpy(xd, bd, (size_t)N * sizeof(suncomplextype));

  /* Solve */
  ok = klu_z_solve(content->Symbolic, content->Numeric, (int)N, 1,
                   (double *)xd, &content->Common);

  if (!ok)
  {
    content->last_flag = 1;
    return SUN_ERR_EXT_FAIL;
  }

  content->last_flag = 0;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSpace_ComplexSparse(SUNLinearSolver S, long int *lenrw,
                                         long int *leniw)
{
  (void)S;
  *lenrw = 0;
  *leniw = 0;  /* KLU manages its own memory */
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolFree_ComplexSparse(SUNLinearSolver S)
{
  if (S == NULL) { return SUN_SUCCESS; }

  if (S->content != NULL)
  {
    SUNLinearSolverContent_ComplexSparse content =
        (SUNLinearSolverContent_ComplexSparse)(S->content);

    if (content->Numeric != NULL)
    {
      klu_z_free_numeric(&content->Numeric, &content->Common);
    }
    if (content->Symbolic != NULL)
    {
      klu_free_symbolic(&content->Symbolic, &content->Common);
    }

    free(content);
    S->content = NULL;
  }

  if (S->ops != NULL) { free(S->ops); S->ops = NULL; }
  free(S);

  return SUN_SUCCESS;
}

sunindextype SUNLinSolLastFlag_ComplexSparse(SUNLinearSolver S)
{
  SUNLinearSolverContent_ComplexSparse content =
      (SUNLinearSolverContent_ComplexSparse)(S->content);
  return content->last_flag;
}
