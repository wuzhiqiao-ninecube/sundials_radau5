/* -----------------------------------------------------------------
 * Complex Band SUNLinearSolver implementation for SUNDIALS.
 *
 * Uses LAPACK zgbtrf/zgbtrs for complex band LU factorization and solve.
 * -----------------------------------------------------------------*/

#include <complex.h>
#include <stdlib.h>
#include <string.h>

#include "sunlinsol_complex_band.h"

/* LAPACK prototypes */
extern void zgbtrf_(sunindextype *m, sunindextype *n, sunindextype *kl,
                    sunindextype *ku, double _Complex *ab, sunindextype *ldab,
                    sunindextype *ipiv, sunindextype *info);
extern void zgbtrs_(const char *trans, sunindextype *n, sunindextype *kl,
                    sunindextype *ku, sunindextype *nrhs, double _Complex *ab,
                    sunindextype *ldab, sunindextype *ipiv, double _Complex *b,
                    sunindextype *ldb, sunindextype *info);

/* ================================================================
 * Constructor
 * ================================================================*/

SUNLinearSolver SUNLinSol_ComplexBand(N_Vector y, SUNMatrix A,
                                      SUNContext sunctx)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_ComplexBand content;
  sunindextype N;

  (void)y;

  N = SM_COLUMNS_ZB(A);

  /* Create generic linear solver shell */
  S = SUNLinSolNewEmpty(sunctx);
  if (S == NULL) { return NULL; }

  /* Attach ops */
  S->ops->gettype    = SUNLinSolGetType_ComplexBand;
  S->ops->getid      = SUNLinSolGetID_ComplexBand;
  S->ops->initialize = SUNLinSolInitialize_ComplexBand;
  S->ops->setup      = SUNLinSolSetup_ComplexBand;
  S->ops->solve      = SUNLinSolSolve_ComplexBand;
  S->ops->space      = SUNLinSolSpace_ComplexBand;
  S->ops->free       = SUNLinSolFree_ComplexBand;
  S->ops->lastflag   = SUNLinSolLastFlag_ComplexBand;

  /* Allocate content */
  content = (SUNLinearSolverContent_ComplexBand)malloc(sizeof *content);
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

SUNLinearSolver_Type SUNLinSolGetType_ComplexBand(SUNLinearSolver S)
{
  (void)S;
  return SUNLINEARSOLVER_DIRECT;
}

SUNLinearSolver_ID SUNLinSolGetID_ComplexBand(SUNLinearSolver S)
{
  (void)S;
  return SUNLINEARSOLVER_CUSTOM;
}

SUNErrCode SUNLinSolInitialize_ComplexBand(SUNLinearSolver S)
{
  SUNLinearSolverContent_ComplexBand content =
      (SUNLinearSolverContent_ComplexBand)(S->content);
  content->last_flag = 0;
  return SUN_SUCCESS;
}

/* LU factorization via zgbtrf */
int SUNLinSolSetup_ComplexBand(SUNLinearSolver S, SUNMatrix A)
{
  SUNLinearSolverContent_ComplexBand content =
      (SUNLinearSolverContent_ComplexBand)(S->content);
  sunindextype N    = content->N;
  sunindextype ml   = SM_LBAND_ZB(A);
  sunindextype mu   = SM_UBAND_ZB(A);
  sunindextype ldab = SM_LDIM_ZB(A);
  sunindextype info;

  zgbtrf_(&N, &N, &ml, &mu, SM_DATA_ZB(A), &ldab, content->pivots, &info);

  content->last_flag = info;

  if (info > 0) { return SUNLS_LUFACT_FAIL; }
  if (info < 0) { return SUN_ERR_EXT_FAIL; }

  return SUN_SUCCESS;
}

/* Solve A*x = b using factored A (from Setup) */
int SUNLinSolSolve_ComplexBand(SUNLinearSolver S, SUNMatrix A,
                               N_Vector x, N_Vector b, sunrealtype tol)
{
  SUNLinearSolverContent_ComplexBand content =
      (SUNLinearSolverContent_ComplexBand)(S->content);
  sunindextype N    = content->N;
  sunindextype ml   = SM_LBAND_ZB(A);
  sunindextype mu   = SM_UBAND_ZB(A);
  sunindextype ldab = SM_LDIM_ZB(A);
  sunindextype nrhs = 1;
  sunindextype info;
  suncomplextype *xd, *bd;

  (void)tol;

  /* Copy b into x (zgbtrs solves in-place) */
  xd = NV_DATA_ZS(x);
  bd = NV_DATA_ZS(b);
  memcpy(xd, bd, (size_t)N * sizeof(suncomplextype));

  zgbtrs_("N", &N, &ml, &mu, &nrhs, SM_DATA_ZB(A), &ldab,
          content->pivots, xd, &N, &info);

  content->last_flag = info;

  if (info != 0) { return SUN_ERR_EXT_FAIL; }

  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSpace_ComplexBand(SUNLinearSolver S, long int *lenrw,
                                      long int *leniw)
{
  SUNLinearSolverContent_ComplexBand content =
      (SUNLinearSolverContent_ComplexBand)(S->content);
  *lenrw = 0;
  *leniw = (long int)content->N;  /* pivot array */
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolFree_ComplexBand(SUNLinearSolver S)
{
  if (S == NULL) { return SUN_SUCCESS; }

  if (S->content != NULL)
  {
    SUNLinearSolverContent_ComplexBand content =
        (SUNLinearSolverContent_ComplexBand)(S->content);
    if (content->pivots != NULL) { free(content->pivots); }
    free(content);
    S->content = NULL;
  }

  if (S->ops != NULL) { free(S->ops); S->ops = NULL; }
  free(S);

  return SUN_SUCCESS;
}

sunindextype SUNLinSolLastFlag_ComplexBand(SUNLinearSolver S)
{
  SUNLinearSolverContent_ComplexBand content =
      (SUNLinearSolverContent_ComplexBand)(S->content);
  return content->last_flag;
}
