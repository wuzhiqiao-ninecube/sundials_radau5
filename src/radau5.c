/* ---------------------------------------------------------------------------
 * radau5.c — Main driver for the RADAU5 solver
 *
 * Implements: Create/Free/Init, SetLinearSolver, optional setters,
 *             statistics getters, and the main Radau5Solve integration loop.
 *
 * Based on the Fortran code by E. Hairer and G. Wanner:
 *   "Solving Ordinary Differential Equations II"
 *   Springer Series in Computational Mathematics 14, 1996.
 * ---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_math.h>
#include <string.h>
#include "radau5_impl.h"
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_band.h>
#include <sunlinsol/sunlinsol_klu.h>
#include "radau5_colgroup.h"

/* radau5_ComputeScal is defined in radau5_linsys.c */

/* ===========================================================================
 * Radau5Create
 * ===========================================================================*/
void* Radau5Create(SUNContext sunctx)
{
  Radau5Mem rmem = (Radau5Mem)calloc(1, sizeof(struct Radau5Mem_));
  if (!rmem) return NULL;

  rmem->sunctx = sunctx;

  /* Defaults matching Fortran radau5.f */
  rmem->safe   = SUN_RCONST(0.9);
  rmem->facl   = SUN_RCONST(5.0);   /* Fortran default: upper bound on quot */
  rmem->facr   = SUN_RCONST(0.125); /* Fortran default: 1/8, lower bound on quot */
  rmem->thet   = SUN_RCONST(0.001);
  rmem->nit    = 7;
  rmem->mxstep = 100000;
  rmem->pred   = 1;
  rmem->quot1  = SUN_RCONST(1.0);
  rmem->quot2  = SUN_RCONST(1.2);  /* Fortran default: 1.2 (NOT 6) */

  /* Scalar tolerances defaults */
  rmem->rtol_s = SUN_RCONST(1.0e-6);
  rmem->atol_s = SUN_RCONST(1.0e-6);
  rmem->itol   = 0;

  /* fnewt default: max(10*uround, min(0.03, sqrt(rtol))) */
  rmem->fnewt  = SUN_RCONST(1.0e-3);

  rmem->first  = 1;
  rmem->reject = 0;
  rmem->last   = 0;
  rmem->caljac = 1;
  rmem->skipdecomp = 0;
  rmem->tol_transformed = 0;
  rmem->sparse_ls_finalized = 0;

  return (void*)rmem;
}

/* ===========================================================================
 * Radau5Free
 * ===========================================================================*/
void Radau5Free(void** radau5_mem)
{
  if (!radau5_mem || !*radau5_mem) return;
  Radau5Mem rmem = RADAU5_MEM(*radau5_mem);

  /* N_Vectors */
  if (rmem->ycur)  { N_VDestroy(rmem->ycur);  rmem->ycur  = NULL; }
  if (rmem->fn)    { N_VDestroy(rmem->fn);     rmem->fn    = NULL; }
  if (rmem->z1)    { N_VDestroy(rmem->z1);     rmem->z1    = NULL; }
  if (rmem->z2)    { N_VDestroy(rmem->z2);     rmem->z2    = NULL; }
  if (rmem->z3)    { N_VDestroy(rmem->z3);     rmem->z3    = NULL; }
  if (rmem->f1)    { N_VDestroy(rmem->f1);     rmem->f1    = NULL; }
  if (rmem->f2)    { N_VDestroy(rmem->f2);     rmem->f2    = NULL; }
  if (rmem->f3)    { N_VDestroy(rmem->f3);     rmem->f3    = NULL; }
  if (rmem->scal)  { N_VDestroy(rmem->scal);   rmem->scal  = NULL; }
  if (rmem->ewt)   { N_VDestroy(rmem->ewt);    rmem->ewt   = NULL; }
  if (rmem->tmp1)  { N_VDestroy(rmem->tmp1);   rmem->tmp1  = NULL; }
  if (rmem->tmp2)  { N_VDestroy(rmem->tmp2);   rmem->tmp2  = NULL; }
  if (rmem->tmp3)  { N_VDestroy(rmem->tmp3);   rmem->tmp3  = NULL; }
  if (rmem->cont1) { N_VDestroy(rmem->cont1);  rmem->cont1 = NULL; }
  if (rmem->cont2) { N_VDestroy(rmem->cont2);  rmem->cont2 = NULL; }
  if (rmem->cont3) { N_VDestroy(rmem->cont3);  rmem->cont3 = NULL; }
  if (rmem->cont4) { N_VDestroy(rmem->cont4);  rmem->cont4 = NULL; }
  if (rmem->rhs2)  { N_VDestroy(rmem->rhs2);   rmem->rhs2  = NULL; }
  if (rmem->sol2)  { N_VDestroy(rmem->sol2);   rmem->sol2  = NULL; }
  if (rmem->y2n)   { N_VDestroy(rmem->y2n);    rmem->y2n   = NULL; }
  if (rmem->id)    { N_VDestroy(rmem->id);      rmem->id    = NULL; }
  if (rmem->rtol_v){ N_VDestroy(rmem->rtol_v); rmem->rtol_v = NULL; }
  if (rmem->atol_v){ N_VDestroy(rmem->atol_v); rmem->atol_v = NULL; }

  /* SUNMatrices */
  if (rmem->J)      { SUNMatDestroy(rmem->J);      rmem->J      = NULL; }
  if (rmem->Jsaved) { SUNMatDestroy(rmem->Jsaved); rmem->Jsaved = NULL; }
  /* M is NOT owned by the solver — it's the user's matrix passed via SetMassFn */
  rmem->M = NULL;
  if (rmem->E1)     { SUNMatDestroy(rmem->E1);     rmem->E1     = NULL; }
  if (rmem->E2)     { SUNMatDestroy(rmem->E2);     rmem->E2     = NULL; }

  /* SUNLinearSolvers */
  if (rmem->LS_E1) { SUNLinSolFree(rmem->LS_E1); rmem->LS_E1 = NULL; }
  if (rmem->LS_E2) { SUNLinSolFree(rmem->LS_E2); rmem->LS_E2 = NULL; }

  /* Column grouping data */
  if (rmem->col_group)     { free(rmem->col_group);     rmem->col_group     = NULL; }
  if (rmem->group_offsets) { free(rmem->group_offsets);  rmem->group_offsets = NULL; }
  if (rmem->group_cols)    { free(rmem->group_cols);     rmem->group_cols    = NULL; }
  if (rmem->sp_colptrs)    { free(rmem->sp_colptrs);    rmem->sp_colptrs    = NULL; }
  if (rmem->sp_rowinds)    { free(rmem->sp_rowinds);    rmem->sp_rowinds    = NULL; }

  /* Rootfinding */
  radau5_root_Free(rmem);

  free(rmem);
  *radau5_mem = NULL;
}

/* ===========================================================================
 * Radau5Init
 * ===========================================================================*/
int Radau5Init(void* radau5_mem, Radau5RhsFn rhs, sunrealtype t0, N_Vector y0)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  Radau5Mem rmem = RADAU5_MEM(radau5_mem);

  if (!rhs) return RADAU5_ILL_INPUT;
  if (!y0)  return RADAU5_ILL_INPUT;

  rmem->rhs = rhs;
  rmem->tn  = t0;
  rmem->n   = N_VGetLength(y0);

  /* Clone y0 into ycur */
  rmem->ycur = N_VClone(y0);
  if (!rmem->ycur) return RADAU5_MEM_FAIL;
  N_VScale(SUN_RCONST(1.0), y0, rmem->ycur);

  /* Allocate all internal N_Vectors by cloning y0 */
  rmem->fn    = N_VClone(y0); if (!rmem->fn)    return RADAU5_MEM_FAIL;
  rmem->z1    = N_VClone(y0); if (!rmem->z1)    return RADAU5_MEM_FAIL;
  rmem->z2    = N_VClone(y0); if (!rmem->z2)    return RADAU5_MEM_FAIL;
  rmem->z3    = N_VClone(y0); if (!rmem->z3)    return RADAU5_MEM_FAIL;
  rmem->f1    = N_VClone(y0); if (!rmem->f1)    return RADAU5_MEM_FAIL;
  rmem->f2    = N_VClone(y0); if (!rmem->f2)    return RADAU5_MEM_FAIL;
  rmem->f3    = N_VClone(y0); if (!rmem->f3)    return RADAU5_MEM_FAIL;
  rmem->scal  = N_VClone(y0); if (!rmem->scal)  return RADAU5_MEM_FAIL;
  rmem->ewt   = N_VClone(y0); if (!rmem->ewt)   return RADAU5_MEM_FAIL;
  rmem->tmp1  = N_VClone(y0); if (!rmem->tmp1)  return RADAU5_MEM_FAIL;
  rmem->tmp2  = N_VClone(y0); if (!rmem->tmp2)  return RADAU5_MEM_FAIL;
  rmem->tmp3  = N_VClone(y0); if (!rmem->tmp3)  return RADAU5_MEM_FAIL;
  rmem->cont1 = N_VClone(y0); if (!rmem->cont1) return RADAU5_MEM_FAIL;
  rmem->cont2 = N_VClone(y0); if (!rmem->cont2) return RADAU5_MEM_FAIL;
  rmem->cont3 = N_VClone(y0); if (!rmem->cont3) return RADAU5_MEM_FAIL;
  rmem->cont4 = N_VClone(y0); if (!rmem->cont4) return RADAU5_MEM_FAIL;

  /* Initialize method constants */
  int ret = radau5_InitConstants(rmem);
  if (ret != RADAU5_SUCCESS) return ret;

  /* Reset counters and flags */
  rmem->nstep  = 0;
  rmem->naccpt = 0;
  rmem->nrejct = 0;
  rmem->nfcn   = 0;
  rmem->njac   = 0;
  rmem->ndec   = 0;
  rmem->nsol   = 0;
  rmem->nnewt  = 0;
  rmem->nsing  = 0;
  rmem->first  = 1;
  rmem->reject = 0;
  rmem->last   = 0;
  rmem->caljac = 1;
  rmem->skipdecomp = 0;
  rmem->tol_transformed = 0;
  rmem->sparse_ls_finalized = 0;
  rmem->hacc   = SUN_RCONST(0.0);
  rmem->erracc = SUN_RCONST(0.0);
  rmem->faccon = SUN_RCONST(1.0);
  rmem->theta  = SUN_RCONST(1.0);
  rmem->hhfac  = SUN_RCONST(1.0);

  rmem->setup_done = 1;
  return RADAU5_SUCCESS;
}

/* ===========================================================================
 * Radau5SetLinearSolver  (dense only for Phase 1)
 * ===========================================================================*/
int Radau5SetLinearSolver(void* radau5_mem, SUNMatrix J)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  Radau5Mem rmem = RADAU5_MEM(radau5_mem);

  if (!rmem->setup_done) return RADAU5_ILL_INPUT;
  if (!J) return RADAU5_ILL_INPUT;

  SUNMatrix_ID mid = SUNMatGetID(J);
  rmem->mat_id = mid;

  if (mid == SUNMATRIX_DENSE) {
    sunindextype n = rmem->n;

    /* Clone J for internal Jacobian storage */
    rmem->J = SUNMatClone(J);
    if (!rmem->J) return RADAU5_MEM_FAIL;

    rmem->Jsaved = SUNMatClone(J);
    if (!rmem->Jsaved) return RADAU5_MEM_FAIL;

    /* E1: n×n dense */
    rmem->E1 = SUNDenseMatrix(n, n, rmem->sunctx);
    if (!rmem->E1) return RADAU5_MEM_FAIL;

    /* E2: 2n×2n dense */
    rmem->E2 = SUNDenseMatrix(2 * n, 2 * n, rmem->sunctx);
    if (!rmem->E2) return RADAU5_MEM_FAIL;

    /* y2n: serial vector of length 2n for LS_E2 */
    rmem->y2n = N_VNew_Serial(2 * n, rmem->sunctx);
    if (!rmem->y2n) return RADAU5_MEM_FAIL;

    /* rhs2, sol2 for the realified complex solve */
    rmem->rhs2 = N_VNew_Serial(2 * n, rmem->sunctx);
    if (!rmem->rhs2) return RADAU5_MEM_FAIL;
    rmem->sol2 = N_VNew_Serial(2 * n, rmem->sunctx);
    if (!rmem->sol2) return RADAU5_MEM_FAIL;

    /* Linear solvers */
    rmem->LS_E1 = SUNLinSol_Dense(rmem->ycur, rmem->E1, rmem->sunctx);
    if (!rmem->LS_E1) return RADAU5_MEM_FAIL;

    rmem->LS_E2 = SUNLinSol_Dense(rmem->y2n, rmem->E2, rmem->sunctx);
    if (!rmem->LS_E2) return RADAU5_MEM_FAIL;

    /* Initialize linear solvers */
    int ret;
    ret = SUNLinSolInitialize(rmem->LS_E1);
    if (ret != 0) return RADAU5_LSETUP_FAIL;
    ret = SUNLinSolInitialize(rmem->LS_E2);
    if (ret != 0) return RADAU5_LSETUP_FAIL;

  } else if (mid == SUNMATRIX_BAND) {
    sunindextype n = rmem->n;

    /* Extract bandwidths */
    rmem->mu = SM_UBAND_B(J);
    rmem->ml = SM_LBAND_B(J);

    /* Clone J for internal Jacobian storage */
    rmem->J = SUNMatClone(J);
    if (!rmem->J) return RADAU5_MEM_FAIL;

    rmem->Jsaved = SUNMatClone(J);
    if (!rmem->Jsaved) return RADAU5_MEM_FAIL;

    /* E1: n×n band (SUNBandMatrix auto-sets s_mu = mu+ml for LU fill) */
    rmem->E1 = SUNBandMatrix(n, rmem->mu, rmem->ml, rmem->sunctx);
    if (!rmem->E1) return RADAU5_MEM_FAIL;

    /* E2: 2n×2n dense — the off-diagonal blocks at offset n make
       the effective bandwidth ~n, so band storage is not beneficial */
    rmem->E2 = SUNDenseMatrix(2 * n, 2 * n, rmem->sunctx);
    if (!rmem->E2) return RADAU5_MEM_FAIL;

    /* y2n, rhs2, sol2 for the realified complex solve (length 2n) */
    rmem->y2n = N_VNew_Serial(2 * n, rmem->sunctx);
    if (!rmem->y2n) return RADAU5_MEM_FAIL;
    rmem->rhs2 = N_VNew_Serial(2 * n, rmem->sunctx);
    if (!rmem->rhs2) return RADAU5_MEM_FAIL;
    rmem->sol2 = N_VNew_Serial(2 * n, rmem->sunctx);
    if (!rmem->sol2) return RADAU5_MEM_FAIL;

    /* Linear solvers: band for E1, dense for E2 */
    rmem->LS_E1 = SUNLinSol_Band(rmem->ycur, rmem->E1, rmem->sunctx);
    if (!rmem->LS_E1) return RADAU5_MEM_FAIL;

    rmem->LS_E2 = SUNLinSol_Dense(rmem->y2n, rmem->E2, rmem->sunctx);
    if (!rmem->LS_E2) return RADAU5_MEM_FAIL;

    int ret;
    ret = SUNLinSolInitialize(rmem->LS_E1);
    if (ret != 0) return RADAU5_LSETUP_FAIL;
    ret = SUNLinSolInitialize(rmem->LS_E2);
    if (ret != 0) return RADAU5_LSETUP_FAIL;

  } else if (mid == SUNMATRIX_SPARSE) {
    sunindextype n = rmem->n;
    sunindextype nnz_J = SM_NNZ_S(J);

    /* Clone J for internal Jacobian storage.
       SUNMatClone_Sparse allocates but does NOT copy the sparsity pattern,
       so we must copy the structure (indexptrs + indexvals) explicitly. */
    rmem->J = SUNMatClone(J);
    if (!rmem->J) return RADAU5_MEM_FAIL;
    SUNMatCopy(J, rmem->J);  /* copies structure + values */

    rmem->Jsaved = SUNMatClone(J);
    if (!rmem->Jsaved) return RADAU5_MEM_FAIL;
    SUNMatCopy(J, rmem->Jsaved);

    /* E1: same structure as J */
    rmem->E1 = SUNMatClone(J);
    if (!rmem->E1) return RADAU5_MEM_FAIL;
    SUNMatCopy(J, rmem->E1);  /* copy structure; values will be overwritten */

    /* E2: 2n×2n sparse CSC.
       NNZ estimate: 2 copies of J pattern + 4n diagonal entries (worst case) */
    sunindextype nnz_E2 = 2 * nnz_J + 4 * n;
    rmem->E2 = SUNSparseMatrix(2 * n, 2 * n, nnz_E2, CSC_MAT, rmem->sunctx);
    if (!rmem->E2) return RADAU5_MEM_FAIL;

    /* Vectors of length 2n */
    rmem->y2n  = N_VNew_Serial(2 * n, rmem->sunctx);
    if (!rmem->y2n) return RADAU5_MEM_FAIL;
    rmem->rhs2 = N_VNew_Serial(2 * n, rmem->sunctx);
    if (!rmem->rhs2) return RADAU5_MEM_FAIL;
    rmem->sol2 = N_VNew_Serial(2 * n, rmem->sunctx);
    if (!rmem->sol2) return RADAU5_MEM_FAIL;

    /* KLU solvers */
    rmem->LS_E1 = SUNLinSol_KLU(rmem->ycur, rmem->E1, rmem->sunctx);
    if (!rmem->LS_E1) return RADAU5_MEM_FAIL;
    rmem->LS_E2 = SUNLinSol_KLU(rmem->y2n, rmem->E2, rmem->sunctx);
    if (!rmem->LS_E2) return RADAU5_MEM_FAIL;

    int ret;
    ret = SUNLinSolInitialize(rmem->LS_E1);
    if (ret != 0) return RADAU5_LSETUP_FAIL;
    ret = SUNLinSolInitialize(rmem->LS_E2);
    if (ret != 0) return RADAU5_LSETUP_FAIL;

  } else {
    return RADAU5_ILL_INPUT;
  }

  rmem->ls_done = 1;
  return RADAU5_SUCCESS;
}

/* ===========================================================================
 * Radau5Solve — main integration loop
 * ===========================================================================*/
int Radau5Solve(void* radau5_mem, sunrealtype tout, N_Vector yout,
                sunrealtype* tret)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  Radau5Mem rmem = RADAU5_MEM(radau5_mem);

  if (!rmem->setup_done) return RADAU5_ILL_INPUT;
  if (!rmem->ls_done)    return RADAU5_ILL_INPUT;
  if (!yout || !tret)    return RADAU5_ILL_INPUT;

  /* Direction of integration */
  sunrealtype posneg = (tout >= rmem->tn) ? SUN_RCONST(1.0) : SUN_RCONST(-1.0);

  /* --- Tolerance transformation (once, on first call) ---
   * Fortran radau5.f lines 440-455:
   *   EXPM = 2/3; QUOT = ATOL/RTOL
   *   RTOL = 0.1 * RTOL^(2/3); ATOL = RTOL * QUOT
   * Also compute fnewt from the transformed rtol. */
  if (!rmem->tol_transformed)
  {
    sunrealtype expm = SUN_RCONST(2.0) / SUN_RCONST(3.0);
    if (rmem->itol == 0)
    {
      sunrealtype quot = rmem->atol_s / rmem->rtol_s;
      rmem->rtol_s = SUN_RCONST(0.1) * SUNRpowerR(rmem->rtol_s, expm);
      rmem->atol_s = rmem->rtol_s * quot;
    }
    else
    {
      sunrealtype *rv = N_VGetArrayPointer(rmem->rtol_v);
      sunrealtype *av = N_VGetArrayPointer(rmem->atol_v);
      for (sunindextype i = 0; i < rmem->n; i++)
      {
        sunrealtype quot = av[i] / rv[i];
        rv[i] = SUN_RCONST(0.1) * SUNRpowerR(rv[i], expm);
        av[i] = rv[i] * quot;
      }
    }
    /* fnewt: max(10*uround/rtol, min(0.03, sqrt(rtol))) */
    sunrealtype uround = SUN_UNIT_ROUNDOFF;
    sunrealtype rtol_val = (rmem->itol == 0) ? rmem->rtol_s
                                              : N_VMin(rmem->rtol_v);
    if (rtol_val <= SUN_RCONST(0.0)) rtol_val = SUN_RCONST(1.0e-6);
    rmem->fnewt = SUNMAX(SUN_RCONST(10.0) * uround / rtol_val,
                         SUNMIN(SUN_RCONST(0.03), sqrt(rtol_val)));
    rmem->tol_transformed = 1;
  }

  /* Compute initial scal vector */
  radau5_ComputeScal(rmem, rmem->ycur);

  /* Evaluate mass matrix if provided and not yet done */
  if (rmem->mas != NULL && !rmem->mass_evaluated)
  {
    int mret = rmem->mas(rmem->tn, rmem->M, rmem->user_data,
                         rmem->tmp1, rmem->tmp2, rmem->tmp3);
    if (mret != 0) return RADAU5_RHSFUNC_FAIL;
    rmem->mass_evaluated = 1;
  }

  /* Rebuild E1/E2 with union(J,M) pattern when M is sparse.
   * This must happen after mass matrix evaluation and before the main loop.
   * Radau5SetLinearSolver runs before Radau5SetMassFn, so M is unknown at
   * that point — we defer the pattern merge to here (once). */
  if (rmem->mat_id == SUNMATRIX_SPARSE && rmem->M != NULL
      && SUNMatGetID(rmem->M) == SUNMATRIX_SPARSE
      && !rmem->sparse_ls_finalized)
  {
    sunindextype n_loc = rmem->n;

    /* E1: union(J, M) pattern */
    SUNMatrix E1_new = radau5_SparseUnion(rmem->J, rmem->M, rmem->sunctx);
    if (!E1_new) return RADAU5_MEM_FAIL;
    SUNMatDestroy(rmem->E1);
    rmem->E1 = E1_new;

    /* E2: 2*nnz_union (diagonal blocks) + 2*nnz_M (off-diagonal blocks) */
    sunindextype nnz_union = SM_INDEXPTRS_S(E1_new)[n_loc];
    sunindextype nnz_M = SM_INDEXPTRS_S(rmem->M)[n_loc];
    sunindextype nnz_E2 = 2 * nnz_union + 2 * nnz_M;
    SUNMatDestroy(rmem->E2);
    rmem->E2 = SUNSparseMatrix(2 * n_loc, 2 * n_loc, nnz_E2, CSC_MAT,
                                rmem->sunctx);
    if (!rmem->E2) return RADAU5_MEM_FAIL;

    /* Recreate KLU solvers for the new matrices */
    SUNLinSolFree(rmem->LS_E1);
    SUNLinSolFree(rmem->LS_E2);
    rmem->LS_E1 = SUNLinSol_KLU(rmem->ycur, rmem->E1, rmem->sunctx);
    if (!rmem->LS_E1) return RADAU5_MEM_FAIL;
    rmem->LS_E2 = SUNLinSol_KLU(rmem->y2n, rmem->E2, rmem->sunctx);
    if (!rmem->LS_E2) return RADAU5_MEM_FAIL;
    SUNLinSolInitialize(rmem->LS_E1);
    SUNLinSolInitialize(rmem->LS_E2);

    rmem->sparse_ls_finalized = 1;
  }

  /* Evaluate f(t0, y0) -> fn */
  int ret = rmem->rhs(rmem->tn, rmem->ycur, rmem->fn, rmem->user_data);
  if (ret != 0) return RADAU5_RHSFUNC_FAIL;
  rmem->nfcn++;

  /* Initialize rootfinding at t0 (once) */
  if (rmem->root_active && !rmem->root_init_done) {
    ret = radau5_root_Check1(rmem);
    if (ret != RADAU5_SUCCESS) return ret;
  }

  /* Post-root re-entry handling */
  if (rmem->root_active && rmem->irfnd) {
    ret = radau5_root_Check2(rmem);
    if (ret != RADAU5_SUCCESS) return ret;
  }

  /* Set initial step size if not already set by user */
  if (rmem->h == SUN_RCONST(0.0))
    rmem->h = SUN_RCONST(1.0e-6) * posneg;

  /* Clamp to hmax if set */
  if (rmem->hmax > SUN_RCONST(0.0)) {
    sunrealtype hmx = posneg * rmem->hmax;
    if (fabs(rmem->h) > fabs(hmx)) rmem->h = hmx;
  }

  /* Store xold for continuous output */
  rmem->xold = rmem->tn;
  rmem->xsol = rmem->tn;
  rmem->hsol = rmem->h;
  rmem->hhfac = rmem->h;  /* Fortran: HHFAC=H before first step */

  /* -----------------------------------------------------------------------
   * Main integration loop
   * -----------------------------------------------------------------------*/
  while (1) {
    /* Check step count */
    if (rmem->nstep >= rmem->mxstep) return RADAU5_TOO_MANY_STEPS;

    /* Remaining distance to tout */
    sunrealtype dist = tout - rmem->tn;

    /* Check if we are already at tout */
    if (fabs(dist) <= SUN_RCONST(1.0e-14) * fabs(rmem->tn)) {
      N_VScale(SUN_RCONST(1.0), rmem->ycur, yout);
      *tret = rmem->tn;
      return RADAU5_SUCCESS;
    }

    /* Adjust h for last step to hit tout exactly */
    rmem->last = 0;
    if (posneg * (rmem->tn + rmem->h - tout) >= SUN_RCONST(0.0)) {
      rmem->h    = dist;
      rmem->last = 1;
      rmem->skipdecomp = 0;  /* h changed — must refactor */
    }

    /* Clamp to hmax */
    if (rmem->hmax > SUN_RCONST(0.0)) {
      sunrealtype hmx = posneg * rmem->hmax;
      if (fabs(rmem->h) > fabs(hmx)) {
        rmem->h = hmx;
        rmem->skipdecomp = 0;
      }
    }

    /* Take one step attempt */
    ret = radau5_Step(rmem);

    if (ret == RADAU5_SUCCESS) {
      /* Step accepted — tn and ycur already updated inside radau5_Step */

      /* Rootfinding check (before solout and tout check) */
      if (rmem->root_active) {
        int rret = radau5_root_Check3(rmem);
        if (rret == RADAU5_ROOT_RETURN) {
          /* Rewind tn and ycur to the root location */
          N_VScale(SUN_RCONST(1.0), rmem->y_root, rmem->ycur);
          rmem->tn = rmem->troot;
          /* Return to user */
          N_VScale(SUN_RCONST(1.0), rmem->ycur, yout);
          *tret = rmem->tn;
          return RADAU5_ROOT_RETURN;
        } else if (rret != RADAU5_SUCCESS) {
          return rret;
        }
      }

      /* Call solution output callback if provided */
      if (rmem->solout) {
        int sret = rmem->solout(rmem->naccpt, rmem->xold, rmem->tn,
                                rmem->ycur, radau5_mem, rmem->user_data);
        if (sret != 0) {
          N_VScale(SUN_RCONST(1.0), rmem->ycur, yout);
          *tret = rmem->tn;
          return RADAU5_SOLOUT_RETURN;
        }
      }

      /* Check if we reached tout */
      if (rmem->last || posneg * (rmem->tn - tout) >= SUN_RCONST(0.0)) {
        N_VScale(SUN_RCONST(1.0), rmem->ycur, yout);
        *tret = rmem->tn;
        return RADAU5_SUCCESS;
      }

    } else if (ret == RADAU5_SINGULAR_MATRIX) {
      return RADAU5_SINGULAR_MATRIX;
    } else if (ret == RADAU5_STEP_TOO_SMALL) {
      return RADAU5_STEP_TOO_SMALL;
    } else if (ret == RADAU5_RHSFUNC_FAIL) {
      return RADAU5_RHSFUNC_FAIL;
    } else if (ret == RADAU5_LSETUP_FAIL) {
      return RADAU5_LSETUP_FAIL;
    } else if (ret == RADAU5_LSOLVE_FAIL) {
      return RADAU5_LSOLVE_FAIL;
    } else {
      /* Step rejected — radau5_Step has already updated h for retry */
      /* Loop continues with the new (smaller) h */
    }
  }
}

/* ===========================================================================
 * Radau5CalcIC
 * ===========================================================================*/
int Radau5CalcIC(void* radau5_mem, N_Vector id)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  Radau5Mem rmem = RADAU5_MEM(radau5_mem);

  if (!rmem->setup_done) return RADAU5_ILL_INPUT;
  if (!id) return RADAU5_ILL_INPUT;

  /* Clone id into rmem->id */
  if (rmem->id) { N_VDestroy(rmem->id); rmem->id = NULL; }
  rmem->id = N_VClone(id);
  if (!rmem->id) return RADAU5_MEM_FAIL;
  N_VScale(SUN_RCONST(1.0), id, rmem->id);

  return radau5_CalcIC(rmem, id);
}

/* ===========================================================================
 * Radau5Contr — continuous output interpolation
 * ===========================================================================*/
sunrealtype Radau5Contr(void* radau5_mem, sunindextype i, sunrealtype t)
{
  if (!radau5_mem) return SUN_RCONST(0.0);
  Radau5Mem rmem = RADAU5_MEM(radau5_mem);

  sunrealtype s = (t - rmem->xsol) / rmem->hsol;

  sunrealtype* c1 = N_VGetArrayPointer(rmem->cont1);
  sunrealtype* c2 = N_VGetArrayPointer(rmem->cont2);
  sunrealtype* c3 = N_VGetArrayPointer(rmem->cont3);
  sunrealtype* c4 = N_VGetArrayPointer(rmem->cont4);

  /* Horner evaluation of degree-4 polynomial */
  return c1[i] + s * (c2[i] + (s - rmem->c2m1) * (c3[i] + (s - rmem->c1m1) * c4[i]));
}

/* ===========================================================================
 * Optional setters
 * ===========================================================================*/

int Radau5SetJacFn(void* radau5_mem, Radau5JacFn jac)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  RADAU5_MEM(radau5_mem)->jac = jac;
  return RADAU5_SUCCESS;
}

int Radau5SetMassFn(void* radau5_mem, Radau5MassFn mas, SUNMatrix M)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  Radau5Mem rmem = RADAU5_MEM(radau5_mem);
  rmem->mas = mas;
  rmem->M   = M;
  return RADAU5_SUCCESS;
}

int Radau5SetSolOutFn(void* radau5_mem, Radau5SolOutFn solout)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  RADAU5_MEM(radau5_mem)->solout = solout;
  return RADAU5_SUCCESS;
}

int Radau5SetUserData(void* radau5_mem, void* user_data)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  RADAU5_MEM(radau5_mem)->user_data = user_data;
  return RADAU5_SUCCESS;
}

int Radau5SetMaxNumSteps(void* radau5_mem, long int mxsteps)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  if (mxsteps <= 0) return RADAU5_ILL_INPUT;
  RADAU5_MEM(radau5_mem)->mxstep = mxsteps;
  return RADAU5_SUCCESS;
}

int Radau5SetMaxNewtonIter(void* radau5_mem, int maxnit)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  if (maxnit <= 0) return RADAU5_ILL_INPUT;
  RADAU5_MEM(radau5_mem)->nit = maxnit;
  return RADAU5_SUCCESS;
}

int Radau5SetInitStep(void* radau5_mem, sunrealtype h0)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  RADAU5_MEM(radau5_mem)->h = h0;
  return RADAU5_SUCCESS;
}

int Radau5SetMaxStep(void* radau5_mem, sunrealtype hmax)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  if (hmax < SUN_RCONST(0.0)) return RADAU5_ILL_INPUT;
  RADAU5_MEM(radau5_mem)->hmax = hmax;
  return RADAU5_SUCCESS;
}

int Radau5SetSafetyFactor(void* radau5_mem, sunrealtype safe)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  if (safe <= SUN_RCONST(0.0) || safe >= SUN_RCONST(1.0)) return RADAU5_ILL_INPUT;
  RADAU5_MEM(radau5_mem)->safe = safe;
  return RADAU5_SUCCESS;
}

int Radau5SetStepSizeController(void* radau5_mem, int pred)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  if (pred != 1 && pred != 2) return RADAU5_ILL_INPUT;
  RADAU5_MEM(radau5_mem)->pred = pred;
  return RADAU5_SUCCESS;
}

int Radau5SetDAEIndex(void* radau5_mem, sunindextype nind1,
                      sunindextype nind2, sunindextype nind3)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  Radau5Mem rmem = RADAU5_MEM(radau5_mem);
  rmem->nind1 = nind1;
  rmem->nind2 = nind2;
  rmem->nind3 = nind3;
  return RADAU5_SUCCESS;
}

int Radau5SetStartNewton(void* radau5_mem, int startn)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  RADAU5_MEM(radau5_mem)->startn = startn;
  return RADAU5_SUCCESS;
}

/* ===========================================================================
 * Radau5ResetForDiscontinuity — reset solver state at a discontinuity
 *
 * Matches Fortran radau5 entry behavior: sets h=h0, FIRST=.TRUE.,
 * REJECT=.FALSE., NSING=0, FACCON=1.0. The solver will use zero starting
 * values for the first Newton iteration (due to FIRST=1) and then switch
 * to extrapolation for subsequent steps (FIRST is cleared after first
 * accepted step).
 * ===========================================================================*/
int Radau5ResetForDiscontinuity(void* radau5_mem, sunrealtype h0)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  Radau5Mem rmem = RADAU5_MEM(radau5_mem);
  rmem->h      = h0;
  rmem->first  = 1;
  rmem->reject = 0;
  rmem->nsing  = 0;
  rmem->faccon = SUN_RCONST(1.0);
  return RADAU5_SUCCESS;
}

int Radau5SetSchurDecomp(void* radau5_mem, int use_schur)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  if (use_schur != 0 && use_schur != 1) return RADAU5_ILL_INPUT;
  Radau5Mem rmem = RADAU5_MEM(radau5_mem);
  rmem->use_schur = use_schur;
  /* Re-run InitConstants to populate US/TS and override T/TI/u1 */
  return radau5_InitConstants(rmem);
}

/* ===========================================================================
 * Radau5SetSparsityPattern — set sparsity pattern for sparse DQ Jacobian
 *
 * Extracts the CSC structure from S, computes column grouping immediately,
 * and stores the result in Radau5Mem for use by radau5_DQJacSparse.
 * ===========================================================================*/
int Radau5SetSparsityPattern(void* radau5_mem, SUNMatrix S)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  Radau5Mem rmem = RADAU5_MEM(radau5_mem);

  if (!rmem->setup_done) return RADAU5_ILL_INPUT;
  if (!S) return RADAU5_ILL_INPUT;

  /* Validate: must be SUNMATRIX_SPARSE in CSC format */
  if (SUNMatGetID(S) != SUNMATRIX_SPARSE) return RADAU5_ILL_INPUT;
  if (SM_SPARSETYPE_S(S) != CSC_MAT) return RADAU5_ILL_INPUT;

  /* Validate dimensions */
  sunindextype n = rmem->n;
  if (SM_COLUMNS_S(S) != n || SM_ROWS_S(S) != n) return RADAU5_ILL_INPUT;

  /* Free any previous grouping data */
  if (rmem->col_group)     { free(rmem->col_group);     rmem->col_group     = NULL; }
  if (rmem->group_offsets) { free(rmem->group_offsets);  rmem->group_offsets = NULL; }
  if (rmem->group_cols)    { free(rmem->group_cols);     rmem->group_cols    = NULL; }
  if (rmem->sp_colptrs)    { free(rmem->sp_colptrs);    rmem->sp_colptrs    = NULL; }
  if (rmem->sp_rowinds)    { free(rmem->sp_rowinds);    rmem->sp_rowinds    = NULL; }

  /* Copy sparsity pattern */
  sunindextype* src_colptrs = SM_INDEXPTRS_S(S);
  sunindextype* src_rowinds = SM_INDEXVALS_S(S);
  sunindextype nnz = src_colptrs[n];

  rmem->sp_colptrs = (sunindextype*)malloc((size_t)(n + 1) * sizeof(sunindextype));
  if (!rmem->sp_colptrs) return RADAU5_MEM_FAIL;
  memcpy(rmem->sp_colptrs, src_colptrs, (size_t)(n + 1) * sizeof(sunindextype));

  rmem->sp_rowinds = (sunindextype*)malloc((size_t)nnz * sizeof(sunindextype));
  if (!rmem->sp_rowinds) {
    free(rmem->sp_colptrs); rmem->sp_colptrs = NULL;
    return RADAU5_MEM_FAIL;
  }
  memcpy(rmem->sp_rowinds, src_rowinds, (size_t)nnz * sizeof(sunindextype));
  rmem->sp_nnz = nnz;

  /* Compute column grouping */
  int ret = radau5_ComputeColGroup(rmem->sp_colptrs, rmem->sp_rowinds,
                                   n, nnz,
                                   &rmem->col_group, &rmem->ngroups);
  if (ret != 0) {
    free(rmem->sp_colptrs); rmem->sp_colptrs = NULL;
    free(rmem->sp_rowinds); rmem->sp_rowinds = NULL;
    return RADAU5_MEM_FAIL;
  }

  /* Build group lookup */
  ret = radau5_BuildGroupLookup(rmem->col_group, n, rmem->ngroups,
                                &rmem->group_offsets, &rmem->group_cols);
  if (ret != 0) {
    free(rmem->sp_colptrs);  rmem->sp_colptrs  = NULL;
    free(rmem->sp_rowinds);  rmem->sp_rowinds  = NULL;
    free(rmem->col_group);   rmem->col_group   = NULL;
    return RADAU5_MEM_FAIL;
  }

  return RADAU5_SUCCESS;
}

int Radau5SStolerances(void* radau5_mem, sunrealtype rtol, sunrealtype atol)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  if (rtol < SUN_RCONST(0.0) || atol < SUN_RCONST(0.0)) return RADAU5_ILL_INPUT;
  Radau5Mem rmem = RADAU5_MEM(radau5_mem);
  rmem->rtol_s = rtol;
  rmem->atol_s = atol;
  rmem->itol   = 0;
  rmem->tol_transformed = 0;  /* re-apply transformation on next Radau5Solve */
  /* Free any previously set vector tolerances */
  if (rmem->rtol_v) { N_VDestroy(rmem->rtol_v); rmem->rtol_v = NULL; }
  if (rmem->atol_v) { N_VDestroy(rmem->atol_v); rmem->atol_v = NULL; }
  return RADAU5_SUCCESS;
}

int Radau5SVtolerances(void* radau5_mem, N_Vector rtol, N_Vector atol)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  if (!rtol || !atol) return RADAU5_ILL_INPUT;
  Radau5Mem rmem = RADAU5_MEM(radau5_mem);

  if (rmem->rtol_v) { N_VDestroy(rmem->rtol_v); rmem->rtol_v = NULL; }
  if (rmem->atol_v) { N_VDestroy(rmem->atol_v); rmem->atol_v = NULL; }

  rmem->rtol_v = N_VClone(rtol);
  if (!rmem->rtol_v) return RADAU5_MEM_FAIL;
  N_VScale(SUN_RCONST(1.0), rtol, rmem->rtol_v);

  rmem->atol_v = N_VClone(atol);
  if (!rmem->atol_v) return RADAU5_MEM_FAIL;
  N_VScale(SUN_RCONST(1.0), atol, rmem->atol_v);

  rmem->itol = 1;
  rmem->tol_transformed = 0;  /* re-apply transformation on next Radau5Solve */
  return RADAU5_SUCCESS;
}

/* ===========================================================================
 * Statistics getters
 * ===========================================================================*/

int Radau5GetNumSteps(void* radau5_mem, long int* nsteps)
{
  if (!radau5_mem || !nsteps) return RADAU5_MEM_NULL;
  *nsteps = RADAU5_MEM(radau5_mem)->nstep;
  return RADAU5_SUCCESS;
}

int Radau5GetNumRhsEvals(void* radau5_mem, long int* nfcn)
{
  if (!radau5_mem || !nfcn) return RADAU5_MEM_NULL;
  *nfcn = RADAU5_MEM(radau5_mem)->nfcn;
  return RADAU5_SUCCESS;
}

int Radau5GetNumJacEvals(void* radau5_mem, long int* njac)
{
  if (!radau5_mem || !njac) return RADAU5_MEM_NULL;
  *njac = RADAU5_MEM(radau5_mem)->njac;
  return RADAU5_SUCCESS;
}

int Radau5GetNumLinSolves(void* radau5_mem, long int* nsol)
{
  if (!radau5_mem || !nsol) return RADAU5_MEM_NULL;
  *nsol = RADAU5_MEM(radau5_mem)->nsol;
  return RADAU5_SUCCESS;
}

int Radau5GetNumDecomps(void* radau5_mem, long int* ndec)
{
  if (!radau5_mem || !ndec) return RADAU5_MEM_NULL;
  *ndec = RADAU5_MEM(radau5_mem)->ndec;
  return RADAU5_SUCCESS;
}

int Radau5GetNumAccSteps(void* radau5_mem, long int* naccpt)
{
  if (!radau5_mem || !naccpt) return RADAU5_MEM_NULL;
  *naccpt = RADAU5_MEM(radau5_mem)->naccpt;
  return RADAU5_SUCCESS;
}

int Radau5GetNumRejSteps(void* radau5_mem, long int* nrejct)
{
  if (!radau5_mem || !nrejct) return RADAU5_MEM_NULL;
  *nrejct = RADAU5_MEM(radau5_mem)->nrejct;
  return RADAU5_SUCCESS;
}

int Radau5GetNumNewtonIters(void* radau5_mem, long int* nnewt)
{
  if (!radau5_mem || !nnewt) return RADAU5_MEM_NULL;
  *nnewt = RADAU5_MEM(radau5_mem)->nnewt;
  return RADAU5_SUCCESS;
}

int Radau5GetCurrentStep(void* radau5_mem, sunrealtype* hcur)
{
  if (!radau5_mem || !hcur) return RADAU5_MEM_NULL;
  *hcur = RADAU5_MEM(radau5_mem)->h;
  return RADAU5_SUCCESS;
}

int Radau5GetCurrentTime(void* radau5_mem, sunrealtype* tcur)
{
  if (!radau5_mem || !tcur) return RADAU5_MEM_NULL;
  *tcur = RADAU5_MEM(radau5_mem)->tn;
  return RADAU5_SUCCESS;
}

/* ===========================================================================
 * Rootfinding API
 * ===========================================================================*/

int Radau5RootInit(void* radau5_mem, int nrtfn, Radau5RootFn g)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  Radau5Mem rmem = RADAU5_MEM(radau5_mem);

  /* If disabling rootfinding */
  if (nrtfn == 0 || g == NULL) {
    radau5_root_Free(rmem);
    return RADAU5_SUCCESS;
  }

  if (nrtfn < 0) return RADAU5_ILL_INPUT;
  if (!rmem->setup_done) return RADAU5_ILL_INPUT;

  /* Free any previous rootfinding data */
  radau5_root_Free(rmem);

  /* Allocate arrays */
  rmem->glo     = (sunrealtype*)calloc(nrtfn, sizeof(sunrealtype));
  rmem->ghi     = (sunrealtype*)calloc(nrtfn, sizeof(sunrealtype));
  rmem->grout   = (sunrealtype*)calloc(nrtfn, sizeof(sunrealtype));
  rmem->iroots  = (int*)calloc(nrtfn, sizeof(int));
  rmem->rootdir = (int*)calloc(nrtfn, sizeof(int));
  rmem->gactive = (int*)malloc(nrtfn * sizeof(int));

  if (!rmem->glo || !rmem->ghi || !rmem->grout ||
      !rmem->iroots || !rmem->rootdir || !rmem->gactive) {
    radau5_root_Free(rmem);
    return RADAU5_MEM_FAIL;
  }

  /* Initialize gactive to all active */
  for (int i = 0; i < nrtfn; i++)
    rmem->gactive[i] = 1;

  /* Clone y_root from ycur */
  rmem->y_root = N_VClone(rmem->ycur);
  if (!rmem->y_root) {
    radau5_root_Free(rmem);
    return RADAU5_MEM_FAIL;
  }

  rmem->gfun = g;
  rmem->nrtfn = nrtfn;
  rmem->nge = 0;
  rmem->troot = SUN_RCONST(0.0);
  rmem->root_active = 1;
  rmem->root_init_done = 0;
  rmem->irfnd = 0;

  return RADAU5_SUCCESS;
}

int Radau5SetRootDirection(void* radau5_mem, int* rootdir)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  Radau5Mem rmem = RADAU5_MEM(radau5_mem);
  if (!rmem->root_active) return RADAU5_ILL_INPUT;
  if (!rootdir) return RADAU5_ILL_INPUT;

  for (int i = 0; i < rmem->nrtfn; i++)
    rmem->rootdir[i] = rootdir[i];

  return RADAU5_SUCCESS;
}

int Radau5GetRootInfo(void* radau5_mem, int* rootsfound)
{
  if (!radau5_mem) return RADAU5_MEM_NULL;
  Radau5Mem rmem = RADAU5_MEM(radau5_mem);
  if (!rmem->root_active) return RADAU5_ILL_INPUT;
  if (!rootsfound) return RADAU5_ILL_INPUT;

  for (int i = 0; i < rmem->nrtfn; i++)
    rootsfound[i] = rmem->iroots[i];

  return RADAU5_SUCCESS;
}

int Radau5GetNumGEvals(void* radau5_mem, long int* ngevals)
{
  if (!radau5_mem || !ngevals) return RADAU5_MEM_NULL;
  *ngevals = RADAU5_MEM(radau5_mem)->nge;
  return RADAU5_SUCCESS;
}
