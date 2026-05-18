/* ---------------------------------------------------------------------------
 * RADAU5 — Error estimation (ESTRAD in original Fortran)
 *
 * Implements radau5_ErrorEstimate() following the Fortran ESTRAD subroutine
 * from Hairer & Wanner.  Handles both identity mass and general mass cases.
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <sundials/sundials_math.h>
#include "radau5_impl.h"

/* ---------------------------------------------------------------------------
 * radau5_ErrorEstimate
 *
 * Computes the weighted RMS error norm following ESTRAD (IJOB=1).
 *
 * Fortran mapping:
 *   CONT    -> rmem->tmp1   (work vector, holds the error correction)
 *   F2      -> rmem->tmp2   (holds sum_k(dd[k]/h)*z[k], possibly after M*)
 *   Y0      -> rmem->fn     (f(tn, ycur) saved at start of step)
 *   Y       -> rmem->ycur   (current accepted solution)
 *   SCAL    -> rmem->scal
 *   scratch -> rmem->tmp3   (linear combo intermediate / N_VDiv result)
 *
 * Returns RADAU5_SUCCESS (0) on success, RADAU5_RHSFUNC_FAIL on RHS error.
 * *err is set to MAX(SQRT(sum/n), 1e-10).
 * ---------------------------------------------------------------------------*/
int radau5_ErrorEstimate(Radau5Mem rmem, sunrealtype* err)
{
  sunindextype n = rmem->n;
  sunrealtype  h = rmem->h;
  int ns = rmem->ns;

  /* ------------------------------------------------------------------
   * First pass: form F2 = sum_k (dd[k]/h) * z[k], then CONT = F2 + fn.
   *
   * F2   -> rmem->tmp2
   * CONT -> rmem->tmp1
   * ------------------------------------------------------------------ */

  /* Build coefficient array for N_VLinearCombination */
  sunrealtype cvals[RADAU5_NS_MAX];
  N_Vector Xvecs[RADAU5_NS_MAX];
  for (int k = 0; k < ns; k++)
  {
    cvals[k] = rmem->dd[k] / h;
    Xvecs[k] = rmem->z[k];
  }

  if (rmem->M != NULL)
  {
    /* General mass: raw combo into tmp3, then M*tmp3 → tmp2 (f2) */
    N_VLinearCombination(ns, cvals, Xvecs, rmem->tmp3);
    radau5_MassMult(rmem, rmem->tmp3, rmem->tmp2);
  }
  else
  {
    /* Identity mass: combo directly into tmp2 (f2) */
    N_VLinearCombination(ns, cvals, Xvecs, rmem->tmp2);
  }

  /* CONT = F2 + fn */
  N_VLinearSum(SUN_RCONST(1.0), rmem->tmp2, SUN_RCONST(1.0), rmem->fn, rmem->tmp1);

  /* Solve E1 * cont = cont  (in-place, b and x share the same vector) */
  SUNLinSolSolve(rmem->LS_E1, rmem->E1, rmem->tmp1, rmem->tmp1, SUN_RCONST(0.0));
  rmem->nsol++;

  /* Weighted RMS norm: err = sqrt( sum((cont/scal)^2) / n ) */
  N_VDiv(rmem->tmp1, rmem->scal, rmem->tmp3);
  sunrealtype dot = N_VDotProd(rmem->tmp3, rmem->tmp3);
  *err = SUNMAX(SUNRsqrt(dot / (sunrealtype)n), SUN_RCONST(1.0e-10));

  if (*err < SUN_RCONST(1.0)) return RADAU5_SUCCESS;

  /* ------------------------------------------------------------------
   * Second pass (only when err >= 1 and this is the first step or a
   * rejected step): re-evaluate f at the corrected point and redo.
   *
   * Fortran:
   *   CONT(I) = Y(I) + CONT(I)          ! corrected y
   *   CALL FCN(N, X, CONT, F1, ...)
   *   CONT(I) = F1(I) + F2(I)
   *   CALL SOL(E1, CONT)
   * ------------------------------------------------------------------ */
  if (rmem->first || rmem->reject)
  {
    /* corrected y -> tmp1:  tmp1 = ycur + tmp1 */
    N_VLinearSum(SUN_RCONST(1.0), rmem->ycur, SUN_RCONST(1.0), rmem->tmp1, rmem->tmp1);

    /* evaluate f at corrected point; result goes into f[0] (scratch) */
    int retval = rmem->rhs(rmem->tn, rmem->tmp1, rmem->f[0], rmem->user_data);
    rmem->nfcn++;
    if (retval < 0) return RADAU5_RHSFUNC_FAIL;
    if (retval > 0) return RADAU5_RHSFUNC_RECVR;

    /* CONT = F1 + F2:  tmp1 = f[0] + tmp2 */
    N_VLinearSum(SUN_RCONST(1.0), rmem->f[0], SUN_RCONST(1.0), rmem->tmp2, rmem->tmp1);

    SUNLinSolSolve(rmem->LS_E1, rmem->E1, rmem->tmp1, rmem->tmp1,
                   SUN_RCONST(0.0));
    rmem->nsol++;

    N_VDiv(rmem->tmp1, rmem->scal, rmem->tmp3);
    dot = N_VDotProd(rmem->tmp3, rmem->tmp3);
    *err = SUNMAX(SUNRsqrt(dot / (sunrealtype)n), SUN_RCONST(1.0e-10));
  }

  return RADAU5_SUCCESS;
}
