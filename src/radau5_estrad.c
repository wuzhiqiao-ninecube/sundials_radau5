/* ---------------------------------------------------------------------------
 * RADAU5 — Error estimation (ESTRAD in original Fortran)
 *
 * Implements radau5_ErrorEstimate() following the Fortran ESTRAD subroutine
 * from Hairer & Wanner.  Only the identity-mass case (IJOB=1/2) is handled
 * here; the general mass-matrix path is left for a later phase.
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
 *   CONT  -> rmem->tmp1   (work vector, holds the error correction)
 *   F2    -> rmem->tmp2   (holds hee1*z1 + hee2*z2 + hee3*z3)
 *   Y0    -> rmem->fn     (f(tn, ycur) saved at start of step)
 *   Y     -> rmem->ycur   (current accepted solution)
 *   SCAL  -> rmem->scal
 *
 * Returns RADAU5_SUCCESS (0) on success, RADAU5_RHSFUNC_FAIL on RHS error.
 * *err is set to MAX(SQRT(sum/n), 1e-10).
 * ---------------------------------------------------------------------------*/
int radau5_ErrorEstimate(Radau5Mem rmem, sunrealtype* err)
{
  sunindextype i, n = rmem->n;
  sunrealtype  h    = rmem->h;

  /* Scaled error coefficients: hee_k = dd_k / h */
  sunrealtype hee1 = rmem->dd[0] / h;
  sunrealtype hee2 = rmem->dd[1] / h;
  sunrealtype hee3 = rmem->dd[2] / h;

  sunrealtype* z1   = N_VGetArrayPointer(rmem->z[0]);
  sunrealtype* z2   = N_VGetArrayPointer(rmem->z[1]);
  sunrealtype* z3   = N_VGetArrayPointer(rmem->z[2]);
  sunrealtype* fn   = N_VGetArrayPointer(rmem->fn);
  sunrealtype* scal = N_VGetArrayPointer(rmem->scal);
  sunrealtype* f2   = N_VGetArrayPointer(rmem->tmp2);   /* F2 in Fortran */
  sunrealtype* cont = N_VGetArrayPointer(rmem->tmp1);   /* CONT in Fortran */

  /* ------------------------------------------------------------------
   * First pass: form F2 and CONT, then solve E1 * CONT = CONT.
   *
   * Identity mass (IJOB=1):
   *   F2(I)   = HEE1*Z1(I) + HEE2*Z2(I) + HEE3*Z3(I)
   *   CONT(I) = F2(I) + Y0(I)
   *
   * General mass (IJOB=3,5):
   *   raw(I)  = HEE1*Z1(I) + HEE2*Z2(I) + HEE3*Z3(I)
   *   F2      = M * raw
   *   CONT(I) = F2(I) + Y0(I)
   * ------------------------------------------------------------------ */
  if (rmem->M != NULL)
  {
    /* General mass: compute raw combo into tmp3, then M*tmp3 → tmp2 (f2) */
    sunrealtype* raw = N_VGetArrayPointer(rmem->tmp3);
    for (i = 0; i < n; i++)
      raw[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];

    radau5_MassMult(rmem, rmem->tmp3, rmem->tmp2);
    /* f2 pointer already points to tmp2 data */

    for (i = 0; i < n; i++)
      cont[i] = f2[i] + fn[i];
  }
  else
  {
    /* Identity mass */
    for (i = 0; i < n; i++)
    {
      f2[i]   = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
      cont[i] = f2[i] + fn[i];
    }
  }

  /* Solve E1 * cont = cont  (in-place, b and x share the same vector) */
  SUNLinSolSolve(rmem->LS_E1, rmem->E1, rmem->tmp1, rmem->tmp1, SUN_RCONST(0.0));
  rmem->nsol++;

  /* Weighted RMS norm */
  sunrealtype sum = SUN_RCONST(0.0);
  for (i = 0; i < n; i++)
  {
    sunrealtype r = cont[i] / scal[i];
    sum += r * r;
  }
  *err = SUNMAX(SUNRsqrt(sum / (sunrealtype)n), SUN_RCONST(1.0e-10));

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
    sunrealtype* ycur = N_VGetArrayPointer(rmem->ycur);
    sunrealtype* f1   = N_VGetArrayPointer(rmem->f[0]);   /* reuse f1 scratch */

    /* corrected y -> tmp1 */
    for (i = 0; i < n; i++)
      cont[i] = ycur[i] + cont[i];

    /* evaluate f at corrected point; result goes into f1 (tmp scratch) */
    int retval = rmem->rhs(rmem->tn, rmem->tmp1, rmem->f[0], rmem->user_data);
    rmem->nfcn++;
    if (retval < 0) return RADAU5_RHSFUNC_FAIL;
    if (retval > 0) return RADAU5_RHSFUNC_RECVR;

    /* CONT = F1 + F2 */
    for (i = 0; i < n; i++)
      cont[i] = f1[i] + f2[i];

    SUNLinSolSolve(rmem->LS_E1, rmem->E1, rmem->tmp1, rmem->tmp1,
                   SUN_RCONST(0.0));
    rmem->nsol++;

    sum = SUN_RCONST(0.0);
    for (i = 0; i < n; i++)
    {
      sunrealtype r = cont[i] / scal[i];
      sum += r * r;
    }
    *err = SUNMAX(SUNRsqrt(sum / (sunrealtype)n), SUN_RCONST(1.0e-10));
  }

  return RADAU5_SUCCESS;
}
