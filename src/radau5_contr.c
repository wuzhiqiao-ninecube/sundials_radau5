/* ---------------------------------------------------------------------------
 * RADAU5 — Continuous output coefficient update (post-accepted step)
 *
 * Implements radau5_UpdateContinuousOutput() following lines 1017-1027 of
 * the original Fortran radau5.f (Hairer & Wanner).
 *
 * After an accepted step the polynomial coefficients stored in
 * cont1..cont4 allow dense output via Radau5Contr():
 *
 *   y(t) ≈ cont1[i] + s*(cont2[i] + (s-c2m1)*(cont3[i] + (s-c1m1)*cont4[i]))
 *
 * where s = (t - xsol) / hsol  (Fortran CONTR5 formula).
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <sundials/sundials_math.h>
#include "radau5_impl.h"

void radau5_UpdateContinuousOutput(Radau5Mem rmem)
{
  sunindextype i, n = rmem->n;

  sunrealtype* z1    = N_VGetArrayPointer(rmem->z1);
  sunrealtype* z2    = N_VGetArrayPointer(rmem->z2);
  sunrealtype* z3    = N_VGetArrayPointer(rmem->z3);
  sunrealtype* ycur  = N_VGetArrayPointer(rmem->ycur);
  sunrealtype* cont1 = N_VGetArrayPointer(rmem->cont1);
  sunrealtype* cont2 = N_VGetArrayPointer(rmem->cont2);
  sunrealtype* cont3 = N_VGetArrayPointer(rmem->cont3);
  sunrealtype* cont4 = N_VGetArrayPointer(rmem->cont4);

  sunrealtype c1    = rmem->c1;
  sunrealtype c2    = rmem->c2;
  sunrealtype c1m1  = rmem->c1m1;   /* c1 - 1 */
  sunrealtype c2m1  = rmem->c2m1;   /* c2 - 1 */
  sunrealtype c1mc2 = rmem->c1mc2;  /* c1 - c2 */

  /* Fortran radau5.f lines 1017-1027:
   *
   *   Z2I        = Z2(I)
   *   Z1I        = Z1(I)
   *   CONT(I+N)  = (Z2I - Z3(I)) / C2M1
   *   AK         = (Z1I - Z2I)   / C1MC2
   *   ACONT3     = Z1I / C1
   *   ACONT3     = (AK - ACONT3) / C2
   *   CONT(I+N2) = (AK - CONT(I+N)) / C1M1
   *   CONT(I+N3) = CONT(I+N2) - ACONT3
   *
   * Mapping to our arrays:
   *   CONT(I)    -> cont1[i]  (= ycur[i], the accepted solution)
   *   CONT(I+N)  -> cont2[i]
   *   CONT(I+N2) -> cont3[i]
   *   CONT(I+N3) -> cont4[i]
   *
   * Note: ycur already holds y + z3 (updated by the step driver before
   * this function is called), matching the Fortran "Y(I)=Y(I)+Z3(I)".
   */
  for (i = 0; i < n; i++)
  {
    sunrealtype z1i    = z1[i];
    sunrealtype z2i    = z2[i];
    sunrealtype z3i    = z3[i];

    sunrealtype ak     = (z1i - z2i) / c1mc2;
    sunrealtype acont3 = (ak - z1i / c1) / c2;

    cont1[i] = ycur[i];                          /* CONT(I)    = Y(I) */
    cont2[i] = (z2i - z3i) / c2m1;              /* CONT(I+N)  */
    cont3[i] = (ak  - cont2[i]) / c1m1;         /* CONT(I+N2) */
    cont4[i] = cont3[i] - acont3;               /* CONT(I+N3) */
  }

  /* Store the time-window metadata used by Radau5Contr().
   * xold is already set by the step driver (radau5_step.c) before calling
   * this function.  We only need to update xsol (= new accepted time, the
   * base for the CONTR5 interpolation) and hsol (= step size used).
   * Fortran: XSOL=X, HSOL=HOLD  (lines 1039, 1045 of radau5.f) */
  rmem->xsol = rmem->tn;    /* current accepted time (Fortran XSOL=X) */
  rmem->hsol = rmem->hold;  /* step size that produced this solution   */
}
