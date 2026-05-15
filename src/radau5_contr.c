/* ---------------------------------------------------------------------------
 * RADAU5 — Continuous output coefficient update (post-accepted step)
 *
 * Implements radau5_UpdateContinuousOutput() using the general Newton
 * divided-difference algorithm from Fortran radau.f (Hairer & Wanner).
 *
 * For ns=3 this reproduces the original radau5.f lines 1097-1107.
 * For ns=5,7 this implements the general algorithm from radau.f lines
 * 1252-1264 / 1445-1457.
 *
 * The polynomial is stored in cont[0..ns]:
 *   cont[0][i] = ycur[i]  (the accepted solution, node s=1)
 *   cont[1..ns][i] = divided-difference coefficients
 *
 * Evaluation via Radau5Contr uses the Horner form on Newton basis:
 *   s = (t - xsol)/hsol + 1
 *   val = cont[ns][i]
 *   for k = ns-1, ..., 0:  val = cont[k][i] + (s - node[k]) * val
 *
 * where nodes are: node[0]=1 (=c[ns-1]), node[k]=c[ns-1-k] for k=1..ns-1,
 * node[ns]=0 (implicit, the step start).
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <sundials/sundials_math.h>
#include "radau5_impl.h"

void radau5_UpdateContinuousOutput(Radau5Mem rmem)
{
  sunindextype i, n = rmem->n;
  int ns = rmem->ns;

  /* --- General divided-difference algorithm (works for all ns) ---
   *
   * Fortran radau.f (ns=5 block, lines 1252-1264; ns=7 block, lines 1445-1457):
   *
   *   CONT(I+NS*N) = ZZ(I) / C(1)
   *   DO K=1,NS-1
   *     FACT = 1/(C(NS-K) - C(NS-K+1))
   *     CONT(I+K*N) = (ZZ(I+(NS-K-1)*N) - ZZ(I+(NS-K)*N)) * FACT
   *   END DO
   *   DO J=2,NS
   *     DO K=NS,J,-1
   *       FACT = 1/(C(NS-K) - C(NS-K+J))
   *       CONT(I+K*N) = (CONT(I+K*N) - CONT(I+(K-1)*N)) * FACT
   *     END DO
   *   END DO
   *
   * Mapping: Fortran C(1..NS) → our c[0..ns-1], Fortran C(0)=0 (implicit).
   *   Fortran C(NS-K)   → c[ns-k-1]  (for k>=1, so index ns-2 down to 0)
   *   Fortran C(NS-K+1) → c[ns-k]    (for k>=1, so index ns-1 down to 1)
   *   Fortran C(NS-K+J) → c[ns-k+j-1] if ns-k+j <= ns, else 0 (=C(0))
   *
   * Our cont layout: cont[0]=ycur, cont[1..ns]=divided-diff coefficients.
   * Fortran CONT(I+K*N) for K=1..NS maps to our cont[K][i].
   * Fortran CONT(I+NS*N) maps to our cont[ns][i].
   */

  sunrealtype* cont_ns = N_VGetArrayPointer(rmem->cont[ns]);

  /* Step 1: cont[ns][i] = z[0][i] / c[0]  (Fortran: CONT(I+NS*N) = ZZ(I)/C(1)) */
  {
    sunrealtype* z0d = N_VGetArrayPointer(rmem->z[0]);
    sunrealtype inv_c0 = SUN_RCONST(1.0) / rmem->c[0];
    for (i = 0; i < n; i++)
      cont_ns[i] = z0d[i] * inv_c0;
  }

  /* Step 2: First divided differences cont[k][i] for k=1..ns-1
   * Fortran: CONT(I+K*N) = (ZZ(I+(NS-K-1)*N) - ZZ(I+(NS-K)*N)) / (C(NS-K) - C(NS-K+1))
   * Our:     cont[k][i]  = (z[ns-k-1][i]     - z[ns-k][i])     / (c[ns-k-1] - c[ns-k])  */
  for (int k = 1; k < ns; k++)
  {
    sunrealtype* cont_k = N_VGetArrayPointer(rmem->cont[k]);
    sunrealtype* z_hi   = N_VGetArrayPointer(rmem->z[ns - k - 1]);
    sunrealtype* z_lo   = N_VGetArrayPointer(rmem->z[ns - k]);
    sunrealtype fact = SUN_RCONST(1.0) / (rmem->c[ns - k - 1] - rmem->c[ns - k]);
    for (i = 0; i < n; i++)
      cont_k[i] = (z_hi[i] - z_lo[i]) * fact;
  }

  /* Step 3: Higher-order divided differences
   * Fortran: DO J=2,NS; DO K=NS,J,-1
   *   FACT = 1/(C(NS-K) - C(NS-K+J))
   *   CONT(I+K*N) = (CONT(I+K*N) - CONT(I+(K-1)*N)) * FACT
   *
   * Node mapping for the denominator:
   *   C(NS-K):   k < ns → c[ns-k-1];  k == ns → C(0) = 0
   *   C(NS-K+J): ns-k+j <= ns → c[ns-k+j-1]; ns-k+j > ns → impossible since k>=j
   *              Actually: ns-k+j, with k>=j and k<=ns, so ns-k+j <= ns. Good.
   *              But when ns-k+j == ns → c[ns-1]. When ns-k == 0 → k==ns → C(0)=0.
   */
  for (int j = 2; j <= ns; j++)
  {
    for (int k = ns; k >= j; k--)
    {
      /* Fortran denominator: C(NS-K) - C(NS-K+J)
       * C(NS-K): if k < ns → c[ns-k-1]; if k == ns → 0 (Fortran C(0)=0)
       * C(NS-K+J): ns-k+j >= 1 and <= ns, so → c[ns-k+j-1] */
      sunrealtype c_left  = (k < ns) ? rmem->c[ns - k - 1] : SUN_RCONST(0.0);
      sunrealtype c_right = rmem->c[ns - k + j - 1];
      sunrealtype fact = SUN_RCONST(1.0) / (c_left - c_right);

      sunrealtype* cont_k  = N_VGetArrayPointer(rmem->cont[k]);
      sunrealtype* cont_km = N_VGetArrayPointer(rmem->cont[k - 1]);
      for (i = 0; i < n; i++)
        cont_k[i] = (cont_k[i] - cont_km[i]) * fact;
    }
  }

  /* Step 4: cont[0][i] = ycur[i] (the accepted solution, node s=1=c[ns-1]) */
  {
    sunrealtype* cont0 = N_VGetArrayPointer(rmem->cont[0]);
    sunrealtype* ycurd = N_VGetArrayPointer(rmem->ycur);
    for (i = 0; i < n; i++)
      cont0[i] = ycurd[i];
  }

  /* Store the time-window metadata used by Radau5Contr().
   * xold is already set by the step driver (radau5_step.c) before calling
   * this function.  We only need to update xsol (= new accepted time, the
   * base for the CONTRA interpolation) and hsol (= step size used).
   * Fortran: XSOL=X, HSOL=HOLD */
  rmem->xsol = rmem->tn;    /* current accepted time (Fortran XSOL=X) */
  rmem->hsol = rmem->hold;  /* step size that produced this solution   */
}
