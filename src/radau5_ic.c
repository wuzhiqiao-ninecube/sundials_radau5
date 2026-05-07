/* ---------------------------------------------------------------------------
 * RADAU5 — Consistent initial conditions for index-1 DAEs
 *
 * Implements radau5_CalcIC() / Radau5CalcIC().
 *
 * Strategy: Newton iteration on the algebraic equations only.
 * The E1 matrix (fac1*I - J) with a large fac1 makes the differential
 * rows nearly identity so their corrections stay near zero, while the
 * algebraic rows (id[i] == 0) are driven to satisfy f_a(t0, y) = 0.
 *
 * Up to MAXITER Newton steps are taken.  Convergence is declared when
 * the algebraic residual RMS norm falls below the tolerance.
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <sundials/sundials_math.h>
#include "radau5_impl.h"

#define MAXITER 10
#define IC_TOL  SUN_RCONST(1.0e-6)
/* Large fac1 makes differential rows of E1 behave like (fac1*I) so the
 * Newton correction for differential variables is ~ rhs/fac1 ≈ 0.          */
#define IC_FAC1 SUN_RCONST(1.0e6)

int radau5_CalcIC(Radau5Mem rmem, N_Vector id)
{
  sunindextype i, n = rmem->n;
  int          iter, retval;

  sunrealtype* id_data   = N_VGetArrayPointer(id);
  sunrealtype* ycur      = N_VGetArrayPointer(rmem->ycur);
  sunrealtype* fn_data   = N_VGetArrayPointer(rmem->fn);
  sunrealtype* rhs_vec   = N_VGetArrayPointer(rmem->tmp1);  /* Newton RHS */
  sunrealtype* delta     = N_VGetArrayPointer(rmem->tmp2);  /* correction  */

  /* -----------------------------------------------------------------------
   * Compute scal, evaluate f(t0, y0), and build the Jacobian so that
   * radau5_BuildE1 has a valid rmem->J to work from.
   * ----------------------------------------------------------------------- */
  radau5_ComputeScal(rmem, rmem->ycur);

  retval = rmem->rhs(rmem->tn, rmem->ycur, rmem->fn, rmem->user_data);
  rmem->nfcn++;
  if (retval != 0) return RADAU5_RHSFUNC_FAIL;

  if (rmem->jac != NULL)
  {
    retval = rmem->jac(rmem->tn, rmem->ycur, rmem->fn, rmem->J,
                       rmem->user_data, rmem->tmp1, rmem->tmp2, rmem->tmp3);
    if (retval != 0) return RADAU5_RHSFUNC_FAIL;
  }
  else
  {
    retval = radau5_DQJacDense(rmem, rmem->tn, rmem->ycur, rmem->fn);
    if (retval != RADAU5_SUCCESS) return retval;
  }
  rmem->njac++;

  /* -----------------------------------------------------------------------
   * Build and factor E1 with a large fac1 so that differential rows are
   * dominated by fac1*I and algebraic rows (where M row = 0) reduce to -J.
   * ----------------------------------------------------------------------- */
  retval = radau5_BuildE1(rmem, IC_FAC1);
  if (retval != RADAU5_SUCCESS) return RADAU5_LSETUP_FAIL;

  retval = radau5_DecompE1(rmem);
  if (retval != RADAU5_SUCCESS) return RADAU5_LSETUP_FAIL;

  /* -----------------------------------------------------------------------
   * Newton loop
   * ----------------------------------------------------------------------- */
  for (iter = 0; iter < MAXITER; iter++)
  {
    /* Evaluate f(t0, ycur) */
    retval = rmem->rhs(rmem->tn, rmem->ycur, rmem->fn, rmem->user_data);
    rmem->nfcn++;
    if (retval != 0) return RADAU5_RHSFUNC_FAIL;

    /* Form Newton RHS:
     *   for algebraic vars (id[i] == 0): rhs[i] = -fn[i]
     *   for differential vars           : rhs[i] =  0
     * The large fac1 in E1 ensures differential corrections ≈ 0.           */
    sunrealtype alg_norm = SUN_RCONST(0.0);
    sunindextype nalg    = 0;
    for (i = 0; i < n; i++)
    {
      if (id_data[i] == SUN_RCONST(0.0))
      {
        rhs_vec[i] = -fn_data[i];
        alg_norm  += fn_data[i] * fn_data[i];
        nalg++;
      }
      else
      {
        rhs_vec[i] = SUN_RCONST(0.0);
      }
    }

    /* Convergence check on algebraic residual */
    if (nalg > 0)
    {
      alg_norm = SUNRsqrt(alg_norm / (sunrealtype)nalg);
      if (alg_norm <= IC_TOL) return RADAU5_SUCCESS;
    }
    else
    {
      /* No algebraic variables — nothing to do */
      return RADAU5_SUCCESS;
    }

    /* Solve E1 * delta = rhs_vec  (tmp1 holds rhs, result overwrites it) */
    retval = SUNLinSolSolve(rmem->LS_E1, rmem->E1, rmem->tmp2, rmem->tmp1,
                            SUN_RCONST(0.0));
    rmem->nsol++;
    if (retval != 0) return RADAU5_LSOLVE_FAIL;

    /* Apply correction to algebraic variables only */
    for (i = 0; i < n; i++)
    {
      if (id_data[i] == SUN_RCONST(0.0))
        ycur[i] += delta[i];
    }
  }

  /* Final residual check after MAXITER iterations */
  retval = rmem->rhs(rmem->tn, rmem->ycur, rmem->fn, rmem->user_data);
  rmem->nfcn++;
  if (retval != 0) return RADAU5_RHSFUNC_FAIL;

  sunrealtype alg_norm = SUN_RCONST(0.0);
  sunindextype nalg    = 0;
  for (i = 0; i < n; i++)
  {
    if (id_data[i] == SUN_RCONST(0.0))
    {
      alg_norm += fn_data[i] * fn_data[i];
      nalg++;
    }
  }
  if (nalg > 0) alg_norm = SUNRsqrt(alg_norm / (sunrealtype)nalg);

  return (alg_norm <= IC_TOL) ? RADAU5_SUCCESS : RADAU5_IC_FAIL;
}

