/* ---------------------------------------------------------------------------
 * RADAU5 — Rootfinding (event detection) implementation
 *
 * Implements sign-change detection and Illinois method root refinement
 * using the RADAU5 continuous output polynomial (Radau5Contr).
 *
 * Algorithm follows SUNDIALS CVODE/ARKODE rootfinding:
 *   - After each accepted step, evaluate user's g(t,y) at step end
 *   - If sign change detected vs step start, refine with Illinois method
 *   - Illinois uses dense output interpolation (no extra RHS evals)
 *   - Solver stops at root, user queries which functions fired
 *
 * Reference: Hiebert & Shampine, SAND80-0180, 1980 (Illinois algorithm)
 * ---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sundials/sundials_math.h>
#include "radau5_impl.h"

/* Maximum iterations for Illinois refinement */
#define RADAU5_ROOT_MAXITER 100

/* ---------------------------------------------------------------------------
 * Helper: evaluate user's root function g(t, y) -> gout
 * ---------------------------------------------------------------------------*/
static int radau5_root_EvalG(Radau5Mem rmem, sunrealtype t, N_Vector y,
                             sunrealtype* gout)
{
  int ret = rmem->gfun(t, y, gout, rmem->user_data);
  rmem->nge++;
  return ret;
}

/* ---------------------------------------------------------------------------
 * Helper: fill y_root by interpolating all components at time t
 * t must be in [xold, xsol] (the last accepted step interval)
 * ---------------------------------------------------------------------------*/
static void radau5_root_InterpY(Radau5Mem rmem, sunrealtype t)
{
  sunrealtype* yd = N_VGetArrayPointer(rmem->y_root);
  for (sunindextype i = 0; i < rmem->n; i++)
    yd[i] = Radau5Contr((void*)rmem, i, t);
}

/* ---------------------------------------------------------------------------
 * Helper: sign function (returns -1, 0, or +1)
 * ---------------------------------------------------------------------------*/
static int radau5_sign(sunrealtype x)
{
  if (x > SUN_RCONST(0.0)) return 1;
  if (x < SUN_RCONST(0.0)) return -1;
  return 0;
}

/* ---------------------------------------------------------------------------
 * radau5_RootFind — Illinois method refinement
 *
 * On entry: glo[] has values at xold, ghi[] has values at tn.
 * A sign change has been detected. Refines the bracket [xold, tn] to
 * locate the earliest root to within ttol = 100*uround*max(|tlo|,|thi|).
 *
 * Algorithm: Single-bracket Illinois method driven by the function whose
 * root estimate is closest to thi (largest gfrac). This matches the
 * SUNDIALS CVODE/ARKODE approach (Hiebert & Shampine, SAND80-0180).
 *
 * On exit: rmem->troot and rmem->y_root hold the root location and state.
 *          rmem->grout[] holds g values at the root.
 *          rmem->iroots[] indicates which functions fired.
 * ---------------------------------------------------------------------------*/
static int radau5_RootFind(Radau5Mem rmem)
{
  int i, nrtfn = rmem->nrtfn;
  sunrealtype tlo = rmem->xold;
  sunrealtype thi = rmem->tn;

  /* Local copies of g values at bracket endpoints (Illinois modifies these) */
  sunrealtype* glo_loc = (sunrealtype*)malloc(nrtfn * sizeof(sunrealtype));
  sunrealtype* ghi_loc = (sunrealtype*)malloc(nrtfn * sizeof(sunrealtype));
  if (!glo_loc || !ghi_loc) {
    free(glo_loc); free(ghi_loc);
    return RADAU5_MEM_FAIL;
  }
  memcpy(glo_loc, rmem->glo, nrtfn * sizeof(sunrealtype));
  memcpy(ghi_loc, rmem->ghi, nrtfn * sizeof(sunrealtype));

  /* Identify which functions have valid sign changes */
  int* sgnchg = (int*)calloc(nrtfn, sizeof(int));
  if (!sgnchg) {
    free(glo_loc); free(ghi_loc);
    return RADAU5_MEM_FAIL;
  }
  for (i = 0; i < nrtfn; i++) {
    if (!rmem->gactive[i]) continue;
    if (radau5_sign(glo_loc[i]) != radau5_sign(ghi_loc[i])) {
      int dir = rmem->rootdir[i];
      if (dir == 0 || (dir > 0 && ghi_loc[i] > glo_loc[i]) ||
          (dir < 0 && ghi_loc[i] < glo_loc[i])) {
        sgnchg[i] = 1;
      }
    }
  }

  /* Find the driving function (imax): the one whose root is estimated
     closest to thi, i.e., largest |ghi/(ghi-glo)| fraction */
  int imax = -1;
  sunrealtype maxfrac = SUN_RCONST(0.0);
  for (i = 0; i < nrtfn; i++) {
    if (!sgnchg[i]) continue;
    sunrealtype dg = fabs(ghi_loc[i] - glo_loc[i]);
    if (dg < SUN_RCONST(1.0e-300)) continue;
    sunrealtype gfrac = fabs(ghi_loc[i]) / dg;
    if (imax < 0 || gfrac > maxfrac) {
      maxfrac = gfrac;
      imax = i;
    }
  }
  if (imax < 0) {
    /* No valid sign change (shouldn't happen) — return midpoint */
    rmem->troot = SUN_RCONST(0.5) * (tlo + thi);
    radau5_root_InterpY(rmem, rmem->troot);
    radau5_root_EvalG(rmem, rmem->troot, rmem->y_root, rmem->grout);
    free(glo_loc); free(ghi_loc); free(sgnchg);
    return RADAU5_ROOT_RETURN;
  }

  /* Illinois iteration — single bracket [tlo, thi] driven by imax */
  int side = 0;  /* 0=initial, 1=thi moved last, 2=tlo moved last */
  sunrealtype alpha = SUN_RCONST(1.0);

  for (int iter = 0; iter < RADAU5_ROOT_MAXITER; iter++) {
    sunrealtype ttol = SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF *
                       SUNMAX(fabs(tlo), fabs(thi));
    if (fabs(thi - tlo) <= ttol) break;

    /* Secant step using the driving function imax, with Illinois weight */
    sunrealtype tmid;
    sunrealtype dg = ghi_loc[imax] - alpha * glo_loc[imax];
    if (fabs(dg) > SUN_RCONST(1.0e-300)) {
      tmid = thi - ghi_loc[imax] * (thi - tlo) / dg;
    } else {
      tmid = SUN_RCONST(0.5) * (tlo + thi);
    }

    /* Clamp to interior */
    sunrealtype ttol2 = SUN_RCONST(0.5) * ttol;
    if (tmid <= tlo + ttol2) tmid = tlo + ttol2;
    if (tmid >= thi - ttol2) tmid = thi - ttol2;

    /* Evaluate g at tmid */
    radau5_root_InterpY(rmem, tmid);
    int gret = radau5_root_EvalG(rmem, tmid, rmem->y_root, rmem->grout);
    if (gret < 0) {
      free(glo_loc); free(ghi_loc); free(sgnchg);
      return RADAU5_ROOTFN_FAIL;
    }
    sunrealtype* gmid = rmem->grout;

    /* Update bracket based on sign of gmid[imax] */
    if (radau5_sign(gmid[imax]) == radau5_sign(glo_loc[imax])) {
      /* Root is in [tmid, thi] — move tlo to tmid */
      tlo = tmid;
      for (i = 0; i < nrtfn; i++)
        glo_loc[i] = gmid[i];
      if (side == 2) { alpha *= SUN_RCONST(2.0); } else { alpha = SUN_RCONST(1.0); }
      side = 2;
    } else {
      /* Root is in [tlo, tmid] — move thi to tmid */
      thi = tmid;
      for (i = 0; i < nrtfn; i++)
        ghi_loc[i] = gmid[i];
      if (side == 1) { alpha *= SUN_RCONST(2.0); } else { alpha = SUN_RCONST(1.0); }
      side = 1;
    }

    /* Re-evaluate imax: which sign-changing function has root closest to thi */
    imax = -1;
    maxfrac = SUN_RCONST(0.0);
    for (i = 0; i < nrtfn; i++) {
      if (!sgnchg[i]) continue;
      /* Check if this function still has a sign change in [tlo, thi] */
      if (radau5_sign(glo_loc[i]) == radau5_sign(ghi_loc[i])) {
        sgnchg[i] = 0;
        continue;
      }
      sunrealtype dg2 = fabs(ghi_loc[i] - glo_loc[i]);
      if (dg2 < SUN_RCONST(1.0e-300)) continue;
      sunrealtype gfrac = fabs(ghi_loc[i]) / dg2;
      if (imax < 0 || gfrac > maxfrac) {
        maxfrac = gfrac;
        imax = i;
      }
    }
    if (imax < 0) break; /* All roots resolved */
  }

  /* Root located: pick the bracket endpoint where |g[imax]| is smallest.
   * The standard midpoint approach fails when the interpolant doesn't cross
   * zero exactly (due to interpolation error) — in that case thi converges
   * to the near-zero point while tlo stays stuck, making the midpoint wrong.
   * Using the endpoint with smallest |g| gives the best root estimate. */
  if (imax >= 0 && fabs(glo_loc[imax]) < fabs(ghi_loc[imax])) {
    rmem->troot = tlo;
  } else {
    rmem->troot = thi;
  }
  radau5_root_InterpY(rmem, rmem->troot);
  int gret = radau5_root_EvalG(rmem, rmem->troot, rmem->y_root, rmem->grout);
  if (gret < 0) {
    free(glo_loc); free(ghi_loc); free(sgnchg);
    return RADAU5_ROOTFN_FAIL;
  }

  /* Fill iroots: which functions crossed zero, and in which direction */
  int any_root = 0;
  for (i = 0; i < nrtfn; i++) {
    rmem->iroots[i] = 0;
    if (!rmem->gactive[i]) continue;
    /* Compare grout against the original glo (stored in rmem->glo).
     * Note: after Illinois convergence, grout may be a tiny residual with
     * the same sign as glo (not exactly zero). Use the original ghi[] to
     * determine if a crossing actually occurred in this step. */
    if (radau5_sign(rmem->ghi[i]) != radau5_sign(rmem->glo[i])) {
      int dir = rmem->rootdir[i];
      if (dir == 0 || (dir > 0 && rmem->ghi[i] > rmem->glo[i]) ||
          (dir < 0 && rmem->ghi[i] < rmem->glo[i])) {
        rmem->iroots[i] = (rmem->ghi[i] > rmem->glo[i]) ? 1 : -1;
        any_root = 1;
      }
    } else if (rmem->grout[i] == SUN_RCONST(0.0) &&
               rmem->glo[i] != SUN_RCONST(0.0)) {
      rmem->iroots[i] = (rmem->glo[i] < SUN_RCONST(0.0)) ? 1 : -1;
      any_root = 1;
    }
  }

  free(glo_loc); free(ghi_loc); free(sgnchg);

  /* If no function actually crossed zero, this was a false alarm */
  if (!any_root) return RADAU5_SUCCESS;

  return RADAU5_ROOT_RETURN;
}

/* ===========================================================================
 * radau5_root_Check1 — Initialize root function values at t0
 *
 * Called once before the first step. Evaluates g(t0, y0) to establish
 * the baseline glo[]. If any g[i]==0 at t0, probes slightly ahead to
 * determine the sign direction; if still zero, marks gactive[i]=0.
 * ===========================================================================*/
int radau5_root_Check1(Radau5Mem rmem)
{
  if (!rmem->root_active) return RADAU5_SUCCESS;
  if (rmem->root_init_done) return RADAU5_SUCCESS;

  int nrtfn = rmem->nrtfn;

  /* Evaluate g(t0, y0) */
  int ret = radau5_root_EvalG(rmem, rmem->tn, rmem->ycur, rmem->glo);
  if (ret < 0) return RADAU5_ROOTFN_FAIL;

  /* Check for any g[i] == 0 at t0 */
  int any_zero = 0;
  for (int i = 0; i < nrtfn; i++) {
    if (rmem->glo[i] == SUN_RCONST(0.0)) { any_zero = 1; break; }
  }

  if (any_zero) {
    /* Probe slightly ahead to establish sign baseline */
    sunrealtype smallh = SUNMAX(
      SUN_RCONST(100.0) * SUN_UNIT_ROUNDOFF * fabs(rmem->tn),
      SUN_RCONST(1.0e-13));
    if (rmem->h < SUN_RCONST(0.0)) smallh = -smallh;

    /* We don't have dense output yet (no step taken), so evaluate g at
       (t0, y0) — the solution hasn't changed. Instead, just mark zero
       functions as inactive. They'll be reactivated if they become nonzero
       after the first step. */
    for (int i = 0; i < nrtfn; i++) {
      if (rmem->glo[i] == SUN_RCONST(0.0)) {
        rmem->gactive[i] = 0;
      }
    }
  }

  rmem->root_init_done = 1;
  return RADAU5_SUCCESS;
}

/* ===========================================================================
 * radau5_root_Check2 — Post-root re-entry handling
 *
 * Called at the start of Radau5Solve when the previous call returned
 * RADAU5_ROOT_RETURN. Re-establishes the g baseline at the root time
 * and avoids re-detecting the same root.
 *
 * Strategy: mark root functions that fired (iroots[i] != 0) as inactive.
 * They will be reactivated in Check3 once their value becomes nonzero
 * (i.e., the solution has moved away from the root).
 * ===========================================================================*/
int radau5_root_Check2(Radau5Mem rmem)
{
  if (!rmem->root_active) return RADAU5_SUCCESS;
  if (!rmem->irfnd) return RADAU5_SUCCESS;

  int nrtfn = rmem->nrtfn;

  /* Re-evaluate g at current time (= troot from last return) */
  int ret = radau5_root_EvalG(rmem, rmem->tn, rmem->ycur, rmem->glo);
  if (ret < 0) return RADAU5_ROOTFN_FAIL;

  /* Mark fired root functions as inactive to prevent re-detection.
   * Also mark any function that is exactly zero as inactive. */
  for (int i = 0; i < nrtfn; i++) {
    if (rmem->iroots[i] != 0 || rmem->glo[i] == SUN_RCONST(0.0)) {
      rmem->gactive[i] = 0;
    }
  }

  rmem->irfnd = 0;
  return RADAU5_SUCCESS;
}

/* ===========================================================================
 * radau5_root_Check3 — Post-step sign-change detection
 *
 * Called after each accepted step. Evaluates g(tn, ycur) and checks for
 * sign changes vs glo[]. If found, calls radau5_RootFind to refine.
 *
 * IMPORTANT: glo[] must represent g values at xold (start of the current
 * accepted step), because RootFind uses Radau5Contr which is only valid
 * on [xold, tn]. We evaluate g at xold using the interpolant to ensure
 * consistency between glo[] and the bracket used by Illinois.
 * ===========================================================================*/
int radau5_root_Check3(Radau5Mem rmem)
{
  if (!rmem->root_active) return RADAU5_SUCCESS;

  int nrtfn = rmem->nrtfn;

  /* Evaluate g at the new step endpoint */
  int ret = radau5_root_EvalG(rmem, rmem->tn, rmem->ycur, rmem->ghi);
  if (ret < 0) return RADAU5_ROOTFN_FAIL;

  /* Evaluate g at xold (start of this step) using the interpolant.
   * This ensures glo[] is consistent with the Radau5Contr bracket [xold, tn].
   * On the very first step after Check1, glo was set at t0 which equals xold,
   * so this is redundant but harmless. On subsequent steps, this correctly
   * updates glo to the step start rather than carrying stale values. */
  radau5_root_InterpY(rmem, rmem->xold);
  ret = radau5_root_EvalG(rmem, rmem->xold, rmem->y_root, rmem->glo);
  if (ret < 0) return RADAU5_ROOTFN_FAIL;

  /* Reactivate any previously inactive functions that now have nonzero values */
  for (int i = 0; i < nrtfn; i++) {
    if (!rmem->gactive[i] && rmem->ghi[i] != SUN_RCONST(0.0)) {
      rmem->gactive[i] = 1;
      rmem->glo[i] = rmem->ghi[i]; /* can't detect crossing this step */
    }
  }

  /* Scan for sign changes */
  int zroot = 0;
  for (int i = 0; i < nrtfn; i++) {
    if (!rmem->gactive[i]) continue;
    if (radau5_sign(rmem->glo[i]) != radau5_sign(rmem->ghi[i])) {
      /* Check direction filter */
      int dir = rmem->rootdir[i];
      if (dir == 0 || (dir > 0 && rmem->ghi[i] > rmem->glo[i]) ||
          (dir < 0 && rmem->ghi[i] < rmem->glo[i])) {
        zroot = 1;
        break;
      }
    }
  }

  if (!zroot) {
    /* No sign change — advance baseline and continue */
    memcpy(rmem->glo, rmem->ghi, nrtfn * sizeof(sunrealtype));
    return RADAU5_SUCCESS;
  }

  /* Sign change detected — refine root location */
  ret = radau5_RootFind(rmem);
  if (ret == RADAU5_ROOT_RETURN) {
    /* Advance glo to grout (root location) for next step */
    memcpy(rmem->glo, rmem->grout, nrtfn * sizeof(sunrealtype));
    rmem->irfnd = 1;
  }
  return ret;
}

/* ===========================================================================
 * radau5_root_Free — Free all rootfinding memory
 * ===========================================================================*/
void radau5_root_Free(Radau5Mem rmem)
{
  if (!rmem) return;
  if (rmem->glo)     { free(rmem->glo);     rmem->glo     = NULL; }
  if (rmem->ghi)     { free(rmem->ghi);     rmem->ghi     = NULL; }
  if (rmem->grout)   { free(rmem->grout);   rmem->grout   = NULL; }
  if (rmem->iroots)  { free(rmem->iroots);  rmem->iroots  = NULL; }
  if (rmem->rootdir) { free(rmem->rootdir); rmem->rootdir = NULL; }
  if (rmem->gactive) { free(rmem->gactive); rmem->gactive = NULL; }
  if (rmem->y_root)  { N_VDestroy(rmem->y_root); rmem->y_root = NULL; }
  rmem->gfun = NULL;
  rmem->nrtfn = 0;
  rmem->root_active = 0;
}
