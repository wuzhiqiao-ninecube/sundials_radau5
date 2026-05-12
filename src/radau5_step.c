/* ---------------------------------------------------------------------------
 * radau5_step.c — Single step attempt for the RADAU5 solver
 *
 * Implements radau5_Step(), which performs ONE complete step attempt,
 * looping internally until the step is accepted or a fatal error occurs.
 *
 * Translates the Fortran RADCOR main loop (labels 10-78) from:
 *   Hairer & Wanner, "Solving ODEs II", radau5.f lines 793-1101.
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <sundials/sundials_math.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include "radau5_impl.h"

int radau5_Step(Radau5Mem rmem)
{
  sunindextype n    = rmem->n;
  sunrealtype  uround = SUN_UNIT_ROUNDOFF;
  int ns = rmem->ns;

  /* Collocation node offsets */
  sunrealtype c1   = rmem->c[0];
  sunrealtype c2   = rmem->c[1];
  sunrealtype c1m1 = (rmem->c[0] - SUN_RCONST(1.0));
  sunrealtype c2m1 = (rmem->c[1] - SUN_RCONST(1.0));
  (void)c1m1; /* used only in extrapolation via cont arrays */
  (void)c2m1;

  /* Eigenvalue parameters */
  sunrealtype u1   = rmem->u1;
  sunrealtype alph = rmem->alph[0];
  sunrealtype beta = rmem->beta_eig[0];

  sunrealtype TI11,TI12,TI13, TI21,TI22,TI23, TI31,TI32,TI33; /* set below based on use_schur */

  /* TI matrix rows (for F1/F2/F3 from Z extrapolation) */
  if (rmem->use_schur) {
    TI11 = rmem->US_mat[0], TI12 = rmem->US_mat[3], TI13 = rmem->US_mat[6];
    TI21 = rmem->US_mat[1], TI22 = rmem->US_mat[4], TI23 = rmem->US_mat[7];
    TI31 = rmem->US_mat[2], TI32 = rmem->US_mat[5], TI33 = rmem->US_mat[8];
  }else {
    TI11 = rmem->TI_mat[0], TI12 = rmem->TI_mat[1], TI13 = rmem->TI_mat[2];
    TI21 = rmem->TI_mat[3], TI22 = rmem->TI_mat[4], TI23 = rmem->TI_mat[5];
    TI31 = rmem->TI_mat[6], TI32 = rmem->TI_mat[7], TI33 = rmem->TI_mat[8];
  }

  /* Step-size control parameters */
  sunrealtype safe  = rmem->safe;
  sunrealtype facl  = rmem->facl;
  sunrealtype facr  = rmem->facr;
  sunrealtype thet  = rmem->thet;
  sunrealtype quot1 = rmem->quot1;
  sunrealtype quot2 = rmem->quot2;

  /* cfac = safe * (1 + 2*nit), used in step-size formula */
  sunrealtype cfac  = safe * (SUN_RCONST(1.0) + SUN_RCONST(2.0) * rmem->nit);

  /* Raw data pointers — set inside loops */
  sunrealtype *z1d, *z2d, *z3d;
  sunrealtype *f1d, *f2d, *f3d;
  sunrealtype *scald, *ycurd;

  int ret;
  int newt;
  sunrealtype err;
  sunrealtype fac1, alphn, betan;
  sunrealtype fac, quot, hnew, qt;
  sunrealtype facgus;

  /* =========================================================================
   * Entry: decide where to start based on caljac/skipdecomp flags.
   * skipdecomp==1: skip both Jacobian and decomposition (goto label30)
   *   — only valid when h hasn't changed (qt in [quot1,quot2])
   * caljac==0: reuse Jacobian, but refactor (goto label20)
   * caljac==1: recompute Jacobian and refactor (goto label10)
   * On the very first call, caljac==1 (set in Radau5Init).
   * =======================================================================*/
  if (rmem->skipdecomp)
  {
    rmem->skipdecomp = 0;
    /* Recompute fac1, alphn, betan for the (unchanged) h */
    fac1  = u1   / rmem->h;
    alphn = alph / rmem->h;
    betan = beta / rmem->h;
    goto label30;
  }
  if (!rmem->caljac)
    goto label20;

  /* =========================================================================
   * label10: Compute Jacobian
   * =======================================================================*/
label10:
  if (rmem->jac != NULL)
  {
    ret = rmem->jac(rmem->tn, rmem->ycur, rmem->fn, rmem->J,
                    rmem->user_data, rmem->tmp1, rmem->tmp2, rmem->tmp3);
    if (ret != 0) return RADAU5_RHSFUNC_FAIL;
  }
  else
  {
    if (rmem->mat_id == SUNMATRIX_SPARSE)
    {
      if (rmem->sp_colptrs == NULL)
        return RADAU5_ILL_INPUT;  /* sparse requires analytic Jac or sparsity pattern */
      ret = radau5_DQJacSparse(rmem, rmem->tn, rmem->ycur, rmem->fn);
      if (ret == RADAU5_RHSFUNC_RECVR) goto label78;
      if (ret != RADAU5_SUCCESS) return ret;
    }
    else if (rmem->mat_id == SUNMATRIX_BAND)
      ret = radau5_DQJacBand(rmem, rmem->tn, rmem->ycur, rmem->fn);
    else
      ret = radau5_DQJacDense(rmem, rmem->tn, rmem->ycur, rmem->fn);
    if (ret == RADAU5_RHSFUNC_RECVR) goto label78;
    if (ret != RADAU5_SUCCESS) return ret;
  }
  rmem->njac++;
  rmem->caljac = 1;

  /* =========================================================================
   * label20: Build and factor E1, E2
   * =======================================================================*/
label20:
  fac1  = u1   / rmem->h;
  alphn = alph / rmem->h;
  betan = beta / rmem->h;

  radau5_BuildE1(rmem, fac1);
  radau5_BuildE2(rmem, alphn, betan);

  ret = radau5_DecompE1(rmem);
  if (ret != RADAU5_SUCCESS) goto label78;

  ret = radau5_DecompE2(rmem);
  if (ret != RADAU5_SUCCESS) goto label78;

  /* =========================================================================
   * label30: Newton starting values, Newton iteration, error estimate
   * =======================================================================*/
label30:
  rmem->nstep++;

  if (rmem->nstep > rmem->mxstep)
    return RADAU5_TOO_MANY_STEPS;

  if (SUN_RCONST(0.1) * fabs(rmem->h) <= fabs(rmem->tn) * uround)
    return RADAU5_STEP_TOO_SMALL;

  /* DAE index-2/3 scaling of scal */
  if (rmem->nind2 > 0)
  {
    scald = N_VGetArrayPointer(rmem->scal);
    for (sunindextype i = rmem->nind1; i < rmem->nind1 + rmem->nind2; i++)
      scald[i] /= rmem->hhfac;
  }
  if (rmem->nind3 > 0)
  {
    scald = N_VGetArrayPointer(rmem->scal);
    for (sunindextype i = rmem->nind1 + rmem->nind2;
         i < rmem->nind1 + rmem->nind2 + rmem->nind3; i++)
      scald[i] /= (rmem->hhfac * rmem->hhfac);
  }

  /* Newton starting values */
  z1d = N_VGetArrayPointer(rmem->z[0]);
  z2d = N_VGetArrayPointer(rmem->z[1]);
  z3d = N_VGetArrayPointer(rmem->z[2]);
  f1d = N_VGetArrayPointer(rmem->f[0]);
  f2d = N_VGetArrayPointer(rmem->f[1]);
  f3d = N_VGetArrayPointer(rmem->f[2]);

  if (rmem->first || rmem->startn)
  {
    for (sunindextype i = 0; i < n; i++)
    {
      z1d[i] = SUN_RCONST(0.0);
      z2d[i] = SUN_RCONST(0.0);
      z3d[i] = SUN_RCONST(0.0);
      f1d[i] = SUN_RCONST(0.0);
      f2d[i] = SUN_RCONST(0.0);
      f3d[i] = SUN_RCONST(0.0);
    }
  }
  else
  {
    /* Extrapolate from previous step using continuous output coefficients */
    sunrealtype *ak1d = N_VGetArrayPointer(rmem->cont[1]);
    sunrealtype *ak2d = N_VGetArrayPointer(rmem->cont[2]);
    sunrealtype *ak3d = N_VGetArrayPointer(rmem->cont[3]);

    sunrealtype c3q = rmem->h / rmem->hold;
    sunrealtype c1q = c1 * c3q;
    sunrealtype c2q = c2 * c3q;

    for (sunindextype i = 0; i < n; i++)
    {
      sunrealtype ak1 = ak1d[i];
      sunrealtype ak2 = ak2d[i];
      sunrealtype ak3 = ak3d[i];

      sunrealtype z1i = c1q * (ak1 + (c1q - c2m1) * (ak2 + (c1q - c1m1) * ak3));
      sunrealtype z2i = c2q * (ak1 + (c2q - c2m1) * (ak2 + (c2q - c1m1) * ak3));
      sunrealtype z3i = c3q * (ak1 + (c3q - c2m1) * (ak2 + (c3q - c1m1) * ak3));

      z1d[i] = z1i;
      z2d[i] = z2i;
      z3d[i] = z3i;

      f1d[i] = TI11 * z1i + TI12 * z2i + TI13 * z3i;
      f2d[i] = TI21 * z1i + TI22 * z2i + TI23 * z3i;
      f3d[i] = TI31 * z1i + TI32 * z2i + TI33 * z3i;
    }
  }

  /* Newton iteration */
  ret = radau5_Newton(rmem, &newt);
  if (ret == RADAU5_NEWT_PREDICT)
  {
    /* Newton predicted slow convergence and already reduced h.
     * Go directly to refactor (like Fortran: IF(CALJAC) GOTO 20; GOTO 10) */
    if (rmem->naccpt >= 1) rmem->nrejct++;
    if (rmem->caljac) goto label20;
    goto label10;
  }
  if (ret != RADAU5_SUCCESS)
    goto label78;

  /* Error estimate */
  ret = radau5_ErrorEstimate(rmem, &err);
  if (ret == RADAU5_RHSFUNC_RECVR) goto label78;
  if (ret != RADAU5_SUCCESS) return ret;

  /* Compute proposed new step size */
  fac  = SUNMIN(safe, cfac / (sunrealtype)(newt + 2 * rmem->nit));
  quot = SUNMAX(facr, SUNMIN(facl, SUNRpowerR(err, SUN_RCONST(1.0) / (sunrealtype)(ns + 1)) / fac));
  hnew = rmem->h / quot;


  /* =========================================================================
   * Accept or reject
   * =======================================================================*/
  if (err < SUN_RCONST(1.0))
  {
    /* --- Step accepted --- */
    rmem->first  = 0;
    rmem->naccpt++;

    /* Gustafsson predictive controller (pred == 1) */
    if (rmem->pred == 1 && rmem->naccpt > 1)
    {
      facgus = (rmem->hacc / rmem->h)
             * SUNRpowerR(err * err / rmem->erracc, SUN_RCONST(1.0) / (sunrealtype)(ns + 1))
             / safe;
      facgus = SUNMAX(facr, SUNMIN(facl, facgus));
      quot   = SUNMAX(quot, facgus);
      hnew   = rmem->h / quot;
    }
    rmem->hacc   = rmem->h;
    rmem->erracc = SUNMAX(SUN_RCONST(0.01), err);

    /* Advance time and solution */
    rmem->xold = rmem->tn;
    rmem->hold = rmem->h;
    rmem->tn  += rmem->h;

    ycurd = N_VGetArrayPointer(rmem->ycur);
    z3d   = N_VGetArrayPointer(rmem->z[2]);
    for (sunindextype i = 0; i < n; i++)
      ycurd[i] += z3d[i];

    /* Update continuous output coefficients */
    radau5_UpdateContinuousOutput(rmem);

    /* Recompute error weights at new y */
    radau5_ComputeScal(rmem, rmem->ycur);

    /* Evaluate f at new (tn, ycur) */
    ret = rmem->rhs(rmem->tn, rmem->ycur, rmem->fn, rmem->user_data);
    if (ret != 0) return RADAU5_RHSFUNC_FAIL;
    rmem->nfcn++;

    /* Bound hnew and record hopt */
    sunrealtype posneg = (rmem->h >= SUN_RCONST(0.0))
                         ? SUN_RCONST(1.0) : SUN_RCONST(-1.0);
    sunrealtype hmaxn  = (rmem->hmax > SUN_RCONST(0.0))
                         ? rmem->hmax : fabs(rmem->h) * SUN_RCONST(1.0e4);
    hnew = posneg * SUNMIN(fabs(hnew), hmaxn);
    rmem->hopt = SUNMIN(fabs(rmem->h), fabs(hnew));

    /* If previous step was rejected, don't grow beyond |h| */
    if (rmem->reject)
      hnew = posneg * SUNMIN(fabs(hnew), fabs(rmem->h));

    rmem->reject = 0;
    rmem->caljac = 0;
    rmem->nsing  = 0;  /* reset singular counter on success */

    /* Decide whether to reuse Jacobian / factorization on NEXT call.
     * We store the decision but always return to Radau5Solve so it can
     * check tout, call solout, etc.  The goto labels are only used for
     * retry loops within a single step attempt (reject / Newton fail). */
    qt = hnew / rmem->h;
    rmem->hhfac = rmem->h;

    if (rmem->theta <= thet && qt >= quot1 && qt <= quot2)
    {
      /* Reuse both Jacobian and factorization (Fortran "goto 30").
       * Set h=hnew — the factorization is slightly stale but qt ∈ [1,6]
       * means hnew ≈ h, so the mismatch is small. */
      rmem->h = hnew;
      rmem->caljac = 0;
      rmem->skipdecomp = 1;
      return RADAU5_SUCCESS;
    }

    rmem->h     = hnew;
    rmem->hhfac = rmem->h;

    if (rmem->theta <= thet)
    {
      rmem->caljac = 0;  /* reuse Jacobian, but refactor next time */
    }
    else
    {
      rmem->caljac = 1;  /* recompute Jacobian next time */
    }

    return RADAU5_SUCCESS;
  }
  else
  {
    /* --- Step rejected --- */
    rmem->reject = 1;
    rmem->last   = 0;

    if (rmem->first)
    {
      rmem->h    *= SUN_RCONST(0.1);
      rmem->hhfac = SUN_RCONST(0.1);
    }
    else
    {
      rmem->hhfac = hnew / rmem->h;
      rmem->h     = hnew;
    }

    if (rmem->naccpt >= 1) rmem->nrejct++;

    if (rmem->caljac) goto label20;
    goto label10;
  }

  /* =========================================================================
   * label78: Unexpected step rejection (Newton failure, singular matrix,
   *          or recoverable RHS error).
   *
   * Fortran radau5.f label 78 (lines 1108-1124):
   *   - IER != 0 (singular matrix): NSING++, if >= 5 → fatal
   *   - IERR != 0 (RHS error):      NSING++, if >= 5 → fatal, reset IERR
   *   - Otherwise (convergence failure): NSING NOT incremented
   *   Then: H *= 0.5, retry.
   *
   * The key insight: only singular matrix and RHS errors count toward the
   * nsing limit. Newton convergence failures just halve h and retry
   * indefinitely (bounded only by step-too-small check).
   * =======================================================================*/
label78:
  if (ret == RADAU5_SINGULAR_MATRIX || ret == RADAU5_LSETUP_FAIL)
  {
    rmem->nsing++;
    if (rmem->nsing >= 5)
      return RADAU5_SINGULAR_MATRIX;
  }
  else if (ret == RADAU5_RHSFUNC_RECVR)
  {
    rmem->nsing++;
    if (rmem->nsing >= 5)
      return RADAU5_RHSFUNC_FAIL;
  }
  /* CONV_FAILURE / NEWT_PREDICT: do NOT increment nsing */

  rmem->h    *= SUN_RCONST(0.5);
  rmem->hhfac = SUN_RCONST(0.5);
  rmem->reject = 1;
  rmem->last   = 0;

  if (rmem->caljac) goto label20;
  goto label10;
}
