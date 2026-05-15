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

  /* Eigenvalue parameters */
  sunrealtype u1   = rmem->u1;
  sunrealtype alph = rmem->alph[0];
  sunrealtype beta = rmem->beta_eig[0];

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
  sunrealtype *scald, *ycurd;

  int ret;
  int newt = 0;
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

    /* Variable-order: ikeep path — if enough steps passed, try order increase */
    if (rmem->variab && rmem->ikeep)
    {
      rmem->ichan++;
      rmem->ikeep = 0;
      if (rmem->ichan >= 10 && rmem->ns < rmem->nsmax)
        goto label20;  /* retry with potential order increase */
    }

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
   * Variable-order selection (Fortran radau.f lines 900-935)
   * Placed before label20 so that the new ns is used for E1/E2 build.
   * =======================================================================*/
label20:
  if (rmem->variab)
  {
    rmem->ichan++;
    sunrealtype hquot = rmem->h / rmem->hold;
    rmem->thetat = SUNMIN(SUN_RCONST(10.0),
                          SUNMAX(rmem->theta, rmem->thetat * SUN_RCONST(0.5)));
    rmem->nsnew = rmem->ns;

    /* Consider increasing order */
    if (rmem->newt_prev > 1 && rmem->thetat <= rmem->vitu
        && hquot < rmem->hhou && hquot > rmem->hhod)
      rmem->nsnew = SUNMIN(rmem->nsmax, rmem->ns + 2);

    /* Consider decreasing order */
    if (rmem->thetat >= rmem->vitd || rmem->unexp)
      rmem->nsnew = SUNMAX(rmem->nsmin, rmem->ns - 2);
    if (rmem->ichan >= 1 && rmem->unexn)
      rmem->nsnew = SUNMAX(rmem->nsmin, rmem->ns - 2);

    /* Don't increase too soon after a change */
    if (rmem->ichan <= 10)
      rmem->nsnew = SUNMIN(rmem->ns, rmem->nsnew);

    rmem->change = (rmem->ns != rmem->nsnew) ? 1 : 0;
    rmem->unexn = 0;
    rmem->unexp = 0;

    if (rmem->change)
    {
      ret = radau5_ChangeOrder(rmem, rmem->nsnew);
      if (ret != RADAU5_SUCCESS) return ret;
      rmem->ichan = 1;
      /* Update local variables that depend on ns */
      ns = rmem->ns;
      u1 = rmem->u1;
      alph = rmem->alph[0];
      beta = rmem->beta_eig[0];
    }
  }

  /* =========================================================================
   * Build and factor E1, E2
   * =======================================================================*/
  fac1  = u1   / rmem->h;
  alphn = alph / rmem->h;
  betan = beta / rmem->h;

  radau5_BuildE1(rmem, fac1);
  for (int pk = 0; pk < rmem->npairs; pk++) {
    sunrealtype alphn_k = rmem->alph[pk] / rmem->h;
    sunrealtype betan_k = rmem->beta_eig[pk] / rmem->h;
    radau5_BuildE2(rmem, pk, alphn_k, betan_k);
  }

  ret = radau5_DecompE1(rmem);
  if (ret != RADAU5_SUCCESS) goto label78;

  for (int pk = 0; pk < rmem->npairs; pk++) {
    ret = radau5_DecompE2(rmem, pk);
    if (ret != RADAU5_SUCCESS) goto label78;
  }
  ret = RADAU5_SUCCESS;

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
  if (rmem->first || rmem->startn || rmem->change)
  {
    for (int k = 0; k < ns; k++)
    {
      sunrealtype* zd = N_VGetArrayPointer(rmem->z[k]);
      sunrealtype* fd = N_VGetArrayPointer(rmem->f[k]);
      for (sunindextype i = 0; i < n; i++)
      {
        zd[i] = SUN_RCONST(0.0);
        fd[i] = SUN_RCONST(0.0);
      }
    }
  }
  else
  {
    /* Extrapolate from previous step using continuous output coefficients.
     *
     * General algorithm from Fortran radau.f (ns=5/7 blocks):
     *   HQUOT = H/HOLD
     *   DO K=1,NS
     *     CCQ = C(K)*HQUOT
     *     DO I=1,N
     *       VAL = CONT(I+NS*N)
     *       DO L=NS-1,1,-1
     *         VAL = CONT(I+L*N) + (CCQ - C(NS-L) + 1) * VAL
     *       END DO
     *       ZZ(I+(K-1)*N) = CCQ * VAL
     *     END DO
     *   END DO
     *
     * Node mapping: Fortran C(NS-L) → c[ns-l-1] for l=1..ns-1.
     * The "+1" comes from the shift s = (t-xsol)/hsol + 1 in CONTRA;
     * here CCQ = c[k]*hquot plays the role of s for the new collocation point.
     */
    sunrealtype hquot = rmem->h / rmem->hold;

    /* Cache cont pointers outside inner loops */
    sunrealtype* contd[RADAU5_NS_MAX + 1];
    for (int l = 0; l <= ns; l++)
      contd[l] = N_VGetArrayPointer(rmem->cont[l]);

    for (int k = 0; k < ns; k++)
    {
      sunrealtype ccq = rmem->c[k] * hquot;
      sunrealtype* zd = N_VGetArrayPointer(rmem->z[k]);

      for (sunindextype i = 0; i < n; i++)
      {
        sunrealtype val = contd[ns][i];
        for (int l = ns - 1; l >= 1; l--)
          val = contd[l][i] + (ccq - rmem->c[ns - l - 1] + SUN_RCONST(1.0)) * val;
        zd[i] = ccq * val;
      }
    }

    /* Compute f = TI * z (forward transform for Newton starting values) */
    sunrealtype *zd_arr[RADAU5_NS_MAX], *fd_arr[RADAU5_NS_MAX];
    for (int k = 0; k < ns; k++) {
      zd_arr[k] = N_VGetArrayPointer(rmem->z[k]);
      fd_arr[k] = N_VGetArrayPointer(rmem->f[k]);
    }

    for (sunindextype i = 0; i < n; i++)
    {
      sunrealtype zvals[RADAU5_NS_MAX];
      for (int k = 0; k < ns; k++)
        zvals[k] = zd_arr[k][i];

      for (int k = 0; k < ns; k++)
      {
        sunrealtype sum = SUN_RCONST(0.0);
        if (rmem->use_schur) {
          for (int j = 0; j < ns; j++)
            sum += rmem->US_mat[j * ns + k] * zvals[j];
        } else {
          for (int j = 0; j < ns; j++)
            sum += rmem->TI_mat[k * ns + j] * zvals[j];
        }
        fd_arr[k][i] = sum;
      }
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
    int order_changed = rmem->change;
    rmem->first  = 0;
    rmem->change = 0;
    rmem->newt_prev = newt;
    rmem->naccpt++;

    /* Gustafsson predictive controller (pred == 1, skip if order just changed) */
    if (rmem->pred == 1 && rmem->naccpt > 1 && !order_changed)
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
    sunrealtype* zlast = N_VGetArrayPointer(rmem->z[ns-1]);
    for (sunindextype i = 0; i < n; i++)
      ycurd[i] += zlast[i];

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
      rmem->ikeep = 1;  /* mark for variable-order ikeep logic */
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

  rmem->unexp = 1;  /* mark unexpected rejection for order selection */
  rmem->h    *= SUN_RCONST(0.5);
  rmem->hhfac = SUN_RCONST(0.5);
  rmem->reject = 1;
  rmem->last   = 0;

  if (rmem->caljac) goto label20;
  goto label10;
}
