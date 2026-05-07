/* ---------------------------------------------------------------------------
 * RADAU5 — simplified Newton iteration
 *
 * Implements radau5_Newton(), which performs the inner Newton loop for the
 * 3-stage Radau IIA collocation system.  E1 and E2 must already be factored
 * before this routine is called.
 *
 * Reference: Hairer & Wanner, "Solving ODEs II", Fortran RADCOR lines 906-985
 *            and subroutine SLVRAD / SLVRAI.
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <sundials/sundials_math.h>
#include <nvector/nvector_serial.h>
#include "radau5_impl.h"

int radau5_Newton(Radau5Mem rmem, int* newt_out)
{
  /* -------------------------------------------------------------------------
   * Unpack frequently used fields
   * -----------------------------------------------------------------------*/
  sunindextype n        = rmem->n;
  sunrealtype  h        = rmem->h;
  sunrealtype  tn       = rmem->tn;
  sunrealtype  uround   = SUN_UNIT_ROUNDOFF;

  /* Collocation nodes */
  sunrealtype c1 = rmem->c1;
  sunrealtype c2 = rmem->c2;

  /* Eigenvalue scalings */
  sunrealtype u1   = rmem->u1;
  sunrealtype alph = rmem->alph;
  sunrealtype beta = rmem->beta;

  /* Transformation matrices */
  sunrealtype TI11 = rmem->TI11, TI12 = rmem->TI12, TI13 = rmem->TI13;
  sunrealtype TI21 = rmem->TI21, TI22 = rmem->TI22, TI23 = rmem->TI23;
  sunrealtype TI31 = rmem->TI31, TI32 = rmem->TI32, TI33 = rmem->TI33;

  sunrealtype T11 = rmem->T11, T12 = rmem->T12, T13 = rmem->T13;
  sunrealtype T21 = rmem->T21, T22 = rmem->T22, T23 = rmem->T23;
  sunrealtype T31 = rmem->T31; /* T32 = 1, T33 = 0 (implicit) */

  /* Newton control */
  int          nit   = rmem->nit;
  sunrealtype  fnewt = rmem->fnewt;

  /* Working vectors */
  N_Vector ycur = rmem->ycur;
  N_Vector z1   = rmem->z1;
  N_Vector z2   = rmem->z2;
  N_Vector z3   = rmem->z3;
  N_Vector f1   = rmem->f1;
  N_Vector f2   = rmem->f2;
  N_Vector f3   = rmem->f3;
  N_Vector tmp1 = rmem->tmp1;
  N_Vector scal = rmem->scal;
  N_Vector rhs2 = rmem->rhs2;
  N_Vector sol2 = rmem->sol2;

  /* Raw data pointers (set inside loop) */
  sunrealtype *z1_data, *z2_data, *z3_data;
  sunrealtype *f1_data, *f2_data, *f3_data;
  sunrealtype *scal_data, *rhs2_data, *sol2_data;

  /* Convergence tracking */
  sunrealtype faccon  = SUNMAX(rmem->faccon, uround);
  faccon = SUNRpowerR(faccon, SUN_RCONST(0.8));
  sunrealtype theta   = SUNRabs(rmem->thet);
  sunrealtype dynold  = SUN_RCONST(1.0);   /* initialised before first use */
  sunrealtype thqold  = SUN_RCONST(0.0);
  sunrealtype dyno    = SUN_RCONST(0.0);

  int newt = 0;

  /* =========================================================================
   * Newton loop (Fortran label 40)
   * =======================================================================*/
  while (1)
  {
    if (newt >= nit)
    {
      *newt_out = newt;
      return RADAU5_CONV_FAILURE;
    }

    /* ------------------------------------------------------------------
     * 1. Evaluate RHS at the three collocation points.
     *    On entry z1,z2,z3 hold the current stage increments Z_i.
     *    After the calls they hold f(t+c_i*h, y+Z_i).
     * ----------------------------------------------------------------*/
    /* z1 <- f(t+c1*h, y+Z1) */
    N_VLinearSum(SUN_RCONST(1.0), ycur, SUN_RCONST(1.0), z1, tmp1);
    {
      int rhsret = rmem->rhs(tn + c1 * h, tmp1, z1, rmem->user_data);
      if (rhsret < 0) return RADAU5_RHSFUNC_FAIL;
      if (rhsret > 0) return RADAU5_RHSFUNC_RECVR;
    }

    /* z2 <- f(t+c2*h, y+Z2) */
    N_VLinearSum(SUN_RCONST(1.0), ycur, SUN_RCONST(1.0), z2, tmp1);
    {
      int rhsret = rmem->rhs(tn + c2 * h, tmp1, z2, rmem->user_data);
      if (rhsret < 0) return RADAU5_RHSFUNC_FAIL;
      if (rhsret > 0) return RADAU5_RHSFUNC_RECVR;
    }

    /* z3 <- f(t+h, y+Z3) */
    N_VLinearSum(SUN_RCONST(1.0), ycur, SUN_RCONST(1.0), z3, tmp1);
    {
      int rhsret = rmem->rhs(tn + h, tmp1, z3, rmem->user_data);
      if (rhsret < 0) return RADAU5_RHSFUNC_FAIL;
      if (rhsret > 0) return RADAU5_RHSFUNC_RECVR;
    }

    rmem->nfcn += 3;

    /* ------------------------------------------------------------------
     * 2. Apply TI transform in-place.
     *    z1,z2,z3 currently hold rhs values a1,a2,a3.
     *    After: z_i = sum_j TI_{ij} * a_j
     * ----------------------------------------------------------------*/
    z1_data = N_VGetArrayPointer(z1);
    z2_data = N_VGetArrayPointer(z2);
    z3_data = N_VGetArrayPointer(z3);

    for (sunindextype i = 0; i < n; i++)
    {
      sunrealtype a1 = z1_data[i];
      sunrealtype a2 = z2_data[i];
      sunrealtype a3 = z3_data[i];
      z1_data[i] = TI11 * a1 + TI12 * a2 + TI13 * a3;
      z2_data[i] = TI21 * a1 + TI22 * a2 + TI23 * a3;
      z3_data[i] = TI31 * a1 + TI32 * a2 + TI33 * a3;
    }

    /* ------------------------------------------------------------------
     * 3. Form linear-system RHS.
     *
     * Eigenvalue mode (use_schur==0):
     *    For identity mass (M==NULL):
     *      z1(i) -= fac1 * f1(i)
     *      z2(i) += (-f2(i))*alphn - (-f3(i))*betan
     *      z3(i) += (-f3(i))*alphn + (-f2(i))*betan
     *
     *    For general mass M (Fortran SLVRAD IJOB=3,5):
     *      S_k = -(M * F_k)  for k=1,2,3
     *      z1(i) += fac1  * S1(i)
     *      z2(i) += alphn * S2(i) - betan * S3(i)
     *      z3(i) += alphn * S3(i) + betan * S2(i)
     *
     * Schur mode (use_schur==1):
     *    rhs1 = w1 - (TS[0][0]/h)*M*F1 - (TS[0][1]/h)*M*F2 - (TS[0][2]/h)*M*F3
     *    rhs2 = w2 - (TS[1][0]/h)*M*F1 - (TS[1][1]/h)*M*F2 - (TS[1][2]/h)*M*F3
     *    rhs3 = w3                                           - (TS[2][2]/h)*M*F3
     *    Then block back-substitution:
     *      Solve E1 * dF3 = rhs3
     *      rhs1 -= (TS[0][2]/h)*M*dF3;  rhs2 -= (TS[1][2]/h)*M*dF3
     *      Solve E2 * [dF1; dF2] = [rhs1; rhs2]
     * ----------------------------------------------------------------*/
    sunrealtype fac1   = u1   / h;
    sunrealtype alphn  = alph / h;
    sunrealtype betan  = beta / h;

    f1_data = N_VGetArrayPointer(f1);
    f2_data = N_VGetArrayPointer(f2);
    f3_data = N_VGetArrayPointer(f3);

    if (rmem->use_schur)
    {
      /* --- Schur mode: full TS row scalings + block back-substitution --- */
      sunrealtype ts00 = rmem->TS[0][0] / h;
      sunrealtype ts01 = rmem->TS[0][1] / h;
      sunrealtype ts02 = rmem->TS[0][2] / h;
      sunrealtype ts10 = rmem->TS[1][0] / h;
      sunrealtype ts11 = rmem->TS[1][1] / h;
      sunrealtype ts12 = rmem->TS[1][2] / h;
      sunrealtype ts22 = rmem->TS[2][2] / h;  /* = fac1 = u1/h */

      if (rmem->M != NULL)
      {
        /* General mass: M * F_k products needed.
         * We compute M*f3, M*f2, M*f1 sequentially, building the RHS. */
        sunrealtype *mfd;

        /* --- rhs3 = w3 - ts22*M*F3 --- */
        radau5_MassMult(rmem, f3, tmp1);
        mfd = N_VGetArrayPointer(tmp1);
        for (sunindextype ii = 0; ii < n; ii++)
          z3_data[ii] -= ts22 * mfd[ii];

        /* Save M*F3 for the ts02, ts12 terms */
        sunrealtype mf3_buf[1024];
        sunrealtype *mf3_ptr = (n <= 1024) ? mf3_buf
                              : (sunrealtype*)malloc((size_t)n * sizeof(sunrealtype));
        for (sunindextype ii = 0; ii < n; ii++)
          mf3_ptr[ii] = mfd[ii];

        /* --- M*F2 -> apply to rhs1, rhs2 --- */
        radau5_MassMult(rmem, f2, tmp1);
        mfd = N_VGetArrayPointer(tmp1);
        sunrealtype mf2_buf[1024];
        sunrealtype *mf2_ptr = (n <= 1024) ? mf2_buf
                              : (sunrealtype*)malloc((size_t)n * sizeof(sunrealtype));
        for (sunindextype ii = 0; ii < n; ii++)
          mf2_ptr[ii] = mfd[ii];

        /* --- M*F1 -> apply to rhs1, rhs2 --- */
        radau5_MassMult(rmem, f1, tmp1);
        mfd = N_VGetArrayPointer(tmp1);

        for (sunindextype ii = 0; ii < n; ii++)
        {
          z1_data[ii] -= ts00 * mfd[ii] + ts01 * mf2_ptr[ii] + ts02 * mf3_ptr[ii];
          z2_data[ii] -= ts10 * mfd[ii] + ts11 * mf2_ptr[ii] + ts12 * mf3_ptr[ii];
        }

        if (n > 1024) { free(mf3_ptr); free(mf2_ptr); }
      }
      else
      {
        /* Identity mass: M = I, so M*F_k = F_k */
        for (sunindextype ii = 0; ii < n; ii++)
        {
          z1_data[ii] -= ts00 * f1_data[ii] + ts01 * f2_data[ii] + ts02 * f3_data[ii];
          z2_data[ii] -= ts10 * f1_data[ii] + ts11 * f2_data[ii] + ts12 * f3_data[ii];
          z3_data[ii] -= ts22 * f3_data[ii];
        }
      }

      /* ------------------------------------------------------------------
       * Block back-substitution:
       * Step A: Solve E1 * dF3 = rhs3 (z3 holds rhs3, result in z3)
       * ----------------------------------------------------------------*/
      if (SUNLinSolSolve(rmem->LS_E1, rmem->E1, z3, z3, SUN_RCONST(0.0)) != 0)
        return RADAU5_LSOLVE_FAIL;

      /* Step B: rhs1 -= ts02*M*dF3;  rhs2 -= ts12*M*dF3
       * z3 now holds dF3 */
      if (rmem->M != NULL)
      {
        radau5_MassMult(rmem, z3, tmp1);
        sunrealtype *mdf3 = N_VGetArrayPointer(tmp1);
        for (sunindextype ii = 0; ii < n; ii++)
        {
          z1_data[ii] -= ts02 * mdf3[ii];
          z2_data[ii] -= ts12 * mdf3[ii];
        }
      }
      else
      {
        for (sunindextype ii = 0; ii < n; ii++)
        {
          z1_data[ii] -= ts02 * z3_data[ii];
          z2_data[ii] -= ts12 * z3_data[ii];
        }
      }

      /* Step C: Solve E2 * [dF1; dF2] = [rhs1; rhs2] */
      rhs2_data = N_VGetArrayPointer(rhs2);
      for (sunindextype ii = 0; ii < n; ii++)
      {
        rhs2_data[ii]     = z1_data[ii];
        rhs2_data[ii + n] = z2_data[ii];
      }

      if (SUNLinSolSolve(rmem->LS_E2, rmem->E2, sol2, rhs2, SUN_RCONST(0.0)) != 0)
        return RADAU5_LSOLVE_FAIL;

      sol2_data = N_VGetArrayPointer(sol2);
      for (sunindextype ii = 0; ii < n; ii++)
      {
        z1_data[ii] = sol2_data[ii];
        z2_data[ii] = sol2_data[ii + n];
      }

      /* z3 already holds dF3 from Step A */

      rmem->nsol++;
      newt++;
    }
    else
    {
      /* --- Eigenvalue mode (original) --- */
      if (rmem->M != NULL)
      {
        /* General mass matrix: compute S_k = -(M * F_k) using tmp vectors. */
        sunrealtype *s1d, *s2d, *s3d;

        /* S1 = -(M * f1) */
        radau5_MassMult(rmem, f1, tmp1);
        s1d = N_VGetArrayPointer(tmp1);

        /* Apply S1 to z1 */
        for (sunindextype ii = 0; ii < n; ii++)
          z1_data[ii] += fac1 * (-s1d[ii]);

        /* S2 = -(M * f2) → store in tmp1 */
        radau5_MassMult(rmem, f2, tmp1);
        s2d = N_VGetArrayPointer(tmp1);

        sunrealtype s2_buf[1024];
        sunrealtype *s2_ptr = (n <= 1024) ? s2_buf
                              : (sunrealtype*)malloc((size_t)n * sizeof(sunrealtype));
        for (sunindextype ii = 0; ii < n; ii++)
          s2_ptr[ii] = -s2d[ii];

        /* S3 = -(M * f3) → store in tmp1 */
        radau5_MassMult(rmem, f3, tmp1);
        s3d = N_VGetArrayPointer(tmp1);

        for (sunindextype ii = 0; ii < n; ii++)
        {
          sunrealtype s2_i = s2_ptr[ii];
          sunrealtype s3_i = -s3d[ii];
          z2_data[ii] += alphn * s2_i - betan * s3_i;
          z3_data[ii] += alphn * s3_i + betan * s2_i;
        }

        if (n > 1024) free(s2_ptr);
      }
      else
      {
        /* Identity mass: S_k = -F_k */
        for (sunindextype ii = 0; ii < n; ii++)
        {
          z1_data[ii] -= fac1 * f1_data[ii];

          sunrealtype s2 = -f2_data[ii];
          sunrealtype s3 = -f3_data[ii];
          z2_data[ii] += s2 * alphn - s3 * betan;
          z3_data[ii] += s3 * alphn + s2 * betan;
        }
      }

      /* ------------------------------------------------------------------
       * 4. Solve E1 * dz1 = z1  (result overwrites z1)
       * ----------------------------------------------------------------*/
      if (SUNLinSolSolve(rmem->LS_E1, rmem->E1, z1, z1, SUN_RCONST(0.0)) != 0)
        return RADAU5_LSOLVE_FAIL;

      /* ------------------------------------------------------------------
       * 5. Solve E2 * [dz2; dz3] = [z2; z3]
       *    Pack into the 2n vector rhs2, solve into sol2, then unpack.
       * ----------------------------------------------------------------*/
      rhs2_data = N_VGetArrayPointer(rhs2);
      for (sunindextype ii = 0; ii < n; ii++)
      {
        rhs2_data[ii]     = z2_data[ii];
        rhs2_data[ii + n] = z3_data[ii];
      }

      if (SUNLinSolSolve(rmem->LS_E2, rmem->E2, sol2, rhs2, SUN_RCONST(0.0)) != 0)
        return RADAU5_LSOLVE_FAIL;

      sol2_data = N_VGetArrayPointer(sol2);
      for (sunindextype ii = 0; ii < n; ii++)
      {
        z2_data[ii] = sol2_data[ii];
        z3_data[ii] = sol2_data[ii + n];
      }

      rmem->nsol++;
      newt++;
    } /* end eigenvalue/schur mode */

    /* ------------------------------------------------------------------
     * 6. Convergence norm  (RMS over all 3n components, scaled)
     * ----------------------------------------------------------------*/
    scal_data = N_VGetArrayPointer(scal);
    dyno = SUN_RCONST(0.0);
    for (sunindextype i = 0; i < n; i++)
    {
      sunrealtype denom = scal_data[i];
      dyno += (z1_data[i] / denom) * (z1_data[i] / denom)
            + (z2_data[i] / denom) * (z2_data[i] / denom)
            + (z3_data[i] / denom) * (z3_data[i] / denom);
    }
    dyno = SUNRsqrt(dyno / (SUN_RCONST(3.0) * (sunrealtype)n));

    /* ------------------------------------------------------------------
     * 7. Convergence monitoring (Fortran lines 949-972)
     * ----------------------------------------------------------------*/
    if (newt > 1 && newt < nit)
    {
      sunrealtype thq = dyno / dynold;
      sunrealtype thqold_prev = thqold;
      thqold = thq;

      if (newt == 2)
        theta = thq;
      else
        theta = SUNRsqrt(thq * thqold_prev);

      if (theta < SUN_RCONST(0.99))
      {
        faccon = theta / (SUN_RCONST(1.0) - theta);
        sunrealtype dyth = faccon * dyno
                         * SUNRpowerR(theta, (sunrealtype)(nit - 1 - newt))
                         / fnewt;
        if (dyth >= SUN_RCONST(1.0))
        {
          /* Predict required step-size reduction and signal retry.
           * The Fortran code reduces h here and goes directly to
           * label 20 or 10 (NOT label 78). We signal this with
           * RADAU5_NEWT_PREDICT so the caller skips the label78
           * h *= 0.5 reduction. */
          sunrealtype qnewt = SUNMAX(SUN_RCONST(1.0e-4),
                                     SUNMIN(SUN_RCONST(20.0), dyth));
          rmem->hhfac = SUN_RCONST(0.8)
                      * SUNRpowerR(qnewt,
                                   -SUN_RCONST(1.0)
                                   / (SUN_RCONST(4.0)
                                      + (sunrealtype)(nit - 1 - newt)));
          rmem->h    *= rmem->hhfac;
          rmem->reject = 1;
          rmem->last   = 0;
          *newt_out    = newt;
          return RADAU5_NEWT_PREDICT;
        }
      }
      else
      {
        /* theta >= 0.99: diverging */
        *newt_out = newt;
        return RADAU5_CONV_FAILURE;
      }
    }

    dynold = SUNMAX(dyno, uround);

    /* ------------------------------------------------------------------
     * 8. Update F1,F2,F3 (TI-space increments) and back-transform to
     *    Z1,Z2,Z3 (physical-space stage increments).
     *
     *    f_i += dz_i   (accumulate Newton correction in TI space)
     *
     *    Eigenvalue mode:
     *      Z1 = T11*f1 + T12*f2 + T13*f3
     *      Z2 = T21*f1 + T22*f2 + T23*f3
     *      Z3 = T31*f1 +      f2            (T32=1, T33=0)
     *
     *    Schur mode (T = US, full 3×3):
     *      Z1 = US[0][0]*f1 + US[0][1]*f2 + US[0][2]*f3
     *      Z2 = US[1][0]*f1 + US[1][1]*f2 + US[1][2]*f3
     *      Z3 = US[2][0]*f1 + US[2][1]*f2 + US[2][2]*f3
     * ----------------------------------------------------------------*/
    if (rmem->use_schur)
    {
      sunrealtype US00 = rmem->US[0][0], US01 = rmem->US[0][1], US02 = rmem->US[0][2];
      sunrealtype US10 = rmem->US[1][0], US11 = rmem->US[1][1], US12 = rmem->US[1][2];
      sunrealtype US20 = rmem->US[2][0], US21 = rmem->US[2][1], US22 = rmem->US[2][2];

      for (sunindextype ii = 0; ii < n; ii++)
      {
        sunrealtype f1i = f1_data[ii] + z1_data[ii];
        sunrealtype f2i = f2_data[ii] + z2_data[ii];
        sunrealtype f3i = f3_data[ii] + z3_data[ii];

        f1_data[ii] = f1i;
        f2_data[ii] = f2i;
        f3_data[ii] = f3i;

        z1_data[ii] = US00 * f1i + US01 * f2i + US02 * f3i;
        z2_data[ii] = US10 * f1i + US11 * f2i + US12 * f3i;
        z3_data[ii] = US20 * f1i + US21 * f2i + US22 * f3i;
      }
    }
    else
    {
      for (sunindextype ii = 0; ii < n; ii++)
      {
        sunrealtype f1i = f1_data[ii] + z1_data[ii];
        sunrealtype f2i = f2_data[ii] + z2_data[ii];
        sunrealtype f3i = f3_data[ii] + z3_data[ii];

        f1_data[ii] = f1i;
        f2_data[ii] = f2i;
        f3_data[ii] = f3i;

        z1_data[ii] = T11 * f1i + T12 * f2i + T13 * f3i;
        z2_data[ii] = T21 * f1i + T22 * f2i + T23 * f3i;
        z3_data[ii] = T31 * f1i +       f2i;  /* T32=1, T33=0 */
      }
    }

    /* ------------------------------------------------------------------
     * 9. Check convergence
     * ----------------------------------------------------------------*/
    if (faccon * dyno <= fnewt)
      break; /* converged */

  } /* end Newton loop */

  /* Store updated state back into rmem */
  rmem->faccon = faccon;
  rmem->theta  = theta;
  *newt_out    = newt;
  rmem->nnewt += newt;

  return RADAU5_SUCCESS;
}
