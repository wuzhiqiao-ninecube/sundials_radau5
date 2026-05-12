/* ---------------------------------------------------------------------------
 * RADAU5 — Generalized Newton iteration for variable-order (ns=3,5,7)
 *
 * Implements radau5_Newton() for arbitrary ns stages. The stage ordering is:
 *   z[0]/f[0] = real eigenvalue component (solved with E1)
 *   z[2k+1]/f[2k+1], z[2k+2]/f[2k+2] = complex pair k (solved with E2[k])
 *
 * This matches the Fortran T/TI matrix layout where row 0 corresponds to
 * the real eigenvalue and rows 1,2 / 3,4 / 5,6 to complex pairs.
 *
 * Reference: Hairer & Wanner, "Solving ODEs II", Fortran radau.f
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <nvector/nvector_serial.h>
#include "radau5_impl.h"

int radau5_Newton(Radau5Mem rmem, int* newt_out)
{
  sunindextype n  = rmem->n;
  sunrealtype  h  = rmem->h;
  sunrealtype  tn = rmem->tn;
  int ns     = rmem->ns;
  int npairs = rmem->npairs;

  /* Newton control */
  int         nit   = rmem->nit;
  sunrealtype fnewt = rmem->fnewt;
  sunrealtype uround = SUN_UNIT_ROUNDOFF;

  N_Vector ycur = rmem->ycur;
  N_Vector tmp1 = rmem->tmp1;
  N_Vector scal = rmem->scal;

  /* Convergence tracking */
  sunrealtype faccon = SUNMAX(rmem->faccon, uround);
  faccon = SUNRpowerR(faccon, SUN_RCONST(0.8));
  sunrealtype theta  = SUNRabs(rmem->thet);
  sunrealtype dynold = SUN_RCONST(1.0);
  sunrealtype thqold = SUN_RCONST(0.0);
  sunrealtype dyno   = SUN_RCONST(0.0);

  int newt = 0;

  /* Temporary buffer for mass-mult results */
  sunrealtype mf_stack[1024];
  sunrealtype *mf_buf = (n <= 1024) ? mf_stack
                       : (sunrealtype*)malloc((size_t)n * sizeof(sunrealtype));

  /* =========================================================================
   * Newton loop
   * =======================================================================*/
  while (1)
  {
    if (newt >= nit)
    {
      *newt_out = newt;
      if (n > 1024) free(mf_buf);
      return RADAU5_CONV_FAILURE;
    }

    /* ------------------------------------------------------------------
     * 1. Evaluate RHS at all ns collocation points.
     *    z[k] holds current stage increment Z_k.
     *    After: z[k] = f(t + c[k]*h, y + Z_k)
     * ----------------------------------------------------------------*/
    for (int k = 0; k < ns; k++)
    {
      N_VLinearSum(SUN_RCONST(1.0), ycur, SUN_RCONST(1.0), rmem->z[k], tmp1);
      sunrealtype tk = tn + rmem->c[k] * h;
      int rhsret = rmem->rhs(tk, tmp1, rmem->z[k], rmem->user_data);
      if (rhsret < 0) { if (n > 1024) free(mf_buf); return RADAU5_RHSFUNC_FAIL; }
      if (rhsret > 0) { if (n > 1024) free(mf_buf); return RADAU5_RHSFUNC_RECVR; }
    }
    rmem->nfcn += ns;

    /* ------------------------------------------------------------------
     * 2. Apply TI forward transform in-place.
     *    z[k] currently holds RHS values a[k].
     *    After: z[k] = sum_j TI[k][j] * a[j]
     * ----------------------------------------------------------------*/
    {
      sunrealtype a_vals[RADAU5_NS_MAX];
      sunrealtype *zd[RADAU5_NS_MAX];
      for (int k = 0; k < ns; k++)
        zd[k] = N_VGetArrayPointer(rmem->z[k]);

      for (sunindextype i = 0; i < n; i++)
      {
        for (int k = 0; k < ns; k++)
          a_vals[k] = zd[k][i];
        for (int k = 0; k < ns; k++)
        {
          sunrealtype sum = SUN_RCONST(0.0);
          for (int j = 0; j < ns; j++)
            sum += rmem->TI_mat[k * ns + j] * a_vals[j];
          zd[k][i] = sum;
        }
      }
    }

    /* ------------------------------------------------------------------
     * 3. Form linear-system RHS (eigenvalue mode only for now).
     *
     * Stage ordering:
     *   z[0] = real eigenvalue: z[0] -= fac1 * M*f[0]
     *   z[2k+1],z[2k+2] = complex pair k:
     *     z[2k+1] += alphn[k]*S[2k+1] - betan[k]*S[2k+2]
     *     z[2k+2] += alphn[k]*S[2k+2] + betan[k]*S[2k+1]
     *   where S[j] = -(M * f[j])  (or -f[j] for identity mass)
     * ----------------------------------------------------------------*/
    sunrealtype fac1 = rmem->u1 / h;

    sunrealtype *fd[RADAU5_NS_MAX], *zd[RADAU5_NS_MAX];
    for (int k = 0; k < ns; k++)
    {
      fd[k] = N_VGetArrayPointer(rmem->f[k]);
      zd[k] = N_VGetArrayPointer(rmem->z[k]);
    }

    if (rmem->use_schur)
    {
      /* --- Schur mode --- */
      /* Form RHS: z[r] -= sum_j (TS[r][j]/h) * M*f[j] for upper-triangular TS */
      sunrealtype *mf_ptrs[RADAU5_NS_MAX];
      sunrealtype *mf_alloc = NULL;

      if (rmem->M != NULL)
      {
        mf_alloc = (sunrealtype*)malloc((size_t)(ns * n) * sizeof(sunrealtype));
        for (int j = 0; j < ns; j++)
        {
          mf_ptrs[j] = mf_alloc + j * n;
          radau5_MassMult(rmem, rmem->f[j], tmp1);
          sunrealtype *td = N_VGetArrayPointer(tmp1);
          for (sunindextype ii = 0; ii < n; ii++)
            mf_ptrs[j][ii] = td[ii];
        }
      }
      else
      {
        for (int j = 0; j < ns; j++)
          mf_ptrs[j] = fd[j];
      }

      /* Apply TS row scalings (TS is quasi-triangular, NOT upper triangular) */
      for (int r = 0; r < ns; r++)
      {
        for (sunindextype ii = 0; ii < n; ii++)
        {
          sunrealtype sum = SUN_RCONST(0.0);
          for (int j = 0; j < ns; j++)
          {
            sunrealtype tsrj = rmem->TS_mat[r * ns + j];
            if (tsrj != SUN_RCONST(0.0))
              sum += (tsrj / h) * mf_ptrs[j][ii];
          }
          zd[r][ii] -= sum;
        }
      }

      /* Block back-substitution:
       * 1. Solve E1 for z[ns-1] (1×1 block at bottom-right) */
      if (SUNLinSolSolve(rmem->LS_E1, rmem->E1, rmem->z[ns-1], rmem->z[ns-1],
                         SUN_RCONST(0.0)) != 0)
      {
        if (mf_alloc) free(mf_alloc);
        if (n > 1024) free(mf_buf);
        return RADAU5_LSOLVE_FAIL;
      }

      /* 2. Back-substitute z[ns-1] into earlier rows */
      {
        sunrealtype *df_last;
        if (rmem->M != NULL)
        {
          radau5_MassMult(rmem, rmem->z[ns-1], tmp1);
          df_last = N_VGetArrayPointer(tmp1);
        }
        else
          df_last = N_VGetArrayPointer(rmem->z[ns-1]);

        for (int r = 0; r < ns - 1; r++)
        {
          sunrealtype coeff = rmem->TS_mat[r * ns + (ns - 1)] / h;
          if (coeff != SUN_RCONST(0.0))
            for (sunindextype ii = 0; ii < n; ii++)
              zd[r][ii] -= coeff * df_last[ii];
        }
      }

      /* 3. Solve each 2×2 block from bottom pair to top */
      for (int pk = npairs - 1; pk >= 0; pk--)
      {
        int r0 = 2 * pk, r1 = r0 + 1;

        sunrealtype *rhs2d = N_VGetArrayPointer(rmem->rhs2[pk]);
        for (sunindextype ii = 0; ii < n; ii++)
        {
          rhs2d[ii]     = zd[r0][ii];
          rhs2d[ii + n] = zd[r1][ii];
        }

        if (SUNLinSolSolve(rmem->LS_E2[pk], rmem->E2[pk], rmem->sol2[pk],
                           rmem->rhs2[pk], SUN_RCONST(0.0)) != 0)
        {
          if (mf_alloc) free(mf_alloc);
          if (n > 1024) free(mf_buf);
          return RADAU5_LSOLVE_FAIL;
        }

        sunrealtype *sol2d = N_VGetArrayPointer(rmem->sol2[pk]);
        for (sunindextype ii = 0; ii < n; ii++)
        {
          zd[r0][ii] = sol2d[ii];
          zd[r1][ii] = sol2d[ii + n];
        }

        /* Back-substitute into earlier rows */
        if (pk > 0)
        {
          for (int sc = r0; sc <= r1; sc++)
          {
            sunrealtype *df_col;
            if (rmem->M != NULL)
            {
              radau5_MassMult(rmem, rmem->z[sc], tmp1);
              df_col = N_VGetArrayPointer(tmp1);
            }
            else
              df_col = zd[sc];

            for (int r = 0; r < r0; r++)
            {
              sunrealtype coeff = rmem->TS_mat[r * ns + sc] / h;
              if (coeff != SUN_RCONST(0.0))
                for (sunindextype ii = 0; ii < n; ii++)
                  zd[r][ii] -= coeff * df_col[ii];
            }
          }
        }
      }

      if (mf_alloc) free(mf_alloc);
      rmem->nsol++;
      newt++;
    }
    else
    {
      /* --- Eigenvalue mode --- */
      if (rmem->M != NULL)
      {
        /* General mass: S[j] = -(M * f[j]) */
        /* Real eigenvalue: z[0] -= fac1 * M*f[0] */
        radau5_MassMult(rmem, rmem->f[0], tmp1);
        sunrealtype *mfd = N_VGetArrayPointer(tmp1);
        for (sunindextype ii = 0; ii < n; ii++)
          zd[0][ii] -= fac1 * mfd[ii];

        /* Complex pairs */
        for (int pk = 0; pk < npairs; pk++)
        {
          int r0 = 2 * pk + 1, r1 = r0 + 1;
          sunrealtype alphn_k = rmem->alph[pk] / h;
          sunrealtype betan_k = rmem->beta_eig[pk] / h;

          /* S[r0] = -(M * f[r0]) */
          radau5_MassMult(rmem, rmem->f[r0], tmp1);
          mfd = N_VGetArrayPointer(tmp1);
          for (sunindextype ii = 0; ii < n; ii++)
            mf_buf[ii] = -mfd[ii];

          /* S[r1] = -(M * f[r1]) */
          radau5_MassMult(rmem, rmem->f[r1], tmp1);
          mfd = N_VGetArrayPointer(tmp1);

          for (sunindextype ii = 0; ii < n; ii++)
          {
            sunrealtype s_r0 = mf_buf[ii];
            sunrealtype s_r1 = -mfd[ii];
            zd[r0][ii] += alphn_k * s_r0 - betan_k * s_r1;
            zd[r1][ii] += alphn_k * s_r1 + betan_k * s_r0;
          }
        }
      }
      else
      {
        /* Identity mass: S[j] = -f[j] */
        /* Real eigenvalue: z[0] -= fac1 * f[0] */
        for (sunindextype ii = 0; ii < n; ii++)
          zd[0][ii] -= fac1 * fd[0][ii];

        /* Complex pairs */
        for (int pk = 0; pk < npairs; pk++)
        {
          int r0 = 2 * pk + 1, r1 = r0 + 1;
          sunrealtype alphn_k = rmem->alph[pk] / h;
          sunrealtype betan_k = rmem->beta_eig[pk] / h;

          for (sunindextype ii = 0; ii < n; ii++)
          {
            sunrealtype s_r0 = -fd[r0][ii];
            sunrealtype s_r1 = -fd[r1][ii];
            zd[r0][ii] += alphn_k * s_r0 - betan_k * s_r1;
            zd[r1][ii] += alphn_k * s_r1 + betan_k * s_r0;
          }
        }
      }

      /* Solve E1 for z[0] (real eigenvalue) */
      if (SUNLinSolSolve(rmem->LS_E1, rmem->E1, rmem->z[0], rmem->z[0],
                         SUN_RCONST(0.0)) != 0)
      {
        if (n > 1024) free(mf_buf);
        return RADAU5_LSOLVE_FAIL;
      }

      /* Solve E2[pk] for each complex pair */
      for (int pk = 0; pk < npairs; pk++)
      {
        int r0 = 2 * pk + 1, r1 = r0 + 1;
        sunrealtype *rhs2d = N_VGetArrayPointer(rmem->rhs2[pk]);
        for (sunindextype ii = 0; ii < n; ii++)
        {
          rhs2d[ii]     = zd[r0][ii];
          rhs2d[ii + n] = zd[r1][ii];
        }

        if (SUNLinSolSolve(rmem->LS_E2[pk], rmem->E2[pk], rmem->sol2[pk],
                           rmem->rhs2[pk], SUN_RCONST(0.0)) != 0)
        {
          if (n > 1024) free(mf_buf);
          return RADAU5_LSOLVE_FAIL;
        }

        sunrealtype *sol2d = N_VGetArrayPointer(rmem->sol2[pk]);
        for (sunindextype ii = 0; ii < n; ii++)
        {
          zd[r0][ii] = sol2d[ii];
          zd[r1][ii] = sol2d[ii + n];
        }
      }

      rmem->nsol++;
      newt++;
    } /* end eigenvalue/schur mode */

    /* ------------------------------------------------------------------
     * 4. Convergence norm (RMS over all ns*n components, scaled)
     * ----------------------------------------------------------------*/
    {
      sunrealtype *scal_data = N_VGetArrayPointer(scal);
      dyno = SUN_RCONST(0.0);
      for (sunindextype i = 0; i < n; i++)
      {
        sunrealtype denom = scal_data[i];
        for (int k = 0; k < ns; k++)
        {
          sunrealtype val = zd[k][i] / denom;
          dyno += val * val;
        }
      }
      dyno = SUNRsqrt(dyno / ((sunrealtype)ns * (sunrealtype)n));
    }

    /* ------------------------------------------------------------------
     * 5. Convergence monitoring
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
          if (n > 1024) free(mf_buf);
          return RADAU5_NEWT_PREDICT;
        }
      }
      else
      {
        *newt_out = newt;
        if (n > 1024) free(mf_buf);
        return RADAU5_CONV_FAILURE;
      }
    }

    dynold = SUNMAX(dyno, uround);

    /* ------------------------------------------------------------------
     * 6. Update f[k] and back-transform to z[k].
     *    f[k] += dz[k]  (accumulate Newton correction in TI space)
     *    Z[r] = sum_j T[r][j] * f[j]  (back-transform to physical space)
     * ----------------------------------------------------------------*/
    {
      sunrealtype f_vals[RADAU5_NS_MAX];
      for (sunindextype ii = 0; ii < n; ii++)
      {
        for (int k = 0; k < ns; k++)
          f_vals[k] = fd[k][ii] + zd[k][ii];

        /* Store updated f values */
        for (int k = 0; k < ns; k++)
          fd[k][ii] = f_vals[k];

        /* Back-transform: Z[r] = sum_j T[r][j] * f[j] */
        for (int r = 0; r < ns; r++)
        {
          sunrealtype sum = SUN_RCONST(0.0);
          for (int j = 0; j < ns; j++)
            sum += rmem->T_mat[r * ns + j] * f_vals[j];
          zd[r][ii] = sum;
        }
      }
    }

    /* ------------------------------------------------------------------
     * 7. Check convergence
     * ----------------------------------------------------------------*/
    if (faccon * dyno <= fnewt)
      break; /* converged */

  } /* end Newton loop */

  /* Store updated state */
  rmem->faccon = faccon;
  rmem->theta  = theta;
  *newt_out    = newt;
  rmem->nnewt += newt;

  if (n > 1024) free(mf_buf);
  return RADAU5_SUCCESS;
}
