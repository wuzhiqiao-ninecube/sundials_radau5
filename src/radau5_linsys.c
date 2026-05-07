/* ---------------------------------------------------------------------------
 * radau5_linsys.c — Linear system infrastructure for RADAU5
 *
 * Implements:
 *   radau5_InitConstants  — Radau IIA method constants (from Fortran radau5.f)
 *   radau5_DQJacDense     — Dense finite-difference Jacobian (CVODE-style)
 *   radau5_BuildE1        — Assemble E1 = fac1*M - J  (n×n real system)
 *   radau5_BuildE2        — Assemble E2 realified 2n×2n complex system
 *   radau5_DecompE1       — Factor E1 via SUNLinSolSetup
 *   radau5_DecompE2       — Factor E2 via SUNLinSolSetup
 *   radau5_ComputeScal    — Error weight vector scal[i] = atol[i] + rtol[i]*|y[i]|
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include "radau5_impl.h"

/* ---------------------------------------------------------------------------
 * radau5_InitConstants
 *
 * Initialize all Radau IIA 3-stage order-5 method constants.
 * Transcribed exactly from Fortran radau5.f lines 706-737.
 * ---------------------------------------------------------------------------*/
int radau5_InitConstants(Radau5Mem rmem)
{
  sunrealtype SQ6, c1, c2, c1m1, c2m1, c1mc2;
  sunrealtype u1, alph, beta, cno;
  sunrealtype uround, rtol;

  SQ6   = sqrt(SUN_RCONST(6.0));
  c1    = (SUN_RCONST(4.0) - SQ6) / SUN_RCONST(10.0);
  c2    = (SUN_RCONST(4.0) + SQ6) / SUN_RCONST(10.0);
  c1m1  = c1 - SUN_RCONST(1.0);
  c2m1  = c2 - SUN_RCONST(1.0);
  c1mc2 = c1 - c2;

  rmem->c1   = c1;
  rmem->c2   = c2;
  rmem->c1m1 = c1m1;
  rmem->c2m1 = c2m1;
  rmem->c1mc2 = c1mc2;

  rmem->dd1 = -(SUN_RCONST(13.0) + SUN_RCONST(7.0) * SQ6) / SUN_RCONST(3.0);
  rmem->dd2 = (-SUN_RCONST(13.0) + SUN_RCONST(7.0) * SQ6) / SUN_RCONST(3.0);
  rmem->dd3 = -SUN_RCONST(1.0) / SUN_RCONST(3.0);

  u1   = (SUN_RCONST(6.0) + pow(SUN_RCONST(81.0), SUN_RCONST(1.0)/SUN_RCONST(3.0))
                           - pow(SUN_RCONST(9.0),  SUN_RCONST(1.0)/SUN_RCONST(3.0)))
         / SUN_RCONST(30.0);
  alph = (SUN_RCONST(12.0) - pow(SUN_RCONST(81.0), SUN_RCONST(1.0)/SUN_RCONST(3.0))
                            + pow(SUN_RCONST(9.0),  SUN_RCONST(1.0)/SUN_RCONST(3.0)))
         / SUN_RCONST(60.0);
  beta = (pow(SUN_RCONST(81.0), SUN_RCONST(1.0)/SUN_RCONST(3.0))
        + pow(SUN_RCONST(9.0),  SUN_RCONST(1.0)/SUN_RCONST(3.0)))
         * sqrt(SUN_RCONST(3.0)) / SUN_RCONST(60.0);

  cno  = alph * alph + beta * beta;

  /* Invert: u1 → 1/u1, alph → alph/cno, beta → beta/cno */
  u1   = SUN_RCONST(1.0) / u1;
  alph = alph / cno;
beta = beta / cno;

  rmem->u1   = u1;
  rmem->alph = alph;
  rmem->beta = beta;

  /* Eigenvector matrix T (hardcoded from Fortran) */
  rmem->T11 =  SUN_RCONST(9.1232394870892942792e-02);
  rmem->T12 = -SUN_RCONST(0.14125529502095420843);
  rmem->T13 = -SUN_RCONST(3.0029194105147424492e-02);
  rmem->T21 =  SUN_RCONST(0.24171793270710701896);
  rmem->T22 =  SUN_RCONST(0.20412935229379993199);
  rmem->T23 =  SUN_RCONST(0.38294211275726193779);
  rmem->T31 =  SUN_RCONST(0.96604818261509293619);
  /* T32 = 1, T33 = 0 are implicit in the Fortran formulation */

  /* Inverse eigenvector matrix TI */
  rmem->TI11 =  SUN_RCONST(4.3255798900631553510);
  rmem->TI12 =  SUN_RCONST(0.33919925181580986954);
  rmem->TI13 =  SUN_RCONST(0.54177053993587487119);
  rmem->TI21 = -SUN_RCONST(4.1787185915519047273);
  rmem->TI22 = -SUN_RCONST(0.32768282076106238708);
  rmem->TI23 =  SUN_RCONST(0.47662355450055045196);
  rmem->TI31 = -SUN_RCONST(0.50287263494578687595);
  rmem->TI32 =  SUN_RCONST(2.5719269498556054292);
  rmem->TI33 = -SUN_RCONST(0.59603920482822492497);

  /* --- Tolerance transformation is deferred to Radau5Solve (first call)
   * because user may call Radau5SStolerances after Radau5Init. --- */

  /* fnewt will be computed after tolerance transformation in Radau5Solve */
  rmem->fnewt = SUN_RCONST(0.0);  /* sentinel: recompute later */

  /* --- Schur decomposition constants ---
   * If use_schur==1, populate US and TS from precomputed values,
   * and set T/TI from US/US^T for the forward/back transforms.
   * Also set u1 = TS[2][2] (the real eigenvalue for the 1×1 block). */
  if (rmem->use_schur)
  {
    /* US: orthogonal Schur vectors (columns of US) */
    rmem->US[0][0] =  SUN_RCONST( 0.138665108751908);
    rmem->US[0][1] =  SUN_RCONST( 0.046278149309488);
    rmem->US[0][2] =  SUN_RCONST( 0.989257459163847);
    rmem->US[1][0] =  SUN_RCONST(-0.229641242351741);
    rmem->US[1][1] =  SUN_RCONST(-0.970178886551833);
    rmem->US[1][2] =  SUN_RCONST( 0.077574660168092);
    rmem->US[2][0] =  SUN_RCONST(-0.963346711950568);
    rmem->US[2][1] =  SUN_RCONST( 0.237931210616713);
    rmem->US[2][2] =  SUN_RCONST( 0.123902589111344);

    /* TS: upper quasi-triangular Schur form */
    rmem->TS[0][0] =  SUN_RCONST( 2.68108287362775);
    rmem->TS[0][1] =  SUN_RCONST(-8.42387579808538);
    rmem->TS[0][2] =  SUN_RCONST(-4.08857781305964);
    rmem->TS[1][0] =  SUN_RCONST( 1.10461319985220);
    rmem->TS[1][1] =  SUN_RCONST( 2.68108287362775);
    rmem->TS[1][2] =  SUN_RCONST( 4.70015200634394);
    rmem->TS[2][0] =  SUN_RCONST( 0.0);
    rmem->TS[2][1] =  SUN_RCONST( 0.0);
    rmem->TS[2][2] =  SUN_RCONST( 3.63783425274450);

    /* Override eigenvalue scalings:
     * u1 = TS[2][2] (real eigenvalue, used as fac1 = u1/h for E1)
     * alph/beta are not used directly — E2 is built from TS[0:2][0:2] block */
    rmem->u1 = rmem->TS[2][2];

    /* Set T = US and TI = US^T for the Newton transforms.
     * Note: with Schur, T32 and T33 are NOT 1 and 0. We use the
     * full 3×3 T matrix. The Newton back-transform must handle this. */
    rmem->T11 = rmem->US[0][0];
    rmem->T12 = rmem->US[0][1];
    rmem->T13 = rmem->US[0][2];
    rmem->T21 = rmem->US[1][0];
    rmem->T22 = rmem->US[1][1];
    rmem->T23 = rmem->US[1][2];
    rmem->T31 = rmem->US[2][0];
    /* T32, T33 not stored as separate fields in the eigenvalue path
     * (T32=1, T33=0 implicit). For Schur, we handle them in the
     * Newton code using US[][] directly. */

    rmem->TI11 = rmem->US[0][0];  /* US^T row 0 = US col 0 */
    rmem->TI12 = rmem->US[1][0];
    rmem->TI13 = rmem->US[2][0];
    rmem->TI21 = rmem->US[0][1];  /* US^T row 1 = US col 1 */
    rmem->TI22 = rmem->US[1][1];
    rmem->TI23 = rmem->US[2][1];
    rmem->TI31 = rmem->US[0][2];  /* US^T row 2 = US col 2 */
    rmem->TI32 = rmem->US[1][2];
    rmem->TI33 = rmem->US[2][2];
  }

  return RADAU5_SUCCESS;
}

/* ---------------------------------------------------------------------------
 * radau5_SparseLookup
 *
 * Binary search for value at (row, col) in a CSC sparse matrix.
 * Returns 0.0 if the entry is not in the sparsity pattern.
 * ---------------------------------------------------------------------------*/
static sunrealtype radau5_SparseLookup(SUNMatrix A, sunindextype row,
                                       sunindextype col)
{
  sunindextype *Ap = SM_INDEXPTRS_S(A);
  sunindextype *Ai = SM_INDEXVALS_S(A);
  sunrealtype  *Ad = SM_DATA_S(A);
  sunindextype lo = Ap[col], hi = Ap[col + 1];
  while (lo < hi)
  {
    sunindextype mid = lo + (hi - lo) / 2;
    if (Ai[mid] == row) return Ad[mid];
    else if (Ai[mid] < row) lo = mid + 1;
    else hi = mid;
  }
  return SUN_RCONST(0.0);
}

/* ---------------------------------------------------------------------------
 * radau5_SparseUnion
 *
 * Create a new n×n CSC sparse matrix whose sparsity pattern is the union
 * of A's and B's patterns. Values are zeroed. Row indices are sorted.
 * ---------------------------------------------------------------------------*/
SUNMatrix radau5_SparseUnion(SUNMatrix A, SUNMatrix B, SUNContext sunctx)
{
  sunindextype n = SM_COLUMNS_S(A);
  sunindextype *Ap = SM_INDEXPTRS_S(A);
  sunindextype *Ai = SM_INDEXVALS_S(A);
  sunindextype *Bp = SM_INDEXPTRS_S(B);
  sunindextype *Bi = SM_INDEXVALS_S(B);

  /* First pass: count union NNZ per column */
  sunindextype total_nnz = 0;
  sunindextype j;
  for (j = 0; j < n; j++)
  {
    sunindextype ka = Ap[j], ka_end = Ap[j + 1];
    sunindextype kb = Bp[j], kb_end = Bp[j + 1];
    while (ka < ka_end && kb < kb_end)
    {
      if (Ai[ka] < Bi[kb])      { total_nnz++; ka++; }
      else if (Ai[ka] > Bi[kb]) { total_nnz++; kb++; }
      else                       { total_nnz++; ka++; kb++; }
    }
    total_nnz += (ka_end - ka) + (kb_end - kb);
  }

  SUNMatrix C = SUNSparseMatrix(n, n, total_nnz, CSC_MAT, sunctx);
  if (!C) return NULL;

  sunindextype *Cp = SM_INDEXPTRS_S(C);
  sunindextype *Ci = SM_INDEXVALS_S(C);

  /* Second pass: fill pattern */
  sunindextype pos = 0;
  for (j = 0; j < n; j++)
  {
    Cp[j] = pos;
    sunindextype ka = Ap[j], ka_end = Ap[j + 1];
    sunindextype kb = Bp[j], kb_end = Bp[j + 1];
    while (ka < ka_end && kb < kb_end)
    {
      if (Ai[ka] < Bi[kb])      { Ci[pos++] = Ai[ka++]; }
      else if (Ai[ka] > Bi[kb]) { Ci[pos++] = Bi[kb++]; }
      else                       { Ci[pos++] = Ai[ka++]; kb++; }
    }
    while (ka < ka_end) Ci[pos++] = Ai[ka++];
    while (kb < kb_end) Ci[pos++] = Bi[kb++];
  }
  Cp[n] = pos;

  return C;
}

/* ---------------------------------------------------------------------------
 * radau5_DQJacDense
 *
 * Column-by-column forward-difference Jacobian for dense matrices.
 * Mirrors CVODE's cvLsDenseDQJac approach.
 * ---------------------------------------------------------------------------*/
int radau5_DQJacDense(Radau5Mem rmem, sunrealtype t, N_Vector y, N_Vector fy)
{
  sunindextype j, i, n;
  sunrealtype  srur, inc, inc_inv, ysave;
  sunrealtype *y_data, *fy_data, *tmp1_data;
  SUNMatrix    J;

  n         = rmem->n;
  J         = rmem->J;
  srur      = SUNRsqrt(SUN_UNIT_ROUNDOFF);

  fy_data   = N_VGetArrayPointer(fy);

  for (j = 0; j < n; j++)
  {
    y_data    = N_VGetArrayPointer(y);
    ysave     = y_data[j];

    /* Fortran radau5.f: DELT=DSQRT(UROUND*MAX(1.D-5,ABS(YSAFE))) */
    inc = SUNRsqrt(SUN_UNIT_ROUNDOFF * SUNMAX(SUN_RCONST(1.0e-5), fabs(ysave)));

    y_data[j] = ysave + inc;
    {
      int rhsret = rmem->rhs(t, y, rmem->tmp1, rmem->user_data);
      if (rhsret < 0) { y_data[j] = ysave; return RADAU5_RHSFUNC_FAIL; }
      if (rhsret > 0) { y_data[j] = ysave; return RADAU5_RHSFUNC_RECVR; }
    }
    y_data[j] = ysave;

    inc_inv   = SUN_RCONST(1.0) / inc;
    tmp1_data = N_VGetArrayPointer(rmem->tmp1);

    for (i = 0; i < n; i++)
      SM_ELEMENT_D(J, i, j) = (tmp1_data[i] - fy_data[i]) * inc_inv;

    rmem->nfcn++;
  }

  return RADAU5_SUCCESS;
}

/* ---------------------------------------------------------------------------
 * radau5_DQJacBand
 *
 * Band finite-difference Jacobian using Curtis-Powell-Reid (CPR) grouped
 * column perturbation.  Mirrors CVODE's cvLsBandDQJac.
 * Total RHS evaluations: min(mu + ml + 1, n).
 * ---------------------------------------------------------------------------*/
int radau5_DQJacBand(Radau5Mem rmem, sunrealtype t, N_Vector y, N_Vector fy)
{
  sunindextype i, j, n, group, i1, i2;
  sunindextype mu, ml, width, ngroups;
  sunrealtype  srur, fnorm, minInc, inc, inc_inv;
  sunrealtype *y_data, *fy_data, *ytemp_data, *ftemp_data, *scal_data;
  sunrealtype *col_j;
  sunrealtype  inc_arr[1024]; /* stack buffer for per-column inc values */
  sunrealtype *incs;
  SUNMatrix    J;

  n    = rmem->n;
  J    = rmem->J;
  mu   = rmem->mu;
  ml   = rmem->ml;
  srur = SUNRsqrt(SUN_UNIT_ROUNDOFF);

  fnorm  = N_VWrmsNorm(fy, rmem->scal);
  minInc = (fnorm != SUN_RCONST(0.0))
           ? (SUN_RCONST(1000.0) * fabs(rmem->h) * SUN_UNIT_ROUNDOFF
              * (sunrealtype)n * fnorm)
           : SUN_RCONST(1.0);

  y_data    = N_VGetArrayPointer(y);
  fy_data   = N_VGetArrayPointer(fy);
  scal_data = N_VGetArrayPointer(rmem->scal);

  /* Use tmp1 as ytemp (copy of y), tmp2 as ftemp */
  N_VScale(SUN_RCONST(1.0), y, rmem->tmp1);
  ytemp_data = N_VGetArrayPointer(rmem->tmp1);
  ftemp_data = N_VGetArrayPointer(rmem->tmp2);

  /* Allocate inc array if n > stack buffer size */
  incs = (n <= 1024) ? inc_arr : (sunrealtype*)malloc((size_t)n * sizeof(sunrealtype));

  width   = ml + mu + 1;
  ngroups = SUNMIN(width, n);

  SUNMatZero(J);

  for (group = 0; group < ngroups; group++)
  {
    /* Perturb all columns in this group simultaneously */
    for (j = group; j < n; j += width)
    {
      inc = SUNMAX(srur * fabs(y_data[j]), minInc / scal_data[j]);
      inc = (y_data[j] + inc) - y_data[j];
      if (inc == SUN_RCONST(0.0)) inc = srur;
      incs[j]        = inc;
      ytemp_data[j] += inc;
    }

    /* One RHS evaluation for the whole group */
    {
      int rhsret = rmem->rhs(t, rmem->tmp1, rmem->tmp2, rmem->user_data);
      if (rhsret < 0) { if (n > 1024) free(incs); return RADAU5_RHSFUNC_FAIL; }
      if (rhsret > 0) {
        for (j = group; j < n; j += width)
          ytemp_data[j] = y_data[j];
        if (n > 1024) free(incs);
        return RADAU5_RHSFUNC_RECVR;
      }
    }
    rmem->nfcn++;

    /* Restore ytemp and fill Jacobian columns */
    for (j = group; j < n; j += width)
    {
      ytemp_data[j] = y_data[j];
      inc_inv = SUN_RCONST(1.0) / incs[j];
      col_j   = SUNBandMatrix_Column(J, j);

      i1 = SUNMAX(0, j - mu);
      i2 = SUNMIN(j + ml, n - 1);
      for (i = i1; i <= i2; i++)
        SM_COLUMN_ELEMENT_B(col_j, i, j) =
          (ftemp_data[i] - fy_data[i]) * inc_inv;
    }
  }

  if (n > 1024) free(incs);

  return RADAU5_SUCCESS;
}

/* ---------------------------------------------------------------------------
 * radau5_DQJacSparse
 *
 * Sparse finite-difference Jacobian using column grouping (CPR technique).
 * Columns in the same group are structurally independent (no shared nonzero
 * rows), so they can be perturbed simultaneously in a single RHS evaluation.
 * Total RHS evaluations: ngroups (typically << n for sparse patterns).
 *
 * Requires that Radau5SetSparsityPattern has been called to compute the
 * column grouping (col_group, group_offsets, group_cols).
 * ---------------------------------------------------------------------------*/
int radau5_DQJacSparse(Radau5Mem rmem, sunrealtype t, N_Vector y, N_Vector fy)
{
  sunindextype j, k, p, n;
  sunrealtype  srur, fnorm, minInc, inc, inc_inv;
  sunrealtype *y_data, *fy_data, *ytemp_data, *ftemp_data, *scal_data;
  sunrealtype  inc_arr[1024];
  sunrealtype *incs;

  n    = rmem->n;
  srur = SUNRsqrt(SUN_UNIT_ROUNDOFF);

  fnorm  = N_VWrmsNorm(fy, rmem->scal);
  minInc = (fnorm != SUN_RCONST(0.0))
           ? (SUN_RCONST(1000.0) * fabs(rmem->h) * SUN_UNIT_ROUNDOFF
              * (sunrealtype)n * fnorm)
           : SUN_RCONST(1.0);

  y_data    = N_VGetArrayPointer(y);
  fy_data   = N_VGetArrayPointer(fy);
  scal_data = N_VGetArrayPointer(rmem->scal);

  /* Use tmp1 as ytemp (copy of y), tmp2 as ftemp */
  N_VScale(SUN_RCONST(1.0), y, rmem->tmp1);
  ytemp_data = N_VGetArrayPointer(rmem->tmp1);
  ftemp_data = N_VGetArrayPointer(rmem->tmp2);

  /* Per-column increment storage */
  incs = (n <= 1024) ? inc_arr : (sunrealtype*)malloc((size_t)n * sizeof(sunrealtype));

  /* Zero J data array (preserve CSC structure) */
  sunrealtype* J_data = SM_DATA_S(rmem->J);
  sunindextype* Jp = SM_INDEXPTRS_S(rmem->J);
  sunindextype nnz_J = Jp[n];
  for (p = 0; p < nnz_J; p++) J_data[p] = SUN_RCONST(0.0);

  /* Sparsity pattern and group lookup from SetSparsityPattern */
  const sunindextype* sp_colptrs    = rmem->sp_colptrs;
  const sunindextype* sp_rowinds    = rmem->sp_rowinds;
  const sunindextype* group_offsets = rmem->group_offsets;
  const sunindextype* group_cols    = rmem->group_cols;

  for (sunindextype group = 0; group < rmem->ngroups; group++)
  {
    /* Perturb all columns in this group simultaneously */
    for (k = group_offsets[group]; k < group_offsets[group + 1]; k++)
    {
      j = group_cols[k];
      inc = SUNMAX(srur * fabs(y_data[j]), minInc / scal_data[j]);
      inc = (y_data[j] + inc) - y_data[j];
      if (inc == SUN_RCONST(0.0)) inc = srur;
      incs[j]        = inc;
      ytemp_data[j] += inc;
    }

    /* One RHS evaluation for the whole group */
    {
      int rhsret = rmem->rhs(t, rmem->tmp1, rmem->tmp2, rmem->user_data);
      if (rhsret < 0) { if (n > 1024) free(incs); return RADAU5_RHSFUNC_FAIL; }
      if (rhsret > 0) {
        for (k = group_offsets[group]; k < group_offsets[group + 1]; k++)
          ytemp_data[group_cols[k]] = y_data[group_cols[k]];
        if (n > 1024) free(incs);
        return RADAU5_RHSFUNC_RECVR;
      }
    }
    rmem->nfcn++;

    /* Extract Jacobian columns and restore ytemp */
    for (k = group_offsets[group]; k < group_offsets[group + 1]; k++)
    {
      j = group_cols[k];
      ytemp_data[j] = y_data[j];
      inc_inv = SUN_RCONST(1.0) / incs[j];

      for (p = sp_colptrs[j]; p < sp_colptrs[j + 1]; p++)
      {
        sunindextype row = sp_rowinds[p];
        J_data[p] = (ftemp_data[row] - fy_data[row]) * inc_inv;
      }
    }
  }

  if (n > 1024) free(incs);

  return RADAU5_SUCCESS;
}

/* ---------------------------------------------------------------------------
 * radau5_BuildE1
 *
 * Build E1 = fac1*M - J  (or fac1*I - J when M is NULL).
 * Dispatches on matrix type: dense, band, or sparse.
 * ---------------------------------------------------------------------------*/
int radau5_BuildE1(Radau5Mem rmem, sunrealtype fac1)
{
  sunindextype i, j, k, n;
  SUNMatrix    J, M, E1;

  n  = rmem->n;
  J  = rmem->J;
  M  = rmem->M;
  E1 = rmem->E1;

  if (rmem->mat_id == SUNMATRIX_SPARSE)
  {
    /* --- Sparse (CSC) case ---
     * When M is sparse, E1 has the union(J,M) sparsity pattern.
     * When M is NULL/dense/band, E1 has J's pattern (cloned).
     * Fill: E1[i,j] = fac1*M[i,j] - J[i,j]  (or fac1*delta(i,j) - J[i,j]) */
    sunindextype *E1p = SM_INDEXPTRS_S(E1);
    sunindextype *E1i = SM_INDEXVALS_S(E1);
    sunrealtype  *E1d = SM_DATA_S(E1);

    if (M != NULL && SUNMatGetID(M) == SUNMATRIX_SPARSE)
    {
      /* E1 pattern = union(J, M). Look up J and M values by binary search. */
      for (j = 0; j < n; j++)
      {
        for (k = E1p[j]; k < E1p[j+1]; k++)
        {
          i = E1i[k];
          sunrealtype jij = radau5_SparseLookup(J, i, j);
          sunrealtype mij = radau5_SparseLookup(M, i, j);
          E1d[k] = fac1 * mij - jij;
        }
      }
    }
    else
    {
      /* E1 pattern = J's pattern. Iterate over J directly. */
      sunindextype *Jp = SM_INDEXPTRS_S(J);
      sunindextype *Ji = SM_INDEXVALS_S(J);
      sunrealtype  *Jd = SM_DATA_S(J);

      for (j = 0; j < n; j++)
      {
        for (k = Jp[j]; k < Jp[j+1]; k++)
        {
          i = Ji[k];
          E1d[k] = -Jd[k];
          if (M == NULL)
          {
            if (i == j) E1d[k] += fac1;
          }
          else if (SUNMatGetID(M) == SUNMATRIX_DENSE)
          {
            E1d[k] += fac1 * SM_ELEMENT_D(M, i, j);
          }
          else if (SUNMatGetID(M) == SUNMATRIX_BAND)
          {
            sunindextype muM = SM_UBAND_B(M);
            sunindextype mlM = SM_LBAND_B(M);
            if (i >= j - muM && i <= j + mlM)
              E1d[k] += fac1 * SM_ELEMENT_B(M, i, j);
          }
        }
      }
    }
  }
  else if (rmem->mat_id == SUNMATRIX_BAND)
  {
    /* --- Band case --- */
    sunindextype mu = rmem->mu;
    sunindextype ml = rmem->ml;
    sunindextype i1, i2;

    SUNMatZero(E1);

    for (j = 0; j < n; j++)
    {
      i1 = SUNMAX(0, j - mu);
      i2 = SUNMIN(j + ml, n - 1);
      for (i = i1; i <= i2; i++)
        SM_ELEMENT_B(E1, i, j) = -SM_ELEMENT_B(J, i, j);

      if (M == NULL)
      {
        SM_ELEMENT_B(E1, j, j) += fac1;
      }
      else
      {
        /* M may be band or dense — check its type */
        if (SUNMatGetID(M) == SUNMATRIX_BAND)
        {
          sunindextype muM = SM_UBAND_B(M);
          sunindextype mlM = SM_LBAND_B(M);
          sunindextype im1 = SUNMAX(0, j - muM);
          sunindextype im2 = SUNMIN(j + mlM, n - 1);
          for (i = im1; i <= im2; i++)
            SM_ELEMENT_B(E1, i, j) += fac1 * SM_ELEMENT_B(M, i, j);
        }
        else
        {
          /* Dense M with band E1 — only copy within band */
          for (i = i1; i <= i2; i++)
            SM_ELEMENT_B(E1, i, j) += fac1 * SM_ELEMENT_D(M, i, j);
        }
      }
    }
  }
  else
  {
    /* --- Dense case --- */
    for (j = 0; j < n; j++)
    {
      for (i = 0; i < n; i++)
        SM_ELEMENT_D(E1, i, j) = -SM_ELEMENT_D(J, i, j);

      if (M == NULL)
      {
        SM_ELEMENT_D(E1, j, j) += fac1;
      }
      else
      {
        for (i = 0; i < n; i++)
          SM_ELEMENT_D(E1, i, j) += fac1 * SM_ELEMENT_D(M, i, j);
      }
    }
  }

  return RADAU5_SUCCESS;
}

/* ---------------------------------------------------------------------------
 * radau5_BuildE2
 *
 * Build the 2n×2n system for the coupled block.
 *
 * Eigenvalue mode (use_schur==0):
 *   E2 = [ alphn*M - J,   -betan*M ]
 *        [ betan*M,     alphn*M - J ]
 *
 * Schur mode (use_schur==1):
 *   E2 = [ TS[0][0]/h*M - J,   TS[0][1]/h*M       ]
 *        [ TS[1][0]/h*M,        TS[1][1]/h*M - J   ]
 *
 * (or with I in place of M when M is NULL)
 * J may be dense, band, or sparse. E2 is dense for dense/band J,
 * sparse (CSC) for sparse J.
 * ---------------------------------------------------------------------------*/
int radau5_BuildE2(Radau5Mem rmem, sunrealtype alphn, sunrealtype betan)
{
  sunindextype i, j, k, n;
  sunrealtype  jij, mij;
  SUNMatrix    J, M, E2;

  /* For Schur mode, the 2×2 block has 4 independent coefficients.
   * For eigenvalue mode: a00 = a11 = alphn, a01 = -betan, a10 = betan. */
  sunrealtype a00, a01, a10, a11;

  n  = rmem->n;
  J  = rmem->J;
  M  = rmem->M;
  E2 = rmem->E2;

  if (rmem->use_schur)
  {
    sunrealtype h = rmem->h;
    a00 = rmem->TS[0][0] / h;
    a01 = rmem->TS[0][1] / h;
    a10 = rmem->TS[1][0] / h;
    a11 = rmem->TS[1][1] / h;
  }
  else
  {
    a00 = alphn;
    a01 = -betan;
    a10 = betan;
    a11 = alphn;
  }

  if (rmem->mat_id == SUNMATRIX_SPARSE)
  {
    /* --- Sparse (CSC) case ---
     * Build the 2n×2n CSC matrix.
     * Column layout of E2 (2n columns):
     *   col j     (j=0..n-1): top-left block + bottom-left block
     *   col j+n   (j=0..n-1): top-right block + bottom-right block
     *
     * When M is sparse, diagonal blocks use union(J,M) pattern and
     * off-diagonal blocks use M's pattern. When M is NULL or dense,
     * the original approach is used.
     */
    sunindextype *Jp = SM_INDEXPTRS_S(J);
    sunindextype *Ji = SM_INDEXVALS_S(J);
    sunrealtype  *Jd = SM_DATA_S(J);

    sunindextype *E2p = SM_INDEXPTRS_S(E2);
    sunindextype *E2i = SM_INDEXVALS_S(E2);
    sunrealtype  *E2d = SM_DATA_S(E2);

    sunindextype pos = 0;
    int m_is_sparse = (M != NULL && SUNMatGetID(M) == SUNMATRIX_SPARSE);

    if (m_is_sparse)
    {
      /* E1 already holds the union(J,M) pattern from lazy init */
      SUNMatrix E1 = rmem->E1;
      sunindextype *Up = SM_INDEXPTRS_S(E1);  /* union pattern */
      sunindextype *Ui = SM_INDEXVALS_S(E1);
      sunindextype *Mp = SM_INDEXPTRS_S(M);
      sunindextype *Mi = SM_INDEXVALS_S(M);

      for (j = 0; j < n; j++)
      {
        E2p[j] = pos;

        /* Top-left block: a00*M[i,j] - J[i,j] over union pattern */
        for (k = Up[j]; k < Up[j+1]; k++)
        {
          i = Ui[k];
          E2i[pos] = i;
          E2d[pos] = a00 * radau5_SparseLookup(M, i, j)
                   - radau5_SparseLookup(J, i, j);
          pos++;
        }

        /* Bottom-left block: a10*M[i,j] over M's pattern */
        for (k = Mp[j]; k < Mp[j+1]; k++)
        {
          E2i[pos] = Mi[k] + n;
          E2d[pos] = a10 * SM_DATA_S(M)[k];
          pos++;
        }
      }

      for (j = 0; j < n; j++)
      {
        E2p[j + n] = pos;

        /* Top-right block: a01*M[i,j] over M's pattern */
        for (k = Mp[j]; k < Mp[j+1]; k++)
        {
          E2i[pos] = Mi[k];
          E2d[pos] = a01 * SM_DATA_S(M)[k];
          pos++;
        }

        /* Bottom-right block: a11*M[i,j] - J[i,j] over union pattern */
        for (k = Up[j]; k < Up[j+1]; k++)
        {
          i = Ui[k];
          E2i[pos] = i + n;
          E2d[pos] = a11 * radau5_SparseLookup(M, i, j)
                   - radau5_SparseLookup(J, i, j);
          pos++;
        }
      }

      E2p[2 * n] = pos;
    }
    else
    {
      /* M is NULL or dense — original approach */
      for (j = 0; j < n; j++)
      {
        E2p[j] = pos;

        /* Top-left block: a00*M[i,j] - J[i,j] for each (i,j) in J's pattern */
        for (k = Jp[j]; k < Jp[j+1]; k++)
        {
          i = Ji[k];
          mij = (M == NULL) ? ((i == j) ? SUN_RCONST(1.0) : SUN_RCONST(0.0))
                            : SM_ELEMENT_D(M, i, j);
          E2i[pos] = i;
          E2d[pos] = a00 * mij - Jd[k];
          pos++;
        }

        /* Bottom-left block: a10*M[i,j] */
        if (M == NULL)
        {
          E2i[pos] = j + n;
          E2d[pos] = a10;
          pos++;
        }
        else
        {
          for (i = 0; i < n; i++)
          {
            mij = SM_ELEMENT_D(M, i, j);
            if (mij != SUN_RCONST(0.0))
            {
              E2i[pos] = i + n;
              E2d[pos] = a10 * mij;
              pos++;
            }
          }
        }
      }

      for (j = 0; j < n; j++)
      {
        E2p[j + n] = pos;

        /* Top-right block: a01*M[i,j] */
        if (M == NULL)
        {
          E2i[pos] = j;
          E2d[pos] = a01;
          pos++;
        }
        else
        {
          for (i = 0; i < n; i++)
          {
            mij = SM_ELEMENT_D(M, i, j);
            if (mij != SUN_RCONST(0.0))
            {
              E2i[pos] = i;
              E2d[pos] = a01 * mij;
              pos++;
            }
          }
        }

        /* Bottom-right block: a11*M[i,j] - J[i,j] */
        for (k = Jp[j]; k < Jp[j+1]; k++)
        {
          i = Ji[k];
          mij = (M == NULL) ? ((i == j) ? SUN_RCONST(1.0) : SUN_RCONST(0.0))
                            : SM_ELEMENT_D(M, i, j);
          E2i[pos] = i + n;
          E2d[pos] = a11 * mij - Jd[k];
          pos++;
        }
      }

      E2p[2 * n] = pos;
    }
  }
  else
  {
    /* --- Dense / Band case (E2 is always dense) --- */
    int j_is_band = (rmem->mat_id == SUNMATRIX_BAND);
    sunindextype mu = rmem->mu;
    sunindextype ml = rmem->ml;

    for (j = 0; j < n; j++)
    {
      for (i = 0; i < n; i++)
      {
        /* Read J(i,j) */
        if (j_is_band)
          jij = (i >= j - mu && i <= j + ml) ? SM_ELEMENT_B(J, i, j)
                                              : SUN_RCONST(0.0);
        else
          jij = SM_ELEMENT_D(J, i, j);

        /* Read M(i,j) or identity */
        if (M == NULL)
        {
          mij = (i == j) ? SUN_RCONST(1.0) : SUN_RCONST(0.0);
        }
        else if (SUNMatGetID(M) == SUNMATRIX_BAND)
        {
          sunindextype muM = SM_UBAND_B(M);
          sunindextype mlM = SM_LBAND_B(M);
          mij = (i >= j - muM && i <= j + mlM) ? SM_ELEMENT_B(M, i, j)
                                                 : SUN_RCONST(0.0);
        }
        else
        {
          mij = SM_ELEMENT_D(M, i, j);
        }

        SM_ELEMENT_D(E2, i,     j    ) = a00 * mij - jij;
        SM_ELEMENT_D(E2, i,     j + n) = a01 * mij;
        SM_ELEMENT_D(E2, i + n, j    ) = a10 * mij;
        SM_ELEMENT_D(E2, i + n, j + n) = a11 * mij - jij;
      }
    }
  }

  return RADAU5_SUCCESS;
}

/* ---------------------------------------------------------------------------
 * radau5_DecompE1
 *
 * Factor E1 via SUNLinSolSetup. Increments ndec.
 * ---------------------------------------------------------------------------*/
int radau5_DecompE1(Radau5Mem rmem)
{
  int retval;

  retval = SUNLinSolSetup(rmem->LS_E1, rmem->E1);
  rmem->ndec++;

  return (retval == 0) ? RADAU5_SUCCESS : RADAU5_SINGULAR_MATRIX;
}

/* ---------------------------------------------------------------------------
 * radau5_DecompE2
 *
 * Factor E2 via SUNLinSolSetup. Increments ndec.
 * ---------------------------------------------------------------------------*/
int radau5_DecompE2(Radau5Mem rmem)
{
  int retval;

  retval = SUNLinSolSetup(rmem->LS_E2, rmem->E2);
  rmem->ndec++;

  return (retval == 0) ? RADAU5_SUCCESS : RADAU5_SINGULAR_MATRIX;
}

/* ---------------------------------------------------------------------------
 * radau5_ComputeScal
 *
 * Compute scal[i] = atol[i] + rtol[i] * |y[i]|
 * Matches Fortran radau5.f lines 780-788.
 * ---------------------------------------------------------------------------*/
void radau5_ComputeScal(Radau5Mem rmem, N_Vector y)
{
  sunindextype  i, n;
  sunrealtype  *y_data, *scal_data;
  sunrealtype  *rtol_v_data, *atol_v_data;

  n         = rmem->n;
  y_data    = N_VGetArrayPointer(y);
  scal_data = N_VGetArrayPointer(rmem->scal);

  if (rmem->itol == 0)
  {
    /* Scalar tolerances */
    for (i = 0; i < n; i++)
      scal_data[i] = rmem->atol_s + rmem->rtol_s * fabs(y_data[i]);
  }
  else
  {
    /* Vector tolerances */
    rtol_v_data = N_VGetArrayPointer(rmem->rtol_v);
    atol_v_data = N_VGetArrayPointer(rmem->atol_v);
    for (i = 0; i < n; i++)
      scal_data[i] = atol_v_data[i] + rtol_v_data[i] * fabs(y_data[i]);
  }
}

/* ---------------------------------------------------------------------------
 * radau5_MassMult
 *
 * Compute result = M * x, dispatching on M's matrix type (dense or band).
 * ---------------------------------------------------------------------------*/
int radau5_MassMult(Radau5Mem rmem, N_Vector x, N_Vector result)
{
  sunindextype i, j, n = rmem->n;
  SUNMatrix M = rmem->M;
  sunrealtype *xd = N_VGetArrayPointer(x);
  sunrealtype *rd = N_VGetArrayPointer(result);

  if (M == NULL)
  {
    /* Identity mass — just copy */
    N_VScale(SUN_RCONST(1.0), x, result);
    return RADAU5_SUCCESS;
  }

  SUNMatrix_ID mid = SUNMatGetID(M);

  if (mid == SUNMATRIX_DENSE)
  {
    for (i = 0; i < n; i++)
    {
      sunrealtype sum = SUN_RCONST(0.0);
      for (j = 0; j < n; j++)
        sum += SM_ELEMENT_D(M, i, j) * xd[j];
      rd[i] = sum;
    }
  }
  else if (mid == SUNMATRIX_BAND)
  {
    sunindextype muM = SM_UBAND_B(M);
    sunindextype mlM = SM_LBAND_B(M);
    for (i = 0; i < n; i++)
    {
      sunrealtype sum = SUN_RCONST(0.0);
      sunindextype j1 = SUNMAX(0, i - mlM);
      sunindextype j2 = SUNMIN(n - 1, i + muM);
      for (j = j1; j <= j2; j++)
        sum += SM_ELEMENT_B(M, i, j) * xd[j];
      rd[i] = sum;
    }
  }
  else
  {
    /* Fallback: use SUNMatMatvec if available */
    SUNErrCode err = SUNMatMatvec(M, x, result);
    if (err != SUN_SUCCESS) return RADAU5_LSOLVE_FAIL;
  }

  return RADAU5_SUCCESS;
}
