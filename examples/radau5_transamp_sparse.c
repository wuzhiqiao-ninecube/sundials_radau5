/* ---------------------------------------------------------------------------
 * radau5_transamp_sparse.c — Transistor Amplifier (DAE index-1, n=8)
 *
 * Same problem as radau5_transamp.c but with BOTH Jacobian and mass matrix
 * stored as sparse CSC matrices (SUNSparseMatrix + KLU).
 *
 * M*y' = f(t,y) where M is tridiagonal (mlmas=1, mumas=1), 14 NNZ.
 * Jacobian is banded (mljac=2, mujac=1), 16 NNZ.
 *
 * t in [0, 0.2],  reference solution from IVPtestset
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include "radau5.h"

#define NEQ   8
#define J_NNZ 16
#define M_NNZ 14

#define UB    6.0
#define UF    0.026
#define ALPHA 0.99
#define BETA  1.0e-6
#define R0    1000.0
#define R1    9000.0
#define R2    9000.0
#define R3    9000.0
#define R4    9000.0
#define R5    9000.0
#define R6    9000.0
#define R7    9000.0
#define R8    9000.0
#define R9    9000.0
#define C1    1.0e-6
#define C2    2.0e-6
#define C3    3.0e-6
#define C4    4.0e-6
#define C5    5.0e-6
#define PI    3.1415926535897931086244

/* ---------------------------------------------------------------------------
 * RHS
 * ---------------------------------------------------------------------------*/
static int rhs_transamp(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)ud;
  sunrealtype* v = N_VGetArrayPointer(y);
  sunrealtype* f = N_VGetArrayPointer(yd);

  sunrealtype uet = 0.1 * sin(200.0 * PI * t);
  sunrealtype e1 = (v[1] - v[2]) / UF;
  sunrealtype e2 = (v[4] - v[5]) / UF;
  if (e1 > 300.0 || e2 > 300.0) return 1;  /* recoverable: exp overflow */

  sunrealtype fac1 = BETA * (exp(e1) - 1.0);
  sunrealtype fac2 = BETA * (exp(e2) - 1.0);

  f[0] = (v[0] - uet) / R0;
  f[1] = v[1] / R1 + (v[1] - UB) / R2 + (1.0 - ALPHA) * fac1;
  f[2] = v[2] / R3 - fac1;
  f[3] = (v[3] - UB) / R4 + ALPHA * fac1;
  f[4] = v[4] / R5 + (v[4] - UB) / R6 + (1.0 - ALPHA) * fac2;
  f[5] = v[5] / R7 - fac2;
  f[6] = (v[6] - UB) / R8 + ALPHA * fac2;
  f[7] = v[7] / R9;
  return 0;
}

/* ---------------------------------------------------------------------------
 * Helper: set a value in a CSC sparse matrix by (row, col) linear search
 * ---------------------------------------------------------------------------*/
static void sparse_set(SUNMatrix A, sunindextype row, sunindextype col,
                       sunrealtype val)
{
  sunindextype* colptrs = SM_INDEXPTRS_S(A);
  sunindextype* rowinds = SM_INDEXVALS_S(A);
  sunrealtype*  data    = SM_DATA_S(A);
  for (sunindextype k = colptrs[col]; k < colptrs[col + 1]; k++) {
    if (rowinds[k] == row) { data[k] = val; return; }
  }
}

/* ---------------------------------------------------------------------------
 * Analytic Jacobian — sparse CSC (16 NNZ)
 * Sparsity: banded mljac=2, mujac=1
 *   col 0: row {0}
 *   col 1: rows {1,2,3}
 *   col 2: rows {1,2,3}
 *   col 3: row {3}
 *   col 4: rows {4,5,6}
 *   col 5: rows {4,5,6}
 *   col 6: row {6}
 *   col 7: row {7}
 * ---------------------------------------------------------------------------*/
static int jac_transamp(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                        void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;
  sunrealtype* v = N_VGetArrayPointer(y);

  /* Zero data only */
  sunrealtype* data = SM_DATA_S(J);
  sunindextype* cp = SM_INDEXPTRS_S(J);
  for (sunindextype k = 0; k < cp[NEQ]; k++) data[k] = 0.0;

  sunrealtype fac1p = BETA * exp((v[1] - v[2]) / UF) / UF;
  sunrealtype fac2p = BETA * exp((v[4] - v[5]) / UF) / UF;

  /* Diagonal */
  sparse_set(J, 0, 0, 1.0 / R0);
  sparse_set(J, 1, 1, 1.0/R1 + 1.0/R2 + (1.0-ALPHA)*fac1p);
  sparse_set(J, 2, 2, 1.0/R3 + fac1p);
  sparse_set(J, 3, 3, 1.0 / R4);
  sparse_set(J, 4, 4, 1.0/R5 + 1.0/R6 + (1.0-ALPHA)*fac2p);
  sparse_set(J, 5, 5, 1.0/R7 + fac2p);
  sparse_set(J, 6, 6, 1.0 / R8);
  sparse_set(J, 7, 7, 1.0 / R9);

  /* Superdiagonal (mujac=1) */
  sparse_set(J, 1, 2, -(1.0-ALPHA)*fac1p);
  sparse_set(J, 4, 5, -(1.0-ALPHA)*fac2p);

  /* First subdiagonal */
  sparse_set(J, 2, 1, -fac1p);
  sparse_set(J, 3, 2, -ALPHA*fac1p);
  sparse_set(J, 5, 4, -fac2p);
  sparse_set(J, 6, 5, -ALPHA*fac2p);

  /* Second subdiagonal (mljac=2) */
  sparse_set(J, 3, 1, ALPHA*fac1p);
  sparse_set(J, 6, 4, ALPHA*fac2p);

  return 0;
}

/* ---------------------------------------------------------------------------
 * Mass matrix — sparse CSC (14 NNZ)
 * Tridiagonal: mlmas=1, mumas=1
 *   col 0: rows {0,1}    col 1: rows {0,1}    col 2: row {2}
 *   col 3: rows {3,4}    col 4: rows {3,4}    col 5: row {5}
 *   col 6: rows {6,7}    col 7: rows {6,7}
 * ---------------------------------------------------------------------------*/
static int mas_transamp(sunrealtype t, SUNMatrix M, void* ud,
                        N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)ud; (void)t1; (void)t2; (void)t3;

  /* Zero data only — preserve CSC structure */
  sunrealtype* data = SM_DATA_S(M);
  sunindextype* cp = SM_INDEXPTRS_S(M);
  for (sunindextype k = 0; k < cp[NEQ]; k++) data[k] = 0.0;

  /* Main diagonal */
  sparse_set(M, 0, 0, -C1);
  sparse_set(M, 1, 1, -C1);
  sparse_set(M, 2, 2, -C2);
  sparse_set(M, 3, 3, -C3);
  sparse_set(M, 4, 4, -C3);
  sparse_set(M, 5, 5, -C4);
  sparse_set(M, 6, 6, -C5);
  sparse_set(M, 7, 7, -C5);

  /* Upper diagonal */
  sparse_set(M, 0, 1, C1);
  sparse_set(M, 3, 4, C3);
  sparse_set(M, 6, 7, C5);

  /* Lower diagonal */
  sparse_set(M, 1, 0, C1);
  sparse_set(M, 4, 3, C3);
  sparse_set(M, 7, 6, C5);

  return 0;
}

/* ---------------------------------------------------------------------------
 * main
 * ---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  sunrealtype rtol = 1.0e-6;
  sunrealtype atol_val = 1.0e-6;
  sunrealtype h0   = 1.0e-6;
  int use_schur    = 0;
  if (argc > 1) rtol      = atof(argv[1]);
  if (argc > 2) atol_val  = atof(argv[2]);
  if (argc > 3) h0        = atof(argv[3]);
  if (argc > 4) use_schur = atoi(argv[4]);

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  int n = NEQ;
  void* mem = Radau5Create(sunctx);

  /* Initial conditions */
  N_Vector y0 = N_VNew_Serial(n, sunctx);
  sunrealtype* y0d = N_VGetArrayPointer(y0);
  y0d[0] = 0.0;
  y0d[1] = UB / (R2 / R1 + 1.0);
  y0d[2] = y0d[1];
  y0d[3] = UB;
  y0d[4] = UB / (R6 / R5 + 1.0);
  y0d[5] = y0d[4];
  y0d[6] = y0d[3];
  y0d[7] = 0.0;

  Radau5Init(mem, rhs_transamp, 0.0, y0);

  /* Sparse Jacobian template (CSC, 8x8, 16 NNZ) */
  SUNMatrix Jt = SUNSparseMatrix(n, n, J_NNZ, CSC_MAT, sunctx);
  static const sunindextype j_colptrs[NEQ + 1] = {
    0, 1, 4, 7, 8, 11, 14, 15, 16
  };
  static const sunindextype j_rowinds[J_NNZ] = {
    0,          /* col 0 */
    1, 2, 3,    /* col 1 */
    1, 2, 3,    /* col 2 */
    3,          /* col 3 */
    4, 5, 6,    /* col 4 */
    4, 5, 6,    /* col 5 */
    6,          /* col 6 */
    7           /* col 7 */
  };

  {
    sunindextype* cp = SM_INDEXPTRS_S(Jt);
    sunindextype* ri = SM_INDEXVALS_S(Jt);
    for (int i = 0; i <= NEQ; i++) cp[i] = j_colptrs[i];
    for (int k = 0; k < J_NNZ; k++) ri[k] = j_rowinds[k];
    sunrealtype* dat = SM_DATA_S(Jt);
    for (int k = 0; k < J_NNZ; k++) dat[k] = 0.0;
  }

  Radau5SetLinearSolver(mem, Jt);
  Radau5SetSchurDecomp(mem, use_schur);
  Radau5SetJacFn(mem, jac_transamp);

  /* Sparse mass matrix template (CSC, 8x8, 14 NNZ) */
  SUNMatrix Mt = SUNSparseMatrix(n, n, M_NNZ, CSC_MAT, sunctx);
  static const sunindextype m_colptrs[NEQ + 1] = {
    0, 2, 4, 5, 7, 9, 10, 12, 14
  };
  static const sunindextype m_rowinds[M_NNZ] = {
    0, 1,       /* col 0 */
    0, 1,       /* col 1 */
    2,          /* col 2 */
    3, 4,       /* col 3 */
    3, 4,       /* col 4 */
    5,          /* col 5 */
    6, 7,       /* col 6 */
    6, 7        /* col 7 */
  };
  {
    sunindextype* cp = SM_INDEXPTRS_S(Mt);
    sunindextype* ri = SM_INDEXVALS_S(Mt);
    for (int i = 0; i <= NEQ; i++) cp[i] = m_colptrs[i];
    for (int k = 0; k < M_NNZ; k++) ri[k] = m_rowinds[k];
    sunrealtype* dat = SM_DATA_S(Mt);
    for (int k = 0; k < M_NNZ; k++) dat[k] = 0.0;
  }

  Radau5SetMassFn(mem, mas_transamp, Mt);

  /* All 8 variables are index-1 */
  Radau5SetDAEIndex(mem, n, 0, 0);

  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);

  N_Vector yout = N_VNew_Serial(n, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 0.2, yout, &tret);

  /* Reference solution at t=0.2 */
  sunrealtype yref[8] = {
    -0.5562145012262709e-02,
     0.3006522471903042e+01,
     0.2849958788608128e+01,
     0.2926422536206241e+01,
     0.2704617865010554e+01,
     0.2761837778393145e+01,
     0.4770927631616772e+01,
     0.1236995868091548e+01
  };
  sunrealtype* yd = N_VGetArrayPointer(yout);

  printf("=== Transistor Amplifier SPARSE "
         "(rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n",
         rtol, atol_val, h0, use_schur);
  printf("ret=%d, tret=%.6e\n", ret, tret);

  sunrealtype maxrelerr = 0.0;
  for (int i = 0; i < n; i++)
  {
    sunrealtype err    = fabs(yd[i] - yref[i]);
    sunrealtype relerr = (fabs(yref[i]) > 1e-15) ? err / fabs(yref[i]) : err;
    if (relerr > maxrelerr) maxrelerr = relerr;
    printf("y[%d] = %16.10e  ref = %16.10e  rel_err = %.3e\n",
           i, yd[i], yref[i], relerr);
  }

  long int nstep, naccpt, nrejct, nfcn, njac, ndec, nsol;
  Radau5GetNumSteps(mem, &nstep);
  Radau5GetNumAccSteps(mem, &naccpt);
  Radau5GetNumRejSteps(mem, &nrejct);
  Radau5GetNumRhsEvals(mem, &nfcn);
  Radau5GetNumJacEvals(mem, &njac);
  Radau5GetNumDecomps(mem, &ndec);
  Radau5GetNumLinSolves(mem, &nsol);
  printf("nstep=%ld naccpt=%ld nrejct=%ld nfcn=%ld njac=%ld ndec=%ld nsol=%ld\n",
         nstep, naccpt, nrejct, nfcn, njac, ndec, nsol);

  int pass = (ret == 0 && maxrelerr < 1.0e-3);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout);
  SUNMatDestroy(Jt); SUNMatDestroy(Mt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
