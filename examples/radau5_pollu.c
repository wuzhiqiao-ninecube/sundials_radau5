/* ---------------------------------------------------------------------------
 * radau5_pollu.c — Pollution problem (IVPtestset)
 *
 * ODE of dimension 20, 25 reaction rates, sparse Jacobian (KLU).
 * t in [0, 60],  reference solution from pollu.f (solut subroutine).
 *
 * Reference: Test Set for IVP Solvers
 *   http://www.dm.uniba.it/~testset/
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include "radau5.h"

#define NEQ 20
#define NNZ 86  /* 82 Jacobian entries + 4 missing diagonal entries (cols 7,11,14,17) */

/* Rate constants */
#define K1   0.35e+00
#define K2   0.266e+02
#define K3   0.123e+05
#define K4   0.86e-03
#define K5   0.82e-03
#define K6   0.15e+05
#define K7   0.13e-03
#define K8   0.24e+05
#define K9   0.165e+05
#define K10  0.9e+04
#define K11  0.22e-01
#define K12  0.12e+05
#define K13  0.188e+01
#define K14  0.163e+05
#define K15  0.48e+07
#define K16  0.35e-03
#define K17  0.175e-01
#define K18  0.1e+09
#define K19  0.444e+12
#define K20  0.124e+04
#define K21  0.21e+01
#define K22  0.578e+01
#define K23  0.474e-01
#define K24  0.178e+04
#define K25  0.312e+01

/* ---------------------------------------------------------------------------
 * RHS
 * ---------------------------------------------------------------------------*/
static int rhs(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)t; (void)ud;
  sunrealtype* yv = N_VGetArrayPointer(y);
  sunrealtype* f  = N_VGetArrayPointer(yd);

  sunrealtype r1  = K1  * yv[0];
  sunrealtype r2  = K2  * yv[1]  * yv[3];
  sunrealtype r3  = K3  * yv[4]  * yv[1];
  sunrealtype r4  = K4  * yv[6];
  sunrealtype r5  = K5  * yv[6];
  sunrealtype r6  = K6  * yv[6]  * yv[5];
  sunrealtype r7  = K7  * yv[8];
  sunrealtype r8  = K8  * yv[8]  * yv[5];
  sunrealtype r9  = K9  * yv[10] * yv[1];
  sunrealtype r10 = K10 * yv[10] * yv[0];
  sunrealtype r11 = K11 * yv[12];
  sunrealtype r12 = K12 * yv[9]  * yv[1];
  sunrealtype r13 = K13 * yv[13];
  sunrealtype r14 = K14 * yv[0]  * yv[5];
  sunrealtype r15 = K15 * yv[2];
  sunrealtype r16 = K16 * yv[3];
  sunrealtype r17 = K17 * yv[3];
  sunrealtype r18 = K18 * yv[15];
  sunrealtype r19 = K19 * yv[15];
  sunrealtype r20 = K20 * yv[16] * yv[5];
  sunrealtype r21 = K21 * yv[18];
  sunrealtype r22 = K22 * yv[18];
  sunrealtype r23 = K23 * yv[0]  * yv[3];
  sunrealtype r24 = K24 * yv[18] * yv[0];
  sunrealtype r25 = K25 * yv[19];

  f[0]  = -r1  - r10 - r14 - r23 - r24 + r2  + r3  + r9  + r11 + r12 + r22 + r25;
  f[1]  = -r2  - r3  - r9  - r12 + r1  + r21;
  f[2]  = -r15 + r1  + r17 + r19 + r22;
  f[3]  = -r2  - r16 - r17 - r23 + r15;
  f[4]  = -r3  + r4  + r4  + r6  + r7  + r13 + r20;
  f[5]  = -r6  - r8  - r14 - r20 + r3  + r18 + r18;
  f[6]  = -r4  - r5  - r6  + r13;
  f[7]  =  r4  + r5  + r6  + r7;
  f[8]  = -r7  - r8;
  f[9]  = -r12 + r7  + r9;
  f[10] = -r9  - r10 + r8  + r11;
  f[11] =  r9;
  f[12] = -r11 + r10;
  f[13] = -r13 + r12;
  f[14] =  r14;
  f[15] = -r18 - r19 + r16;
  f[16] = -r20;
  f[17] =  r20;
  f[18] = -r21 - r22 - r24 + r23 + r25;
  f[19] = -r25 + r24;

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
 * Analytic Jacobian — fills the pre-allocated CSC sparse matrix
 * ---------------------------------------------------------------------------*/
static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;
  sunrealtype* yv = N_VGetArrayPointer(y);

  /* Zero data only — do NOT call SUNMatZero which destroys CSC structure */
  sunrealtype* data = SM_DATA_S(J);
  sunindextype* colptrs = SM_INDEXPTRS_S(J);
  sunindextype nnz_used = colptrs[NEQ];
  for (sunindextype k = 0; k < nnz_used; k++) data[k] = 0.0;

  /* col 0 (y1) */
  sparse_set(J,  0, 0, -K1 - K10*yv[10] - K14*yv[5] - K23*yv[3] - K24*yv[18]);
  sparse_set(J,  1, 0,  K1);
  sparse_set(J,  2, 0,  K1);
  sparse_set(J,  3, 0, -K23*yv[3]);
  sparse_set(J,  5, 0, -K14*yv[5]);
  sparse_set(J, 10, 0, -K10*yv[10]);
  sparse_set(J, 12, 0,  K10*yv[10]);
  sparse_set(J, 14, 0,  K14*yv[5]);
  sparse_set(J, 18, 0, -K24*yv[18] + K23*yv[3]);
  sparse_set(J, 19, 0,  K24*yv[18]);

  /* col 1 (y2) */
  sparse_set(J,  0, 1,  K2*yv[3] + K9*yv[10] + K3*yv[4] + K12*yv[9]);
  sparse_set(J,  1, 1, -K2*yv[3] - K3*yv[4] - K9*yv[10] - K12*yv[9]);
  sparse_set(J,  3, 1, -K2*yv[3]);
  sparse_set(J,  4, 1, -K3*yv[4]);
  sparse_set(J,  5, 1,  K3*yv[4]);
  sparse_set(J,  9, 1, -K12*yv[9] + K9*yv[10]);
  sparse_set(J, 10, 1, -K9*yv[10]);
  sparse_set(J, 11, 1,  K9*yv[10]);
  sparse_set(J, 13, 1,  K12*yv[9]);

  /* col 2 (y3) */
  sparse_set(J,  2, 2, -K15);
  sparse_set(J,  3, 2,  K15);

  /* col 3 (y4) */
  sparse_set(J,  0, 3, -K23*yv[0] + K2*yv[1]);
  sparse_set(J,  1, 3, -K2*yv[1]);
  sparse_set(J,  2, 3,  K17);
  sparse_set(J,  3, 3, -K2*yv[1] - K16 - K17 - K23*yv[0]);
  sparse_set(J, 15, 3,  K16);
  sparse_set(J, 18, 3,  K23*yv[0]);

  /* col 4 (y5) */
  sparse_set(J,  0, 4,  K3*yv[1]);
  sparse_set(J,  1, 4, -K3*yv[1]);
  sparse_set(J,  4, 4, -K3*yv[1]);
  sparse_set(J,  5, 4,  K3*yv[1]);

  /* col 5 (y6) */
  sparse_set(J,  0, 5, -K14*yv[0]);
  sparse_set(J,  4, 5,  K6*yv[6] + K20*yv[16]);
  sparse_set(J,  5, 5, -K6*yv[6] - K8*yv[8] - K14*yv[0] - K20*yv[16]);
  sparse_set(J,  6, 5, -K6*yv[6]);
  sparse_set(J,  7, 5,  K6*yv[6]);
  sparse_set(J,  8, 5, -K8*yv[8]);
  sparse_set(J, 10, 5,  K8*yv[8]);
  sparse_set(J, 14, 5,  K14*yv[0]);
  sparse_set(J, 16, 5, -K20*yv[16]);
  sparse_set(J, 17, 5,  K20*yv[16]);

  /* col 6 (y7) */
  sparse_set(J,  4, 6,  2.0*K4 + K6*yv[5]);
  sparse_set(J,  5, 6, -K6*yv[5]);
  sparse_set(J,  6, 6, -K4 - K5 - K6*yv[5]);
  sparse_set(J,  7, 6,  K4 + K5 + K6*yv[5]);

  /* col 7 (y8): no nonzeros */

  /* col 8 (y9) */
  sparse_set(J,  4, 8,  K7);
  sparse_set(J,  5, 8, -K8*yv[5]);
  sparse_set(J,  7, 8,  K7);
  sparse_set(J,  8, 8, -K7 - K8*yv[5]);
  sparse_set(J,  9, 8,  K7);
  sparse_set(J, 10, 8,  K8*yv[5]);

  /* col 9 (y10) */
  sparse_set(J,  0, 9,  K12*yv[1]);
  sparse_set(J,  1, 9, -K12*yv[1]);
  sparse_set(J,  9, 9, -K12*yv[1]);
  sparse_set(J, 13, 9,  K12*yv[1]);

  /* col 10 (y11) */
  sparse_set(J,  0, 10, -K10*yv[0] + K9*yv[1]);
  sparse_set(J,  1, 10, -K9*yv[1]);
  sparse_set(J,  9, 10,  K9*yv[1]);
  sparse_set(J, 10, 10, -K9*yv[1] - K10*yv[0]);
  sparse_set(J, 11, 10,  K9*yv[1]);
  sparse_set(J, 12, 10,  K10*yv[0]);

  /* col 11 (y12): no nonzeros */

  /* col 12 (y13) */
  sparse_set(J,  0, 12,  K11);
  sparse_set(J, 10, 12,  K11);
  sparse_set(J, 12, 12, -K11);

  /* col 13 (y14) */
  sparse_set(J,  4, 13,  K13);
  sparse_set(J,  6, 13,  K13);
  sparse_set(J, 13, 13, -K13);

  /* col 14 (y15): no nonzeros */

  /* col 15 (y16) */
  sparse_set(J,  2, 15,  K19);
  sparse_set(J,  5, 15,  2.0*K18);
  sparse_set(J, 15, 15, -K18 - K19);

  /* col 16 (y17) */
  sparse_set(J,  4, 16,  K20*yv[5]);
  sparse_set(J,  5, 16, -K20*yv[5]);
  sparse_set(J, 16, 16, -K20*yv[5]);
  sparse_set(J, 17, 16,  K20*yv[5]);

  /* col 17 (y18): no nonzeros */

  /* col 18 (y19) */
  sparse_set(J,  0, 18, -K24*yv[0] + K22);
  sparse_set(J,  1, 18,  K21);
  sparse_set(J,  2, 18,  K22);
  sparse_set(J, 18, 18, -K21 - K22 - K24*yv[0]);
  sparse_set(J, 19, 18,  K24*yv[0]);

  /* col 19 (y20) */
  sparse_set(J,  0, 19,  K25);
  sparse_set(J, 18, 19,  K25);
  sparse_set(J, 19, 19, -K25);

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
  int use_dq       = 0;  /* 1 = sparse DQ Jacobian via column grouping */
  if (argc > 1) rtol      = atof(argv[1]);
  if (argc > 2) atol_val  = atof(argv[2]);
  if (argc > 3) h0        = atof(argv[3]);
  if (argc > 4) use_schur = atoi(argv[4]);
  if (argc > 5) use_dq    = atoi(argv[5]);

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  /* initial conditions (0-based, from pollu.f init) */
  N_Vector y0 = N_VNew_Serial(NEQ, sunctx);
  N_VConst(0.0, y0);
  sunrealtype* y0v = N_VGetArrayPointer(y0);
  y0v[1]  = 0.2;
  y0v[3]  = 0.04;
  y0v[6]  = 0.1;
  y0v[7]  = 0.3;
  y0v[8]  = 0.01;
  y0v[16] = 0.007;

  /* sparse Jacobian template: CSC, 20x20, NNZ=82 */
  SUNMatrix Jt = SUNSparseMatrix(NEQ, NEQ, NNZ, CSC_MAT, sunctx);

  /* column pointers (21 values) — includes diagonal entries for cols 7,11,14,17 */
  static const sunindextype colptrs[NEQ + 1] = {
    0, 10, 19, 21, 27, 31, 41, 45, 46, 52, 56, 62, 63, 66, 69, 70, 73, 77, 78, 83, 86
  };
  /* row indices (86 values) — sorted within each column */
  static const sunindextype rowinds[NNZ] = {
     0,  1,  2,  3,  5, 10, 12, 14, 18, 19,  /* col  0: 10 */
     0,  1,  3,  4,  5,  9, 10, 11, 13,       /* col  1:  9 */
     2,  3,                                    /* col  2:  2 */
     0,  1,  2,  3, 15, 18,                   /* col  3:  6 */
     0,  1,  4,  5,                            /* col  4:  4 */
     0,  4,  5,  6,  7,  8, 10, 14, 16, 17,  /* col  5: 10 */
     4,  5,  6,  7,                            /* col  6:  4 */
     7,                                        /* col  7:  1 (diagonal only) */
     4,  5,  7,  8,  9, 10,                   /* col  8:  6 */
     0,  1,  9, 13,                            /* col  9:  4 */
     0,  1,  9, 10, 11, 12,                   /* col 10:  6 */
    11,                                        /* col 11:  1 (diagonal only) */
     0, 10, 12,                               /* col 12:  3 */
     4,  6, 13,                               /* col 13:  3 */
    14,                                        /* col 14:  1 (diagonal only) */
     2,  5, 15,                               /* col 15:  3 */
     4,  5, 16, 17,                           /* col 16:  4 */
    17,                                        /* col 17:  1 (diagonal only) */
     0,  1,  2, 18, 19,                       /* col 18:  5 */
     0, 18, 19                                /* col 19:  3 */
  };

  sunindextype* cp = SM_INDEXPTRS_S(Jt);
  sunindextype* ri = SM_INDEXVALS_S(Jt);
  for (int i = 0; i <= NEQ; i++) cp[i] = colptrs[i];
  for (int k = 0; k < NNZ;  k++) ri[k] = rowinds[k];
  sunrealtype* dat = SM_DATA_S(Jt);
  for (int k = 0; k < NNZ;  k++) dat[k] = 0.0;

  /* solver setup */
  void* mem = Radau5Create(sunctx);
  Radau5Init(mem, rhs, 0.0, y0);
  Radau5SetLinearSolver(mem, Jt);
  Radau5SetSchurDecomp(mem, use_schur);
  if (use_dq) {
    /* Sparse DQ Jacobian: register sparsity pattern, skip analytic jac */
    Radau5SetSparsityPattern(mem, Jt);
  } else {
    Radau5SetJacFn(mem, jac);
  }
  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);

  /* solve to t=60 */
  N_Vector yout = N_VNew_Serial(NEQ, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 60.0, yout, &tret);

  /* reference solution from pollu.f solut */
  static const sunrealtype yref[NEQ] = {
    0.5646255480022769e-01,
    0.1342484130422339e+00,
    0.4139734331099427e-08,
    0.5523140207484359e-02,
    0.2018977262302196e-06,
    0.1464541863493966e-06,
    0.7784249118997964e-01,
    0.3245075353396018e+00,
    0.7494013383880406e-02,
    0.1622293157301561e-07,
    0.1135863833257075e-07,
    0.2230505975721359e-02,
    0.2087162882798630e-03,
    0.1396921016840158e-04,
    0.8964884856898295e-02,
    0.4352846369330103e-17,
    0.6899219696263405e-02,
    0.1007803037365946e-03,
    0.1772146513969984e-05,
    0.5682943292316392e-04
  };

  sunrealtype* yd = N_VGetArrayPointer(yout);
  printf("=== Pollution problem (rtol=%.1e atol=%.1e h0=%.1e schur=%d dq=%d) ===\n", rtol, atol_val, h0, use_schur, use_dq);
  printf("ret=%d  tret=%.6e\n\n", ret, tret);
  printf("%-4s  %-22s  %-22s  %-10s\n", "i", "y[i]", "ref[i]", "rel_err");

  int pass = (ret == 0);
  for (int i = 0; i < NEQ; i++) {
    sunrealtype err    = fabs(yd[i] - yref[i]);
    sunrealtype relerr = (fabs(yref[i]) > 0.0) ? err / fabs(yref[i]) : err;
    printf("[%2d]  %22.14e  %22.14e  %.3e\n", i, yd[i], yref[i], relerr);
    if (relerr > 1.0e-3) pass = 0;
  }

  long int nstep, naccpt, nrejct, nfcn, njac, ndec;
  Radau5GetNumSteps(mem, &nstep);
  Radau5GetNumAccSteps(mem, &naccpt);
  Radau5GetNumRejSteps(mem, &nrejct);
  Radau5GetNumRhsEvals(mem, &nfcn);
  Radau5GetNumJacEvals(mem, &njac);
  Radau5GetNumDecomps(mem, &ndec);
  printf("\nnstep=%ld  naccpt=%ld  nrejct=%ld  nfcn=%ld  njac=%ld  ndec=%ld\n",
         nstep, naccpt, nrejct, nfcn, njac, ndec);
  printf("%s\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout); SUNMatDestroy(Jt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
