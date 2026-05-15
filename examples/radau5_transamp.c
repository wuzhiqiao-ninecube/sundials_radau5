/* ---------------------------------------------------------------------------
 * radau5_transamp.c — Transistor Amplifier (DAE index-1, n=8)
 *
 * M*y' = f(t,y) where M is banded (mlmas=1, mumas=1), stored dense here.
 * Banded Jacobian (mljac=2, mujac=1), stored dense here (n=8 is small).
 *
 * t in [0, 0.2],  reference solution from IVPtestset
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

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

static int rhs_transamp(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)ud;
  sunrealtype* v = N_VGetArrayPointer(y);
  sunrealtype* f = N_VGetArrayPointer(yd);

  sunrealtype uet = 0.1 * sin(200.0 * PI * t);

  sunrealtype e1 = (v[1] - v[2]) / UF;
  sunrealtype e2 = (v[4] - v[5]) / UF;

  /* Overflow protection */
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

static int jac_transamp(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                        void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;
  sunrealtype* v = N_VGetArrayPointer(y);
  sunindextype i, j;

  /* Zero the Jacobian */
  for (j = 0; j < 8; j++)
    for (i = 0; i < 8; i++)
      SM_ELEMENT_D(J, i, j) = 0.0;

  sunrealtype fac1p = BETA * exp((v[1] - v[2]) / UF) / UF;
  sunrealtype fac2p = BETA * exp((v[4] - v[5]) / UF) / UF;

  /* Main diagonal */
  SM_ELEMENT_D(J, 0, 0) = 1.0 / R0;
  SM_ELEMENT_D(J, 1, 1) = 1.0 / R1 + 1.0 / R2 + (1.0 - ALPHA) * fac1p;
  SM_ELEMENT_D(J, 2, 2) = 1.0 / R3 + fac1p;
  SM_ELEMENT_D(J, 3, 3) = 1.0 / R4;
  SM_ELEMENT_D(J, 4, 4) = 1.0 / R5 + 1.0 / R6 + (1.0 - ALPHA) * fac2p;
  SM_ELEMENT_D(J, 5, 5) = 1.0 / R7 + fac2p;
  SM_ELEMENT_D(J, 6, 6) = 1.0 / R8;
  SM_ELEMENT_D(J, 7, 7) = 1.0 / R9;

  /* Superdiagonal (mujac=1): dfdy(1,j) → J[j-2, j-1] 0-indexed */
  SM_ELEMENT_D(J, 1, 2) = -(1.0 - ALPHA) * fac1p;  /* dfdy(1,3) */
  SM_ELEMENT_D(J, 4, 5) = -(1.0 - ALPHA) * fac2p;  /* dfdy(1,6) */

  /* First subdiagonal (mljac>=1): dfdy(3,j) → J[j, j-1] 0-indexed */
  SM_ELEMENT_D(J, 2, 1) = -fac1p;                   /* dfdy(3,2) */
  SM_ELEMENT_D(J, 3, 2) = -ALPHA * fac1p;            /* dfdy(3,3) */
  SM_ELEMENT_D(J, 5, 4) = -fac2p;                   /* dfdy(3,5) */
  SM_ELEMENT_D(J, 6, 5) = -ALPHA * fac2p;            /* dfdy(3,6) */

  /* Second subdiagonal (mljac=2): dfdy(4,j) → J[j+1, j-1] 0-indexed */
  SM_ELEMENT_D(J, 3, 1) = ALPHA * fac1p;             /* dfdy(4,2) */
  SM_ELEMENT_D(J, 6, 4) = ALPHA * fac2p;             /* dfdy(4,5) */

  return 0;
}

static int mas_transamp(sunrealtype t, SUNMatrix M, void* ud,
                        N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)ud; (void)t1; (void)t2; (void)t3;
  SUNMatZero(M);

  /* Main diagonal */
  SM_ELEMENT_D(M, 0, 0) = -C1;
  SM_ELEMENT_D(M, 1, 1) = -C1;
  SM_ELEMENT_D(M, 2, 2) = -C2;
  SM_ELEMENT_D(M, 3, 3) = -C3;
  SM_ELEMENT_D(M, 4, 4) = -C3;
  SM_ELEMENT_D(M, 5, 5) = -C4;
  SM_ELEMENT_D(M, 6, 6) = -C5;
  SM_ELEMENT_D(M, 7, 7) = -C5;

  /* Upper diagonal: M[i-1, i] */
  SM_ELEMENT_D(M, 0, 1) = C1;   /* dfddy(1,2) */
  SM_ELEMENT_D(M, 3, 4) = C3;   /* dfddy(1,5) */
  SM_ELEMENT_D(M, 6, 7) = C5;   /* dfddy(1,8) */

  /* Lower diagonal: M[i+1, i] */
  SM_ELEMENT_D(M, 1, 0) = C1;   /* dfddy(3,1) */
  SM_ELEMENT_D(M, 4, 3) = C3;   /* dfddy(3,4) */
  SM_ELEMENT_D(M, 7, 6) = C5;   /* dfddy(3,7) */

  return 0;
}

int main(int argc, char* argv[])
{
  sunrealtype rtol = 1.0e-6;
  sunrealtype atol_val = 1.0e-6;
  sunrealtype h0   = 1.0e-6;
  int use_schur    = 0;
  int nsmin        = 3;
  int nsmax        = 7;
  if (argc > 1) rtol      = atof(argv[1]);
  if (argc > 2) atol_val  = atof(argv[2]);
  if (argc > 3) h0        = atof(argv[3]);
  if (argc > 4) use_schur = atoi(argv[4]);
  if (argc > 5) nsmin     = atoi(argv[5]);
  if (argc > 6) nsmax     = atoi(argv[6]);

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  int n = 8;
  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, nsmin, nsmax);

  /* Initial conditions */
  N_Vector y0 = N_VNew_Serial(n, sunctx);
  sunrealtype* y0d = N_VGetArrayPointer(y0);
  y0d[0] = 0.0;
  y0d[1] = UB / (R2 / R1 + 1.0);   /* 3.0 */
  y0d[2] = y0d[1];                  /* 3.0 */
  y0d[3] = UB;                      /* 6.0 */
  y0d[4] = UB / (R6 / R5 + 1.0);   /* 3.0 */
  y0d[5] = y0d[4];                  /* 3.0 */
  y0d[6] = y0d[3];                  /* 6.0 */
  y0d[7] = 0.0;

  Radau5Init(mem, rhs_transamp, 0.0, y0);

  SUNMatrix Jt = SUNDenseMatrix(n, n, sunctx);
  Radau5SetLinearSolver(mem, Jt, NULL);
  Radau5SetSchurDecomp(mem, use_schur);
  Radau5SetJacFn(mem, jac_transamp);

  /* Mass matrix */
  SUNMatrix Mt = SUNDenseMatrix(n, n, sunctx);
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

  printf("=== Transistor Amplifier (rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n", rtol, atol_val, h0, use_schur);
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
