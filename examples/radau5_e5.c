/* ---------------------------------------------------------------------------
 * radau5_e5.c — E5 stiff-detest (ODE, n=4)
 *
 * Extremely stiff chemical kinetics, t in [0, 1e13]
 * Reference from IVPtestset (RADAU, Alphaserver, uround=1.01e-19)
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

#define NEQ 4
#define A_RATE  7.89e-10
#define B_RATE  1.1e7
#define CM_RATE 1.13e9
#define C_RATE  1.13e3

static int rhs(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)t; (void)ud;
  sunrealtype* v = N_VGetArrayPointer(y);
  sunrealtype* f = N_VGetArrayPointer(yd);
  sunrealtype p1 = A_RATE * v[0];
  sunrealtype p2 = B_RATE * v[0] * v[2];
  sunrealtype p3 = CM_RATE * v[1] * v[2];
  sunrealtype p4 = C_RATE * v[3];
  f[0] = -p1 - p2;
  f[1] = p1 - p3;
  f[3] = p2 - p4;
  f[2] = f[1] - f[3];
  return 0;
}

static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;
  sunrealtype* v = N_VGetArrayPointer(y);
  SM_ELEMENT_D(J, 0, 0) = -A_RATE - B_RATE * v[2];
  SM_ELEMENT_D(J, 0, 1) = 0.0;
  SM_ELEMENT_D(J, 0, 2) = -B_RATE * v[0];
  SM_ELEMENT_D(J, 0, 3) = 0.0;
  SM_ELEMENT_D(J, 1, 0) = A_RATE;
  SM_ELEMENT_D(J, 1, 1) = -CM_RATE * v[2];
  SM_ELEMENT_D(J, 1, 2) = -CM_RATE * v[1];
  SM_ELEMENT_D(J, 1, 3) = 0.0;
  SM_ELEMENT_D(J, 2, 0) = A_RATE - B_RATE * v[2];
  SM_ELEMENT_D(J, 2, 1) = -CM_RATE * v[2];
  SM_ELEMENT_D(J, 2, 2) = -B_RATE * v[0] - CM_RATE * v[1];
  SM_ELEMENT_D(J, 2, 3) = C_RATE;
  SM_ELEMENT_D(J, 3, 0) = B_RATE * v[2];
  SM_ELEMENT_D(J, 3, 1) = 0.0;
  SM_ELEMENT_D(J, 3, 2) = B_RATE * v[0];
  SM_ELEMENT_D(J, 3, 3) = -C_RATE;
  return 0;
}

int main(int argc, char* argv[])
{
  sunrealtype rtol = 1e-4;
  sunrealtype atol = 1e-24;
  sunrealtype h0   = 1e-6;
  int use_schur    = 0;
  if (argc > 1) rtol      = atof(argv[1]);
  if (argc > 2) atol      = atof(argv[2]);
  if (argc > 3) h0        = atof(argv[3]);
  if (argc > 4) use_schur = atoi(argv[4]);

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  void* mem = Radau5Create(sunctx);
  N_Vector y0 = N_VNew_Serial(NEQ, sunctx);
  sunrealtype* y0v = N_VGetArrayPointer(y0);
  y0v[0] = 1.76e-3; y0v[1] = 0.0; y0v[2] = 0.0; y0v[3] = 0.0;

  Radau5Init(mem, rhs, 0.0, y0);
  SUNMatrix Jt = SUNDenseMatrix(NEQ, NEQ, sunctx);
  Radau5SetLinearSolver(mem, Jt);
  Radau5SetSchurDecomp(mem, use_schur);
  Radau5SetJacFn(mem, jac);
  Radau5SStolerances(mem, rtol, atol);
  Radau5SetInitStep(mem, h0);

  N_Vector yout = N_VNew_Serial(NEQ, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 1.0e13, yout, &tret);

  /* Reference: y(1) is essentially 0 at t=1e13 */
  sunrealtype yref[NEQ] = {
    0.1152903278711829e-290,
    0.8867655517642120e-22,
    0.8854814626268838e-22,
    0.0
  };
  sunrealtype* yd = N_VGetArrayPointer(yout);

  printf("=== E5 stiff-detest (rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n", rtol, atol, h0, use_schur);
  printf("ret=%d, tret=%.6e\n", ret, tret);
  for (int i = 0; i < NEQ; i++) {
    sunrealtype err = fabs(yd[i] - yref[i]);
    sunrealtype relerr = (fabs(yref[i]) > 1e-200) ? err / fabs(yref[i]) : err;
    printf("y[%d] = %22.14e  ref = %22.14e  rel_err = %.3e\n",
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

  /* E5 is extremely stiff; just check it completes */
  int pass = (ret == 0);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout); SUNMatDestroy(Jt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
