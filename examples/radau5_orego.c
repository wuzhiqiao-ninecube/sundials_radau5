/* ---------------------------------------------------------------------------
 * radau5_orego.c — Oregonator (stiff ODE, n=3)
 *
 * f(1) = 77.27*(y2 + y1*(1 - 8.375e-6*y1 - y2))
 * f(2) = (y3 - (1+y1)*y2) / 77.27
 * f(3) = 0.161*(y1 - y3)
 *
 * t in [0, 360],  y0 = [1, 2, 3]
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

#define NEQ 3

static int rhs(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)t; (void)ud;
  sunrealtype* v = N_VGetArrayPointer(y);
  sunrealtype* f = N_VGetArrayPointer(yd);
  f[0] = 77.27 * (v[1] + v[0] * (1.0 - 8.375e-6 * v[0] - v[1]));
  f[1] = (v[2] - (1.0 + v[0]) * v[1]) / 77.27;
  f[2] = 0.161 * (v[0] - v[2]);
  return 0;
}

static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;
  sunrealtype* v = N_VGetArrayPointer(y);
  SM_ELEMENT_D(J, 0, 0) = 77.27 * (1.0 - 2.0 * 8.375e-6 * v[0] - v[1]);
  SM_ELEMENT_D(J, 0, 1) = 77.27 * (1.0 - v[0]);
  SM_ELEMENT_D(J, 0, 2) = 0.0;
  SM_ELEMENT_D(J, 1, 0) = -v[1] / 77.27;
  SM_ELEMENT_D(J, 1, 1) = -(1.0 + v[0]) / 77.27;
  SM_ELEMENT_D(J, 1, 2) = 1.0 / 77.27;
  SM_ELEMENT_D(J, 2, 0) = 0.161;
  SM_ELEMENT_D(J, 2, 1) = 0.0;
  SM_ELEMENT_D(J, 2, 2) = -0.161;
  return 0;
}

int main(int argc, char* argv[])
{
  sunrealtype rtol = 1e-6;
  sunrealtype atol = 1e-6;
  sunrealtype h0   = 1e-6;
  int use_schur    = 0;
  int nsmin        = 3;
  int nsmax        = 7;
  if (argc > 1) rtol      = atof(argv[1]);
  if (argc > 2) atol      = atof(argv[2]);
  if (argc > 3) h0        = atof(argv[3]);
  if (argc > 4) use_schur = atoi(argv[4]);
  if (argc > 5) nsmin     = atoi(argv[5]);
  if (argc > 6) nsmax     = atoi(argv[6]);

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, nsmin, nsmax);
  N_Vector y0 = N_VNew_Serial(NEQ, sunctx);
  sunrealtype* y0v = N_VGetArrayPointer(y0);
  y0v[0] = 1.0; y0v[1] = 2.0; y0v[2] = 3.0;

  Radau5Init(mem, rhs, 0.0, y0);
  SUNMatrix Jt = SUNDenseMatrix(NEQ, NEQ, sunctx);
  Radau5SetLinearSolver(mem, Jt, NULL);
  Radau5SetSchurDecomp(mem, use_schur);
  Radau5SetJacFn(mem, jac);
  Radau5SStolerances(mem, rtol, atol);
  Radau5SetInitStep(mem, h0);

  N_Vector yout = N_VNew_Serial(NEQ, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 360.0, yout, &tret);

  sunrealtype yref[NEQ] = {
    0.1000814870318523e+01,
    0.1228178521549917e+04,
    0.1320554942846706e+03
  };
  sunrealtype* yd = N_VGetArrayPointer(yout);

  printf("=== Oregonator (rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n", rtol, atol, h0, use_schur);
  printf("ret=%d, tret=%.6e\n", ret, tret);
  sunrealtype maxrelerr = 0.0;
  for (int i = 0; i < NEQ; i++) {
    sunrealtype err = fabs(yd[i] - yref[i]);
    sunrealtype relerr = (fabs(yref[i]) > 0) ? err / fabs(yref[i]) : err;
    if (relerr > maxrelerr) maxrelerr = relerr;
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

  int pass = (ret == 0 && maxrelerr < 1.0e-3);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout); SUNMatDestroy(Jt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
