/* ---------------------------------------------------------------------------
 * radau5_vdpol.c — Van der Pol oscillator (eps=1e-6, very stiff)
 *
 * y1' = y2
 * y2' = ((1-y1^2)*y2 - y1) / eps,   eps = 1e-6
 *
 * t in [0, 2],  y0 = [2, 0]
 * Reference: y(2) = [1.7061677321..., -0.8928097010...]
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

#define EPS 1.0e-6

static int rhs(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)t; (void)ud;
  sunrealtype* yv  = N_VGetArrayPointer(y);
  sunrealtype* ydv = N_VGetArrayPointer(yd);
  ydv[0] = yv[1];
  ydv[1] = ((1.0 - yv[0] * yv[0]) * yv[1] - yv[0]) / EPS;
  return 0;
}

static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;
  sunrealtype* yv = N_VGetArrayPointer(y);
  SM_ELEMENT_D(J, 0, 0) = 0.0;
  SM_ELEMENT_D(J, 0, 1) = 1.0;
  SM_ELEMENT_D(J, 1, 0) = (-2.0 * yv[0] * yv[1] - 1.0) / EPS;
  SM_ELEMENT_D(J, 1, 1) = (1.0 - yv[0] * yv[0]) / EPS;
  return 0;
}

int main(int argc, char* argv[])
{
  sunrealtype rtol = 1e-6;
  sunrealtype atol = 1e-6;
  sunrealtype h0   = 1e-6;
  int use_schur    = 0;
  if (argc > 1) rtol      = atof(argv[1]);
  if (argc > 2) atol      = atof(argv[2]);
  if (argc > 3) h0        = atof(argv[3]);
  if (argc > 4) use_schur = atoi(argv[4]);

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  void* mem = Radau5Create(sunctx);
  N_Vector y0 = N_VNew_Serial(2, sunctx);
  N_VGetArrayPointer(y0)[0] = 2.0;
  N_VGetArrayPointer(y0)[1] = 0.0;

  Radau5Init(mem, rhs, 0.0, y0);
  SUNMatrix Jt = SUNDenseMatrix(2, 2, sunctx);
  Radau5SetLinearSolver(mem, Jt);
  Radau5SetSchurDecomp(mem, use_schur);
  Radau5SetJacFn(mem, jac);
  Radau5SStolerances(mem, rtol, atol);
  Radau5SetInitStep(mem, h0);

  N_Vector yout = N_VNew_Serial(2, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 2.0, yout, &tret);

  sunrealtype yref[2] = {0.1706167732170483e+01, -0.8928097010247975e+00};
  sunrealtype* yd = N_VGetArrayPointer(yout);

  printf("=== Van der Pol (rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n", rtol, atol, h0, use_schur);
  printf("ret=%d, tret=%.6e\n", ret, tret);
  printf("y[0] = %20.14e  ref = %20.14e  err = %.3e\n",
         yd[0], yref[0], fabs(yd[0] - yref[0]));
  printf("y[1] = %20.14e  ref = %20.14e  err = %.3e\n",
         yd[1], yref[1], fabs(yd[1] - yref[1]));

  long int nstep, naccpt, nrejct, nfcn, njac, ndec;
  Radau5GetNumSteps(mem, &nstep);
  Radau5GetNumAccSteps(mem, &naccpt);
  Radau5GetNumRejSteps(mem, &nrejct);
  Radau5GetNumRhsEvals(mem, &nfcn);
  Radau5GetNumJacEvals(mem, &njac);
  Radau5GetNumDecomps(mem, &ndec);
  printf("nstep=%ld naccpt=%ld nrejct=%ld nfcn=%ld njac=%ld ndec=%ld\n",
         nstep, naccpt, nrejct, nfcn, njac, ndec);

  int pass = (ret == 0 && fabs(yd[0] - yref[0]) < 1.0e-4
                       && fabs(yd[1] - yref[1]) < 1.0e-4);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout); SUNMatDestroy(Jt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
