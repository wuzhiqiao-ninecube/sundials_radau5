/* ---------------------------------------------------------------------------
 * radau5_rober.c — Robertson chemical kinetics (very stiff ODE)
 *
 * y1' = -0.04*y1 + 1e4*y2*y3
 * y2' =  0.04*y1 - 1e4*y2*y3 - 3e7*y2^2
 * y3' =  3e7*y2^2
 *
 * t in [0, 1e11],  y0 = [1, 0, 0]
 * Reference: y(1e11) = [2.0833e-8, 8.3334e-13, 1.0]
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

static int rhs(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)t; (void)ud;
  sunrealtype* yv  = N_VGetArrayPointer(y);
  sunrealtype* ydv = N_VGetArrayPointer(yd);
  ydv[0] = -0.04 * yv[0] + 1.0e4 * yv[1] * yv[2];
  ydv[1] =  0.04 * yv[0] - 1.0e4 * yv[1] * yv[2] - 3.0e7 * yv[1] * yv[1];
  ydv[2] =  3.0e7 * yv[1] * yv[1];
  return 0;
}

static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;
  sunrealtype* yv = N_VGetArrayPointer(y);
  SM_ELEMENT_D(J, 0, 0) = -0.04;
  SM_ELEMENT_D(J, 0, 1) =  1.0e4 * yv[2];
  SM_ELEMENT_D(J, 0, 2) =  1.0e4 * yv[1];
  SM_ELEMENT_D(J, 1, 0) =  0.04;
  SM_ELEMENT_D(J, 1, 1) = -1.0e4 * yv[2] - 6.0e7 * yv[1];
  SM_ELEMENT_D(J, 1, 2) = -1.0e4 * yv[1];
  SM_ELEMENT_D(J, 2, 0) =  0.0;
  SM_ELEMENT_D(J, 2, 1) =  6.0e7 * yv[1];
  SM_ELEMENT_D(J, 2, 2) =  0.0;
  return 0;
}

int main(int argc, char* argv[])
{
  sunrealtype rtol = 1e-10;
  sunrealtype atol = 1e-14;
  sunrealtype h0   = 1e-12;
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
  N_Vector y0 = N_VNew_Serial(3, sunctx);
  N_VGetArrayPointer(y0)[0] = 1.0;
  N_VGetArrayPointer(y0)[1] = 0.0;
  N_VGetArrayPointer(y0)[2] = 0.0;

  Radau5Init(mem, rhs, 0.0, y0);
  SUNMatrix Jt = SUNDenseMatrix(3, 3, sunctx);
  Radau5SetLinearSolver(mem, Jt, NULL);
  Radau5SetSchurDecomp(mem, use_schur);
  Radau5SetJacFn(mem, jac);
  Radau5SStolerances(mem, rtol, atol);
  Radau5SetInitStep(mem, h0);

  N_Vector yout = N_VNew_Serial(3, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 1.0e11, yout, &tret);

  sunrealtype yref[3] = {0.2083340149701255e-07,
                          0.8333360770334713e-13,
                          0.9999999791665050e+00};
  sunrealtype* yd = N_VGetArrayPointer(yout);

  printf("=== Robertson (rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n", rtol, atol, h0, use_schur);
  printf("ret=%d, tret=%.6e\n", ret, tret);
  for (int i = 0; i < 3; i++) {
    sunrealtype err = fabs(yd[i] - yref[i]);
    sunrealtype relerr = (fabs(yref[i]) > 0) ? err / fabs(yref[i]) : err;
    printf("y[%d] = %20.14e  ref = %20.14e  rel_err = %.3e\n",
           i, yd[i], yref[i], relerr);
  }

  long int nstep, naccpt, nrejct, nfcn, njac, ndec;
  Radau5GetNumSteps(mem, &nstep);
  Radau5GetNumAccSteps(mem, &naccpt);
  Radau5GetNumRejSteps(mem, &nrejct);
  Radau5GetNumRhsEvals(mem, &nfcn);
  Radau5GetNumJacEvals(mem, &njac);
  Radau5GetNumDecomps(mem, &ndec);
  printf("nstep=%ld naccpt=%ld nrejct=%ld nfcn=%ld njac=%ld ndec=%ld\n",
         nstep, naccpt, nrejct, nfcn, njac, ndec);

  /* Check relative error on y[0] and y[2] (y[1] is tiny) */
  sunrealtype relerr0 = fabs(yd[0] - yref[0]) / fabs(yref[0]);
  sunrealtype relerr2 = fabs(yd[2] - yref[2]) / fabs(yref[2]);
  int pass = (ret == 0 && relerr0 < 1.0e-3 && relerr2 < 1.0e-6);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout); SUNMatDestroy(Jt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
