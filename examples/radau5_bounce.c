/* ---------------------------------------------------------------------------
 * RADAU5 — Bouncing ball example with event detection
 *
 * Demonstrates rootfinding: a ball is dropped from height y0=10 with zero
 * initial velocity. The ODE is:
 *   y1' = y2       (height)
 *   y2' = -9.81    (velocity, gravity)
 *
 * Event function: g = y1 (ball hits ground when height = 0, falling)
 * Direction: -1 (only detect falling crossings)
 *
 * On each bounce, velocity is reversed with coefficient of restitution = 0.9.
 * Integration continues until 10 bounces or t = 10.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

/* RHS: y1' = y2, y2' = -9.81 */
static int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* yd  = N_VGetArrayPointer(y);
  sunrealtype* ydd = N_VGetArrayPointer(ydot);
  (void)t; (void)user_data;

  ydd[0] = yd[1];        /* dy1/dt = y2 */
  ydd[1] = -9.8;        /* dy2/dt = -g */
  return 0;
}

/* Jacobian (constant) */
static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  (void)t; (void)y; (void)fy; (void)user_data;
  (void)tmp1; (void)tmp2; (void)tmp3;

  SM_ELEMENT_D(J, 0, 0) = 0.0;
  SM_ELEMENT_D(J, 0, 1) = 1.0;
  SM_ELEMENT_D(J, 1, 0) = 0.0;
  SM_ELEMENT_D(J, 1, 1) = 0.0;
  return 0;
}

/* Root function: g = y1 (height) */
static int rootfn(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data)
{
  (void)t; (void)user_data;
  gout[0] = N_VGetArrayPointer(y)[0];
  return 0;
}

int main(int argc, char* argv[])
{
  int nsmin = 3;
  int nsmax = 7;
  if (argc > 1) nsmin = atoi(argv[1]);
  if (argc > 2) nsmax = atoi(argv[2]);

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  /* Create solver */
  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, nsmin, nsmax);

  /* Initial conditions: height=10, velocity=0 */
  N_Vector y = N_VNew_Serial(2, sunctx);
  sunrealtype* yd = N_VGetArrayPointer(y);
  yd[0] = 0.0;  /* height */
  yd[1] = 20.0;   /* velocity */

  Radau5Init(mem, rhs, 0.0, y);

  /* Linear solver (dense 2x2) */
  SUNMatrix J = SUNDenseMatrix(2, 2, sunctx);
  Radau5SetLinearSolver(mem, J, NULL);
  Radau5SetJacFn(mem, jac);

  /* Tolerances */
  Radau5SStolerances(mem, 1.0e-10, 1.0e-12);
  Radau5SetInitStep(mem, 1.0e-4);

  /* Rootfinding: detect g = y1 = 0, falling direction only */
  Radau5RootInit(mem, 1, rootfn);
  int rootdir[1] = {-1};
  Radau5SetRootDirection(mem, rootdir);

  /* Integration loop */
  sunrealtype t = 0.0;
  sunrealtype tout = 30.0;
  sunrealtype cor = 0.9;  /* coefficient of restitution */
  int nbounce = 0;
  int max_bounces = 10;

  N_Vector yout = N_VClone(y);

  printf("Bouncing ball with event detection (RADAU5)\n");
  printf("  Initial height = 10.0, g = 9.81, restitution = %.2f\n", cor);
  printf("  %-4s  %12s  %12s  %12s\n", "Bounce", "Time", "Height", "Velocity");
  printf("  %-4s  %12s  %12s  %12s\n", "------", "--------", "--------", "--------");

  while (nbounce < max_bounces) {
    int ret = Radau5Solve(mem, tout, yout, &t);

    if (ret == RADAU5_ROOT_RETURN) {
      nbounce++;
      sunrealtype* ydata = N_VGetArrayPointer(yout);
      printf("  %-4d  %12.8f  %12.4e  %12.6f\n",
             nbounce, t, ydata[0], ydata[1]);

      /* Bounce: reverse velocity with restitution */
      ydata[1] = -cor * ydata[1];

      /* Reset solver at discontinuity (velocity jump) */
      N_VScale(1.0, yout, y);
      Radau5Init(mem, rhs, t, y);
      Radau5SetLinearSolver(mem, J, NULL);
      Radau5SetJacFn(mem, jac);
      Radau5SStolerances(mem, 1.0e-10, 1.0e-12);
      Radau5SetInitStep(mem, 1.0e-4);
      Radau5RootInit(mem, 1, rootfn);
      Radau5SetRootDirection(mem, rootdir);

    } else if (ret == RADAU5_SUCCESS) {
      printf("  Reached tout = %.2f\n", tout);
      break;
    } else {
      printf("  Solver error: %d\n", ret);
      break;
    }
  }

  /* Print statistics */
  long int nsteps, nfcn, njac, nge;
  Radau5GetNumSteps(mem, &nsteps);
  Radau5GetNumRhsEvals(mem, &nfcn);
  Radau5GetNumJacEvals(mem, &njac);
  Radau5GetNumGEvals(mem, &nge);
  printf("\n  Statistics (final segment):\n");
  printf("    Steps: %ld, RHS evals: %ld, Jac evals: %ld, G evals: %ld\n",
         nsteps, nfcn, njac, nge);
  printf("    Total bounces: %d\n", nbounce);

  /* Analytical check for first bounce:
     y = 10 - 0.5*9.81*t^2 = 0 => t = sqrt(20/9.81) ≈ 1.42784... */
  printf("\n  Analytical first bounce time: %.10f\n", sqrt(20.0/9.81));

  /* Cleanup */
  N_VDestroy(y);
  N_VDestroy(yout);
  SUNMatDestroy(J);
  Radau5Free(&mem);
  SUNContext_Free(&sunctx);

  return 0;
}
