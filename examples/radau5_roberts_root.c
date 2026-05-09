/* ---------------------------------------------------------------------------
 * RADAU5 — Robertson chemical kinetics with event detection
 *
 * Ported from SUNDIALS cvRoberts_dns.c. The Robertson problem is a classic
 * stiff chemical kinetics system:
 *
 *   dy1/dt = -0.04*y1 + 1e4*y2*y3
 *   dy2/dt =  0.04*y1 - 1e4*y2*y3 - 3e7*y2^2
 *   dy3/dt =  3e7*y2^2
 *
 * y(0) = [1, 0, 0], t in [0, 4e10]
 *
 * Two root functions (non-terminal events):
 *   g0 = y1 - 1e-4   (species 1 drops to 1e-4)
 *   g1 = y3 - 0.01   (species 3 reaches 0.01)
 *
 * Since RADAU5 rootfinding is terminal (solver stops at each root),
 * we simply resume integration after each root return.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

#define NEQ 3

/* RHS function */
static int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* yd  = N_VGetArrayPointer(y);
  sunrealtype* ydd = N_VGetArrayPointer(ydot);
  (void)t; (void)user_data;

  sunrealtype y1 = yd[0], y2 = yd[1], y3 = yd[2];

  ydd[0] = -0.04 * y1 + 1.0e4 * y2 * y3;
  ydd[2] = 3.0e7 * y2 * y2;
  ydd[1] = -ydd[0] - ydd[2];

  return 0;
}

/* Jacobian */
static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  sunrealtype* yd = N_VGetArrayPointer(y);
  (void)t; (void)fy; (void)user_data; (void)tmp1; (void)tmp2; (void)tmp3;

  sunrealtype y2 = yd[1], y3 = yd[2];

  SM_ELEMENT_D(J, 0, 0) = -0.04;
  SM_ELEMENT_D(J, 0, 1) = 1.0e4 * y3;
  SM_ELEMENT_D(J, 0, 2) = 1.0e4 * y2;

  SM_ELEMENT_D(J, 1, 0) = 0.04;
  SM_ELEMENT_D(J, 1, 1) = -1.0e4 * y3 - 6.0e7 * y2;
  SM_ELEMENT_D(J, 1, 2) = -1.0e4 * y2;

  SM_ELEMENT_D(J, 2, 0) = 0.0;
  SM_ELEMENT_D(J, 2, 1) = 6.0e7 * y2;
  SM_ELEMENT_D(J, 2, 2) = 0.0;

  return 0;
}

/* Root function: g0 = y1 - 1e-4, g1 = y3 - 0.01 */
static int rootfn(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data)
{
  sunrealtype* yd = N_VGetArrayPointer(y);
  (void)t; (void)user_data;

  gout[0] = yd[0] - 1.0e-4;
  gout[1] = yd[2] - 0.01;

  return 0;
}

int main(void)
{
  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  void* mem = Radau5Create(sunctx);

  /* Initial conditions */
  N_Vector y = N_VNew_Serial(NEQ, sunctx);
  sunrealtype* yd = N_VGetArrayPointer(y);
  yd[0] = 1.0;
  yd[1] = 0.0;
  yd[2] = 0.0;

  Radau5Init(mem, rhs, 0.0, y);

  /* Linear solver */
  SUNMatrix J = SUNDenseMatrix(NEQ, NEQ, sunctx);
  Radau5SetLinearSolver(mem, J);
  Radau5SetJacFn(mem, jac);

  /* Tolerances (matching CVODE example) */
  Radau5SStolerances(mem, 1.0e-4, 1.0e-8);
  Radau5SetInitStep(mem, 1.0e-4);

  /* Rootfinding: 2 root functions */
  Radau5RootInit(mem, 2, rootfn);

  /* Integration: output at decades from 0.4 to 4e10 */
  N_Vector yout = N_VClone(y);
  sunrealtype t = 0.0;
  sunrealtype tout = 0.4;
  int iout = 0;
  int nout = 12;

  printf("Robertson chemical kinetics with event detection (RADAU5)\n");
  printf("  Root functions: g0 = y1 - 1e-4, g1 = y3 - 0.01\n\n");

  while (iout < nout) {
    int ret = Radau5Solve(mem, tout, yout, &t);

    if (ret == RADAU5_ROOT_RETURN) {
      int rootsfound[2];
      Radau5GetRootInfo(mem, rootsfound);
      sunrealtype* yo = N_VGetArrayPointer(yout);
      printf("  Root at t = %12.4e   y = %14.6e %14.6e %14.6e\n",
             t, yo[0], yo[1], yo[2]);
      printf("    rootsfound[] = %3d %3d\n", rootsfound[0], rootsfound[1]);
      /* Continue integration (non-terminal behavior) */
    } else if (ret == RADAU5_SUCCESS) {
      sunrealtype* yo = N_VGetArrayPointer(yout);
      printf("  t = %12.4e   y = %14.6e %14.6e %14.6e\n",
             t, yo[0], yo[1], yo[2]);
      iout++;
      tout *= 10.0;
    } else {
      printf("  Solver error: %d at t = %e\n", ret, t);
      break;
    }
  }

  /* Statistics */
  long int nsteps, nfcn, njac, nge;
  Radau5GetNumSteps(mem, &nsteps);
  Radau5GetNumRhsEvals(mem, &nfcn);
  Radau5GetNumJacEvals(mem, &njac);
  Radau5GetNumGEvals(mem, &nge);
  printf("\n  Statistics:\n");
  printf("    Steps: %ld, RHS evals: %ld, Jac evals: %ld, G evals: %ld\n",
         nsteps, nfcn, njac, nge);

  /* Cleanup */
  N_VDestroy(y);
  N_VDestroy(yout);
  SUNMatDestroy(J);
  Radau5Free(&mem);
  SUNContext_Free(&sunctx);

  return 0;
}
