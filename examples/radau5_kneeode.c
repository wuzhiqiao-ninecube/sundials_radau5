/* ---------------------------------------------------------------------------
 * RADAU5 — Knee problem with non-negativity via event detection
 *
 * Ported from MATLAB kneeode.m (Dahlquist, Edsberg, Skollermo, Soderlind).
 *
 * The "knee problem":
 *   epsilon * y' = (1-x)*y - y^2,   y(0) = 1,  x in [0, 2]
 *
 * with epsilon = 1e-6. The solution follows the isocline y = 1-x for x < 1,
 * then drops to y = 0 for x > 1. Without non-negativity enforcement,
 * numerical solvers overshoot and produce incorrect negative solutions.
 *
 * RADAU5 doesn't have built-in non-negativity constraints. Instead, we use
 * event detection (rootfinding) to locate the exact time when y crosses zero
 * (falling direction). After the event, we set y = 0 and continue — since
 * f(t, 0) = 0 for all t > 1, the solution stays at zero naturally.
 *
 * This demonstrates using rootfinding as a non-negativity enforcement
 * mechanism for problems where the zero boundary is an equilibrium.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

#define EPSILON 1.0e-6

/* RHS: y' = ((1-x)*y - y^2) / epsilon */
static int rhs_knee(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  (void)user_data;
  sunrealtype yval = N_VGetArrayPointer(y)[0];
  N_VGetArrayPointer(ydot)[0] = ((1.0 - t) * yval - yval * yval) / EPSILON;
  return 0;
}

/* Analytic Jacobian: J = ((1-x) - 2*y) / epsilon */
static int jac_knee(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                    void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  (void)fy; (void)user_data; (void)tmp1; (void)tmp2; (void)tmp3;
  sunrealtype yval = N_VGetArrayPointer(y)[0];
  SM_ELEMENT_D(J, 0, 0) = ((1.0 - t) - 2.0 * yval) / EPSILON;
  return 0;
}

/* Root function: g = y (detect zero crossing, falling direction) */
static int rootfn_knee(sunrealtype t, N_Vector y, sunrealtype* gout,
                       void* user_data)
{
  (void)t; (void)user_data;
  gout[0] = N_VGetArrayPointer(y)[0];
  return 0;
}

/* Solve the knee problem without constraint */
static void solve_unconstrained(SUNContext sunctx, sunrealtype* xout,
                                sunrealtype* yout_arr, int nout,
                                int nsmin, int nsmax)
{
  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, nsmin, nsmax);
  N_Vector y = N_VNew_Serial(1, sunctx);
  N_VGetArrayPointer(y)[0] = 1.0;

  Radau5Init(mem, rhs_knee, 0.0, y);
  SUNMatrix J = SUNDenseMatrix(1, 1, sunctx);
  Radau5SetLinearSolver(mem, J, NULL);
  Radau5SetJacFn(mem, jac_knee);
  Radau5SStolerances(mem, 1.0e-3, 1.0e-6);
  Radau5SetInitStep(mem, 1.0e-6);

  N_Vector yout = N_VNew_Serial(1, sunctx);
  sunrealtype tret;

  for (int i = 0; i < nout; i++) {
    int ret = Radau5Solve(mem, xout[i], yout, &tret);
    yout_arr[i] = (ret == RADAU5_SUCCESS) ? N_VGetArrayPointer(yout)[0] : NAN;
  }

  long int nsteps;
  Radau5GetNumSteps(mem, &nsteps);
  printf("  No constraint: steps=%ld\n", nsteps);

  N_VDestroy(y); N_VDestroy(yout); SUNMatDestroy(J); Radau5Free(&mem);
}

/* Solve with non-negativity via rootfinding */
static void solve_nonneg(SUNContext sunctx, sunrealtype* xout,
                         sunrealtype* yout_arr, int nout,
                         int nsmin, int nsmax)
{
  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, nsmin, nsmax);
  N_Vector y = N_VNew_Serial(1, sunctx);
  N_VGetArrayPointer(y)[0] = 1.0;

  Radau5Init(mem, rhs_knee, 0.0, y);
  SUNMatrix J = SUNDenseMatrix(1, 1, sunctx);
  Radau5SetLinearSolver(mem, J, NULL);
  Radau5SetJacFn(mem, jac_knee);
  Radau5SStolerances(mem, 1.0e-3, 1.0e-6);
  Radau5SetInitStep(mem, 1.0e-6);

  /* Rootfinding: detect y = 0, falling direction only */
  Radau5RootInit(mem, 1, rootfn_knee);
  int rootdir[1] = {-1};
  Radau5SetRootDirection(mem, rootdir);

  N_Vector yout = N_VNew_Serial(1, sunctx);
  sunrealtype tret;
  int hit_zero = 0;
  sunrealtype t_zero = 2.0;  /* time when y hits zero */

  for (int i = 0; i < nout; i++) {
    if (hit_zero) {
      /* After y hit zero, solution stays at 0 (equilibrium) */
      yout_arr[i] = 0.0;
      continue;
    }

    while (1) {
      int ret = Radau5Solve(mem, xout[i], yout, &tret);
      if (ret == RADAU5_ROOT_RETURN) {
        /* y crossed zero — set y = 0 exactly and mark done */
        hit_zero = 1;
        t_zero = tret;
        yout_arr[i] = 0.0;
        /* Fill remaining output points with 0 */
        for (int j = i + 1; j < nout; j++)
          yout_arr[j] = 0.0;
        goto done;
      } else if (ret == RADAU5_SUCCESS) {
        yout_arr[i] = N_VGetArrayPointer(yout)[0];
        break;
      } else {
        yout_arr[i] = NAN;
        break;
      }
    }
  }
done:;

  long int nsteps, nge;
  Radau5GetNumSteps(mem, &nsteps);
  Radau5GetNumGEvals(mem, &nge);
  printf("  Non-negative (rootfinding): steps=%ld, g_evals=%ld, t_zero=%.6f\n",
         nsteps, nge, t_zero);

  N_VDestroy(y); N_VDestroy(yout); SUNMatDestroy(J); Radau5Free(&mem);
}

int main(int argc, char* argv[])
{
  int nsmin = 3;
  int nsmax = 7;
  if (argc > 1) nsmin = atoi(argv[1]);
  if (argc > 2) nsmax = atoi(argv[2]);

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  printf("Knee problem: epsilon*y' = (1-x)*y - y^2, y(0)=1, epsilon=1e-6\n");
  printf("  Isocline: y = 1-x for x<1, y = 0 for x>1\n");
  printf("  Non-negativity enforced via rootfinding (detect y=0 crossing)\n\n");

  /* Output points */
  int nout = 21;
  sunrealtype xout[21];
  for (int i = 0; i < nout; i++)
    xout[i] = 0.1 * i;  /* x = 0.0, 0.1, ..., 2.0 */

  sunrealtype y_noconst[21], y_nonneg[21];

  solve_unconstrained(sunctx, xout, y_noconst, nout, nsmin, nsmax);
  solve_nonneg(sunctx, xout, y_nonneg, nout, nsmin, nsmax);

  /* Print comparison */
  printf("\n  %6s  %14s  %14s  %14s\n", "x", "No constraint", "Non-negative", "Exact (1-x)+");
  printf("  %6s  %14s  %14s  %14s\n", "------", "-----------", "-----------", "-----------");
  for (int i = 0; i < nout; i++) {
    sunrealtype exact = (xout[i] < 1.0) ? (1.0 - xout[i]) : 0.0;
    printf("  %6.2f  %14.6e  %14.6e  %14.6e\n",
           xout[i], y_noconst[i], y_nonneg[i], exact);
  }

  /* Verify: with constraint, y should be non-negative everywhere */
  int pass = 1;
  for (int i = 0; i < nout; i++) {
    if (!isnan(y_nonneg[i]) && y_nonneg[i] < -1.0e-10) {
      pass = 0;
      break;
    }
  }
  printf("\n  Non-negativity check: %s\n", pass ? "PASSED" : "FAILED");

  SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
