/* ---------------------------------------------------------------------------
 * RADAU5 — Rocket ascent with multi-phase event detection
 *
 * Ported from SUNDIALS cvRocket_dns.c. A rocket ascends vertically with
 * decreasing mass (fuel burn). The system is:
 *
 *   dH/dt = v                    (height)
 *   dv/dt = a(t,v)               (velocity)
 *
 * where acceleration depends on engine state:
 *   Engine ON:  a = F/(M_r + M_f0 - r*t) - D*v - g
 *   Engine OFF: a = -D*v - g
 *
 * Phase 1 (engine on): Two root functions
 *   g0 = M_f0 - r*t             (fuel exhaustion)
 *   g1 = H - H_cutoff           (height cutoff at 4000 ft)
 *
 * Phase 2 (engine off): One root function
 *   g0 = v                      (maximum height when v=0)
 *
 * This demonstrates:
 *   - Terminal events with solver restart
 *   - Changing the root function set between phases
 *   - User data to switch RHS behavior
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

/* Problem constants */
#define FORCE   2200.0    /* engine thrust (lbf) */
#define MASSR   10.0      /* rocket mass without fuel (slugs) */
#define MASSF0  1.0       /* initial fuel mass (slugs) */
#define BRATE   0.1       /* fuel burn rate (slugs/s) */
#define DRAG    0.3       /* drag coefficient (1/s) */
#define GRAV    32.0      /* gravitational acceleration (ft/s^2) */
#define HCUT    4000.0    /* height cutoff for engine (ft) */

/* User data: engine state */
typedef struct {
  int engine_on;
} RocketData;

/* RHS */
static int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  RocketData* rd = (RocketData*)user_data;
  sunrealtype* yd  = N_VGetArrayPointer(y);
  sunrealtype* ydd = N_VGetArrayPointer(ydot);

  sunrealtype v = yd[1];
  ydd[0] = v;  /* dH/dt = v */

  sunrealtype acc = 0.0;
  if (rd->engine_on) {
    acc = FORCE / (MASSR + MASSF0 - BRATE * t);
  }
  ydd[1] = acc - DRAG * v - GRAV;  /* dv/dt */

  return 0;
}

/* Jacobian */
static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  (void)t; (void)y; (void)fy; (void)user_data;
  (void)tmp1; (void)tmp2; (void)tmp3;

  SM_ELEMENT_D(J, 0, 0) = 0.0;
  SM_ELEMENT_D(J, 0, 1) = 1.0;
  SM_ELEMENT_D(J, 1, 0) = 0.0;
  SM_ELEMENT_D(J, 1, 1) = -DRAG;

  return 0;
}

/* Root function — Phase 1 (engine on): g0 = fuel mass, g1 = H - Hcut */
static int rootfn_phase1(sunrealtype t, N_Vector y, sunrealtype* gout,
                         void* user_data)
{
  sunrealtype* yd = N_VGetArrayPointer(y);
  (void)user_data;

  gout[0] = MASSF0 - BRATE * t;   /* fuel exhaustion */
  gout[1] = yd[0] - HCUT;         /* height cutoff */

  return 0;
}

/* Root function — Phase 2 (engine off): g0 = v (max height) */
static int rootfn_phase2(sunrealtype t, N_Vector y, sunrealtype* gout,
                         void* user_data)
{
  sunrealtype* yd = N_VGetArrayPointer(y);
  (void)t; (void)user_data;

  gout[0] = yd[1];  /* velocity = 0 → max height */

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

  RocketData rd = { .engine_on = 1 };

  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, nsmin, nsmax);

  /* Initial conditions: H=0, v=0 */
  N_Vector y = N_VNew_Serial(2, sunctx);
  sunrealtype* yd = N_VGetArrayPointer(y);
  yd[0] = 0.0;  /* height */
  yd[1] = 0.0;  /* velocity */

  Radau5Init(mem, rhs, 0.0, y);
  SUNMatrix J = SUNDenseMatrix(2, 2, sunctx);
  Radau5SetLinearSolver(mem, J, NULL);
  Radau5SetJacFn(mem, jac);
  Radau5SetUserData(mem, &rd);
  Radau5SStolerances(mem, 1.0e-5, 1.0e-2);
  Radau5SetInitStep(mem, 1.0e-3);

  /* Phase 1: engine on, 2 root functions */
  Radau5RootInit(mem, 2, rootfn_phase1);

  N_Vector yout = N_VClone(y);
  sunrealtype t = 0.0;
  sunrealtype tout = 1.0;
  int iout = 0;
  int nout = 70;

  printf("Rocket ascent with multi-phase event detection (RADAU5)\n");
  printf("  Phase 1: engine ON, roots = fuel exhaustion OR height cutoff\n");
  printf("  Phase 2: engine OFF, root = max height (v=0)\n\n");
  printf("  %12s  %12s  %12s\n", "Time", "Height", "Velocity");
  printf("  %12s  %12s  %12s\n", "--------", "--------", "--------");

  int phase = 1;
  int engine_cutoff_done = 0;
  int max_height_found = 0;

  while (iout < nout) {
    int ret = Radau5Solve(mem, tout, yout, &t);
    sunrealtype* yo = N_VGetArrayPointer(yout);

    if (ret == RADAU5_ROOT_RETURN) {
      if (phase == 1) {
        /* Engine cutoff event */
        int rootsfound[2];
        Radau5GetRootInfo(mem, rootsfound);
        printf("  ** ENGINE CUTOFF at t = %.6f\n", t);
        printf("     H = %.2f ft, v = %.2f ft/s\n", yo[0], yo[1]);
        if (rootsfound[0] != 0) printf("     Cause: fuel exhaustion\n");
        if (rootsfound[1] != 0) printf("     Cause: height cutoff\n");

        /* Switch to phase 2: engine off */
        rd.engine_on = 0;
        phase = 2;
        engine_cutoff_done = 1;

        /* Reinitialize with new root function (1 component: v=0) */
        N_VScale(1.0, yout, y);
        Radau5Init(mem, rhs, t, y);
        Radau5SetLinearSolver(mem, J, NULL);
        Radau5SetJacFn(mem, jac);
        Radau5SetUserData(mem, &rd);
        Radau5SStolerances(mem, 1.0e-5, 1.0e-2);
        Radau5SetInitStep(mem, 0.1);
        Radau5RootInit(mem, 1, rootfn_phase2);
        /* Only detect falling velocity (direction = -1) to avoid t=0 */
        int dir[1] = {-1};
        Radau5SetRootDirection(mem, dir);

      } else {
        /* Max height event — report once, then disable rootfinding */
        printf("  ** MAX HEIGHT at t = %.6f\n", t);
        printf("     H = %.2f ft, v = %.6f ft/s\n", yo[0], yo[1]);
        max_height_found = 1;
        /* Disable rootfinding for the descent phase */
        Radau5RootInit(mem, 0, NULL);
      }
    } else if (ret == RADAU5_SUCCESS) {
      printf("  %12.4f  %12.2f  %12.4f\n", t, yo[0], yo[1]);
      iout++;
      tout += 1.0;
      /* Stop if rocket has landed (H < 0) */
      if (yo[0] < 0.0) {
        printf("  ** LANDED (H < 0)\n");
        break;
      }
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
  printf("\n  Statistics (final segment):\n");
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
