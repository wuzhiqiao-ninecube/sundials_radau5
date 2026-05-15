/* ---------------------------------------------------------------------------
 * RADAU5 — Restricted three-body orbit with direction-sensitive events
 *
 * Ported from MATLAB orbitode.m (Shampine & Gordon, p.246).
 * The restricted three-body problem: a body of infinitesimal mass orbits
 * in the gravitational field of two massive bodies (e.g., Earth-Moon).
 *
 * State: y = [x, y_pos, vx, vy] (position + velocity in rotating frame)
 *
 * Parameters:
 *   mu = 1/82.45 (mass ratio)
 *   mustar = 1 - mu
 *
 * ODE system (rotating frame):
 *   dx/dt  = vx
 *   dy/dt  = vy
 *   dvx/dt = 2*vy + x - mustar*(x+mu)/r13 - mu*(x-mustar)/r23
 *   dvy/dt = -2*vx + y - mustar*y/r13 - mu*y/r23
 *
 * where r13 = ((x+mu)^2 + y^2)^1.5, r23 = ((x-mustar)^2 + y^2)^1.5
 *
 * Event function: dDSQ/dt = 2*((x-x0)*vx + (y-y0)*vy)
 *   This is the time derivative of the squared distance from the initial point.
 *   Two events on the SAME function value, distinguished by direction:
 *     g0: dDSQ/dt, direction = +1 → local minimum of distance (return to start)
 *     g1: dDSQ/dt, direction = -1 → local maximum of distance (farthest point)
 *
 *   g0 is terminal (stop when orbit returns to start).
 *   g1 is non-terminal (just record the farthest point).
 *
 * Since RADAU5 rootfinding is always terminal, we handle the non-terminal
 * event (g1) by simply resuming integration after it fires.
 *
 * Initial conditions: y0 = [1.2, 0, 0, -1.04935750983031990726]
 * Integration: t in [0, 7]
 * Tolerances: rtol=1e-5, atol=1e-4
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

#define NEQ 4

/* Problem parameters */
static const sunrealtype mu_param = 1.0 / 82.45;

/* Initial conditions (shared with event function) */
static const sunrealtype x0 = 1.2;
static const sunrealtype y0_pos = 0.0;

/* RHS function */
static int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* yd  = N_VGetArrayPointer(y);
  sunrealtype* ydd = N_VGetArrayPointer(ydot);
  (void)t; (void)user_data;

  sunrealtype mu = mu_param;
  sunrealtype mustar = 1.0 - mu;

  sunrealtype x  = yd[0], yp = yd[1];
  sunrealtype vx = yd[2], vy = yd[3];

  sunrealtype r13 = pow((x + mu) * (x + mu) + yp * yp, 1.5);
  sunrealtype r23 = pow((x - mustar) * (x - mustar) + yp * yp, 1.5);

  ydd[0] = vx;
  ydd[1] = vy;
  ydd[2] = 2.0 * vy + x - mustar * (x + mu) / r13 - mu * (x - mustar) / r23;
  ydd[3] = -2.0 * vx + yp - mustar * yp / r13 - mu * yp / r23;

  return 0;
}

/* Root function: two events on the same value dDSQ/dt, different directions.
 * g0 = dDSQ/dt with direction +1 (local minimum of distance → return to start)
 * g1 = dDSQ/dt with direction -1 (local maximum of distance → farthest point) */
static int rootfn(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data)
{
  sunrealtype* yd = N_VGetArrayPointer(y);
  (void)t; (void)user_data;

  /* dDSQ/dt = 2 * ((x-x0)*vx + (y-y0)*vy) */
  sunrealtype dDSQdt = 2.0 * ((yd[0] - x0) * yd[2] + (yd[1] - y0_pos) * yd[3]);

  gout[0] = dDSQdt;  /* direction +1: local minimum (return to start) */
  gout[1] = dDSQdt;  /* direction -1: local maximum (farthest point) */

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

  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, nsmin, nsmax);

  /* Initial conditions */
  N_Vector y = N_VNew_Serial(NEQ, sunctx);
  sunrealtype* yd = N_VGetArrayPointer(y);
  yd[0] = 1.2;
  yd[1] = 0.0;
  yd[2] = 0.0;
  yd[3] = -1.04935750983031990726;

  Radau5Init(mem, rhs, 0.0, y);

  /* Linear solver (dense 4x4, DQ Jacobian) */
  SUNMatrix J = SUNDenseMatrix(NEQ, NEQ, sunctx);
  Radau5SetLinearSolver(mem, J, NULL);

  /* Tolerances (matching MATLAB orbitode) */
  Radau5SStolerances(mem, 1.0e-5, 1.0e-4);
  Radau5SetInitStep(mem, 1.0e-3);

  /* Rootfinding: 2 root functions, directions +1 and -1 */
  Radau5RootInit(mem, 2, rootfn);
  int rootdir[2] = {+1, -1};
  Radau5SetRootDirection(mem, rootdir);

  /* Integration */
  N_Vector yout = N_VClone(y);
  sunrealtype t = 0.0;
  sunrealtype tout = 7.0;

  printf("Restricted three-body orbit with direction-sensitive events (RADAU5)\n");
  printf("  mu = 1/82.45, y0 = [1.2, 0, 0, -1.0494]\n");
  printf("  Event 0: dDSQ/dt=0, direction=+1 (local min distance, TERMINAL)\n");
  printf("  Event 1: dDSQ/dt=0, direction=-1 (local max distance, non-terminal)\n\n");

  int n_min = 0, n_max = 0;
  int done = 0;

  while (!done) {
    int ret = Radau5Solve(mem, tout, yout, &t);
    sunrealtype* yo = N_VGetArrayPointer(yout);

    if (ret == RADAU5_ROOT_RETURN) {
      int rootsfound[2];
      Radau5GetRootInfo(mem, rootsfound);

      /* Compute distance from initial point */
      sunrealtype dx = yo[0] - x0;
      sunrealtype dy = yo[1] - y0_pos;
      sunrealtype dist = sqrt(dx * dx + dy * dy);

      if (rootsfound[0] != 0) {
        /* Event 0: local minimum of distance (return to start) — TERMINAL */
        n_min++;
        printf("  LOCAL MIN (return) at t = %.10f\n", t);
        printf("    pos = (%.10f, %.10f), dist = %.6e\n", yo[0], yo[1], dist);
        done = 1;  /* Stop — this is the terminal event */
      }
      if (rootsfound[1] != 0) {
        /* Event 1: local maximum of distance (farthest point) — non-terminal */
        n_max++;
        printf("  LOCAL MAX (farthest) at t = %.10f\n", t);
        printf("    pos = (%.10f, %.10f), dist = %.6f\n", yo[0], yo[1], dist);
        /* Continue integration (non-terminal) */
      }
    } else if (ret == RADAU5_SUCCESS) {
      printf("  Reached tout = %.2f without terminal event\n", tout);
      done = 1;
    } else {
      printf("  Solver error: %d at t = %e\n", ret, t);
      done = 1;
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
  printf("    Local minima (returns): %d, Local maxima (farthest): %d\n",
         n_min, n_max);

  /* Reference: the orbit period is approximately T ≈ 6.19 */
  printf("\n  Expected: orbit returns to start near t ≈ 6.19\n");

  /* Cleanup */
  N_VDestroy(y);
  N_VDestroy(yout);
  SUNMatDestroy(J);
  Radau5Free(&mem);
  SUNContext_Free(&sunctx);

  return 0;
}
