/* ---------------------------------------------------------------------------
 * RADAU5 — Kepler two-body orbit with half-orbit counting via events
 *
 * Ported from SUNDIALS ark_kepler.c. A body orbits in a central force field.
 * The root function g = q2 (y-coordinate) fires at each half-orbit crossing
 * of the x-axis.
 *
 * State: y = [q1, q2, p1, p2] (position + momentum, 4 components)
 *
 * Parameters: eccentricity e = 0.6
 *
 * ODE system (Hamiltonian):
 *   dq1/dt = p1
 *   dq2/dt = p2
 *   dp1/dt = -q1 / (q1^2 + q2^2)^(3/2)
 *   dp2/dt = -q2 / (q1^2 + q2^2)^(3/2)
 *
 * Initial conditions (perihelion):
 *   q1(0) = 1 - e = 0.4
 *   q2(0) = 0
 *   p1(0) = 0
 *   p2(0) = sqrt((1+e)/(1-e)) = 2.0
 *
 * Event: g = q2 (y-position crosses zero → half-orbit)
 *   No direction filter — counts both upward and downward crossings.
 *   Non-terminal: integration continues to tf.
 *
 * The orbital period is T = 2*pi for this normalization.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

#define NEQ 4
#define ECC 0.6  /* eccentricity */

/* RHS function */
static int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  sunrealtype* yd  = N_VGetArrayPointer(y);
  sunrealtype* ydd = N_VGetArrayPointer(ydot);
  (void)t; (void)user_data;

  sunrealtype q1 = yd[0], q2 = yd[1];
  sunrealtype r3 = pow(q1 * q1 + q2 * q2, 1.5);

  ydd[0] = yd[2];       /* dq1/dt = p1 */
  ydd[1] = yd[3];       /* dq2/dt = p2 */
  ydd[2] = -q1 / r3;   /* dp1/dt */
  ydd[3] = -q2 / r3;   /* dp2/dt */

  return 0;
}

/* Root function: g = q2 (y-position) */
static int rootfn(sunrealtype t, N_Vector y, sunrealtype* gout, void* user_data)
{
  (void)t; (void)user_data;
  gout[0] = N_VGetArrayPointer(y)[1];  /* q2 */
  return 0;
}

int main(void)
{
  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  void* mem = Radau5Create(sunctx);

  /* Initial conditions at perihelion */
  N_Vector y = N_VNew_Serial(NEQ, sunctx);
  sunrealtype* yd = N_VGetArrayPointer(y);
  yd[0] = 1.0 - ECC;                    /* q1 = 0.4 */
  yd[1] = 0.0;                           /* q2 = 0 */
  yd[2] = 0.0;                           /* p1 = 0 */
  yd[3] = sqrt((1.0 + ECC) / (1.0 - ECC));  /* p2 = 2.0 */

  Radau5Init(mem, rhs, 0.0, y);

  /* Linear solver (dense 4x4, DQ Jacobian) */
  SUNMatrix J = SUNDenseMatrix(NEQ, NEQ, sunctx);
  Radau5SetLinearSolver(mem, J);

  /* Tolerances */
  Radau5SStolerances(mem, 1.0e-8, 1.0e-10);
  Radau5SetInitStep(mem, 1.0e-3);

  /* Rootfinding: g = q2, no direction filter (both crossings) */
  Radau5RootInit(mem, 1, rootfn);

  /* Integrate for 5 full orbits (T = 2*pi ≈ 6.2832) */
  N_Vector yout = N_VClone(y);
  sunrealtype t = 0.0;
  sunrealtype tout =  32 * M_PI; /* ~31.416 */
  int n_crossings = 0;

  printf("Kepler two-body orbit with half-orbit event counting (RADAU5)\n");
  printf("  Eccentricity e = %.1f, period T = 2*pi = %.6f\n", ECC, 2.0*M_PI);
  printf("  Event: g = q2 (x-axis crossing), no direction filter\n");
  printf("  Integrating for 31 orbits (t_final = %.4f)\n\n", tout);
  printf("  %4s  %12s  %12s  %12s  %12s\n",
         "Half", "Time", "q1", "q2", "Energy err");
  printf("  %4s  %12s  %12s  %12s  %12s\n",
         "----", "--------", "--------", "--------", "----------");

  /* Initial energy: H = 0.5*(p1^2+p2^2) - 1/sqrt(q1^2+q2^2) */
  sunrealtype H0 = 0.5 * (yd[2]*yd[2] + yd[3]*yd[3])
                 - 1.0 / sqrt(yd[0]*yd[0] + yd[1]*yd[1]);

  while (1) {
    int ret = Radau5Solve(mem, tout, yout, &t);
    sunrealtype* yo = N_VGetArrayPointer(yout);

    if (ret == RADAU5_ROOT_RETURN) {
      n_crossings++;
      /* Compute energy error */
      sunrealtype H = 0.5 * (yo[2]*yo[2] + yo[3]*yo[3])
                    - 1.0 / sqrt(yo[0]*yo[0] + yo[1]*yo[1]);
      printf("  %4d  %12.8f  %12.8f  %12.4e  %12.4e\n",
             n_crossings, t, yo[0], yo[1], H - H0);
    } else if (ret == RADAU5_SUCCESS) {
      break;
    } else {
      printf("  Solver error: %d at t = %e\n", ret, t);
      break;
    }
  }

  /* Final state */
  sunrealtype* yo = N_VGetArrayPointer(yout);
  sunrealtype Hf = 0.5 * (yo[2]*yo[2] + yo[3]*yo[3])
                 - 1.0 / sqrt(yo[0]*yo[0] + yo[1]*yo[1]);

  /* Statistics */
  long int nsteps, nfcn, njac, nge;
  Radau5GetNumSteps(mem, &nsteps);
  Radau5GetNumRhsEvals(mem, &nfcn);
  Radau5GetNumJacEvals(mem, &njac);
  Radau5GetNumGEvals(mem, &nge);
  printf("\n  Statistics:\n");
  printf("    Steps: %ld, RHS evals: %ld, Jac evals: %ld, G evals: %ld\n",
         nsteps, nfcn, njac, nge);
  printf("    Half-orbit crossings: %d (expected: 31 for 16 orbits)\n", n_crossings);
  printf("    Final energy error: %.6e\n", Hf - H0);
  printf("    Period error per orbit: %.6e\n",
         n_crossings > 0 ? (t / (n_crossings / 2.0) - 2.0 * M_PI) : 0.0);

  /* Cleanup */
  N_VDestroy(y);
  N_VDestroy(yout);
  SUNMatDestroy(J);
  Radau5Free(&mem);
  SUNContext_Free(&sunctx);

  return 0;
}
