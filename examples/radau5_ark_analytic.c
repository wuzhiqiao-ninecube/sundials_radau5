/* ---------------------------------------------------------------------------
 * radau5_ark_analytic.c — Analytical ODE test problem (from ARKODE examples)
 *
 * dy/dt = lambda*y + 1/(1+t^2) - lambda*atan(t)
 *
 * for t in [0, 10], y(0) = 0, lambda = -100.
 * Exact solution: y(t) = atan(t).
 *
 * Solved with RADAU5 (3-stage order-5 implicit Runge-Kutta, Radau IIA).
 * Dense Jacobian, user-supplied.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

static sunrealtype LAMBDA = -1000000.0;

/* RHS: f(t,y) = lambda*y + 1/(1+t^2) - lambda*atan(t) */
static int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void* ud)
{
  (void)ud;
  sunrealtype u = N_VGetArrayPointer(y)[0];
  N_VGetArrayPointer(ydot)[0] =
      LAMBDA * u + 1.0 / (1.0 + t * t) - LAMBDA * atan(t);
  return 0;
}

/* Jacobian: J = df/dy = lambda */
static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)y; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;
  SM_ELEMENT_D(J, 0, 0) = LAMBDA;
  return 0;
}

int main(int argc, char* argv[])
{
  sunrealtype rtol = 1.0e-5;
  sunrealtype atol_val = 1.0e-10;
  sunrealtype h0 = 1.0e-4;
  int use_schur = 0;
  if (argc > 1) rtol      = atof(argv[1]);
  if (argc > 2) atol_val  = atof(argv[2]);
  if (argc > 3) h0        = atof(argv[3]);
  if (argc > 4) use_schur = atoi(argv[4]);

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  int n = 1;
  void* mem = Radau5Create(sunctx);

  /* Initial condition: y(0) = 0 */
  N_Vector y0 = N_VNew_Serial(n, sunctx);
  N_VConst(0.0, y0);

  Radau5Init(mem, rhs, 0.0, y0);

  SUNMatrix Jt = SUNDenseMatrix(n, n, sunctx);
  Radau5SetLinearSolver(mem, Jt);
  Radau5SetSchurDecomp(mem, use_schur);
  Radau5SetJacFn(mem, jac);
  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);

  N_Vector yout = N_VNew_Serial(n, sunctx);
  sunrealtype tret;

  printf("\nAnalytical ODE test problem (RADAU5):\n");
  printf("   lambda = %.1f\n", LAMBDA);
  printf("   rtol = %.1e, atol = %.1e, h0 = %.1e, schur = %d\n\n",
         rtol, atol_val, h0, use_schur);
  printf("        t           u          u_exact        abs_err\n");
  printf("   --------------------------------------------------------\n");

  /* Integrate to t=1,2,...,10, printing at each output point */
  sunrealtype tout;
  int ret = 0;
  for (tout = 1.0; tout <= 10.0; tout += 1.0)
  {
    ret = Radau5Solve(mem, tout, yout, &tret);
    if (ret < 0) {
      printf("   Solver failure at t=%.1f, ret=%d\n", tout, ret);
      break;
    }
    sunrealtype u = N_VGetArrayPointer(yout)[0];
    sunrealtype u_exact = atan(tret);
    sunrealtype err = fabs(u - u_exact);
    printf("  %10.6f  %14.10f  %14.10f  %.3e\n", tret, u, u_exact, err);
  }
  printf("   --------------------------------------------------------\n");

  /* Final error check */
  sunrealtype u_final = N_VGetArrayPointer(yout)[0];
  sunrealtype u_exact = atan(10.0);
  sunrealtype final_err = fabs(u_final - u_exact);
  sunrealtype ewt = 1.0 / (rtol * fabs(u_exact) + atol_val);
  sunrealtype werr = ewt * final_err;

  /* Statistics */
  long int nstep, naccpt, nrejct, nfcn, njac, ndec, nsol;
  Radau5GetNumSteps(mem, &nstep);
  Radau5GetNumAccSteps(mem, &naccpt);
  Radau5GetNumRejSteps(mem, &nrejct);
  Radau5GetNumRhsEvals(mem, &nfcn);
  Radau5GetNumJacEvals(mem, &njac);
  Radau5GetNumDecomps(mem, &ndec);
  Radau5GetNumLinSolves(mem, &nsol);

  printf("\nFinal Solver Statistics:\n");
  printf("   Internal solver steps = %ld (accepted = %ld, rejected = %ld)\n",
         nstep, naccpt, nrejct);
  printf("   Total RHS evals = %ld\n", nfcn);
  printf("   Total Jacobian evals = %ld\n", njac);
  printf("   Total LU decompositions = %ld\n", ndec);
  printf("   Total linear solves = %ld\n", nsol);
  printf("\n   Final solution: y(10) = %.15e\n", u_final);
  printf("   Exact solution: atan(10) = %.15e\n", u_exact);
  printf("   Absolute error: %.3e\n", final_err);
  printf("   Weighted error: %.3e\n\n", werr);

  int pass = (ret >= 0 && werr < 1.5);
  printf("%s\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout);
  SUNMatDestroy(Jt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
