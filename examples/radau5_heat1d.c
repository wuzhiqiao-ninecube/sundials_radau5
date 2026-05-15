/* ---------------------------------------------------------------------------
 * radau5_heat1d.c — 1D heat equation with band Jacobian
 *
 * u_t = u_xx on [0,1], u(0)=u(1)=0, u(x,0)=sin(pi*x)
 * Exact solution: u(x,t) = exp(-pi^2*t)*sin(pi*x)
 *
 * Discretize with n interior points: dx = 1/(n+1)
 * y_i' = (y_{i-1} - 2*y_i + y_{i+1}) / dx^2,  y_0 = y_{n+1} = 0
 *
 * Jacobian is tridiagonal: mu=1, ml=1
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_band.h>
#include "radau5.h"

#define N_INTERIOR 100
#define PI 3.14159265358979323846

typedef struct { sunrealtype dx2_inv; } HeatData;

static int rhs_heat(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)t;
  HeatData* data = (HeatData*)ud;
  sunrealtype dx2_inv = data->dx2_inv;
  sunindextype n = N_VGetLength(y);
  sunrealtype* yv  = N_VGetArrayPointer(y);
  sunrealtype* ydv = N_VGetArrayPointer(yd);

  /* Interior points with Dirichlet BCs y_0 = y_{n+1} = 0 */
  ydv[0] = dx2_inv * (0.0 - 2.0 * yv[0] + yv[1]);
  for (sunindextype i = 1; i < n - 1; i++)
    ydv[i] = dx2_inv * (yv[i-1] - 2.0 * yv[i] + yv[i+1]);
  ydv[n-1] = dx2_inv * (yv[n-2] - 2.0 * yv[n-1] + 0.0);

  return 0;
}

int main(int argc, char* argv[])
{
  sunrealtype rtol = 1.0e-10;
  sunrealtype atol_val = 1.0e-12;
  sunrealtype h0   = 1.0e-6;
  int use_schur    = 0;
  int nsmin        = 3;
  int nsmax        = 7;
  if (argc > 1) rtol      = atof(argv[1]);
  if (argc > 2) atol_val  = atof(argv[2]);
  if (argc > 3) h0        = atof(argv[3]);
  if (argc > 4) use_schur = atoi(argv[4]);
  if (argc > 5) nsmin     = atoi(argv[5]);
  if (argc > 6) nsmax     = atoi(argv[6]);

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  sunindextype n = N_INTERIOR;
  sunrealtype dx = 1.0 / (n + 1);
  HeatData hdata;
  hdata.dx2_inv = 1.0 / (dx * dx);

  /* Initial condition: u(x,0) = sin(pi*x) */
  N_Vector y0 = N_VNew_Serial(n, sunctx);
  sunrealtype* y0d = N_VGetArrayPointer(y0);
  for (sunindextype i = 0; i < n; i++)
    y0d[i] = sin(PI * (i + 1) * dx);

  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, nsmin, nsmax);
  Radau5Init(mem, rhs_heat, 0.0, y0);
  Radau5SetUserData(mem, &hdata);

  /* Band Jacobian template: tridiagonal (mu=1, ml=1) */
  SUNMatrix Jt = SUNBandMatrix(n, 1, 1, sunctx);
  Radau5SetLinearSolver(mem, Jt, NULL);
  Radau5SetSchurDecomp(mem, use_schur);

  /* Use DQ Jacobian (no analytic Jacobian) */
  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);

  /* Solve to t=0.1 */
  sunrealtype tend = 0.1;
  N_Vector yout = N_VNew_Serial(n, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, tend, yout, &tret);

  /* Exact solution: u(x,t) = exp(-pi^2*t)*sin(pi*x) */
  sunrealtype decay = exp(-PI * PI * tend);
  sunrealtype* yd = N_VGetArrayPointer(yout);
  sunrealtype maxerr = 0.0;
  for (sunindextype i = 0; i < n; i++) {
    sunrealtype exact = decay * sin(PI * (i + 1) * dx);
    sunrealtype err = fabs(yd[i] - exact);
    if (err > maxerr) maxerr = err;
  }

  /* Also compute spatial discretization error for reference */
  sunrealtype disc_err = PI * PI * dx * dx / 12.0 * decay; /* O(dx^2) */

  printf("=== 1D Heat Equation (n=%d, rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n", (int)n, rtol, atol_val, h0, use_schur);
  printf("ret=%d, tret=%.10e\n", ret, tret);
  printf("max |y - exact_discrete| = %.3e\n", maxerr);
  printf("spatial discretization error ~ %.3e\n", disc_err);

  long int nstep, naccpt, nrejct, nfcn, njac, ndec;
  Radau5GetNumSteps(mem, &nstep);
  Radau5GetNumAccSteps(mem, &naccpt);
  Radau5GetNumRejSteps(mem, &nrejct);
  Radau5GetNumRhsEvals(mem, &nfcn);
  Radau5GetNumJacEvals(mem, &njac);
  Radau5GetNumDecomps(mem, &ndec);
  printf("nstep=%ld naccpt=%ld nrejct=%ld nfcn=%ld njac=%ld ndec=%ld\n",
         nstep, naccpt, nrejct, nfcn, njac, ndec);

  /* The error is dominated by spatial discretization O(dx^2) ~ 3e-5,
     not time integration. Check that total error is within 2x of disc error. */
  int pass = (ret == 0 && maxerr < 2.0 * disc_err + 1.0e-8);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout); SUNMatDestroy(Jt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
