/* ---------------------------------------------------------------------------
 * radau5_hires.c — HIRES problem (stiff ODE, n=8)
 *
 * t in [0, 321.8122],  y0 = [1,0,0,0,0,0,0,0.0057]
 * Reference solution from IVPtestset (Cray C90, RADAU5, tol=1.1e-18)
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
  sunrealtype* v  = N_VGetArrayPointer(y);
  sunrealtype* f  = N_VGetArrayPointer(yd);
  f[0] = -1.71*v[0] + 0.43*v[1] + 8.32*v[2] + 0.0007;
  f[1] =  1.71*v[0] - 8.75*v[1];
  f[2] = -10.03*v[2] + 0.43*v[3] + 0.035*v[4];
  f[3] =  8.32*v[1] + 1.71*v[2] - 1.12*v[3];
  f[4] = -1.745*v[4] + 0.43*(v[5] + v[6]);
  f[5] = -280.0*v[5]*v[7] + 0.69*v[3] + 1.71*v[4]
         - 0.43*v[5] + 0.69*v[6];
  f[6] =  280.0*v[5]*v[7] - 1.81*v[6];
  f[7] = -f[6];
  return 0;
}

int main(int argc, char* argv[])
{
  sunrealtype rtol = 1.0e-6;
  sunrealtype atol_val = 1.0e-6;
  sunrealtype h0   = 1.0e-6;
  int use_schur    = 0;
  if (argc > 1) rtol      = atof(argv[1]);
  if (argc > 2) atol_val  = atof(argv[2]);
  if (argc > 3) h0        = atof(argv[3]);
  if (argc > 4) use_schur = atoi(argv[4]);

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  int n = 8;
  void* mem = Radau5Create(sunctx);
  N_Vector y0 = N_VNew_Serial(n, sunctx);
  sunrealtype* y0d = N_VGetArrayPointer(y0);
  y0d[0] = 1.0; y0d[1] = 0.0; y0d[2] = 0.0; y0d[3] = 0.0;
  y0d[4] = 0.0; y0d[5] = 0.0; y0d[6] = 0.0; y0d[7] = 0.0057;

  Radau5Init(mem, rhs, 0.0, y0);
  SUNMatrix Jt = SUNDenseMatrix(n, n, sunctx);
  Radau5SetLinearSolver(mem, Jt);
  Radau5SetSchurDecomp(mem, use_schur);
  /* Use DQ Jacobian (no analytic Jacobian provided) */
  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);

  N_Vector yout = N_VNew_Serial(n, sunctx);
  sunrealtype tret;
  sunrealtype tend = 321.8122;
  int ret = Radau5Solve(mem, tend, yout, &tret);

  sunrealtype yref[8] = {
    0.7371312573325668e-3, 0.1442485726316185e-3,
    0.5888729740967575e-4, 0.1175651343283149e-2,
    0.2386356198831331e-2, 0.6238968252742796e-2,
    0.2849998395185769e-2, 0.2850001604814231e-2
  };
  sunrealtype* yd = N_VGetArrayPointer(yout);

  printf("=== HIRES (rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n", rtol, atol_val, h0, use_schur);
  printf("ret=%d, tret=%.6e\n", ret, tret);

  sunrealtype maxrelerr = 0.0;
  for (int i = 0; i < n; i++) {
    sunrealtype err = fabs(yd[i] - yref[i]);
    sunrealtype relerr = (fabs(yref[i]) > 1e-15) ? err / fabs(yref[i]) : err;
    if (relerr > maxrelerr) maxrelerr = relerr;
    printf("y[%d] = %16.10e  ref = %16.10e  rel_err = %.3e\n",
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

  int pass = (ret == 0 && maxrelerr < 1.0e-3);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout); SUNMatDestroy(Jt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
