/* ---------------------------------------------------------------------------
 * radau5_medakzo.c — Medical Akzo Nobel problem (ODE, n=400, band ml=2 mu=2)
 *
 * From IVPtestset_2.4: reaction-diffusion system on [0,1] with N=200 grid
 * points, yielding 400 ODEs with pentadiagonal Jacobian.
 *
 * t in [0, 20],  y(2j-1)=0, y(2j)=1 for j=1..200
 * Discontinuity at t=5 (phi switches from 2 to 0).
 *
 * Reference solution from IVPtestset (PSIDE on Cray C90, tol=1e-10).
 * We check 8 representative components at t=20.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_band.h>
#include "radau5.h"

#define NEQN 400
#define NHALF 200
#define K_RATE 100.0
#define C_PARAM 4.0

static int rhs_medakzo(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)ud;
  sunrealtype* yv = N_VGetArrayPointer(y);
  sunrealtype* f  = N_VGetArrayPointer(yd);
  sunindextype N = NHALF;
  sunrealtype dzeta = 1.0 / (sunrealtype)N;
  sunrealtype dzeta2 = dzeta * dzeta;
  sunrealtype phi = (t <= 5.0) ? 2.0 : 0.0;

  /* j=1 (first grid point) */
  {
    sunrealtype zeta = dzeta;
    sunrealtype dum = (zeta - 1.0) * (zeta - 1.0) / C_PARAM;
    sunrealtype alpha = 2.0 * (zeta - 1.0) * dum / C_PARAM;
    sunrealtype beta = dum * dum;

    f[0] = (phi - 2.0 * yv[0] + yv[2]) * beta / dzeta2
         + alpha * (yv[2] - phi) / (2.0 * dzeta)
         - K_RATE * yv[0] * yv[1];
    f[1] = -K_RATE * yv[0] * yv[1];
  }

  /* j=2..N-1 (interior grid points) */
  for (sunindextype j = 2; j <= N - 1; j++)
  {
    sunindextype i = 2 * j - 2; /* 0-based index for y(2j-1) */
    sunrealtype zeta = j * dzeta;
    sunrealtype dum = (zeta - 1.0) * (zeta - 1.0) / C_PARAM;
    sunrealtype alpha = 2.0 * (zeta - 1.0) * dum / C_PARAM;
    sunrealtype beta = dum * dum;

    f[i] = (yv[i-2] - 2.0 * yv[i] + yv[i+2]) * beta / dzeta2
         + alpha * (yv[i+2] - yv[i-2]) / (2.0 * dzeta)
         - K_RATE * yv[i] * yv[i+1];
    f[i+1] = -K_RATE * yv[i] * yv[i+1];
  }

  /* j=N (last grid point, Neumann-like BC: y(2N+1) not present) */
  {
    sunindextype i = 2 * N - 2;
    f[i]   = -K_RATE * yv[i] * yv[i+1];
    f[i+1] = -K_RATE * yv[i] * yv[i+1];
  }

  return 0;
}

int main(int argc, char* argv[])
{
  sunrealtype rtol = 1.0e-7;
  sunrealtype atol_val = 1.0e-7;
  sunrealtype h0   = 1.0e-6;
  int use_schur    = 0;
  if (argc > 1) rtol      = atof(argv[1]);
  if (argc > 2) atol_val  = atof(argv[2]);
  if (argc > 3) h0        = atof(argv[3]);
  if (argc > 4) use_schur = atoi(argv[4]);

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  sunindextype n = NEQN;

  /* Initial conditions: y(2j-1)=0, y(2j)=1 */
  N_Vector y0 = N_VNew_Serial(n, sunctx);
  sunrealtype* y0d = N_VGetArrayPointer(y0);
  for (sunindextype j = 0; j < NHALF; j++)
  {
    y0d[2*j]   = 0.0;
    y0d[2*j+1] = 1.0;
  }

  void* mem = Radau5Create(sunctx);
  Radau5Init(mem, rhs_medakzo, 0.0, y0);

  /* Band Jacobian template: ml=2, mu=2 */
  SUNMatrix Jt = SUNBandMatrix(n, 2, 2, sunctx);
  Radau5SetLinearSolver(mem, Jt);
  Radau5SetSchurDecomp(mem, use_schur);

  /* Use DQ Jacobian (analytic Jacobian available but let's test DQ band) */
  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);

  /* Solve in two phases due to discontinuity at t=5 */
  N_Vector yout = N_VNew_Serial(n, sunctx);
  sunrealtype tret;
  int ret;

  /* Phase 1: t=0 to t=5 */
  ret = Radau5Solve(mem, 5.0, yout, &tret);
  if (ret != 0) {
    printf("=== Medical Akzo Nobel ===\n");
    printf("Phase 1 FAILED: ret=%d at tret=%.6e\n", ret, tret);
    goto cleanup;
  }

  /* Phase 2: t=5 to t=20 (phi switches to 0) */
  ret = Radau5Solve(mem, 20.0, yout, &tret);

  /* Reference solution at t=20 (8 representative components) */
  struct { int idx; sunrealtype ref; } checks[] = {
    {  78, 0.2339942217046434e-03 },  /* y(79) */
    {  79, -0.1127916494884468e-141 }, /* y(80) — essentially 0 */
    { 132, 0.3576835958481664e-03 },  /* y(133) */
    { 133, 0.5931668289615909e-108 }, /* y(134) — essentially 0 */
    { 170, 0.3085949832350532e-03 },  /* y(171) */
    { 171, 0.6193143746524996e-046 }, /* y(172) — essentially 0 */
    { 198, 0.1173741304462833e-03 },  /* y(199) */
    { 199, 0.6190822732534586e-05 },  /* y(200) */
  };
  int nchecks = 8;

  sunrealtype* yd = N_VGetArrayPointer(yout);

  printf("=== Medical Akzo Nobel (rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n", rtol, atol_val, h0, use_schur);
  printf("ret=%d, tret=%.6e\n", ret, tret);

  sunrealtype maxrelerr = 0.0;
  int npass = 0;
  for (int k = 0; k < nchecks; k++)
  {
    int idx = checks[k].idx;
    sunrealtype ref = checks[k].ref;
    sunrealtype val = yd[idx];
    sunrealtype err = fabs(val - ref);
    sunrealtype relerr;

    /* For components that are essentially zero, use absolute error */
    if (fabs(ref) < 1.0e-10) {
      relerr = err; /* absolute */
      printf("y[%3d] = %12.6e  ref = %12.6e  abs_err = %.3e %s\n",
             idx+1, val, ref, err, (err < 1.0e-6) ? "OK" : "!");
      if (err < 1.0e-3) npass++;
    } else {
      relerr = err / fabs(ref);
      if (relerr > maxrelerr) maxrelerr = relerr;
      printf("y[%3d] = %12.6e  ref = %12.6e  rel_err = %.3e %s\n",
             idx+1, val, ref, relerr, (relerr < 1.0e-2) ? "OK" : "!");
      if (relerr < 1.0e-2) npass++;
    }
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

  int pass = (ret == 0 && npass == nchecks);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

cleanup:
  N_VDestroy(y0);
  N_VDestroy(yout);
  SUNMatDestroy(Jt);
  Radau5Free(&mem);
  SUNContext_Free(&sunctx);
  return (ret == 0) ? 0 : 1;
}
