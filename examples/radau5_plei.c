/* ---------------------------------------------------------------------------
 * radau5_plei.c — Pleiades 7-body gravitational problem (IVPtestset)
 *
 * ODE, n=28, dense Jacobian (analytic).
 * t in [0, 3],  y0 = positions and velocities of 7 bodies with masses m_j = j.
 *
 * Reference: Test Set for IVP Solvers
 *   http://www.dm.uniba.it/~testset/
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

#define NEQ 28

static int rhs(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)t; (void)ud;
  sunrealtype* v = N_VGetArrayPointer(y);
  sunrealtype* f = N_VGetArrayPointer(yd);

  /* positions' = velocities */
  for (int i = 0; i < 14; i++) f[i] = v[i + 14];

  /* accelerations */
  for (int i = 0; i < 7; i++) {
    sunrealtype ax = 0.0, ay = 0.0;
    for (int j = 0; j < 7; j++) {
      if (j == i) continue;
      sunrealtype dx  = v[j]     - v[i];
      sunrealtype dy  = v[j + 7] - v[i + 7];
      sunrealtype rij = dx * dx + dy * dy;
      sunrealtype r32 = rij * sqrt(rij);
      sunrealtype mj  = (sunrealtype)(j + 1);
      ax += mj * dx / r32;
      ay += mj * dy / r32;
    }
    f[i + 14] = ax;
    f[i + 21] = ay;
  }
  return 0;
}

static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;
  sunrealtype* v = N_VGetArrayPointer(y);

  /* zero the matrix */
  SUNMatZero(J);

  /* d(pos_i)/d(vel_i) = 1 */
  for (int i = 0; i < 14; i++) SM_ELEMENT_D(J, i, i + 14) = 1.0;

  /* off-diagonal acceleration blocks */
  for (int i = 1; i < 7; i++) {
    for (int j = 0; j < i; j++) {
      sunrealtype mi  = (sunrealtype)(i + 1);
      sunrealtype mj  = (sunrealtype)(j + 1);
      sunrealtype dx  = v[j]     - v[i];
      sunrealtype dy  = v[j + 7] - v[i + 7];
      sunrealtype rij = dx * dx + dy * dy;
      sunrealtype r32 = rij * sqrt(rij);       /* rij^1.5 */
      sunrealtype r52 = rij * rij * sqrt(rij); /* rij^2.5 */

      /* d(ax_i)/d(x_j), d(ax_j)/d(x_i) */
      sunrealtype fjh_xx = (1.0 - 3.0 * dx * dx / rij) / r32;
      SM_ELEMENT_D(J, i + 14, j)     = mj * fjh_xx;
      SM_ELEMENT_D(J, j + 14, i)     = mi * fjh_xx;

      /* d(ay_i)/d(y_j), d(ay_j)/d(y_i) */
      sunrealtype fjh_yy = (1.0 - 3.0 * dy * dy / rij) / r32;
      SM_ELEMENT_D(J, i + 21, j + 7) = mj * fjh_yy;
      SM_ELEMENT_D(J, j + 21, i + 7) = mi * fjh_yy;

      /* cross terms: d(ax)/d(y) and d(ay)/d(x) */
      sunrealtype fjh_xy = -3.0 * dx * dy / r52;
      SM_ELEMENT_D(J, i + 14, j + 7) = mj * fjh_xy;
      SM_ELEMENT_D(J, j + 14, i + 7) = mi * fjh_xy;
      SM_ELEMENT_D(J, i + 21, j)     = mj * fjh_xy;
      SM_ELEMENT_D(J, j + 21, i)     = mi * fjh_xy;
    }
  }

  /* diagonal entries: negative sum of off-diagonal entries in same row */
  for (int i = 0; i < 7; i++) {
    sunrealtype sumxx = 0.0, sumxy = 0.0, sumyy = 0.0;
    for (int j = 0; j < 7; j++) {
      if (j == i) continue;
      sumxx += SM_ELEMENT_D(J, i + 14, j);
      sumxy += SM_ELEMENT_D(J, i + 14, j + 7);
      sumyy += SM_ELEMENT_D(J, i + 21, j + 7);
    }
    SM_ELEMENT_D(J, i + 14, i)     = -sumxx;
    SM_ELEMENT_D(J, i + 14, i + 7) = -sumxy;
    SM_ELEMENT_D(J, i + 21, i)     = -sumxy;
    SM_ELEMENT_D(J, i + 21, i + 7) = -sumyy;
  }

  return 0;
}

int main(int argc, char* argv[])
{
  sunrealtype rtol = 1.0e-6;
  sunrealtype atol_val = 1.0e-6;
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

  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, nsmin, nsmax);
  N_Vector y0 = N_VNew_Serial(NEQ, sunctx);
  sunrealtype* y0d = N_VGetArrayPointer(y0);

  /* x-positions */
  y0d[0]  =  3.0;  y0d[1]  =  3.0;  y0d[2]  = -1.0;  y0d[3]  = -3.0;
  y0d[4]  =  2.0;  y0d[5]  = -2.0;  y0d[6]  =  2.0;
  /* y-positions */
  y0d[7]  =  3.0;  y0d[8]  = -3.0;  y0d[9]  =  2.0;  y0d[10] =  0.0;
  y0d[11] =  0.0;  y0d[12] = -4.0;  y0d[13] =  4.0;
  /* x-velocities */
  y0d[14] =  0.0;  y0d[15] =  0.0;  y0d[16] =  0.0;  y0d[17] =  0.0;
  y0d[18] =  0.0;  y0d[19] =  1.75; y0d[20] = -1.5;
  /* y-velocities */
  y0d[21] =  0.0;  y0d[22] =  0.0;  y0d[23] =  0.0;  y0d[24] = -1.25;
  y0d[25] =  1.0;  y0d[26] =  0.0;  y0d[27] =  0.0;

  Radau5Init(mem, rhs, 0.0, y0);
  SUNMatrix Jt = SUNDenseMatrix(NEQ, NEQ, sunctx);
  Radau5SetLinearSolver(mem, Jt, NULL);
  Radau5SetSchurDecomp(mem, use_schur);
  Radau5SetJacFn(mem, jac);
  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);

  N_Vector yout = N_VNew_Serial(NEQ, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 3.0, yout, &tret);

  sunrealtype yref[NEQ] = {
    /* x-positions */
     0.3706139143970502e+00,  0.3237284092057233e+01, -0.3222559032418324e+01,
     0.6597091455775310e+00,  0.3425581707156584e+00,  0.1562172101400631e+01,
    -0.7003092922212495e+00,
    /* y-positions */
    -0.3943437585517392e+01, -0.3271380973972550e+01,  0.5225081843456543e+01,
    -0.2590612434977470e+01,  0.1198213693392275e+01, -0.2429682344935824e+00,
     0.1091449240428980e+01,
    /* x-velocities */
     0.3417003806314313e+01,  0.1354584501625501e+01, -0.2590065597810775e+01,
     0.2025053734714242e+01, -0.1155815100160448e+01, -0.8072988170223021e+00,
     0.5952396354208710e+00,
    /* y-velocities */
    -0.3741244961234010e+01,  0.3773459685750630e+00,  0.9386858869551073e+00,
     0.3667922227200571e+00, -0.3474046353808490e+00,  0.2344915448180937e+01,
    -0.1947020434263292e+01
  };

  sunrealtype* yd = N_VGetArrayPointer(yout);

  printf("=== Pleiades (rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n", rtol, atol_val, h0, use_schur);
  printf("ret=%d, tret=%.6e\n", ret, tret);

  sunrealtype maxrelerr = 0.0;
  for (int i = 0; i < NEQ; i++) {
    sunrealtype err    = fabs(yd[i] - yref[i]);
    sunrealtype relerr = (fabs(yref[i]) > 1e-15) ? err / fabs(yref[i]) : err;
    if (i < 14 && relerr > maxrelerr) maxrelerr = relerr;
    printf("y[%2d] = %18.12e  ref = %18.12e  rel_err = %.3e\n",
           i, yd[i], yref[i], relerr);
  }

  long int nstep, naccpt, nrejct, nfcn, njac, ndec, nsol;
  Radau5GetNumSteps(mem, &nstep);
  Radau5GetNumAccSteps(mem, &naccpt);
  Radau5GetNumRejSteps(mem, &nrejct);
  Radau5GetNumRhsEvals(mem, &nfcn);
  Radau5GetNumJacEvals(mem, &njac);
  Radau5GetNumDecomps(mem, &ndec);
  Radau5GetNumLinSolves(mem, &nsol);
  printf("nstep=%ld naccpt=%ld nrejct=%ld nfcn=%ld njac=%ld ndec=%ld nsol=%ld\n",
         nstep, naccpt, nrejct, nfcn, njac, ndec, nsol);

  int pass = (ret == 0 && maxrelerr < 1.0e-2);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout); SUNMatDestroy(Jt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
