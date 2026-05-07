/* ---------------------------------------------------------------------------
 * radau5_chemakzo.c — Chemical Akzo Nobel (DAE index-1, n=6)
 *
 * M*y' = f(t,y) where M = diag(1,1,1,1,1,0)
 *
 * 5 ODEs + 1 algebraic constraint: y6 = ks*y1*y4
 *
 * t in [0, 180],  reference solution from IVPtestset (PSIDE, tol=1e-19)
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

#define K1  18.7
#define K2  0.58
#define K3  0.09
#define K4  0.42
#define KBIG 34.4
#define KLA 3.3
#define KS  115.83
#define PO2 0.9
#define HEN 737.0

static int rhs_chemakzo(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)t; (void)ud;
  sunrealtype* v = N_VGetArrayPointer(y);
  sunrealtype* f = N_VGetArrayPointer(yd);

  sunrealtype r1  = K1 * pow(v[0], 4) * sqrt(fabs(v[1]));
  sunrealtype r2  = K2 * v[2] * v[3];
  sunrealtype r3  = (K2 / KBIG) * v[0] * v[4];
  sunrealtype r4  = K3 * v[0] * v[3] * v[3];
  sunrealtype r5  = K4 * v[5] * v[5] * sqrt(fabs(v[1]));
  sunrealtype fin = KLA * (PO2 / HEN - v[1]);

  f[0] = -2.0 * r1 + r2 - r3 - r4;
  f[1] = -0.5 * r1 - r4 - 0.5 * r5 + fin;
  f[2] = r1 - r2 + r3;
  f[3] = -r2 + r3 - 2.0 * r4;
  f[4] = r2 - r3 + r5;
  f[5] = KS * v[0] * v[3] - v[5];  /* algebraic constraint */

  return 0;
}

static int jac_chemakzo(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                        void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;
  sunrealtype* v = N_VGetArrayPointer(y);
  sunindextype i, j;

  /* Zero the Jacobian */
  for (j = 0; j < 6; j++)
    for (i = 0; i < 6; i++)
      SM_ELEMENT_D(J, i, j) = 0.0;

  sunrealtype sq_y2 = sqrt(fabs(v[1]));
  sunrealtype y1_3  = v[0] * v[0] * v[0];

  /* dr1/dy */
  sunrealtype dr1_dy1 = 4.0 * K1 * y1_3 * sq_y2;
  sunrealtype dr1_dy2 = (v[1] > 0.0) ? 0.5 * K1 * pow(v[0], 4) / sq_y2 : 0.0;

  /* dr2/dy */
  sunrealtype dr2_dy3 = K2 * v[3];
  sunrealtype dr2_dy4 = K2 * v[2];

  /* dr3/dy */
  sunrealtype dr3_dy1 = (K2 / KBIG) * v[4];
  sunrealtype dr3_dy5 = (K2 / KBIG) * v[0];

  /* dr4/dy */
  sunrealtype dr4_dy1 = K3 * v[3] * v[3];
  sunrealtype dr4_dy4 = 2.0 * K3 * v[0] * v[3];

  /* dr5/dy */
  sunrealtype dr5_dy2 = (v[1] > 0.0) ? 0.5 * K4 * v[5] * v[5] / sq_y2 : 0.0;
  sunrealtype dr5_dy6 = 2.0 * K4 * v[5] * sq_y2;

  /* f[0] = -2*r1 + r2 - r3 - r4 */
  SM_ELEMENT_D(J, 0, 0) = -2.0 * dr1_dy1 - dr3_dy1 - dr4_dy1;
  SM_ELEMENT_D(J, 0, 1) = -2.0 * dr1_dy2;
  SM_ELEMENT_D(J, 0, 2) = dr2_dy3;
  SM_ELEMENT_D(J, 0, 3) = dr2_dy4 - dr4_dy4;
  SM_ELEMENT_D(J, 0, 4) = -dr3_dy5;

  /* f[1] = -0.5*r1 - r4 - 0.5*r5 + fin */
  SM_ELEMENT_D(J, 1, 0) = -0.5 * dr1_dy1 - dr4_dy1;
  SM_ELEMENT_D(J, 1, 1) = -0.5 * dr1_dy2 - 0.5 * dr5_dy2 - KLA;
  SM_ELEMENT_D(J, 1, 3) = -dr4_dy4;
  SM_ELEMENT_D(J, 1, 5) = -0.5 * dr5_dy6;

  /* f[2] = r1 - r2 + r3 */
  SM_ELEMENT_D(J, 2, 0) = dr1_dy1 + dr3_dy1;
  SM_ELEMENT_D(J, 2, 1) = dr1_dy2;
  SM_ELEMENT_D(J, 2, 2) = -dr2_dy3;
  SM_ELEMENT_D(J, 2, 3) = -dr2_dy4;
  SM_ELEMENT_D(J, 2, 4) = dr3_dy5;

  /* f[3] = -r2 + r3 - 2*r4 */
  SM_ELEMENT_D(J, 3, 0) = dr3_dy1 - 2.0 * dr4_dy1;
  SM_ELEMENT_D(J, 3, 2) = -dr2_dy3;
  SM_ELEMENT_D(J, 3, 3) = -dr2_dy4 - 2.0 * dr4_dy4;
  SM_ELEMENT_D(J, 3, 4) = dr3_dy5;

  /* f[4] = r2 - r3 + r5 */
  SM_ELEMENT_D(J, 4, 0) = -dr3_dy1;
  SM_ELEMENT_D(J, 4, 1) = dr5_dy2;
  SM_ELEMENT_D(J, 4, 2) = dr2_dy3;
  SM_ELEMENT_D(J, 4, 3) = dr2_dy4;
  SM_ELEMENT_D(J, 4, 4) = -dr3_dy5;
  SM_ELEMENT_D(J, 4, 5) = dr5_dy6;

  /* f[5] = ks*y1*y4 - y6 (algebraic) */
  SM_ELEMENT_D(J, 5, 0) = KS * v[3];
  SM_ELEMENT_D(J, 5, 3) = KS * v[0];
  SM_ELEMENT_D(J, 5, 5) = -1.0;

  return 0;
}

static int mas_chemakzo(sunrealtype t, SUNMatrix M, void* ud,
                        N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)ud; (void)t1; (void)t2; (void)t3;
  SUNMatZero(M);
  SM_ELEMENT_D(M, 0, 0) = 1.0;
  SM_ELEMENT_D(M, 1, 1) = 1.0;
  SM_ELEMENT_D(M, 2, 2) = 1.0;
  SM_ELEMENT_D(M, 3, 3) = 1.0;
  SM_ELEMENT_D(M, 4, 4) = 1.0;
  /* M(5,5) = 0 — algebraic equation */
  return 0;
}

int main(int argc, char* argv[])
{
  sunrealtype rtol = 1.0e-7;
  sunrealtype atol_val = 1.0e-7;
  sunrealtype h0   = 1.0e-10;
  int use_schur    = 0;
  if (argc > 1) rtol      = atof(argv[1]);
  if (argc > 2) atol_val  = atof(argv[2]);
  if (argc > 3) h0        = atof(argv[3]);
  if (argc > 4) use_schur = atoi(argv[4]);

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  int n = 6;
  void* mem = Radau5Create(sunctx);

  /* Initial conditions */
  N_Vector y0 = N_VNew_Serial(n, sunctx);
  sunrealtype* y0d = N_VGetArrayPointer(y0);
  y0d[0] = 0.444;
  y0d[1] = 0.00123;
  y0d[2] = 0.0;
  y0d[3] = 0.007;
  y0d[4] = 0.0;
  y0d[5] = KS * y0d[0] * y0d[3];  /* consistent IC */

  Radau5Init(mem, rhs_chemakzo, 0.0, y0);

  SUNMatrix Jt = SUNDenseMatrix(n, n, sunctx);
  Radau5SetLinearSolver(mem, Jt);
  Radau5SetSchurDecomp(mem, use_schur);
  Radau5SetJacFn(mem, jac_chemakzo);

  /* Mass matrix */
  SUNMatrix Mt = SUNDenseMatrix(n, n, sunctx);
  Radau5SetMassFn(mem, mas_chemakzo, Mt);

  /* DAE index: 5 index-1 vars + 1 index-1 algebraic (all counted as nind1) */
  Radau5SetDAEIndex(mem, n, 0, 0);

  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);

  N_Vector yout = N_VNew_Serial(n, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 180.0, yout, &tret);

  /* Reference solution at t=180 */
  sunrealtype yref[6] = {
    0.1150794920661702e+00,
    0.1203831471567715e-02,
    0.1611562887407974e+00,
    0.3656156421249283e-03,
    0.1708010885264404e-01,
    0.4873531310307455e-02
  };
  sunrealtype* yd = N_VGetArrayPointer(yout);

  printf("=== Chemical Akzo Nobel (rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n", rtol, atol_val, h0, use_schur);
  printf("ret=%d, tret=%.6e\n", ret, tret);

  sunrealtype maxrelerr = 0.0;
  for (int i = 0; i < n; i++)
  {
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

  /* Check: algebraic constraint y6 = ks*y1*y4 should hold */
  sunrealtype alg_err = fabs(yd[5] - KS * yd[0] * yd[3]);
  printf("algebraic constraint |y6 - ks*y1*y4| = %.3e\n", alg_err);

  int pass = (ret == 0 && maxrelerr < 1.0e-3 && alg_err < 1.0e-6);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout);
  SUNMatDestroy(Jt); SUNMatDestroy(Mt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
