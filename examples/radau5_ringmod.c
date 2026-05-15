/* ---------------------------------------------------------------------------
 * radau5_ringmod.c — Ring Modulator problem (IVPtestset)
 *
 * ODE, n=15, dense Jacobian (analytic).
 * t in [0, 1e-3],  y0 = 0 (all components)
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

#define NEQ 15

/* Physical parameters */
#define C    1.6e-8
#define CS   2e-12
#define CP   1e-8
#define R    25e3
#define RP   50.0
#define LH   4.45
#define LS1  2e-3
#define LS2  5e-4
#define LS3  5e-4
#define RG1  36.3
#define RG2  17.3
#define RG3  17.3
#define RI   50.0
#define RC   600.0

#define GAMMA  40.67286402e-9
#define DELTA  17.7493332
#define PI     3.141592653589793238462643383

/* ---------------------------------------------------------------------------
 * RHS
 * ---------------------------------------------------------------------------*/
static int rhs(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)ud;
  sunrealtype* v = N_VGetArrayPointer(y);
  sunrealtype* f = N_VGetArrayPointer(yd);

  sunrealtype uin1 = 0.5 * sin(2000.0 * PI * t);
  sunrealtype uin2 = 2.0 * sin(20000.0 * PI * t);

  sunrealtype ud1 =  v[2] - v[4] - v[6] - uin2;
  sunrealtype ud2 = -v[3] + v[5] - v[6] - uin2;
  sunrealtype ud3 =  v[3] + v[4] + v[6] + uin2;
  sunrealtype ud4 = -v[2] - v[5] + v[6] + uin2;

  /* Overflow protection */
  sunrealtype udmax = ud1;
  if (ud2 > udmax) udmax = ud2;
  if (ud3 > udmax) udmax = ud3;
  if (ud4 > udmax) udmax = ud4;
  if (DELTA * udmax > 300.0) return 1;  /* recoverable: exp overflow */

  sunrealtype qud1 = GAMMA * (exp(DELTA * ud1) - 1.0);
  sunrealtype qud2 = GAMMA * (exp(DELTA * ud2) - 1.0);
  sunrealtype qud3 = GAMMA * (exp(DELTA * ud3) - 1.0);
  sunrealtype qud4 = GAMMA * (exp(DELTA * ud4) - 1.0);

  f[0]  = (v[7] - 0.5*v[9] + 0.5*v[10] + v[13] - v[0]/R) / C;
  f[1]  = (v[8] - 0.5*v[11] + 0.5*v[12] + v[14] - v[1]/R) / C;
  f[2]  = (v[9] - qud1 + qud4) / CS;
  f[3]  = (-v[10] + qud2 - qud3) / CS;
  f[4]  = (v[11] + qud1 - qud3) / CS;
  f[5]  = (-v[12] - qud2 + qud4) / CS;
  f[6]  = (-v[6]/RP + qud1 + qud2 - qud3 - qud4) / CP;
  f[7]  = -v[0] / LH;
  f[8]  = -v[1] / LH;
  f[9]  = (0.5*v[0] - v[2] - RG2*v[9]) / LS2;
  f[10] = (-0.5*v[0] + v[3] - RG3*v[10]) / LS3;
  f[11] = (0.5*v[1] - v[4] - RG2*v[11]) / LS2;
  f[12] = (-0.5*v[1] + v[5] - RG3*v[12]) / LS3;
  f[13] = (-v[0] + uin1 - (RI + RG1)*v[13]) / LS1;
  f[14] = (-v[1] - (RC + RG1)*v[14]) / LS1;

  return 0;
}

/* ---------------------------------------------------------------------------
 * Analytic Jacobian
 * ---------------------------------------------------------------------------*/
static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;
  sunrealtype* v = N_VGetArrayPointer(y);

  /* Zero the entire matrix first */
  for (int i = 0; i < NEQ; i++)
    for (int j = 0; j < NEQ; j++)
      SM_ELEMENT_D(J, i, j) = 0.0;

  sunrealtype uin2 = 2.0 * sin(20000.0 * PI * t);

  sunrealtype ud1 =  v[2] - v[4] - v[6] - uin2;
  sunrealtype ud2 = -v[3] + v[5] - v[6] - uin2;
  sunrealtype ud3 =  v[3] + v[4] + v[6] + uin2;
  sunrealtype ud4 = -v[2] - v[5] + v[6] + uin2;

  sunrealtype qpud1 = GAMMA * DELTA * exp(DELTA * ud1);
  sunrealtype qpud2 = GAMMA * DELTA * exp(DELTA * ud2);
  sunrealtype qpud3 = GAMMA * DELTA * exp(DELTA * ud3);
  sunrealtype qpud4 = GAMMA * DELTA * exp(DELTA * ud4);

  /* Row 0 */
  SM_ELEMENT_D(J, 0, 0)  = -1.0 / (C * R);
  SM_ELEMENT_D(J, 0, 7)  =  1.0 / C;
  SM_ELEMENT_D(J, 0, 9)  = -0.5 / C;
  SM_ELEMENT_D(J, 0, 10) =  0.5 / C;
  SM_ELEMENT_D(J, 0, 13) =  1.0 / C;

  /* Row 1 */
  SM_ELEMENT_D(J, 1, 1)  = -1.0 / (C * R);
  SM_ELEMENT_D(J, 1, 8)  =  1.0 / C;
  SM_ELEMENT_D(J, 1, 11) = -0.5 / C;
  SM_ELEMENT_D(J, 1, 12) =  0.5 / C;
  SM_ELEMENT_D(J, 1, 14) =  1.0 / C;

  /* Row 2 */
  SM_ELEMENT_D(J, 2, 2) = (-qpud1 - qpud4) / CS;
  SM_ELEMENT_D(J, 2, 4) =  qpud1 / CS;
  SM_ELEMENT_D(J, 2, 5) = -qpud4 / CS;
  SM_ELEMENT_D(J, 2, 6) =  (qpud1 + qpud4) / CS;
  SM_ELEMENT_D(J, 2, 9) =  1.0 / CS;

  /* Row 3 */
  SM_ELEMENT_D(J, 3, 3)  = (-qpud2 - qpud3) / CS;
  SM_ELEMENT_D(J, 3, 4)  = -qpud3 / CS;
  SM_ELEMENT_D(J, 3, 5)  =  qpud2 / CS;
  SM_ELEMENT_D(J, 3, 6)  = (-qpud2 - qpud3) / CS;
  SM_ELEMENT_D(J, 3, 10) = -1.0 / CS;

  /* Row 4 */
  SM_ELEMENT_D(J, 4, 2)  =  qpud1 / CS;
  SM_ELEMENT_D(J, 4, 3)  = -qpud3 / CS;
  SM_ELEMENT_D(J, 4, 4)  = (-qpud1 - qpud3) / CS;
  SM_ELEMENT_D(J, 4, 6)  = (-qpud1 - qpud3) / CS;
  SM_ELEMENT_D(J, 4, 11) =  1.0 / CS;

  /* Row 5 */
  SM_ELEMENT_D(J, 5, 2)  = -qpud4 / CS;
  SM_ELEMENT_D(J, 5, 3)  =  qpud2 / CS;
  SM_ELEMENT_D(J, 5, 5)  = (-qpud2 - qpud4) / CS;
  SM_ELEMENT_D(J, 5, 6)  =  (qpud2 + qpud4) / CS;
  SM_ELEMENT_D(J, 5, 12) = -1.0 / CS;

  /* Row 6 */
  SM_ELEMENT_D(J, 6, 2) =  (qpud1 + qpud4) / CP;
  SM_ELEMENT_D(J, 6, 3) = (-qpud2 - qpud3) / CP;
  SM_ELEMENT_D(J, 6, 4) = (-qpud1 - qpud3) / CP;
  SM_ELEMENT_D(J, 6, 5) =  (qpud2 + qpud4) / CP;
  SM_ELEMENT_D(J, 6, 6) = (-qpud1 - qpud2 - qpud3 - qpud4 - 1.0/RP) / CP;

  /* Row 7 */
  SM_ELEMENT_D(J, 7, 0) = -1.0 / LH;

  /* Row 8 */
  SM_ELEMENT_D(J, 8, 1) = -1.0 / LH;

  /* Row 9 */
  SM_ELEMENT_D(J, 9, 0) =  0.5 / LS2;
  SM_ELEMENT_D(J, 9, 2) = -1.0 / LS2;
  SM_ELEMENT_D(J, 9, 9) = -RG2 / LS2;

  /* Row 10 */
  SM_ELEMENT_D(J, 10, 0)  = -0.5 / LS3;
  SM_ELEMENT_D(J, 10, 3)  =  1.0 / LS3;
  SM_ELEMENT_D(J, 10, 10) = -RG3 / LS3;

  /* Row 11 */
  SM_ELEMENT_D(J, 11, 1)  =  0.5 / LS2;
  SM_ELEMENT_D(J, 11, 4)  = -1.0 / LS2;
  SM_ELEMENT_D(J, 11, 11) = -RG2 / LS2;

  /* Row 12 */
  SM_ELEMENT_D(J, 12, 1)  = -0.5 / LS3;
  SM_ELEMENT_D(J, 12, 5)  =  1.0 / LS3;
  SM_ELEMENT_D(J, 12, 12) = -RG3 / LS3;

  /* Row 13 */
  SM_ELEMENT_D(J, 13, 0)  = -1.0 / LS1;
  SM_ELEMENT_D(J, 13, 13) = -(RI + RG1) / LS1;

  /* Row 14 */
  SM_ELEMENT_D(J, 14, 1)  = -1.0 / LS1;
  SM_ELEMENT_D(J, 14, 14) = -(RC + RG1) / LS1;

  return 0;
}

/* ---------------------------------------------------------------------------
 * main
 * ---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  sunrealtype rtol = 1.0e-7;
  sunrealtype atol_val = 1.0e-7;
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

  N_Vector y0 = N_VNew_Serial(NEQ, sunctx);
  N_VConst(0.0, y0);

  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, nsmin, nsmax);
  Radau5Init(mem, rhs, 0.0, y0);

  SUNMatrix Jt = SUNDenseMatrix(NEQ, NEQ, sunctx);
  Radau5SetLinearSolver(mem, Jt, NULL);
  Radau5SetSchurDecomp(mem, use_schur);
  Radau5SetJacFn(mem, jac);
  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);
  Radau5SetMaxNumSteps(mem, 500000);

  N_Vector yout = N_VNew_Serial(NEQ, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 1.0e-3, yout, &tret);

  static const sunrealtype yref[NEQ] = {
    -0.2339057358486745e-01,
    -0.7367485485540825e-02,
     0.2582956709291169e+00,
    -0.4064465721283450e+00,
    -0.4039455665149794e+00,
     0.2607966765422943e+00,
     0.1106761861269975e+00,
     0.2939904342435596e-06,
    -0.2840029933642329e-07,
     0.7267198267264553e-03,
     0.7929487196960840e-03,
    -0.7255283495698965e-03,
    -0.7941401968526521e-03,
     0.7088495416976114e-04,
     0.2390059075236570e-04
  };

  sunrealtype* yd = N_VGetArrayPointer(yout);

  printf("=== Ring Modulator (rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n", rtol, atol_val, h0, use_schur);
  printf("ret=%d  tret=%.6e\n\n", ret, tret);
  printf("%-4s  %-22s  %-22s  %-10s\n", "i", "y[i]", "ref[i]", "rel_err");

  sunrealtype maxrelerr = 0.0;
  for (int i = 0; i < NEQ; i++) {
    sunrealtype err    = fabs(yd[i] - yref[i]);
    sunrealtype relerr = (fabs(yref[i]) > 1e-15) ? err / fabs(yref[i]) : err;
    if (relerr > maxrelerr) maxrelerr = relerr;
    printf("[%2d]  %22.14e  %22.14e  %.3e\n", i, yd[i], yref[i], relerr);
  }

  long int nstep, naccpt, nrejct, nfcn, njac, ndec, nsol;
  Radau5GetNumSteps(mem, &nstep);
  Radau5GetNumAccSteps(mem, &naccpt);
  Radau5GetNumRejSteps(mem, &nrejct);
  Radau5GetNumRhsEvals(mem, &nfcn);
  Radau5GetNumJacEvals(mem, &njac);
  Radau5GetNumDecomps(mem, &ndec);
  Radau5GetNumLinSolves(mem, &nsol);
  printf("\nnstep=%ld  naccpt=%ld  nrejct=%ld  nfcn=%ld  njac=%ld  ndec=%ld  nsol=%ld\n",
         nstep, naccpt, nrejct, nfcn, njac, ndec, nsol);

  int pass = (ret == 0 && maxrelerr < 1.0e-1);
  printf("max_rel_err=%.3e\n", maxrelerr);
  printf("%s\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout); SUNMatDestroy(Jt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
