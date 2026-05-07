/* ---------------------------------------------------------------------------
 * radau5_caraxis.c — Car Axis problem (DAE index-3, n=10)
 *
 * M*y' = f(t,y) where M = diag(1,1,1,1, m*eps^2/2, ..., m*eps^2/2, 0, 0)
 *
 * 4 positions + 4 velocities + 2 Lagrange multipliers
 * DAE index: nind1=4, nind2=4, nind3=2
 *
 * t in [0, 3],  reference solution from IVPtestset (GAMD, quad precision)
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

#define NEQ 10

/* Physical parameters */
#define MM   10.0    /* mass */
#define EPS  0.01    /* epsilon */
#define LL   1.0     /* L */
#define L0   0.5     /* L0 */
#define RR   0.1     /* r */
#define WW   10.0    /* omega */
#define GG   1.0     /* gravity */

/* Mass effective coefficient: M*eps^2/2 */
#define MEFF (MM * EPS * EPS / 2.0)

/* ---------------------------------------------------------------------------
 * RHS: M*y' = f(t, y)
 * ---------------------------------------------------------------------------*/
static int rhs_caraxis(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)ud;
  sunrealtype* v = N_VGetArrayPointer(y);
  sunrealtype* f = N_VGetArrayPointer(yd);

  sunrealtype xl   = v[0];
  sunrealtype yl   = v[1];
  sunrealtype xr   = v[2];
  sunrealtype yr   = v[3];
  sunrealtype lam1 = v[8];
  sunrealtype lam2 = v[9];

  /* Moving pivot */
  sunrealtype yb = RR * sin(WW * t);
  sunrealtype xb = sqrt(LL * LL - yb * yb);

  /* Spring lengths */
  sunrealtype Ll = sqrt(xl * xl + yl * yl);
  sunrealtype Lr = sqrt((xr - xb) * (xr - xb) + (yr - yb) * (yr - yb));

  /* f(1..4) = velocities */
  f[0] = v[4];
  f[1] = v[5];
  f[2] = v[6];
  f[3] = v[7];

  /* f(5..8) = forces (Fortran stores M*eps^2/2 * y'' on LHS, so RHS has raw forces) */
  f[4] = (L0 - Ll) * xl / Ll + lam1 * xb + 2.0 * lam2 * (xl - xr);
  f[5] = (L0 - Ll) * yl / Ll + lam1 * yb + 2.0 * lam2 * (yl - yr)
         - MM * EPS * EPS * GG / 2.0;
  f[6] = (L0 - Lr) * (xr - xb) / Lr - 2.0 * lam2 * (xl - xr);
  f[7] = (L0 - Lr) * (yr - yb) / Lr - 2.0 * lam2 * (yl - yr)
         - MM * EPS * EPS * GG / 2.0;

  /* f(9..10) = algebraic constraints */
  f[8] = xb * xl + yb * yl;
  f[9] = (xl - xr) * (xl - xr) + (yl - yr) * (yl - yr) - LL * LL;

  return 0;
}

/* ---------------------------------------------------------------------------
 * Mass matrix: M = diag(1,1,1,1, MEFF,MEFF,MEFF,MEFF, 0,0)
 * ---------------------------------------------------------------------------*/
static int mas_caraxis(sunrealtype t, SUNMatrix M, void* ud,
                       N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)ud; (void)t1; (void)t2; (void)t3;
  SUNMatZero(M);
  /* Positions: M[0..3] = 1 */
  SM_ELEMENT_D(M, 0, 0) = 1.0;
  SM_ELEMENT_D(M, 1, 1) = 1.0;
  SM_ELEMENT_D(M, 2, 2) = 1.0;
  SM_ELEMENT_D(M, 3, 3) = 1.0;
  /* Velocities: M[4..7] = M*eps^2/2 = 5e-4 */
  SM_ELEMENT_D(M, 4, 4) = MEFF;
  SM_ELEMENT_D(M, 5, 5) = MEFF;
  SM_ELEMENT_D(M, 6, 6) = MEFF;
  SM_ELEMENT_D(M, 7, 7) = MEFF;
  /* Multipliers: M[8..9] = 0 (algebraic) */
  return 0;
}

int main(int argc, char* argv[])
{
  sunrealtype rtol = 1.0e-6;
  sunrealtype atol_val = 1.0e-6;
  sunrealtype h0   = 1.0e-10;
  int use_schur    = 0;
  if (argc > 1) rtol      = atof(argv[1]);
  if (argc > 2) atol_val  = atof(argv[2]);
  if (argc > 3) h0        = atof(argv[3]);
  if (argc > 4) use_schur = atoi(argv[4]);

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  void* mem = Radau5Create(sunctx);

  /* Initial conditions (0-based, from Fortran init subroutine) */
  N_Vector y0 = N_VNew_Serial(NEQ, sunctx);
  sunrealtype* y0v = N_VGetArrayPointer(y0);
  y0v[0] = 0.0;          /* xl */
  y0v[1] = 0.5;          /* yl = L0 */
  y0v[2] = 1.0;          /* xr = L */
  y0v[3] = 0.5;          /* yr = L0 */
  y0v[4] = -0.5;         /* xla = -L0/L */
  y0v[5] = 0.0;          /* yla */
  y0v[6] = -0.5;         /* xra = -L0/L */
  y0v[7] = 0.0;          /* yra */
  y0v[8] = 0.0;          /* lam1 */
  y0v[9] = 0.0;          /* lam2 */

  Radau5Init(mem, rhs_caraxis, 0.0, y0);

  SUNMatrix Jt = SUNDenseMatrix(NEQ, NEQ, sunctx);
  Radau5SetLinearSolver(mem, Jt);
  Radau5SetSchurDecomp(mem, use_schur);
  /* No Radau5SetJacFn — use DQ Jacobian */

  SUNMatrix Mt = SUNDenseMatrix(NEQ, NEQ, sunctx);
  Radau5SetMassFn(mem, mas_caraxis, Mt);
  Radau5SetDAEIndex(mem, 4, 4, 2);

  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);

  N_Vector yout = N_VNew_Serial(NEQ, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 3.0, yout, &tret);

  /* Reference solution at t=3 (GAMD, quad precision, from IVPtestset) */
  static const sunrealtype yref[NEQ] = {
    0.493455784275402809122e-1,
    0.496989460230171153861e+0,
    0.104174252488542151681e+1,
    0.373911027265361256927e+0,
   -0.770583684040972357970e-1,
    0.744686658723778553466e-2,
    0.175568157537232222276e-1,
    0.770341043779251976443e+0,
   -0.473688659084893324729e-2,
   -0.110468033125734368808e-2
  };
  sunrealtype* yd = N_VGetArrayPointer(yout);
  printf("=== Car Axis (rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n",
         rtol, atol_val, h0, use_schur);
  printf("ret=%d, tret=%.6e\n", ret, tret);

  sunrealtype maxrelerr_pos = 0.0;  /* positions only (index-1) */
  sunrealtype maxrelerr_all = 0.0;
  for (int i = 0; i < NEQ; i++) {
    sunrealtype err = fabs(yd[i] - yref[i]);
    sunrealtype relerr = (fabs(yref[i]) > 1e-15) ? err / fabs(yref[i]) : err;
    if (relerr > maxrelerr_all) maxrelerr_all = relerr;
    if (i < 4 && relerr > maxrelerr_pos) maxrelerr_pos = relerr;
    printf("y[%2d] = %22.14e  ref = %22.14e  rel_err = %.3e\n",
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

  /* Index-3 DAE: positions (ind1) converge as O(tol), velocities (ind2) as
     O(tol^{2/3}), multipliers (ind3) as O(tol^{1/3}).  Check positions < 1e-2. */
  int pass = (ret == 0 && maxrelerr_pos < 1.0e-2);
  printf("max_pos_rel_err=%.3e  max_all_rel_err=%.3e\n",
         maxrelerr_pos, maxrelerr_all);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout);
  SUNMatDestroy(Jt); SUNMatDestroy(Mt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
