/* ---------------------------------------------------------------------------
 * radau5_beam.c — Beam problem (stiff ODE, n=80)
 *
 * Elastic beam modeled as N=40 segments. State: th(1..N) angles,
 * th(N+1..2N) angular velocities. RHS involves solving a tridiagonal
 * system at each evaluation.
 *
 * t in [0, 5],  y0 = all zeros
 * Reference solution from IVPtestset (Alphaserver DS20E, RADAU5, tol=1.1e-18)
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

#define BEAM_N     40
#define BEAM_NN    (2 * BEAM_N)
#define BEAM_NSQ   (BEAM_N * BEAM_N)
#define BEAM_NQUAT (BEAM_NSQ * BEAM_NSQ)

static const sunrealtype PI_VAL = 3.14159265358979324;

static int rhs(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)ud;
  const int N = BEAM_N;
  const int NSQ = BEAM_NSQ;
  const long long NQUATR = (long long)BEAM_NQUAT;

  sunrealtype* th = N_VGetArrayPointer(y);
  sunrealtype* df = N_VGetArrayPointer(yd);

  sunrealtype sth[BEAM_N + 1], cth[BEAM_N + 1];
  sunrealtype V[BEAM_N], W[BEAM_N], U[BEAM_N];
  sunrealtype alpha[BEAM_N], beta[BEAM_N];

  int i;
  sunrealtype thdiff, term1, term2, q;

  /* --- sin/cos of angle differences (indices 1..N-1 map to i=1..N-1) --- */
  for (i = 1; i < N; i++) {
    thdiff = th[i] - th[i - 1];
    sth[i + 1] = sin(thdiff);   /* sth[2..N] */
    cth[i + 1] = cos(thdiff);   /* cth[2..N] */
  }

  /* --- Compute V = stiffness + external force --- */
  if (t > PI_VAL) {
    /* t > pi: no external force */
    V[0] = (-3.0 * th[0] + th[1]) * NQUATR;
    for (i = 1; i < N - 1; i++)
      V[i] = (th[i - 1] - 2.0 * th[i] + th[i + 1]) * NQUATR;
    V[N - 1] = (th[N - 2] - th[N - 1]) * NQUATR;
  } else {
    /* t <= pi: external force Fx = -Fy = -1.5*sin^2(t) */
    sunrealtype fabs_val = 1.5 * sin(t) * sin(t);
    sunrealtype Fx = -fabs_val;
    sunrealtype Fy =  fabs_val;
    V[0] = (-3.0 * th[0] + th[1]) * NQUATR
           + NSQ * (Fy * cos(th[0]) - Fx * sin(th[0]));
    for (i = 1; i < N - 1; i++) {
      term1 = (th[i - 1] - 2.0 * th[i] + th[i + 1]) * NQUATR;
      term2 = NSQ * (Fy * cos(th[i]) - Fx * sin(th[i]));
      V[i] = term1 + term2;
    }
    V[N - 1] = (th[N - 2] - th[N - 1]) * NQUATR
               + NSQ * (Fy * cos(th[N - 1]) - Fx * sin(th[N - 1]));
  }

  /* --- Compute W = D*V --- */
  W[0] = sth[2] * V[1];
  for (i = 1; i < N - 1; i++)
    W[i] = -sth[i + 1] * V[i - 1] + sth[i + 2] * V[i + 1];
  W[N - 1] = -sth[N] * V[N - 2];

  /* --- Add centrifugal term: W(i) += th(N+i)^2 --- */
  for (i = 0; i < N; i++)
    W[i] += th[N + i] * th[N + i];

  /* --- Solve tridiagonal system C*W = W (Thomas algorithm) --- */
  alpha[0] = 1.0;
  for (i = 1; i < N; i++) {
    alpha[i] = 2.0;
    beta[i - 1] = -cth[i + 1];
  }
  alpha[N - 1] = 3.0;

  /* Backward elimination */
  for (i = N - 2; i >= 0; i--) {
    q = beta[i] / alpha[i + 1];
    W[i] = W[i] - W[i + 1] * q;
    alpha[i] = alpha[i] - beta[i] * q;
  }
  /* Forward substitution */
  W[0] = W[0] / alpha[0];
  for (i = 1; i < N; i++)
    W[i] = (W[i] - beta[i - 1] * W[i - 1]) / alpha[i];

  /* --- Compute U = C*V + D*W --- */
  U[0] = V[0] - cth[2] * V[1] + sth[2] * W[1];
  for (i = 1; i < N - 1; i++)
    U[i] = 2.0 * V[i] - cth[i + 1] * V[i - 1] - cth[i + 2] * V[i + 1]
           - sth[i + 1] * W[i - 1] + sth[i + 2] * W[i + 1];
  U[N - 1] = 3.0 * V[N - 1] - cth[N] * V[N - 2] - sth[N] * W[N - 2];

  /* --- Output: df(i) = th(N+i), df(N+i) = U(i) --- */
  for (i = 0; i < N; i++) {
    df[i]     = th[N + i];
    df[N + i] = U[i];
  }

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

  int n = BEAM_NN;
  void* mem = Radau5Create(sunctx);
  N_Vector y0 = N_VNew_Serial(n, sunctx);
  sunrealtype* y0d = N_VGetArrayPointer(y0);
  for (int i = 0; i < n; i++) y0d[i] = 0.0;

  Radau5Init(mem, rhs, 0.0, y0);
  SUNMatrix Jt = SUNDenseMatrix(n, n, sunctx);
  Radau5SetLinearSolver(mem, Jt);
  Radau5SetSchurDecomp(mem, use_schur);
  /* Use DQ Jacobian (no analytic Jacobian provided) */
  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);

  N_Vector yout = N_VNew_Serial(n, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 5.0, yout, &tret);

  sunrealtype* yd = N_VGetArrayPointer(yout);

  /* Reference solution at t=5 (from IVPtestset solut) */
  struct { int idx; sunrealtype val; } ref[] = {
    {  0, -0.5792366591285007e-02 },
    {  9, -0.9119134654647947e-01 },
    { 19, -0.1483595472463012e+00 },
    { 39, -0.1767201761075488e+00 },
    { 40,  0.3747362681329794e-01 },
    { 49,  0.5987609702624270e+00 },
    { 59,  0.9868747821728363e+00 },
    { 79,  0.1186724615113034e+01 },
  };
  int nref = (int)(sizeof(ref) / sizeof(ref[0]));

  printf("=== Beam (rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n",
         rtol, atol_val, h0, use_schur);
  printf("ret=%d, tret=%.6e\n", ret, tret);

  sunrealtype maxrelerr = 0.0;
  for (int k = 0; k < nref; k++) {
    int idx = ref[k].idx;
    sunrealtype rval = ref[k].val;
    sunrealtype err = fabs(yd[idx] - rval);
    sunrealtype relerr = (fabs(rval) > 1e-15) ? err / fabs(rval) : err;
    if (relerr > maxrelerr) maxrelerr = relerr;
    printf("y[%2d] = %20.14e  ref = %20.14e  rel_err = %.3e\n",
           idx, yd[idx], rval, relerr);
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

  int pass = (ret == 0 && maxrelerr < 1.0e-2);
  printf("max_rel_err = %.3e\n", maxrelerr);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout); SUNMatDestroy(Jt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
