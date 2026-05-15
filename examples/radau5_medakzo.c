/* ---------------------------------------------------------------------------
 * radau5_medakzo.c — Medical Akzo Nobel problem (ODE, n=400, sparse Jacobian)
 *
 * From IVPtestset_2.4: reaction-diffusion system on [0,1] with N=200 grid
 * points, yielding 400 ODEs with pentadiagonal Jacobian (ml=2, mu=2).
 * Now using sparse CSC matrix with KLU and analytic Jacobian from Fortran.
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
#include <sunmatrix/sunmatrix_sparse.h>
#include "radau5.h"

#define NEQN 400
#define NHALF 200
#define K_RATE 100.0
#define C_PARAM 4.0

/* NNZ for the sparse Jacobian (pentadiagonal, computed from structure) */
#define NNZ 1197

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

/* ---------------------------------------------------------------------------
 * Helper: set a value in a CSC sparse matrix by (row, col) linear search
 * ---------------------------------------------------------------------------*/
static void sparse_set(SUNMatrix A, sunindextype row, sunindextype col,
                       sunrealtype val)
{
  sunindextype* colptrs = SM_INDEXPTRS_S(A);
  sunindextype* rowinds = SM_INDEXVALS_S(A);
  sunrealtype*  data    = SM_DATA_S(A);
  for (sunindextype k = colptrs[col]; k < colptrs[col + 1]; k++) {
    if (rowinds[k] == row) { data[k] = val; return; }
  }
}

/* ---------------------------------------------------------------------------
 * Build the CSC sparsity pattern for the medakzo Jacobian.
 * The pattern is pentadiagonal (ml=mu=2) but sparse — we only store
 * structurally nonzero entries from the analytic Jacobian.
 * ---------------------------------------------------------------------------*/
static void build_sparsity_pattern(SUNMatrix J)
{
  sunindextype* colptrs = SM_INDEXPTRS_S(J);
  sunindextype* rowinds = SM_INDEXVALS_S(J);
  sunindextype N = NHALF;
  sunindextype idx = 0;

  /* Col 0: rows {0, 1, 2} */
  colptrs[0] = 0;
  rowinds[idx++] = 0;
  rowinds[idx++] = 1;
  rowinds[idx++] = 2;

  /* Col 1: rows {0, 1} */
  colptrs[1] = idx;
  rowinds[idx++] = 0;
  rowinds[idx++] = 1;

  /* Interior odd cols: col = 2j-2 (0-based) for j=2..N-1 */
  for (sunindextype j = 2; j <= N - 1; j++)
  {
    sunindextype c = 2 * j - 2;
    colptrs[c] = idx;
    rowinds[idx++] = c - 2;  /* from previous j's dfdy(1,...) or boundary */
    rowinds[idx++] = c;      /* diagonal */
    rowinds[idx++] = c + 1;  /* from even dfdy(4,...) */
    if (j < N - 1) {
      rowinds[idx++] = c + 2;  /* from next j's dfdy(5,...) */
    }

    /* Even col: col = 2j-1 (0-based), rows {c_even-1, c_even} */
    sunindextype ce = 2 * j - 1;
    colptrs[ce] = idx;
    rowinds[idx++] = ce - 1;  /* dfdy(2,...) */
    rowinds[idx++] = ce;      /* diagonal */
  }

  /* Col 398 (= 2N-2): rows {396, 398, 399} */
  sunindextype c398 = 2 * N - 2;
  colptrs[c398] = idx;
  rowinds[idx++] = c398 - 2;  /* 396 */
  rowinds[idx++] = c398;      /* 398 */
  rowinds[idx++] = c398 + 1;  /* 399 */

  /* Col 399 (= 2N-1): rows {398, 399} */
  sunindextype c399 = 2 * N - 1;
  colptrs[c399] = idx;
  rowinds[idx++] = c399 - 1;  /* 398 */
  rowinds[idx++] = c399;      /* 399 */

  colptrs[NEQN] = idx;
}

/* ---------------------------------------------------------------------------
 * Analytic Jacobian — translated from Fortran medakzo.f jeval
 * ---------------------------------------------------------------------------*/
static int jac_medakzo(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                       void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;
  sunrealtype* yv = N_VGetArrayPointer(y);
  sunindextype N = NHALF;
  sunrealtype dzeta = 1.0 / (sunrealtype)N;
  sunrealtype dzeta2 = dzeta * dzeta;

  sunrealtype* data = SM_DATA_S(J);
  sunindextype* colptrs = SM_INDEXPTRS_S(J);
  sunindextype nnz_used = colptrs[NEQN];
  for (sunindextype k = 0; k < nnz_used; k++) data[k] = 0.0;

  /* j=1 boundary */
  {
    sunrealtype dum = (dzeta - 1.0) * (dzeta - 1.0) / C_PARAM;
    sunrealtype alpha = 2.0 * (dzeta - 1.0) * dum / C_PARAM;
    sunrealtype beta = dum * dum;
    sunrealtype bz = beta / dzeta2;

    /* Col 0: diagonal and below */
    sparse_set(J, 0, 0, -2.0 * bz - K_RATE * yv[1]);
    sparse_set(J, 1, 0, -K_RATE * yv[1]);

    /* Col 1 */
    sparse_set(J, 0, 1, -K_RATE * yv[0]);
    sparse_set(J, 1, 1, -K_RATE * yv[0]);

    /* Col 2 gets entry from boundary: dfdy(1,3) → (0, 2) */
    sparse_set(J, 0, 2, bz + alpha / (2.0 * dzeta));
  }

  /* j=2..N-1 interior */
  for (sunindextype j = 2; j <= N - 1; j++)
  {
    sunindextype i = 2 * j - 2;  /* 0-based odd index */
    sunrealtype zeta = j * dzeta;
    sunrealtype dum = (zeta - 1.0) * (zeta - 1.0) / C_PARAM;
    sunrealtype alpha = 2.0 * (zeta - 1.0) * dum / C_PARAM;
    sunrealtype beta = dum * dum;
    sunrealtype bz = beta / dzeta2;

    /* dfdy(5, i_odd-2) → (i, i-2): contributes to col i-2, row i */
    sparse_set(J, i, i - 2, bz - alpha / (2.0 * dzeta));

    /* dfdy(3, i_odd) → (i, i): diagonal of odd col */
    sparse_set(J, i, i, -2.0 * bz - K_RATE * yv[i + 1]);

    /* dfdy(1, i_odd+2) → (i, i+2): contributes to col i+2, row i */
    sparse_set(J, i, i + 2, bz + alpha / (2.0 * dzeta));

    /* dfdy(2, i_odd+1) → (i, i+1): contributes to col i+1, row i */
    sparse_set(J, i, i + 1, -K_RATE * yv[i]);

    /* dfdy(4, i_even-1) → (i+1, i) */
    sparse_set(J, i + 1, i, -K_RATE * yv[i + 1]);

    /* dfdy(3, i_even) → (i+1, i+1) */
    sparse_set(J, i + 1, i + 1, -K_RATE * yv[i]);
  }

  /* j=N boundary */
  {
    sunindextype i = 2 * N - 2;
    sparse_set(J, i, i, -K_RATE * yv[i + 1]);
    sparse_set(J, i, i + 1, -K_RATE * yv[i]);
    sparse_set(J, i + 1, i, -K_RATE * yv[i + 1]);
    sparse_set(J, i + 1, i + 1, -K_RATE * yv[i]);
  }

  return 0;
}

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
  Radau5SetOrderLimits(mem, nsmin, nsmax);
  Radau5Init(mem, rhs_medakzo, 0.0, y0);

  /* Sparse Jacobian template: CSC, 400x400, NNZ entries */
  SUNMatrix Jt = SUNSparseMatrix(n, n, NNZ, CSC_MAT, sunctx);
  build_sparsity_pattern(Jt);

  Radau5SetLinearSolver(mem, Jt, NULL);
  Radau5SetSchurDecomp(mem, use_schur);
  Radau5SetJacFn(mem, jac_medakzo);

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

  printf("=== Medical Akzo Nobel (rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n",
         rtol, atol_val, h0, use_schur);
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
