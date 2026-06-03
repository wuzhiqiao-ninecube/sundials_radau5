/* ---------------------------------------------------------------------------
 * radau5_plate.c — Vibrating plate problem (IVPtestset)
 *
 * 2D biharmonic operator on NX*NY grid, second-order ODE reformulated as
 * first-order system of dimension N = 2*NX*NY = 80.
 *
 * y[k]'       = y[k+NDEMI]                         (displacement -> velocity)
 * y[k+NDEMI]' = -omega*y[k+NDEMI] - FAC*UC + force (plate dynamics)
 *
 * where UC is the 13-point biharmonic stencil and force is a time/space-
 * dependent Gaussian loading on grid rows NACHS1=2 and NACHS2=4.
 *
 * Parameters: NX=8, NY=5, omega=1000, stiffn=100, weight=200
 *             delx=2/(NX+1), fac=stiffn/delx^4
 * t in [0, 7], y(0) = 0
 *
 * Sparse Jacobian (CSC) solved via KLU.
 * Reference: IVPtestset 2.4, plate problem
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include "radau5.h"

/* Grid and physical parameters (matching Fortran COMMON/TRANS/) */
#define MX     8
#define MY     5
#define NDEMI  (MX * MY)   /* 40 */
#define NEQ    (2 * NDEMI) /* 80 */
#define NACHS1 2
#define NACHS2 4

static const sunrealtype OMEGA  = 1000.0;
static const sunrealtype STIFFN = 100.0;
static const sunrealtype WEIGHT = 200.0;

/* Derived constants (computed once in main, stored globally) */
static sunrealtype DELX;
static sunrealtype FAC;

/* ---------------------------------------------------------------------------
 * RHS: translates Fortran FPLATE exactly
 * ---------------------------------------------------------------------------*/
static int rhs_plate(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)ud;
  sunrealtype* yv  = N_VGetArrayPointer(y);
  sunrealtype* f   = N_VGetArrayPointer(yd);

  int nx = MX, ny = MY, nxm1 = MX - 1, nym1 = MY - 1;

  for (int j = 1; j <= ny; j++) {
    for (int i = 1; i <= nx; i++) {
      int k = (i - 1) + nx * (j - 1);

      f[k] = yv[k + NDEMI];

      /* Biharmonic stencil UC (13-point) */
      sunrealtype uc = 16.0 * yv[k];
      if (i > 1) { uc += yv[k]; uc -= 8.0 * yv[k - 1]; }
      if (i < nx) { uc += yv[k]; uc -= 8.0 * yv[k + 1]; }
      if (j > 1) { uc += yv[k]; uc -= 8.0 * yv[k - nx]; }
      if (j < ny) { uc += yv[k]; uc -= 8.0 * yv[k + nx]; }
      if (i > 1 && j > 1)   uc += 2.0 * yv[k - nx - 1];
      if (i < nx && j > 1)  uc += 2.0 * yv[k - nx + 1];
      if (i > 1 && j < ny)  uc += 2.0 * yv[k + nx - 1];
      if (i < nx && j < ny) uc += 2.0 * yv[k + nx + 1];
      if (i > 2)    uc += yv[k - 2];
      if (i < nxm1) uc += yv[k + 2];
      if (j > 2)    uc += yv[k - 2 * nx];
      if (j < nym1) uc += yv[k + 2 * nx];

      /* External force on rows NACHS1 and NACHS2 */
      sunrealtype force = 0.0;
      if (j == NACHS1 || j == NACHS2) {
        sunrealtype xi = i * DELX;
        force = exp(-5.0 * (t - xi - 2.0) * (t - xi - 2.0))
              + exp(-5.0 * (t - xi - 5.0) * (t - xi - 5.0));
      }

      /* Second equation: plate dynamics */
      f[k + NDEMI] = -OMEGA * yv[k + NDEMI] - FAC * uc + force * WEIGHT;
    }
  }
  return 0;
}

/* ---------------------------------------------------------------------------
 * Sparse CSC helpers: operate on raw pointers to avoid repeated macro calls
 * ---------------------------------------------------------------------------*/
static void sparse_set(sunrealtype* data, sunindextype* colptrs,
                       sunindextype* rowinds, sunindextype row,
                       sunindextype col, sunrealtype val)
{
  for (sunindextype p = colptrs[col]; p < colptrs[col + 1]; p++) {
    if (rowinds[p] == row) { data[p] = val; return; }
  }
}

static void sparse_add(sunrealtype* data, sunindextype* colptrs,
                       sunindextype* rowinds, sunindextype row,
                       sunindextype col, sunrealtype val)
{
  for (sunindextype p = colptrs[col]; p < colptrs[col + 1]; p++) {
    if (rowinds[p] == row) { data[p] += val; return; }
  }
}

/* ---------------------------------------------------------------------------
 * Analytic Jacobian: translates Fortran JPLATE
 *
 * Structure (80x80):
 *   Row k (k<NDEMI):       dF[k]/dy[k+NDEMI] = 1  (identity in upper-right)
 *   Row k+NDEMI (k<NDEMI): biharmonic stencil + damping
 * ---------------------------------------------------------------------------*/
static int jac_plate(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                     void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)y; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;

  sunrealtype*  data    = SM_DATA_S(J);
  sunindextype* colptrs = SM_INDEXPTRS_S(J);
  sunindextype* rowinds = SM_INDEXVALS_S(J);

  /* Zero all data */
  sunindextype nnz_total = colptrs[NEQ];
  for (sunindextype p = 0; p < nnz_total; p++) data[p] = 0.0;

  int nx = MX, ny = MY, nxm1 = MX - 1, nym1 = MY - 1;

  for (int j = 1; j <= ny; j++) {
    for (int i = 1; i <= nx; i++) {
      int k = (i - 1) + nx * (j - 1);

      /* Row k: dF[k]/dy[k+NDEMI] = 1 */
      sparse_set(data, colptrs, rowinds, k, k + NDEMI, 1.0);

      /* Row k+NDEMI: damping term */
      sparse_set(data, colptrs, rowinds, k + NDEMI, k + NDEMI, -OMEGA);

      /* Biharmonic stencil w.r.t. displacement columns */
      sparse_add(data, colptrs, rowinds, k + NDEMI, k, -FAC * 16.0);
      if (i > 1) {
        sparse_add(data, colptrs, rowinds, k + NDEMI, k, -FAC);
        sparse_set(data, colptrs, rowinds, k + NDEMI, k - 1, FAC * 8.0);
      }
      if (i < nx) {
        sparse_add(data, colptrs, rowinds, k + NDEMI, k, -FAC);
        sparse_set(data, colptrs, rowinds, k + NDEMI, k + 1, FAC * 8.0);
      }
      if (j > 1) {
        sparse_add(data, colptrs, rowinds, k + NDEMI, k, -FAC);
        sparse_set(data, colptrs, rowinds, k + NDEMI, k - nx, FAC * 8.0);
      }
      if (j < ny) {
        sparse_add(data, colptrs, rowinds, k + NDEMI, k, -FAC);
        sparse_set(data, colptrs, rowinds, k + NDEMI, k + nx, FAC * 8.0);
      }
      if (i > 1 && j > 1)   sparse_set(data, colptrs, rowinds, k + NDEMI, k - nx - 1, -FAC * 2.0);
      if (i < nx && j > 1)  sparse_set(data, colptrs, rowinds, k + NDEMI, k - nx + 1, -FAC * 2.0);
      if (i > 1 && j < ny)  sparse_set(data, colptrs, rowinds, k + NDEMI, k + nx - 1, -FAC * 2.0);
      if (i < nx && j < ny) sparse_set(data, colptrs, rowinds, k + NDEMI, k + nx + 1, -FAC * 2.0);
      if (i > 2)    sparse_set(data, colptrs, rowinds, k + NDEMI, k - 2, -FAC);
      if (i < nxm1) sparse_set(data, colptrs, rowinds, k + NDEMI, k + 2, -FAC);
      if (j > 2)    sparse_set(data, colptrs, rowinds, k + NDEMI, k - 2 * nx, -FAC);
      if (j < nym1) sparse_set(data, colptrs, rowinds, k + NDEMI, k + 2 * nx, -FAC);
    }
  }
  return 0;
}

/* ---------------------------------------------------------------------------
 * Build the CSC sparsity pattern for the 80x80 Jacobian.
 * Returns the number of nonzeros and fills colptrs[NEQ+1] and rowinds[].
 * ---------------------------------------------------------------------------*/
static sunindextype build_sparsity(sunindextype* colptrs, sunindextype* rowinds)
{
  int nx = MX, ny = MY, nxm1 = MX - 1, nym1 = MY - 1;

  /* First pass: count nonzeros per column.
   * IMPORTANT: include diagonal entries for ALL columns — BuildE1 with M=NULL
   * adds fac1 only to existing diagonal positions in the CSC pattern. */
  sunindextype count[NEQ];
  for (int c = 0; c < NEQ; c++) count[c] = 1; /* diagonal entry for every column */

  for (int j = 1; j <= ny; j++) {
    for (int i = 1; i <= nx; i++) {
      int k = (i - 1) + nx * (j - 1);

      /* Row k depends on column k+NDEMI (velocity) */
      count[k + NDEMI]++;

      /* Row k+NDEMI depends on column k+NDEMI: already counted via diagonal init */

      /* Row k+NDEMI depends on displacement columns via biharmonic stencil */
      count[k]++; /* central point (row k+NDEMI in col k, in addition to diagonal) */
      if (i > 1)  count[k - 1]++;
      if (i < nx) count[k + 1]++;
      if (j > 1)  count[k - nx]++;
      if (j < ny) count[k + nx]++;
      if (i > 1 && j > 1)   count[k - nx - 1]++;
      if (i < nx && j > 1)  count[k - nx + 1]++;
      if (i > 1 && j < ny)  count[k + nx - 1]++;
      if (i < nx && j < ny) count[k + nx + 1]++;
      if (i > 2)    count[k - 2]++;
      if (i < nxm1) count[k + 2]++;
      if (j > 2)    count[k - 2 * nx]++;
      if (j < nym1) count[k + 2 * nx]++;
    }
  }

  /* Build column pointers */
  colptrs[0] = 0;
  for (int c = 0; c < NEQ; c++)
    colptrs[c + 1] = colptrs[c] + count[c];
  sunindextype nnz = colptrs[NEQ];

  /* Second pass: fill row indices */
  sunindextype pos[NEQ];
  for (int c = 0; c < NEQ; c++) pos[c] = colptrs[c];

  /* Insert diagonal entries for all columns first */
  for (int c = 0; c < NEQ; c++)
    rowinds[pos[c]++] = c;

  for (int j = 1; j <= ny; j++) {
    for (int i = 1; i <= nx; i++) {
      int k = (i - 1) + nx * (j - 1);

      /* Row k in column k+NDEMI (off-diagonal: k != k+NDEMI) */
      rowinds[pos[k + NDEMI]++] = k;

      /* Row k+NDEMI in displacement columns */
      rowinds[pos[k]++] = k + NDEMI; /* central */
      if (i > 1)  rowinds[pos[k - 1]++] = k + NDEMI;
      if (i < nx) rowinds[pos[k + 1]++] = k + NDEMI;
      if (j > 1)  rowinds[pos[k - nx]++] = k + NDEMI;
      if (j < ny) rowinds[pos[k + nx]++] = k + NDEMI;
      if (i > 1 && j > 1)   rowinds[pos[k - nx - 1]++] = k + NDEMI;
      if (i < nx && j > 1)  rowinds[pos[k - nx + 1]++] = k + NDEMI;
      if (i > 1 && j < ny)  rowinds[pos[k + nx - 1]++] = k + NDEMI;
      if (i < nx && j < ny) rowinds[pos[k + nx + 1]++] = k + NDEMI;
      if (i > 2)    rowinds[pos[k - 2]++] = k + NDEMI;
      if (i < nxm1) rowinds[pos[k + 2]++] = k + NDEMI;
      if (j > 2)    rowinds[pos[k - 2 * nx]++] = k + NDEMI;
      if (j < nym1) rowinds[pos[k + 2 * nx]++] = k + NDEMI;
    }
  }

  /* Sort row indices within each column (insertion sort — small columns) */
  for (int c = 0; c < NEQ; c++) {
    sunindextype start = colptrs[c];
    sunindextype end   = colptrs[c + 1];
    for (sunindextype a = start + 1; a < end; a++) {
      sunindextype tmp = rowinds[a];
      sunindextype b = a - 1;
      while (b >= start && rowinds[b] > tmp) {
        rowinds[b + 1] = rowinds[b];
        b--;
      }
      rowinds[b + 1] = tmp;
    }
  }

  return nnz;
}

/* ---------------------------------------------------------------------------
 * Reference solution at t=7 (from res_exact_pic.txt)
 * ---------------------------------------------------------------------------*/
static const sunrealtype yref[NEQ] = {
  0.000490143813851336, 0.000980081485560611, 0.001462893811482190,
  0.001915822464411935, 0.002285152533727002, 0.002461353376688549,
  0.002254597413097122, 0.001438312591933600, 0.000849025149228402,
  0.001697885005625757, 0.002535239886068847, 0.003323989552181772,
  0.003977902193560667, 0.004320231736082990, 0.004025679955083897,
  0.002643206356123840, 0.000980287627702671, 0.001960162971121222,
  0.002925787622964379, 0.003831644928823870, 0.004570305067454005,
  0.004922706753377098, 0.004509194826194244, 0.002876625183867201,
  0.000849025149228402, 0.001697885005625757, 0.002535239886068847,
  0.003323989552181772, 0.003977902193560667, 0.004320231736082990,
  0.004025679955083897, 0.002643206356123840, 0.000490143813851336,
  0.000980081485560611, 0.001462893811482190, 0.001915822464411935,
  0.002285152533727002, 0.002461353376688549, 0.002254597413097122,
  0.001438312591933600,
 -0.001177590304545409, -0.002409005827992214, -0.003722140831656533,
 -0.005078780056048207, -0.006302661811097914, -0.006973399942926759,
 -0.006394575120415784, -0.003960464551310118, -0.002040148244040460,
 -0.004174829877953482, -0.006456510337516159, -0.008832503276738242,
 -0.011029624807177369, -0.012352389570141255, -0.011524177328690540,
 -0.007253301886026949, -0.002355180609090818, -0.004818011655984428,
 -0.007444281663313065, -0.010157560112096413, -0.012605323622195828,
 -0.013946799885853517, -0.012789150240831569, -0.007920929102620235,
 -0.002040148244040460, -0.004174829877953482, -0.006456510337516159,
 -0.008832503276738242, -0.011029624807177369, -0.012352389570141255,
 -0.011524177328690540, -0.007253301886026949, -0.001177590304545409,
 -0.002409005827992214, -0.003722140831656533, -0.005078780056048207,
 -0.006302661811097914, -0.006973399942926759, -0.006394575120415784,
 -0.003960464551310118
};

/* ---------------------------------------------------------------------------
 * main
 * ---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  sunrealtype rtol = 1.0e-10;
  sunrealtype atol_val = 1.0e-10;
  sunrealtype h0   = 1.0e-2;
  int use_schur    = 0;
  int nsmin        = 3;
  int nsmax        = 13;
  if (argc > 1) rtol      = atof(argv[1]);
  if (argc > 2) atol_val  = atof(argv[2]);
  if (argc > 3) h0        = atof(argv[3]);
  if (argc > 4) use_schur = atoi(argv[4]);
  if (argc > 5) nsmin     = atoi(argv[5]);
  if (argc > 6) nsmax     = atoi(argv[6]);

  /* Compute derived constants */
  sunrealtype denom = MX + 1;
  DELX = 2.0 / denom;
  sunrealtype ush4 = 1.0 / (DELX * DELX * DELX * DELX);
  FAC = STIFFN * ush4;

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  /* Initial conditions: y(0) = 0 */
  N_Vector y0 = N_VNew_Serial(NEQ, sunctx);
  N_VConst(0.0, y0);

  /* Build sparsity pattern */
  sunindextype colptrs[NEQ + 1];
  sunindextype rowinds_buf[NEQ + NDEMI + NDEMI * 13]; /* upper bound on nnz */
  sunindextype nnz = build_sparsity(colptrs, rowinds_buf);

  /* Create sparse Jacobian template */
  SUNMatrix Jt = SUNSparseMatrix(NEQ, NEQ, nnz, CSC_MAT, sunctx);
  sunindextype* cp = SM_INDEXPTRS_S(Jt);
  sunindextype* ri = SM_INDEXVALS_S(Jt);
  for (int c = 0; c <= NEQ; c++) cp[c] = colptrs[c];
  for (sunindextype k = 0; k < nnz; k++) ri[k] = rowinds_buf[k];
  sunrealtype* dat = SM_DATA_S(Jt);
  for (sunindextype k = 0; k < nnz; k++) dat[k] = 0.0;

  /* Solver setup */
  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, nsmin, nsmax);
  Radau5Init(mem, rhs_plate, 0.0, y0);
  Radau5SetLinearSolver(mem, Jt, NULL);
  Radau5SetSchurDecomp(mem, use_schur);
  Radau5SetJacFn(mem, jac_plate);
  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);

  /* Solve to t=7 */
  N_Vector yout = N_VNew_Serial(NEQ, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 7.0, yout, &tret);

  /* Compare with reference */
  sunrealtype* yd = N_VGetArrayPointer(yout);
  sunrealtype maxerr = 0.0;
  for (int i = 0; i < NEQ; i++) {
    sunrealtype err = fabs(yd[i] - yref[i]);
    if (err > maxerr) maxerr = err;
  }

  printf("=== Plate problem (n=%d, rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n",
         NEQ, rtol, atol_val, h0, use_schur);
  printf("ret=%d, tret=%.10e\n", ret, tret);
  printf("max |y - yref| = %.6e\n", maxerr);

  /* Print first few and last few components */
  printf("\nDisplacement (first 8):\n");
  for (int i = 0; i < 8; i++)
    printf("  y[%2d] = %20.14e  ref = %20.14e  err = %.3e\n",
           i, yd[i], yref[i], fabs(yd[i] - yref[i]));
  printf("Velocity (first 8):\n");
  for (int i = NDEMI; i < NDEMI + 8; i++)
    printf("  y[%2d] = %20.14e  ref = %20.14e  err = %.3e\n",
           i, yd[i], yref[i], fabs(yd[i] - yref[i]));

  long int nstep, naccpt, nrejct, nfcn, njac, ndec;
  Radau5GetNumSteps(mem, &nstep);
  Radau5GetNumAccSteps(mem, &naccpt);
  Radau5GetNumRejSteps(mem, &nrejct);
  Radau5GetNumRhsEvals(mem, &nfcn);
  Radau5GetNumJacEvals(mem, &njac);
  Radau5GetNumDecomps(mem, &ndec);
  printf("\nnstep=%ld naccpt=%ld nrejct=%ld nfcn=%ld njac=%ld ndec=%ld\n",
         nstep, naccpt, nrejct, nfcn, njac, ndec);

  /* Pass if max error < 1e-6 (tight tolerance used) */
  int pass = (ret == 0 && maxerr < 1.0e-6);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout); SUNMatDestroy(Jt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
