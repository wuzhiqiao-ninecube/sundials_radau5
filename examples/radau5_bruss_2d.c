/* ---------------------------------------------------------------------------
 * radau5_bruss_2d.c — 2D Brusselator with diffusion (IVPtestset)
 *
 * Reaction-diffusion on [0,1]^2 with periodic BCs, ns*ns grid:
 *   u_ij' = 1 + u_ij^2*v_ij - 4.4*u_ij + alpha*ns^2*(Laplacian u)
 *   v_ij' = 3.4*u_ij - u_ij^2*v_ij + alpha*ns^2*(Laplacian v)
 *
 * Plus inhomogeneous forcing: +5 on u where (x-0.3)^2+(y-0.6)^2 <= 0.01
 *   (only active for t >= 1.1)
 *
 * Parameters: ns (default 128, input arg), alpha=0.1, tend=11.5
 * System dimension: N = 2*ns^2
 * Ordering: u components [0..ns^2-1], v components [ns^2..2*ns^2-1]
 *
 * Jacobian: sparse (CSC) via KLU — 5-point Laplacian stencil + reaction.
 * Uses DQ Jacobian via sparsity pattern (column grouping).
 *
 * Reference: IVPtestset 2.4, bruss_2d problem (NS=128)
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include "radau5.h"

static const sunrealtype ALPHA = 0.1;

/* Problem dimensions (set in main from command-line arg) */
static int g_ns;
static int g_nssq;
static int g_neq;

/* ---------------------------------------------------------------------------
 * RHS: translates Fortran FBRUS exactly (periodic BCs, 5-point Laplacian)
 * ---------------------------------------------------------------------------*/
static int rhs_bruss2d(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)ud;
  sunrealtype* yv = N_VGetArrayPointer(y);
  sunrealtype* f  = N_VGetArrayPointer(yd);

  int ns = g_ns;
  int nssq = g_nssq;
  sunrealtype fac = ALPHA * nssq;
  sunrealtype radsq = 0.01;  /* 0.1^2 */
  sunrealtype bet = (t >= 1.1) ? 5.0 : 0.0;
  sunrealtype ans = (sunrealtype)ns;

  for (int idx = 0; idx < nssq; idx++) {
    int ileft  = (idx % ns == 0) ? idx + ns - 1 : idx - 1;
    int iright = ((idx + 1) % ns == 0) ? idx - ns + 1 : idx + 1;
    int ilow   = (idx < ns) ? idx + ns * (ns - 1) : idx - ns;
    int iup    = (idx >= ns * (ns - 1)) ? idx - ns * (ns - 1) : idx + ns;

    sunrealtype uij = yv[idx];
    sunrealtype vij = yv[idx + nssq];

    sunrealtype ulap = yv[ileft] + yv[iright] + yv[ilow] + yv[iup] - 4.0 * uij;
    sunrealtype vlap = yv[ileft + nssq] + yv[iright + nssq]
                     + yv[ilow + nssq] + yv[iup + nssq] - 4.0 * vij;

    f[idx] = 1.0 + uij * uij * vij - 4.4 * uij + fac * ulap;
    f[idx + nssq] = 3.4 * uij - uij * uij * vij + fac * vlap;

    /* Inhomogeneous forcing */
    int iy = idx / ns;
    int ix = idx - iy * ns;
    sunrealtype yy = (iy + 1) / ans;
    sunrealtype xx = (ix + 1) / ans;
    if ((xx - 0.3) * (xx - 0.3) + (yy - 0.6) * (yy - 0.6) <= radsq)
      f[idx] += bet;
  }

  return 0;
}

/* ---------------------------------------------------------------------------
 * Build the CSC sparsity pattern for the neq x neq Jacobian.
 * Each equation has 6 nonzeros. Total nnz = 2*nssq*6 = 12*nssq.
 * ---------------------------------------------------------------------------*/
static sunindextype build_sparsity(int ns, int nssq, int neq,
                                   sunindextype* colptrs, sunindextype* rowinds)
{
  /* Count nonzeros per column */
  sunindextype* count = (sunindextype*)calloc(neq, sizeof(sunindextype));

  for (int idx = 0; idx < nssq; idx++) {
    int ileft  = (idx % ns == 0) ? idx + ns - 1 : idx - 1;
    int iright = ((idx + 1) % ns == 0) ? idx - ns + 1 : idx + 1;
    int ilow   = (idx < ns) ? idx + ns * (ns - 1) : idx - ns;
    int iup    = (idx >= ns * (ns - 1)) ? idx - ns * (ns - 1) : idx + ns;

    count[idx]++;          count[idx + nssq]++;
    count[ileft]++;        count[iright]++;
    count[ilow]++;         count[iup]++;

    count[idx]++;          count[idx + nssq]++;
    count[ileft + nssq]++; count[iright + nssq]++;
    count[ilow + nssq]++;  count[iup + nssq]++;
  }

  /* Build column pointers */
  colptrs[0] = 0;
  for (int c = 0; c < neq; c++)
    colptrs[c + 1] = colptrs[c] + count[c];
  sunindextype nnz = colptrs[neq];

  /* Fill row indices */
  sunindextype* pos = (sunindextype*)malloc(neq * sizeof(sunindextype));
  for (int c = 0; c < neq; c++) pos[c] = colptrs[c];

  for (int idx = 0; idx < nssq; idx++) {
    int ileft  = (idx % ns == 0) ? idx + ns - 1 : idx - 1;
    int iright = ((idx + 1) % ns == 0) ? idx - ns + 1 : idx + 1;
    int ilow   = (idx < ns) ? idx + ns * (ns - 1) : idx - ns;
    int iup    = (idx >= ns * (ns - 1)) ? idx - ns * (ns - 1) : idx + ns;

    /* Column idx (u[idx]): rows idx and idx+nssq */
    rowinds[pos[idx]++] = idx;
    rowinds[pos[idx]++] = idx + nssq;

    /* Column idx+nssq (v[idx]): rows idx and idx+nssq */
    rowinds[pos[idx + nssq]++] = idx;
    rowinds[pos[idx + nssq]++] = idx + nssq;

    /* u-diffusion neighbors: row idx */
    rowinds[pos[ileft]++] = idx;
    rowinds[pos[iright]++] = idx;
    rowinds[pos[ilow]++] = idx;
    rowinds[pos[iup]++] = idx;

    /* v-diffusion neighbors: row idx+nssq */
    rowinds[pos[ileft + nssq]++] = idx + nssq;
    rowinds[pos[iright + nssq]++] = idx + nssq;
    rowinds[pos[ilow + nssq]++] = idx + nssq;
    rowinds[pos[iup + nssq]++] = idx + nssq;
  }

  /* Sort row indices within each column */
  for (int c = 0; c < neq; c++) {
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

  free(count);
  free(pos);
  return nnz;
}

/* ---------------------------------------------------------------------------
 * Reference solution for NS=128 (from res_exact_pic.txt)
 * The Fortran driver outputs y(1..N step 267) at t=1.5 (123 values)
 * then continues to t=11.5 and outputs again (123 values). Total 246 values.
 * Indices (0-based): 0, 267, 534, ..., 32574
 * First 62 are in u-block (<16384), next 61 in v-block (>=16384).
 * ---------------------------------------------------------------------------*/
#define NREF_NS128 123

static const sunrealtype yref_1p5_ns128[NREF_NS128] = {
  1.2007578117315e+00, 1.2071204234193e+00, 1.2127863884331e+00,
  1.2157617220458e+00, 1.2151093442108e+00, 1.2112905524236e+00,
  1.2056859410479e+00, 1.1998898328063e+00, 1.1952587894741e+00,
  1.1927956372396e+00, 1.1932396870857e+00, 1.1971503541783e+00,
  1.2055616497079e+00, 1.2169427714320e+00, 1.2295833293359e+00,
  1.2393063578306e+00, 1.2415908122779e+00, 1.2352778236779e+00,
  1.2238349046615e+00, 1.2122739363694e+00, 1.2042812695349e+00,
  1.2020684394042e+00, 1.2075151367927e+00, 1.2231292616291e+00,
  1.2547489060414e+00, 1.3016316028197e+00, 1.3582643581483e+00,
  1.3894456743653e+00, 1.3601498433737e+00, 1.3044776128510e+00,
  1.2581840281292e+00, 1.2286351913707e+00, 1.2143888573917e+00,
  1.2141325350876e+00, 1.2293197777817e+00, 1.2659356621240e+00,
  1.3356244413895e+00, 1.4674720333067e+00, 1.6053381251621e+00,
  1.5591470207945e+00, 1.3931224076493e+00, 1.2981426542882e+00,
  1.2464115103246e+00, 1.2196716965141e+00, 1.2094950385771e+00,
  1.2123225243568e+00, 1.2274555038523e+00, 1.2526076530490e+00,
  1.2866192927743e+00, 1.3143375219986e+00, 1.3153445935744e+00,
  1.2892222713513e+00, 1.2557010484650e+00, 1.2281973090240e+00,
  1.2099023942812e+00, 1.1998903217010e+00, 1.1965820263759e+00,
  1.1986949902357e+00, 1.2049686869748e+00, 1.2121925833509e+00,
  1.2194685996132e+00, 1.2229651212161e+00, 1.9700044196047e+00,
  1.9716505598087e+00, 1.9743223329114e+00, 1.9772588831630e+00,
  1.9798195166071e+00, 1.9815258237722e+00, 1.9820407096824e+00,
  1.9811753388835e+00, 1.9788000316493e+00, 1.9754532290811e+00,
  1.9717658000974e+00, 1.9687864667444e+00, 1.9675341211224e+00,
  1.9684449724633e+00, 1.9710748217853e+00, 1.9743718128204e+00,
  1.9771924556130e+00, 1.9785885155543e+00, 1.9778152391372e+00,
  1.9742920555027e+00, 1.9671207583570e+00, 1.9577053428037e+00,
  1.9478238004438e+00, 1.9416020656759e+00, 1.9426539587413e+00,
  1.9500277235397e+00, 1.9595045256959e+00, 1.9676425466124e+00,
  1.9727247599386e+00, 1.9739305111626e+00, 1.9705739025632e+00,
  1.9618403341202e+00, 1.9466785179803e+00, 1.9265546345509e+00,
  1.9080465450988e+00, 1.9065284293651e+00, 1.9231333660498e+00,
  1.9429661244960e+00, 1.9585068007395e+00, 1.9686174344299e+00,
  1.9735768477365e+00, 1.9737106567676e+00, 1.9690950951913e+00,
  1.9604235865377e+00, 1.9483827805512e+00, 1.9372501579054e+00,
  1.9334891165869e+00, 1.9397353621287e+00, 1.9512050086826e+00,
  1.9625007946002e+00, 1.9711642274952e+00, 1.9766459210976e+00,
  1.9790857246758e+00, 1.9787796010252e+00, 1.9761627798577e+00,
  1.9725896098149e+00, 1.9685877682870e+00, 1.9661032633447e+00,
  1.9661299554931e+00, 1.9684734345736e+00, 1.9721013899568e+00
};

static const sunrealtype yref_11p5_ns128[NREF_NS128] = {
  5.8626859760337e-01, 6.0614584136584e-01, 6.2274093987758e-01,
  6.3034801854856e-01, 6.2638049355771e-01, 6.1232572769272e-01,
  5.9268847004576e-01, 5.7302547784275e-01, 5.5829678427123e-01,
  5.5202331231966e-01, 5.5612011347714e-01, 5.7096226765928e-01,
  5.9722970804949e-01, 6.2856172496642e-01, 6.5942200300378e-01,
  6.8018926301181e-01, 6.8225786014605e-01, 6.6458061144073e-01,
  6.3543285485798e-01, 6.0664071695989e-01, 5.8777126589944e-01,
  5.8489359323551e-01, 6.0195287910865e-01, 6.4209212125657e-01,
  7.1303463456556e-01, 8.0461271610449e-01, 9.0166508808861e-01,
  9.4966574666909e-01, 9.0214958180332e-01, 8.0587890808635e-01,
  7.1709184558778e-01, 6.5574383667880e-01, 6.2575966420098e-01,
  6.2779892346122e-01, 6.6410744910453e-01, 7.4132317958896e-01,
  8.6792000780301e-01, 1.0642575793649e+00, 1.2444371049107e+00,
  1.1819249666518e+00, 9.5574674941750e-01, 7.9720829284789e-01,
  6.9600612893850e-01, 6.3832489370098e-01, 6.1586234298976e-01,
  6.2378770665632e-01, 6.5921960763540e-01, 7.1318117176321e-01,
  7.8026599821705e-01, 8.3063429159613e-01, 8.3143399187363e-01,
  7.8178061460561e-01, 7.1302229760172e-01, 6.5090706932807e-01,
  6.0593261032429e-01, 5.8013427047908e-01, 5.7227858112794e-01,
  5.7975092020609e-01, 5.9848875057030e-01, 6.1906071996385e-01,
  6.3883992038819e-01, 6.4772762967967e-01, 4.8324588657870e+00,
  4.8461284995064e+00, 4.8650839305965e+00, 4.8840991508353e+00,
  4.8988335187177e+00, 4.9064425154573e+00, 4.9056187082969e+00,
  4.8965079146691e+00, 4.8800616905729e+00, 4.8602320808578e+00,
  4.8413299375847e+00, 4.8287518984708e+00, 4.8267820185441e+00,
  4.8361781758601e+00, 4.8532552822269e+00, 4.8716354914919e+00,
  4.8849993526290e+00, 4.8886635121590e+00, 4.8797443334033e+00,
  4.8569240651317e+00, 4.8174663531700e+00, 4.7708116423158e+00,
  4.7256263053705e+00, 4.6998807992396e+00, 4.7078408448258e+00,
  4.7442411811947e+00, 4.7898080606581e+00, 4.8279906813859e+00,
  4.8498600348162e+00, 4.8513163992105e+00, 4.8300485026328e+00,
  4.7843438940209e+00, 4.7120894759109e+00, 4.6242511153373e+00,
  4.5495213295546e+00, 4.5459687708281e+00, 4.6169813888437e+00,
  4.7050991028430e+00, 4.7768974942999e+00, 4.8238579789964e+00,
  4.8452504556597e+00, 4.8422356350317e+00, 4.8162834461999e+00,
  4.7728138642143e+00, 4.7167411227554e+00, 4.6683467565610e+00,
  4.6546316267232e+00, 4.6857425230702e+00, 4.7415380012651e+00,
  4.7977477210612e+00, 4.8414538966210e+00, 4.8684502468829e+00,
  4.8786658925749e+00, 4.8738918540265e+00, 4.8574328660616e+00,
  4.8377318854693e+00, 4.8174306634359e+00, 4.8068583109958e+00,
  4.8106439164646e+00, 4.8272740052129e+00, 4.8507305111147e+00
};

/* ---------------------------------------------------------------------------
 * main
 * Usage: radau5_bruss_2d [rtol atol h0 use_schur nsmin nsmax [ns_grid]]
 * ---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  sunrealtype rtol = 1.0e-6;
  sunrealtype atol_val = 1.0e-6;
  sunrealtype h0   = 1.0e-4;
  int use_schur    = 0;
  int nsmin        = 3;
  int nsmax        = 13;
  int ns_grid      = 128;  /* default NS=128 to match Fortran reference */
  if (argc > 1) rtol      = atof(argv[1]);
  if (argc > 2) atol_val  = atof(argv[2]);
  if (argc > 3) h0        = atof(argv[3]);
  if (argc > 4) use_schur = atoi(argv[4]);
  if (argc > 5) nsmin     = atoi(argv[5]);
  if (argc > 6) nsmax     = atoi(argv[6]);
  if (argc > 7) ns_grid   = atoi(argv[7]);

  /* Set global dimensions */
  g_ns   = ns_grid;
  g_nssq = ns_grid * ns_grid;
  g_neq  = 2 * g_nssq;

  int ns   = g_ns;
  int nssq = g_nssq;
  int neq  = g_neq;

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  /* Initial conditions (matching Fortran driver) */
  N_Vector y0 = N_VNew_Serial(neq, sunctx);
  sunrealtype* y0v = N_VGetArrayPointer(y0);
  sunrealtype ans = (sunrealtype)ns;

  /* u-component: u(i,j) = 22*yy*(1-yy)^1.5, yy = j/ns (j=0..ns-1) */
  for (int j = 0; j < ns; j++) {
    sunrealtype yy = (sunrealtype)j / ans;
    sunrealtype uval = 22.0 * yy * pow(1.0 - yy, 1.5);
    for (int i = 0; i < ns; i++)
      y0v[j * ns + i] = uval;
  }
  /* v-component: v(i,j) = 27*xx*(1-xx)^1.5, xx = i/ns (i=0..ns-1) */
  for (int i = 0; i < ns; i++) {
    sunrealtype xx = (sunrealtype)i / ans;
    sunrealtype vval = 27.0 * xx * pow(1.0 - xx, 1.5);
    for (int j = 0; j < ns; j++)
      y0v[j * ns + i + nssq] = vval;
  }

  /* Build sparsity pattern */
  sunindextype* colptrs = (sunindextype*)malloc((neq + 1) * sizeof(sunindextype));
  sunindextype* rowinds_buf = (sunindextype*)malloc(
      (sunindextype)(neq * 6) * sizeof(sunindextype));
  sunindextype nnz = build_sparsity(ns, nssq, neq, colptrs, rowinds_buf);

  /* Create sparse Jacobian template */
  SUNMatrix Jt = SUNSparseMatrix(neq, neq, nnz, CSC_MAT, sunctx);
  sunindextype* cp = SM_INDEXPTRS_S(Jt);
  sunindextype* ri = SM_INDEXVALS_S(Jt);
  for (int c = 0; c <= neq; c++) cp[c] = colptrs[c];
  for (sunindextype k = 0; k < nnz; k++) ri[k] = rowinds_buf[k];
  sunrealtype* dat = SM_DATA_S(Jt);
  for (sunindextype k = 0; k < nnz; k++) dat[k] = 0.0;
  free(colptrs);
  free(rowinds_buf);

  /* Solver setup */
  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, nsmin, nsmax);
  Radau5Init(mem, rhs_bruss2d, 0.0, y0);
  Radau5SetLinearSolver(mem, Jt, NULL);
  Radau5SetSchurDecomp(mem, use_schur);
  Radau5SetSparsityPattern(mem, Jt);
  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);

  /* Phase 1: solve to t=1.5 */
  N_Vector yout = N_VNew_Serial(neq, sunctx);
  sunrealtype tret;
  int ret1 = Radau5Solve(mem, 1.5, yout, &tret);

  sunrealtype* yd = N_VGetArrayPointer(yout);
  sunrealtype maxerr1 = 0.0;

  /* Validate phase 1 against reference if NS=128 */
  if (ns == 128) {
    for (int r = 0; r < NREF_NS128; r++) {
      int idx = r * 267;
      if (idx >= neq) break;
      sunrealtype err = fabs(yd[idx] - yref_1p5_ns128[r]);
      if (err > maxerr1) maxerr1 = err;
    }
  }

  printf("=== Brusselator 2D (ns=%d, n=%d, rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n",
         ns, neq, rtol, atol_val, h0, use_schur);
  printf("Phase 1 (t=1.5): ret=%d", ret1);
  if (ns == 128) printf(", max err = %.6e", maxerr1);
  printf("\n");

  /* Phase 2: continue to t=11.5 */
  int ret2 = Radau5Solve(mem, 11.5, yout, &tret);
  yd = N_VGetArrayPointer(yout);
  sunrealtype maxerr2 = 0.0;

  if (ns == 128) {
    for (int r = 0; r < NREF_NS128; r++) {
      int idx = r * 267;
      if (idx >= neq) break;
      sunrealtype err = fabs(yd[idx] - yref_11p5_ns128[r]);
      if (err > maxerr2) maxerr2 = err;
    }
  }

  printf("Phase 2 (t=11.5): ret=%d", ret2);
  if (ns == 128) printf(", max err = %.6e", maxerr2);
  printf("\n");

  /* Print a few sampled values */
  printf("\nSampled solution at t=11.5:\n");
  int stride = (neq > 1000) ? neq / 8 : 1;
  for (int idx = 0; idx < neq && idx < 8 * stride; idx += stride) {
    const char* comp = (idx < nssq) ? "u" : "v";
    int idx2 = (idx < nssq) ? idx : idx - nssq;
    int gi = idx2 % ns, gj = idx2 / ns;
    printf("  y[%5d] = %s[%3d,%3d] = %18.12e\n", idx, comp, gi, gj, yd[idx]);
  }

  long int nstep, naccpt, nrejct, nfcn, njac, ndec;
  Radau5GetNumSteps(mem, &nstep);
  Radau5GetNumAccSteps(mem, &naccpt);
  Radau5GetNumRejSteps(mem, &nrejct);
  Radau5GetNumRhsEvals(mem, &nfcn);
  Radau5GetNumJacEvals(mem, &njac);
  Radau5GetNumDecomps(mem, &ndec);
  printf("\nnstep=%ld naccpt=%ld nrejct=%ld nfcn=%ld njac=%ld ndec=%ld\n",
         nstep, naccpt, nrejct, nfcn, njac, ndec);

  /* Pass criterion */
  int pass;
  if (ns == 128)
    pass = (ret1 == 0 && ret2 == 0 && maxerr1 < 5.0e-2 && maxerr2 < 5.0e-2);
  else
    pass = (ret1 == 0 && ret2 == 0);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout); SUNMatDestroy(Jt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}