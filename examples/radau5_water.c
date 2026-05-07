/* ---------------------------------------------------------------------------
 * radau5_water.c — Water Tube System (DAE index-2, n=49)
 *
 * M*y' = f(t,y) where M = diag(v,..,v, 0,..,0, c,c, 0,..,0)
 *   y[0..17]:  pipe flows (phi)        — M[i]=v (differential)
 *   y[18..35]: friction factors (lambda) — M[i]=0 (algebraic)
 *   y[36..37]: tank pressures           — M[i]=c (differential)
 *   y[38..48]: node pressures           — M[i]=0 (algebraic)
 *
 * DAE index: nind1=38, nind2=11, nind3=0
 * DQ Jacobian (no analytic Jac)
 *
 * t in [0, 61200] (17 hours in seconds)
 * Reference solution from IVPtestset (PSIDE, tol=1e-14)
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_types.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

#define NEQ    49
#define NNODES 13

/* Physical parameters (from water.f) */
#define NU     1.31e-6           /* kinematic viscosity [m^2/s] */
#define GG     9.8               /* gravity [m/s^2] */
#define RHO    1.0e3             /* density [kg/m^3] */
#define RCRIT  2.3e3             /* critical Reynolds number */
#define LENGTH 1.0e3             /* pipe length [m] */
#define KK     2.0e-4            /* roughness [m] */
#define DD     1.0               /* diameter [m] */
#define BB     2.0e2             /* tank parameter */
#define PI_VAL 3.141592653589793238462643383

/* Derived constants */
#define AA     (PI_VAL * DD * DD / 4.0)       /* cross-section area */
#define MU     (NU * RHO)                     /* dynamic viscosity */
#define VV     (RHO * LENGTH / AA)            /* mass diagonal for flows */
#define CC     (BB / (RHO * GG))              /* mass diagonal for tanks */

/* Pipe network topology: 18 pipes connecting 13 nodes.
 * Each pipe k has endpoints (from_node[k], to_node[k]) (1-based node indices).
 * phi(from,to) = y[k], lambda(from,to) = y[k+18], pressures p(node) mapped below.
 */
/* Pipe topology: pipe k connects from_node[k] -> to_node[k] (0-based node index) */
static const int from_node[18] = {
  0, 1, 1, 2, 2, 3, 4, 5, 6, 6, 7, 7, 8, 10, 10, 11, 11, 12
};
static const int to_node[18] = {
  1, 2, 5, 3, 4, 4, 9, 4, 3, 7, 4, 9, 7, 8, 11, 6, 7, 10
};

/* Pressure mapping: p(node) = y[p_idx[node]] (0-based) */
/* Nodes 0-based: 0..12.  Tanks are nodes 4 and 7 (y[36],y[37]).
 * Remaining nodes map to y[38..48]. */
static const int p_idx[NNODES] = {
  /* node 0 */ 38,
  /* node 1 */ 39,
  /* node 2 */ 40,
  /* node 3 */ 41,
  /* node 4 */ 36,  /* tank */
  /* node 5 */ 42,
  /* node 6 */ 43,
  /* node 7 */ 37,  /* tank */
  /* node 8 */ 44,
  /* node 9 */ 45,
  /* node10 */ 46,
  /* node11 */ 47,
  /* node12 */ 48
};

/* ---------------------------------------------------------------------------
 * RHS: M*y' = f(t, y)
 *
 * Faithfully translated from water.f feval subroutine.
 * Uses 2D arrays phi[i][j], lambda[i][j], fdba[i][j], rghres[i][j],
 * netflo[n] indexed by 0-based node indices.
 * ---------------------------------------------------------------------------*/
static int rhs_water(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)ud;
  sunrealtype* yv = N_VGetArrayPointer(y);
  sunrealtype* f  = N_VGetArrayPointer(yd);

  const sunrealtype a = PI_VAL * DD * DD / 4.0;
  const sunrealtype mu = NU * RHO;

  int i, j, n;

  /* Local 2D arrays (0-based, NNODES x NNODES) */
  sunrealtype phi[NNODES][NNODES];
  sunrealtype lam[NNODES][NNODES];
  sunrealtype p[NNODES];
  sunrealtype fdba[NNODES][NNODES];
  sunrealtype netflo[NNODES];
  sunrealtype ein[NNODES], eout[NNODES];
  sunrealtype rghres[NNODES][NNODES];

  for (j = 0; j < NNODES; j++) {
    ein[j]  = 0.0;
    eout[j] = 0.0;
  }

  /* Time-dependent inflows/outflows */
  sunrealtype that  = t / 3600.0;
  sunrealtype that2 = that * that;

  ein[0]  = (1.0 - cos(exp(-that) - 1.0)) / 200.0;   /* node 1 (0-based: 0) */
  ein[12] = (1.0 - cos(exp(-that) - 1.0)) / 80.0;    /* node 13 (0-based: 12) */
  eout[9] = that2 * (3.0 * that2 - 92.0 * that + 720.0) / 1.0e6; /* node 10 (0-based: 9) */

  /* Initialize phi and lambda arrays */
  for (j = 0; j < NNODES; j++)
    for (i = 0; i < NNODES; i++) {
      phi[i][j] = 0.0;
      lam[i][j] = 1.0;
    }

  /* Map state vector into phi(i,j) — Fortran 1-based to C 0-based */
  phi[0][1]   = yv[0];
  phi[1][2]   = yv[1];
  phi[1][5]   = yv[2];
  phi[2][3]   = yv[3];
  phi[2][4]   = yv[4];
  phi[3][4]   = yv[5];
  phi[4][9]   = yv[6];
  phi[5][4]   = yv[7];
  phi[6][3]   = yv[8];
  phi[6][7]   = yv[9];
  phi[7][4]   = yv[10];
  phi[7][9]   = yv[11];
  phi[8][7]   = yv[12];
  phi[10][8]  = yv[13];
  phi[10][11] = yv[14];
  phi[11][6]  = yv[15];
  phi[11][7]  = yv[16];
  phi[12][10] = yv[17];

  /* Map state vector into lambda(i,j) */
  lam[0][1]   = yv[18];
  lam[1][2]   = yv[19];
  lam[1][5]   = yv[20];
  lam[2][3]   = yv[21];
  lam[2][4]   = yv[22];
  lam[3][4]   = yv[23];
  lam[4][9]   = yv[24];
  lam[5][4]   = yv[25];
  lam[6][3]   = yv[26];
  lam[6][7]   = yv[27];
  lam[7][4]   = yv[28];
  lam[7][9]   = yv[29];
  lam[8][7]   = yv[30];
  lam[10][8]  = yv[31];
  lam[10][11] = yv[32];
  lam[11][6]  = yv[33];
  lam[11][7]  = yv[34];
  lam[12][10] = yv[35];

  /* Map pressures: p(node) from state vector */
  p[4]  = yv[36];
  p[7]  = yv[37];
  p[0]  = yv[38];
  p[1]  = yv[39];
  p[2]  = yv[40];
  p[3]  = yv[41];
  p[5]  = yv[42];
  p[6]  = yv[43];
  p[8]  = yv[44];
  p[9]  = yv[45];
  p[10] = yv[46];
  p[11] = yv[47];
  p[12] = yv[48];

  /* Compute fdba and rghres for all node pairs */
  for (j = 0; j < NNODES; j++) {
    for (i = 0; i < NNODES; i++) {
      if (lam[i][j] < 0.0) return 1;  /* recoverable: sqrt of negative */
      sunrealtype rtla = sqrt(lam[i][j]);
      sunrealtype r = fabs(phi[i][j] * DD / (NU * a));
      if (r > RCRIT) {
        rghres[i][j] = 1.0 / rtla - 1.74
                       + 2.0 * log10(2.0 * KK / DD + 18.7 / (r * rtla));
        fdba[i][j] = p[i] - p[j]
                     - lam[i][j] * RHO * LENGTH * phi[i][j] * phi[i][j]
                       / (a * a * DD);
      } else {
        rghres[i][j] = 1.0 / rtla - 1.74
                       + 2.0 * log10(2.0 * KK / DD + 18.7 / (RCRIT * rtla));
        fdba[i][j] = p[i] - p[j]
                     - 32.0 * mu * LENGTH * phi[i][j] / (a * DD * DD);
      }
    }
  }

  /* Compute net flow at each node */
  for (n = 0; n < NNODES; n++) {
    netflo[n] = ein[n] - eout[n];
    for (i = 0; i < NNODES; i++)
      netflo[n] += phi[i][n];   /* inflow to node n */
    for (j = 0; j < NNODES; j++)
      netflo[n] -= phi[n][j];   /* outflow from node n */
  }

  /* f[0..17]: flow equations (pressure diff - friction loss) */
  f[0]  = fdba[0][1];
  f[1]  = fdba[1][2];
  f[2]  = fdba[1][5];
  f[3]  = fdba[2][3];
  f[4]  = fdba[2][4];
  f[5]  = fdba[3][4];
  f[6]  = fdba[4][9];
  f[7]  = fdba[5][4];
  f[8]  = fdba[6][3];
  f[9]  = fdba[6][7];
  f[10] = fdba[7][4];
  f[11] = fdba[7][9];
  f[12] = fdba[8][7];
  f[13] = fdba[10][8];
  f[14] = fdba[10][11];
  f[15] = fdba[11][6];
  f[16] = fdba[11][7];
  f[17] = fdba[12][10];

  /* f[18..35]: friction factor equations (Colebrook-White) */
  f[18] = rghres[0][1];
  f[19] = rghres[1][2];
  f[20] = rghres[1][5];
  f[21] = rghres[2][3];
  f[22] = rghres[2][4];
  f[23] = rghres[3][4];
  f[24] = rghres[4][9];
  f[25] = rghres[5][4];
  f[26] = rghres[6][3];
  f[27] = rghres[6][7];
  f[28] = rghres[7][4];
  f[29] = rghres[7][9];
  f[30] = rghres[8][7];
  f[31] = rghres[10][8];
  f[32] = rghres[10][11];
  f[33] = rghres[11][6];
  f[34] = rghres[11][7];
  f[35] = rghres[12][10];

  /* f[36..48]: pressure/flow balance at nodes */
  f[36] = netflo[4];   /* tank node 5 (0-based: 4) */
  f[37] = netflo[7];   /* tank node 8 (0-based: 7) */
  f[38] = netflo[0];
  f[39] = netflo[1];
  f[40] = netflo[2];
  f[41] = netflo[3];
  f[42] = netflo[5];
  f[43] = netflo[6];
  f[44] = netflo[8];
  f[45] = netflo[9];
  f[46] = netflo[10];
  f[47] = netflo[11];
  f[48] = netflo[12];

  return 0;
}

/* ---------------------------------------------------------------------------
 * Mass matrix: M = diag(v,..,v, 0,..,0, c,c, 0,..,0)
 *   M[0..17]  = v = rho*length/a  (pipe flow inertia)
 *   M[18..35] = 0                  (algebraic: friction factors)
 *   M[36..37] = c = b/(rho*g)     (tank capacitance)
 *   M[38..48] = 0                  (algebraic: pressure/flow balance)
 * ---------------------------------------------------------------------------*/
static int mas_water(sunrealtype t, SUNMatrix M, void* ud,
                     N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)ud; (void)t1; (void)t2; (void)t3;
  SUNMatZero(M);

  /* Pipe flows: M[0..17] = v */
  for (int i = 0; i < 18; i++)
    SM_ELEMENT_D(M, i, i) = VV;

  /* Friction factors: M[18..35] = 0 (already zero) */

  /* Tank pressures: M[36..37] = c */
  SM_ELEMENT_D(M, 36, 36) = CC;
  SM_ELEMENT_D(M, 37, 37) = CC;

  /* Node pressures: M[38..48] = 0 (already zero) */

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

  void* mem = Radau5Create(sunctx);

  /* Initial conditions (from water.f init subroutine, 0-based) */
  N_Vector y0 = N_VNew_Serial(NEQ, sunctx);
  sunrealtype* y0v = N_VGetArrayPointer(y0);
  int i;

  /* y[0..17] = 0 (flows) */
  for (i = 0; i < 18; i++) y0v[i] = 0.0;
  /* y[18..35] = 0.47519404529185289807e-1 (friction factors) */
  for (i = 18; i < 36; i++) y0v[i] = 0.47519404529185289807e-1;
  /* y[36..48] = 109800 (pressures) */
  for (i = 36; i < 49; i++) y0v[i] = 109800.0;

  Radau5Init(mem, rhs_water, 0.0, y0);

  SUNMatrix Jt = SUNDenseMatrix(NEQ, NEQ, sunctx);
  Radau5SetLinearSolver(mem, Jt);
  Radau5SetSchurDecomp(mem, use_schur);
  /* No Radau5SetJacFn — use DQ Jacobian */

  /* Mass matrix */
  SUNMatrix Mt = SUNDenseMatrix(NEQ, NEQ, sunctx);
  Radau5SetMassFn(mem, mas_water, Mt);

  /* DAE index: nind1=38, nind2=11, nind3=0 */
  Radau5SetDAEIndex(mem, 38, 11, 0);

  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);

  N_Vector yout = N_VNew_Serial(NEQ, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 61200.0, yout, &tret);

  /* Reference solution at t=61200 (PSIDE, tol=1e-14, from water.f solut) */
  static const sunrealtype yref[NEQ] = {
    0.2298488296477430e-02,
    0.1188984650746585e-02,
    0.1109503645730845e-02,
    0.1589620100314825e-03,
    0.1030022640715102e-02,
    0.8710606306836165e-03,
    0.3243571480903489e-02,
    0.1109503645730845e-02,
    0.7120986206521341e-03,
    0.6414613963833099e-03,
    0.9416978549524347e-03,
    0.3403428519096511e-02,
    0.2397639310739395e-02,
    0.2397639310739395e-02,
    0.3348581430454180e-02,
    0.1353560017035444e-02,
    0.1995021413418736e-02,
    0.5746220741193575e-02,
    0.4751940452918529e-01,
    0.4751940452918529e-01,
    0.4751940452918529e-01,
    0.4751940452918529e-01,
    0.4751940452918529e-01,
    0.4751940452918529e-01,
    0.4311196778792902e-01,
    0.4751940452918529e-01,
    0.4751940452918529e-01,
    0.4751940452918529e-01,
    0.4751940452918529e-01,
    0.4249217433601160e-01,
    0.4732336439609648e-01,
    0.4732336439609648e-01,
    0.4270002118868241e-01,
    0.4751940452918529e-01,
    0.4751940452918529e-01,
    0.3651427026675656e-01,
    0.1111268591478108e+06,
    0.1111270045592387e+06,
    0.1111271078730254e+06,
    0.1111269851929858e+06,
    0.1111269255355337e+06,
    0.1111269322658045e+06,
    0.1111269221703983e+06,
    0.1111270121140691e+06,
    0.1111274419515807e+06,
    0.1111255158881087e+06,
    0.1111278793439227e+06,
    0.1111270995171642e+06,
    0.1111298338971779e+06
  };

  sunrealtype* yd = N_VGetArrayPointer(yout);
  printf("=== Water Tube System (rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n",
         rtol, atol_val, h0, use_schur);
  printf("ret=%d, tret=%.6e\n", ret, tret);

  sunrealtype maxrelerr = 0.0;
  for (i = 0; i < NEQ; i++) {
    sunrealtype err = fabs(yd[i] - yref[i]);
    sunrealtype relerr = (fabs(yref[i]) > 1e-15) ? err / fabs(yref[i]) : err;
    if (relerr > maxrelerr) maxrelerr = relerr;
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

  /* Pass criterion: max relative error < 1e-2 */
  int pass = (ret == 0 && maxrelerr < 1.0e-2);
  printf("max_rel_err=%.3e\n", maxrelerr);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout);
  SUNMatDestroy(Jt); SUNMatDestroy(Mt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
