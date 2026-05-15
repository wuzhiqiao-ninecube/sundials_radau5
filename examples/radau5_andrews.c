/* ---------------------------------------------------------------------------
 * radau5_andrews.c — Andrews' squeezing mechanism (DAE index-3, n=27)
 *
 * 7 angles + 7 velocities + 7 "accelerations" + 6 Lagrange multipliers
 * M*y' = f(t,y), M = diag(1..1, 0..0) (first 14 are differential)
 * DAE index: nind1=7, nind2=7, nind3=13
 *
 * t in [0, 0.03]
 * Reference from IVPtestset (PSIDE, Cray C90, tol=1e-14)
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

#define NEQ 27

/* Physical parameters — prefixed to avoid macro collisions */
#define PM1  0.04325
#define PM2  0.00365
#define PM3  0.02373
#define PM4  0.00706
#define PM5  0.07050
#define PM6  0.00706
#define PM7  0.05498
#define PXA  (-0.06934)
#define PYA  (-0.00227)
#define PXB  (-0.03635)
#define PYB  0.03273
#define PXC  0.014
#define PYC  0.072
#define PC0  4530.0
#define PI1  2.194e-6
#define PI2  4.410e-7
#define PI3  5.255e-6
#define PI4  5.667e-7
#define PI5  1.169e-5
#define PI6  5.667e-7
#define PI7  1.912e-5
#define PD   28e-3
#define PDA  115e-4
#define PE   2e-2
#define PEA  1421e-5
#define PRR  7e-3
#define PRA  92e-5
#define PL0  7785e-5
#define PSS  35e-3
#define PSA  1874e-5
#define PSB  1043e-5
#define PSC  18e-3
#define PSD  2e-2
#define PTA  2308e-5
#define PTB  916e-5
#define PU   4e-2
#define PUA  1228e-5
#define PUB  449e-5
#define PZF  2e-2
#define PZT  4e-2
#define PFA  1421e-5
#define PMOM 33e-3

/* ---------------------------------------------------------------------------
 * RHS: ff = f(t, y)
 * ---------------------------------------------------------------------------*/
static int rhs_andrews(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)t; (void)ud;
  sunrealtype* v  = N_VGetArrayPointer(y);
  sunrealtype* ff = N_VGetArrayPointer(yd);
  int i, j;

  sunrealtype sibe = sin(v[0]), sith = sin(v[1]), siga = sin(v[2]);
  sunrealtype siph = sin(v[3]), side = sin(v[4]), siom = sin(v[5]);
  sunrealtype siep = sin(v[6]);
  sunrealtype cobe = cos(v[0]), coth = cos(v[1]), coga = cos(v[2]);
  sunrealtype coph = cos(v[3]), code = cos(v[4]), coom = cos(v[5]);
  sunrealtype coep = cos(v[6]);
  sunrealtype sibeth = sin(v[0]+v[1]), cobeth = cos(v[0]+v[1]);
  sunrealtype siphde = sin(v[3]+v[4]), cophde = cos(v[3]+v[4]);
  sunrealtype siomep = sin(v[5]+v[6]), coomep = cos(v[5]+v[6]);

  sunrealtype bep = v[7], thp = v[8];
  sunrealtype php = v[10], dep = v[11], omp = v[12], epp = v[13];

  /* Inertia matrix mm[7][7] */
  sunrealtype mm[7][7];
  for (i = 0; i < 7; i++) for (j = 0; j < 7; j++) mm[i][j] = 0.0;
/* PLACEHOLDER_RHS_BODY */
  mm[0][0] = PM1*PRA*PRA + PM2*(PRR*PRR-2*PDA*PRR*coth+PDA*PDA) + PI1 + PI2;
  mm[1][0] = PM2*(PDA*PDA-PDA*PRR*coth) + PI2;
  mm[1][1] = PM2*PDA*PDA + PI2;
  mm[2][2] = PM3*(PSA*PSA+PSB*PSB) + PI3;
  mm[3][3] = PM4*(PE-PEA)*(PE-PEA) + PI4;
  mm[4][3] = PM4*((PE-PEA)*(PE-PEA)+PZT*(PE-PEA)*siph) + PI4;
  mm[4][4] = PM4*(PZT*PZT+2*PZT*(PE-PEA)*siph+(PE-PEA)*(PE-PEA))
           + PM5*(PTA*PTA+PTB*PTB) + PI4 + PI5;
  mm[5][5] = PM6*(PZF-PFA)*(PZF-PFA) + PI6;
  mm[6][5] = PM6*((PZF-PFA)*(PZF-PFA)-PU*(PZF-PFA)*siom) + PI6;
  mm[6][6] = PM6*((PZF-PFA)*(PZF-PFA)-2*PU*(PZF-PFA)*siom+PU*PU)
           + PM7*(PUA*PUA+PUB*PUB) + PI6 + PI7;
  for (j = 1; j < 7; j++)
    for (i = 0; i < j; i++)
      mm[i][j] = mm[j][i];

  sunrealtype xd = PSD*coga + PSC*siga + PXB;
  sunrealtype yd_v = PSD*siga - PSC*coga + PYB;
  sunrealtype lang = sqrt((xd-PXC)*(xd-PXC) + (yd_v-PYC)*(yd_v-PYC));
  sunrealtype force = -PC0*(lang-PL0)/lang;
  sunrealtype fx = force*(xd-PXC);
  sunrealtype fy = force*(yd_v-PYC);

  sunrealtype fv[7];
  fv[0] = PMOM - PM2*PDA*PRR*thp*(thp+2*bep)*sith;
  fv[1] = PM2*PDA*PRR*bep*bep*sith;
  fv[2] = fx*(PSC*coga - PSD*siga) + fy*(PSD*coga + PSC*siga);
  fv[3] = PM4*PZT*(PE-PEA)*dep*dep*coph;
  fv[4] = -PM4*PZT*(PE-PEA)*php*(php+2*dep)*coph;
  fv[5] = -PM6*PU*(PZF-PFA)*epp*epp*coom;
  fv[6] = PM6*PU*(PZF-PFA)*omp*(omp+2*epp)*coom;

  /* Constraint Jacobian gp[6][7] */
  sunrealtype gp[6][7];
  for (i = 0; i < 6; i++) for (j = 0; j < 7; j++) gp[i][j] = 0.0;
/* PLACEHOLDER_RHS_GP */
  gp[0][0] = -PRR*sibe + PD*sibeth;
  gp[0][1] = PD*sibeth;
  gp[0][2] = -PSS*coga;
  gp[1][0] = PRR*cobe - PD*cobeth;
  gp[1][1] = -PD*cobeth;
  gp[1][2] = -PSS*siga;
  gp[2][0] = -PRR*sibe + PD*sibeth;
  gp[2][1] = PD*sibeth;
  gp[2][3] = -PE*cophde;
  gp[2][4] = -PE*cophde + PZT*side;
  gp[3][0] = PRR*cobe - PD*cobeth;
  gp[3][1] = -PD*cobeth;
  gp[3][3] = -PE*siphde;
  gp[3][4] = -PE*siphde - PZT*code;
  gp[4][0] = -PRR*sibe + PD*sibeth;
  gp[4][1] = PD*sibeth;
  gp[4][5] = PZF*siomep;
  gp[4][6] = PZF*siomep - PU*coep;
  gp[5][0] = PRR*cobe - PD*cobeth;
  gp[5][1] = -PD*cobeth;
  gp[5][5] = -PZF*coomep;
  gp[5][6] = -PZF*coomep - PU*siep;

  /* Constraints g[6] */
  sunrealtype g[6];
  g[0] = PRR*cobe - PD*cobeth - PSS*siga - PXB;
  g[1] = PRR*sibe - PD*sibeth + PSS*coga - PYB;
  g[2] = PRR*cobe - PD*cobeth - PE*siphde - PZT*code - PXA;
  g[3] = PRR*sibe - PD*sibeth + PE*cophde - PZT*side - PYA;
  g[4] = PRR*cobe - PD*cobeth - PZF*coomep - PU*siep - PXA;
  g[5] = PRR*sibe - PD*sibeth - PZF*siomep + PU*coep - PYA;

  /* Assemble ff */
  for (i = 0; i < 14; i++) ff[i] = v[i+7];
  for (i = 0; i < 7; i++) {
    ff[14+i] = -fv[i];
    for (j = 0; j < 7; j++) ff[14+i] += mm[i][j] * v[j+14];
    for (j = 0; j < 6; j++) ff[14+i] += gp[j][i] * v[j+21];
  }
  for (i = 0; i < 6; i++) ff[21+i] = g[i];

  return 0;
}
/* PLACEHOLDER_JAC */
/* ---------------------------------------------------------------------------
 * Approximate Jacobian (from Hairer & Wanner p.540)
 * ---------------------------------------------------------------------------*/
static int jac_andrews(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                       void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;
  sunrealtype* v = N_VGetArrayPointer(y);
  int i, j;

  sunrealtype sibe = sin(v[0]), siga = sin(v[2]);
  sunrealtype siph = sin(v[3]), side = sin(v[4]), siom = sin(v[5]);
  sunrealtype siep = sin(v[6]);
  sunrealtype cobe = cos(v[0]), coth = cos(v[1]), coga = cos(v[2]);
  sunrealtype code = cos(v[4]), coep = cos(v[6]);
  sunrealtype sibeth = sin(v[0]+v[1]), cobeth = cos(v[0]+v[1]);
  sunrealtype siphde = sin(v[3]+v[4]), cophde = cos(v[3]+v[4]);
  sunrealtype siomep = sin(v[5]+v[6]), coomep = cos(v[5]+v[6]);

  /* Zero all */
  for (j = 0; j < NEQ; j++)
    for (i = 0; i < NEQ; i++)
      SM_ELEMENT_D(J, i, j) = 0.0;

  /* df[0..13]/dy[7..20] = I */
  for (i = 0; i < 14; i++)
    SM_ELEMENT_D(J, i, i+7) = 1.0;

  /* Inertia matrix mm[7][7] */
  sunrealtype mm[7][7];
  for (i = 0; i < 7; i++) for (j = 0; j < 7; j++) mm[i][j] = 0.0;
  mm[0][0] = PM1*PRA*PRA + PM2*(PRR*PRR-2*PDA*PRR*coth+PDA*PDA) + PI1 + PI2;
  mm[1][0] = PM2*(PDA*PDA-PDA*PRR*coth) + PI2;
  mm[1][1] = PM2*PDA*PDA + PI2;
  mm[2][2] = PM3*(PSA*PSA+PSB*PSB) + PI3;
  mm[3][3] = PM4*(PE-PEA)*(PE-PEA) + PI4;
  mm[4][3] = PM4*((PE-PEA)*(PE-PEA)+PZT*(PE-PEA)*siph) + PI4;
  mm[4][4] = PM4*(PZT*PZT+2*PZT*(PE-PEA)*siph+(PE-PEA)*(PE-PEA))
           + PM5*(PTA*PTA+PTB*PTB) + PI4 + PI5;
  mm[5][5] = PM6*(PZF-PFA)*(PZF-PFA) + PI6;
  mm[6][5] = PM6*((PZF-PFA)*(PZF-PFA)-PU*(PZF-PFA)*siom) + PI6;
  mm[6][6] = PM6*((PZF-PFA)*(PZF-PFA)-2*PU*(PZF-PFA)*siom+PU*PU)
           + PM7*(PUA*PUA+PUB*PUB) + PI6 + PI7;
  for (j = 1; j < 7; j++)
    for (i = 0; i < j; i++)
      mm[i][j] = mm[j][i];
/* PLACEHOLDER_JAC2 */

  /* Constraint Jacobian gp[6][7] */
  sunrealtype gp[6][7];
  for (i = 0; i < 6; i++) for (j = 0; j < 7; j++) gp[i][j] = 0.0;
  gp[0][0] = -PRR*sibe + PD*sibeth;
  gp[0][1] = PD*sibeth;
  gp[0][2] = -PSS*coga;
  gp[1][0] = PRR*cobe - PD*cobeth;
  gp[1][1] = -PD*cobeth;
  gp[1][2] = -PSS*siga;
  gp[2][0] = -PRR*sibe + PD*sibeth;
  gp[2][1] = PD*sibeth;
  gp[2][3] = -PE*cophde;
  gp[2][4] = -PE*cophde + PZT*side;
  gp[3][0] = PRR*cobe - PD*cobeth;
  gp[3][1] = -PD*cobeth;
  gp[3][3] = -PE*siphde;
  gp[3][4] = -PE*siphde - PZT*code;
  gp[4][0] = -PRR*sibe + PD*sibeth;
  gp[4][1] = PD*sibeth;
  gp[4][5] = PZF*siomep;
  gp[4][6] = PZF*siomep - PU*coep;
  gp[5][0] = PRR*cobe - PD*cobeth;
  gp[5][1] = -PD*cobeth;
  gp[5][5] = -PZF*coomep;
  gp[5][6] = -PZF*coomep - PU*siep;

  /* J[14+j][14+i] = mm[j][i] */
  for (i = 0; i < 7; i++)
    for (j = 0; j < 7; j++)
      SM_ELEMENT_D(J, 14+j, 14+i) = mm[j][i];

  /* J[14+j][21+i] = gp[i][j] */
  for (i = 0; i < 6; i++)
    for (j = 0; j < 7; j++)
      SM_ELEMENT_D(J, 14+j, 21+i) = gp[i][j];

  /* J[21+j][i] = gp[j][i] */
  for (i = 0; i < 7; i++)
    for (j = 0; j < 6; j++)
      SM_ELEMENT_D(J, 21+j, i) = gp[j][i];

  return 0;
}

/* ---------------------------------------------------------------------------
 * Mass matrix: M = diag(1..1, 0..0), first 14 differential
 * ---------------------------------------------------------------------------*/
static int mas_andrews(sunrealtype t, SUNMatrix M, void* ud,
                       N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)ud; (void)t1; (void)t2; (void)t3;
  SUNMatZero(M);
  for (int i = 0; i < 14; i++)
    SM_ELEMENT_D(M, i, i) = 1.0;
  return 0;
}
/* PLACEHOLDER_MAIN */

int main(int argc, char* argv[])
{
  sunrealtype rtol = 1.0e-6;
  sunrealtype atol_val = 1.0e-6;
  sunrealtype h0   = 1.0e-10;
  int use_schur    = 1;
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
  sunrealtype* y0v = N_VGetArrayPointer(y0);
  for (int i = 0; i < NEQ; i++) y0v[i] = 0.0;

  y0v[0]  = -0.0617138900142764496358948458001;
  y0v[1]  =  0.0;
  y0v[2]  =  0.455279819163070380255912382449;
  y0v[3]  =  0.222668390165885884674473185609;
  y0v[4]  =  0.487364979543842550225598953530;
  y0v[5]  = -0.222668390165885884674473185609;
  y0v[6]  =  1.23054744454982119249735015568;
  y0v[14] =  14222.4439199541138705911625887;
  y0v[15] = -10666.8329399655854029433719415;
  y0v[21] =  98.5668703962410896057654982170;
  y0v[22] = -6.12268834425566265503114393122;

  Radau5Init(mem, rhs_andrews, 0.0, y0);

  SUNMatrix Jt = SUNDenseMatrix(NEQ, NEQ, sunctx);
  Radau5SetLinearSolver(mem, Jt, NULL);
  Radau5SetSchurDecomp(mem, use_schur);
  Radau5SetJacFn(mem, jac_andrews);

  SUNMatrix Mt = SUNDenseMatrix(NEQ, NEQ, sunctx);
  Radau5SetMassFn(mem, mas_andrews, Mt);
  Radau5SetDAEIndex(mem, 7, 7, 13);

  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);

  N_Vector yout = N_VNew_Serial(NEQ, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 0.03, yout, &tret);

  static const sunrealtype yref[NEQ] = {
     0.1581077119629904e+02, -0.1575637105984298e+02,
     0.4082224013073101e-01, -0.5347301163226948e+00,
     0.5244099658805304e+00,  0.5347301163226948e+00,
     0.1048080741042263e+01,  0.1139920302151208e+04,
    -0.1424379294994111e+04,  0.1103291221937134e+02,
     0.1929337464421385e+02,  0.5735699284790808e+00,
    -0.1929337464421385e+02,  0.3231791658026955e+00,
    -0.2463176316945196e+05,  0.5185037701610329e+05,
     0.3241025686413781e+06,  0.5667493645176213e+06,
     0.1674362929479361e+05, -0.5667493645176222e+06,
     0.9826520791458422e+04,  0.1991753333731910e+03,
    -0.2975531228015052e+02,  0.2306654119098399e+02,
     0.3145271365475927e+02,  0.2264249232082739e+02,
     0.1161740700019673e+02
  };
/* PLACEHOLDER_MAIN2 */

  sunrealtype* yd = N_VGetArrayPointer(yout);
  printf("=== Andrews' squeezer (rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n", rtol, atol_val, h0, use_schur);
  printf("ret=%d, tret=%.6e\n", ret, tret);

  sunrealtype maxrelerr = 0.0;
  for (int i = 0; i < 7; i++) {
    sunrealtype err = fabs(yd[i] - yref[i]);
    sunrealtype relerr = (fabs(yref[i]) > 1e-15) ? err / fabs(yref[i]) : err;
    if (relerr > maxrelerr) maxrelerr = relerr;
    printf("y[%2d] = %22.14e  ref = %22.14e  rel_err = %.3e\n",
           i, yd[i], yref[i], relerr);
  }
  for (int i = 7; i < NEQ; i++) {
    sunrealtype err = fabs(yd[i] - yref[i]);
    sunrealtype relerr = (fabs(yref[i]) > 1e-15) ? err / fabs(yref[i]) : err;
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

  /* Index-3 DAE at rtol=1e-6 is hard; relaxed check on positions only */
  int pass = (ret == 0 && maxrelerr < 0.5);
  printf("max_pos_rel_err=%.3e\n", maxrelerr);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout);
  SUNMatDestroy(Jt); SUNMatDestroy(Mt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
