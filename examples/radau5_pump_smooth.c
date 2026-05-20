/*
 * radau5_pump_smooth.c — Charge pump DAE-2 with smoothed charge functions.
 *
 * The original charge functions QG, QS, QD have C0-but-not-C1 transitions
 * at the saturation/linear boundary (max(ugd-vte,0)) and the depletion/
 * inversion boundary (ugs=vte). This version replaces those hard switches
 * with softplus/sigmoid smoothing to produce C-infinity right-hand sides
 * and a consistent analytic Jacobian.
 *
 * Based on pumpdae_smooth.m (MATLAB) and radau5_pump_test.c (C original).
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_types.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

/* Physical parameters */
#define VT0   0.20
#define GAMMA 0.035
#define PHI   1.01
#define COX   4.0e-12
#define CAPD  0.40e-12
#define CAPS  1.60e-12
#define VHIGH  20.0
#define DELTAT 120.0e-9
#define T1_PLS 50.0e-9
#define T2_PLS 60.0e-9
#define T3_PLS 110.0e-9

/* Smoothing parameter (voltage scale for transition width) */
#define EPS_SMOOTH 1e-4

/* =========================================================================
 * Smoothing helper functions
 * ========================================================================= */

static sunrealtype softplus(sunrealtype x, sunrealtype eps)
{
  /* softplus(x, eps) = eps * log(1 + exp(x/eps)), numerically stable */
  sunrealtype z = x / eps;
  if (z > 30.0) return x;
  if (z < -30.0) return 0.0;
  return eps * log(1.0 + exp(z));
}

static sunrealtype dsoftplus(sunrealtype x, sunrealtype eps)
{
  /* Derivative of softplus: sigmoid(x/eps) */
  sunrealtype z = x / eps;
  if (z > 30.0) return 1.0;
  if (z < -30.0) return 0.0;
  return 1.0 / (1.0 + exp(-z));
}

static sunrealtype sigmoid_s(sunrealtype x, sunrealtype eps)
{
  /* Smooth step: 1/(1+exp(-x/eps)) */
  sunrealtype z = x / eps;
  if (z > 30.0) return 1.0;
  if (z < -30.0) return 0.0;
  return 1.0 / (1.0 + exp(-z));
}

static sunrealtype dsigmoid_s(sunrealtype x, sunrealtype eps)
{
  /* Derivative of sigmoid: sig*(1-sig)/eps */
  sunrealtype s = sigmoid_s(x, eps);
  return s * (1.0 - s) / eps;
}

/* =========================================================================
 * Smoothed charge functions
 * ========================================================================= */

static sunrealtype fn_qgate_smooth(sunrealtype vgb, sunrealtype vgs,
                                   sunrealtype vgd)
{
  sunrealtype ugs, ugd, ugb, ubs, vfb, vte;
  sunrealtype q_A, q_B, q_C, arg_B, ugst, ugdt, S, H, sigma_BC;

  /* Swap so ugs >= ugd */
  if ((vgs - vgd) <= 0.0) { ugs = vgd; ugd = vgs; }
  else { ugs = vgs; ugd = vgd; }
  ugb = vgb;
  ubs = ugs - ugb;

  vfb = VT0 - GAMMA * sqrt(PHI) - PHI;
  vte = VT0 + GAMMA * (sqrt(PHI - ubs) - sqrt(PHI));

  /* Region A */
  q_A = COX * (ugb - vfb);

  /* Region B */
  arg_B = (GAMMA / 2.0) * (GAMMA / 2.0) + ugb - vfb;
  if (arg_B < 1e-30) arg_B = 1e-30;
  q_B = COX * GAMMA * (sqrt(arg_B) - GAMMA / 2.0);

  /* Region C (smoothed) */
  ugst = softplus(ugs - vte, EPS_SMOOTH);
  ugdt = softplus(ugd - vte, EPS_SMOOTH);
  S = ugst + ugdt;
  if (S > 1e-30)
    H = ugdt + ugst - (ugdt * ugst) / S;
  else
    H = 0.0;
  q_C = COX * ((2.0 / 3.0) * H + GAMMA * sqrt(PHI - ubs));

  /* Smooth blending at B/C boundary */
  sigma_BC = sigmoid_s(ugs - vte, EPS_SMOOTH);

  if (ugb <= vfb)
    return q_A;
  else
    return (1.0 - sigma_BC) * q_B + sigma_BC * q_C;
}

static sunrealtype fn_qsrc_smooth(sunrealtype vgb, sunrealtype vgs,
                                  sunrealtype vgd)
{
  sunrealtype ugs, ugd, ugb, ubs, vfb, vte;
  sunrealtype ugst, ugdt, S, H, q_C, sigma_BC;

  if ((vgs - vgd) <= 0.0) { ugs = vgd; ugd = vgs; }
  else { ugs = vgs; ugd = vgd; }
  ugb = vgb;
  ubs = ugs - ugb;

  vfb = VT0 - GAMMA * sqrt(PHI) - PHI;
  vte = VT0 + GAMMA * (sqrt(PHI - ubs) - sqrt(PHI));

  /* Region C (smoothed) */
  ugst = softplus(ugs - vte, EPS_SMOOTH);
  ugdt = softplus(ugd - vte, EPS_SMOOTH);
  S = ugst + ugdt;
  if (S > 1e-30)
    H = ugdt + ugst - (ugdt * ugst) / S;
  else
    H = 0.0;
  q_C = -COX * (1.0 / 3.0) * H;

  /* Smooth blending: q = sigma * q_C (regions A,B have q=0) */
  sigma_BC = sigmoid_s(ugs - vte, EPS_SMOOTH);
  return sigma_BC * q_C;
}

static sunrealtype fn_qdrain_smooth(sunrealtype vgb, sunrealtype vgs,
                                    sunrealtype vgd)
{
  /* qdrain is identical to qsrc in this model */
  return fn_qsrc_smooth(vgb, vgs, vgd);
}

/* =========================================================================
 * Analytic Jacobian of smoothed charge functions
 * ========================================================================= */

static void fn_dqgate_smooth(sunrealtype vgb, sunrealtype vgs, sunrealtype vgd,
                             sunrealtype *dq_dvgb, sunrealtype *dq_dvgs,
                             sunrealtype *dq_dvgd)
{
  sunrealtype ugs, ugd, ugb, ubs, vfb, vte, sqrt_phi_ubs, dvte_dubs;
  sunrealtype arg_B, q_B, q_C;
  sunrealtype ugst, ugdt, S, H, dsp_s, dsp_d;
  sunrealtype dH_dugst, dH_dugdt;
  sunrealtype dugst_dugs, dugst_dugb, dugdt_dugd, dugdt_dugs, dugdt_dugb;
  sunrealtype dsqrt_dugs, dsqrt_dugb;
  sunrealtype dqC_dugs, dqC_dugb, dqC_dugd;
  sunrealtype dqB_dugb;
  sunrealtype sigma_BC, dsigma_BC, dsigma_dugs, dsigma_dugb, dsigma_dugd;
  sunrealtype dq_dugb, dq_dugs, dq_dugd;
  int swapped;

  swapped = (vgs - vgd) <= 0.0;
  if (swapped) { ugs = vgd; ugd = vgs; }
  else { ugs = vgs; ugd = vgd; }
  ugb = vgb;
  ubs = ugs - ugb;

  vfb = VT0 - GAMMA * sqrt(PHI) - PHI;
  sqrt_phi_ubs = sqrt(PHI - ubs);
  vte = VT0 + GAMMA * (sqrt_phi_ubs - sqrt(PHI));
  dvte_dubs = -GAMMA / (2.0 * sqrt_phi_ubs);

  /* Region B derivatives */
  arg_B = (GAMMA / 2.0) * (GAMMA / 2.0) + ugb - vfb;
  if (arg_B < 1e-30) arg_B = 1e-30;
  dqB_dugb = COX * GAMMA / (2.0 * sqrt(arg_B));

  /* Region C derivatives (smoothed ugst, ugdt) */
  ugst = softplus(ugs - vte, EPS_SMOOTH);
  dsp_s = dsoftplus(ugs - vte, EPS_SMOOTH);
  ugdt = softplus(ugd - vte, EPS_SMOOTH);
  dsp_d = dsoftplus(ugd - vte, EPS_SMOOTH);

  S = ugst + ugdt;
  if (S > 1e-30) {
    H = ugdt + ugst - (ugdt * ugst) / S;
    dH_dugst = 1.0 - (ugdt * ugdt) / (S * S);
    dH_dugdt = 1.0 - (ugst * ugst) / (S * S);
  } else {
    H = 0.0;
    dH_dugst = 1.0;
    dH_dugdt = 0.0;
  }

  /* Chain rules: d(ugs-vte)/dugs = 1 - dvte_dubs, d(ugs-vte)/dugb = dvte_dubs */
  dugst_dugs = dsp_s * (1.0 - dvte_dubs);
  dugst_dugb = dsp_s * dvte_dubs;
  dugdt_dugd = dsp_d * 1.0;
  dugdt_dugs = dsp_d * (-dvte_dubs);
  dugdt_dugb = dsp_d * dvte_dubs;

  dsqrt_dugs = -GAMMA / (2.0 * sqrt_phi_ubs);
  dsqrt_dugb =  GAMMA / (2.0 * sqrt_phi_ubs);

  dqC_dugs = COX * ((2.0/3.0)*(dH_dugst*dugst_dugs + dH_dugdt*dugdt_dugs) + dsqrt_dugs);
  dqC_dugb = COX * ((2.0/3.0)*(dH_dugst*dugst_dugb + dH_dugdt*dugdt_dugb) + dsqrt_dugb);
  dqC_dugd = COX * (2.0/3.0) * dH_dugdt * dugdt_dugd;

  /* Sigmoid blending derivatives */
  sigma_BC = sigmoid_s(ugs - vte, EPS_SMOOTH);
  dsigma_BC = dsigmoid_s(ugs - vte, EPS_SMOOTH);
  dsigma_dugs = dsigma_BC * (1.0 - dvte_dubs);
  dsigma_dugb = dsigma_BC * dvte_dubs;
  dsigma_dugd = 0.0;

  if (ugb <= vfb) {
    /* Region A (hard boundary, already C1) */
    dq_dugb = COX;
    dq_dugs = 0.0;
    dq_dugd = 0.0;
  } else {
    /* Blended: q = (1-sigma)*q_B + sigma*q_C */
    q_B = COX * GAMMA * (sqrt(arg_B) - GAMMA / 2.0);
    q_C = COX * ((2.0/3.0)*H + GAMMA*sqrt_phi_ubs);

    dq_dugb = dsigma_dugb*(q_C - q_B) + (1.0-sigma_BC)*dqB_dugb + sigma_BC*dqC_dugb;
    dq_dugs = dsigma_dugs*(q_C - q_B) + sigma_BC*dqC_dugs;
    dq_dugd = dsigma_dugd*(q_C - q_B) + sigma_BC*dqC_dugd;
  }

  /* Map back from (ugb, ugs, ugd) to (vgb, vgs, vgd) */
  *dq_dvgb = dq_dugb;
  if (swapped) { *dq_dvgs = dq_dugd; *dq_dvgd = dq_dugs; }
  else { *dq_dvgs = dq_dugs; *dq_dvgd = dq_dugd; }
}

static void fn_dqsrc_smooth(sunrealtype vgb, sunrealtype vgs, sunrealtype vgd,
                            sunrealtype *dq_dvgb, sunrealtype *dq_dvgs,
                            sunrealtype *dq_dvgd)
{
  sunrealtype ugs, ugd, ugb, ubs, vfb, vte, sqrt_phi_ubs, dvte_dubs;
  sunrealtype ugst, ugdt, S, H, dsp_s, dsp_d;
  sunrealtype dH_dugst, dH_dugdt;
  sunrealtype dugst_dugs, dugst_dugb, dugdt_dugd, dugdt_dugs, dugdt_dugb;
  sunrealtype dqC_dugs, dqC_dugb, dqC_dugd, q_C;
  sunrealtype sigma_BC, dsigma_BC, dsigma_dugs, dsigma_dugb, dsigma_dugd;
  sunrealtype dq_dugb, dq_dugs, dq_dugd;
  int swapped;

  swapped = (vgs - vgd) <= 0.0;
  if (swapped) { ugs = vgd; ugd = vgs; }
  else { ugs = vgs; ugd = vgd; }
  ugb = vgb;
  ubs = ugs - ugb;

  vfb = VT0 - GAMMA * sqrt(PHI) - PHI;
  sqrt_phi_ubs = sqrt(PHI - ubs);
  vte = VT0 + GAMMA * (sqrt_phi_ubs - sqrt(PHI));
  dvte_dubs = -GAMMA / (2.0 * sqrt_phi_ubs);

  /* Region C derivatives */
  ugst = softplus(ugs - vte, EPS_SMOOTH);
  dsp_s = dsoftplus(ugs - vte, EPS_SMOOTH);
  ugdt = softplus(ugd - vte, EPS_SMOOTH);
  dsp_d = dsoftplus(ugd - vte, EPS_SMOOTH);

  S = ugst + ugdt;
  if (S > 1e-30) {
    H = ugdt + ugst - (ugdt * ugst) / S;
    dH_dugst = 1.0 - (ugdt * ugdt) / (S * S);
    dH_dugdt = 1.0 - (ugst * ugst) / (S * S);
  } else {
    H = 0.0;
    dH_dugst = 1.0;
    dH_dugdt = 0.0;
  }

  dugst_dugs = dsp_s * (1.0 - dvte_dubs);
  dugst_dugb = dsp_s * dvte_dubs;
  dugdt_dugd = dsp_d * 1.0;
  dugdt_dugs = dsp_d * (-dvte_dubs);
  dugdt_dugb = dsp_d * dvte_dubs;

  dqC_dugs = -COX * (1.0/3.0) * (dH_dugst*dugst_dugs + dH_dugdt*dugdt_dugs);
  dqC_dugb = -COX * (1.0/3.0) * (dH_dugst*dugst_dugb + dH_dugdt*dugdt_dugb);
  dqC_dugd = -COX * (1.0/3.0) * dH_dugdt * dugdt_dugd;

  q_C = -COX * (1.0/3.0) * H;

  /* Sigmoid blending: q = sigma * q_C */
  sigma_BC = sigmoid_s(ugs - vte, EPS_SMOOTH);
  dsigma_BC = dsigmoid_s(ugs - vte, EPS_SMOOTH);
  dsigma_dugs = dsigma_BC * (1.0 - dvte_dubs);
  dsigma_dugb = dsigma_BC * dvte_dubs;
  dsigma_dugd = 0.0;

  /* d/dx [sigma * q_C] = dsigma * q_C + sigma * dqC */
  dq_dugb = dsigma_dugb * q_C + sigma_BC * dqC_dugb;
  dq_dugs = dsigma_dugs * q_C + sigma_BC * dqC_dugs;
  dq_dugd = dsigma_dugd * q_C + sigma_BC * dqC_dugd;

  /* Map back */
  *dq_dvgb = dq_dugb;
  if (swapped) { *dq_dvgs = dq_dugd; *dq_dvgd = dq_dugs; }
  else { *dq_dvgs = dq_dugs; *dq_dvgd = dq_dugd; }
}

static void fn_dqdrain_smooth(sunrealtype vgb, sunrealtype vgs, sunrealtype vgd,
                              sunrealtype *dq_dvgb, sunrealtype *dq_dvgs,
                              sunrealtype *dq_dvgd)
{
  fn_dqsrc_smooth(vgb, vgs, vgd, dq_dvgb, dq_dvgs, dq_dvgd);
}

/* =========================================================================
 * Input voltage (piecewise linear, unchanged)
 * ========================================================================= */

static sunrealtype vin(sunrealtype t)
{
  sunrealtype dummy = fmod(t, DELTAT);
  if (dummy < T1_PLS) return 0.0;
  else if (dummy < T2_PLS) return (dummy - T1_PLS) * 0.10e+9 * VHIGH;
  else if (dummy < T3_PLS) return VHIGH;
  else return (DELTAT - dummy) * 0.10e+9 * VHIGH;
}

/* =========================================================================
 * RHS function
 * ========================================================================= */

static int rhs_pump(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)ud;
  sunrealtype* v = N_VGetArrayPointer(y);
  sunrealtype* f = N_VGetArrayPointer(yd);
  sunrealtype vgb = v[5], vgs = v[5] - v[6], vgd = v[5] - v[7];

  f[0] = -v[8];
  f[1] = 0.0;
  f[2] = 0.0;
  f[3] = -v[5] + vin(t);
  f[4] = v[0] - fn_qgate_smooth(vgb, vgs, vgd);
  f[5] = v[1] - CAPS * v[6];
  f[6] = v[2] - fn_qsrc_smooth(vgb, vgs, vgd);
  f[7] = v[3] - CAPD * v[7];
  f[8] = v[4] - fn_qdrain_smooth(vgb, vgs, vgd);
  return 0;
}

/* =========================================================================
 * Analytic Jacobian
 * ========================================================================= */

static int jac_pump(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                    void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;
  sunrealtype* v = N_VGetArrayPointer(y);
  sunrealtype vgb = v[5], vgs = v[5] - v[6], vgd = v[5] - v[7];
  sunrealtype dq_dvgb, dq_dvgs, dq_dvgd;

  SUNMatZero(J);

  /* Trivial entries */
  SM_ELEMENT_D(J, 0, 8) = -1.0;           /* df[0]/dy[8] */
  SM_ELEMENT_D(J, 3, 5) = -1.0;           /* df[3]/dy[5] */
  SM_ELEMENT_D(J, 4, 0) =  1.0;           /* df[4]/dy[0] */
  SM_ELEMENT_D(J, 5, 1) =  1.0;           /* df[5]/dy[1] */
  SM_ELEMENT_D(J, 5, 6) = -CAPS;          /* df[5]/dy[6] */
  SM_ELEMENT_D(J, 6, 2) =  1.0;           /* df[6]/dy[2] */
  SM_ELEMENT_D(J, 7, 3) =  1.0;           /* df[7]/dy[3] */
  SM_ELEMENT_D(J, 7, 7) = -CAPD;          /* df[7]/dy[7] */
  SM_ELEMENT_D(J, 8, 4) =  1.0;           /* df[8]/dy[4] */

  /* Row 4: f[4] = v[0] - qgate(vgb, vgs, vgd) */
  /* dy[5]: vgb+=eps, vgs+=eps, vgd+=eps => dq/dy5 = dq_dvgb + dq_dvgs + dq_dvgd */
  /* dy[6]: vgs-=eps => dq/dy6 = -dq_dvgs */
  /* dy[7]: vgd-=eps => dq/dy7 = -dq_dvgd */
  fn_dqgate_smooth(vgb, vgs, vgd, &dq_dvgb, &dq_dvgs, &dq_dvgd);
  SM_ELEMENT_D(J, 4, 5) = -(dq_dvgb + dq_dvgs + dq_dvgd);
  SM_ELEMENT_D(J, 4, 6) = dq_dvgs;
  SM_ELEMENT_D(J, 4, 7) = dq_dvgd;

  /* Row 6: f[6] = v[2] - qsrc(vgb, vgs, vgd) */
  fn_dqsrc_smooth(vgb, vgs, vgd, &dq_dvgb, &dq_dvgs, &dq_dvgd);
  SM_ELEMENT_D(J, 6, 5) = -(dq_dvgb + dq_dvgs + dq_dvgd);
  SM_ELEMENT_D(J, 6, 6) = dq_dvgs;
  SM_ELEMENT_D(J, 6, 7) = dq_dvgd;

  /* Row 8: f[8] = v[4] - qdrain(vgb, vgs, vgd) */
  fn_dqdrain_smooth(vgb, vgs, vgd, &dq_dvgb, &dq_dvgs, &dq_dvgd);
  SM_ELEMENT_D(J, 8, 5) = -(dq_dvgb + dq_dvgs + dq_dvgd);
  SM_ELEMENT_D(J, 8, 6) = dq_dvgs;
  SM_ELEMENT_D(J, 8, 7) = dq_dvgd;

  return 0;
}

/* =========================================================================
 * Mass matrix (constant, singular)
 * ========================================================================= */

static int mas_pump(sunrealtype t, SUNMatrix M, void* ud,
                    N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)ud; (void)t1; (void)t2; (void)t3;
  SUNMatZero(M);
  SM_ELEMENT_D(M, 0, 0) = 1.0;
  SM_ELEMENT_D(M, 1, 1) = 1.0;
  SM_ELEMENT_D(M, 1, 2) = 1.0;
  SM_ELEMENT_D(M, 2, 3) = 1.0;
  SM_ELEMENT_D(M, 2, 4) = 1.0;
  return 0;
}

/* =========================================================================
 * Main program
 * ========================================================================= */

int main(int argc, char* argv[])
{
  sunrealtype rtol_val = 1e-7, atol_val = 1e-7, h0 = 1e-3;
  int use_schur = 0, nsmin = 3, nsmax = 13;
  if (argc > 1) rtol_val  = atof(argv[1]);
  if (argc > 2) atol_val  = atof(argv[2]);
  if (argc > 3) h0        = atof(argv[3]);
  if (argc > 4) use_schur = atoi(argv[4]);
  if (argc > 5) nsmin     = atoi(argv[5]);
  if (argc > 6) nsmax     = atoi(argv[6]);

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);
  int n = 9;
  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, nsmin, nsmax);

  /* Initial conditions */
  N_Vector y0 = N_VNew_Serial(n, sunctx);
  sunrealtype* y0d = N_VGetArrayPointer(y0);
  y0d[0] = fn_qgate_smooth(0.0, 0.0, 0.0);
  y0d[1] = 0.0;
  y0d[2] = fn_qsrc_smooth(0.0, 0.0, 0.0);
  y0d[3] = 0.0;
  y0d[4] = fn_qdrain_smooth(0.0, 0.0, 0.0);
  for (int i = 5; i < 9; i++) y0d[i] = 0.0;

  Radau5Init(mem, rhs_pump, 0.0, y0);

  /* Linear solver and Jacobian */
  SUNMatrix Jt = SUNDenseMatrix(n, n, sunctx);
  Radau5SetLinearSolver(mem, Jt, NULL);
  Radau5SetSchurDecomp(mem, use_schur);
  Radau5SetJacFn(mem, jac_pump);

  /* Mass matrix */
  SUNMatrix Mt = SUNDenseMatrix(n, n, sunctx);
  Radau5SetMassFn(mem, mas_pump, Mt);
  Radau5SetDAEIndex(mem, 8, 1, 0);

  /* Component-wise tolerances */
  N_Vector rtol_v = N_VNew_Serial(n, sunctx);
  N_Vector atol_v = N_VNew_Serial(n, sunctx);
  sunrealtype* rtol_d = N_VGetArrayPointer(rtol_v);
  sunrealtype* atol_d = N_VGetArrayPointer(atol_v);
  for (int i = 0; i < 5; i++) { rtol_d[i] = rtol_val; atol_d[i] = atol_val*1e-6; }
  for (int i = 5; i < 8; i++) { rtol_d[i] = rtol_val; atol_d[i] = atol_val; }
  rtol_d[8] = rtol_val; atol_d[8] = atol_val;
  Radau5SVtolerances(mem, rtol_v, atol_v);

  Radau5SetInitStep(mem, h0);
  Radau5SetMaxNumSteps(mem, 1000000);

  N_Vector yout = N_VNew_Serial(n, sunctx);
  sunrealtype tret;

  /* Segmented integration at Vin discontinuity times */
  int ndisc = 39;
  sunrealtype disc[41];
  disc[0] = 0.0;
  disc[1] = 50.0e-9;
  disc[2] = 60.0e-9;
  disc[3] = 110.0e-9;
  for (int i = 4; i <= 39; i++)
    disc[i] = disc[i % 4] + (sunrealtype)(i / 4) * 120.0e-9;
  disc[40] = 1200.0e-9;

  printf("Charge pump (smoothed, eps=%.1e) nsmin=%d nsmax=%d %s\n",
         EPS_SMOOTH, nsmin, nsmax, use_schur ? "Schur" : "Eigen");

  int ret = 0;
  for (int seg = 0; seg <= ndisc; seg++) {
    ret = Radau5Solve(mem, disc[seg+1], yout, &tret);
    if (seg <= 2 || ret < 0) {
      sunrealtype* yd = N_VGetArrayPointer(yout);
      printf("[seg %2d] t_end=%.4e ret=%d tret=%.6e y[5]=%.6e y[8]=%.6e\n",
             seg, disc[seg+1], ret, tret, yd[5], yd[8]);
    }
    if (ret < 0) break;
    if (seg < ndisc) Radau5ResetForDiscontinuity(mem, 1e-12);
  }

  /* Final results */
  printf("\nFinal: ret=%d\n", ret);
  int pass = 0;
  if (ret >= 0) {
    sunrealtype* yd = N_VGetArrayPointer(yout);
    sunrealtype ref0 = 1.26280042987676e-13;
    sunrealtype ref8 = 1.52255686815578e-04;
    for (int i = 0; i < 9; i++)
      printf("  y[%d] = %.14e\n", i, yd[i]);
    printf("  ref y[0] = %.14e\n", ref0);
    printf("  ref y[8] = %.14e\n", ref8);
    sunrealtype err0 = fabs(yd[0] - ref0) / fabs(ref0);
    sunrealtype err8 = fabs(yd[8] - ref8) / fabs(ref8);
    printf("  rel_err y[0] = %.3e\n", err0);
    printf("  rel_err y[8] = %.3e\n", err8);
    sunrealtype maxerr = err0 > err8 ? err0 : err8;
    pass = (maxerr < 1e-2);
    printf("  %s (max_rel_err = %.3e)\n", pass ? "PASSED" : "FAILED", maxerr);
  }

  /* Statistics */
  long int nstep, naccpt, nfcn, njac, ndec;
  Radau5GetNumSteps(mem, &nstep);
  Radau5GetNumAccSteps(mem, &naccpt);
  Radau5GetNumRhsEvals(mem, &nfcn);
  Radau5GetNumJacEvals(mem, &njac);
  Radau5GetNumDecomps(mem, &ndec);
  printf("nstep=%ld naccpt=%ld nfcn=%ld njac=%ld ndec=%ld\n",
         nstep, naccpt, nfcn, njac, ndec);

  /* Cleanup */
  N_VDestroy(y0); N_VDestroy(yout);
  N_VDestroy(rtol_v); N_VDestroy(atol_v);
  SUNMatDestroy(Jt); SUNMatDestroy(Mt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
