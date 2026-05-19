/* ---------------------------------------------------------------------------
 * radau5_pump.c — Charge Pump (DAE index-2, n=9)
 *
 * Faithful C translation of the Fortran source:
 *   IVPtestset_2.4/src/problems/pump.f
 *
 * M*y' = f(t,y) where M is a 9×9 singular mass matrix (5 nonzeros)
 *
 * MOS transistor charge pump circuit simulation.
 * State vector (0-based):
 *   y[0..4] = charges (qgate, caps*V_s, qsrc, capd*V_d, qdrain)
 *   y[5..7] = node voltages (Vg, Vs, Vd)
 *   y[8]    = current
 *
 * DAE index: ind(1..8)=1, ind(9)=2 → nind1=8, nind2=1
 *
 * t in [0, 1200e-9] with 39 discontinuities (pulsed input vin(t))
 * Reference solution from IVPtestset (GAMD quadruple precision)
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_types.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

/* MOS transistor parameters */
#define VT0   0.20
#define GAMMA 0.035
#define PHI   1.01
#define COX   4.0e-12
#define CAPD  0.40e-12
#define CAPS  1.60e-12

/* Input voltage parameters */
#define VHIGH  20.0
#define DELTAT 120.0
#define T1_PLS 50.0
#define T2_PLS 60.0
#define T3_PLS 110.0

/* =========================================================================
 * Pulsed input voltage (trapezoidal, period 120ns)
 * =========================================================================*/
static sunrealtype vin(sunrealtype t)
{
  sunrealtype dummy = fmod(t, DELTAT);
  if (dummy < T1_PLS)
    return 0.0;
  else if (dummy < T2_PLS)
    return (dummy - T1_PLS) * 0.10 * VHIGH;
  else if (dummy < T3_PLS)
    return VHIGH;
  else
    return (DELTAT - dummy) * 0.10 * VHIGH;
}

/* =========================================================================
 * Charge functions (MOS transistor model)
 * =========================================================================*/
static sunrealtype fn_qgate(sunrealtype vgb, sunrealtype vgs, sunrealtype vgd)
{
  sunrealtype ugs, ugd, ugb, ubs, vfb, vte, ugst, ugdt;

  if ((vgs - vgd) <= 0.0) {
    ugs = vgd;
    ugd = vgs;
  } else {
    ugs = vgs;
    ugd = vgd;
  }

  ugb = vgb;
  ubs = ugs - ugb;

  vfb = VT0 - GAMMA * sqrt(PHI) - PHI;
  vte = VT0 + GAMMA * (sqrt(PHI - ubs) - sqrt(PHI));

  if (ugb <= vfb) {
    return COX * (ugb - vfb);
  } else if (ugb > vfb && ugs <= vte) {
    return COX * GAMMA *
           (sqrt((GAMMA / 2.0) * (GAMMA / 2.0) + ugb - vfb) - GAMMA / 2.0);
  } else {
    ugst = ugs - vte;
    ugdt = (ugd > vte) ? ugd - vte : 0.0;
    return COX * ((2.0 / 3.0) * (ugdt + ugst -
           ((ugdt * ugst) / (ugdt + ugst))) +
           GAMMA * sqrt(PHI - ubs));
  }
}

static sunrealtype fn_qsrc(sunrealtype vgb, sunrealtype vgs, sunrealtype vgd)
{
  sunrealtype ugs, ugd, ugb, ubs, vfb, vte, ugst, ugdt;

  if ((vgs - vgd) <= 0.0) {
    ugs = vgd;
    ugd = vgs;
  } else {
    ugs = vgs;
    ugd = vgd;
  }

  ugb = vgb;
  ubs = ugs - ugb;

  vfb = VT0 - GAMMA * sqrt(PHI) - PHI;
  vte = VT0 + GAMMA * (sqrt(PHI - ubs) - sqrt(PHI));

  if (ugb <= vfb) {
    return 0.0;
  } else if (ugb > vfb && ugs <= vte) {
    return 0.0;
  } else {
    ugst = ugs - vte;
    ugdt = (ugd >= vte) ? ugd - vte : 0.0;
    return -COX * (1.0 / 3.0) * (ugdt + ugst -
           ((ugdt * ugst) / (ugdt + ugst)));
  }
}

static sunrealtype fn_qdrain(sunrealtype vgb, sunrealtype vgs, sunrealtype vgd)
{
  sunrealtype ugs, ugd, ugb, ubs, vfb, vte, ugst, ugdt;

  if ((vgs - vgd) <= 0.0) {
    ugs = vgd;
    ugd = vgs;
  } else {
    ugs = vgs;
    ugd = vgd;
  }

  ugb = vgb;
  ubs = ugs - ugb;

  vfb = VT0 - GAMMA * sqrt(PHI) - PHI;
  vte = VT0 + GAMMA * (sqrt(PHI - ubs) - sqrt(PHI));

  if (ugb <= vfb) {
    return 0.0;
  } else if (ugb > vfb && ugs <= vte) {
    return 0.0;
  } else {
    ugst = ugs - vte;
    ugdt = (ugd >= vte) ? ugd - vte : 0.0;
    return -COX * (1.0 / 3.0) * (ugdt + ugst -
           ((ugdt * ugst) / (ugdt + ugst)));
  }
}

/* =========================================================================
 * Derivatives of charge functions w.r.t. (vgb, vgs, vgd)
 *
 * Each charge function Q(vgb, vgs, vgd) has:
 *   Step 1: direction swap based on sign(vgs - vgd)
 *   Step 2: three operating regions (accumulation, depletion, inversion)
 *
 * These derivative functions compute ∂Q/∂vgb, ∂Q/∂vgs, ∂Q/∂vgd analytically.
 * =========================================================================*/

static void fn_dqgate(sunrealtype vgb, sunrealtype vgs, sunrealtype vgd,
                      sunrealtype* dq_dvgb, sunrealtype* dq_dvgs, sunrealtype* dq_dvgd)
{
  sunrealtype ugs, ugd, ugb, ubs, vfb, vte;
  int swapped = ((vgs - vgd) <= 0.0);

  if (swapped) { ugs = vgd; ugd = vgs; }
  else         { ugs = vgs; ugd = vgd; }

  ugb = vgb;
  ubs = ugs - ugb;

  vfb = VT0 - GAMMA * sqrt(PHI) - PHI;
  vte = VT0 + GAMMA * (sqrt(PHI - ubs) - sqrt(PHI));

  /* Derivatives w.r.t. intermediate variables (ugs, ugd, ugb) */
  sunrealtype dq_dugb = 0.0, dq_dugs = 0.0, dq_dugd = 0.0;

  if (ugb <= vfb)
  {
    /* Region A: qgate = COX * (ugb - vfb) */
    dq_dugb = COX;
    dq_dugs = 0.0;
    dq_dugd = 0.0;
  }
  else if (ugs <= vte)
  {
    /* Region B: qgate = COX*GAMMA*(sqrt(G^2/4 + ugb - vfb) - G/2) */
    sunrealtype arg = (GAMMA/2.0)*(GAMMA/2.0) + ugb - vfb;
    dq_dugb = COX * GAMMA / (2.0 * sqrt(arg));
    dq_dugs = 0.0;
    dq_dugd = 0.0;
  }
  else
  {
    /* Region C (inversion):
     * qgate = COX * ((2/3)*H(ugst, ugdt) + GAMMA*sqrt(PHI - ubs))
     * where H = ugst + ugdt - ugst*ugdt/(ugst+ugdt)
     *       ugst = ugs - vte,  ugdt = max(ugd - vte, 0)
     *       vte = VT0 + GAMMA*(sqrt(PHI-ubs) - sqrt(PHI))
     *       ubs = ugs - ugb
     */
    sunrealtype sqrt_phi_ubs = sqrt(PHI - ubs);
    sunrealtype dvte_dubs = -GAMMA / (2.0 * sqrt_phi_ubs);

    sunrealtype ugst = ugs - vte;
    int ugd_above_vte = (ugd > vte);
    sunrealtype ugdt = ugd_above_vte ? (ugd - vte) : 0.0;

    sunrealtype S = ugst + ugdt;
    /* ∂H/∂ugst = 1 - ugdt^2/S^2,  ∂H/∂ugdt = 1 - ugst^2/S^2 */
    sunrealtype dH_dugst = (S > 0.0) ? 1.0 - (ugdt*ugdt)/(S*S) : 1.0;
    sunrealtype dH_dugdt = (S > 0.0) ? 1.0 - (ugst*ugst)/(S*S) : 0.0;

    /* dugst/dugs = 1 - dvte/dubs * (dubs/dugs=1) = 1 - dvte_dubs */
    sunrealtype dugst_dugs = 1.0 - dvte_dubs;
    /* dugst/dugb = 0 - dvte/dubs * (dubs/dugb=-1) = dvte_dubs */
    sunrealtype dugst_dugb = dvte_dubs;

    /* dugdt/dugd: if ugd > vte, dugdt/dugd = 1; else 0 */
    sunrealtype dugdt_dugd = ugd_above_vte ? 1.0 : 0.0;
    /* dugdt/dugs: if ugd > vte, dugdt/dugs = -dvte/dubs; else 0 */
    sunrealtype dugdt_dugs = ugd_above_vte ? (-dvte_dubs) : 0.0;
    /* dugdt/dugb: if ugd > vte, dugdt/dugb = dvte_dubs; else 0 */
    sunrealtype dugdt_dugb = ugd_above_vte ? dvte_dubs : 0.0;

    /* d(GAMMA*sqrt(PHI-ubs))/dugs = GAMMA * (-1)/(2*sqrt(PHI-ubs)) */
    sunrealtype dsqrt_dugs = -GAMMA / (2.0 * sqrt_phi_ubs);
    /* d(GAMMA*sqrt(PHI-ubs))/dugb = GAMMA * (1)/(2*sqrt(PHI-ubs)) */
    sunrealtype dsqrt_dugb = GAMMA / (2.0 * sqrt_phi_ubs);

    /* dqgate/dugs */
    dq_dugs = COX * ((2.0/3.0) * (dH_dugst * dugst_dugs + dH_dugdt * dugdt_dugs)
                     + dsqrt_dugs);
    /* dqgate/dugb */
    dq_dugb = COX * ((2.0/3.0) * (dH_dugst * dugst_dugb + dH_dugdt * dugdt_dugb)
                     + dsqrt_dugb);
    /* dqgate/dugd */
    dq_dugd = COX * (2.0/3.0) * dH_dugdt * dugdt_dugd;
  }

  /* Map from (ugb, ugs, ugd) derivatives back to (vgb, vgs, vgd) */
  /* ugb = vgb always, so dq/dvgb = dq/dugb */
  *dq_dvgb = dq_dugb;
  if (swapped)
  {
    /* ugs = vgd, ugd = vgs */
    *dq_dvgs = dq_dugd;
    *dq_dvgd = dq_dugs;
  }
  else
  {
    /* ugs = vgs, ugd = vgd */
    *dq_dvgs = dq_dugs;
    *dq_dvgd = dq_dugd;
  }
}

static void fn_dqsrc(sunrealtype vgb, sunrealtype vgs, sunrealtype vgd,
                     sunrealtype* dq_dvgb, sunrealtype* dq_dvgs, sunrealtype* dq_dvgd)
{
  sunrealtype ugs, ugd, ugb, ubs, vfb, vte;
  int swapped = ((vgs - vgd) <= 0.0);

  if (swapped) { ugs = vgd; ugd = vgs; }
  else         { ugs = vgs; ugd = vgd; }

  ugb = vgb;
  ubs = ugs - ugb;

  vfb = VT0 - GAMMA * sqrt(PHI) - PHI;
  vte = VT0 + GAMMA * (sqrt(PHI - ubs) - sqrt(PHI));

  if (ugb <= vfb || ugs <= vte)
  {
    /* Regions A and B: qsrc = 0 */
    *dq_dvgb = 0.0; *dq_dvgs = 0.0; *dq_dvgd = 0.0;
    return;
  }

  /* Region C: qsrc = -COX*(1/3)*(ugst + ugdt - ugst*ugdt/(ugst+ugdt)) */
  sunrealtype sqrt_phi_ubs = sqrt(PHI - ubs);
  sunrealtype dvte_dubs = -GAMMA / (2.0 * sqrt_phi_ubs);

  sunrealtype ugst = ugs - vte;
  int ugd_above_vte = (ugd >= vte);
  sunrealtype ugdt = ugd_above_vte ? (ugd - vte) : 0.0;

  sunrealtype S = ugst + ugdt;
  sunrealtype dH_dugst = (S > 0.0) ? 1.0 - (ugdt*ugdt)/(S*S) : 1.0;
  sunrealtype dH_dugdt = (S > 0.0) ? 1.0 - (ugst*ugst)/(S*S) : 0.0;

  sunrealtype dugst_dugs = 1.0 - dvte_dubs;
  sunrealtype dugst_dugb = dvte_dubs;
  sunrealtype dugdt_dugd = ugd_above_vte ? 1.0 : 0.0;
  sunrealtype dugdt_dugs = ugd_above_vte ? (-dvte_dubs) : 0.0;
  sunrealtype dugdt_dugb = ugd_above_vte ? dvte_dubs : 0.0;

  sunrealtype dq_dugs = -COX * (1.0/3.0) * (dH_dugst * dugst_dugs + dH_dugdt * dugdt_dugs);
  sunrealtype dq_dugb = -COX * (1.0/3.0) * (dH_dugst * dugst_dugb + dH_dugdt * dugdt_dugb);
  sunrealtype dq_dugd = -COX * (1.0/3.0) * dH_dugdt * dugdt_dugd;

  *dq_dvgb = dq_dugb;
  if (swapped) { *dq_dvgs = dq_dugd; *dq_dvgd = dq_dugs; }
  else         { *dq_dvgs = dq_dugs; *dq_dvgd = dq_dugd; }
}

/* fn_qdrain is identical to fn_qsrc in this model */
static void fn_dqdrain(sunrealtype vgb, sunrealtype vgs, sunrealtype vgd,
                       sunrealtype* dq_dvgb, sunrealtype* dq_dvgs, sunrealtype* dq_dvgd)
{
  fn_dqsrc(vgb, vgs, vgd, dq_dvgb, dq_dvgs, dq_dvgd);
}

/* =========================================================================
 * Analytic Jacobian: J[i][j] = ∂f[i]/∂y[j]
 *
 * Chain rule for charge functions called as Q(y5, y5-y6, y5-y7):
 *   ∂Q/∂y5 = ∂Q/∂vgb + ∂Q/∂vgs + ∂Q/∂vgd
 *   ∂Q/∂y6 = -∂Q/∂vgs
 *   ∂Q/∂y7 = -∂Q/∂vgd
 * =========================================================================*/
static int jac_pump(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                    void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;
  sunrealtype* v = N_VGetArrayPointer(y);

  SUNMatZero(J);

  /* Row 0: f[0] = -1e-9 * y[8] */
  SM_ELEMENT_D(J, 0, 8) = -1.0e-9;

  /* Row 1: f[1] = 0 — all zero */
  /* Row 2: f[2] = 0 — all zero */

  /* Row 3: f[3] = 1e-9 * (-y[5] + vin(t)) */
  SM_ELEMENT_D(J, 3, 5) = -1.0e-9;

  /* Row 4: f[4] = 1e-9 * (y[0] - qgate(y5, y5-y6, y5-y7)) */
  SM_ELEMENT_D(J, 4, 0) = 1.0e-9;
  {
    sunrealtype dq_dvgb, dq_dvgs, dq_dvgd;
    fn_dqgate(v[5], v[5]-v[6], v[5]-v[7], &dq_dvgb, &dq_dvgs, &dq_dvgd);
    SM_ELEMENT_D(J, 4, 5) = -1.0e-9 * (dq_dvgb + dq_dvgs + dq_dvgd);
    SM_ELEMENT_D(J, 4, 6) = -1.0e-9 * (-dq_dvgs);
    SM_ELEMENT_D(J, 4, 7) = -1.0e-9 * (-dq_dvgd);
  }

  /* Row 5: f[5] = 1e-9 * (y[1] - CAPS * y[6]) */
  SM_ELEMENT_D(J, 5, 1) = 1.0e-9;
  SM_ELEMENT_D(J, 5, 6) = -1.0e-9 * CAPS;

  /* Row 6: f[6] = 1e-9 * (y[2] - qsrc(y5, y5-y6, y5-y7)) */
  SM_ELEMENT_D(J, 6, 2) = 1.0e-9;
  {
    sunrealtype dq_dvgb, dq_dvgs, dq_dvgd;
    fn_dqsrc(v[5], v[5]-v[6], v[5]-v[7], &dq_dvgb, &dq_dvgs, &dq_dvgd);
    SM_ELEMENT_D(J, 6, 5) = -1.0e-9 * (dq_dvgb + dq_dvgs + dq_dvgd);
    SM_ELEMENT_D(J, 6, 6) = -1.0e-9 * (-dq_dvgs);
    SM_ELEMENT_D(J, 6, 7) = -1.0e-9 * (-dq_dvgd);
  }

  /* Row 7: f[7] = 1e-9 * (y[3] - CAPD * y[7]) */
  SM_ELEMENT_D(J, 7, 3) = 1.0e-9;
  SM_ELEMENT_D(J, 7, 7) = -1.0e-9 * CAPD;

  /* Row 8: f[8] = 1e-9 * (y[4] - qdrain(y5, y5-y6, y5-y7)) */
  SM_ELEMENT_D(J, 8, 4) = 1.0e-9;
  {
    sunrealtype dq_dvgb, dq_dvgs, dq_dvgd;
    fn_dqdrain(v[5], v[5]-v[6], v[5]-v[7], &dq_dvgb, &dq_dvgs, &dq_dvgd);
    SM_ELEMENT_D(J, 8, 5) = -1.0e-9 * (dq_dvgb + dq_dvgs + dq_dvgd);
    SM_ELEMENT_D(J, 8, 6) = -1.0e-9 * (-dq_dvgs);
    SM_ELEMENT_D(J, 8, 7) = -1.0e-9 * (-dq_dvgd);
  }

  return 0;
}

/* =========================================================================
 * RHS: M*y' = f(t, y)
 * =========================================================================*/
static int rhs_pump(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)ud;
  sunrealtype* v = N_VGetArrayPointer(y);
  sunrealtype* f = N_VGetArrayPointer(yd);

  f[0] = -1e-9*v[8];
  f[1] = 0.0;
  f[2] = 0.0;
  f[3] = 1e-9*(-v[5] + vin(t));
  f[4] = 1e-9*(v[0] - fn_qgate(v[5], v[5] - v[6], v[5] - v[7]));
  f[5] = 1e-9*(v[1] - CAPS * v[6]);
  f[6] = 1e-9*(v[2] - fn_qsrc(v[5], v[5] - v[6], v[5] - v[7]));
  f[7] = 1e-9*(v[3] - CAPD * v[7]);
  f[8] = 1e-9*(v[4] - fn_qdrain(v[5], v[5] - v[6], v[5] - v[7]));

  return 0;
}

/* =========================================================================
 * Mass matrix (dense 9×9, only 5 nonzero entries)
 * Decoded from Fortran banded storage (mlmas=0, mumas=2):
 *   M[0,0]=1, M[1,1]=1, M[1,2]=1, M[2,3]=1, M[2,4]=1
 * =========================================================================*/
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

int main(int argc, char* argv[])
{
  sunrealtype rtol = 1.0e-7;
  sunrealtype atol_val = 1.0e-7;
  sunrealtype h0   = 0.001;
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

  int n = 9;


  /* Initial conditions (from pump.f init subroutine) */
  N_Vector y0 = N_VNew_Serial(n, sunctx);
  sunrealtype* y0d = N_VGetArrayPointer(y0);
  y0d[0] = fn_qgate(0.0, 0.0, 0.0);
  y0d[1] = 0.0;
  y0d[2] = fn_qsrc(0.0, 0.0, 0.0);
  y0d[3] = 0.0;
  y0d[4] = fn_qdrain(0.0, 0.0, 0.0);
  y0d[5] = 0.0;
  y0d[6] = 0.0;
  y0d[7] = 0.0;
  y0d[8] = 0.0;



  /* Component-dependent tolerances (from pump.f settolerances):
   *   Charges (y[0..4]) are much smaller than voltages/current:
   *     atol(i) = Tol * 1e-6  for i=1..5 (0-based: 0..4)
   *     atol(i) = Tol         for i=6..8 (0-based: 5..7)
   *     rtol(i) = Tol         for i=1..8 (0-based: 0..7)
   *     atol(9) = rtol(9) = Tol (0-based: index 8)
   */
  N_Vector rtol_v = N_VNew_Serial(n, sunctx);
  N_Vector atol_v = N_VNew_Serial(n, sunctx);
  sunrealtype* rtol_d = N_VGetArrayPointer(rtol_v);
  sunrealtype* atol_d = N_VGetArrayPointer(atol_v);
  int i;
  for (i = 0; i < 5; i++) {
    rtol_d[i] = rtol;
    atol_d[i] = atol_val * 1.0e-6;
  }
  for (i = 5; i < 8; i++) {
    rtol_d[i] = rtol;
    atol_d[i] = atol_val;
  }
  rtol_d[8] = rtol;
  atol_d[8] = atol_val;


  N_Vector yout = N_VNew_Serial(n, sunctx);
  sunrealtype tret;
  int ret = RADAU5_SUCCESS;

  /* -----------------------------------------------------------------------
   * Segmented integration: 39 discontinuities from pulsed input vin(t)
   * Discontinuity times from Fortran:
   *   t(1)=50e-9, t(2)=60e-9, t(3)=110e-9
   *   t(i) = t(mod(i,4)) + floor(i/4)*120e-9  for i=4..39
   *   t(40) = 1200e-9 (final time)
   * ----------------------------------------------------------------------- */
  int ndisc = 39;
  sunrealtype disc[41];
  disc[0] = 0.0;
  disc[1] = 50.0;
  disc[2] = 60.0;
  disc[3] = 110.0;
  for (int i = 4; i <= 39; i++)
    disc[i] = disc[i % 4] + (sunrealtype)(i / 4) * 120.0;
  disc[40] = 1200.0;

  for (int seg = 0; seg < ndisc + 1; seg++) {
    sunrealtype t_seg_end = disc[seg + 1];

    void* mem = Radau5Create(sunctx);
    Radau5SetOrderLimits(mem, nsmin, nsmax);

    Radau5Init(mem, rhs_pump, disc[seg], y0);

    SUNMatrix Jt = SUNDenseMatrix(n, n, sunctx);
    Radau5SetLinearSolver(mem, Jt, NULL);
    Radau5SetSchurDecomp(mem, use_schur);
    /* Radau5SetJacFn(mem, jac_pump); — temporarily use DQ to compare */

    /* Mass matrix */
    SUNMatrix Mt = SUNDenseMatrix(n, n, sunctx);
    Radau5SetMassFn(mem, mas_pump, Mt);

    /* DAE index: ind(1..8)=1, ind(9)=2 */
    Radau5SetDAEIndex(mem, 8, 1, 0);
    Radau5SVtolerances(mem, rtol_v, atol_v);

    Radau5SetInitStep(mem, h0);

    ret = Radau5Solve(mem, t_seg_end, yout, &tret);

    long int nstep, naccpt, nrejct, nfcn, njac, ndec;
    Radau5GetNumSteps(mem, &nstep);
    Radau5GetNumAccSteps(mem, &naccpt);
    Radau5GetNumRejSteps(mem, &nrejct);
    Radau5GetNumRhsEvals(mem, &nfcn);
    Radau5GetNumJacEvals(mem, &njac);
    Radau5GetNumDecomps(mem, &ndec);
    printf("nstep=%ld naccpt=%ld nrejct=%ld nfcn=%ld njac=%ld ndec=%ld\n",
           nstep, naccpt, nrejct, nfcn, njac, ndec);
    if (ret < 0) {
      printf("Radau5Solve failed at segment %d, tret=%.6e, ret=%d\n",
             seg, tret, ret);
      break;
    }
    N_VScale(1.0, yout, y0);  /* update initial conditions for next segment */
    SUNMatDestroy(Jt); SUNMatDestroy(Mt);
    Radau5Free(&mem);
  }

  /* Reference solution at t=1200e-9 (GAMD quadruple precision) */
  sunrealtype yref[9];
  yref[0] = 0.126280042987675933170893e-12;
  for (i = 1; i <= 7; i++) yref[i] = 0.0;
  yref[8] = 0.152255686815577679043511e-03;

  sunrealtype* yd = N_VGetArrayPointer(yout);

  printf("=== Charge Pump (rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n",
         rtol, atol_val, h0, use_schur);
  printf("ret=%d, tret=%.6e\n", ret, tret);

  sunrealtype maxrelerr = 0.0;
  for (i = 0; i < n; i++) {
    sunrealtype err = fabs(yd[i] - yref[i]);
    sunrealtype relerr = (fabs(yref[i]) > 1e-15) ? err / fabs(yref[i]) : err;
    if (relerr > maxrelerr) maxrelerr = relerr;
    printf("y[%d] = %22.14e  ref = %22.14e  rel_err = %.3e\n",
           i, yd[i], yref[i], relerr);
  }



  /* Pass criterion: check only y[0] and y[8] (y[1..7] are ~0) */
  sunrealtype rel0 = fabs(yd[0] - yref[0]) / fabs(yref[0]);
  sunrealtype rel8 = fabs(yd[8] - yref[8]) / fabs(yref[8]);
  sunrealtype maxrel = (rel0 > rel8) ? rel0 : rel8;
  printf("max_rel_err (y[0],y[8]) = %.3e\n", maxrel);

  int pass = (ret >= 0 && maxrel < 1.0e-2);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout);
  N_VDestroy(rtol_v); N_VDestroy(atol_v);

  SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
