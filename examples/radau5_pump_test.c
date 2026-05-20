#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_types.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

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

static sunrealtype vin(sunrealtype t) {
  sunrealtype dummy = fmod(t, DELTAT);
  if (dummy < T1_PLS) return 0.0;
  else if (dummy < T2_PLS) return (dummy - T1_PLS) * 0.10e+9 * VHIGH;
  else if (dummy < T3_PLS) return VHIGH;
  else return (DELTAT - dummy) * 0.10e+9 * VHIGH;
}

static sunrealtype fn_qgate(sunrealtype vgb, sunrealtype vgs, sunrealtype vgd) {
  sunrealtype ugs, ugd, ugb, ubs, vfb, vte;
  if ((vgs - vgd) <= 0.0) { ugs = vgd; ugd = vgs; }
  else { ugs = vgs; ugd = vgd; }
  ugb = vgb; ubs = ugs - ugb;
  vfb = VT0 - GAMMA * sqrt(PHI) - PHI;
  vte = VT0 + GAMMA * (sqrt(PHI - ubs) - sqrt(PHI));
  if (ugb <= vfb) return COX * (ugb - vfb);
  else if (ugb > vfb && ugs <= vte)
    return COX * GAMMA * (sqrt((GAMMA/2.0)*(GAMMA/2.0) + ugb - vfb) - GAMMA/2.0);
  else {
    sunrealtype ugst = ugs - vte;
    sunrealtype ugdt = (ugd > vte) ? ugd - vte : 0.0;
    return COX * ((2.0/3.0) * (ugdt + ugst - ((ugdt*ugst)/(ugdt+ugst))) + GAMMA * sqrt(PHI-ubs));
  }
}
static sunrealtype fn_qsrc(sunrealtype vgb, sunrealtype vgs, sunrealtype vgd) {
  sunrealtype ugs, ugd, ugb, ubs, vfb, vte;
  if ((vgs - vgd) <= 0.0) { ugs = vgd; ugd = vgs; }
  else { ugs = vgs; ugd = vgd; }
  ugb = vgb; ubs = ugs - ugb;
  vfb = VT0 - GAMMA * sqrt(PHI) - PHI;
  vte = VT0 + GAMMA * (sqrt(PHI - ubs) - sqrt(PHI));
  if (ugb <= vfb) return 0.0;
  else if (ugb > vfb && ugs <= vte) return 0.0;
  else {
    sunrealtype ugst = ugs - vte;
    sunrealtype ugdt = (ugd >= vte) ? ugd - vte : 0.0;
    return -COX * (1.0/3.0) * (ugdt + ugst - ((ugdt*ugst)/(ugdt+ugst)));
  }
}
static sunrealtype fn_qdrain(sunrealtype vgb, sunrealtype vgs, sunrealtype vgd) {
  sunrealtype ugs, ugd, ugb, ubs, vfb, vte;
  if ((vgs - vgd) <= 0.0) { ugs = vgd; ugd = vgs; }
  else { ugs = vgs; ugd = vgd; }
  ugb = vgb; ubs = ugs - ugb;
  vfb = VT0 - GAMMA * sqrt(PHI) - PHI;
  vte = VT0 + GAMMA * (sqrt(PHI - ubs) - sqrt(PHI));
  if (ugb <= vfb) return 0.0;
  else if (ugb > vfb && ugs <= vte) return 0.0;
  else {
    sunrealtype ugst = ugs - vte;
    sunrealtype ugdt = (ugd >= vte) ? ugd - vte : 0.0;
    return -COX * (1.0/3.0) * (ugdt + ugst - ((ugdt*ugst)/(ugdt+ugst)));
  }
}

static int rhs_pump(sunrealtype t, N_Vector y, N_Vector yd, void* ud) {
  (void)ud;
  sunrealtype* v = N_VGetArrayPointer(y);
  sunrealtype* f = N_VGetArrayPointer(yd);
  f[0] = -v[8];
  f[1] = 0.0;
  f[2] = 0.0;
  f[3] = -v[5] + vin(t);
  f[4] = v[0] - fn_qgate(v[5], v[5]-v[6], v[5]-v[7]);
  f[5] = v[1] - CAPS*v[6];
  f[6] = v[2] - fn_qsrc(v[5], v[5]-v[6], v[5]-v[7]);
  f[7] = v[3] - CAPD*v[7];
  f[8] = v[4] - fn_qdrain(v[5], v[5]-v[6], v[5]-v[7]);
  return 0;
}

/* Analytic Jacobian using centered finite differences for charge derivatives */
static int jac_pump(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                    void* ud, N_Vector t1, N_Vector t2, N_Vector t3) {
  (void)t; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;
  sunrealtype* v = N_VGetArrayPointer(y);
  
  /* Zero the Jacobian */
  SUNMatZero(J);
  
  /* Trivial entries */
  SM_ELEMENT_D(J, 0, 8) = -1.0;           /* df[0]/dy[8] */
  SM_ELEMENT_D(J, 3, 5) = -1.0;           /* df[3]/dy[5] */
  SM_ELEMENT_D(J, 4, 0) = 1.0;            /* df[4]/dy[0] */
  SM_ELEMENT_D(J, 5, 1) = 1.0;            /* df[5]/dy[1] */
  SM_ELEMENT_D(J, 5, 6) = -CAPS;          /* df[5]/dy[6] */
  SM_ELEMENT_D(J, 6, 2) = 1.0;            /* df[6]/dy[2] */
  SM_ELEMENT_D(J, 7, 3) = 1.0;            /* df[7]/dy[3] */
  SM_ELEMENT_D(J, 7, 7) = -CAPD;          /* df[7]/dy[7] */
  SM_ELEMENT_D(J, 8, 4) = 1.0;            /* df[8]/dy[4] */
  
  /* Charge function derivatives w.r.t. voltages y[5], y[6], y[7] */
  /* Use centered finite differences with small increment */
  sunrealtype eps = 1.0e-7;
  sunrealtype vgb = v[5], vgs = v[5] - v[6], vgd = v[5] - v[7];
  
  /* d/dy[5]: vgb += eps, vgs += eps, vgd += eps */
  sunrealtype qg_p = fn_qgate(vgb+eps, vgs+eps, vgd+eps);
  sunrealtype qg_m = fn_qgate(vgb-eps, vgs-eps, vgd-eps);
  sunrealtype qs_p = fn_qsrc(vgb+eps, vgs+eps, vgd+eps);
  sunrealtype qs_m = fn_qsrc(vgb-eps, vgs-eps, vgd-eps);
  sunrealtype qd_p = fn_qdrain(vgb+eps, vgs+eps, vgd+eps);
  sunrealtype qd_m = fn_qdrain(vgb-eps, vgs-eps, vgd-eps);
  SM_ELEMENT_D(J, 4, 5) = -(qg_p - qg_m) / (2.0*eps);
  SM_ELEMENT_D(J, 6, 5) = -(qs_p - qs_m) / (2.0*eps);
  SM_ELEMENT_D(J, 8, 5) = -(qd_p - qd_m) / (2.0*eps);
  
  /* d/dy[6]: vgs -= eps (since vgs = y[5]-y[6]) */
  qg_p = fn_qgate(vgb, vgs-eps, vgd);
  qg_m = fn_qgate(vgb, vgs+eps, vgd);
  qs_p = fn_qsrc(vgb, vgs-eps, vgd);
  qs_m = fn_qsrc(vgb, vgs+eps, vgd);
  qd_p = fn_qdrain(vgb, vgs-eps, vgd);
  qd_m = fn_qdrain(vgb, vgs+eps, vgd);
  SM_ELEMENT_D(J, 4, 6) = -(qg_p - qg_m) / (2.0*eps);
  SM_ELEMENT_D(J, 6, 6) = -(qs_p - qs_m) / (2.0*eps);
  SM_ELEMENT_D(J, 8, 6) = -(qd_p - qd_m) / (2.0*eps);
  
  /* d/dy[7]: vgd -= eps (since vgd = y[5]-y[7]) */
  qg_p = fn_qgate(vgb, vgs, vgd-eps);
  qg_m = fn_qgate(vgb, vgs, vgd+eps);
  qs_p = fn_qsrc(vgb, vgs, vgd-eps);
  qs_m = fn_qsrc(vgb, vgs, vgd+eps);
  qd_p = fn_qdrain(vgb, vgs, vgd-eps);
  qd_m = fn_qdrain(vgb, vgs, vgd+eps);
  SM_ELEMENT_D(J, 4, 7) = -(qg_p - qg_m) / (2.0*eps);
  SM_ELEMENT_D(J, 6, 7) = -(qs_p - qs_m) / (2.0*eps);
  SM_ELEMENT_D(J, 8, 7) = -(qd_p - qd_m) / (2.0*eps);
  
  return 0;
}

static int mas_pump(sunrealtype t, SUNMatrix M, void* ud,
                    N_Vector t1, N_Vector t2, N_Vector t3) {
  (void)t; (void)ud; (void)t1; (void)t2; (void)t3;
  SUNMatZero(M);
  SM_ELEMENT_D(M, 0, 0) = 1.0;
  SM_ELEMENT_D(M, 1, 1) = 1.0;
  SM_ELEMENT_D(M, 1, 2) = 1.0;
  SM_ELEMENT_D(M, 2, 3) = 1.0;
  SM_ELEMENT_D(M, 2, 4) = 1.0;
  return 0;
}

int main() {
  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);
  int n = 9;
  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, 3, 3);

  N_Vector y0 = N_VNew_Serial(n, sunctx);
  sunrealtype* y0d = N_VGetArrayPointer(y0);
  y0d[0] = fn_qgate(0.0, 0.0, 0.0);
  for (int i = 1; i < 9; i++) y0d[i] = 0.0;

  Radau5Init(mem, rhs_pump, 0.0, y0);
  SUNMatrix Jt = SUNDenseMatrix(n, n, sunctx);
  Radau5SetLinearSolver(mem, Jt, NULL);
  Radau5SetJacFn(mem, jac_pump);  // Use analytic Jacobian!
  SUNMatrix Mt = SUNDenseMatrix(n, n, sunctx);
  Radau5SetMassFn(mem, mas_pump, Mt);
  Radau5SetDAEIndex(mem, 8, 1, 0);

  N_Vector rtol_v = N_VNew_Serial(n, sunctx);
  N_Vector atol_v = N_VNew_Serial(n, sunctx);
  sunrealtype* rtol_d = N_VGetArrayPointer(rtol_v);
  sunrealtype* atol_d = N_VGetArrayPointer(atol_v);
  sunrealtype tol = 1e-7;
  for (int i = 0; i < 5; i++) { rtol_d[i] = tol; atol_d[i] = tol*1e-6; }
  for (int i = 5; i < 8; i++) { rtol_d[i] = tol; atol_d[i] = tol; }
  rtol_d[8] = tol; atol_d[8] = tol;
  Radau5SVtolerances(mem, rtol_v, atol_v);
  
  Radau5SetInitStep(mem, 1e-12);
  Radau5SetMaxNumSteps(mem, 1000000);

  N_Vector yout = N_VNew_Serial(n, sunctx);
  sunrealtype tret;

  // Segmented integration
  int ndisc = 39;
  sunrealtype disc[41];
  disc[0] = 0.0;
  disc[1] = 50.0e-9;
  disc[2] = 60.0e-9;
  disc[3] = 110.0e-9;
  for (int i = 4; i <= 39; i++)
    disc[i] = disc[i % 4] + (sunrealtype)(i / 4) * 120.0e-9;
  disc[40] = 1200.0e-9;

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
  
  printf("\nFinal: ret=%d\n", ret);
  if (ret >= 0) {
    sunrealtype* yd = N_VGetArrayPointer(yout);
    for (int i = 0; i < 9; i++)
      printf("  y[%d] = %.14e\n", i, yd[i]);
    printf("  ref y[0] = 1.26280042987676e-13\n");
    printf("  ref y[8] = 1.52255686815578e-04\n");
  }

  long int nstep, naccpt, nfcn, njac, ndec;
  Radau5GetNumSteps(mem, &nstep);
  Radau5GetNumAccSteps(mem, &naccpt);
  Radau5GetNumRhsEvals(mem, &nfcn);
  Radau5GetNumJacEvals(mem, &njac);
  Radau5GetNumDecomps(mem, &ndec);
  printf("nstep=%ld naccpt=%ld nfcn=%ld njac=%ld ndec=%ld\n",
         nstep, naccpt, nfcn, njac, ndec);

  N_VDestroy(y0); N_VDestroy(yout);
  N_VDestroy(rtol_v); N_VDestroy(atol_v);
  SUNMatDestroy(Jt); SUNMatDestroy(Mt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return 0;
}
