/* ---------------------------------------------------------------------------
 * radau5_tba2.c -- Two Bit Adding Unit (DAE index-1, n=350)
 *
 * Faithful C translation of the Fortran source:
 *   IVPtestset_2.4/src/problems/tba.f
 *
 * M*y' = f(t,y) where M = diag(1,...,1,0,...,0)  (175 ones, 175 zeros)
 *
 * MOSFET circuit simulation of a two-bit adder computing:
 *   A1*2+A0 + B1*2+B0 + CIN = C*4 + S1*2 + S0
 *
 * State vector (0-based):
 *   y[0..174]   = charges g(U)
 *   y[175..349] = node voltages U
 *
 * t in [0, 315] with 63 discontinuities at t=5*k, k=1..63
 * Reference solution from IVPtestset (RADAU5, rtol=atol=1e-5, h0=4e-5)
 *
 * Output signals at t=315:
 *   y[223] (node 49)  ~ 0.2040  (S0)
 *   y[304] (node 130) ~ 4.9972  (S1)
 *   y[322] (node 148) ~ 0.2039  (C)
 * ---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_context.h>
#include <sundials/sundials_types.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

/* =========================================================================
 * Circuit parameters (Fortran COMMON /CONST/)
 * =========================================================================*/
typedef struct
{
  sunrealtype RGS, RGD, RBS, RBD;
  sunrealtype CGS, CGD, CBD, CBS;
  sunrealtype DELTA, CTIME, STIFF;
  sunrealtype CURIS, VTH, VDD, VBB;
  sunrealtype CLOAD, COUT;
} TbaParams;

/* =========================================================================
 * MOSFET model functions (Shichman-Hodges)
 * =========================================================================*/

/* Voltage-dependent bulk capacitance CBDBS */
static sunrealtype CBDBS(sunrealtype V, const TbaParams* p)
{
  sunrealtype PHIB = 0.87;
  if (V <= 0.0)
    return p->CBD / sqrt(1.0 - V / PHIB);
  else
    return p->CBD * (1.0 + V / (2.0 * PHIB));
}

/* Bulk-source junction current IBS */
static sunrealtype fn_IBS(sunrealtype VBS, const TbaParams* p)
{
  if (VBS <= 0.0)
    return -p->CURIS * (exp(VBS / p->VTH) - 1.0);
  else
    return 0.0;
}

/* Bulk-drain junction current IBD */
static sunrealtype fn_IBD(sunrealtype VBD, const TbaParams* p)
{
  if (VBD <= 0.0)
    return -p->CURIS * (exp(VBD / p->VTH) - 1.0);
  else
    return 0.0;
}

/* Drain current for VDS > 0: GDSP */
static sunrealtype GDSP(int NED, sunrealtype VDS, sunrealtype VGS,
                         sunrealtype VBS, int* ierr, const TbaParams* p)
{
  sunrealtype VT0, CGAMMA, PHI, BETA, VTE;

  if (NED == 0) {
    VT0 = -2.43; CGAMMA = 0.2; PHI = 1.28;
    BETA = 53.5e-6 * p->CTIME * p->STIFF;
  } else if (NED == 1) {
    VT0 = 0.2; CGAMMA = 0.035; PHI = 1.01;
    BETA = 4.0 * 43.7e-6 * p->CTIME * p->STIFF;
  } else if (NED == 2) {
    VT0 = 0.2; CGAMMA = 0.035; PHI = 1.01;
    BETA = 8.0 * 43.7e-6 * p->CTIME * p->STIFF;
  } else {
    VT0 = 0.2; CGAMMA = 0.035; PHI = 1.01;
    BETA = 12.0 * 43.7e-6 * p->CTIME * p->STIFF;
  }

  if (PHI - VBS < 0.0 || PHI < 0.0) { *ierr = -1; return 0.0; }

  VTE = VT0 + CGAMMA * (sqrt(PHI - VBS) - sqrt(PHI));

  if (VGS - VTE <= 0.0)
    return 0.0;
  else if (VGS - VTE <= VDS)
    return -BETA * (VGS - VTE) * (VGS - VTE) * (1.0 + p->DELTA * VDS);
  else
    return -BETA * VDS * (2.0 * (VGS - VTE) - VDS) * (1.0 + p->DELTA * VDS);
}

/* Drain current for VDS < 0: GDSM */
static sunrealtype GDSM(int NED, sunrealtype VDS, sunrealtype VGD,
                         sunrealtype VBD, int* ierr, const TbaParams* p)
{
  sunrealtype VT0, CGAMMA, PHI, BETA, VTE;

  if (NED == 0) {
    VT0 = -2.43; CGAMMA = 0.2; PHI = 1.28;
    BETA = 53.5e-6 * p->CTIME * p->STIFF;
  } else if (NED == 1) {
    VT0 = 0.2; CGAMMA = 0.035; PHI = 1.01;
    BETA = 4.0 * 43.7e-6 * p->CTIME * p->STIFF;
  } else if (NED == 2) {
    VT0 = 0.2; CGAMMA = 0.035; PHI = 1.01;
    BETA = 8.0 * 43.7e-6 * p->CTIME * p->STIFF;
  } else {
    VT0 = 0.2; CGAMMA = 0.035; PHI = 1.01;
    BETA = 12.0 * 43.7e-6 * p->CTIME * p->STIFF;
  }

  if (PHI - VBD < 0.0 || PHI < 0.0) { *ierr = -1; return 0.0; }

  VTE = VT0 + CGAMMA * (sqrt(PHI - VBD) - sqrt(PHI));

  if (VGD - VTE <= 0.0)
    return 0.0;
  else if (VGD - VTE <= -VDS)
    return BETA * (VGD - VTE) * (VGD - VTE) * (1.0 - p->DELTA * VDS);
  else
    return -BETA * VDS * (2.0 * (VGD - VTE) + VDS) * (1.0 - p->DELTA * VDS);
}
/* Drain-source current IDS (dispatches on sign of VDS) */
static sunrealtype fn_IDS(int NED, sunrealtype VDS, sunrealtype VGS,
                           sunrealtype VBS, sunrealtype VGD, sunrealtype VBD,
                           int* ierr, const TbaParams* p)
{
  sunrealtype val;
  if (VDS > 0.0) {
    val = GDSP(NED, VDS, VGS, VBS, ierr, p);
  } else if (VDS == 0.0) {
    val = 0.0;
  } else {
    val = GDSM(NED, VDS, VGD, VBD, ierr, p);
  }
  if (*ierr == -1) return val;
  return val;
}

/* =========================================================================
 * Input signal generation (PULSE subroutine)
 * =========================================================================*/
static void PULSE(sunrealtype X, sunrealtype* VIN, sunrealtype* VIND,
                  sunrealtype LOW, sunrealtype HIGH, sunrealtype DELAY,
                  sunrealtype T1, sunrealtype T2, sunrealtype T3,
                  sunrealtype PERIOD)
{
  sunrealtype TIME = fmod(X, PERIOD);

  if (TIME > (DELAY + T1 + T2 + T3)) {
    *VIN = LOW; *VIND = 0.0;
  } else if (TIME > (DELAY + T1 + T2)) {
    *VIN = ((HIGH - LOW) / T3) * (DELAY + T1 + T2 + T3 - TIME) + LOW;
    *VIND = -((HIGH - LOW) / T3);
  } else if (TIME > (DELAY + T1)) {
    *VIN = HIGH; *VIND = 0.0;
  } else if (TIME > DELAY) {
    *VIN = ((HIGH - LOW) / T1) * (TIME - DELAY) + LOW;
    *VIND = ((HIGH - LOW) / T1);
  } else {
    *VIN = LOW; *VIND = 0.0;
  }
}

/* =========================================================================
 * Gate subroutines for FCN (circuit equations)
 * All use 0-based indexing: Fortran I -> C I-1 passed as argument.
 * Inside gates, Y[I], Y[I+1], etc. with I already 0-based.
 * =========================================================================*/

/* NOR gate: 13 nodes I..I+12 (0-based) */
static void gate_NOR(int N, int I, sunrealtype U1, sunrealtype U2,
                     sunrealtype U1D, sunrealtype U2D,
                     sunrealtype* Y, sunrealtype* F, int* ierr,
                     const TbaParams* p)
{
  F[I] = -(Y[I] - Y[I+4]) / p->RGS
         - fn_IDS(0, Y[I+1]-Y[I], Y[I+4]-Y[I],
                  Y[I+2]-Y[I+4], Y[I+4]-Y[I+1], Y[I+3]-p->VDD, ierr, p);
  F[I+1] = -(Y[I+1] - p->VDD) / p->RGD
           + fn_IDS(0, Y[I+1]-Y[I], Y[I+4]-Y[I],
                    Y[I+2]-Y[I+4], Y[I+4]-Y[I+1], Y[I+3]-p->VDD, ierr, p);
  F[I+2] = -(Y[I+2] - p->VBB) / p->RBS + fn_IBS(Y[I+2]-Y[I+4], p);
  F[I+3] = -(Y[I+3] - p->VBB) / p->RBD + fn_IBD(Y[I+3]-p->VDD, p);

  /* Result node I+4 */
  F[I+4] = -(Y[I+4]-Y[I]) / p->RGS - fn_IBS(Y[I+2]-Y[I+4], p)
           -(Y[I+4]-Y[I+6]) / p->RGD - fn_IBD(Y[I+8]-Y[I+4], p)
           -(Y[I+4]-Y[I+10]) / p->RGD - fn_IBD(Y[I+12]-Y[I+4], p);

  F[I+5] = p->CGS*U1D - Y[I+5]/p->RGS
           - fn_IDS(1, Y[I+6]-Y[I+5], U1-Y[I+5],
                    Y[I+7], U1-Y[I+6], Y[I+8]-Y[I+4], ierr, p);
  F[I+6] = p->CGD*U1D - (Y[I+6]-Y[I+4])/p->RGD
           + fn_IDS(1, Y[I+6]-Y[I+5], U1-Y[I+5],
                    Y[I+7], U1-Y[I+6], Y[I+8]-Y[I+4], ierr, p);
  F[I+7] = -(Y[I+7] - p->VBB) / p->RBS + fn_IBS(Y[I+7], p);
  F[I+8] = -(Y[I+8] - p->VBB) / p->RBD + fn_IBD(Y[I+8]-Y[I+4], p);

  F[I+9] = p->CGS*U2D - Y[I+9]/p->RGS
           - fn_IDS(1, Y[I+10]-Y[I+9], U2-Y[I+9],
                    Y[I+11], U2-Y[I+10], Y[I+12]-Y[I+4], ierr, p);
  F[I+10] = p->CGD*U2D - (Y[I+10]-Y[I+4])/p->RGD
            + fn_IDS(1, Y[I+10]-Y[I+9], U2-Y[I+9],
                     Y[I+11], U2-Y[I+10], Y[I+12]-Y[I+4], ierr, p);
  F[I+11] = -(Y[I+11] - p->VBB) / p->RBS + fn_IBS(Y[I+11], p);
  F[I+12] = -(Y[I+12] - p->VBB) / p->RBD + fn_IBD(Y[I+12]-Y[I+4], p);

  if (*ierr == -1) return;
}
/* ANDOI gate: 18 nodes I..I+17 (0-based) */
static void gate_ANDOI(int N, int I, sunrealtype U1, sunrealtype U2,
                       sunrealtype U3, sunrealtype U1D, sunrealtype U2D,
                       sunrealtype U3D, sunrealtype* Y, sunrealtype* F,
                       int* ierr, const TbaParams* p)
{
  F[I] = -(Y[I]-Y[I+4])/p->RGS
         - fn_IDS(0, Y[I+1]-Y[I], Y[I+4]-Y[I],
                  Y[I+2]-Y[I+4], Y[I+4]-Y[I+1], Y[I+3]-p->VDD, ierr, p);
  F[I+1] = -(Y[I+1]-p->VDD)/p->RGD
           + fn_IDS(0, Y[I+1]-Y[I], Y[I+4]-Y[I],
                    Y[I+2]-Y[I+4], Y[I+4]-Y[I+1], Y[I+3]-p->VDD, ierr, p);
  F[I+2] = -(Y[I+2]-p->VBB)/p->RBS + fn_IBS(Y[I+2]-Y[I+4], p);
  F[I+3] = -(Y[I+3]-p->VBB)/p->RBD + fn_IBD(Y[I+3]-p->VDD, p);

  /* Result node I+4 */
  F[I+4] = -(Y[I+4]-Y[I])/p->RGS - fn_IBS(Y[I+2]-Y[I+4], p)
           -(Y[I+4]-Y[I+6])/p->RGD - fn_IBD(Y[I+8]-Y[I+4], p)
           -(Y[I+4]-Y[I+10])/p->RGD - fn_IBD(Y[I+12]-Y[I+4], p);

  F[I+5] = p->CGS*U1D - Y[I+5]/p->RGS
           - fn_IDS(1, Y[I+6]-Y[I+5], U1-Y[I+5],
                    Y[I+7], U1-Y[I+6], Y[I+8]-Y[I+4], ierr, p);
  F[I+6] = p->CGD*U1D - (Y[I+6]-Y[I+4])/p->RGD
           + fn_IDS(1, Y[I+6]-Y[I+5], U1-Y[I+5],
                    Y[I+7], U1-Y[I+6], Y[I+8]-Y[I+4], ierr, p);
  F[I+7] = -(Y[I+7]-p->VBB)/p->RBS + fn_IBS(Y[I+7], p);
  F[I+8] = -(Y[I+8]-p->VBB)/p->RBD + fn_IBD(Y[I+8]-Y[I+4], p);

  F[I+9] = p->CGS*U2D - (Y[I+9]-Y[I+13])/p->RGS
           - fn_IDS(2, Y[I+10]-Y[I+9], U2-Y[I+9],
                    Y[I+11]-Y[I+13], U2-Y[I+10], Y[I+12]-Y[I+4], ierr, p);
  F[I+10] = p->CGD*U2D - (Y[I+10]-Y[I+4])/p->RGD
            + fn_IDS(2, Y[I+10]-Y[I+9], U2-Y[I+9],
                     Y[I+11]-Y[I+13], U2-Y[I+10], Y[I+12]-Y[I+4], ierr, p);
  F[I+11] = -(Y[I+11]-p->VBB)/p->RBS + fn_IBS(Y[I+11]-Y[I+13], p);
  F[I+12] = -(Y[I+12]-p->VBB)/p->RBD + fn_IBD(Y[I+12]-Y[I+4], p);

  /* Coupling node I+13 */
  F[I+13] = -(Y[I+13]-Y[I+9])/p->RGS - fn_IBS(Y[I+11]-Y[I+13], p)
            -(Y[I+13]-Y[I+15])/p->RGD - fn_IBD(Y[I+17]-Y[I+13], p);

  F[I+14] = p->CGS*U3D - Y[I+14]/p->RGS
            - fn_IDS(2, Y[I+15]-Y[I+14], U3-Y[I+14],
                     Y[I+16], U3-Y[I+15], Y[I+17]-Y[I+13], ierr, p);
  F[I+15] = p->CGD*U3D - (Y[I+15]-Y[I+13])/p->RGD
            + fn_IDS(2, Y[I+15]-Y[I+14], U3-Y[I+14],
                     Y[I+16], U3-Y[I+15], Y[I+17]-Y[I+13], ierr, p);
  F[I+16] = -(Y[I+16]-p->VBB)/p->RBS + fn_IBS(Y[I+16], p);
  F[I+17] = -(Y[I+17]-p->VBB)/p->RBD + fn_IBD(Y[I+17]-Y[I+13], p);

  if (*ierr == -1) return;
}
/* NAND gate: 14 nodes I..I+13 (0-based) */
static void gate_NAND(int N, int I, sunrealtype U1, sunrealtype U2,
                      sunrealtype U1D, sunrealtype U2D,
                      sunrealtype* Y, sunrealtype* F, int* ierr,
                      const TbaParams* p)
{
  F[I] = -(Y[I]-Y[I+4])/p->RGS
         - fn_IDS(0, Y[I+1]-Y[I], Y[I+4]-Y[I],
                  Y[I+2]-Y[I+4], Y[I+4]-Y[I+1], Y[I+3]-p->VDD, ierr, p);
  F[I+1] = -(Y[I+1]-p->VDD)/p->RGD
           + fn_IDS(0, Y[I+1]-Y[I], Y[I+4]-Y[I],
                    Y[I+2]-Y[I+4], Y[I+4]-Y[I+1], Y[I+3]-p->VDD, ierr, p);
  F[I+2] = -(Y[I+2]-p->VBB)/p->RBS + fn_IBS(Y[I+2]-Y[I+4], p);
  F[I+3] = -(Y[I+3]-p->VBB)/p->RBD + fn_IBD(Y[I+3]-p->VDD, p);

  /* Result node I+4 */
  F[I+4] = -(Y[I+4]-Y[I])/p->RGS - fn_IBS(Y[I+2]-Y[I+4], p)
           -(Y[I+4]-Y[I+6])/p->RGD - fn_IBD(Y[I+8]-Y[I+4], p);

  F[I+5] = p->CGS*U1D - (Y[I+5]-Y[I+9])/p->RGS
           - fn_IDS(2, Y[I+6]-Y[I+5], U1-Y[I+5],
                    Y[I+7]-Y[I+9], U1-Y[I+6], Y[I+8]-Y[I+4], ierr, p);
  F[I+6] = p->CGD*U1D - (Y[I+6]-Y[I+4])/p->RGD
           + fn_IDS(2, Y[I+6]-Y[I+5], U1-Y[I+5],
                    Y[I+7]-Y[I+9], U1-Y[I+6], Y[I+8]-Y[I+4], ierr, p);
  F[I+7] = -(Y[I+7]-p->VBB)/p->RBS + fn_IBS(Y[I+7]-Y[I+9], p);
  F[I+8] = -(Y[I+8]-p->VBB)/p->RBD + fn_IBD(Y[I+8]-Y[I+4], p);

  /* Coupling node I+9 */
  F[I+9] = -(Y[I+9]-Y[I+5])/p->RGS - fn_IBS(Y[I+7]-Y[I+9], p)
           -(Y[I+9]-Y[I+11])/p->RGD - fn_IBD(Y[I+13]-Y[I+9], p);

  F[I+10] = p->CGS*U2D - Y[I+10]/p->RGS
            - fn_IDS(2, Y[I+11]-Y[I+10], U2-Y[I+10],
                     Y[I+12], U2-Y[I+11], Y[I+13]-Y[I+9], ierr, p);
  F[I+11] = p->CGD*U2D - (Y[I+11]-Y[I+9])/p->RGD
            + fn_IDS(2, Y[I+11]-Y[I+10], U2-Y[I+10],
                     Y[I+12], U2-Y[I+11], Y[I+13]-Y[I+9], ierr, p);
  F[I+12] = -(Y[I+12]-p->VBB)/p->RBS + fn_IBS(Y[I+12], p);
  F[I+13] = -(Y[I+13]-p->VBB)/p->RBD + fn_IBD(Y[I+13]-Y[I+9], p);

  if (*ierr == -1) return;
}
/* ORANI gate: 18 nodes I..I+17 (0-based) */
static void gate_ORANI(int N, int I, sunrealtype U1, sunrealtype U2,
                       sunrealtype U3, sunrealtype U1D, sunrealtype U2D,
                       sunrealtype U3D, sunrealtype* Y, sunrealtype* F,
                       int* ierr, const TbaParams* p)
{
  F[I] = -(Y[I]-Y[I+4])/p->RGS
         - fn_IDS(0, Y[I+1]-Y[I], Y[I+4]-Y[I],
                  Y[I+2]-Y[I+4], Y[I+4]-Y[I+1], Y[I+3]-p->VDD, ierr, p);
  F[I+1] = -(Y[I+1]-p->VDD)/p->RGD
           + fn_IDS(0, Y[I+1]-Y[I], Y[I+4]-Y[I],
                    Y[I+2]-Y[I+4], Y[I+4]-Y[I+1], Y[I+3]-p->VDD, ierr, p);
  F[I+2] = -(Y[I+2]-p->VBB)/p->RBS + fn_IBS(Y[I+2]-Y[I+4], p);
  F[I+3] = -(Y[I+3]-p->VBB)/p->RBD + fn_IBD(Y[I+3]-p->VDD, p);

  /* Result node I+4 */
  F[I+4] = -(Y[I+4]-Y[I])/p->RGS - fn_IBS(Y[I+2]-Y[I+4], p)
           -(Y[I+4]-Y[I+6])/p->RGD - fn_IBD(Y[I+8]-Y[I+4], p);

  F[I+5] = p->CGS*U1D - (Y[I+5]-Y[I+9])/p->RGS
           - fn_IDS(2, Y[I+6]-Y[I+5], U1-Y[I+5],
                    Y[I+7]-Y[I+9], U1-Y[I+6], Y[I+8]-Y[I+4], ierr, p);
  F[I+6] = p->CGD*U1D - (Y[I+6]-Y[I+4])/p->RGD
           + fn_IDS(2, Y[I+6]-Y[I+5], U1-Y[I+5],
                    Y[I+7]-Y[I+9], U1-Y[I+6], Y[I+8]-Y[I+4], ierr, p);
  F[I+7] = -(Y[I+7]-p->VBB)/p->RBS + fn_IBS(Y[I+7]-Y[I+9], p);
  F[I+8] = -(Y[I+8]-p->VBB)/p->RBD + fn_IBD(Y[I+8]-Y[I+4], p);

  /* Coupling node I+9 */
  F[I+9] = -(Y[I+9]-Y[I+5])/p->RGS - fn_IBS(Y[I+7]-Y[I+9], p)
           -(Y[I+9]-Y[I+11])/p->RGD - fn_IBD(Y[I+13]-Y[I+9], p)
           -(Y[I+9]-Y[I+15])/p->RGD - fn_IBD(Y[I+17]-Y[I+9], p);

  F[I+10] = p->CGS*U2D - Y[I+10]/p->RGS
            - fn_IDS(2, Y[I+11]-Y[I+10], U2-Y[I+10],
                     Y[I+12], U2-Y[I+11], Y[I+13]-Y[I+9], ierr, p);
  F[I+11] = p->CGD*U2D - (Y[I+11]-Y[I+9])/p->RGD
            + fn_IDS(2, Y[I+11]-Y[I+10], U2-Y[I+10],
                     Y[I+12], U2-Y[I+11], Y[I+13]-Y[I+9], ierr, p);
  F[I+12] = -(Y[I+12]-p->VBB)/p->RBS + fn_IBS(Y[I+12], p);
  F[I+13] = -(Y[I+13]-p->VBB)/p->RBD + fn_IBD(Y[I+13]-Y[I+9], p);

  F[I+14] = p->CGS*U3D - Y[I+14]/p->RGS
            - fn_IDS(2, Y[I+15]-Y[I+14], U3-Y[I+14],
                     Y[I+16], U3-Y[I+15], Y[I+17]-Y[I+9], ierr, p);
  F[I+15] = p->CGD*U3D - (Y[I+15]-Y[I+9])/p->RGD
            + fn_IDS(2, Y[I+15]-Y[I+14], U3-Y[I+14],
                     Y[I+16], U3-Y[I+15], Y[I+17]-Y[I+9], ierr, p);
  F[I+16] = -(Y[I+16]-p->VBB)/p->RBS + fn_IBS(Y[I+16], p);
  F[I+17] = -(Y[I+17]-p->VBB)/p->RBD + fn_IBD(Y[I+17]-Y[I+9], p);

  if (*ierr == -1) return;
}
/* ANDOIP gate (ANDOI with capacitive coupling at result node):
 * 18 nodes I..I+17 (0-based).  Result node I+4 has extra coupling
 * to Fortran nodes 163,165 -> C indices 162,164. */
static void gate_ANDOIP(int N, int I, sunrealtype U1, sunrealtype U2,
                        sunrealtype U3, sunrealtype U1D, sunrealtype U2D,
                        sunrealtype U3D, sunrealtype* Y, sunrealtype* F,
                        int* ierr, const TbaParams* p)
{
  F[I] = -(Y[I]-Y[I+4])/p->RGS
         - fn_IDS(0, Y[I+1]-Y[I], Y[I+4]-Y[I],
                  Y[I+2]-Y[I+4], Y[I+4]-Y[I+1], Y[I+3]-p->VDD, ierr, p);
  F[I+1] = -(Y[I+1]-p->VDD)/p->RGD
           + fn_IDS(0, Y[I+1]-Y[I], Y[I+4]-Y[I],
                    Y[I+2]-Y[I+4], Y[I+4]-Y[I+1], Y[I+3]-p->VDD, ierr, p);
  F[I+2] = -(Y[I+2]-p->VBB)/p->RBS + fn_IBS(Y[I+2]-Y[I+4], p);
  F[I+3] = -(Y[I+3]-p->VBB)/p->RBD + fn_IBD(Y[I+3]-p->VDD, p);

  /* Result node I+4 -- extra coupling to Fortran Y(163),Y(165) = C Y[162],Y[164] */
  F[I+4] = -(Y[I+4]-Y[I])/p->RGS - fn_IBS(Y[I+2]-Y[I+4], p)
           -(Y[I+4]-Y[I+6])/p->RGD - fn_IBD(Y[I+8]-Y[I+4], p)
           -(Y[I+4]-Y[I+10])/p->RGD - fn_IBD(Y[I+12]-Y[I+4], p)
           -(Y[I+4]-Y[162])/p->RGD - fn_IBD(Y[164]-Y[I+4], p);

  F[I+5] = p->CGS*U1D - Y[I+5]/p->RGS
           - fn_IDS(1, Y[I+6]-Y[I+5], U1-Y[I+5],
                    Y[I+7], U1-Y[I+6], Y[I+8]-Y[I+4], ierr, p);
  F[I+6] = p->CGD*U1D - (Y[I+6]-Y[I+4])/p->RGD
           + fn_IDS(1, Y[I+6]-Y[I+5], U1-Y[I+5],
                    Y[I+7], U1-Y[I+6], Y[I+8]-Y[I+4], ierr, p);
  F[I+7] = -(Y[I+7]-p->VBB)/p->RBS + fn_IBS(Y[I+7], p);
  F[I+8] = -(Y[I+8]-p->VBB)/p->RBD + fn_IBD(Y[I+8]-Y[I+4], p);

  F[I+9] = p->CGS*U2D - (Y[I+9]-Y[I+13])/p->RGS
           - fn_IDS(2, Y[I+10]-Y[I+9], U2-Y[I+9],
                    Y[I+11]-Y[I+13], U2-Y[I+10], Y[I+12]-Y[I+4], ierr, p);
  F[I+10] = p->CGD*U2D - (Y[I+10]-Y[I+4])/p->RGD
            + fn_IDS(2, Y[I+10]-Y[I+9], U2-Y[I+9],
                     Y[I+11]-Y[I+13], U2-Y[I+10], Y[I+12]-Y[I+4], ierr, p);
  F[I+11] = -(Y[I+11]-p->VBB)/p->RBS + fn_IBS(Y[I+11]-Y[I+13], p);
  F[I+12] = -(Y[I+12]-p->VBB)/p->RBD + fn_IBD(Y[I+12]-Y[I+4], p);

  /* Coupling node I+13 */
  F[I+13] = -(Y[I+13]-Y[I+9])/p->RGS - fn_IBS(Y[I+11]-Y[I+13], p)
            -(Y[I+13]-Y[I+15])/p->RGD - fn_IBD(Y[I+17]-Y[I+13], p);

  F[I+14] = p->CGS*U3D - Y[I+14]/p->RGS
            - fn_IDS(2, Y[I+15]-Y[I+14], U3-Y[I+14],
                     Y[I+16], U3-Y[I+15], Y[I+17]-Y[I+13], ierr, p);
  F[I+15] = p->CGD*U3D - (Y[I+15]-Y[I+13])/p->RGD
            + fn_IDS(2, Y[I+15]-Y[I+14], U3-Y[I+14],
                     Y[I+16], U3-Y[I+15], Y[I+17]-Y[I+13], ierr, p);
  F[I+16] = -(Y[I+16]-p->VBB)/p->RBS + fn_IBS(Y[I+16], p);
  F[I+17] = -(Y[I+17]-p->VBB)/p->RBD + fn_IBD(Y[I+17]-Y[I+13], p);

  if (*ierr == -1) return;
}
/* =========================================================================
 * FCN: Circuit equations for 175 voltage nodes (0-based Y[0..174])
 * Translated from Fortran FCN subroutine.
 * =========================================================================*/
static void FCN(int N, sunrealtype X, sunrealtype* Y, sunrealtype* F,
                int* ierr, const TbaParams* p)
{
  sunrealtype V1, V1D, V2, V2D, V3, V3D, V4, V4D, CIN, CIND;

  /* Input signals */
  PULSE(X, &V1, &V1D, 0.0, 5.0, 0.0, 5.0, 5.0, 5.0, 20.0);
  PULSE(X, &V2, &V2D, 0.0, 5.0, 10.0, 5.0, 15.0, 5.0, 40.0);
  PULSE(X, &V3, &V3D, 0.0, 5.0, 30.0, 5.0, 35.0, 5.0, 80.0);
  PULSE(X, &V4, &V4D, 0.0, 5.0, 70.0, 5.0, 75.0, 5.0, 160.0);
  PULSE(X, &CIN, &CIND, 0.0, 5.0, 150.0, 5.0, 155.0, 5.0, 320.0);

  /* NOR-gate 1: Fortran nodes 1--13 -> C 0--12 */
  gate_NOR(N, 0, V1, V2, V1D, V2D, Y, F, ierr, p);

  /* ANDOI-gate 1: Fortran nodes 14--31 -> C 13--30
   * Fortran: CALL ANDOI(N,14,Y(5),V2,V1,0,V2D,V1D,...) */
  gate_ANDOI(N, 13, Y[4], V2, V1, 0.0, V2D, V1D, Y, F, ierr, p);

  /* NOR-gate 2: Fortran nodes 32--44 -> C 31--43
   * Fortran: CALL NOR(N,32,Y(18),CIN,0,CIND,...) */
  gate_NOR(N, 31, Y[17], CIN, 0.0, CIND, Y, F, ierr, p);

  /* ANDOI-gate 2: Fortran nodes 45--62 -> C 44--61
   * Fortran: CALL ANDOI(N,45,Y(36),CIN,Y(18),0,CIND,0,...) */
  gate_ANDOI(N, 44, Y[35], CIN, Y[17], 0.0, CIND, 0.0, Y, F, ierr, p);

  /* ANDOI-gate 3: Fortran nodes 63--80 -> C 62--79
   * Fortran: CALL ANDOI(N,63,Y(5),CIN,Y(18),0,CIND,0,...) */
  gate_ANDOI(N, 62, Y[4], CIN, Y[17], 0.0, CIND, 0.0, Y, F, ierr, p);

  /* NOR-gate 3: Fortran nodes 81--93 -> C 80--92 */
  gate_NOR(N, 80, V3, V4, V3D, V4D, Y, F, ierr, p);

  /* ANDOI-gate 4: Fortran nodes 94--111 -> C 93--110
   * Fortran: CALL ANDOI(N,94,Y(85),V4,V3,0,V4D,V3D,...) */
  gate_ANDOI(N, 93, Y[84], V4, V3, 0.0, V4D, V3D, Y, F, ierr, p);

  /* NAND-gate: Fortran nodes 112--125 -> C 111--124
   * Fortran: CALL NAND(N,112,Y(67),Y(98),0,0,...) */
  gate_NAND(N, 111, Y[66], Y[97], 0.0, 0.0, Y, F, ierr, p);

  /* ORANI-gate 1: Fortran nodes 126--143 -> C 125--142
   * Fortran: CALL ORANI(N,126,Y(116),Y(67),Y(98),0,0,0,...) */
  gate_ORANI(N, 125, Y[115], Y[66], Y[97], 0.0, 0.0, 0.0, Y, F, ierr, p);

  /* ANDOIP-gate 5: Fortran nodes 144--161 -> C 143--160
   * Fortran: CALL ANDOIP(N,144,Y(85),Y(5),Y(98),0,0,0,...) */
  gate_ANDOIP(N, 143, Y[84], Y[4], Y[97], 0.0, 0.0, 0.0, Y, F, ierr, p);
  /* Three additional enhancement transistors in series
   * Fortran nodes 162--175 -> C indices 161--174
   *
   * First transistor:  DRAIN=node148(C:147), GATE=node98(C:97), SOURCE=node166(C:165)
   * Second transistor: DRAIN=node166(C:165), GATE=node18(C:17), SOURCE=node171(C:170)
   * Third transistor:  DRAIN=node171(C:170), GATE=CIN, SOURCE=MASS(0)
   */

  /* Fortran F(162) -> C F[161] */
  F[161] = -(Y[161]-Y[165])/p->RGS
           - fn_IDS(3, Y[162]-Y[161], Y[97]-Y[161],
                    Y[163]-Y[165], Y[97]-Y[162], Y[164]-Y[147], ierr, p);
  /* Fortran F(163) -> C F[162] */
  F[162] = -(Y[162]-Y[147])/p->RGD
           + fn_IDS(3, Y[162]-Y[161], Y[97]-Y[161],
                    Y[163]-Y[165], Y[97]-Y[162], Y[164]-Y[147], ierr, p);
  /* Fortran F(164) -> C F[163] */
  F[163] = -(Y[163]-p->VBB)/p->RBS + fn_IBS(Y[163]-Y[165], p);
  /* Fortran F(165) -> C F[164] */
  F[164] = -(Y[164]-p->VBB)/p->RBD + fn_IBD(Y[164]-Y[147], p);

  /* Fortran F(166) -> C F[165] */
  F[165] = -fn_IBS(Y[163]-Y[165], p) - (Y[165]-Y[161])/p->RGS
           -fn_IBD(Y[169]-Y[165], p) - (Y[165]-Y[167])/p->RGD;

  /* Fortran F(167) -> C F[166] */
  F[166] = -(Y[166]-Y[170])/p->RGS
           - fn_IDS(3, Y[167]-Y[166], Y[17]-Y[166],
                    Y[168]-Y[170], Y[17]-Y[167], Y[169]-Y[165], ierr, p);
  /* Fortran F(168) -> C F[167] */
  F[167] = -(Y[167]-Y[165])/p->RGD
           + fn_IDS(3, Y[167]-Y[166], Y[17]-Y[166],
                    Y[168]-Y[170], Y[17]-Y[167], Y[169]-Y[165], ierr, p);
  /* Fortran F(169) -> C F[168] */
  F[168] = -(Y[168]-p->VBB)/p->RBS + fn_IBS(Y[168]-Y[170], p);
  /* Fortran F(170) -> C F[169] */
  F[169] = -(Y[169]-p->VBB)/p->RBD + fn_IBD(Y[169]-Y[165], p);

  /* Fortran F(171) -> C F[170] */
  F[170] = -fn_IBS(Y[168]-Y[170], p) - (Y[170]-Y[166])/p->RGS
           -fn_IBD(Y[174]-Y[170], p) - (Y[170]-Y[172])/p->RGD;

  /* Fortran F(172) -> C F[171] */
  F[171] = p->CGS*CIND - Y[171]/p->RGS
           - fn_IDS(3, Y[172]-Y[171], CIN-Y[171],
                    Y[173], CIN-Y[172], Y[174]-Y[170], ierr, p);
  /* Fortran F(173) -> C F[172] */
  F[172] = p->CGD*CIND - (Y[172]-Y[170])/p->RGD
           + fn_IDS(3, Y[172]-Y[171], CIN-Y[171],
                    Y[173], CIN-Y[172], Y[174]-Y[170], ierr, p);
  /* Fortran F(174) -> C F[173] */
  F[173] = -(Y[173]-p->VBB)/p->RBS + fn_IBS(Y[173], p);
  /* Fortran F(175) -> C F[174] */
  F[174] = -(Y[174]-p->VBB)/p->RBD + fn_IBD(Y[174]-Y[170], p);

  if (*ierr == -1) return;
}
/* =========================================================================
 * Charge gate subroutines for GCN (0-based indexing)
 * =========================================================================*/

/* Charge function for NOR gate: 13 nodes I..I+12 */
static void charge_NOR(int N, const sunrealtype* U, int I, sunrealtype* G,
                       const TbaParams* p)
{
  G[I]   += p->CGS * (U[I] - U[I+4]);
  G[I+1] += p->CGD * (U[I+1] - U[I+4]);
  G[I+2] += CBDBS(U[I+2]-U[I+4], p) * (U[I+2]-U[I+4]);
  G[I+3] += CBDBS(U[I+3]-p->VDD, p) * U[I+3];
  G[I+4] += p->CGS * (U[I+4]-U[I]) + p->CGD * (U[I+4]-U[I+1])
          + CBDBS(U[I+2]-U[I+4], p) * (U[I+4]-U[I+2])
          + CBDBS(U[I+8]-U[I+4], p) * (U[I+4]-U[I+8])
          + CBDBS(U[I+12]-U[I+4], p) * (U[I+4]-U[I+12])
          + p->CLOAD * U[I+4];
  G[I+5]  += p->CGS * U[I+5];
  G[I+6]  += p->CGD * U[I+6];
  G[I+7]  += CBDBS(U[I+7], p) * U[I+7];
  G[I+8]  += CBDBS(U[I+8]-U[I+4], p) * (U[I+8]-U[I+4]);
  G[I+9]  += p->CGS * U[I+9];
  G[I+10] += p->CGD * U[I+10];
  G[I+11] += CBDBS(U[I+11], p) * U[I+11];
  /* Fortran: G(I+12)=CBDBS(U(I+12)-U(I+4))*(U(I+12)-U(I+14)) */
  G[I+12] += CBDBS(U[I+12]-U[I+4], p) * (U[I+12]-U[I+14]);
}
/* Charge function for ANDOI gate: 18 nodes I..I+17 */
static void charge_ANDOI(int N, const sunrealtype* U, int I, sunrealtype* G,
                         const TbaParams* p)
{
  G[I]   += p->CGS * (U[I] - U[I+4]);
  G[I+1] += p->CGD * (U[I+1] - U[I+4]);
  G[I+2] += CBDBS(U[I+2]-U[I+4], p) * (U[I+2]-U[I+4]);
  G[I+3] += CBDBS(U[I+3]-p->VDD, p) * U[I+3];
  G[I+4] += p->CGS * (U[I+4]-U[I]) + p->CGD * (U[I+4]-U[I+1])
          + CBDBS(U[I+2]-U[I+4], p) * (U[I+4]-U[I+2])
          + CBDBS(U[I+8]-U[I+4], p) * (U[I+4]-U[I+8])
          + CBDBS(U[I+12]-U[I+4], p) * (U[I+4]-U[I+12])
          + p->CLOAD * U[I+4];
  G[I+5]  += p->CGS * U[I+5];
  G[I+6]  += p->CGD * U[I+6];
  G[I+7]  += CBDBS(U[I+7], p) * U[I+7];
  G[I+8]  += CBDBS(U[I+8]-U[I+4], p) * (U[I+8]-U[I+4]);
  G[I+9]  += p->CGS * U[I+9];
  G[I+10] += p->CGD * U[I+10];
  G[I+11] += CBDBS(U[I+11]-U[I+13], p) * (U[I+11]-U[I+13]);
  G[I+12] += CBDBS(U[I+12]-U[I+4], p) * (U[I+12]-U[I+4]);
  G[I+13] += CBDBS(U[I+11]-U[I+13], p) * (U[I+13]-U[I+11])
           + CBDBS(U[I+17]-U[I+13], p) * (U[I+13]-U[I+17])
           + p->CLOAD * U[I+13];
  G[I+14] += p->CGS * U[I+14];
  G[I+15] += p->CGD * U[I+15];
  G[I+16] += CBDBS(U[I+16], p) * U[I+16];
  G[I+17] += CBDBS(U[I+17]-U[I+13], p) * (U[I+17]-U[I+13]);
}
/* Charge function for NAND gate: 14 nodes I..I+13 */
static void charge_NAND(int N, const sunrealtype* U, int I, sunrealtype* G,
                        const TbaParams* p)
{
  G[I]   += p->CGS * (U[I] - U[I+4]);
  G[I+1] += p->CGD * (U[I+1] - U[I+4]);
  G[I+2] += CBDBS(U[I+2]-U[I+4], p) * (U[I+2]-U[I+4]);
  G[I+3] += CBDBS(U[I+3]-p->VDD, p) * U[I+3];
  G[I+4] += p->CGS * (U[I+4]-U[I]) + p->CGD * (U[I+4]-U[I+1])
          + CBDBS(U[I+2]-U[I+4], p) * (U[I+4]-U[I+2])
          + CBDBS(U[I+8]-U[I+4], p) * (U[I+4]-U[I+8])
          + p->CLOAD * U[I+4];
  G[I+5]  += p->CGS * U[I+5];
  G[I+6]  += p->CGD * U[I+6];
  G[I+7]  += CBDBS(U[I+7]-U[I+9], p) * (U[I+7]-U[I+9]);
  G[I+8]  += CBDBS(U[I+8]-U[I+4], p) * (U[I+8]-U[I+4]);
  G[I+9]  += CBDBS(U[I+7]-U[I+9], p) * (U[I+9]-U[I+7])
           + CBDBS(U[I+13]-U[I+9], p) * (U[I+9]-U[I+13])
           + p->CLOAD * U[I+9];
  G[I+10] += p->CGS * U[I+10];
  G[I+11] += p->CGD * U[I+11];
  G[I+12] += CBDBS(U[I+12], p) * U[I+12];
  G[I+13] += CBDBS(U[I+13]-U[I+9], p) * (U[I+13]-U[I+9]);
}
/* Charge function for ORANI gate: 18 nodes I..I+17 */
static void charge_ORANI(int N, const sunrealtype* U, int I, sunrealtype* G,
                         const TbaParams* p)
{
  G[I]   += p->CGS * (U[I] - U[I+4]);
  G[I+1] += p->CGD * (U[I+1] - U[I+4]);
  G[I+2] += CBDBS(U[I+2]-U[I+4], p) * (U[I+2]-U[I+4]);
  G[I+3] += CBDBS(U[I+3]-p->VDD, p) * U[I+3];
  G[I+4] += p->CGS * (U[I+4]-U[I]) + p->CGD * (U[I+4]-U[I+1])
          + CBDBS(U[I+2]-U[I+4], p) * (U[I+4]-U[I+2])
          + CBDBS(U[I+8]-U[I+4], p) * (U[I+4]-U[I+8])
          + p->CLOAD * U[I+4];
  G[I+5]  += p->CGS * U[I+5];
  G[I+6]  += p->CGD * U[I+6];
  G[I+7]  += CBDBS(U[I+7]-U[I+9], p) * (U[I+7]-U[I+9]);
  G[I+8]  += CBDBS(U[I+8]-U[I+4], p) * (U[I+8]-U[I+4]);
  G[I+9]  += CBDBS(U[I+7]-U[I+9], p) * (U[I+9]-U[I+7])
           + CBDBS(U[I+13]-U[I+9], p) * (U[I+9]-U[I+13])
           + CBDBS(U[I+17]-U[I+9], p) * (U[I+9]-U[I+17])
           + p->CLOAD * U[I+9];
  G[I+10] += p->CGS * U[I+10];
  G[I+11] += p->CGD * U[I+11];
  G[I+12] += CBDBS(U[I+12], p) * U[I+12];
  G[I+13] += CBDBS(U[I+13]-U[I+9], p) * (U[I+13]-U[I+9]);
  G[I+14] += p->CGS * U[I+14];
  G[I+15] += p->CGD * U[I+15];
  G[I+16] += CBDBS(U[I+16], p) * U[I+16];
  G[I+17] += CBDBS(U[I+17]-U[I+9], p) * (U[I+17]-U[I+9]);
}
/* =========================================================================
 * GCN: Charge function for 175 voltage nodes (0-based U[0..174], G[0..174])
 * Translated from Fortran GCN subroutine.
 * =========================================================================*/
static void GCN(int N, const sunrealtype* U, sunrealtype* G,
                const TbaParams* p)
{
  int i;
  for (i = 0; i < N; i++) G[i] = 0.0;

  /* Ten logical subcircuits (Fortran 1-based -> C 0-based) */
  charge_NOR(N, U, 0, G, p);       /* NOR-gate 1: Fortran 1 -> C 0 */
  charge_ANDOI(N, U, 13, G, p);    /* ANDOI-gate 1: Fortran 14 -> C 13 */
  charge_NOR(N, U, 31, G, p);      /* NOR-gate 2: Fortran 32 -> C 31 */
  charge_ANDOI(N, U, 44, G, p);    /* ANDOI-gate 2: Fortran 45 -> C 44 */
  charge_ANDOI(N, U, 62, G, p);    /* ANDOI-gate 3: Fortran 63 -> C 62 */
  charge_NOR(N, U, 80, G, p);      /* NOR-gate 3: Fortran 81 -> C 80 */
  charge_ANDOI(N, U, 93, G, p);    /* ANDOI-gate 4: Fortran 94 -> C 93 */
  charge_NAND(N, U, 111, G, p);    /* NAND-gate: Fortran 112 -> C 111 */
  charge_ORANI(N, U, 125, G, p);   /* ORANI-gate 1: Fortran 126 -> C 125 */
  charge_ANDOI(N, U, 143, G, p);   /* ANDOI-gate 5: Fortran 144 -> C 143 */

  /* Capacitive coupling: result node NOR-gate 1
   * Fortran G(5) -> C G[4], Fortran U(19)->C U[18], etc. */
  G[4] += p->CGS*(U[4]-U[18]) + p->CGD*(U[4]-U[19])
        + p->CGS*(U[4]-U[67]) + p->CGD*(U[4]-U[68])
        + p->CGS*(U[4]-U[152]) + p->CGD*(U[4]-U[153]);
  G[18]  -= p->CGS * U[4];
  G[19]  -= p->CGD * U[4];
  G[67]  -= p->CGS * U[4];
  G[68]  -= p->CGD * U[4];
  G[152] -= p->CGS * U[4];
  G[153] -= p->CGD * U[4];

  /* Capacitive coupling: result node ANDOI-gate 1
   * Fortran G(18) -> C G[17] */
  G[17] += p->CGS*(U[17]-U[36]) + p->CGD*(U[17]-U[37])
         + p->CGS*(U[17]-U[58]) + p->CGD*(U[17]-U[59])
         + p->CGS*(U[17]-U[76]) + p->CGD*(U[17]-U[77])
         + p->CGS*(U[17]-U[166]) + p->CGD*(U[17]-U[167]);
  G[36] -= p->CGS * U[17];
  G[37] -= p->CGD * U[17];
  G[58] -= p->CGS * U[17];
  G[59] -= p->CGD * U[17];
  G[76] -= p->CGS * U[17];
  G[77] -= p->CGD * U[17];

  /* Capacitive coupling: result node NOR-gate 2
   * Fortran G(36) -> C G[35] */
  G[35] += p->CGS*(U[35]-U[49]) + p->CGD*(U[35]-U[50]);
  G[49] -= p->CGS * U[35];
  G[50] -= p->CGD * U[35];

  /* Capacitive coupling: result node ANDOI-gate 2 === S0
   * Fortran G(49) -> C G[48] */
  G[48] += p->COUT * U[48];

  /* Capacitive coupling: result node ANDOI-gate 3
   * Fortran G(67) -> C G[66] */
  G[66] += p->CGS*(U[66]-U[116]) + p->CGD*(U[66]-U[117])
         + p->CGS*(U[66]-U[135]) + p->CGD*(U[66]-U[136]);
  G[116] -= p->CGS * U[66];
  G[117] -= p->CGD * U[66];
  G[135] -= p->CGS * U[66];
  G[136] -= p->CGD * U[66];

  /* Capacitive coupling: result node NOR-gate 3
   * Fortran G(85) -> C G[84] */
  G[84] += p->CGS*(U[84]-U[98]) + p->CGD*(U[84]-U[99])
         + p->CGS*(U[84]-U[148]) + p->CGD*(U[84]-U[149]);
  G[98]  -= p->CGS * U[84];
  G[99]  -= p->CGD * U[84];
  G[148] -= p->CGS * U[84];
  G[149] -= p->CGD * U[84];

  /* Capacitive coupling: result node ANDOI-gate 4
   * Fortran G(98) -> C G[97] */
  G[97] += p->CGS*(U[97]-U[121]) + p->CGD*(U[97]-U[122])
         + p->CGS*(U[97]-U[139]) + p->CGD*(U[97]-U[140])
         + p->CGS*(U[97]-U[157]) + p->CGD*(U[97]-U[158])
         + p->CGS*(U[97]-U[161]) + p->CGD*(U[97]-U[162]);
  G[121] -= p->CGS * U[97];
  G[122] -= p->CGD * U[97];
  G[139] -= p->CGS * U[97];
  G[140] -= p->CGD * U[97];
  G[157] -= p->CGS * U[97];
  G[158] -= p->CGD * U[97];

  /* Capacitive coupling: result NAND-gate
   * Fortran G(116) -> C G[115] */
  G[115] += p->CGS*(U[115]-U[130]) + p->CGD*(U[115]-U[131]);
  G[130] -= p->CGS * U[115];
  G[131] -= p->CGD * U[115];

  /* Capacitive coupling: result node ORANI-gate === S1
   * Fortran G(130) -> C G[129] */
  G[129] += p->COUT * U[129];

  /* Capacitive coupling: result ANDOI-gate 5 === Cinvers
   * Fortran G(148) -> C G[147] */
  G[147] += CBDBS(U[164]-U[147], p) * (U[147]-U[164]) + p->COUT * U[147];

  /* Charge function of three additional transistors
   * Fortran nodes 162--175 -> C 161--174 */
  G[161] += p->CGS * (U[161] - U[97]);
  G[162] += p->CGD * (U[162] - U[97]);
  G[163] += CBDBS(U[163]-U[165], p) * (U[163]-U[165]);
  G[164] += CBDBS(U[164]-U[147], p) * (U[164]-U[147]);
  G[165] += CBDBS(U[163]-U[165], p) * (U[165]-U[163])
          + CBDBS(U[169]-U[165], p) * (U[165]-U[169])
          + p->CLOAD * U[165];
  G[166] += p->CGS * (U[166] - U[17]);
  G[167] += p->CGD * (U[167] - U[17]);
  G[168] += CBDBS(U[168]-U[170], p) * (U[168]-U[170]);
  G[169] += CBDBS(U[169]-U[165], p) * (U[169]-U[165]);
  G[170] += CBDBS(U[168]-U[170], p) * (U[170]-U[168])
          + CBDBS(U[174]-U[170], p) * (U[170]-U[174])
          + p->CLOAD * U[170];
  G[171] += p->CGS * U[171];
  G[172] += p->CGD * U[172];
  G[173] += CBDBS(U[173], p) * U[173];
  G[174] += CBDBS(U[174]-U[170], p) * (U[174]-U[170]);
}
/* =========================================================================
 * RHS wrapper: f(t, y) for the full 350-dimensional system
 * y[0..174] = charges, y[175..349] = voltages
 *
 * f[0..174]   = FCN(175, t, U, res)   (current equations)
 * f[175..349] = y[i] - GCN(U)[i-175]  (algebraic: stored charge = g(U))
 * =========================================================================*/
static int rhs_tba(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  TbaParams* p = (TbaParams*)ud;
  sunrealtype* yv = N_VGetArrayPointer(y);
  sunrealtype* f  = N_VGetArrayPointer(yd);

  sunrealtype x[175], res[175], gres[175];
  int i, ierr = 0;

  /* Extract voltages: Fortran x(i)=y(i+175) -> C x[i]=yv[i+175] */
  for (i = 0; i < 175; i++) x[i] = yv[i + 175];

  /* Compute circuit currents */
  FCN(175, t, x, res, &ierr, p);
  if (ierr == -1) return 1;  /* recoverable */

  for (i = 0; i < 175; i++) f[i] = res[i];

  /* Compute charges from voltages */
  GCN(175, x, gres, p);

  /* Algebraic constraints: y[i] - g(U)[i-175] = 0 */
  for (i = 0; i < 175; i++) f[i + 175] = yv[i] - gres[i];

  return 0;
}

/* Mass matrix: M = diag(1,...,1,0,...,0) -- 175 ones, 175 zeros */
static int mas_tba(sunrealtype t, SUNMatrix M, void* ud,
                   N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)ud; (void)t1; (void)t2; (void)t3;
  SUNMatZero(M);
  for (int i = 0; i < 175; i++)
    SM_ELEMENT_D(M, i, i) = 1.0;
  return 0;
}
/* =========================================================================
 * Initialize parameters and initial conditions
 * =========================================================================*/
static void init_params(TbaParams* p)
{
  p->CTIME = 1.0e4;
  p->STIFF = 5.0;
  p->RGS   = 40.0 / (p->CTIME * p->STIFF);
  p->RGD   = 40.0 / (p->CTIME * p->STIFF);
  p->RBS   = 100.0 / (p->CTIME * p->STIFF);
  p->RBD   = 100.0 / (p->CTIME * p->STIFF);
  p->CGS   = 0.6e-4 * p->CTIME;
  p->CGD   = 0.6e-4 * p->CTIME;
  p->CBD   = 2.4e-5 * p->CTIME;
  p->CBS   = 2.4e-5 * p->CTIME;
  p->DELTA = 0.02;
  p->CURIS = 1.0e-15 * p->CTIME * p->STIFF;
  p->VTH   = 25.85;
  p->VDD   = 5.0;
  p->VBB   = -2.5;
  p->CLOAD = 0.0;
  p->COUT  = 2.0e-4 * p->CTIME - p->CLOAD;
}

static void init_voltages(sunrealtype* U)
{
  /* 175 initial voltage values (Fortran U(1)..U(175) -> C U[0]..U[174]) */
  U[  0] =  4.999999999996544e+00;
  U[  1] =  4.999999999999970e+00;
  U[  2] = -2.499999999999975e+00;
  U[  3] = -2.499999999999975e+00;
  U[  4] =  4.999999999996514e+00;
  U[  5] =  0.000000000000000e+00;
  U[  6] =  4.999999999996514e+00;
  U[  7] = -2.499999999999991e+00;
  U[  8] = -2.499999999999975e+00;
  U[  9] =  0.000000000000000e+00;
  U[ 10] =  4.999999999996514e+00;
  U[ 11] = -2.499999999999991e+00;
  U[ 12] = -2.499999999999975e+00;
  U[ 13] =  0.215858486765796e+00;
  U[ 14] =  4.988182208251953e+00;
  U[ 15] = -2.499999999999990e+00;
  U[ 16] = -2.499999999999975e+00;
  U[ 17] =  0.204040695017748e+00;
  U[ 18] =  0.011817791748026e+00;
  U[ 19] =  0.192222903269723e+00;
  U[ 20] = -2.499999999999991e+00;
  U[ 21] = -2.499999999999990e+00;
  U[ 22] = -0.228160951881239e+00;
  U[ 23] =  0.204040695017748e+00;
  U[ 24] = -2.499999999999992e+00;
  U[ 25] = -2.499999999999990e+00;
  U[ 26] = -0.228160951881241e+00;
  U[ 27] =  0.000000000000000e+00;
  U[ 28] = -0.228160951881239e+00;
  U[ 29] = -2.499999999999991e+00;
  U[ 30] = -2.499999999999992e+00;
  U[ 31] =  4.999999999996547e+00;
  U[ 32] =  4.999999999999970e+00;
  U[ 33] = -2.499999999999975e+00;
  U[ 34] = -2.499999999999975e+00;
  U[ 35] =  4.999999999996517e+00;
  U[ 36] =  0.000000000000000e+00;
  U[ 37] =  4.999999999996517e+00;
  U[ 38] = -2.499999999999991e+00;
  U[ 39] = -2.499999999999975e+00;
  U[ 40] =  0.000000000000000e+00;
  U[ 41] =  4.999999999996517e+00;
  U[ 42] = -2.499999999999991e+00;
  U[ 43] = -2.499999999999975e+00;
  U[ 44] =  0.215858484247529e+00;
  U[ 45] =  4.988182208251953e+00;
  U[ 46] = -2.499999999999990e+00;
  U[ 47] = -2.499999999999975e+00;
  U[ 48] =  0.204040692499482e+00;
  U[ 49] =  0.011817791748035e+00;
  U[ 50] =  0.192222900751447e+00;
  U[ 51] = -2.499999999999991e+00;
  U[ 52] = -2.499999999999990e+00;
  U[ 53] = -0.026041071738432e+00;
  U[ 54] =  0.204040692499482e+00;
  U[ 55] = -2.499999999999992e+00;
  U[ 56] = -2.499999999999990e+00;
  U[ 57] = -0.026041071738434e+00;
  U[ 58] =  0.000000000000000e+00;
  U[ 59] = -0.026041071738432e+00;
  U[ 60] = -2.499999999999991e+00;
  U[ 61] = -2.499999999999992e+00;
  U[ 62] =  0.215858484880918e+00;
  U[ 63] =  4.988182208251953e+00;
  U[ 64] = -2.499999999999990e+00;
  U[ 65] = -2.499999999999975e+00;
  U[ 66] =  0.204040693132870e+00;
  U[ 67] =  0.011817791748026e+00;
  U[ 68] =  0.192222901384845e+00;
  U[ 69] = -2.499999999999991e+00;
  U[ 70] = -2.499999999999990e+00;
  U[ 71] = -0.026041071737961e+00;
  U[ 72] =  0.204040693132870e+00;
  U[ 73] = -2.499999999999992e+00;
  U[ 74] = -2.499999999999990e+00;
  U[ 75] = -0.026041071737963e+00;
  U[ 76] =  0.000000000000000e+00;
  U[ 77] = -0.026041071737961e+00;
  U[ 78] = -2.499999999999991e+00;
  U[ 79] = -2.499999999999992e+00;
  U[ 80] =  4.999999999996546e+00;
  U[ 81] =  4.999999999999970e+00;
  U[ 82] = -2.499999999999975e+00;
  U[ 83] = -2.499999999999975e+00;
  U[ 84] =  4.999999999996516e+00;
  U[ 85] =  0.000000000000000e+00;
  U[ 86] =  4.999999999996516e+00;
  U[ 87] = -2.499999999999991e+00;
  U[ 88] = -2.499999999999975e+00;
  U[ 89] =  0.000000000000000e+00;
  U[ 90] =  4.999999999996516e+00;
  U[ 91] = -2.499999999999991e+00;
  U[ 92] = -2.499999999999975e+00;
  U[ 93] =  0.215858481060569e+00;
  U[ 94] =  4.988182208251953e+00;
  U[ 95] = -2.499999999999990e+00;
  U[ 96] = -2.499999999999975e+00;
  U[ 97] =  0.204040689312522e+00;
  U[ 98] =  0.011817791748023e+00;
  U[ 99] =  0.192222897564498e+00;
  U[100] = -2.499999999999991e+00;
  U[101] = -2.499999999999990e+00;
  U[102] =  4.734672533390068e+00;
  U[103] =  0.204040689312522e+00;
  U[104] = -2.499999999999977e+00;
  U[105] = -2.499999999999990e+00;
  U[106] =  4.734672533390062e+00;
  U[107] =  0.000000000000000e+00;
  U[108] =  4.734672533390068e+00;
  U[109] = -2.499999999999991e+00;
  U[110] = -2.499999999999977e+00;
  U[111] =  4.999999999996870e+00;
  U[112] =  4.999999999999972e+00;
  U[113] = -2.499999999999975e+00;
  U[114] = -2.499999999999975e+00;
  U[115] =  4.999999999996843e+00;
  U[116] = -0.025968303070038e+00;
  U[117] =  4.999999999996843e+00;
  U[118] = -2.499999999999992e+00;
  U[119] = -2.499999999999975e+00;
  U[120] = -0.025968303070040e+00;
  U[121] =  0.000000000000000e+00;
  U[122] = -0.025968303070038e+00;
  U[123] = -2.499999999999991e+00;
  U[124] = -2.499999999999992e+00;
  U[125] =  4.999999999997699e+00;
  U[126] =  4.999999999999980e+00;
  U[127] = -2.499999999999975e+00;
  U[128] = -2.499999999999975e+00;
  U[129] =  4.999999999997678e+00;
  U[130] =  4.744923533081106e+00;
  U[131] =  4.999999999997678e+00;
  U[132] = -2.499999999999977e+00;
  U[133] = -2.499999999999975e+00;
  U[134] =  4.744923533081098e+00;
  U[135] =  0.000000000000000e+00;
  U[136] =  4.744923533081106e+00;
  U[137] = -2.499999999999991e+00;
  U[138] = -2.499999999999977e+00;
  U[139] =  0.000000000000000e+00;
  U[140] =  4.744923533081106e+00;
  U[141] = -2.499999999999991e+00;
  U[142] = -2.499999999999977e+00;
  U[143] =  0.215858484844162e+00;
  U[144] =  4.988182208251953e+00;
  U[145] = -2.499999999999990e+00;
  U[146] = -2.499999999999975e+00;
  U[147] =  0.204040693096114e+00;
  U[148] =  0.011817791748023e+00;
  U[149] =  0.192222901348091e+00;
  U[150] = -2.499999999999991e+00;
  U[151] = -2.499999999999990e+00;
  U[152] =  0.204040693096045e+00;
  U[153] =  0.204040693096107e+00;
  U[154] = -2.499999999999990e+00;
  U[155] = -2.499999999999990e+00;
  U[156] =  0.204040693096037e+00;
  U[157] =  0.000000000000000e+00;
  U[158] =  0.204040693096037e+00;
  U[159] = -2.499999999999991e+00;
  U[160] = -2.499999999999990e+00;
  U[161] = -0.026017361873565e+00;
  U[162] =  0.204040693096114e+00;
  U[163] = -2.499999999999992e+00;
  U[164] = -2.499999999999990e+00;
  U[165] = -0.026017361873568e+00;
  U[166] = -0.026017590106916e+00;
  U[167] = -0.026017361873565e+00;
  U[168] = -2.499999999999992e+00;
  U[169] = -2.499999999999992e+00;
  U[170] = -0.026017590106918e+00;
  U[171] =  0.000000000000000e+00;
  U[172] = -0.026017590106916e+00;
  U[173] = -2.499999999999991e+00;
  U[174] = -2.499999999999992e+00;
}
/* =========================================================================
 * Main
 * =========================================================================*/
int main(int argc, char* argv[])
{
  sunrealtype rtol = 1.0e-5;
  sunrealtype atol_val = 1.0e-5;
  sunrealtype h0 = 4.0e-5;
  int use_schur = 0;
  int nsmin        = 3;
  int nsmax        = 7;
  if (argc > 1) rtol     = atof(argv[1]);
  if (argc > 2) atol_val = atof(argv[2]);
  if (argc > 3) h0       = atof(argv[3]);
  if (argc > 4) use_schur = atoi(argv[4]);
  if (argc > 5) nsmin     = atoi(argv[5]);
  if (argc > 6) nsmax     = atoi(argv[6]);

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  int n = 350;
  TbaParams params;
  init_params(&params);

  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, nsmin, nsmax);

  /* Initial conditions */
  N_Vector y0 = N_VNew_Serial(n, sunctx);
  sunrealtype* y0d = N_VGetArrayPointer(y0);

  sunrealtype Uinit[175], Ginit[175];
  init_voltages(Uinit);
  GCN(175, Uinit, Ginit, &params);

  /* y[0..174] = charges g(U), y[175..349] = voltages U */
  int i;
  for (i = 0; i < 175; i++) y0d[i] = Ginit[i];
  for (i = 0; i < 175; i++) y0d[i + 175] = Uinit[i];

  Radau5Init(mem, rhs_tba, 0.0, y0);

  SUNMatrix Jt = SUNDenseMatrix(n, n, sunctx);
  Radau5SetLinearSolver(mem, Jt, NULL);
  Radau5SetSchurDecomp(mem, use_schur);
  /* No analytic Jacobian -- use DQ */

  /* Mass matrix */
  SUNMatrix Mt = SUNDenseMatrix(n, n, sunctx);
  Radau5SetMassFn(mem, mas_tba, Mt);

  /* DAE index: all 350 variables are index-1 */
  Radau5SetDAEIndex(mem, n, 0, 0);

  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);
  Radau5SetUserData(mem, &params);
  Radau5SetMaxNumSteps(mem, 1000000);

  N_Vector yout = N_VNew_Serial(n, sunctx);
  sunrealtype tret;
  int ret;

  /* Segmented solve at discontinuity points t=5,10,...,315
   * (matching Fortran driver which calls radau5 64 times with h=h0 reset).
   * The input signals have sharp transitions at t=5*k, so the solver must
   * restart with a small step at each boundary. */
  ret = 0;
  for (int seg = 1; seg <= 63; seg++) {
    sunrealtype t_disc = 5.0 * seg;
    /* Reset solver state at discontinuity, matching Fortran radau5 entry:
     * h=h0, FIRST=.TRUE., REJECT=.FALSE., NSING=0, FACCON=1.0
     * FIRST is automatically cleared after the first accepted step,
     * so extrapolation is used for subsequent steps within the segment. */
    Radau5ResetForDiscontinuity(mem, h0);
    ret = Radau5Solve(mem, t_disc, yout, &tret);
    if (ret < 0) {
      printf("Radau5Solve failed at segment %d (t=%.1f), ret=%d\n",
             seg, t_disc, ret);
      break;
    }
  }

  /* Reference solution at t=315: check 3 output signals
   * Fortran indices: y(224), y(305), y(323) -> C indices: 223, 304, 322 */
  sunrealtype yref_224 = 0.2040419147264534e+00;  /* S0: node 49 */
  sunrealtype yref_305 = 0.4997238455712048e+01;  /* S1: node 130 */
  sunrealtype yref_323 = 0.2038985905095614e+00;  /* C:  node 148 */

  sunrealtype* yd = N_VGetArrayPointer(yout);

  printf("=== Two Bit Adding Unit (rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n",
         rtol, atol_val, h0, use_schur);
  printf("ret=%d, tret=%.6e\n", ret, tret);

  sunrealtype err1 = fabs(yd[223] - yref_224);
  sunrealtype rel1 = err1 / fabs(yref_224);
  sunrealtype err2 = fabs(yd[304] - yref_305);
  sunrealtype rel2 = err2 / fabs(yref_305);
  sunrealtype err3 = fabs(yd[322] - yref_323);
  sunrealtype rel3 = err3 / fabs(yref_323);

  printf("y[223] (S0) = %16.10e  ref = %16.10e  rel_err = %.3e\n",
         yd[223], yref_224, rel1);
  printf("y[304] (S1) = %16.10e  ref = %16.10e  rel_err = %.3e\n",
         yd[304], yref_305, rel2);
  printf("y[322] (C)  = %16.10e  ref = %16.10e  rel_err = %.3e\n",
         yd[322], yref_323, rel3);

  sunrealtype maxrelerr = rel1;
  if (rel2 > maxrelerr) maxrelerr = rel2;
  if (rel3 > maxrelerr) maxrelerr = rel3;
  printf("max_rel_err = %.3e\n", maxrelerr);

  long int nstep, naccpt, nrejct, nfcn, njac, ndec;
  Radau5GetNumSteps(mem, &nstep);
  Radau5GetNumAccSteps(mem, &naccpt);
  Radau5GetNumRejSteps(mem, &nrejct);
  Radau5GetNumRhsEvals(mem, &nfcn);
  Radau5GetNumJacEvals(mem, &njac);
  Radau5GetNumDecomps(mem, &ndec);
  printf("nstep=%ld naccpt=%ld nrejct=%ld nfcn=%ld njac=%ld ndec=%ld\n",
         nstep, naccpt, nrejct, nfcn, njac, ndec);

  int pass = (ret >= 0 && maxrelerr < 0.1);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout);
  SUNMatDestroy(Jt); SUNMatDestroy(Mt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
