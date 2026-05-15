/* ---------------------------------------------------------------------------
 * radau5_crank.c — Slider Crank mechanism (DAE index-2, n=24)
 *
 * Flexible slider crank with 3 rigid + 4 elastic DOFs.
 * M*y' = f(t,y), M = diag(1..1, 0..0) (first 14 differential)
 * DAE index: nind1=14, nind2=10, nind3=0
 *
 * State vector (0-based):
 *   y[0..6]:   positions (phi1, phi2, x3, q1..q4)
 *   y[7..13]:  velocities
 *   y[14..20]: accelerations
 *   y[21..23]: Lagrange multipliers
 *
 * t in [0, 0.1]
 * Reference: Simeon, B.: Modelling a flexible slider crank mechanism
 *            by a mixed system of DAEs and PDEs,
 *            Math. Modelling of Systems 2, 1-18 (1996)
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

#define NEQ 24

/* Physical parameters */
#define COMEGA 150.0
#define CL1    0.15
#define CL2    0.30
#define CM1    0.36
#define CM2    0.151104
#define CM3    0.075552
#define CJ1    0.002727
#define CJ2    0.0045339259
#define CEE    2.0e11
#define CBB    0.008
#define CHH    0.008
#define CRHO   7870.0
#define CGRAV  0.0
#define CPI    3.1415927

#define CNQ 4
#define CNP 7

/* ---------------------------------------------------------------------------
 * FE matrices — computed once at first call (mirrors Fortran SAVE)
 * Stored in row-major C convention: mat[row][col]
 * Fortran MQ(J,I) column-major → C feMQ[j][i] means same element
 * ---------------------------------------------------------------------------*/
static int fe_init = 0;
static sunrealtype feMQ[CNQ][CNQ], feKQ[CNQ][CNQ];
static sunrealtype feBQ[CNQ][CNQ], feDQ[CNQ][CNQ];
static sunrealtype fec1[CNQ], fec2[CNQ], fec12[CNQ], fec21[CNQ];

static void init_fe(void)
{
  int i, j;
  sunrealtype FACM = CRHO * CBB * CHH * CL2;
  sunrealtype FACK = CEE * CBB * CHH / CL2;
  sunrealtype FACB = CBB * CHH * CL2;
  sunrealtype PI2  = CPI * CPI;
  sunrealtype PI3  = PI2 * CPI;
  sunrealtype PI4  = PI2 * PI2;
  sunrealtype HoL  = CHH / CL2;

  for (i = 0; i < CNQ; i++) {
    for (j = 0; j < CNQ; j++) {
      feMQ[i][j] = 0.0; feKQ[i][j] = 0.0;
      feBQ[i][j] = 0.0; feDQ[i][j] = 0.0;
    }
    fec1[i] = 0.0; fec2[i] = 0.0;
    fec12[i] = 0.0; fec21[i] = 0.0;
  }

  /* Mass matrix MQ (Fortran indices 1-based → C 0-based) */
  feMQ[0][0] = FACM * 0.5;
  feMQ[1][1] = FACM * 0.5;
  feMQ[2][2] = FACM * 8.0;
  feMQ[2][3] = FACM * 1.0;
  feMQ[3][2] = FACM * 1.0;
  feMQ[3][3] = FACM * 2.0;

  /* Stiffness matrix KQ */
  feKQ[0][0] = FACK * PI4 / 24.0 * (HoL * HoL);
  feKQ[1][1] = FACK * PI4 * 2.0 / 3.0 * (HoL * HoL);
  feKQ[2][2] = FACK * 16.0 / 3.0;
  feKQ[2][3] = -FACK * 8.0 / 3.0;
  feKQ[3][2] = -FACK * 8.0 / 3.0;
  feKQ[3][3] = FACK * 7.0 / 3.0;
  /* Coupling matrix BQ
   * Fortran BQ(i,j) col-major stored as BQ(NQMAX,NQMAX)
   * In Fortran: BQ(1,3), BQ(1,4), BQ(2,4), BQ(3,1), BQ(4,1), BQ(4,2)
   * Fortran BQ(row,col) → C feBQ[row-1][col-1] */
  feBQ[0][2] = -FACB * 16.0 / PI3;
  feBQ[0][3] = FACB * (8.0 / PI3 - 1.0 / CPI);
  feBQ[1][3] = FACB * 0.5 / CPI;
  feBQ[2][0] = FACB * 16.0 / PI3;
  feBQ[3][0] = -FACB * (8.0 / PI3 - 1.0 / CPI);
  feBQ[3][1] = -FACB * 0.5 / CPI;

  /* Vectors c1, c2, c12, c21 (Fortran 1-based → C 0-based) */
  fec1[2]  = FACB * 2.0 / 3.0;
  fec1[3]  = FACB * 1.0 / 6.0;
  fec2[0]  = FACB * 2.0 / CPI;
  fec12[2] = CL2 * FACB * 1.0 / 3.0;
  fec12[3] = CL2 * FACB * 1.0 / 6.0;
  fec21[0] = CL2 * FACB * 1.0 / CPI;
  fec21[1] = -CL2 * FACB * 0.5 / CPI;

  /* No damping (ipar(2)=0) */
  fe_init = 1;
}

/* Simple dot product */
static sunrealtype dot(int n, const sunrealtype* a, const sunrealtype* b)
{
  sunrealtype s = 0.0;
  for (int i = 0; i < n; i++) s += a[i] * b[i];
  return s;
}
/* ---------------------------------------------------------------------------
 * Column extraction helper: get column 'col' of a CNQ x CNQ matrix
 * ---------------------------------------------------------------------------*/
static void get_col(const sunrealtype mat[CNQ][CNQ], int col, sunrealtype out[CNQ])
{
  for (int i = 0; i < CNQ; i++) out[i] = mat[i][col];
}

/* ---------------------------------------------------------------------------
 * RHS function: implements the feval wrapper from crank.f
 *
 * feval calls RESMBS(ityp=1, iequa=0, icall=0, t, y, dy=0, f, ...)
 * then negates f[0..13].
 *
 * RESMBS with dy=0:
 *   delta[i]     = 0 - x[np+i]       = -velocity        (i=0..6)
 *   delta[np+i]  = 0 - x[2*np+i]     = -acceleration    (i=0..6)
 *   delta[2*np+i] = AM_col_i . w - F[i] + GP^T . lambda (i=0..6)
 *   delta[nx-3..nx-1] = VLC[0..2]    (index-2: velocity constraints)
 *
 * After negation of delta[0..13]:
 *   f[i]    = velocity[i]             (i=0..6)
 *   f[np+i] = acceleration[i]        (i=0..6)
 *   f[2*np+i] = AM_col_i . w - F[i] + GP^T . lambda  (unchanged)
 *   f[nx-3..nx-1] = VLC              (unchanged)
 * ---------------------------------------------------------------------------*/
static int rhs_crank(sunrealtype t, N_Vector yv, N_Vector ydv, void* ud)
{
  (void)ud;
  if (!fe_init) init_fe();

  sunrealtype* x = N_VGetArrayPointer(yv);
  sunrealtype* delta = N_VGetArrayPointer(ydv);
  int i, j;

  /* Trig functions */
  sunrealtype cosp1  = cos(x[0]);
  sunrealtype cosp2  = cos(x[1]);
  sunrealtype sinp1  = sin(x[0]);
  sunrealtype sinp2  = sin(x[1]);
  sunrealtype cosp12 = cos(x[0] - x[1]);
  sunrealtype sinp12 = sin(x[0] - x[1]);

  /* Velocities of phi1, phi2 */
  sunrealtype v1 = x[CNP];      /* x[7] = phi1_dot */
  sunrealtype v2 = x[CNP + 1];  /* x[8] = phi2_dot */
  /* Elastic coordinates q[0..3] and their velocities qd[0..3]
   * Fortran: Q(I) = X(3+I), QD(I) = X(NP+3+I) for I=1..NQ
   * C: q[i] = x[3+i], qd[i] = x[CNP+3+i] for i=0..3 */
  sunrealtype q[CNQ], qd[CNQ];
  for (i = 0; i < CNQ; i++) {
    q[i]  = x[3 + i];       /* x[3..6] */
    qd[i] = x[CNP + 3 + i]; /* x[10..13] */
  }

  /* Scalar products with FE vectors */
  sunrealtype c1TQ   = dot(CNQ, fec1, q);
  sunrealtype c1TQD  = dot(CNQ, fec1, qd);
  sunrealtype c2TQ   = dot(CNQ, fec2, q);
  sunrealtype c2TQD  = dot(CNQ, fec2, qd);
  sunrealtype c12TQ  = dot(CNQ, fec12, q);
  sunrealtype c12TQD = dot(CNQ, fec12, qd);

  /* Matrix-vector products:
   * Fortran MQQ(I) = PDOT(NQ, MQ(1,I), 1, Q, 1) = column I of MQ dotted with Q
   * In Fortran col-major: MQ(1,I)..MQ(NQ,I) is column I
   * Our feMQ[row][col]: column I = feMQ[0][I], feMQ[1][I], ...
   * So MQQ[i] = sum_j feMQ[j][i] * q[j] */
  sunrealtype MQQ[CNQ], KQQ[CNQ], DQQD[CNQ], QTBQ[CNQ], BQQD[CNQ];
  sunrealtype col_tmp[CNQ];

  for (i = 0; i < CNQ; i++) {
    get_col(feMQ, i, col_tmp);
    MQQ[i] = dot(CNQ, col_tmp, q);

    get_col(feKQ, i, col_tmp);
    KQQ[i] = dot(CNQ, col_tmp, q);

    get_col(feDQ, i, col_tmp);
    DQQD[i] = dot(CNQ, col_tmp, qd);

    /* QTBQ(I) = PDOT(NQ, Q, 1, BQ(1,I), 1) = Q^T * column I of BQ */
    get_col(feBQ, i, col_tmp);
    QTBQ[i] = dot(CNQ, q, col_tmp);

    /* BQQD(I) = PDOT(NQ, BQ(I,1), NQMAX, QD, 1)
     * This is row I of BQ dotted with QD (stride NQMAX in Fortran = row stride)
     * In C row-major: feBQ[i][0..3] is row i */
    BQQD[i] = dot(CNQ, feBQ[i], qd);
  }

  sunrealtype QTMQQ   = dot(CNQ, q, MQQ);
  sunrealtype QDTMQQ  = dot(CNQ, qd, MQQ);
  sunrealtype QDTBQQD = dot(CNQ, qd, BQQD);
  /* Kinematic equations (with dy=0):
   * delta[i]     = 0 - x[CNP+i]   = -velocity     (i=0..6)
   * delta[CNP+i] = 0 - x[2*CNP+i] = -acceleration (i=0..6) */
  for (i = 0; i < CNP; i++) {
    delta[i]       = -x[CNP + i];
    delta[CNP + i] = -x[2 * CNP + i];
  }

  /* Compute inertia matrix AM[CNP][CNP] (0-based)
   * Fortran AM(I,J) col-major → C AM[i][j] row-major, same indexing */
  sunrealtype AM[CNP][CNP];
  for (i = 0; i < CNP; i++)
    for (j = 0; j < CNP; j++)
      AM[i][j] = 0.0;

  AM[0][0] = CJ1 + CM2 * CL1 * CL1;
  AM[0][1] = 0.5 * CL1 * CL2 * CM2 * cosp12;
  AM[1][1] = CJ2;
  AM[0][2] = 0.0;
  AM[1][2] = 0.0;
  AM[2][0] = 0.0;
  AM[2][1] = 0.0;
  AM[2][2] = CM3;

  /* Flexible coupling to AM(1,2) and AM(2,2) */
  AM[0][1] += CRHO * CL1 * (sinp12 * c2TQ + cosp12 * c1TQ);
  AM[1][1] += QTMQQ + 2.0 * CRHO * c12TQ;

  /* AM(1, 3+I) and AM(2, 3+I) for I=0..3 (Fortran I=1..NQ) */
  for (i = 0; i < CNQ; i++) {
    AM[0][3 + i] = CRHO * CL1 * (-sinp12 * fec1[i] + cosp12 * fec2[i]);
    AM[1][3 + i] = CRHO * fec21[i] + CRHO * QTBQ[i];
    AM[2][3 + i] = 0.0;
  }

  /* MQ block: AM(3+J, 3+I) = MQ(J,I) for J=1..I, I=1..NQ
   * Fortran: upper triangle of MQ placed in AM(3+J, 3+I)
   * C 0-based: AM[3+j][3+i] = feMQ[j][i] for j <= i */
  for (i = 0; i < CNQ; i++)
    for (j = 0; j <= i; j++)
      AM[3 + j][3 + i] = feMQ[j][i];

  /* Symmetrize: AM[j][i] = AM[i][j] for j > i */
  for (i = 0; i < CNP; i++)
    for (j = i + 1; j < CNP; j++)
      AM[j][i] = AM[i][j];
  /* Compute constraint matrix GP[3][CNP]
   * Fortran GP(constraint, dof) → C GP[constraint][dof]
   * KU=4 in Fortran (1-based) → Q index 3 (0-based), so QKU = q[3]
   * KV=0 in Fortran → QKV = 0 */
  sunrealtype QKU = q[3];   /* Fortran KU=4 → C index 3 */
  sunrealtype QKV = 0.0;    /* Fortran KV=0 → no lateral tip */
  sunrealtype GP[3][CNP];

  GP[0][0] = CL1 * cosp1;
  GP[0][1] = CL2 * cosp2 + QKU * cosp2 - QKV * sinp2;
  GP[0][2] = 0.0;
  GP[1][0] = CL1 * sinp1;
  GP[1][1] = CL2 * sinp2 + QKU * sinp2 + QKV * cosp2;
  GP[1][2] = 1.0;
  GP[2][0] = 1.0;
  GP[2][1] = 0.0;
  GP[2][2] = 0.0;

  for (i = 0; i < CNQ; i++) {
    GP[0][3 + i] = 0.0;
    GP[1][3 + i] = 0.0;
    GP[2][3 + i] = 0.0;
  }
  /* KU != 0: GP(1, 3+KU) = sinp2, GP(2, 3+KU) = -cosp2
   * Fortran KU=4 → GP(*,3+4)=GP(*,7) → C GP[*][6] since dof index 3+3=6 */
  GP[0][6] = sinp2;
  GP[1][6] = -cosp2;
  /* KV=0: no contribution */

  /* Forces — rigid motion entries (Fortran F(1..NP), C F[0..6]) */
  sunrealtype F[CNP];
  F[0] = -0.5 * CL1 * CGRAV * (CM1 + 2.0 * CM2) * cosp1
         - 0.5 * CL1 * CL2 * CM2 * v2 * v2 * sinp12;
  F[1] = -0.5 * CL2 * CGRAV * CM2 * cosp2
         + 0.5 * CL1 * CL2 * CM2 * v1 * v1 * sinp12;
  F[2] = 0.0;

  /* Superposition of flexible motion (term f^e) */
  F[0] += CRHO * CL1 * v2 * v2 * (-sinp12 * c1TQ + cosp12 * c2TQ)
        - 2.0 * CRHO * CL1 * v2 * (cosp12 * c1TQD + sinp12 * c2TQD);
  F[1] += CRHO * CL1 * v1 * v1 * (sinp12 * c1TQ - cosp12 * c2TQ)
        - 2.0 * CRHO * v2 * c12TQD - 2.0 * v2 * QDTMQQ
        - CRHO * QDTBQQD - CRHO * CGRAV * (cosp2 * c1TQ - sinp2 * c2TQ);

  /* Coriolis and gravity terms for elastic DOFs (Gamma) */
  for (i = 0; i < CNQ; i++) {
    F[3 + i] = v2 * v2 * MQQ[i]
      + CRHO * (v2 * v2 * fec12[i]
                + CL1 * v1 * v1 * (cosp12 * fec1[i] + sinp12 * fec2[i])
                + 2.0 * v2 * BQQD[i])
      - CRHO * CGRAV * (sinp2 * fec1[i] + cosp2 * fec2[i]);
  }

  /* Stiffness + damping: F(3+I) -= KQQ(I) + DQQD(I) */
  for (i = 0; i < CNQ; i++)
    F[3 + i] -= KQQ[i] + DQQD[i];

  /* ipar(1)=0: no nonlinear stiffness term */

  /* Dynamics part II: delta[2*CNP+i] = AM[*][i] . w - F[i] + GP^T . lambda
   * Fortran: DELTA(2*NP+I) = PDOT(NP, AM(1,I), 1, X(2*NP+1), 1)
   *          - F(I) + GP(1,I)*X(NX-2) + GP(2,I)*X(NX-1) + GP(3,I)*X(NX)
   * Fortran AM(1,I) with stride 1 = column I of AM
   * C: AM[*][i] = column i → we need dot of column i with w[0..6]
   * w = x[2*CNP .. 2*CNP+6] = x[14..20]
   * lambda = x[21..23] */
  sunrealtype* w = &x[2 * CNP];       /* x[14..20] */
  sunrealtype lam0 = x[21], lam1 = x[22], lam2 = x[23];

  for (i = 0; i < CNP; i++) {
    sunrealtype am_dot_w = 0.0;
    for (j = 0; j < CNP; j++)
      am_dot_w += AM[j][i] * w[j];
    delta[2 * CNP + i] = am_dot_w - F[i]
                        + GP[0][i] * lam0 + GP[1][i] * lam1 + GP[2][i] * lam2;
  }

  /* Velocity level constraints (index-2: ITYP=1)
   * VLC(k) = sum_i GP(k,i) * x[CNP+i] for k=0..2, plus VLC(3) -= OMEGA
   * Fortran: VLC(1)=0, VLC(2)=0, VLC(3)=-OMEGA
   *          VLC(k) += sum_{I=1}^{NP} GP(k,I)*X(NP+I) */
  sunrealtype VLC[3] = {0.0, 0.0, -COMEGA};
  for (i = 0; i < CNP; i++) {
    VLC[0] += GP[0][i] * x[CNP + i];
    VLC[1] += GP[1][i] * x[CNP + i];
    VLC[2] += GP[2][i] * x[CNP + i];
  }

  delta[21] = VLC[0];
  delta[22] = VLC[1];
  delta[23] = VLC[2];

  /* feval negation: f(i) = -f(i) for i=1..14 (Fortran) → i=0..13 (C) */
  for (i = 0; i < 14; i++)
    delta[i] = -delta[i];

  return 0;
}
/* ---------------------------------------------------------------------------
 * Mass matrix: M = diag(1..1, 0..0), first 14 differential
 * ---------------------------------------------------------------------------*/
static int mas_crank(sunrealtype t, SUNMatrix M, void* ud,
                     N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)ud; (void)t1; (void)t2; (void)t3;
  SUNMatZero(M);
  for (int i = 0; i < 14; i++)
    SM_ELEMENT_D(M, i, i) = 1.0;
  return 0;
}

int main(int argc, char* argv[])
{
  sunrealtype rtol = 1.0e-6;
  sunrealtype atol_val = 1.0e-6;
  sunrealtype h0   = 1.0e-10;
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

  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, nsmin, nsmax);
  N_Vector y0 = N_VNew_Serial(NEQ, sunctx);
  sunrealtype* y0v = N_VGetArrayPointer(y0);

  /* Initial conditions from crank.f init() */
  /* Position variables: phi1, phi2, x3, q1..q4 */
  y0v[0]  =  0.0;
  y0v[1]  =  0.0;
  y0v[2]  =  4.500169330000000e-01;
  y0v[3]  =  0.0;
  y0v[4]  =  0.0;
  y0v[5]  =  1.033398630000000e-05;
  y0v[6]  =  1.693279690000000e-05;
  /* Velocity variables */
  y0v[7]  =  1.500000000000000e+02;
  y0v[8]  = -7.499576703969453e+01;
  y0v[9]  = -2.689386719979040e-06;
  y0v[10] =  4.448961125815990e-01;
  y0v[11] =  4.634339319238670e-03;
  y0v[12] = -1.785910760000550e-06;
  y0v[13] = -2.689386719979040e-06;
  /* Acceleration variables */
  y0v[14] =  0.0;
  y0v[15] = -1.344541576709835e-03;
  y0v[16] = -5.062194924490193e+03;
  y0v[17] = -6.829725665986310e-05;
  y0v[18] =  1.813207639590617e-20;
  y0v[19] = -4.268463266810281e+00;
  y0v[20] =  2.098339029337557e-01;
  /* Lagrange multipliers */
  y0v[21] = -6.552727150584648e-08;
  y0v[22] =  3.824589509350831e+02;
  y0v[23] = -4.635908708561371e-09;

  Radau5Init(mem, rhs_crank, 0.0, y0);

  SUNMatrix Jt = SUNDenseMatrix(NEQ, NEQ, sunctx);
  Radau5SetLinearSolver(mem, Jt, NULL);
  Radau5SetSchurDecomp(mem, use_schur);
  /* DQ Jacobian (no analytic Jac) */

  SUNMatrix Mt = SUNDenseMatrix(NEQ, NEQ, sunctx);
  Radau5SetMassFn(mem, mas_crank, Mt);
  Radau5SetDAEIndex(mem, 14, 10, 0);

  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);

  N_Vector yout = N_VNew_Serial(NEQ, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 0.1, yout, &tret);

  /* Reference solution at t=0.1 from crank.f solut() */
  static const sunrealtype yref[NEQ] = {
     0.1500000000000104e+02, -0.3311734988256260e+00,
     0.1697373328427860e+00,  0.1893192899613509e-03,
     0.2375751249879174e-04, -0.5323896770569702e-05,
    -0.8363313279112129e-05,  0.1500000000000000e+03,
     0.6025346755138369e+02, -0.8753116326670527e+01,
    -0.3005541400289738e-01, -0.5500431812571696e-02,
     0.4974111734266989e-03,  0.1105560003626645e-02,
     0.0,                     0.6488737541276957e+04,
     0.2167938629509884e+04,  0.3391137060286523e+02,
     0.1715134772216488e+00, -0.1422449408912512e+01,
     0.1003946428124810e+01, -0.6232935833287916e+02,
    -0.1637920993367306e+03,  0.2529857947066878e+02
  };

  sunrealtype* yd = N_VGetArrayPointer(yout);
  printf("=== Slider Crank (rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n",
         rtol, atol_val, h0, use_schur);
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

  /* Pass criterion: positions y[0..6] with max relative error < 1e-2 */
  int pass = (ret == 0 && maxrelerr < 1.0e-2);
  printf("max_pos_rel_err=%.3e\n", maxrelerr);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout);
  SUNMatDestroy(Jt); SUNMatDestroy(Mt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
