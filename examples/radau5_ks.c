/* ---------------------------------------------------------------------------
 * radau5_ks.c — Kuramoto-Sivashinsky equation in Fourier space (IVPtestset)
 *
 * d/dt u = -(dx^2 + dx^4)u - u*u'
 *
 * Formulated as Fourier modes:
 *   dt u_n = (n^2*q^2 - n^4*q^4)*u_n - (i*q*n/2) sum_{n1+n2=n} u_n1*u_n2
 *
 * Parameters: MMM=9, NH=512, N=1024, ND=1022 (modes 1..NH-1, real+imag)
 *             q=0.025, tend=100
 *
 * Jacobian is diagonal (only the linear stiffness term).
 * Band matrix with mu=0, ml=0. DQ Jacobian would also work but analytic
 * is trivial here.
 *
 * Reference: IVPtestset 2.4, KS problem
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_band.h>
#include "radau5.h"

#define MMM   9
#define NH    (1 << MMM)    /* 512 */
#define NFFT  (2 * NH)      /* 1024 */
#define NEQ   (NFFT - 2)    /* 1022 */

static const sunrealtype QQ = 0.025;

/* Global: the zero-mode value (constant through integration) */
static sunrealtype UZERO;

/* FFT twiddle factor storage */
#define MAXN_FFT 8192
static double cstore[2 * MAXN_FFT];
static int ntbl = 0;

/* ---------------------------------------------------------------------------
 * FFT: Cooley-Tukey radix-2 (translated from Fortran fft subroutine)
 * x is a complex array stored as real/imag pairs: x[2*n]
 * isign >= 0: forward, isign < 0: backward
 * ---------------------------------------------------------------------------*/
static void fft(double* x, int n, int isign)
{
  int i, j, k, m, jj, ii, iir, ijr;
  double pi, temp, vr1, vr2, wr1, wr2, xr1, xr2;

  if (n > MAXN_FFT) { fprintf(stderr, "n too big in fft\n"); exit(1); }

  /* Build twiddle table if needed */
  if (n > ntbl) {
    ntbl = n;
    pi = 4.0 * atan(1.0);
    j = 1;
    int icnt = 0;
    while (j < n) {
      for (k = 0; k < j; k++) {
        icnt++;
        cstore[2 * icnt - 2] = cos(k * (pi / j));
        cstore[2 * icnt - 1] = sin(k * (pi / j));
      }
      j *= 2;
    }
  }

  /* Bit reversal */
  j = 1;
  for (i = 1; i <= n; i++) {
    if (i < j) {
      jj = 2 * j; ii = 2 * i;
      temp = x[jj - 2]; x[jj - 2] = x[ii - 2]; x[ii - 2] = temp;
      temp = x[jj - 1]; x[jj - 1] = x[ii - 1]; x[ii - 1] = temp;
    }
    m = n / 2;
    while (j > m && m >= 1) { j -= m; m /= 2; }
    j += m;
  }

  /* Butterfly passes */
  j = 1;
  int icnt = 0;
  while (j < n) {
    jj = 2 * j;
    for (k = 1; k <= j; k++) {
      icnt++;
      wr1 = cstore[2 * icnt - 2];
      wr2 = cstore[2 * icnt - 1];
      if (isign < 0) wr2 = -wr2;

      for (i = k; i <= n; i += jj) {
        iir = 2 * i;
        ijr = 2 * (i + j);
        xr1 = x[ijr - 2];
        xr2 = x[ijr - 1];
        vr1 = wr1 * xr1 - wr2 * xr2;
        temp = x[iir - 2];
        x[iir - 2] = temp + vr1;
        x[ijr - 2] = temp - vr1;
        vr2 = wr1 * xr2 + wr2 * xr1;
        temp = x[iir - 1];
        x[ijr - 1] = temp - vr2;
        x[iir - 1] = temp + vr2;
      }
    }
    j = jj;
  }
}

/* ---------------------------------------------------------------------------
 * REALFT: Real FFT wrapper (translated from Fortran realft subroutine)
 * isign=+1: forward (real -> half-complex), isign=-1: backward
 * ---------------------------------------------------------------------------*/
static void realft(double* x, int n, int isign)
{
  double c1 = 0.5;
  double pi = 4.0 * atan(1.0);
  double theta = pi / n;
  double c2, wpr, wpi, wr, wi, wtemp;
  double h1r, h2i, uu, vv, h2r, h1i;
  double x1, x3, x2, x4, wr0, wr1v, wr2v;
  int i, i1, i3, n2p3;

  if (isign == 1) {
    c2 = -0.5;
    fft(x, n, 1);
  } else {
    c2 = 0.5;
    theta = -theta;
  }

  wtemp = sin(0.5 * theta);
  wpr = -2.0 * wtemp * wtemp;
  wpi = sin(theta);
  wr = 1.0 + wpr;
  wi = wpi;
  n2p3 = 2 * n + 3;

  for (i = 2; i <= n / 2; i++) {
    /* Optimized path (ppp=1 in Fortran) */
    i1 = 2 * i - 1 - 1;  /* 0-based index for x[i1] = Fortran x(2*i-1) */
    i3 = (n2p3 - 1) - (2 * i - 1) - 1; /* 0-based for x(i3) = Fortran x(n2p3-2*i) */

    x1 = x[i1];
    x3 = x[i3];
    h2i = c2 * (x1 - x3);
    uu = wi * h2i;
    vv = wr * h2i;
    h1r = c1 * (x1 + x3);
    x2 = x[i1 + 1];
    x4 = x[i3 + 1];
    h2r = -c2 * (x2 + x4);
    uu = wr * h2r - uu;
    vv = wi * h2r + vv;
    x[i1] = h1r + uu;
    x[i3] = h1r - uu;
    h1i = c1 * (x2 - x4);
    x[i1 + 1] = h1i + vv;
    x[i3 + 1] = -h1i + vv;
    wr0 = wi * wpi;
    wr1v = (1.0 + wpr);
    wi = wi * wr1v;
    wtemp = wr;
    wr2v = wtemp * wr1v;
    wi = wtemp * wpi + wi;
    wr = wr2v - wr0;
  }

  if (isign == 1) {
    h1r = x[0];
    x[0] = h1r + x[1];
    x[1] = h1r - x[1];
  } else {
    h1r = x[0];
    x[0] = c1 * (h1r + x[1]);
    x[1] = c1 * (h1r - x[1]);
    fft(x, n, -1);
  }
}

/* ---------------------------------------------------------------------------
 * RHS: translates Fortran FKS exactly
 *
 * The nonlinear term u*u' is computed via:
 *   1. Reconstruct full Fourier array U (with UZERO and 0 for mode 0)
 *   2. Inverse FFT to physical space, multiply by 2
 *   3. Square: U[i] = U[i]^2 / 2
 *   4. Forward FFT back, normalize by N
 *   5. Linear diagonal + nonlinear coupling
 * ---------------------------------------------------------------------------*/
static int rhs_ks(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)t; (void)ud;
  sunrealtype* yv = N_VGetArrayPointer(y);
  sunrealtype* f  = N_VGetArrayPointer(yd);

  double U[NFFT];

  /* Copy y to U: U[0]=UZERO, U[1]=0, U[2..N-1]=y[0..ND-1] */
  U[0] = UZERO;
  U[1] = 0.0;
  for (int i = 2; i < NFFT; i++)
    U[i] = yv[i - 2];

  /* Inverse FFT and multiply by 2 */
  realft(U, NH, -1);
  for (int i = 0; i < NFFT; i++)
    U[i] = 2.0 * U[i];

  /* Square: u^2/2 */
  for (int i = 0; i < NFFT; i++)
    U[i] = U[i] * U[i] / 2.0;

  /* Forward FFT and normalize */
  realft(U, NH, +1);
  double inv_an = 1.0 / (double)NFFT;
  for (int i = 0; i < NFFT; i++)
    U[i] = U[i] * inv_an;

  /* Compute RHS for each mode j=1..NH-1 */
  for (int j = 1; j < NH; j++) {
    double qj = QQ * j;
    double diag = qj * qj * (1.0 - qj * qj);
    f[2 * j - 2] = diag * yv[2 * j - 2] + QQ * j * U[2 * j + 1];
    f[2 * j - 1] = diag * yv[2 * j - 1] - QQ * j * U[2 * j];
  }

  return 0;
}

/* ---------------------------------------------------------------------------
 * Analytic Jacobian: diagonal only (linear stiffness)
 * Band matrix with mu=0, ml=0 (diagonal band)
 * ---------------------------------------------------------------------------*/
static int jac_ks(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                  void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)y; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;

  for (int j = 1; j < NH; j++) {
    double qj = QQ * j;
    double diag = qj * qj * (1.0 - qj * qj);
    SM_ELEMENT_B(J, 2 * j - 2, 2 * j - 2) = diag;
    SM_ELEMENT_B(J, 2 * j - 1, 2 * j - 1) = diag;
  }

  return 0;
}

/* ---------------------------------------------------------------------------
 * Reference solution at t=100 (from res_exact_pic.txt)
 * 47 values: y[0..9], y[10..99 step 5], y[100..1021 step 50]
 * ---------------------------------------------------------------------------*/
static const sunrealtype yref[] = {
   0.0387294818017374,
   0.0011594920936697,
   0.0161135569739015,
  -0.0497378990194494,
   0.0161388094706610,
  -0.0790294186417376,
  -0.0562945063709753,
  -0.0194501629216551,
  -0.0317390967305054,
  -0.0256301653146683,
  -0.0896908152487471,
   0.0452570608746347,
   0.0110762324013310,
   0.0095676261122614,
   0.1089375008451597,
  -0.1201229482736044,
  -0.1354295218545231,
  -0.0146538030311935,
  -0.0092223615970260,
  -0.2760611800802050,
  -0.0174686339809433,
   0.1982033914963349,
  -0.1768843940011891,
   0.3090308065592430,
  -0.0952958329765616,
  -0.1063846133920493,
   0.0390151203721951,
   0.0451950368314852,
   0.0477521089259518,
   0.0112038077853544,
   0.0003656188545739,
  -0.0001443415867781,
  -0.0000090034942150,
   0.0000016423881092,
   0.0000001490634819,
   0.0000000334389501,
  -0.0000000020970060,
  -0.0000000017295829,
   0.0000000000353463,
   0.0000000000036143,
   0.0000000000004182,
  -0.0000000000000002,
  -0.0000000000000070,
  -0.0000000000000012,
   0.0000000000000001,
   0.0000000000000000,
   0.0000000000000000
};
#define NREF 47

/* Reference indices matching Fortran output pattern:
   y[0..9], y[10..99 step 5], y[100..1021 step 50] */
static void get_ref_indices(int* indices)
{
  int cnt = 0;
  for (int i = 0; i < 10; i++)
    indices[cnt++] = i;
  for (int i = 10; i < 100; i += 5)
    indices[cnt++] = i;
  for (int i = 100; i < NEQ; i += 50)
    indices[cnt++] = i;
}

/* ---------------------------------------------------------------------------
 * main
 * ---------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  sunrealtype rtol = 1.0e-14;
  sunrealtype atol_val = 1.0e-14;
  sunrealtype h0   = 1.0e-6;
  int use_schur    = 1;
  int nsmin        = 13;
  int nsmax        = 13;
  if (argc > 1) rtol      = atof(argv[1]);
  if (argc > 2) atol_val  = atof(argv[2]);
  if (argc > 3) h0        = atof(argv[3]);
  if (argc > 4) use_schur = atoi(argv[4]);
  if (argc > 5) nsmin     = atoi(argv[5]);
  if (argc > 6) nsmax     = atoi(argv[6]);

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  /* Compute initial condition in physical space, then FFT to Fourier modes */
  double U[NFFT];
  double an = (double)NFFT;
  double delx = 1.0 / an;
  for (int i = 0; i < NFFT; i++) {
    double x = delx * i;
    double u1 = (x < 0.05) ? x : 0.1 - x;
    if (x > 0.1) u1 = 0.0;
    double u2 = (x >= 0.2 && x <= 0.3) ? 20.0 * (x - 0.2) * (0.3 - x) : 0.0;
    double u3 = 0.0;
    if (x >= 0.6 && x <= 0.7) u3 = (x - 0.6 < 0.7 - x) ? x - 0.6 : 0.7 - x;
    double u4 = 0.0;
    if (x >= 0.9 && x <= 1.0) u4 = (x - 0.9 < 1.0 - x) ? x - 0.9 : 1.0 - x;
    double val = u1;
    if (u2 > val) val = u2;
    if (u3 > val) val = u3;
    if (u4 > val) val = u4;
    if (val < 0.0) val = 0.0;
    U[i] = 16.0 * val;
  }

  /* Forward FFT and normalize */
  realft(U, NH, +1);
  for (int i = 0; i < NFFT; i++)
    U[i] = U[i] / an;

  /* Extract initial conditions */
  UZERO = U[0];

  N_Vector y0 = N_VNew_Serial(NEQ, sunctx);
  sunrealtype* y0v = N_VGetArrayPointer(y0);
  for (int i = 0; i < NEQ; i++)
    y0v[i] = U[i + 2];

  /* Solver setup */
  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, nsmin, nsmax);
  Radau5Init(mem, rhs_ks, 0.0, y0);

  /* Diagonal Jacobian: band matrix with mu=0, ml=0 */
  SUNMatrix Jt = SUNBandMatrix(NEQ, 0, 0, sunctx);
  Radau5SetLinearSolver(mem, Jt, NULL);
  Radau5SetSchurDecomp(mem, use_schur);
  Radau5SetJacFn(mem, jac_ks);
  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);

  /* Solve to t=100 */
  N_Vector yout = N_VNew_Serial(NEQ, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 100.0, yout, &tret);

  /* Compare with reference at selected indices */
  sunrealtype* yd = N_VGetArrayPointer(yout);
  int ref_indices[NREF];
  get_ref_indices(ref_indices);

  sunrealtype maxerr = 0.0;
  for (int r = 0; r < NREF; r++) {
    int idx = ref_indices[r];
    sunrealtype err = fabs(yd[idx] - yref[r]);
    if (err > maxerr) maxerr = err;
  }

  printf("=== KS problem (n=%d, rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n",
         NEQ, rtol, atol_val, h0, use_schur);
  printf("ret=%d, tret=%.10e\n", ret, tret);
  printf("max |y - yref| at ref points = %.6e\n", maxerr);

  /* Print first 10 components */
  printf("\nFirst 10 modes:\n");
  for (int i = 0; i < 10; i++)
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

  /* Pass criterion: relative error < 1e-3 for significant modes */
  int pass = (ret == 0 && maxerr < 1.0e-3);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout); SUNMatDestroy(Jt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
