/* ---------------------------------------------------------------------------
 * radau5_fekete.c — Fekete problem (DAE index-2, n=160)
 *
 * 20 charged particles on the unit sphere seeking equilibrium under
 * Coulomb repulsion with holonomic constraints |p_i|^2 = 1.
 *
 * State vector (0-based, nart=20):
 *   y[0..59]    positions  p(i,k), i=0..19, k=0..2
 *   y[60..119]  velocities q(i,k)
 *   y[120..139] Lagrange multipliers lambda(i)
 *   y[140..159] Lagrange multipliers mu(i)
 *
 * M*y' = f(t,y), M = diag(1..1, 0..0) (first 120 differential)
 * DAE index: nind1=120, nind2=40, nind3=0
 *
 * t in [0, 1000]
 * Reference from IVPtestset (RADAU5, Cray C90, tol=1e-12)
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

#define NEQ  160
#define NART 20

#ifndef M_PI
#define M_PI 3.141592653589793238462643383
#endif

/* ---------------------------------------------------------------------------
 * RHS: f(t, y) — Coulomb forces on unit sphere with stabilization
 *
 * Translated from feval in fekete.f:
 *   pp(i,k) = q(i,k) + 2*mu(i)*p(i,k)
 *   qp(i,k) = -alpha*q(i,k) + 2*lam(i)*p(i,k) + sum_j f(i,j,k)
 *   phi(i)  = sum_k p(i,k)^2 - 1
 *   gpq(i)  = sum_k 2*p(i,k)*q(i,k)
 * where f(i,j,k) = (p(i,k)-p(j,k)) / |p_i - p_j|^2  (2D Coulomb on sphere)
 * ---------------------------------------------------------------------------*/
static int rhs_fekete(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)t; (void)ud;
  sunrealtype* yv = N_VGetArrayPointer(y);
  sunrealtype* dy = N_VGetArrayPointer(yd);

  const sunrealtype alpha = 0.5;

  /* Extract p, q, lam, mu from state vector (0-based) */
  /* p(i,k) = yv[3*i + k],  i=0..19, k=0..2 */
  /* q(i,k) = yv[3*NART + 3*i + k] */
  /* lam(i) = yv[6*NART + i] */
  /* mu(i)  = yv[7*NART + i] */

  /* Compute pairwise forces f(i,j,k) = (p_i - p_j) / |p_i - p_j|^2 */
  sunrealtype fij[NART][NART][3];
  int i, j, k;
  for (i = 0; i < NART; i++) {
    for (j = 0; j < NART; j++) {
      if (i == j) {
        fij[i][j][0] = 0.0;
        fij[i][j][1] = 0.0;
        fij[i][j][2] = 0.0;
      } else {
        sunrealtype rn = 0.0;
        for (k = 0; k < 3; k++) {
          sunrealtype d = yv[3*i + k] - yv[3*j + k];
          rn += d * d;
        }
        for (k = 0; k < 3; k++) {
          fij[i][j][k] = (yv[3*i + k] - yv[3*j + k]) / rn;
        }
      }
    }
  }

  /* Compute derivatives */
  for (i = 0; i < NART; i++) {
    sunrealtype lam_i = yv[6*NART + i];
    sunrealtype mu_i  = yv[7*NART + i];
    for (k = 0; k < 3; k++) {
      sunrealtype p_ik = yv[3*i + k];
      sunrealtype q_ik = yv[3*NART + 3*i + k];
      /* pp(i,k) = q(i,k) + 2*mu(i)*p(i,k) */
      dy[3*i + k] = q_ik + 2.0 * mu_i * p_ik;

      /* qp(i,k) = -alpha*q(i,k) + 2*lam(i)*p(i,k) + sum_j f(i,j,k) */
      sunrealtype qp_ik = -alpha * q_ik + 2.0 * lam_i * p_ik;
      for (j = 0; j < NART; j++) {
        qp_ik += fij[i][j][k];
      }
      dy[3*NART + 3*i + k] = qp_ik;
    }
  }

  /* Constraint equations */
  for (i = 0; i < NART; i++) {
    sunrealtype phi_i = -1.0;
    sunrealtype gpq_i = 0.0;
    for (k = 0; k < 3; k++) {
      sunrealtype p_ik = yv[3*i + k];
      sunrealtype q_ik = yv[3*NART + 3*i + k];
      phi_i += p_ik * p_ik;
      gpq_i += 2.0 * p_ik * q_ik;
    }
    dy[6*NART + i] = phi_i;
    dy[7*NART + i] = gpq_i;
  }

  return 0;
}

/* ---------------------------------------------------------------------------
 * Mass matrix: M = diag(1..1, 0..0), first 120 differential
 * ---------------------------------------------------------------------------*/
static int mas_fekete(sunrealtype t, SUNMatrix M, void* ud,
                      N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)ud; (void)t1; (void)t2; (void)t3;
  SUNMatZero(M);
  for (int i = 0; i < 6*NART; i++)
    SM_ELEMENT_D(M, i, i) = 1.0;
  return 0;
}

/* ---------------------------------------------------------------------------
 * Raw feval on plain arrays (used by fekete_init for consistent IC)
 * Same logic as rhs_fekete but operates on sunrealtype* directly.
 * ---------------------------------------------------------------------------*/
static void fekete_feval_raw(const sunrealtype* yv, sunrealtype* dy)
{
  const sunrealtype alpha = 0.5;
  int i, j, k;

  /* Pairwise forces */
  sunrealtype fij[NART][NART][3];
  for (i = 0; i < NART; i++) {
    for (j = 0; j < NART; j++) {
      if (i == j) {
        fij[i][j][0] = fij[i][j][1] = fij[i][j][2] = 0.0;
      } else {
        sunrealtype rn = 0.0;
        for (k = 0; k < 3; k++) {
          sunrealtype d = yv[3*i + k] - yv[3*j + k];
          rn += d * d;
        }
        for (k = 0; k < 3; k++)
          fij[i][j][k] = (yv[3*i + k] - yv[3*j + k]) / rn;
      }
    }
  }

  for (i = 0; i < NART; i++) {
    sunrealtype lam_i = yv[6*NART + i];
    sunrealtype mu_i  = yv[7*NART + i];
    for (k = 0; k < 3; k++) {
      sunrealtype p_ik = yv[3*i + k];
      sunrealtype q_ik = yv[3*NART + 3*i + k];
      dy[3*i + k] = q_ik + 2.0 * mu_i * p_ik;
      sunrealtype qp_ik = -alpha * q_ik + 2.0 * lam_i * p_ik;
      for (j = 0; j < NART; j++) qp_ik += fij[i][j][k];
      dy[3*NART + 3*i + k] = qp_ik;
    }
  }

  for (i = 0; i < NART; i++) {
    sunrealtype phi_i = -1.0, gpq_i = 0.0;
    for (k = 0; k < 3; k++) {
      sunrealtype p_ik = yv[3*i + k];
      sunrealtype q_ik = yv[3*NART + 3*i + k];
      phi_i += p_ik * p_ik;
      gpq_i += 2.0 * p_ik * q_ik;
    }
    dy[6*NART + i] = phi_i;
    dy[7*NART + i] = gpq_i;
  }
}

/* ---------------------------------------------------------------------------
 * Initial conditions (from init in fekete.f)
 * ---------------------------------------------------------------------------*/
static void fekete_init(sunrealtype* yv)
{
  int i, j, k;
  sunrealtype a, b;

  /* Band 1: particles 0..2 (Fortran i=1..3) */
  for (i = 0; i < 3; i++) {
    a = 2.0 * M_PI * (double)(i + 1) / 3.0 + M_PI / 13.0;
    b = 3.0 * M_PI / 8.0;
    yv[3*i + 0] = cos(a) * cos(b);
    yv[3*i + 1] = sin(a) * cos(b);
    yv[3*i + 2] = sin(b);
  }
  /* Band 2: particles 3..9 (Fortran i=4..10) */
  for (i = 3; i < 10; i++) {
    a = 2.0 * M_PI * (double)(i - 2) / 7.0 + M_PI / 29.0;
    b = M_PI / 8.0;
    yv[3*i + 0] = cos(a) * cos(b);
    yv[3*i + 1] = sin(a) * cos(b);
    yv[3*i + 2] = sin(b);
  }
  /* Band 3: particles 10..15 (Fortran i=11..16) */
  for (i = 10; i < 16; i++) {
    a = 2.0 * M_PI * (double)(i - 9) / 6.0 + M_PI / 7.0;
    b = -2.0 * M_PI / 15.0;
    yv[3*i + 0] = cos(a) * cos(b);
    yv[3*i + 1] = sin(a) * cos(b);
    yv[3*i + 2] = sin(b);
  }
  /* Band 4: particles 16..19 (Fortran i=17..20) */
  for (i = 16; i < 20; i++) {
    a = 2.0 * M_PI * (double)(i - 16) / 4.0 + M_PI / 17.0;
    b = -3.0 * M_PI / 10.0;
    yv[3*i + 0] = cos(a) * cos(b);
    yv[3*i + 1] = sin(a) * cos(b);
    yv[3*i + 2] = sin(b);
  }

  /* Velocities = 0 */
  for (i = 3*NART; i < 6*NART; i++) yv[i] = 0.0;
  /* Multipliers = 0 */
  for (i = 6*NART; i < 8*NART; i++) yv[i] = 0.0;

  /* Compute consistent initial lambda values.
   * Replicate feval logic on raw arrays (avoids N_Vector dependency). */
  sunrealtype dy[NEQ];
  fekete_feval_raw(yv, dy);

  /* lam(i) = -sum_k p(i,k)*qp(i,k) / 2 */
  for (i = 0; i < NART; i++) {
    sunrealtype s = 0.0;
    for (k = 0; k < 3; k++) {
      s += yv[3*i + k] * dy[3*NART + 3*i + k];
    }
    yv[6*NART + i] = -s / 2.0;
  }

  /* Call feval again with updated lambda (makes yprime consistent) */
  fekete_feval_raw(yv, dy);
}

int main(int argc, char* argv[])
{
  sunrealtype rtol = 1.0e-6;
  sunrealtype atol_val = 1.0e-6;
  sunrealtype h0   = 1.0e-6;
  int use_schur    = 1;
  if (argc > 1) rtol      = atof(argv[1]);
  if (argc > 2) atol_val  = atof(argv[2]);
  if (argc > 3) h0        = atof(argv[3]);
  if (argc > 4) use_schur = atoi(argv[4]);

  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  void* mem = Radau5Create(sunctx);

  /* Initial conditions */
  N_Vector y0 = N_VNew_Serial(NEQ, sunctx);
  sunrealtype* y0v = N_VGetArrayPointer(y0);
  fekete_init(y0v);

  Radau5Init(mem, rhs_fekete, 0.0, y0);

  SUNMatrix Jt = SUNDenseMatrix(NEQ, NEQ, sunctx);
  Radau5SetLinearSolver(mem, Jt);
  Radau5SetSchurDecomp(mem, use_schur);
  /* No Radau5SetJacFn — use DQ Jacobian */

  SUNMatrix Mt = SUNDenseMatrix(NEQ, NEQ, sunctx);
  Radau5SetMassFn(mem, mas_fekete, Mt);
  Radau5SetDAEIndex(mem, 120, 40, 0);

  Radau5SStolerances(mem, rtol, atol_val);
  Radau5SetInitStep(mem, h0);

  N_Vector yout = N_VNew_Serial(NEQ, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 1000.0, yout, &tret);

  /* Reference solution at t=1000 (from solut in fekete.f, first 6 components) */
  static const sunrealtype yref[6] = {
    -0.4070263380333202e+00,
     0.3463758772791802e+00,
     0.8451942450030429e+00,
     0.7752934752521549e-01,
    -0.2628662719972299e+00,
     0.9617122871829146e+00
  };

  sunrealtype* yd = N_VGetArrayPointer(yout);
  printf("=== Fekete problem (rtol=%.1e atol=%.1e h0=%.1e schur=%d) ===\n",
         rtol, atol_val, h0, use_schur);
  printf("ret=%d, tret=%.6e\n", ret, tret);

  sunrealtype maxrelerr = 0.0;
  for (int i = 0; i < 6; i++) {
    sunrealtype err = fabs(yd[i] - yref[i]);
    sunrealtype relerr = (fabs(yref[i]) > 1e-15)
                           ? err / fabs(yref[i]) : err;
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

  int pass = (ret == 0 && maxrelerr < 1.0e-2);
  printf("max_pos_rel_err=%.3e\n", maxrelerr);
  printf("%s\n\n", pass ? "PASSED" : "FAILED");

  N_VDestroy(y0); N_VDestroy(yout);
  SUNMatDestroy(Jt); SUNMatDestroy(Mt);
  Radau5Free(&mem); SUNContext_Free(&sunctx);
  return pass ? 0 : 1;
}
