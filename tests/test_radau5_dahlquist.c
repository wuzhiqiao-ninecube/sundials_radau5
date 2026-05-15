/* ---------------------------------------------------------------------------
 * test_radau5_dahlquist.c — Order verification on y' = lambda*y
 *
 * Solves y' = -y from t=0 to t=1 (exact: y = e^{-t}) at multiple
 * tolerances and verifies the global error scales as O(tol).
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

static int nfail = 0;

#define CHECK(cond, msg) do { \
  if (!(cond)) { fprintf(stderr, "FAIL: %s\n", msg); nfail++; } \
} while(0)

/* RHS: y' = -y */
static int rhs_exp(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  sunrealtype* ydata  = N_VGetArrayPointer(y);
  sunrealtype* yddata = N_VGetArrayPointer(yd);
  yddata[0] = -ydata[0];
  return 0;
}

/* Analytic Jacobian: J = -1 */
static int jac_exp(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                   void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  SM_ELEMENT_D(J, 0, 0) = -1.0;
  return 0;
}

static sunrealtype run_one(SUNContext sunctx, sunrealtype rtol, sunrealtype atol,
                           int use_analytic_jac)
{
  void* mem = Radau5Create(sunctx);
  N_Vector y0 = N_VNew_Serial(1, sunctx);
  N_VConst(1.0, y0);

  Radau5SetOrderLimits(mem, 3, 3);  /* fixed ns=3 for convergence order test */
  Radau5Init(mem, rhs_exp, 0.0, y0);

  SUNMatrix J = SUNDenseMatrix(1, 1, sunctx);
  Radau5SetLinearSolver(mem, J, NULL);

  if (use_analytic_jac)
    Radau5SetJacFn(mem, jac_exp);

  Radau5SStolerances(mem, rtol, atol);
  Radau5SetInitStep(mem, rtol); /* scale initial step with tolerance */

  N_Vector yout = N_VNew_Serial(1, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 1.0, yout, &tret);

  sunrealtype err = 0.0;
  if (ret == RADAU5_SUCCESS) {
    sunrealtype exact = exp(-1.0);
    sunrealtype* yd = N_VGetArrayPointer(yout);
    err = fabs(yd[0] - exact);
  } else {
    fprintf(stderr, "  Solve failed with ret=%d at rtol=%.1e\n", ret, rtol);
    err = 1.0; /* signal failure */
  }

  N_VDestroy(y0);
  N_VDestroy(yout);
  SUNMatDestroy(J);
  Radau5Free(&mem);

  return err;
}

int main(void)
{
  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  printf("--- Dahlquist y'=-y order test (DQ Jacobian) ---\n");

  sunrealtype tols[]  = { 1.0e-4, 1.0e-6, 1.0e-8, 1.0e-10 };
  sunrealtype errs[4];
  int ntols = 4;

  for (int i = 0; i < ntols; i++) {
    errs[i] = run_one(sunctx, tols[i], tols[i], 0);
    printf("  rtol=%.1e  err=%.3e\n", tols[i], errs[i]);
  }

  /* Check that error decreases with tolerance */
  for (int i = 1; i < ntols; i++) {
    if (errs[i-1] > 1.0e-15 && errs[i] > 1.0e-15) {
      sunrealtype ratio = log(errs[i-1] / errs[i]) / log(tols[i-1] / tols[i]);
      printf("  order estimate [%d->%d]: %.2f\n", i-1, i, ratio);
      /* For a well-behaved order-5 method, global error ~ C*tol,
         so ratio should be approximately 1.0 (within some margin) */
      CHECK(ratio > 0.5, "order ratio should be > 0.5");
      CHECK(ratio < 2.0, "order ratio should be < 2.0");
    }
  }

  /* Verify error is actually small at tight tolerance */
  CHECK(errs[ntols-1] < 1.0e-8, "error at tol=1e-10 should be < 1e-8");

  /* --- Repeat with analytic Jacobian --- */
  printf("--- Dahlquist y'=-y order test (analytic Jacobian) ---\n");

  sunrealtype errs_aj[4];
  for (int i = 0; i < ntols; i++) {
    errs_aj[i] = run_one(sunctx, tols[i], tols[i], 1);
    printf("  rtol=%.1e  err=%.3e\n", tols[i], errs_aj[i]);
  }

  CHECK(errs_aj[ntols-1] < 1.0e-8,
        "error at tol=1e-10 (analytic J) should be < 1e-8");

  SUNContext_Free(&sunctx);

  if (nfail == 0)
    printf("test_radau5_dahlquist: PASSED\n");
  else
    printf("test_radau5_dahlquist: %d FAILURES\n", nfail);

  return nfail;
}
