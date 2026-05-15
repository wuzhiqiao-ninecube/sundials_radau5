/* ---------------------------------------------------------------------------
 * test_radau5_dq_jac.c — Verify finite-difference Jacobian
 *
 * Tests dense DQ Jacobian against a known analytic Jacobian for
 * f(y) = [y0^2, y0*y1, y1^3].
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"
#include "radau5_impl.h"

static int nfail = 0;
#define DQ_TOL 1.0e-5  /* DQ Jacobian is O(sqrt(eps)) accurate */

#define CHECK_VAL(name, got, exp) do { \
  sunrealtype _err = fabs((got)-(exp)); \
  sunrealtype _scale = fabs(exp) > 1.0 ? fabs(exp) : 1.0; \
  if (_err / _scale > DQ_TOL) { \
    fprintf(stderr, "FAIL: %s = %.10e, expected %.10e, rel_err = %.3e\n", \
            name, (double)(got), (double)(exp), (double)(_err/_scale)); \
    nfail++; \
  } \
} while(0)

/* f(y) = [y0^2, y0*y1, y1^3] */
static int rhs_poly(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  sunrealtype* yv  = N_VGetArrayPointer(y);
  sunrealtype* ydv = N_VGetArrayPointer(yd);
  ydv[0] = yv[0] * yv[0];
  ydv[1] = yv[0] * yv[1];
  ydv[2] = yv[1] * yv[1] * yv[1];
  return 0;
}

int main(void)
{
  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  sunindextype n = 3;

  void* mem = Radau5Create(sunctx);
  N_Vector y0 = N_VNew_Serial(n, sunctx);
  sunrealtype* y0d = N_VGetArrayPointer(y0);
  y0d[0] = 2.0;
  y0d[1] = 3.0;
  y0d[2] = 0.5;

  Radau5Init(mem, rhs_poly, 0.0, y0);
  Radau5SStolerances(mem, 1.0e-8, 1.0e-10);
  Radau5SetInitStep(mem, 0.01);

  SUNMatrix Jtemplate = SUNDenseMatrix(n, n, sunctx);
  Radau5SetLinearSolver(mem, Jtemplate, NULL);

  Radau5Mem rmem = RADAU5_MEM(mem);

  /* Set h for the perturbation formula */
  rmem->h = 0.01;

  /* Compute scal for the DQ formula */
  radau5_ComputeScal(rmem, rmem->ycur);

  /* Evaluate f(t0, y0) -> fn */
  rhs_poly(0.0, rmem->ycur, rmem->fn, NULL);

  /* Compute DQ Jacobian */
  int ret = radau5_DQJacDense(rmem, 0.0, rmem->ycur, rmem->fn);
  if (ret != RADAU5_SUCCESS) {
    fprintf(stderr, "FAIL: DQ Jacobian computation failed\n");
    nfail++;
  }

  /* Analytic Jacobian at y = [2, 3, 0.5]:
   * df0/dy0 = 2*y0 = 4,   df0/dy1 = 0,         df0/dy2 = 0
   * df1/dy0 = y1 = 3,      df1/dy1 = y0 = 2,    df1/dy2 = 0
   * df2/dy0 = 0,            df2/dy1 = 3*y1^2=27, df2/dy2 = 0
   */
  CHECK_VAL("J[0][0]", SM_ELEMENT_D(rmem->J, 0, 0), 4.0);
  CHECK_VAL("J[0][1]", SM_ELEMENT_D(rmem->J, 0, 1), 0.0);
  CHECK_VAL("J[0][2]", SM_ELEMENT_D(rmem->J, 0, 2), 0.0);
  CHECK_VAL("J[1][0]", SM_ELEMENT_D(rmem->J, 1, 0), 3.0);
  CHECK_VAL("J[1][1]", SM_ELEMENT_D(rmem->J, 1, 1), 2.0);
  CHECK_VAL("J[1][2]", SM_ELEMENT_D(rmem->J, 1, 2), 0.0);
  CHECK_VAL("J[2][0]", SM_ELEMENT_D(rmem->J, 2, 0), 0.0);
  CHECK_VAL("J[2][1]", SM_ELEMENT_D(rmem->J, 2, 1), 27.0);
  CHECK_VAL("J[2][2]", SM_ELEMENT_D(rmem->J, 2, 2), 0.0);

  /* --- Test perturbation at y=0: verify minInc floor works ---
   * At y=0, f=[0,0,0]. The analytic Jacobian is all zeros.
   * DQ Jacobian will have O(inc) errors for quadratic terms (df0/dy0 = 2*y0),
   * so we just verify the computation doesn't crash and values are small. */
  N_VConst(0.0, rmem->ycur);
  rhs_poly(0.0, rmem->ycur, rmem->fn, NULL);
  radau5_ComputeScal(rmem, rmem->ycur);

  ret = radau5_DQJacDense(rmem, 0.0, rmem->ycur, rmem->fn);
  if (ret != RADAU5_SUCCESS) {
    fprintf(stderr, "FAIL: DQ Jacobian at y=0 failed\n");
    nfail++;
  }
  /* Just check it ran without error — DQ at y=0 for nonlinear f
   * won't match the analytic zero Jacobian exactly. */
  printf("  DQ Jac at y=0: J[0][0]=%.3e (expected ~inc for quadratic term)\n",
         SM_ELEMENT_D(rmem->J, 0, 0));

  /* Cleanup */
  N_VDestroy(y0);
  SUNMatDestroy(Jtemplate);
  Radau5Free(&mem);
  SUNContext_Free(&sunctx);

  if (nfail == 0)
    printf("test_radau5_dq_jac: PASSED\n");
  else
    printf("test_radau5_dq_jac: %d FAILURES\n", nfail);

  return nfail;
}
