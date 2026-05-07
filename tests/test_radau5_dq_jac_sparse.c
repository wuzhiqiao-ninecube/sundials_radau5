/* ---------------------------------------------------------------------------
 * test_radau5_dq_jac_sparse.c — Verify sparse DQ Jacobian via column grouping
 *
 * Uses f(y) = [y0^2, y0*y1, y1^3] with a known sparse pattern and compares
 * radau5_DQJacSparse output against the analytic Jacobian.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include "radau5.h"
#include "radau5_impl.h"

static int nfail = 0;
#define DQ_TOL 1.0e-5

#define CHECK_VAL(name, got, exp) do { \
  sunrealtype _err = fabs((got)-(exp)); \
  sunrealtype _scale = fabs(exp) > 1.0 ? fabs(exp) : 1.0; \
  if (_err / _scale > DQ_TOL) { \
    fprintf(stderr, "FAIL: %s = %.10e, expected %.10e, rel_err = %.3e\n", \
            name, (double)(got), (double)(exp), (double)(_err/_scale)); \
    nfail++; \
  } \
} while(0)

/* f(y) = [y0^2, y0*y1, y1^3]
 * Analytic Jacobian (CSC, only nonzeros):
 *   col 0: row 0 = 2*y0, row 1 = y1
 *   col 1: row 1 = y0,   row 2 = 3*y1^2
 *   col 2: (empty)
 */
static int rhs_poly(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)t; (void)ud;
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

  sunindextype n   = 3;
  sunindextype nnz = 4; /* nonzeros: (0,0),(1,0),(1,1),(2,1) */

  /* Build sparse J template: CSC 3x3 */
  SUNMatrix Jt = SUNSparseMatrix(n, n, nnz, CSC_MAT, sunctx);

  sunindextype* cp = SM_INDEXPTRS_S(Jt);
  sunindextype* ri = SM_INDEXVALS_S(Jt);
  sunrealtype*  dv = SM_DATA_S(Jt);

  /* col 0: rows 0,1 */
  cp[0] = 0; ri[0] = 0; ri[1] = 1;
  /* col 1: rows 1,2 */
  cp[1] = 2; ri[2] = 1; ri[3] = 2;
  /* col 2: empty */
  cp[2] = 4;
  cp[3] = 4;
  dv[0] = dv[1] = dv[2] = dv[3] = SUN_RCONST(0.0);

  /* Set up solver */
  N_Vector y0 = N_VNew_Serial(n, sunctx);
  sunrealtype* y0d = N_VGetArrayPointer(y0);
  y0d[0] = 2.0; y0d[1] = 3.0; y0d[2] = 0.5;

  void* mem = Radau5Create(sunctx);
  Radau5Init(mem, rhs_poly, 0.0, y0);
  Radau5SStolerances(mem, 1.0e-8, 1.0e-10);
  Radau5SetInitStep(mem, 0.01);
  Radau5SetLinearSolver(mem, Jt);

  /* Register sparsity pattern — triggers column grouping */
  int ret = Radau5SetSparsityPattern(mem, Jt);
  if (ret != RADAU5_SUCCESS) {
    fprintf(stderr, "FAIL: Radau5SetSparsityPattern returned %d\n", ret);
    nfail++;
  }

  Radau5Mem rmem = RADAU5_MEM(mem);

  /* Verify grouping: cols 0 and 1 share row 1, so they must be in different
   * groups. Col 2 is empty (group -1). ngroups must be 2. */
  if (rmem->ngroups != 2) {
    fprintf(stderr, "FAIL: ngroups = %ld, expected 2\n", (long)rmem->ngroups);
    nfail++;
  }
  if (rmem->col_group[0] == rmem->col_group[1]) {
    fprintf(stderr, "FAIL: cols 0 and 1 assigned to same group %ld\n",
            (long)rmem->col_group[0]);
    nfail++;
  }
  if (rmem->col_group[2] != -1) {
    fprintf(stderr, "FAIL: empty col 2 should have group -1, got %ld\n",
            (long)rmem->col_group[2]);
    nfail++;
  }

  /* Set h and compute scal */
  rmem->h = 0.01;
  radau5_ComputeScal(rmem, rmem->ycur);

  /* Evaluate f(t0, y0) -> fn */
  rhs_poly(0.0, rmem->ycur, rmem->fn, NULL);

  /* Compute sparse DQ Jacobian */
  ret = radau5_DQJacSparse(rmem, 0.0, rmem->ycur, rmem->fn);
  if (ret != RADAU5_SUCCESS) {
    fprintf(stderr, "FAIL: radau5_DQJacSparse returned %d\n", ret);
    nfail++;
  }

  /* Analytic Jacobian at y=[2,3,0.5]:
   *   J[0,0] = 2*y0 = 4
   *   J[1,0] = y1   = 3
   *   J[1,1] = y0   = 2
   *   J[2,1] = 3*y1^2 = 27
   */
  sunrealtype* Jd = SM_DATA_S(rmem->J);
  /* CSC order: Jd[0]=(0,0), Jd[1]=(1,0), Jd[2]=(1,1), Jd[3]=(2,1) */
  CHECK_VAL("J[0,0]", Jd[0],  4.0);
  CHECK_VAL("J[1,0]", Jd[1],  3.0);
  CHECK_VAL("J[1,1]", Jd[2],  2.0);
  CHECK_VAL("J[2,1]", Jd[3], 27.0);

  /* Cleanup */
  N_VDestroy(y0);
  SUNMatDestroy(Jt);
  Radau5Free(&mem);
  SUNContext_Free(&sunctx);

  if (nfail == 0)
    printf("test_radau5_dq_jac_sparse: PASSED\n");
  else
    printf("test_radau5_dq_jac_sparse: %d FAILURES\n", nfail);

  return nfail;
}
