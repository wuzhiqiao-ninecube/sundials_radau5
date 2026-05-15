/* ---------------------------------------------------------------------------
 * test_radau5_build_e1.c — Verify E1 matrix assembly
 *
 * For a known 3×3 Jacobian, builds E1 = fac1*I - J and checks each element.
 * Also tests the mass matrix case E1 = fac1*M - J.
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
#define TOL 1.0e-14

#define CHECK_VAL(name, got, exp) do { \
  if (fabs((got)-(exp)) > TOL) { \
    fprintf(stderr, "FAIL: %s = %.17e, expected %.17e\n", name, (double)(got), (double)(exp)); \
    nfail++; \
  } \
} while(0)

static int dummy_rhs(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)t; (void)y; (void)ud;
  N_VConst(0.0, yd);
  return 0;
}

int main(void)
{
  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  sunindextype n = 3;
  sunrealtype fac1 = 5.0;

  /* Create solver and initialize minimally */
  void* mem = Radau5Create(sunctx);
  N_Vector y0 = N_VNew_Serial(n, sunctx);
  N_VConst(1.0, y0);

  Radau5Init(mem, dummy_rhs, 0.0, y0);

  SUNMatrix Jtemplate = SUNDenseMatrix(n, n, sunctx);
  Radau5SetLinearSolver(mem, Jtemplate, NULL);

  Radau5Mem rmem = RADAU5_MEM(mem);

  /* Fill J with known values:
   * J = [ 1  2  3 ]
   *     [ 4  5  6 ]
   *     [ 7  8  9 ] */
  for (sunindextype j = 0; j < n; j++)
    for (sunindextype i = 0; i < n; i++)
      SM_ELEMENT_D(rmem->J, i, j) = (sunrealtype)(i * n + j + 1);

  /* --- Test 1: Identity mass (M == NULL) ---
   * E1 = fac1*I - J
   * E1[i][j] = -J[i][j] + fac1*(i==j) */
  radau5_BuildE1(rmem, fac1);

  for (sunindextype j = 0; j < n; j++) {
    for (sunindextype i = 0; i < n; i++) {
      sunrealtype expected = -(sunrealtype)(i * n + j + 1);
      if (i == j) expected += fac1;
      char buf[64];
      sprintf(buf, "E1_ident[%d][%d]", (int)i, (int)j);
      CHECK_VAL(buf, SM_ELEMENT_D(rmem->E1, i, j), expected);
    }
  }

  /* --- Test 2: General mass matrix ---
   * M = [ 2  0  0 ]
   *     [ 0  3  0 ]
   *     [ 0  0  1 ]
   * E1 = fac1*M - J */
  rmem->M = SUNDenseMatrix(n, n, sunctx);
  SUNMatZero(rmem->M);
  SM_ELEMENT_D(rmem->M, 0, 0) = 2.0;
  SM_ELEMENT_D(rmem->M, 1, 1) = 3.0;
  SM_ELEMENT_D(rmem->M, 2, 2) = 1.0;

  radau5_BuildE1(rmem, fac1);

  /* E1[i][j] = fac1*M[i][j] - J[i][j] */
  sunrealtype Mvals[3][3] = {{2,0,0},{0,3,0},{0,0,1}};
  for (sunindextype j = 0; j < n; j++) {
    for (sunindextype i = 0; i < n; i++) {
      sunrealtype expected = fac1 * Mvals[i][j] - (sunrealtype)(i * n + j + 1);
      char buf[64];
      sprintf(buf, "E1_mass[%d][%d]", (int)i, (int)j);
      CHECK_VAL(buf, SM_ELEMENT_D(rmem->E1, i, j), expected);
    }
  }

  /* Cleanup — set M back to NULL before Free (it's not owned by the solver in the normal path) */
  SUNMatDestroy(rmem->M);
  rmem->M = NULL;

  N_VDestroy(y0);
  SUNMatDestroy(Jtemplate);
  Radau5Free(&mem);
  SUNContext_Free(&sunctx);

  if (nfail == 0)
    printf("test_radau5_build_e1: PASSED\n");
  else
    printf("test_radau5_build_e1: %d FAILURES\n", nfail);

  return nfail;
}
