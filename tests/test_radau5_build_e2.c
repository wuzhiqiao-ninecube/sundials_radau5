/* ---------------------------------------------------------------------------
 * test_radau5_build_e2.c — Verify E2 realified 2n×2n matrix assembly
 *
 * For a known 2×2 Jacobian, builds the realified complex matrix E2 and
 * checks the 2×2 block structure.
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

  sunindextype n = 2;
  sunrealtype alphn = 3.0;
  sunrealtype betan = 2.0;

  /* Create solver and initialize minimally */
  void* mem = Radau5Create(sunctx);
  N_Vector y0 = N_VNew_Serial(n, sunctx);
  N_VConst(1.0, y0);

  Radau5Init(mem, dummy_rhs, 0.0, y0);

  SUNMatrix Jtemplate = SUNDenseMatrix(n, n, sunctx);
  Radau5SetLinearSolver(mem, Jtemplate);

  Radau5Mem rmem = RADAU5_MEM(mem);

  /* Fill J = [ 1  2 ]
   *          [ 3  4 ] */
  SM_ELEMENT_D(rmem->J, 0, 0) = 1.0;
  SM_ELEMENT_D(rmem->J, 0, 1) = 2.0;
  SM_ELEMENT_D(rmem->J, 1, 0) = 3.0;
  SM_ELEMENT_D(rmem->J, 1, 1) = 4.0;

  /* --- Test 1: Identity mass (M == NULL) ---
   * E2 = [ alphn*I - J,   -betan*I ]
   *      [ betan*I,     alphn*I - J ]
   *
   * Top-left:     alphn*delta(i,j) - J(i,j)
   * Top-right:   -betan*delta(i,j)
   * Bottom-left:  betan*delta(i,j)
   * Bottom-right: alphn*delta(i,j) - J(i,j)
   */
  radau5_BuildE2(rmem, alphn, betan);

  sunrealtype Jv[2][2] = {{1,2},{3,4}};

  for (sunindextype j = 0; j < n; j++) {
    for (sunindextype i = 0; i < n; i++) {
      sunrealtype delta_ij = (i == j) ? 1.0 : 0.0;
      sunrealtype e_tl = alphn * delta_ij - Jv[i][j];
      sunrealtype e_tr = -betan * delta_ij;
      sunrealtype e_bl = betan * delta_ij;
      sunrealtype e_br = alphn * delta_ij - Jv[i][j];

      char buf[64];
      sprintf(buf, "E2_TL[%d][%d]", (int)i, (int)j);
      CHECK_VAL(buf, SM_ELEMENT_D(rmem->E2, i, j), e_tl);

      sprintf(buf, "E2_TR[%d][%d+n]", (int)i, (int)j);
      CHECK_VAL(buf, SM_ELEMENT_D(rmem->E2, i, j + n), e_tr);

      sprintf(buf, "E2_BL[%d+n][%d]", (int)i, (int)j);
      CHECK_VAL(buf, SM_ELEMENT_D(rmem->E2, i + n, j), e_bl);

      sprintf(buf, "E2_BR[%d+n][%d+n]", (int)i, (int)j);
      CHECK_VAL(buf, SM_ELEMENT_D(rmem->E2, i + n, j + n), e_br);
    }
  }

  /* --- Test 2: General mass matrix ---
   * M = [ 2  0 ]
   *     [ 0  3 ] */
  rmem->M = SUNDenseMatrix(n, n, sunctx);
  SUNMatZero(rmem->M);
  SM_ELEMENT_D(rmem->M, 0, 0) = 2.0;
  SM_ELEMENT_D(rmem->M, 1, 1) = 3.0;

  radau5_BuildE2(rmem, alphn, betan);

  sunrealtype Mv[2][2] = {{2,0},{0,3}};

  for (sunindextype j = 0; j < n; j++) {
    for (sunindextype i = 0; i < n; i++) {
      sunrealtype e_tl = alphn * Mv[i][j] - Jv[i][j];
      sunrealtype e_tr = -betan * Mv[i][j];
      sunrealtype e_bl = betan * Mv[i][j];
      sunrealtype e_br = alphn * Mv[i][j] - Jv[i][j];

      char buf[64];
      sprintf(buf, "E2_mass_TL[%d][%d]", (int)i, (int)j);
      CHECK_VAL(buf, SM_ELEMENT_D(rmem->E2, i, j), e_tl);

      sprintf(buf, "E2_mass_TR[%d][%d+n]", (int)i, (int)j);
      CHECK_VAL(buf, SM_ELEMENT_D(rmem->E2, i, j + n), e_tr);

      sprintf(buf, "E2_mass_BL[%d+n][%d]", (int)i, (int)j);
      CHECK_VAL(buf, SM_ELEMENT_D(rmem->E2, i + n, j), e_bl);

      sprintf(buf, "E2_mass_BR[%d+n][%d+n]", (int)i, (int)j);
      CHECK_VAL(buf, SM_ELEMENT_D(rmem->E2, i + n, j + n), e_br);
    }
  }

  /* Cleanup */
  SUNMatDestroy(rmem->M);
  rmem->M = NULL;

  N_VDestroy(y0);
  SUNMatDestroy(Jtemplate);
  Radau5Free(&mem);
  SUNContext_Free(&sunctx);

  if (nfail == 0)
    printf("test_radau5_build_e2: PASSED\n");
  else
    printf("test_radau5_build_e2: %d FAILURES\n", nfail);

  return nfail;
}
