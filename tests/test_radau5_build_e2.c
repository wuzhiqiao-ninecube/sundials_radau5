/* ---------------------------------------------------------------------------
 * test_radau5_build_e2.c — Verify complex n×n E2 matrix assembly and solve
 *
 * For a known 2×2 Jacobian, builds the complex E2 matrix and verifies
 * that the complex solve produces the same result as the analytical solution.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"
#include "radau5_impl.h"

static int nfail = 0;
#define TOL 1.0e-12

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
  Radau5SetLinearSolver(mem, Jtemplate, NULL);

  Radau5Mem rmem = RADAU5_MEM(mem);

  /* Fill J = [ 1  2 ]
   *          [ 3  4 ] */
  SM_ELEMENT_D(rmem->J, 0, 0) = 1.0;
  SM_ELEMENT_D(rmem->J, 0, 1) = 2.0;
  SM_ELEMENT_D(rmem->J, 1, 0) = 3.0;
  SM_ELEMENT_D(rmem->J, 1, 1) = 4.0;

  /* --- Test 1: Identity mass (M == NULL), eigenvalue mode ---
   * Complex matrix K = (alphn + i*betan)*I - J
   * K[0,0] = (3+2i) - 1 = 2+2i
   * K[0,1] = 0 - 2 = -2
   * K[1,0] = 0 - 3 = -3
   * K[1,1] = (3+2i) - 4 = -1+2i
   */
  rmem->h = 1.0;  /* so alphn/h = alphn */
  radau5_BuildE2c(rmem, 0, alphn, betan);

  /* Check omega = 1 for eigenvalue mode (a01=-betan, a10=betan) */
  CHECK_VAL("omega_eigen", rmem->e2c_omega[0], 1.0);

  /* Check E2c_data: complex column-major via SM_ELEMENT_ZD
   * Col 0: (0,0)=(2,2), (1,0)=(-3,0)
   * Col 1: (0,1)=(-2,0), (1,1)=(-1,2) */
  suncomplextype *E2d = SM_DATA_ZD(rmem->E2c_mat[0]);
  CHECK_VAL("E2c[0,0].re", creal(E2d[0]), 2.0);
  CHECK_VAL("E2c[0,0].im", cimag(E2d[0]), 2.0);
  CHECK_VAL("E2c[1,0].re", creal(E2d[1]), -3.0);
  CHECK_VAL("E2c[1,0].im", cimag(E2d[1]), 0.0);
  CHECK_VAL("E2c[0,1].re", creal(E2d[2]), -2.0);
  CHECK_VAL("E2c[0,1].im", cimag(E2d[2]), 0.0);
  CHECK_VAL("E2c[1,1].re", creal(E2d[3]), -1.0);
  CHECK_VAL("E2c[1,1].im", cimag(E2d[3]), 2.0);

  /* --- Test 2: Factor and solve ---
   * Solve K*z = d where d = (1+0i, 0+0i)
   * i.e. rhs_re = [1, 0], rhs_im = [0, 0] */
  int ret = radau5_DecompE2c(rmem, 0);
  CHECK_VAL("decomp_ret", (double)ret, 0.0);

  sunrealtype rhs_re[2] = {1.0, 0.0};
  sunrealtype rhs_im[2] = {0.0, 0.0};
  radau5_SolveE2c(rmem, 0, rhs_re, rhs_im);

  /* Verify: K * x = rhs
   * x = K^{-1} * [1, 0]^T
   * K = [[2+2i, -2], [-3, -1+2i]]
   * det(K) = (2+2i)(-1+2i) - (-2)(-3) = (-2+4i-2i+4i^2) - 6
   *        = (-2+2i-4) - 6 = -12+2i
   * K^{-1} = 1/det * [[-1+2i, 2], [3, 2+2i]]
   * x[0] = (-1+2i)/(-12+2i) = (-1+2i)(-12-2i)/(144+4) = (12+2i-24i-4i^2)/(148)
   *       = (12+2i-24i+4)/148 = (16-22i)/148 = (8-11i)/74
   * x[1] = 3/(-12+2i) = 3*(-12-2i)/148 = (-36-6i)/148 = (-18-3i)/74
   */
  sunrealtype x0_re_exp = 8.0 / 74.0;
  sunrealtype x0_im_exp = -11.0 / 74.0;
  sunrealtype x1_re_exp = -18.0 / 74.0;
  sunrealtype x1_im_exp = -3.0 / 74.0;

  CHECK_VAL("solve_x0_re", rhs_re[0], x0_re_exp);
  CHECK_VAL("solve_x0_im", rhs_im[0], x0_im_exp);
  CHECK_VAL("solve_x1_re", rhs_re[1], x1_re_exp);
  CHECK_VAL("solve_x1_im", rhs_im[1], x1_im_exp);

  /* --- Test 3: General mass matrix ---
   * M = [ 2  0 ]
   *     [ 0  3 ]
   * K = (3+2i)*M - J = [(6+4i)-1, -2] = [5+4i, -2]
   *                     [-3, (9+6i)-4]   [-3, 5+6i]
   */
  rmem->M = SUNDenseMatrix(n, n, sunctx);
  SUNMatZero(rmem->M);
  SM_ELEMENT_D(rmem->M, 0, 0) = 2.0;
  SM_ELEMENT_D(rmem->M, 1, 1) = 3.0;

  radau5_BuildE2c(rmem, 0, alphn, betan);

  E2d = SM_DATA_ZD(rmem->E2c_mat[0]);
  CHECK_VAL("E2c_mass[0,0].re", creal(E2d[0]), 5.0);
  CHECK_VAL("E2c_mass[0,0].im", cimag(E2d[0]), 4.0);
  CHECK_VAL("E2c_mass[1,0].re", creal(E2d[1]), -3.0);
  CHECK_VAL("E2c_mass[1,0].im", cimag(E2d[1]), 0.0);
  CHECK_VAL("E2c_mass[0,1].re", creal(E2d[2]), -2.0);
  CHECK_VAL("E2c_mass[0,1].im", cimag(E2d[2]), 0.0);
  CHECK_VAL("E2c_mass[1,1].re", creal(E2d[3]), 5.0);
  CHECK_VAL("E2c_mass[1,1].im", cimag(E2d[3]), 6.0);

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
