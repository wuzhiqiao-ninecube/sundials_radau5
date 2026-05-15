/* ---------------------------------------------------------------------------
 * test_radau5_constants.c — Verify Radau IIA method constants
 *
 * Checks that the hardcoded constants (collocation nodes, eigenvalues,
 * eigenvector matrices T/TI, error coefficients) match the Fortran values
 * and satisfy the algebraic identity TI * T = I.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include "radau5.h"
#include "radau5_impl.h"

#define TOL 1.0e-14

static int rhs_zero(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)t; (void)y; (void)ud;
  N_VConst(0.0, yd);
  return 0;
}

static int check(const char* name, sunrealtype got, sunrealtype expected)
{
  sunrealtype err = fabs(got - expected);
  if (err > TOL) {
    fprintf(stderr, "FAIL: %s = %.17e, expected %.17e, diff = %.3e\n",
            name, got, expected, err);
    return 1;
  }
  return 0;
}

int main(void)
{
  int failures = 0;
  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  /* Create solver just to get constants initialized */
  void* mem = Radau5Create(sunctx);
  if (!mem) { fprintf(stderr, "Radau5Create failed\n"); return 1; }

  /* Need a dummy y0 for Init */
  N_Vector y0 = N_VNew_Serial(1, sunctx);
  N_VConst(1.0, y0);

  int ret = Radau5Init(mem, rhs_zero, 0.0, y0);
  if (ret != RADAU5_SUCCESS) {
    fprintf(stderr, "Radau5Init failed: %d\n", ret);
    return 1;
  }

  Radau5Mem rmem = RADAU5_MEM(mem);

  /* --- Collocation nodes (ns=3) --- */
  sunrealtype sq6 = sqrt(6.0);
  failures += check("c[0]", rmem->c[0], (4.0 - sq6) / 10.0);
  failures += check("c[1]", rmem->c[1], (4.0 + sq6) / 10.0);
  failures += check("c[2]", rmem->c[2], 1.0);

  /* --- Eigenvalues of A^{-1} (after inversion) --- */
  sunrealtype u1_raw = (6.0 + pow(81.0, 1.0/3.0) - pow(9.0, 1.0/3.0)) / 30.0;
  sunrealtype alph_raw = (12.0 - pow(81.0, 1.0/3.0) + pow(9.0, 1.0/3.0)) / 60.0;
  sunrealtype beta_raw = (pow(81.0, 1.0/3.0) + pow(9.0, 1.0/3.0)) * sqrt(3.0) / 60.0;
  sunrealtype cno = alph_raw * alph_raw + beta_raw * beta_raw;

  failures += check("u1",       rmem->u1,          1.0 / u1_raw);
  failures += check("alph[0]",  rmem->alph[0],     alph_raw / cno);
  failures += check("beta[0]",  rmem->beta_eig[0], beta_raw / cno);

  /* --- Error coefficients --- */
  failures += check("dd[0]", rmem->dd[0], -(13.0 + 7.0 * sq6) / 3.0);
  failures += check("dd[1]", rmem->dd[1], (-13.0 + 7.0 * sq6) / 3.0);
  failures += check("dd[2]", rmem->dd[2], -1.0 / 3.0);

  /* --- Eigenvector matrix T (row-major T_mat[ns*ns]) --- */
  failures += check("T[0][0]", rmem->T_mat[0],  9.1232394870892942792e-02);
  failures += check("T[0][1]", rmem->T_mat[1], -0.14125529502095420843);
  failures += check("T[0][2]", rmem->T_mat[2], -3.0029194105147424492e-02);
  failures += check("T[1][0]", rmem->T_mat[3],  0.24171793270710701896);
  failures += check("T[1][1]", rmem->T_mat[4],  0.20412935229379993199);
  failures += check("T[1][2]", rmem->T_mat[5],  0.38294211275726193779);
  failures += check("T[2][0]", rmem->T_mat[6],  0.96604818261509293619);

  /* --- Inverse eigenvector matrix TI (row-major TI_mat[ns*ns]) --- */
  failures += check("TI[0][0]", rmem->TI_mat[0],  4.3255798900631553510);
  failures += check("TI[0][1]", rmem->TI_mat[1],  0.33919925181580986954);
  failures += check("TI[0][2]", rmem->TI_mat[2],  0.54177053993587487119);
  failures += check("TI[1][0]", rmem->TI_mat[3], -4.1787185915519047273);
  failures += check("TI[1][1]", rmem->TI_mat[4], -0.32768282076106238708);
  failures += check("TI[1][2]", rmem->TI_mat[5],  0.47662355450055045196);
  failures += check("TI[2][0]", rmem->TI_mat[6], -0.50287263494578687595);
  failures += check("TI[2][1]", rmem->TI_mat[7],  2.5719269498556054292);
  failures += check("TI[2][2]", rmem->TI_mat[8], -0.59603920482822492497);

  /* --- Verify TI * T = I (3x3 identity) ---
   * T has T[2][1]=1, T[2][2]=0 implicit in the Fortran. */
  int ns = rmem->ns;
  sunrealtype T[3][3], TI[3][3];
  for (int r = 0; r < ns; r++)
    for (int c = 0; c < ns; c++) {
      T[r][c]  = rmem->T_mat[r * ns + c];
      TI[r][c] = rmem->TI_mat[r * ns + c];
    }

  for (int r = 0; r < ns; r++) {
    for (int c = 0; c < ns; c++) {
      sunrealtype sum = 0.0;
      for (int k = 0; k < ns; k++)
        sum += TI[r][k] * T[k][c];
      sunrealtype expected = (r == c) ? 1.0 : 0.0;
      char buf[64];
      sprintf(buf, "TI*T[%d][%d]", r, c);
      failures += check(buf, sum, expected);
    }
  }

  /* Cleanup */
  N_VDestroy(y0);
  Radau5Free(&mem);
  SUNContext_Free(&sunctx);

  if (failures == 0)
    printf("test_radau5_constants: PASSED\n");
  else
    printf("test_radau5_constants: %d FAILURES\n", failures);

  return failures;
}