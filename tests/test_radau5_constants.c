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

  /* --- Collocation nodes --- */
  sunrealtype sq6 = sqrt(6.0);
  failures += check("c1", rmem->c1, (4.0 - sq6) / 10.0);
  failures += check("c2", rmem->c2, (4.0 + sq6) / 10.0);
  failures += check("c1m1", rmem->c1m1, rmem->c1 - 1.0);
  failures += check("c2m1", rmem->c2m1, rmem->c2 - 1.0);
  failures += check("c1mc2", rmem->c1mc2, rmem->c1 - rmem->c2);

  /* --- Eigenvalues of A^{-1} (after inversion) --- */
  sunrealtype u1_raw = (6.0 + pow(81.0, 1.0/3.0) - pow(9.0, 1.0/3.0)) / 30.0;
  sunrealtype alph_raw = (12.0 - pow(81.0, 1.0/3.0) + pow(9.0, 1.0/3.0)) / 60.0;
  sunrealtype beta_raw = (pow(81.0, 1.0/3.0) + pow(9.0, 1.0/3.0)) * sqrt(3.0) / 60.0;
  sunrealtype cno = alph_raw * alph_raw + beta_raw * beta_raw;

  failures += check("u1",   rmem->u1,   1.0 / u1_raw);
  failures += check("alph", rmem->alph, alph_raw / cno);
  failures += check("beta", rmem->beta, beta_raw / cno);

  /* --- Error coefficients --- */
  failures += check("dd1", rmem->dd1, -(13.0 + 7.0 * sq6) / 3.0);
  failures += check("dd2", rmem->dd2, (-13.0 + 7.0 * sq6) / 3.0);
  failures += check("dd3", rmem->dd3, -1.0 / 3.0);

  /* --- Eigenvector matrix T (Fortran hardcoded values) --- */
  failures += check("T11", rmem->T11,  9.1232394870892942792e-02);
  failures += check("T12", rmem->T12, -0.14125529502095420843);
  failures += check("T13", rmem->T13, -3.0029194105147424492e-02);
  failures += check("T21", rmem->T21,  0.24171793270710701896);
  failures += check("T22", rmem->T22,  0.20412935229379993199);
  failures += check("T23", rmem->T23,  0.38294211275726193779);
  failures += check("T31", rmem->T31,  0.96604818261509293619);

  /* --- Inverse eigenvector matrix TI --- */
  failures += check("TI11", rmem->TI11,  4.3255798900631553510);
  failures += check("TI12", rmem->TI12,  0.33919925181580986954);
  failures += check("TI13", rmem->TI13,  0.54177053993587487119);
  failures += check("TI21", rmem->TI21, -4.1787185915519047273);
  failures += check("TI22", rmem->TI22, -0.32768282076106238708);
  failures += check("TI23", rmem->TI23,  0.47662355450055045196);
  failures += check("TI31", rmem->TI31, -0.50287263494578687595);
  failures += check("TI32", rmem->TI32,  2.5719269498556054292);
  failures += check("TI33", rmem->TI33, -0.59603920482822492497);

  /* --- Verify TI * T = I (3x3 identity) ---
   * T has T32=1, T33=0 implicit in the Fortran.
   * T = [ T11  T12  T13 ]    TI = [ TI11  TI12  TI13 ]
   *     [ T21  T22  T23 ]         [ TI21  TI22  TI23 ]
   *     [ T31   1    0  ]         [ TI31  TI32  TI33 ]
   */
  sunrealtype T[3][3] = {
    { rmem->T11, rmem->T12, rmem->T13 },
    { rmem->T21, rmem->T22, rmem->T23 },
    { rmem->T31, 1.0,       0.0       }
  };
  sunrealtype TI[3][3] = {
    { rmem->TI11, rmem->TI12, rmem->TI13 },
    { rmem->TI21, rmem->TI22, rmem->TI23 },
    { rmem->TI31, rmem->TI32, rmem->TI33 }
  };

  for (int r = 0; r < 3; r++) {
    for (int c = 0; c < 3; c++) {
      sunrealtype sum = 0.0;
      for (int k = 0; k < 3; k++)
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
