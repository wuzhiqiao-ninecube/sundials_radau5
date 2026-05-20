/* ---------------------------------------------------------------------------
 * test_radau5_constants_high.c — Verify Radau IIA method constants for ns=9,11,13
 *
 * For each ns in {9, 11, 13}, checks:
 *   1. Collocation nodes c[i] sum to known value
 *   2. Last collocation node c[ns-1] == 1.0
 *   3. Error coefficients dd[ns-1] == -1/ns
 *   4. TI * T = I (identity matrix)
 *   5. Eigenvalue consistency: u1 matches expected value
 *   6. ALPH/BETA array sizes match npairs = (ns-1)/2
 *   7. Schur mode: US * US^T = I (orthogonality)
 *   8. Schur mode: TS diagonal matches u1
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include "radau5.h"
#include "radau5_impl.h"

#define TOL 1.0e-10
#define TOL_ORTHO 1.0e-9

static int rhs_zero(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)t; (void)y; (void)ud;
  N_VConst(0.0, yd);
  return 0;
}

static int check(const char* name, sunrealtype got, sunrealtype expected, sunrealtype tol)
{
  sunrealtype err = fabs(got - expected);
  sunrealtype scale = fabs(expected) > 1.0 ? fabs(expected) : 1.0;
  if (err / scale > tol) {
    fprintf(stderr, "FAIL: %s = %.17e, expected %.17e, rel_err = %.3e\n",
            name, got, expected, err / scale);
    return 1;
  }
  return 0;
}

/* Test TI * T = I for given ns */
static int test_TI_T_identity(Radau5Mem rmem, int ns)
{
  int failures = 0;
  for (int r = 0; r < ns; r++) {
    for (int c = 0; c < ns; c++) {
      sunrealtype sum = 0.0;
      for (int k = 0; k < ns; k++)
        sum += rmem->TI_mat[r * ns + k] * rmem->T_mat[k * ns + c];
      sunrealtype expected = (r == c) ? 1.0 : 0.0;
      char buf[64];
      sprintf(buf, "ns=%d TI*T[%d][%d]", ns, r, c);
      failures += check(buf, sum, expected, TOL);
    }
  }
  return failures;
}

/* Test US * US^T = I (orthogonality) for given ns */
static int test_US_orthogonality(Radau5Mem rmem, int ns)
{
  int failures = 0;
  for (int r = 0; r < ns; r++) {
    for (int c = 0; c < ns; c++) {
      sunrealtype sum = 0.0;
      for (int k = 0; k < ns; k++)
        sum += rmem->US_mat[r * ns + k] * rmem->US_mat[c * ns + k];
      sunrealtype expected = (r == c) ? 1.0 : 0.0;
      char buf[64];
      sprintf(buf, "ns=%d US*US^T[%d][%d]", ns, r, c);
      failures += check(buf, sum, expected, TOL_ORTHO);
    }
  }
  return failures;
}

/* Test constants for a given ns value in eigenvalue mode */
static int test_ns_eigen(int ns_val, SUNContext sunctx)
{
  int failures = 0;
  printf("  Testing ns=%d (eigenvalue mode)...\n", ns_val);

  void* mem = Radau5Create(sunctx);
  if (!mem) { fprintf(stderr, "Radau5Create failed\n"); return 1; }

  Radau5SetNumStages(mem, ns_val);

  N_Vector y0 = N_VNew_Serial(1, sunctx);
  N_VConst(1.0, y0);
  int ret = Radau5Init(mem, rhs_zero, 0.0, y0);
  if (ret != RADAU5_SUCCESS) {
    fprintf(stderr, "Radau5Init failed for ns=%d: %d\n", ns_val, ret);
    N_VDestroy(y0);
    Radau5Free(&mem);
    return 1;
  }

  Radau5Mem rmem = RADAU5_MEM(mem);
  int npairs = (ns_val - 1) / 2;

  /* Check c[ns-1] == 1.0 */
  {
    char buf[64];
    sprintf(buf, "ns=%d c[%d]", ns_val, ns_val - 1);
    failures += check(buf, rmem->c[ns_val - 1], 1.0, 1e-15);
  }

  /* Check dd[ns-1] == -1/ns */
  {
    char buf[64];
    sprintf(buf, "ns=%d dd[%d]", ns_val, ns_val - 1);
    failures += check(buf, rmem->dd[ns_val - 1], -1.0 / (sunrealtype)ns_val, 1e-14);
  }

  /* Check collocation nodes are in (0, 1] and increasing */
  for (int i = 0; i < ns_val; i++) {
    if (rmem->c[i] <= 0.0 || rmem->c[i] > 1.0) {
      fprintf(stderr, "FAIL: ns=%d c[%d] = %.17e not in (0,1]\n",
              ns_val, i, rmem->c[i]);
      failures++;
    }
    if (i > 0 && rmem->c[i] <= rmem->c[i-1]) {
      fprintf(stderr, "FAIL: ns=%d c[%d] <= c[%d]\n", ns_val, i, i-1);
      failures++;
    }
  }

  /* Check u1 > 0 */
  if (rmem->u1 <= 0.0) {
    fprintf(stderr, "FAIL: ns=%d u1 = %.17e <= 0\n", ns_val, rmem->u1);
    failures++;
  }

  /* Check alph[k] > 0 and beta[k] > 0 */
  for (int k = 0; k < npairs; k++) {
    if (rmem->alph[k] <= 0.0) {
      fprintf(stderr, "FAIL: ns=%d alph[%d] = %.17e <= 0\n", ns_val, k, rmem->alph[k]);
      failures++;
    }
    if (rmem->beta_eig[k] <= 0.0) {
      fprintf(stderr, "FAIL: ns=%d beta[%d] = %.17e <= 0\n", ns_val, k, rmem->beta_eig[k]);
      failures++;
    }
  }

  /* Check TI * T = I */
  failures += test_TI_T_identity(rmem, ns_val);

  /* Check T last row pattern: T[ns-1][0]=val, T[ns-1][1]=1, T[ns-1][2]=0, ... */
  for (int j = 1; j < ns_val; j++) {
    sunrealtype expected = (j % 2 == 1) ? 1.0 : 0.0;
    char buf[64];
    sprintf(buf, "ns=%d T[%d][%d]", ns_val, ns_val - 1, j);
    failures += check(buf, rmem->T_mat[(ns_val - 1) * ns_val + j], expected, 1e-15);
  }

  N_VDestroy(y0);
  Radau5Free(&mem);
  return failures;
}

/* Test constants for a given ns value in Schur mode */
static int test_ns_schur(int ns_val, SUNContext sunctx)
{
  int failures = 0;
  printf("  Testing ns=%d (Schur mode)...\n", ns_val);

  void* mem = Radau5Create(sunctx);
  if (!mem) { fprintf(stderr, "Radau5Create failed\n"); return 1; }

  Radau5SetNumStages(mem, ns_val);
  Radau5SetSchurDecomp(mem, 1);

  N_Vector y0 = N_VNew_Serial(1, sunctx);
  N_VConst(1.0, y0);
  int ret = Radau5Init(mem, rhs_zero, 0.0, y0);
  if (ret != RADAU5_SUCCESS) {
    fprintf(stderr, "Radau5Init failed for ns=%d schur: %d\n", ns_val, ret);
    N_VDestroy(y0);
    Radau5Free(&mem);
    return 1;
  }

  Radau5Mem rmem = RADAU5_MEM(mem);

  /* Check US orthogonality: US * US^T = I */
  failures += test_US_orthogonality(rmem, ns_val);

  /* Check T = US (in Schur mode) */
  for (int i = 0; i < ns_val * ns_val; i++) {
    if (fabs(rmem->T_mat[i] - rmem->US_mat[i]) > 1e-15) {
      fprintf(stderr, "FAIL: ns=%d T_mat[%d] != US_mat[%d]\n", ns_val, i, i);
      failures++;
      break;
    }
  }

  /* Check TI = US^T (transpose) */
  for (int r = 0; r < ns_val; r++) {
    for (int c = 0; c < ns_val; c++) {
      sunrealtype ti_val = rmem->TI_mat[r * ns_val + c];
      sunrealtype us_t_val = rmem->US_mat[c * ns_val + r];
      if (fabs(ti_val - us_t_val) > 1e-15) {
        fprintf(stderr, "FAIL: ns=%d TI[%d][%d] != US^T[%d][%d]\n",
                ns_val, r, c, r, c);
        failures++;
        break;
      }
    }
    if (failures) break;
  }

  /* Check u1 = TS[ns-1][ns-1] */
  {
    sunrealtype ts_diag = rmem->TS_mat[(ns_val - 1) * ns_val + (ns_val - 1)];
    char buf[64];
    sprintf(buf, "ns=%d u1==TS[%d][%d]", ns_val, ns_val - 1, ns_val - 1);
    failures += check(buf, rmem->u1, ts_diag, 1e-15);
  }

  /* Check TS is upper quasi-triangular: below the 2x2 blocks should be zero */
  /* The structure is: pairs of 2x2 blocks on diagonal, then 1x1 block at bottom-right */
  int npairs = (ns_val - 1) / 2;
  /* Check zeros below the block structure */
  for (int r = 2; r < ns_val; r++) {
    int block_start = (r < ns_val - 1) ? (r / 2) * 2 : ns_val - 1;
    for (int c = 0; c < block_start - 1; c++) {
      sunrealtype val = rmem->TS_mat[r * ns_val + c];
      if (fabs(val) > 1e-12) {
        fprintf(stderr, "FAIL: ns=%d TS[%d][%d] = %.3e (expected 0)\n",
                ns_val, r, c, val);
        failures++;
      }
    }
  }

  N_VDestroy(y0);
  Radau5Free(&mem);
  return failures;
}

/* Test specific known values for ns=9 */
static int test_ns9_values(SUNContext sunctx)
{
  int failures = 0;
  printf("  Checking ns=9 specific values...\n");

  void* mem = Radau5Create(sunctx);
  Radau5SetNumStages(mem, 9);
  N_Vector y0 = N_VNew_Serial(1, sunctx);
  N_VConst(1.0, y0);
  Radau5Init(mem, rhs_zero, 0.0, y0);
  Radau5Mem rmem = RADAU5_MEM(mem);

  failures += check("ns9 c[0]", rmem->c[0], 0.177799151473634518132051010376790613e-01, 1e-14);
  failures += check("ns9 c[7]", rmem->c[7], 0.955366044710030149266878978141692238e+00, 1e-14);
  failures += check("ns9 u1", rmem->u1, 0.115873509212862784070515753302441063e+02, 1e-13);
  failures += check("ns9 alph[0]", rmem->alph[0], 0.496612926068677791126709775546702135e+01, 1e-13);
  failures += check("ns9 beta[3]", rmem->beta_eig[3], 0.332134053152182946011157633342576940e+01, 1e-13);

  N_VDestroy(y0);
  Radau5Free(&mem);
  return failures;
}

/* Test specific known values for ns=11 */
static int test_ns11_values(SUNContext sunctx)
{
  int failures = 0;
  printf("  Checking ns=11 specific values...\n");

  void* mem = Radau5Create(sunctx);
  Radau5SetNumStages(mem, 11);
  N_Vector y0 = N_VNew_Serial(1, sunctx);
  N_VConst(1.0, y0);
  Radau5Init(mem, rhs_zero, 0.0, y0);
  Radau5Mem rmem = RADAU5_MEM(mem);

  failures += check("ns11 c[0]", rmem->c[0], 0.119176134324155969097455869589859760e-01, 1e-14);
  failures += check("ns11 c[9]", rmem->c[9], 0.969970967838513502956935642365592059e+00, 1e-14);
  failures += check("ns11 u1", rmem->u1, 0.142380399544621108935030412346331511e+02, 1e-13);
  failures += check("ns11 alph[4]", rmem->alph[4], 0.139626435483485822564042721435663498e+02, 1e-13);
  failures += check("ns11 beta[0]", rmem->beta_eig[0], 0.176032980318069105751033443311784989e+02, 1e-13);

  N_VDestroy(y0);
  Radau5Free(&mem);
  return failures;
}

/* Test specific known values for ns=13 */
static int test_ns13_values(SUNContext sunctx)
{
  int failures = 0;
  printf("  Checking ns=13 specific values...\n");

  void* mem = Radau5Create(sunctx);
  Radau5SetNumStages(mem, 13);
  N_Vector y0 = N_VNew_Serial(1, sunctx);
  N_VConst(1.0, y0);
  Radau5Init(mem, rhs_zero, 0.0, y0);
  Radau5Mem rmem = RADAU5_MEM(mem);

  failures += check("ns13 c[0]", rmem->c[0], 0.853905498842741936866446087783980280e-02, 1e-14);
  failures += check("ns13 c[11]", rmem->c[11], 0.978437936834149639091906916916996038e+00, 1e-14);
  failures += check("ns13 u1", rmem->u1, 0.168888189439781927912425829292606156e+02, 1e-13);
  failures += check("ns13 alph[5]", rmem->alph[5], 0.166544961771492438512633627168438846e+02, 1e-13);
  failures += check("ns13 beta[5]", rmem->beta_eig[5], 0.336581446671059695895140239139227454e+01, 1e-13);

  N_VDestroy(y0);
  Radau5Free(&mem);
  return failures;
}

int main(void)
{
  int failures = 0;
  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  printf("test_radau5_constants_high: ns=9, 11, 13 constants\n");

  /* Eigenvalue mode tests */
  failures += test_ns_eigen(9, sunctx);
  failures += test_ns_eigen(11, sunctx);
  failures += test_ns_eigen(13, sunctx);

  /* Schur mode tests */
  failures += test_ns_schur(9, sunctx);
  failures += test_ns_schur(11, sunctx);
  failures += test_ns_schur(13, sunctx);

  /* Specific value checks */
  failures += test_ns9_values(sunctx);
  failures += test_ns11_values(sunctx);
  failures += test_ns13_values(sunctx);

  SUNContext_Free(&sunctx);

  if (failures == 0)
    printf("test_radau5_constants_high: PASSED\n");
  else
    printf("test_radau5_constants_high: %d FAILURES\n", failures);

  return failures;
}
