/* Unit test for SUNMatrix_ComplexDense and SUNLinSol_ComplexDense */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_nvector.h>
#include "nvector_complex_serial.h"
#include "sunmatrix_complex_dense.h"
#include "sunlinsol_complex_dense.h"

#define TOL 1.0e-12

static int check_val(const char *name, double computed, double expected)
{
  double err = fabs(computed - expected);
  if (err > TOL)
  {
    printf("FAIL %s: computed=%.16e expected=%.16e err=%.2e\n",
           name, computed, expected, err);
    return 1;
  }
  return 0;
}

int main(void)
{
  SUNContext sunctx;
  int fails = 0;

  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  /* ---- Test SUNMatrix_ComplexDense ---- */
  {
    sunindextype M = 3, N = 3;
    SUNMatrix A = SUNMatNew_ComplexDense(M, N, sunctx);
    if (A == NULL) { printf("FAIL: SUNMatNew_ComplexDense\n"); return 1; }

    /* Check dimensions */
    if (SM_ROWS_ZD(A) != 3 || SM_COLUMNS_ZD(A) != 3)
    {
      printf("FAIL: dimensions\n"); fails++;
    }

    /* Fill A with known values:
     * A = [ 1+i,   2+i,   3    ]
     *     [ 4,     5+2i,  6-i  ]
     *     [ 7-i,   8,     9+3i ]
     */
    SM_ELEMENT_ZD(A, 0, 0) = 1.0 + 1.0*I;
    SM_ELEMENT_ZD(A, 1, 0) = 4.0;
    SM_ELEMENT_ZD(A, 2, 0) = 7.0 - 1.0*I;
    SM_ELEMENT_ZD(A, 0, 1) = 2.0 + 1.0*I;
    SM_ELEMENT_ZD(A, 1, 1) = 5.0 + 2.0*I;
    SM_ELEMENT_ZD(A, 2, 1) = 8.0;
    SM_ELEMENT_ZD(A, 0, 2) = 3.0;
    SM_ELEMENT_ZD(A, 1, 2) = 6.0 - 1.0*I;
    SM_ELEMENT_ZD(A, 2, 2) = 9.0 + 3.0*I;

    /* Test Matvec: y = A*x */
    N_Vector x = N_VNew_ComplexSerial(3, sunctx);
    N_Vector y = N_VNew_ComplexSerial(3, sunctx);
    suncomplextype *xd = NV_DATA_ZS(x);
    suncomplextype *yd;

    xd[0] = 1.0;
    xd[1] = 0.0 + 1.0*I;
    xd[2] = -1.0 + 0.5*I;

    SUNMatMatvec(A, x, y);
    yd = NV_DATA_ZS(y);

    /* Verify y = A*x manually:
     * y[0] = (1+i)*1 + (2+i)*(i) + 3*(-1+0.5i)
     *       = 1+i + 2i+i^2 + -3+1.5i = 1+i + 2i-1 + -3+1.5i = -3 + 4.5i
     */
    suncomplextype y0_exp = (1.0+1.0*I)*1.0 + (2.0+1.0*I)*(1.0*I) + 3.0*(-1.0+0.5*I);
    suncomplextype y1_exp = 4.0*1.0 + (5.0+2.0*I)*(1.0*I) + (6.0-1.0*I)*(-1.0+0.5*I);
    suncomplextype y2_exp = (7.0-1.0*I)*1.0 + 8.0*(1.0*I) + (9.0+3.0*I)*(-1.0+0.5*I);

    if (cabs(yd[0] - y0_exp) > TOL)
    {
      printf("FAIL Matvec[0]: got (%.6f,%.6f) exp (%.6f,%.6f)\n",
             creal(yd[0]), cimag(yd[0]), creal(y0_exp), cimag(y0_exp));
      fails++;
    }
    if (cabs(yd[1] - y1_exp) > TOL)
    {
      printf("FAIL Matvec[1]: got (%.6f,%.6f) exp (%.6f,%.6f)\n",
             creal(yd[1]), cimag(yd[1]), creal(y1_exp), cimag(y1_exp));
      fails++;
    }
    if (cabs(yd[2] - y2_exp) > TOL)
    {
      printf("FAIL Matvec[2]: got (%.6f,%.6f) exp (%.6f,%.6f)\n",
             creal(yd[2]), cimag(yd[2]), creal(y2_exp), cimag(y2_exp));
      fails++;
    }

    /* Test HermitianTransposeVec: y = A^H * x */
    SUNMatHermitianTransposeVec_ComplexDense(A, x, y);
    yd = NV_DATA_ZS(y);

    suncomplextype yh0_exp = conj(1.0+1.0*I)*1.0 + conj(4.0)*(1.0*I) + conj(7.0-1.0*I)*(-1.0+0.5*I);
    suncomplextype yh1_exp = conj(2.0+1.0*I)*1.0 + conj(5.0+2.0*I)*(1.0*I) + conj(8.0)*(-1.0+0.5*I);
    suncomplextype yh2_exp = conj(3.0)*1.0 + conj(6.0-1.0*I)*(1.0*I) + conj(9.0+3.0*I)*(-1.0+0.5*I);

    if (cabs(yd[0] - yh0_exp) > TOL)
    {
      printf("FAIL HermTransVec[0]: err=%.2e\n", cabs(yd[0] - yh0_exp));
      fails++;
    }
    if (cabs(yd[1] - yh1_exp) > TOL)
    {
      printf("FAIL HermTransVec[1]: err=%.2e\n", cabs(yd[1] - yh1_exp));
      fails++;
    }
    if (cabs(yd[2] - yh2_exp) > TOL)
    {
      printf("FAIL HermTransVec[2]: err=%.2e\n", cabs(yd[2] - yh2_exp));
      fails++;
    }

    /* Test Zero */
    SUNMatZero(A);
    suncomplextype *Ad = SM_DATA_ZD(A);
    for (sunindextype k = 0; k < 9; k++)
    {
      if (cabs(Ad[k]) > 0.0)
      {
        printf("FAIL Zero[%lld]\n", (long long)k);
        fails++;
      }
    }

    /* Test ScaleAddI: A = 2*A + I (A is zero, so result is I) */
    SUNMatScaleAddI(2.0, A);
    for (sunindextype i = 0; i < 3; i++)
    {
      for (sunindextype j = 0; j < 3; j++)
      {
        suncomplextype expected = (i == j) ? 1.0 : 0.0;
        if (cabs(SM_ELEMENT_ZD(A, i, j) - expected) > TOL)
        {
          printf("FAIL ScaleAddI[%lld,%lld]\n", (long long)i, (long long)j);
          fails++;
        }
      }
    }

    N_VDestroy(x);
    N_VDestroy(y);
    SUNMatDestroy(A);
  }

  /* ---- Test SUNLinSol_ComplexDense (solve A*x = b) ---- */
  {
    sunindextype n = 4;
    SUNMatrix A = SUNMatNew_ComplexDense(n, n, sunctx);
    N_Vector x = N_VNew_ComplexSerial(n, sunctx);
    N_Vector b = N_VNew_ComplexSerial(n, sunctx);
    N_Vector x_exact = N_VNew_ComplexSerial(n, sunctx);

    /* Set known solution */
    suncomplextype *xed = NV_DATA_ZS(x_exact);
    xed[0] = 1.0 + 2.0*I;
    xed[1] = -1.0 + 0.5*I;
    xed[2] = 3.0 - 1.0*I;
    xed[3] = 0.5 + 0.5*I;

    /* Fill A with a well-conditioned complex matrix */
    SM_ELEMENT_ZD(A, 0, 0) = 4.0 + 1.0*I;
    SM_ELEMENT_ZD(A, 0, 1) = 1.0 - 1.0*I;
    SM_ELEMENT_ZD(A, 0, 2) = 0.5;
    SM_ELEMENT_ZD(A, 0, 3) = 0.0 + 0.5*I;
    SM_ELEMENT_ZD(A, 1, 0) = 1.0 + 1.0*I;
    SM_ELEMENT_ZD(A, 1, 1) = 5.0 - 2.0*I;
    SM_ELEMENT_ZD(A, 1, 2) = 0.0 + 1.0*I;
    SM_ELEMENT_ZD(A, 1, 3) = 1.0;
    SM_ELEMENT_ZD(A, 2, 0) = 0.5 - 0.5*I;
    SM_ELEMENT_ZD(A, 2, 1) = 0.0 + 1.0*I;
    SM_ELEMENT_ZD(A, 2, 2) = 6.0 + 1.0*I;
    SM_ELEMENT_ZD(A, 2, 3) = 1.0 - 1.0*I;
    SM_ELEMENT_ZD(A, 3, 0) = 0.0 - 0.5*I;
    SM_ELEMENT_ZD(A, 3, 1) = 1.0;
    SM_ELEMENT_ZD(A, 3, 2) = 1.0 + 1.0*I;
    SM_ELEMENT_ZD(A, 3, 3) = 7.0 + 3.0*I;

    /* Compute b = A * x_exact */
    SUNMatMatvec(A, x_exact, b);

    /* Create linear solver */
    SUNLinearSolver LS = SUNLinSol_ComplexDense(x, A, sunctx);
    if (LS == NULL) { printf("FAIL: SUNLinSol_ComplexDense\n"); return 1; }

    /* Setup (factorize) */
    int ret = SUNLinSolSetup(LS, A);
    if (ret != 0)
    {
      printf("FAIL: SUNLinSolSetup returned %d\n", ret);
      fails++;
    }
    else
    {
      /* Solve */
      ret = SUNLinSolSolve(LS, A, x, b, 0.0);
      if (ret != 0)
      {
        printf("FAIL: SUNLinSolSolve returned %d\n", ret);
        fails++;
      }
      else
      {
        /* Check solution */
        suncomplextype *xd = NV_DATA_ZS(x);
        for (sunindextype i = 0; i < n; i++)
        {
          double err = cabs(xd[i] - xed[i]);
          if (err > 1.0e-10)
          {
            printf("FAIL Solve[%lld]: got (%.10f,%.10f) exp (%.10f,%.10f) err=%.2e\n",
                   (long long)i, creal(xd[i]), cimag(xd[i]),
                   creal(xed[i]), cimag(xed[i]), err);
            fails++;
          }
        }
      }
    }

    SUNLinSolFree(LS);
    N_VDestroy(x);
    N_VDestroy(b);
    N_VDestroy(x_exact);
    SUNMatDestroy(A);
  }

  SUNContext_Free(&sunctx);

  if (fails == 0)
  {
    printf("ALL TESTS PASSED\n");
  }
  else
  {
    printf("%d TEST(S) FAILED\n", fails);
  }

  return fails;
}
