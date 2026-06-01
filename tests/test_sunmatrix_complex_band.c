/* Unit test for SUNMatrix_ComplexBand and SUNLinSol_ComplexBand */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_nvector.h>
#include "nvector_complex_serial.h"
#include "sunmatrix_complex_band.h"
#include "sunlinsol_complex_band.h"

#define TOL 1.0e-10

int main(void)
{
  SUNContext sunctx;
  int fails = 0;

  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  /* Test with a 5x5 tridiagonal matrix (mu=1, ml=1) */
  {
    sunindextype n = 5, mu = 1, ml = 1;
    SUNMatrix A = SUNMatNew_ComplexBand(n, mu, ml, sunctx);
    N_Vector x = N_VNew_ComplexSerial(n, sunctx);
    N_Vector b = N_VNew_ComplexSerial(n, sunctx);
    N_Vector x_exact = N_VNew_ComplexSerial(n, sunctx);

    if (A == NULL) { printf("FAIL: SUNMatNew_ComplexBand\n"); return 1; }

    /* Fill tridiagonal: diagonal = 4+i, sub/super = -1+0.5i */
    sunindextype j;
    for (j = 0; j < n; j++)
    {
      SM_ELEMENT_ZB(A, j, j) = 4.0 + 1.0*I;
      if (j > 0)   { SM_ELEMENT_ZB(A, j, j-1) = -1.0 + 0.5*I; }
      if (j < n-1) { SM_ELEMENT_ZB(A, j, j+1) = -1.0 + 0.5*I; }
    }

    /* Set known solution */
    suncomplextype *xed = NV_DATA_ZS(x_exact);
    xed[0] = 1.0 + 2.0*I;
    xed[1] = -1.0 + 0.5*I;
    xed[2] = 3.0 - 1.0*I;
    xed[3] = 0.5 + 1.5*I;
    xed[4] = -2.0 + 0.0*I;

    /* Test Matvec: b = A * x_exact */
    SUNMatMatvec(A, x_exact, b);

    /* Verify matvec manually for element 0:
     * b[0] = (4+i)*x[0] + (-1+0.5i)*x[1]
     */
    suncomplextype b0_exp = (4.0+1.0*I)*xed[0] + (-1.0+0.5*I)*xed[1];
    suncomplextype *bd = NV_DATA_ZS(b);
    if (cabs(bd[0] - b0_exp) > TOL)
    {
      printf("FAIL Matvec[0]: err=%.2e\n", cabs(bd[0] - b0_exp));
      fails++;
    }

    /* Create linear solver and solve */
    SUNLinearSolver LS = SUNLinSol_ComplexBand(x, A, sunctx);
    if (LS == NULL) { printf("FAIL: SUNLinSol_ComplexBand\n"); return 1; }

    int ret = SUNLinSolSetup(LS, A);
    if (ret != 0)
    {
      printf("FAIL: SUNLinSolSetup returned %d\n", ret);
      fails++;
    }
    else
    {
      ret = SUNLinSolSolve(LS, A, x, b, 0.0);
      if (ret != 0)
      {
        printf("FAIL: SUNLinSolSolve returned %d\n", ret);
        fails++;
      }
      else
      {
        suncomplextype *xd = NV_DATA_ZS(x);
        for (sunindextype i = 0; i < n; i++)
        {
          double err = cabs(xd[i] - xed[i]);
          if (err > TOL)
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

  /* Test with a wider band (mu=2, ml=1) */
  {
    sunindextype n = 4, mu = 2, ml = 1;
    SUNMatrix A = SUNMatNew_ComplexBand(n, mu, ml, sunctx);
    N_Vector x = N_VNew_ComplexSerial(n, sunctx);
    N_Vector b = N_VNew_ComplexSerial(n, sunctx);
    N_Vector x_exact = N_VNew_ComplexSerial(n, sunctx);

    /* Fill band matrix */
    SM_ELEMENT_ZB(A, 0, 0) = 5.0 + 1.0*I;
    SM_ELEMENT_ZB(A, 0, 1) = 1.0 - 0.5*I;
    SM_ELEMENT_ZB(A, 0, 2) = 0.5 + 0.2*I;
    SM_ELEMENT_ZB(A, 1, 0) = -1.0 + 0.3*I;
    SM_ELEMENT_ZB(A, 1, 1) = 6.0 - 2.0*I;
    SM_ELEMENT_ZB(A, 1, 2) = 1.0;
    SM_ELEMENT_ZB(A, 1, 3) = 0.3 + 0.1*I;
    SM_ELEMENT_ZB(A, 2, 1) = -0.5 + 0.1*I;
    SM_ELEMENT_ZB(A, 2, 2) = 7.0 + 0.5*I;
    SM_ELEMENT_ZB(A, 2, 3) = 1.0 - 0.5*I;
    SM_ELEMENT_ZB(A, 3, 2) = -1.0 + 0.2*I;
    SM_ELEMENT_ZB(A, 3, 3) = 4.0 + 2.0*I;

    suncomplextype *xed = NV_DATA_ZS(x_exact);
    xed[0] = 1.0 + 1.0*I;
    xed[1] = 2.0 - 0.5*I;
    xed[2] = -1.0 + 2.0*I;
    xed[3] = 0.5 - 1.0*I;

    SUNMatMatvec(A, x_exact, b);

    SUNLinearSolver LS = SUNLinSol_ComplexBand(x, A, sunctx);
    int ret = SUNLinSolSetup(LS, A);
    if (ret != 0)
    {
      printf("FAIL: Band(2,1) Setup returned %d\n", ret);
      fails++;
    }
    else
    {
      ret = SUNLinSolSolve(LS, A, x, b, 0.0);
      if (ret != 0)
      {
        printf("FAIL: Band(2,1) Solve returned %d\n", ret);
        fails++;
      }
      else
      {
        suncomplextype *xd = NV_DATA_ZS(x);
        for (sunindextype i = 0; i < n; i++)
        {
          double err = cabs(xd[i] - xed[i]);
          if (err > TOL)
          {
            printf("FAIL Band(2,1) Solve[%lld]: err=%.2e\n", (long long)i, err);
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

  /* Test ScaleAddI */
  {
    sunindextype n = 3, mu = 1, ml = 1;
    SUNMatrix A = SUNMatNew_ComplexBand(n, mu, ml, sunctx);

    /* Set A to zero, then ScaleAddI(2, A) should give I */
    SUNMatZero(A);
    SUNMatScaleAddI(2.0, A);

    for (sunindextype j = 0; j < n; j++)
    {
      if (cabs(SM_ELEMENT_ZB(A, j, j) - 1.0) > TOL)
      {
        printf("FAIL ScaleAddI diag[%lld]\n", (long long)j);
        fails++;
      }
    }

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
