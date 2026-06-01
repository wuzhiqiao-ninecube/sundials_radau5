/* Unit test for SUNMatrix_ComplexSparse and SUNLinSol_ComplexSparse (KLU) */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_nvector.h>
#include "nvector_complex_serial.h"
#include "sunmatrix_complex_sparse.h"
#include "sunlinsol_complex_sparse.h"

#define TOL 1.0e-10

int main(void)
{
  SUNContext sunctx;
  int fails = 0;

  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  /* Test with a 4x4 sparse complex matrix (tridiagonal + corner) */
  {
    sunindextype n = 4;
    /* CSC pattern for a tridiagonal + (0,3) and (3,0) entries:
     * Row 0: cols 0, 1, 3
     * Row 1: cols 0, 1, 2
     * Row 2: cols 1, 2, 3
     * Row 3: cols 0, 2, 3
     *
     * CSC representation:
     * col 0: rows 0, 1, 3  (3 entries)
     * col 1: rows 0, 1, 2  (3 entries)
     * col 2: rows 1, 2, 3  (3 entries)
     * col 3: rows 0, 2, 3  (3 entries)
     */
    sunindextype colptrs[] = {0, 3, 6, 9, 12};
    sunindextype rowinds[] = {0, 1, 3,   /* col 0 */
                              0, 1, 2,   /* col 1 */
                              1, 2, 3,   /* col 2 */
                              0, 2, 3};  /* col 3 */
    sunindextype NNZ = 12;

    SUNMatrix A = SUNMatNewFromPattern_ComplexSparse(n, n, colptrs, rowinds,
                                                     sunctx);
    if (A == NULL) { printf("FAIL: SUNMatNewFromPattern\n"); return 1; }

    /* Fill matrix values (diagonally dominant for stability) */
    suncomplextype *Ad = SM_DATA_ZS(A);
    /* col 0: rows 0, 1, 3 */
    Ad[0] = 5.0 + 1.0*I;    /* A(0,0) */
    Ad[1] = -1.0 + 0.5*I;   /* A(1,0) */
    Ad[2] = 0.2 - 0.1*I;    /* A(3,0) */
    /* col 1: rows 0, 1, 2 */
    Ad[3] = -1.0 + 0.3*I;   /* A(0,1) */
    Ad[4] = 6.0 - 2.0*I;    /* A(1,1) */
    Ad[5] = -0.5 + 0.2*I;   /* A(2,1) */
    /* col 2: rows 1, 2, 3 */
    Ad[6] = -0.5 + 0.1*I;   /* A(1,2) */
    Ad[7] = 7.0 + 0.5*I;    /* A(2,2) */
    Ad[8] = -1.0 + 0.3*I;   /* A(3,2) */
    /* col 3: rows 0, 2, 3 */
    Ad[9]  = 0.3 - 0.2*I;   /* A(0,3) */
    Ad[10] = -1.0 + 0.5*I;  /* A(2,3) */
    Ad[11] = 4.0 + 2.0*I;   /* A(3,3) */

    /* Test Matvec */
    N_Vector x_exact = N_VNew_ComplexSerial(n, sunctx);
    N_Vector b = N_VNew_ComplexSerial(n, sunctx);
    N_Vector x = N_VNew_ComplexSerial(n, sunctx);

    suncomplextype *xed = NV_DATA_ZS(x_exact);
    xed[0] = 1.0 + 2.0*I;
    xed[1] = -1.0 + 0.5*I;
    xed[2] = 3.0 - 1.0*I;
    xed[3] = 0.5 + 1.5*I;

    /* b = A * x_exact */
    SUNMatMatvec(A, x_exact, b);

    /* Verify matvec for element 0:
     * b[0] = A(0,0)*x[0] + A(0,1)*x[1] + A(0,3)*x[3]
     */
    suncomplextype b0_exp = Ad[0]*xed[0] + Ad[3]*xed[1] + Ad[9]*xed[3];
    suncomplextype *bd = NV_DATA_ZS(b);
    if (cabs(bd[0] - b0_exp) > TOL)
    {
      printf("FAIL Matvec[0]: err=%.2e\n", cabs(bd[0] - b0_exp));
      fails++;
    }

    /* Create linear solver and solve A*x = b */
    SUNLinearSolver LS = SUNLinSol_ComplexSparse(x, A, sunctx);
    if (LS == NULL) { printf("FAIL: SUNLinSol_ComplexSparse\n"); return 1; }

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

    /* Test refactorization: modify values, re-setup, re-solve */
    Ad[0] = 8.0 + 2.0*I;  /* change A(0,0) */
    Ad[4] = 9.0 - 1.0*I;  /* change A(1,1) */

    /* Recompute b with modified A */
    SUNMatMatvec(A, x_exact, b);

    ret = SUNLinSolSetup(LS, A);
    if (ret != 0)
    {
      printf("FAIL: Refactor Setup returned %d\n", ret);
      fails++;
    }
    else
    {
      ret = SUNLinSolSolve(LS, A, x, b, 0.0);
      if (ret != 0)
      {
        printf("FAIL: Refactor Solve returned %d\n", ret);
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
            printf("FAIL Refactor Solve[%lld]: err=%.2e\n", (long long)i, err);
            fails++;
          }
        }
      }
    }

    /* Test ScaleAddI */
    SUNMatZero(A);
    /* After zero + ScaleAddI(2, A): A should be identity (since 2*0 + I = I) */
    /* But we need diagonal entries in the pattern. They are at indices 0,4,7,11 */
    SUNMatScaleAddI(2.0, A);
    /* Check diagonal */
    if (cabs(Ad[0] - 1.0) > TOL) { printf("FAIL ScaleAddI diag[0]\n"); fails++; }
    if (cabs(Ad[4] - 1.0) > TOL) { printf("FAIL ScaleAddI diag[1]\n"); fails++; }
    if (cabs(Ad[7] - 1.0) > TOL) { printf("FAIL ScaleAddI diag[2]\n"); fails++; }
    if (cabs(Ad[11] - 1.0) > TOL) { printf("FAIL ScaleAddI diag[3]\n"); fails++; }
    /* Check off-diagonal should be 0 */
    if (cabs(Ad[1]) > TOL) { printf("FAIL ScaleAddI off-diag[1]\n"); fails++; }
    if (cabs(Ad[3]) > TOL) { printf("FAIL ScaleAddI off-diag[3]\n"); fails++; }

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
