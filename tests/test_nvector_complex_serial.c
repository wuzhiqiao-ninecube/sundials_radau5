/* Unit test for N_Vector_ComplexSerial */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_nvector.h>
#include <nvector/nvector_serial.h>
#include "nvector_complex_serial.h"

#define N 5
#define TOL 1.0e-14

static int check_val(const char *name, sunrealtype computed, sunrealtype expected)
{
  sunrealtype err = fabs(computed - expected);
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
  N_Vector x, y, z;
  suncomplextype *xd, *yd, *zd;
  int fails = 0;
  sunindextype i;

  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  /* Create vectors */
  x = N_VNew_ComplexSerial(N, sunctx);
  y = N_VNew_ComplexSerial(N, sunctx);
  z = N_VNew_ComplexSerial(N, sunctx);

  if (x == NULL || y == NULL || z == NULL)
  {
    printf("FAIL: N_VNew_ComplexSerial returned NULL\n");
    return 1;
  }

  /* Check length */
  if (N_VGetLength(x) != N)
  {
    printf("FAIL: N_VGetLength returned %lld, expected %d\n",
           (long long)N_VGetLength(x), N);
    fails++;
  }

  /* Fill vectors: x[i] = (i+1) + i*i, y[i] = 1/(i+1) - i*(i+1) */
  xd = NV_DATA_ZS(x);
  yd = NV_DATA_ZS(y);
  for (i = 0; i < N; i++)
  {
    xd[i] = (i + 1.0) + _Complex_I * (double)i;
    yd[i] = 1.0 / (i + 1.0) - _Complex_I * (i + 1.0);
  }

  /* Test N_VLinearSum: z = 2*x + 3*y */
  N_VLinearSum(2.0, x, 3.0, y, z);
  zd = NV_DATA_ZS(z);
  for (i = 0; i < N; i++)
  {
    suncomplextype expected = 2.0 * xd[i] + 3.0 * yd[i];
    sunrealtype err = cabs(zd[i] - expected);
    if (err > TOL)
    {
      printf("FAIL N_VLinearSum[%lld]: err=%.2e\n", (long long)i, err);
      fails++;
    }
  }

  /* Test N_VScale */
  N_VScale(-0.5, x, z);
  zd = NV_DATA_ZS(z);
  for (i = 0; i < N; i++)
  {
    suncomplextype expected = -0.5 * xd[i];
    sunrealtype err = cabs(zd[i] - expected);
    if (err > TOL)
    {
      printf("FAIL N_VScale[%lld]: err=%.2e\n", (long long)i, err);
      fails++;
    }
  }

  /* Test N_VDotProd: Re(sum(conj(x)*y)) */
  {
    sunrealtype dp = N_VDotProd(x, y);
    suncomplextype sum = 0.0;
    for (i = 0; i < N; i++) { sum += conj(xd[i]) * yd[i]; }
    fails += check_val("N_VDotProd", dp, creal(sum));
  }

  /* Test N_VMaxNorm */
  {
    sunrealtype mn = N_VMaxNorm(x);
    sunrealtype expected = 0.0;
    for (i = 0; i < N; i++)
    {
      sunrealtype a = cabs(xd[i]);
      if (a > expected) { expected = a; }
    }
    fails += check_val("N_VMaxNorm", mn, expected);
  }

  /* Test N_VL1Norm */
  {
    sunrealtype l1 = N_VL1Norm(x);
    sunrealtype expected = 0.0;
    for (i = 0; i < N; i++) { expected += cabs(xd[i]); }
    fails += check_val("N_VL1Norm", l1, expected);
  }

  /* Test N_VProd (element-wise complex multiply) */
  N_VProd(x, y, z);
  zd = NV_DATA_ZS(z);
  for (i = 0; i < N; i++)
  {
    suncomplextype expected = xd[i] * yd[i];
    sunrealtype err = cabs(zd[i] - expected);
    if (err > TOL)
    {
      printf("FAIL N_VProd[%lld]: err=%.2e\n", (long long)i, err);
      fails++;
    }
  }

  /* Test N_VComplexDotProd (full complex result) */
  {
    suncomplextype cdp = N_VComplexDotProd_ComplexSerial(x, y);
    suncomplextype expected = 0.0;
    for (i = 0; i < N; i++) { expected += conj(xd[i]) * yd[i]; }
    sunrealtype err = cabs(cdp - expected);
    if (err > TOL)
    {
      printf("FAIL N_VComplexDotProd: err=%.2e\n", err);
      fails++;
    }
  }

  /* Test N_VComplexScale */
  {
    suncomplextype c = 2.0 + 3.0 * _Complex_I;
    N_VComplexScale_ComplexSerial(c, x, z);
    zd = NV_DATA_ZS(z);
    for (i = 0; i < N; i++)
    {
      suncomplextype expected = c * xd[i];
      sunrealtype err = cabs(zd[i] - expected);
      if (err > TOL)
      {
        printf("FAIL N_VComplexScale[%lld]: err=%.2e\n", (long long)i, err);
        fails++;
      }
    }
  }

  /* Test Pack/Unpack */
  {
    N_Vector re = N_VNew_Serial(N, sunctx);
    N_Vector im = N_VNew_Serial(N, sunctx);
    N_Vector re2 = N_VNew_Serial(N, sunctx);
    N_Vector im2 = N_VNew_Serial(N, sunctx);
    sunrealtype *red = N_VGetArrayPointer(re);
    sunrealtype *imd = N_VGetArrayPointer(im);
    sunrealtype omega = 2.5;

    for (i = 0; i < N; i++)
    {
      red[i] = (double)(i + 1);
      imd[i] = (double)(i * 2 + 1);
    }

    /* Pack: z[i] = re[i] + I * omega * im[i] */
    N_VComplexPack_ComplexSerial(re, omega, im, z);
    zd = NV_DATA_ZS(z);
    for (i = 0; i < N; i++)
    {
      suncomplextype expected = red[i] + _Complex_I * (omega * imd[i]);
      sunrealtype err = cabs(zd[i] - expected);
      if (err > TOL)
      {
        printf("FAIL Pack[%lld]: err=%.2e\n", (long long)i, err);
        fails++;
      }
    }

    /* Unpack: re2[i] = creal(z[i]), im2[i] = cimag(z[i]) / omega */
    N_VComplexUnpack_ComplexSerial(z, re2, 1.0 / omega, im2);
    sunrealtype *re2d = N_VGetArrayPointer(re2);
    sunrealtype *im2d = N_VGetArrayPointer(im2);
    for (i = 0; i < N; i++)
    {
      fails += check_val("Unpack_re", re2d[i], red[i]);
      fails += check_val("Unpack_im", im2d[i], imd[i]);
    }

    N_VDestroy(re);
    N_VDestroy(im);
    N_VDestroy(re2);
    N_VDestroy(im2);
  }

  /* Test Clone */
  {
    N_Vector w = N_VClone(x);
    if (w == NULL)
    {
      printf("FAIL: N_VClone returned NULL\n");
      fails++;
    }
    else
    {
      /* Clone should have same length but independent data */
      if (N_VGetLength(w) != N)
      {
        printf("FAIL: Clone length mismatch\n");
        fails++;
      }
      N_VDestroy(w);
    }
  }

  /* Test WrmsNorm */
  {
    /* Set weight vector: all 1.0 */
    suncomplextype *wd = NV_DATA_ZS(y);
    for (i = 0; i < N; i++) { wd[i] = 1.0 + 0.0 * _Complex_I; }

    sunrealtype wrms = N_VWrmsNorm(x, y);
    sunrealtype expected = 0.0;
    for (i = 0; i < N; i++) { expected += cabs(xd[i]) * cabs(xd[i]); }
    expected = sqrt(expected / N);
    fails += check_val("N_VWrmsNorm", wrms, expected);
  }

  N_VDestroy(x);
  N_VDestroy(y);
  N_VDestroy(z);
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
