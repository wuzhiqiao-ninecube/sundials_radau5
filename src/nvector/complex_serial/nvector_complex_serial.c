/* -----------------------------------------------------------------
 * Complex Serial N_Vector implementation for SUNDIALS.
 *
 * Uses C99 double _Complex as the base element type.
 * Follows SUNDIALS N_Vector interface conventions.
 * -----------------------------------------------------------------*/

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "nvector_complex_serial.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/* -----------------------------------------------------------------
 * Private helper prototypes
 * -----------------------------------------------------------------*/

static void VCopy_ComplexSerial(N_Vector x, N_Vector z);
static void VScaleBy_ComplexSerial(sunrealtype c, N_Vector x);
static void VNeg_ComplexSerial(N_Vector x, N_Vector z);

/* ================================================================
 * Constructors / Destructors
 * ================================================================*/

N_Vector N_VNewEmpty_ComplexSerial(sunindextype length, SUNContext sunctx)
{
  N_Vector v;
  N_VectorContent_ComplexSerial content;

  if (length < 0) { return NULL; }

  /* Create the generic N_Vector shell */
  v = N_VNewEmpty(sunctx);
  if (v == NULL) { return NULL; }

  /* Attach ops */
  v->ops->nvgetvectorid  = N_VGetVectorID_ComplexSerial;
  v->ops->nvclone        = N_VClone_ComplexSerial;
  v->ops->nvcloneempty   = N_VCloneEmpty_ComplexSerial;
  v->ops->nvdestroy      = N_VDestroy_ComplexSerial;
  v->ops->nvspace        = N_VSpace_ComplexSerial;
  v->ops->nvgetlength    = N_VGetLength_ComplexSerial;

  /* nvgetarraypointer returns NULL — use N_VGetComplexArrayPointer instead */
  v->ops->nvgetarraypointer = NULL;
  v->ops->nvsetarraypointer = NULL;

  /* Standard vector operations */
  v->ops->nvlinearsum    = N_VLinearSum_ComplexSerial;
  v->ops->nvconst        = N_VConst_ComplexSerial;
  v->ops->nvprod         = N_VProd_ComplexSerial;
  v->ops->nvdiv          = N_VDiv_ComplexSerial;
  v->ops->nvscale        = N_VScale_ComplexSerial;
  v->ops->nvabs          = N_VAbs_ComplexSerial;
  v->ops->nvinv          = N_VInv_ComplexSerial;
  v->ops->nvaddconst     = N_VAddConst_ComplexSerial;
  v->ops->nvdotprod      = N_VDotProd_ComplexSerial;
  v->ops->nvmaxnorm      = N_VMaxNorm_ComplexSerial;
  v->ops->nvwrmsnorm     = N_VWrmsNorm_ComplexSerial;
  v->ops->nvmin          = N_VMin_ComplexSerial;
  v->ops->nvl1norm       = N_VL1Norm_ComplexSerial;
  v->ops->nvcompare      = N_VCompare_ComplexSerial;
  v->ops->nvinvtest      = N_VInvTest_ComplexSerial;
  v->ops->nvminquotient  = N_VMinQuotient_ComplexSerial;

  /* Fused operations */
  v->ops->nvlinearcombination = N_VLinearCombination_ComplexSerial;

  /* Local reductions */
  v->ops->nvwsqrsumlocal = N_VWSqrSumLocal_ComplexSerial;

  /* Print */
  v->ops->nvprint     = N_VPrint_ComplexSerial;
  v->ops->nvprintfile = N_VPrintFile_ComplexSerial;

  /* Allocate content */
  content = (N_VectorContent_ComplexSerial)malloc(sizeof *content);
  if (content == NULL) { N_VDestroy(v); return NULL; }

  content->length   = length;
  content->own_data = SUNFALSE;
  content->data     = NULL;

  v->content = content;

  return v;
}

N_Vector N_VNew_ComplexSerial(sunindextype length, SUNContext sunctx)
{
  N_Vector v;
  suncomplextype *data;

  v = N_VNewEmpty_ComplexSerial(length, sunctx);
  if (v == NULL) { return NULL; }

  if (length > 0)
  {
    data = (suncomplextype *)calloc((size_t)length, sizeof(suncomplextype));
    if (data == NULL) { N_VDestroy(v); return NULL; }
    NV_OWN_DATA_ZS(v) = SUNTRUE;
    NV_DATA_ZS(v)     = data;
  }

  return v;
}

N_Vector N_VMake_ComplexSerial(sunindextype length, suncomplextype *data,
                               SUNContext sunctx)
{
  N_Vector v;

  v = N_VNewEmpty_ComplexSerial(length, sunctx);
  if (v == NULL) { return NULL; }

  NV_OWN_DATA_ZS(v) = SUNFALSE;
  NV_DATA_ZS(v)     = data;

  return v;
}

N_Vector N_VCloneEmpty_ComplexSerial(N_Vector w)
{
  N_Vector v;
  N_VectorContent_ComplexSerial content;

  if (w == NULL) { return NULL; }

  SUNContext sunctx = w->sunctx;

  v = N_VNewEmpty(sunctx);
  if (v == NULL) { return NULL; }

  /* Copy ops from w */
  if (N_VCopyOps(w, v)) { N_VDestroy(v); return NULL; }

  /* Allocate content */
  content = (N_VectorContent_ComplexSerial)malloc(sizeof *content);
  if (content == NULL) { N_VDestroy(v); return NULL; }

  content->length   = NV_LENGTH_ZS(w);
  content->own_data = SUNFALSE;
  content->data     = NULL;

  v->content = content;

  return v;
}

N_Vector N_VClone_ComplexSerial(N_Vector w)
{
  N_Vector v;
  sunindextype length;

  v = N_VCloneEmpty_ComplexSerial(w);
  if (v == NULL) { return NULL; }

  length = NV_LENGTH_ZS(w);
  if (length > 0)
  {
    suncomplextype *data =
        (suncomplextype *)calloc((size_t)length, sizeof(suncomplextype));
    if (data == NULL) { N_VDestroy(v); return NULL; }
    NV_OWN_DATA_ZS(v) = SUNTRUE;
    NV_DATA_ZS(v)     = data;
  }

  return v;
}

void N_VDestroy_ComplexSerial(N_Vector v)
{
  if (v == NULL) { return; }

  if (v->content != NULL)
  {
    N_VectorContent_ComplexSerial content = NV_CONTENT_ZS(v);
    if (content->own_data && content->data != NULL)
    {
      free(content->data);
      content->data = NULL;
    }
    free(content);
    v->content = NULL;
  }

  /* Free ops and vector shell */
  if (v->ops != NULL) { free(v->ops); v->ops = NULL; }
  free(v);
}

/* ================================================================
 * Accessor functions
 * ================================================================*/

sunindextype N_VGetLength_ComplexSerial(N_Vector v)
{
  return NV_LENGTH_ZS(v);
}

suncomplextype *N_VGetComplexArrayPointer_ComplexSerial(N_Vector v)
{
  return NV_DATA_ZS(v);
}

void N_VSetComplexArrayPointer_ComplexSerial(suncomplextype *data, N_Vector v)
{
  NV_DATA_ZS(v) = data;
}

N_Vector_ID N_VGetVectorID_ComplexSerial(N_Vector v)
{
  (void)v;
  return SUNDIALS_NVEC_CUSTOM;
}

void N_VPrint_ComplexSerial(N_Vector v)
{
  N_VPrintFile_ComplexSerial(v, stdout);
}

void N_VPrintFile_ComplexSerial(N_Vector v, FILE *outfile)
{
  sunindextype i, n;
  suncomplextype *vd;

  n  = NV_LENGTH_ZS(v);
  vd = NV_DATA_ZS(v);

  for (i = 0; i < n; i++)
  {
    fprintf(outfile, "  (%22.15e, %22.15e)\n", creal(vd[i]), cimag(vd[i]));
  }
  fprintf(outfile, "\n");
}

void N_VSpace_ComplexSerial(N_Vector v, sunindextype *lrw, sunindextype *liw)
{
  /* Each complex element uses 2 reals */
  *lrw = 2 * NV_LENGTH_ZS(v);
  *liw = 1;
}

/* ================================================================
 * Standard vector operations
 * ================================================================*/

/* z = a*x + b*y (a, b are real scalars) */
void N_VLinearSum_ComplexSerial(sunrealtype a, N_Vector x, sunrealtype b,
                                N_Vector y, N_Vector z)
{
  sunindextype i, n;
  suncomplextype *xd, *yd, *zd;

  n  = NV_LENGTH_ZS(x);
  xd = NV_DATA_ZS(x);
  yd = NV_DATA_ZS(y);
  zd = NV_DATA_ZS(z);

  /* Optimize common cases */
  if (a == ONE && b == ONE)
  {
    for (i = 0; i < n; i++) { zd[i] = xd[i] + yd[i]; }
    return;
  }
  if (a == ONE && b == -ONE)
  {
    for (i = 0; i < n; i++) { zd[i] = xd[i] - yd[i]; }
    return;
  }
  if (a == -ONE && b == ONE)
  {
    for (i = 0; i < n; i++) { zd[i] = yd[i] - xd[i]; }
    return;
  }
  if (a == ONE)
  {
    for (i = 0; i < n; i++) { zd[i] = xd[i] + b * yd[i]; }
    return;
  }
  if (b == ONE)
  {
    for (i = 0; i < n; i++) { zd[i] = a * xd[i] + yd[i]; }
    return;
  }

  /* General case */
  for (i = 0; i < n; i++) { zd[i] = a * xd[i] + b * yd[i]; }
}

/* z[i] = c + 0i for all i */
void N_VConst_ComplexSerial(sunrealtype c, N_Vector z)
{
  sunindextype i, n;
  suncomplextype *zd;

  n  = NV_LENGTH_ZS(z);
  zd = NV_DATA_ZS(z);

  for (i = 0; i < n; i++) { zd[i] = c + 0.0 * _Complex_I; }
}

/* z[i] = x[i] * y[i] (element-wise complex multiply) */
void N_VProd_ComplexSerial(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, n;
  suncomplextype *xd, *yd, *zd;

  n  = NV_LENGTH_ZS(x);
  xd = NV_DATA_ZS(x);
  yd = NV_DATA_ZS(y);
  zd = NV_DATA_ZS(z);

  for (i = 0; i < n; i++) { zd[i] = xd[i] * yd[i]; }
}

/* z[i] = x[i] / y[i] */
void N_VDiv_ComplexSerial(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, n;
  suncomplextype *xd, *yd, *zd;

  n  = NV_LENGTH_ZS(x);
  xd = NV_DATA_ZS(x);
  yd = NV_DATA_ZS(y);
  zd = NV_DATA_ZS(z);

  for (i = 0; i < n; i++) { zd[i] = xd[i] / yd[i]; }
}

/* z[i] = c * x[i] (real scalar) */
void N_VScale_ComplexSerial(sunrealtype c, N_Vector x, N_Vector z)
{
  sunindextype i, n;
  suncomplextype *xd, *zd;

  n  = NV_LENGTH_ZS(x);
  xd = NV_DATA_ZS(x);
  zd = NV_DATA_ZS(z);

  if (z == x)
  {
    VScaleBy_ComplexSerial(c, x);
    return;
  }
  if (c == ONE)  { VCopy_ComplexSerial(x, z); return; }
  if (c == -ONE) { VNeg_ComplexSerial(x, z); return; }

  for (i = 0; i < n; i++) { zd[i] = c * xd[i]; }
}

/* z[i] = |x[i]| (modulus, stored as real in complex slot) */
void N_VAbs_ComplexSerial(N_Vector x, N_Vector z)
{
  sunindextype i, n;
  suncomplextype *xd, *zd;

  n  = NV_LENGTH_ZS(x);
  xd = NV_DATA_ZS(x);
  zd = NV_DATA_ZS(z);

  for (i = 0; i < n; i++) { zd[i] = cabs(xd[i]); }
}

/* z[i] = 1.0 / x[i] */
void N_VInv_ComplexSerial(N_Vector x, N_Vector z)
{
  sunindextype i, n;
  suncomplextype *xd, *zd;

  n  = NV_LENGTH_ZS(x);
  xd = NV_DATA_ZS(x);
  zd = NV_DATA_ZS(z);

  for (i = 0; i < n; i++) { zd[i] = 1.0 / xd[i]; }
}

/* z[i] = x[i] + b (b is real) */
void N_VAddConst_ComplexSerial(N_Vector x, sunrealtype b, N_Vector z)
{
  sunindextype i, n;
  suncomplextype *xd, *zd;

  n  = NV_LENGTH_ZS(x);
  xd = NV_DATA_ZS(x);
  zd = NV_DATA_ZS(z);

  for (i = 0; i < n; i++) { zd[i] = xd[i] + b; }
}

/* Returns Re(sum(conj(x[i]) * y[i])) */
sunrealtype N_VDotProd_ComplexSerial(N_Vector x, N_Vector y)
{
  sunindextype i, n;
  suncomplextype *xd, *yd;
  suncomplextype sum = 0.0;

  n  = NV_LENGTH_ZS(x);
  xd = NV_DATA_ZS(x);
  yd = NV_DATA_ZS(y);

  for (i = 0; i < n; i++) { sum += conj(xd[i]) * yd[i]; }

  return creal(sum);
}

/* max_i |x[i]| */
sunrealtype N_VMaxNorm_ComplexSerial(N_Vector x)
{
  sunindextype i, n;
  suncomplextype *xd;
  sunrealtype mx = ZERO;

  n  = NV_LENGTH_ZS(x);
  xd = NV_DATA_ZS(x);

  for (i = 0; i < n; i++)
  {
    sunrealtype av = cabs(xd[i]);
    if (av > mx) { mx = av; }
  }

  return mx;
}

/* sqrt(sum(|x[i]|^2 * creal(w[i])^2) / n) */
sunrealtype N_VWrmsNorm_ComplexSerial(N_Vector x, N_Vector w)
{
  sunrealtype sum;
  sunindextype n = NV_LENGTH_ZS(x);

  sum = N_VWSqrSumLocal_ComplexSerial(x, w);

  return sqrt(sum / (sunrealtype)n);
}

/* sum(|x[i]|^2 * creal(w[i])^2) — local kernel */
sunrealtype N_VWSqrSumLocal_ComplexSerial(N_Vector x, N_Vector w)
{
  sunindextype i, n;
  suncomplextype *xd, *wd;
  sunrealtype sum = ZERO;

  n  = NV_LENGTH_ZS(x);
  xd = NV_DATA_ZS(x);
  wd = NV_DATA_ZS(w);

  for (i = 0; i < n; i++)
  {
    sunrealtype ax = cabs(xd[i]);
    sunrealtype wi = creal(wd[i]);  /* weight is real part */
    sum += ax * ax * wi * wi;
  }

  return sum;
}

/* min_i creal(x[i]) — only meaningful if vector is real-valued */
sunrealtype N_VMin_ComplexSerial(N_Vector x)
{
  sunindextype i, n;
  suncomplextype *xd;
  sunrealtype mn;

  n  = NV_LENGTH_ZS(x);
  xd = NV_DATA_ZS(x);

  mn = creal(xd[0]);
  for (i = 1; i < n; i++)
  {
    sunrealtype val = creal(xd[i]);
    if (val < mn) { mn = val; }
  }

  return mn;
}

/* sum_i |x[i]| */
sunrealtype N_VL1Norm_ComplexSerial(N_Vector x)
{
  sunindextype i, n;
  suncomplextype *xd;
  sunrealtype sum = ZERO;

  n  = NV_LENGTH_ZS(x);
  xd = NV_DATA_ZS(x);

  for (i = 0; i < n; i++) { sum += cabs(xd[i]); }

  return sum;
}

/* z[i] = (|x[i]| >= c) ? 1.0 : 0.0 */
void N_VCompare_ComplexSerial(sunrealtype c, N_Vector x, N_Vector z)
{
  sunindextype i, n;
  suncomplextype *xd, *zd;

  n  = NV_LENGTH_ZS(x);
  xd = NV_DATA_ZS(x);
  zd = NV_DATA_ZS(z);

  for (i = 0; i < n; i++)
  {
    zd[i] = (cabs(xd[i]) >= c) ? 1.0 : 0.0;
  }
}

/* z[i] = 1/x[i], returns SUNFALSE if any x[i] == 0 */
sunbooleantype N_VInvTest_ComplexSerial(N_Vector x, N_Vector z)
{
  sunindextype i, n;
  suncomplextype *xd, *zd;

  n  = NV_LENGTH_ZS(x);
  xd = NV_DATA_ZS(x);
  zd = NV_DATA_ZS(z);

  for (i = 0; i < n; i++)
  {
    if (cabs(xd[i]) == ZERO) { return SUNFALSE; }
    zd[i] = 1.0 / xd[i];
  }

  return SUNTRUE;
}

/* min_i (creal(num[i]) / creal(denom[i])) where creal(denom[i]) != 0 */
sunrealtype N_VMinQuotient_ComplexSerial(N_Vector num, N_Vector denom)
{
  sunindextype i, n;
  suncomplextype *nd, *dd;
  sunrealtype mn = SUN_BIG_REAL;

  n  = NV_LENGTH_ZS(num);
  nd = NV_DATA_ZS(num);
  dd = NV_DATA_ZS(denom);

  for (i = 0; i < n; i++)
  {
    sunrealtype dr = creal(dd[i]);
    if (dr != ZERO)
    {
      sunrealtype q = creal(nd[i]) / dr;
      if (q < mn) { mn = q; }
    }
  }

  return mn;
}

/* z = c[0]*V[0] + c[1]*V[1] + ... + c[nvec-1]*V[nvec-1] */
SUNErrCode N_VLinearCombination_ComplexSerial(int nvec, sunrealtype *c,
                                              N_Vector *V, N_Vector z)
{
  sunindextype i, n;
  int j;
  suncomplextype *zd;

  n  = NV_LENGTH_ZS(z);
  zd = NV_DATA_ZS(z);

  /* Zero z first */
  memset(zd, 0, (size_t)n * sizeof(suncomplextype));

  for (j = 0; j < nvec; j++)
  {
    suncomplextype *vd = NV_DATA_ZS(V[j]);
    sunrealtype cj = c[j];
    for (i = 0; i < n; i++) { zd[i] += cj * vd[i]; }
  }

  return SUN_SUCCESS;
}

/* ================================================================
 * Complex-specific operations
 * ================================================================*/

/* Full complex dot product: sum(conj(x[i]) * y[i]) */
suncomplextype N_VComplexDotProd_ComplexSerial(N_Vector x, N_Vector y)
{
  sunindextype i, n;
  suncomplextype *xd, *yd;
  suncomplextype sum = 0.0;

  n  = NV_LENGTH_ZS(x);
  xd = NV_DATA_ZS(x);
  yd = NV_DATA_ZS(y);

  for (i = 0; i < n; i++) { sum += conj(xd[i]) * yd[i]; }

  return sum;
}

/* z[i] = c * x[i] (complex scalar) */
void N_VComplexScale_ComplexSerial(suncomplextype c, N_Vector x, N_Vector z)
{
  sunindextype i, n;
  suncomplextype *xd, *zd;

  n  = NV_LENGTH_ZS(x);
  xd = NV_DATA_ZS(x);
  zd = NV_DATA_ZS(z);

  for (i = 0; i < n; i++) { zd[i] = c * xd[i]; }
}

/* z[i] = a*x[i] + b*y[i] (complex scalars) */
void N_VComplexLinearSum_ComplexSerial(suncomplextype a, N_Vector x,
                                       suncomplextype b, N_Vector y,
                                       N_Vector z)
{
  sunindextype i, n;
  suncomplextype *xd, *yd, *zd;

  n  = NV_LENGTH_ZS(x);
  xd = NV_DATA_ZS(x);
  yd = NV_DATA_ZS(y);
  zd = NV_DATA_ZS(z);

  for (i = 0; i < n; i++) { zd[i] = a * xd[i] + b * yd[i]; }
}

/* z[i] = c (complex constant) */
void N_VComplexConst_ComplexSerial(suncomplextype c, N_Vector z)
{
  sunindextype i, n;
  suncomplextype *zd;

  n  = NV_LENGTH_ZS(z);
  zd = NV_DATA_ZS(z);

  for (i = 0; i < n; i++) { zd[i] = c; }
}

/* Pack: z[i] = re[i] + I * scale_im * im[i]
 * re, im are real N_Vectors (any implementation with GetArrayPointer) */
void N_VComplexPack_ComplexSerial(N_Vector re, sunrealtype scale_im,
                                   N_Vector im, N_Vector z)
{
  sunindextype i, n;
  sunrealtype *red, *imd;
  suncomplextype *zd;

  n   = NV_LENGTH_ZS(z);
  red = N_VGetArrayPointer(re);
  imd = N_VGetArrayPointer(im);
  zd  = NV_DATA_ZS(z);

  for (i = 0; i < n; i++)
  {
    zd[i] = red[i] + _Complex_I * (scale_im * imd[i]);
  }
}

/* Unpack: re[i] = creal(z[i]),  im[i] = cimag(z[i]) * scale_im */
void N_VComplexUnpack_ComplexSerial(N_Vector z, N_Vector re,
                                     sunrealtype scale_im, N_Vector im)
{
  sunindextype i, n;
  sunrealtype *red, *imd;
  suncomplextype *zd;

  n   = NV_LENGTH_ZS(z);
  zd  = NV_DATA_ZS(z);
  red = N_VGetArrayPointer(re);
  imd = N_VGetArrayPointer(im);

  for (i = 0; i < n; i++)
  {
    red[i] = creal(zd[i]);
    imd[i] = cimag(zd[i]) * scale_im;
  }
}

/* ================================================================
 * Private helper functions
 * ================================================================*/

static void VCopy_ComplexSerial(N_Vector x, N_Vector z)
{
  sunindextype n = NV_LENGTH_ZS(x);
  memcpy(NV_DATA_ZS(z), NV_DATA_ZS(x), (size_t)n * sizeof(suncomplextype));
}

static void VScaleBy_ComplexSerial(sunrealtype c, N_Vector x)
{
  sunindextype i, n;
  suncomplextype *xd;

  n  = NV_LENGTH_ZS(x);
  xd = NV_DATA_ZS(x);

  for (i = 0; i < n; i++) { xd[i] *= c; }
}

static void VNeg_ComplexSerial(N_Vector x, N_Vector z)
{
  sunindextype i, n;
  suncomplextype *xd, *zd;

  n  = NV_LENGTH_ZS(x);
  xd = NV_DATA_ZS(x);
  zd = NV_DATA_ZS(z);

  for (i = 0; i < n; i++) { zd[i] = -xd[i]; }
}
