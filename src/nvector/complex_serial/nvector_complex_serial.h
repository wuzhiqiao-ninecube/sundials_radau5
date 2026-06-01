/* -----------------------------------------------------------------
 * Complex Serial N_Vector implementation for SUNDIALS.
 *
 * Uses C99 double _Complex as the base element type.
 * Follows SUNDIALS N_Vector interface conventions.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_COMPLEX_SERIAL_H
#define _NVECTOR_COMPLEX_SERIAL_H

#include <complex.h>
#include <stdio.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Complex type matched to SUNDIALS precision configuration */
#if defined(SUNDIALS_DOUBLE_PRECISION)
typedef double _Complex suncomplextype;
#elif defined(SUNDIALS_SINGLE_PRECISION)
typedef float _Complex suncomplextype;
#elif defined(SUNDIALS_EXTENDED_PRECISION)
typedef long double _Complex suncomplextype;
#elif defined(SUNDIALS_FLOAT128_PRECISION)
typedef __float128 _Complex suncomplextype;
#else
typedef double _Complex suncomplextype;
#endif

/* -----------------------------------------------------------------
 * Content structure
 * -----------------------------------------------------------------*/

struct _N_VectorContent_ComplexSerial
{
  sunindextype length;       /* number of complex elements */
  sunbooleantype own_data;   /* data ownership flag */
  suncomplextype *data;      /* complex data array */
};

typedef struct _N_VectorContent_ComplexSerial *N_VectorContent_ComplexSerial;

/* -----------------------------------------------------------------
 * Access macros
 * -----------------------------------------------------------------*/

#define NV_CONTENT_ZS(v)   ((N_VectorContent_ComplexSerial)(v->content))
#define NV_LENGTH_ZS(v)    (NV_CONTENT_ZS(v)->length)
#define NV_OWN_DATA_ZS(v)  (NV_CONTENT_ZS(v)->own_data)
#define NV_DATA_ZS(v)      (NV_CONTENT_ZS(v)->data)
#define NV_Ith_ZS(v, i)    (NV_DATA_ZS(v)[i])

/* -----------------------------------------------------------------
 * Constructors / destructors
 * -----------------------------------------------------------------*/

N_Vector N_VNewEmpty_ComplexSerial(sunindextype length, SUNContext sunctx);
N_Vector N_VNew_ComplexSerial(sunindextype length, SUNContext sunctx);
N_Vector N_VMake_ComplexSerial(sunindextype length, suncomplextype *data,
                               SUNContext sunctx);

N_Vector N_VCloneEmpty_ComplexSerial(N_Vector w);
N_Vector N_VClone_ComplexSerial(N_Vector w);
void N_VDestroy_ComplexSerial(N_Vector v);

/* -----------------------------------------------------------------
 * Accessor functions
 * -----------------------------------------------------------------*/

sunindextype N_VGetLength_ComplexSerial(N_Vector v);

/* Returns the complex data pointer (use this instead of N_VGetArrayPointer) */
suncomplextype *N_VGetComplexArrayPointer_ComplexSerial(N_Vector v);
void N_VSetComplexArrayPointer_ComplexSerial(suncomplextype *data, N_Vector v);

void N_VPrint_ComplexSerial(N_Vector v);
void N_VPrintFile_ComplexSerial(N_Vector v, FILE *outfile);

/* -----------------------------------------------------------------
 * Standard N_Vector ops (real scalar versions)
 * -----------------------------------------------------------------*/

N_Vector_ID N_VGetVectorID_ComplexSerial(N_Vector v);
void N_VLinearSum_ComplexSerial(sunrealtype a, N_Vector x, sunrealtype b,
                                N_Vector y, N_Vector z);
void N_VConst_ComplexSerial(sunrealtype c, N_Vector z);
void N_VProd_ComplexSerial(N_Vector x, N_Vector y, N_Vector z);
void N_VDiv_ComplexSerial(N_Vector x, N_Vector y, N_Vector z);
void N_VScale_ComplexSerial(sunrealtype c, N_Vector x, N_Vector z);
void N_VAbs_ComplexSerial(N_Vector x, N_Vector z);
void N_VInv_ComplexSerial(N_Vector x, N_Vector z);
void N_VAddConst_ComplexSerial(N_Vector x, sunrealtype b, N_Vector z);
sunrealtype N_VDotProd_ComplexSerial(N_Vector x, N_Vector y);
sunrealtype N_VMaxNorm_ComplexSerial(N_Vector x);
sunrealtype N_VWrmsNorm_ComplexSerial(N_Vector x, N_Vector w);
sunrealtype N_VMin_ComplexSerial(N_Vector x);
sunrealtype N_VL1Norm_ComplexSerial(N_Vector x);
void N_VCompare_ComplexSerial(sunrealtype c, N_Vector x, N_Vector z);
sunbooleantype N_VInvTest_ComplexSerial(N_Vector x, N_Vector z);
sunrealtype N_VMinQuotient_ComplexSerial(N_Vector num, N_Vector denom);
void N_VSpace_ComplexSerial(N_Vector v, sunindextype *lrw, sunindextype *liw);

/* Fused operations */
SUNErrCode N_VLinearCombination_ComplexSerial(int nvec, sunrealtype *c,
                                              N_Vector *V, N_Vector z);

/* Local reduction (for compatibility) */
sunrealtype N_VWSqrSumLocal_ComplexSerial(N_Vector x, N_Vector w);

/* -----------------------------------------------------------------
 * Complex-specific operations (beyond standard ops table)
 * -----------------------------------------------------------------*/

/* Full complex dot product: sum(conj(x[i]) * y[i]) */
suncomplextype N_VComplexDotProd_ComplexSerial(N_Vector x, N_Vector y);

/* Complex scalar operations */
void N_VComplexScale_ComplexSerial(suncomplextype c, N_Vector x, N_Vector z);
void N_VComplexLinearSum_ComplexSerial(suncomplextype a, N_Vector x,
                                       suncomplextype b, N_Vector y,
                                       N_Vector z);
void N_VComplexConst_ComplexSerial(suncomplextype c, N_Vector z);

/* Pack/Unpack: convert between real vector pair and complex vector
 *   Pack:   z[i] = re[i] + I * scale_im * im[i]
 *   Unpack: re[i] = creal(z[i]),  im[i] = cimag(z[i]) * scale_im
 * re and im are standard real N_Vectors (e.g., N_Vector_Serial).
 * z is a ComplexSerial N_Vector.
 */
void N_VComplexPack_ComplexSerial(N_Vector re, sunrealtype scale_im,
                                   N_Vector im, N_Vector z);
void N_VComplexUnpack_ComplexSerial(N_Vector z, N_Vector re,
                                     sunrealtype scale_im, N_Vector im);

#ifdef __cplusplus
}
#endif

#endif /* _NVECTOR_COMPLEX_SERIAL_H */
