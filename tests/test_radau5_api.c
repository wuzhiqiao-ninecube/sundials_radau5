/* ---------------------------------------------------------------------------
 * test_radau5_api.c — Public API validation
 *
 * Tests Create/Init/Free, parameter setters, and basic error handling.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

static int nfail = 0;

#define CHECK(cond, msg) do { \
  if (!(cond)) { fprintf(stderr, "FAIL: %s\n", msg); nfail++; } \
} while(0)

static int rhs_zero(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  N_VConst(SUN_RCONST(0.0), yd);
  return 0;
}

int main(void)
{
  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  /* --- Radau5Create --- */
  void* mem = Radau5Create(sunctx);
  CHECK(mem != NULL, "Radau5Create returned NULL");

  /* --- Radau5Init with NULL rhs should fail --- */
  N_Vector y0 = N_VNew_Serial(2, sunctx);
  N_VConst(SUN_RCONST(1.0), y0);

  int ret = Radau5Init(mem, NULL, 0.0, y0);
  CHECK(ret == RADAU5_ILL_INPUT, "Init with NULL rhs should return ILL_INPUT");

  /* --- Radau5Init with valid inputs --- */
  ret = Radau5Init(mem, rhs_zero, 0.0, y0);
  CHECK(ret == RADAU5_SUCCESS, "Init with valid inputs should succeed");

  /* --- Radau5SetLinearSolver with NULL J should fail --- */
  ret = Radau5SetLinearSolver(mem, NULL, NULL);
  CHECK(ret == RADAU5_ILL_INPUT, "SetLinearSolver with NULL J should fail");

  /* --- Radau5SetLinearSolver with valid dense J --- */
  SUNMatrix J = SUNDenseMatrix(2, 2, sunctx);
  ret = Radau5SetLinearSolver(mem, J, NULL);
  CHECK(ret == RADAU5_SUCCESS, "SetLinearSolver with dense J should succeed");

  /* --- Tolerance setters --- */
  ret = Radau5SStolerances(mem, 1.0e-6, 1.0e-8);
  CHECK(ret == RADAU5_SUCCESS, "SStolerances should succeed");

  ret = Radau5SStolerances(mem, -1.0, 1.0e-8);
  CHECK(ret == RADAU5_ILL_INPUT, "SStolerances with negative rtol should fail");

  /* --- Parameter setters --- */
  ret = Radau5SetMaxNumSteps(mem, 50000);
  CHECK(ret == RADAU5_SUCCESS, "SetMaxNumSteps should succeed");

  ret = Radau5SetMaxNumSteps(mem, -1);
  CHECK(ret == RADAU5_ILL_INPUT, "SetMaxNumSteps with negative should fail");

  ret = Radau5SetMaxNewtonIter(mem, 10);
  CHECK(ret == RADAU5_SUCCESS, "SetMaxNewtonIter should succeed");

  ret = Radau5SetInitStep(mem, 1.0e-4);
  CHECK(ret == RADAU5_SUCCESS, "SetInitStep should succeed");

  ret = Radau5SetMaxStep(mem, 0.1);
  CHECK(ret == RADAU5_SUCCESS, "SetMaxStep should succeed");

  ret = Radau5SetSafetyFactor(mem, 0.9);
  CHECK(ret == RADAU5_SUCCESS, "SetSafetyFactor should succeed");

  ret = Radau5SetSafetyFactor(mem, 1.5);
  CHECK(ret == RADAU5_ILL_INPUT, "SetSafetyFactor >= 1 should fail");

  ret = Radau5SetStepSizeController(mem, 1);
  CHECK(ret == RADAU5_SUCCESS, "SetStepSizeController(1) should succeed");

  ret = Radau5SetStepSizeController(mem, 3);
  CHECK(ret == RADAU5_ILL_INPUT, "SetStepSizeController(3) should fail");

  ret = Radau5SetDAEIndex(mem, 2, 0, 0);
  CHECK(ret == RADAU5_SUCCESS, "SetDAEIndex should succeed");

  /* --- Statistics getters should return zero after Init --- */
  long int val;
  Radau5GetNumSteps(mem, &val);
  CHECK(val == 0, "nsteps should be 0 after Init");

  Radau5GetNumRhsEvals(mem, &val);
  CHECK(val == 0, "nfcn should be 0 after Init");

  Radau5GetNumJacEvals(mem, &val);
  CHECK(val == 0, "njac should be 0 after Init");

  /* --- Radau5Solve on y'=0 from t=0 to t=1 --- */
  Radau5SStolerances(mem, 1.0e-8, 1.0e-10);
  Radau5SetInitStep(mem, 0.01);

  N_Vector yout = N_VNew_Serial(2, sunctx);
  sunrealtype tret;
  ret = Radau5Solve(mem, 1.0, yout, &tret);
  CHECK(ret == RADAU5_SUCCESS, "Solve y'=0 should succeed");
  CHECK(fabs(tret - 1.0) < 1.0e-12, "tret should be 1.0");

  sunrealtype* yd = N_VGetArrayPointer(yout);
  CHECK(fabs(yd[0] - 1.0) < 1.0e-10, "y[0] should remain 1.0 for y'=0");
  CHECK(fabs(yd[1] - 1.0) < 1.0e-10, "y[1] should remain 1.0 for y'=0");

  /* Statistics should be nonzero after solve */
  Radau5GetNumSteps(mem, &val);
  CHECK(val > 0, "nsteps should be > 0 after Solve");

  Radau5GetNumRhsEvals(mem, &val);
  CHECK(val > 0, "nfcn should be > 0 after Solve");

  /* --- Radau5Free sets pointer to NULL --- */
  Radau5Free(&mem);
  CHECK(mem == NULL, "Radau5Free should set pointer to NULL");

  /* --- Calling with NULL mem --- */
  ret = Radau5Init(NULL, rhs_zero, 0.0, y0);
  CHECK(ret == RADAU5_MEM_NULL, "Init with NULL mem should return MEM_NULL");

  /* Cleanup */
  N_VDestroy(y0);
  N_VDestroy(yout);
  SUNMatDestroy(J);
  SUNContext_Free(&sunctx);

  if (nfail == 0)
    printf("test_radau5_api: PASSED\n");
  else
    printf("test_radau5_api: %d FAILURES\n", nfail);

  return nfail;
}
