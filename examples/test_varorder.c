/* ---------------------------------------------------------------------------
 * test_varorder.c — Test variable-order mode and fixed ns=5/7
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"
#include "radau5_impl.h"

#define EPS 1.0e-6

static int vdpol_rhs(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)t; (void)ud;
  sunrealtype* yv  = N_VGetArrayPointer(y);
  sunrealtype* ydv = N_VGetArrayPointer(yd);
  ydv[0] = yv[1];
  ydv[1] = ((1.0 - yv[0] * yv[0]) * yv[1] - yv[0]) / EPS;
  return 0;
}

static int vdpol_jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                     void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;
  sunrealtype* yv = N_VGetArrayPointer(y);
  SM_ELEMENT_D(J, 0, 0) = 0.0;
  SM_ELEMENT_D(J, 0, 1) = 1.0;
  SM_ELEMENT_D(J, 1, 0) = (-2.0 * yv[0] * yv[1] - 1.0) / EPS;
  SM_ELEMENT_D(J, 1, 1) = (1.0 - yv[0] * yv[0]) / EPS;
  return 0;
}

static int rober_rhs(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)t; (void)ud;
  sunrealtype* yv = N_VGetArrayPointer(y);
  sunrealtype* ydv = N_VGetArrayPointer(yd);
  ydv[0] = -0.04*yv[0] + 1e4*yv[1]*yv[2];
  ydv[1] =  0.04*yv[0] - 1e4*yv[1]*yv[2] - 3e7*yv[1]*yv[1];
  ydv[2] =  3e7*yv[1]*yv[1];
  return 0;
}

static int rober_jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                     void* ud, N_Vector t1, N_Vector t2, N_Vector t3)
{
  (void)t; (void)fy; (void)ud; (void)t1; (void)t2; (void)t3;
  sunrealtype* yv = N_VGetArrayPointer(y);
  SM_ELEMENT_D(J,0,0) = -0.04;
  SM_ELEMENT_D(J,0,1) = 1e4*yv[2];
  SM_ELEMENT_D(J,0,2) = 1e4*yv[1];
  SM_ELEMENT_D(J,1,0) = 0.04;
  SM_ELEMENT_D(J,1,1) = -1e4*yv[2] - 6e7*yv[1];
  SM_ELEMENT_D(J,1,2) = -1e4*yv[1];
  SM_ELEMENT_D(J,2,0) = 0.0;
  SM_ELEMENT_D(J,2,1) = 6e7*yv[1];
  SM_ELEMENT_D(J,2,2) = 0.0;
  return 0;
}

static int test_vdpol(SUNContext sunctx, int ns, int use_schur, const char* label)
{
  void* mem = Radau5Create(sunctx);
  Radau5SetNumStages(mem, ns);
  N_Vector y0 = N_VNew_Serial(2, sunctx);
  N_VGetArrayPointer(y0)[0] = 2.0;
  N_VGetArrayPointer(y0)[1] = 0.0;
  Radau5Init(mem, vdpol_rhs, 0.0, y0);
  SUNMatrix Jt = SUNDenseMatrix(2, 2, sunctx);
  Radau5SetLinearSolver(mem, Jt, NULL);
  Radau5SetSchurDecomp(mem, use_schur);
  Radau5SetJacFn(mem, vdpol_jac);
  Radau5SStolerances(mem, 1e-6, 1e-6);
  Radau5SetInitStep(mem, 1e-6);
  N_Vector yout = N_VNew_Serial(2, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 2.0, yout, &tret);
  sunrealtype yref[2] = {1.70616773217048e+00, -8.92809701024798e-01};
  sunrealtype* yd = N_VGetArrayPointer(yout);
  sunrealtype err0 = fabs(yd[0] - yref[0]);
  sunrealtype err1 = fabs(yd[1] - yref[1]);
  long int naccpt;
  Radau5GetNumAccSteps(mem, &naccpt);
  int pass = (ret == 0 && err0 < 1e-4 && err1 < 1e-4);
  printf("  %s: ret=%d err=[%.2e,%.2e] naccpt=%ld %s\n",
         label, ret, err0, err1, naccpt, pass ? "PASS" : "FAIL");
  N_VDestroy(y0); N_VDestroy(yout); SUNMatDestroy(Jt);
  Radau5Free(&mem);
  return pass;
}

static int test_rober_varorder(SUNContext sunctx)
{
  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, 3, 7);
  Radau5Mem rmem = RADAU5_MEM(mem);
  rmem->vitu = 0.02;

  N_Vector y0 = N_VNew_Serial(3, sunctx);
  N_VGetArrayPointer(y0)[0] = 1.0;
  N_VGetArrayPointer(y0)[1] = 0.0;
  N_VGetArrayPointer(y0)[2] = 0.0;

  Radau5Init(mem, rober_rhs, 0.0, y0);
  SUNMatrix Jt = SUNDenseMatrix(3, 3, sunctx);
  Radau5SetLinearSolver(mem, Jt, NULL);
  Radau5SetJacFn(mem, rober_jac);
  Radau5SStolerances(mem, 1e-7, 1e-7);
  Radau5SetInitStep(mem, 1e-6);

  N_Vector yout = N_VNew_Serial(3, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 1e11, yout, &tret);

  sunrealtype* yd = N_VGetArrayPointer(yout);
  sunrealtype err2 = fabs(yd[2] - 0.9999999791665);

  long int naccpt;
  Radau5GetNumAccSteps(mem, &naccpt);

  int pass = (ret == 0 && err2 < 1e-5);
  printf("  Robertson varorder (3->7, vitu=0.02): ret=%d err=%.2e naccpt=%ld final_ns=%d %s\n",
         ret, err2, naccpt, rmem->ns, pass ? "PASS" : "FAIL");

  N_VDestroy(y0); N_VDestroy(yout); SUNMatDestroy(Jt);
  Radau5Free(&mem);
  return pass;
}

static int test_vdpol_varorder(SUNContext sunctx)
{
  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, 3, 7);

  N_Vector y0 = N_VNew_Serial(2, sunctx);
  N_VGetArrayPointer(y0)[0] = 2.0;
  N_VGetArrayPointer(y0)[1] = 0.0;

  Radau5Init(mem, vdpol_rhs, 0.0, y0);
  SUNMatrix Jt = SUNDenseMatrix(2, 2, sunctx);
  Radau5SetLinearSolver(mem, Jt, NULL);
  Radau5SetJacFn(mem, vdpol_jac);
  Radau5SStolerances(mem, 1e-6, 1e-6);
  Radau5SetInitStep(mem, 1e-6);

  N_Vector yout = N_VNew_Serial(2, sunctx);
  sunrealtype tret;
  int ret = Radau5Solve(mem, 2.0, yout, &tret);

  sunrealtype yref[2] = {1.70616773217048e+00, -8.92809701024798e-01};
  sunrealtype* yd = N_VGetArrayPointer(yout);
  sunrealtype err0 = fabs(yd[0] - yref[0]);
  sunrealtype err1 = fabs(yd[1] - yref[1]);

  int pass = (ret == 0 && err0 < 1e-4 && err1 < 1e-4);
  printf("  VdPol varorder (3->7, default vitu): ret=%d err=[%.2e,%.2e] %s\n",
         ret, err0, err1, pass ? "PASS" : "FAIL");

  N_VDestroy(y0); N_VDestroy(yout); SUNMatDestroy(Jt);
  Radau5Free(&mem);
  return pass;
}

int main(void)
{
  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  int all_pass = 1;
  printf("=== Variable-Order Tests ===\n");

  if (!test_vdpol(sunctx, 5, 0, "VdPol ns=5 eigen")) all_pass = 0;
  if (!test_vdpol(sunctx, 7, 0, "VdPol ns=7 eigen")) all_pass = 0;
  if (!test_vdpol(sunctx, 5, 1, "VdPol ns=5 schur")) all_pass = 0;
  if (!test_vdpol(sunctx, 7, 1, "VdPol ns=7 schur")) all_pass = 0;
  if (!test_rober_varorder(sunctx)) all_pass = 0;
  if (!test_vdpol_varorder(sunctx)) all_pass = 0;

  printf("\n%s\n", all_pass ? "ALL PASSED" : "SOME FAILED");

  SUNContext_Free(&sunctx);
  return all_pass ? 0 : 1;
}