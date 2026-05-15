/* ---------------------------------------------------------------------------
 * RADAU5 — 3-stage, order-5 implicit Runge-Kutta solver (Radau IIA)
 * Rewritten in C using SUNDIALS abstractions for large-scale sparse support.
 *
 * Based on the Fortran code by E. Hairer and G. Wanner:
 *   "Solving Ordinary Differential Equations II"
 *   Springer Series in Computational Mathematics 14, 1996.
 * ---------------------------------------------------------------------------*/

#ifndef RADAU5_H_
#define RADAU5_H_

#include <sundials/sundials_context.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------------------------------------------------------
 * Return codes
 * ---------------------------------------------------------------------------*/
#define RADAU5_SUCCESS        0
#define RADAU5_TSTOP_RETURN   1
#define RADAU5_SOLOUT_RETURN  2
#define RADAU5_MEM_NULL      -1
#define RADAU5_ILL_INPUT     -2
#define RADAU5_TOO_MANY_STEPS -3
#define RADAU5_STEP_TOO_SMALL -4
#define RADAU5_SINGULAR_MATRIX -5
#define RADAU5_CONV_FAILURE  -6
#define RADAU5_NEWT_PREDICT  -60  /* Newton predicted slow convergence, h already reduced */
#define RADAU5_LSETUP_FAIL   -7
#define RADAU5_LSOLVE_FAIL   -8
#define RADAU5_RHSFUNC_FAIL  -9
#define RADAU5_RHSFUNC_RECVR -90  /* internal: recoverable RHS failure, triggers label78 */
#define RADAU5_MEM_FAIL      -10
#define RADAU5_IC_FAIL       -11
#define RADAU5_ROOTFN_FAIL   -12
#define RADAU5_ROOT_RETURN    3

/* ---------------------------------------------------------------------------
 * User-supplied function types
 * ---------------------------------------------------------------------------*/

/* RHS function: ydot = f(t, y)
 * Return 0 on success, > 0 for recoverable error (solver will reduce
 * step size and retry), < 0 for unrecoverable error (solver halts). */
typedef int (*Radau5RhsFn)(sunrealtype t, N_Vector y, N_Vector ydot,
                           void* user_data);

/* Jacobian function: Jac = df/dy at (t, y). fy = f(t,y) already evaluated. */
typedef int (*Radau5JacFn)(sunrealtype t, N_Vector y, N_Vector fy,
                           SUNMatrix Jac, void* user_data,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Mass matrix function: fills M. */
typedef int (*Radau5MassFn)(sunrealtype t, SUNMatrix M, void* user_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Solution output callback, called after each accepted step. */
typedef int (*Radau5SolOutFn)(long int nr, sunrealtype told, sunrealtype t,
                              N_Vector y, void* radau5_mem, void* user_data);

/* Root function for event detection: evaluates nrtfn scalar functions g_i(t,y).
 * Return 0 on success, < 0 for unrecoverable error. */
typedef int (*Radau5RootFn)(sunrealtype t, N_Vector y, sunrealtype* gout,
                            void* user_data);

/* ---------------------------------------------------------------------------
 * Create / Destroy
 * ---------------------------------------------------------------------------*/
void* Radau5Create(SUNContext sunctx);
void  Radau5Free(void** radau5_mem);

/* ---------------------------------------------------------------------------
 * Initialize
 * ---------------------------------------------------------------------------*/
int Radau5Init(void* radau5_mem, Radau5RhsFn rhs, sunrealtype t0,
               N_Vector y0);

/* ---------------------------------------------------------------------------
 * Linear solver setup — user provides J template and optional M template.
 * M is only needed when both J and M are sparse (for union sparsity pattern).
 * Pass NULL for M when J is dense/band or when there is no mass matrix.
 * ---------------------------------------------------------------------------*/
int Radau5SetLinearSolver(void* radau5_mem, SUNMatrix J, SUNMatrix M);

/* ---------------------------------------------------------------------------
 * Optional setters
 * ---------------------------------------------------------------------------*/
int Radau5SetJacFn(void* radau5_mem, Radau5JacFn jac);
int Radau5SetMassFn(void* radau5_mem, Radau5MassFn mas, SUNMatrix M);
int Radau5SetSolOutFn(void* radau5_mem, Radau5SolOutFn solout);
int Radau5SetUserData(void* radau5_mem, void* user_data);
int Radau5SetMaxNumSteps(void* radau5_mem, long int mxsteps);
int Radau5SetMaxNewtonIter(void* radau5_mem, int maxnit);
int Radau5SetInitStep(void* radau5_mem, sunrealtype h0);
int Radau5SetMaxStep(void* radau5_mem, sunrealtype hmax);
int Radau5SetSafetyFactor(void* radau5_mem, sunrealtype safe);
int Radau5SetStepSizeController(void* radau5_mem, int pred);
int Radau5SetDAEIndex(void* radau5_mem, sunindextype nind1,
                      sunindextype nind2, sunindextype nind3);
int Radau5SetStartNewton(void* radau5_mem, int startn);
int Radau5SetSchurDecomp(void* radau5_mem, int use_schur);
int Radau5ResetForDiscontinuity(void* radau5_mem, sunrealtype h0);
int Radau5SetSparsityPattern(void* radau5_mem, SUNMatrix S);
int Radau5SetNumStages(void* radau5_mem, int ns);
int Radau5SetOrderLimits(void* radau5_mem, int nsmin, int nsmax);
int Radau5SStolerances(void* radau5_mem, sunrealtype rtol, sunrealtype atol);
int Radau5SVtolerances(void* radau5_mem, N_Vector rtol, N_Vector atol);

/* ---------------------------------------------------------------------------
 * Main solver
 * ---------------------------------------------------------------------------*/
int Radau5Solve(void* radau5_mem, sunrealtype tout, N_Vector yout,
                sunrealtype* tret);

/* ---------------------------------------------------------------------------
 * Consistent initial conditions (index-1 DAEs only)
 * ---------------------------------------------------------------------------*/
int Radau5CalcIC(void* radau5_mem, N_Vector id);

/* ---------------------------------------------------------------------------
 * Continuous output — callable from SolOut callback
 * ---------------------------------------------------------------------------*/
sunrealtype Radau5Contr(void* radau5_mem, sunindextype i, sunrealtype t);

/* ---------------------------------------------------------------------------
 * Statistics getters
 * ---------------------------------------------------------------------------*/
int Radau5GetNumSteps(void* radau5_mem, long int* nsteps);
int Radau5GetNumRhsEvals(void* radau5_mem, long int* nfcn);
int Radau5GetNumJacEvals(void* radau5_mem, long int* njac);
int Radau5GetNumLinSolves(void* radau5_mem, long int* nsol);
int Radau5GetNumDecomps(void* radau5_mem, long int* ndec);
int Radau5GetNumAccSteps(void* radau5_mem, long int* naccpt);
int Radau5GetNumRejSteps(void* radau5_mem, long int* nrejct);
int Radau5GetNumNewtonIters(void* radau5_mem, long int* nnewt);
int Radau5GetCurrentStep(void* radau5_mem, sunrealtype* hcur);
int Radau5GetCurrentTime(void* radau5_mem, sunrealtype* tcur);

/* ---------------------------------------------------------------------------
 * Rootfinding (event detection)
 * ---------------------------------------------------------------------------*/
int Radau5RootInit(void* radau5_mem, int nrtfn, Radau5RootFn g);
int Radau5SetRootDirection(void* radau5_mem, int* rootdir);
int Radau5GetRootInfo(void* radau5_mem, int* rootsfound);
int Radau5GetNumGEvals(void* radau5_mem, long int* ngevals);

#ifdef __cplusplus
}
#endif

#endif /* RADAU5_H_ */
