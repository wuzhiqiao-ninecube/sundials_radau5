# RADAU5 User Guide and API Reference

## 1. Overview

RADAU5 is a C implementation of the 3-stage, order-5 implicit Runge-Kutta method
from the Radau IIA family, based on the Fortran code by E. Hairer and G. Wanner
(*Solving Ordinary Differential Equations II*, Springer, 1996). It is designed
for stiff ODEs and differential-algebraic equations (DAEs) of index up to 3.

The solver is built on SUNDIALS abstractions (`N_Vector`, `SUNMatrix`,
`SUNLinearSolver`, `SUNContext`) and supports dense, banded, and sparse (KLU)
linear algebra backends. Two transform modes are available for decoupling the
Newton system: the classical eigenvalue decomposition (default) and an optional
Schur decomposition using orthogonal transforms.

Key features:

- Implicit 3-stage Radau IIA method (A-stable, L-stable, stiffly accurate)
- Adaptive step size control (Gustafsson or classical)
- Simplified Newton iteration with Jacobian reuse
- Dense, band, and sparse (KLU) linear solvers
- Analytic or difference-quotient Jacobian
- Mass matrix support for DAEs
- Continuous output via collocation polynomial
- Consistent initial condition computation for index-1 DAEs

## 2. Header Files and Libraries

### Header File

A single header provides the complete public API:

```c
#include "radau5.h"
```

This header includes the necessary SUNDIALS headers internally:
`sundials_context.h`, `sundials_nvector.h`, `sundials_matrix.h`,
`sundials_linearsolver.h`, and `sundials_types.h`.

### SUNDIALS Dependencies

Depending on the linear algebra backend, include the appropriate SUNDIALS
headers in your application:

```c
/* Always needed */
#include <nvector/nvector_serial.h>

/* Dense problems */
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

/* Band problems */
#include <sunmatrix/sunmatrix_band.h>
#include <sunlinsol/sunlinsol_band.h>

/* Sparse problems (KLU) */
#include <sunmatrix/sunmatrix_sparse.h>
#include <sunlinsol/sunlinsol_klu.h>
```

### Linking

Link your application against the RADAU5 library and the required SUNDIALS
libraries:

```
-lradau5 -lsundials_nvecserial -lsundials_sunmatrixdense \
-lsundials_sunmatrixband -lsundials_sunmatrixsparse \
-lsundials_sunlinsoldense -lsundials_sunlinsolband \
-lsundials_sunlinsolklu -lsundials_core -lm
```

Only the libraries corresponding to the matrix/solver types actually used need
to be linked. For sparse support, the KLU library (`-lklu -lamd -lcolamd
-lbtf -lsuitesparseconfig`) must also be available.

## 3. Skeleton Program

The following numbered steps outline the typical calling sequence for a RADAU5
application. Steps marked *(optional)* may be omitted when defaults are
acceptable.

1. **Create the SUNDIALS context**

   ```c
   SUNContext sunctx;
   SUNContext_Create(SUN_COMM_NULL, &sunctx);
   ```

2. **Create initial condition vector**

   ```c
   N_Vector y0 = N_VNew_Serial(n, sunctx);
   sunrealtype* y0data = N_VGetArrayPointer(y0);
   /* fill y0data[0..n-1] with initial values */
   ```

3. **Create the RADAU5 solver object**

   ```c
   void* mem = Radau5Create(sunctx);
   ```

4. **Initialize the solver**

   ```c
   Radau5Init(mem, rhs, t0, y0);
   ```

5. **Create a template Jacobian matrix and attach the linear solver**

   ```c
   /* Dense */
   SUNMatrix J = SUNDenseMatrix(n, n, sunctx);
   Radau5SetLinearSolver(mem, J);

   /* Band (mu = upper bandwidth, ml = lower bandwidth) */
   SUNMatrix J = SUNBandMatrix(n, mu, ml, sunctx);
   Radau5SetLinearSolver(mem, J);

   /* Sparse (KLU) — analytic Jacobian or DQ via sparsity pattern */
   SUNMatrix J = SUNSparseMatrix(n, n, nnz, CSC_MAT, sunctx);
   Radau5SetLinearSolver(mem, J);
   ```

6. **Set tolerances**

   ```c
   Radau5SStolerances(mem, rtol, atol);
   ```

7. **Set optional inputs** *(optional)*

   ```c
   Radau5SetJacFn(mem, jac);           /* analytic Jacobian */
   Radau5SetSparsityPattern(mem, S);   /* sparse DQ Jacobian (alternative to SetJacFn) */
   Radau5SetMassFn(mem, mas, M);       /* mass matrix (DAEs) */
   Radau5SetInitStep(mem, h0);         /* initial step size */
   Radau5SetMaxStep(mem, hmax);        /* maximum step size */
   Radau5SetSchurDecomp(mem, 1);       /* use Schur transform */
   Radau5SetUserData(mem, user_data);  /* user data pointer */
   ```

8. **Compute consistent initial conditions** *(optional, index-1 DAEs)*

   ```c
   Radau5CalcIC(mem, id);
   ```

9. **Advance the solution**

   ```c
   N_Vector yout = N_VNew_Serial(n, sunctx);
   sunrealtype tret;
   int flag = Radau5Solve(mem, tout, yout, &tret);
   ```

10. **Get optional outputs** *(optional)*

    ```c
    long int nstep, naccpt, nfcn, njac;
    Radau5GetNumSteps(mem, &nstep);
    Radau5GetNumAccSteps(mem, &naccpt);
    Radau5GetNumRhsEvals(mem, &nfcn);
    Radau5GetNumJacEvals(mem, &njac);
    ```

11. **Free memory**

    ```c
    N_VDestroy(y0);
    N_VDestroy(yout);
    SUNMatDestroy(J);
    Radau5Free(&mem);
    SUNContext_Free(&sunctx);
    ```

    Note: `Radau5Free` frees all solver-internal memory but does **not** destroy
    user-owned objects (the template Jacobian `J`, the mass matrix `M`, initial
    condition vectors, or the `SUNContext`).

### Minimal Complete Example — Van der Pol Oscillator

```c
#include <stdio.h>
#include <math.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include "radau5.h"

#define EPS 1.0e-6

static int rhs(sunrealtype t, N_Vector y, N_Vector yd, void* ud)
{
  (void)t; (void)ud;
  sunrealtype* yv  = N_VGetArrayPointer(y);
  sunrealtype* ydv = N_VGetArrayPointer(yd);
  ydv[0] = yv[1];
  ydv[1] = ((1.0 - yv[0] * yv[0]) * yv[1] - yv[0]) / EPS;
  return 0;
}

static int jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
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

int main(void)
{
  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  void* mem = Radau5Create(sunctx);
  N_Vector y0 = N_VNew_Serial(2, sunctx);
  N_VGetArrayPointer(y0)[0] = 2.0;
  N_VGetArrayPointer(y0)[1] = 0.0;

  Radau5Init(mem, rhs, 0.0, y0);

  SUNMatrix J = SUNDenseMatrix(2, 2, sunctx);
  Radau5SetLinearSolver(mem, J);
  Radau5SetJacFn(mem, jac);
  Radau5SStolerances(mem, 1e-6, 1e-6);
  Radau5SetInitStep(mem, 1e-6);

  N_Vector yout = N_VNew_Serial(2, sunctx);
  sunrealtype tret;
  int flag = Radau5Solve(mem, 2.0, yout, &tret);

  sunrealtype* yd = N_VGetArrayPointer(yout);
  printf("flag=%d  t=%.6e\n", flag, tret);
  printf("y[0]=%.14e  y[1]=%.14e\n", yd[0], yd[1]);

  long int nstep, naccpt, nfcn;
  Radau5GetNumSteps(mem, &nstep);
  Radau5GetNumAccSteps(mem, &naccpt);
  Radau5GetNumRhsEvals(mem, &nfcn);
  printf("steps=%ld  accepted=%ld  rhs_evals=%ld\n", nstep, naccpt, nfcn);

  N_VDestroy(y0);
  N_VDestroy(yout);
  SUNMatDestroy(J);
  Radau5Free(&mem);
  SUNContext_Free(&sunctx);
  return 0;
}
```

## 4. User-Supplied Function Types

All user-supplied functions return `int`. A return value of `0` indicates
success; a negative value indicates an unrecoverable error that will cause the
solver to halt.

### Radau5RhsFn

Right-hand side function for the ODE system M y' = f(t, y).

```c
typedef int (*Radau5RhsFn)(sunrealtype t, N_Vector y, N_Vector ydot,
                           void* user_data);
```

| Argument    | Description                                      |
|-------------|--------------------------------------------------|
| `t`         | Current value of the independent variable.       |
| `y`         | Current value of the dependent variable vector.  |
| `ydot`      | On output, must contain f(t, y).                 |
| `user_data` | Pointer to user data set via `Radau5SetUserData`. |

### Radau5JacFn

Jacobian function computing J = df/dy.

```c
typedef int (*Radau5JacFn)(sunrealtype t, N_Vector y, N_Vector fy,
                           SUNMatrix Jac, void* user_data,
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
```

| Argument    | Description                                                  |
|-------------|--------------------------------------------------------------|
| `t`         | Current value of the independent variable.                   |
| `y`         | Current value of the dependent variable vector.              |
| `fy`        | Current value of f(t, y), already evaluated by the solver.   |
| `Jac`       | Output Jacobian matrix. Fill entries of df/dy.               |
| `user_data` | Pointer to user data.                                        |
| `tmp1`–`tmp3` | Scratch vectors of length n available for use.             |

If no Jacobian function is provided (the default), the solver uses
difference-quotient (DQ) approximation for dense and band matrices.
For sparse matrices, the solver uses DQ approximation if a sparsity pattern
has been set via `Radau5SetSparsityPattern`; otherwise, a user-supplied
analytic Jacobian is required.

### Radau5MassFn

Mass matrix function for DAE systems M(t) y' = f(t, y).

```c
typedef int (*Radau5MassFn)(sunrealtype t, SUNMatrix M, void* user_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
```

| Argument    | Description                                          |
|-------------|------------------------------------------------------|
| `t`         | Current value of the independent variable.           |
| `M`         | Output mass matrix. Fill entries of M(t).            |
| `user_data` | Pointer to user data.                                |
| `tmp1`–`tmp3` | Scratch vectors of length n available for use.     |

If no mass function is provided, the solver assumes M = I (identity),
i.e., a standard ODE system y' = f(t, y).

### Radau5SolOutFn

Solution output callback, invoked after each accepted step.

```c
typedef int (*Radau5SolOutFn)(long int nr, sunrealtype told, sunrealtype t,
                              N_Vector y, void* radau5_mem, void* user_data);
```

| Argument     | Description                                                    |
|--------------|----------------------------------------------------------------|
| `nr`         | Step number (1, 2, 3, ...).                                    |
| `told`       | Time at the beginning of the step.                             |
| `t`          | Time at the end of the step (current time).                    |
| `y`          | Solution vector at time `t`.                                   |
| `radau5_mem` | Opaque solver memory; pass to `Radau5Contr` for interpolation. |
| `user_data`  | Pointer to user data.                                          |

Return value:

- `0` — continue integration.
- Negative — halt integration; `Radau5Solve` returns `RADAU5_SOLOUT_RETURN`.

Within this callback, `Radau5Contr(radau5_mem, i, s)` may be called to obtain
the i-th component of the solution at any time `s` in `[told, t]`.

## 5. Solver Creation and Destruction

### Radau5Create

```c
void* Radau5Create(SUNContext sunctx);
```

Creates and returns an opaque pointer to the RADAU5 solver memory block.
All internal fields are initialized to safe defaults. Returns `NULL` on failure
(memory allocation error or `NULL` context).

| Argument  | Description                          |
|-----------|--------------------------------------|
| `sunctx`  | A valid SUNDIALS context object.     |

| Return    | Description                          |
|-----------|--------------------------------------|
| `void*`   | Opaque solver memory, or `NULL` on failure. |

### Radau5Free

```c
void Radau5Free(void** radau5_mem);
```

Frees all memory allocated internally by the solver, including cloned vectors,
internal matrices (`E1`, `E2`, `J`, `Jsaved`), linear solvers (`LS_E1`,
`LS_E2`), and scratch vectors. Sets `*radau5_mem` to `NULL`.

Does **not** free user-owned objects: the template Jacobian matrix passed to
`Radau5SetLinearSolver`, the mass matrix `M` passed to `Radau5SetMassFn`,
initial condition vectors, output vectors, or the `SUNContext`.

| Argument     | Description                                    |
|--------------|------------------------------------------------|
| `radau5_mem` | Pointer to the opaque solver memory pointer.   |

## 6. Solver Initialization

### Radau5Init

```c
int Radau5Init(void* radau5_mem, Radau5RhsFn rhs, sunrealtype t0,
               N_Vector y0);
```

Initializes the solver. Clones `y0` to create the internal solution vector and
all scratch/stage vectors, records the problem dimension `n`, stores the RHS
function pointer, and sets the initial time. Must be called before
`Radau5SetLinearSolver` and before `Radau5Solve`.

| Argument     | Description                                      |
|--------------|--------------------------------------------------|
| `radau5_mem` | Opaque solver memory from `Radau5Create`.        |
| `rhs`        | User-supplied RHS function (must not be `NULL`).  |
| `t0`         | Initial value of the independent variable.       |
| `y0`         | Initial condition vector. Cloned internally.     |

| Return                | Description                              |
|-----------------------|------------------------------------------|
| `RADAU5_SUCCESS`      | Initialization succeeded.                |
| `RADAU5_MEM_NULL`     | `radau5_mem` is `NULL`.                  |
| `RADAU5_ILL_INPUT`    | `rhs` or `y0` is `NULL`, or `n <= 0`.   |
| `RADAU5_MEM_FAIL`     | Memory allocation failed.                |

## 7. Linear Solver Setup

### Radau5SetLinearSolver

```c
int Radau5SetLinearSolver(void* radau5_mem, SUNMatrix J);
```

Configures the internal linear algebra based on the type of the template
Jacobian matrix `J`. The solver clones `J` internally and creates the
appropriate linear solvers for the two Newton subsystems:

- **E1** (n x n real): `fac1 * M - J`
- **E2** (2n x 2n real): the realified complex system coupling two n x n blocks

The behavior depends on the matrix type of `J`:

| Matrix type of `J`       | E1 solver       | E2 solver       | Notes                          |
|--------------------------|-----------------|-----------------|--------------------------------|
| `SUNDenseMatrix`         | Dense (n x n)   | Dense (2n x 2n) | DQ Jacobian available          |
| `SUNBandMatrix`          | Band (n x n)    | Dense (2n x 2n) | DQ Jacobian available          |
| `SUNSparseMatrix` (CSC)  | KLU (n x n)     | KLU (2n x 2n)   | Analytic Jacobian or DQ via `SetSparsityPattern` |

Must be called after `Radau5Init` and before `Radau5Solve`.

| Argument     | Description                                                  |
|--------------|--------------------------------------------------------------|
| `radau5_mem` | Opaque solver memory.                                        |
| `J`          | Template Jacobian matrix. Cloned internally; user retains ownership. |

| Return                | Description                                      |
|-----------------------|--------------------------------------------------|
| `RADAU5_SUCCESS`      | Linear solver setup succeeded.                   |
| `RADAU5_MEM_NULL`     | `radau5_mem` is `NULL`.                          |
| `RADAU5_ILL_INPUT`    | `J` is `NULL` or has unsupported type.           |
| `RADAU5_MEM_FAIL`     | Memory allocation or linear solver creation failed. |

For band matrices, the upper and lower bandwidths are extracted from `J` and
stored internally. The E2 system is always dense (2n x 2n) in the band case
because the block coupling destroys the band structure.

For sparse matrices, the E2 system is also sparse (CSC format). The sparsity
pattern of E2 is constructed from the Jacobian and mass matrix patterns. For
sparse matrices, the user may either provide an analytic Jacobian via
`Radau5SetJacFn`, or provide a sparsity pattern via `Radau5SetSparsityPattern`
to enable automatic difference-quotient Jacobian computation using column
grouping (Curtis-Powell-Reid technique).

## 8. Tolerance Specification

### Radau5SStolerances

```c
int Radau5SStolerances(void* radau5_mem, sunrealtype rtol, sunrealtype atol);
```

Sets scalar relative and absolute tolerances. The error weight vector is
computed as `scal_i = atol + rtol * |y_i|` and updated at each step.

| Argument     | Description                          |
|--------------|--------------------------------------|
| `radau5_mem` | Opaque solver memory.                |
| `rtol`       | Scalar relative tolerance (> 0).     |
| `atol`       | Scalar absolute tolerance (> 0).     |

| Return             | Description                              |
|--------------------|------------------------------------------|
| `RADAU5_SUCCESS`   | Tolerances set successfully.             |
| `RADAU5_MEM_NULL`  | `radau5_mem` is `NULL`.                  |
| `RADAU5_ILL_INPUT` | Negative tolerance value.                |

### Radau5SVtolerances

```c
int Radau5SVtolerances(void* radau5_mem, N_Vector rtol, N_Vector atol);
```

Sets vector relative and absolute tolerances for component-wise error control.
Both vectors must have length n.

| Argument     | Description                                    |
|--------------|------------------------------------------------|
| `radau5_mem` | Opaque solver memory.                          |
| `rtol`       | Vector of relative tolerances (all entries > 0). |
| `atol`       | Vector of absolute tolerances (all entries > 0). |

| Return             | Description                              |
|--------------------|------------------------------------------|
| `RADAU5_SUCCESS`   | Tolerances set successfully.             |
| `RADAU5_MEM_NULL`  | `radau5_mem` is `NULL`.                  |
| `RADAU5_ILL_INPUT` | `NULL` vector or negative entries.       |

**Note on tolerance transformation:** On the first call to `Radau5Solve`, the
solver internally transforms the relative tolerance as
`rtol_internal = 0.1 * rtol^(2/3)`, following the original Fortran RADAU5
convention. This transformation modifies the stored tolerances in-place (the
original values are not preserved). The transformation is applied once and
flagged so it is not repeated on subsequent `Radau5Solve` calls. If tolerances
are re-set via `Radau5SStolerances` or `Radau5SVtolerances` between
`Radau5Solve` calls, the flag is reset and the new tolerances will be
transformed on the next call.

## 9. Optional Input Functions

The following table summarizes all optional input functions. Each returns
`RADAU5_SUCCESS` on success, `RADAU5_MEM_NULL` if the solver memory is `NULL`,
or `RADAU5_ILL_INPUT` for invalid arguments.

| Function                      | Key Argument(s)          | Default        | Description                                    |
|-------------------------------|--------------------------|----------------|------------------------------------------------|
| `Radau5SetJacFn`              | `jac`                    | `NULL` (DQ)    | Analytic Jacobian callback                     |
| `Radau5SetMassFn`             | `mas`, `M`               | `NULL` (M=I)   | Mass matrix callback and template matrix       |
| `Radau5SetSolOutFn`           | `solout`                 | `NULL`         | Solution output callback                       |
| `Radau5SetUserData`           | `user_data`              | `NULL`         | User data pointer passed to all callbacks      |
| `Radau5SetMaxNumSteps`        | `mxsteps`                | 100000         | Maximum number of step attempts                |
| `Radau5SetMaxNewtonIter`      | `maxnit`                 | 7              | Maximum Newton iterations per step             |
| `Radau5SetInitStep`           | `h0`                     | 1e-6           | Initial step size                              |
| `Radau5SetMaxStep`            | `hmax`                   | 0.0 (no limit) | Maximum allowed step size                      |
| `Radau5SetSafetyFactor`       | `safe`                   | 0.9            | Step size safety factor                        |
| `Radau5SetStepSizeController` | `pred`                   | 1              | 1 = Gustafsson, 2 = classical                  |
| `Radau5SetDAEIndex`           | `nind1`, `nind2`, `nind3`| n, 0, 0        | DAE index-component partition                  |
| `Radau5SetStartNewton`        | `startn`                 | 0              | 0 = extrapolate, 1 = zero starting values      |
| `Radau5SetSchurDecomp`        | `use_schur`              | 0              | 0 = eigenvalue decomp, 1 = Schur decomp        |
| `Radau5SetSparsityPattern`    | `S`                      | `NULL`         | Sparsity pattern for sparse DQ Jacobian (column grouping) |

### Radau5SetJacFn

```c
int Radau5SetJacFn(void* radau5_mem, Radau5JacFn jac);
```

Sets the user-supplied Jacobian function. If `jac` is `NULL` (the default), the
solver uses internal difference-quotient approximation (dense and band only).

| Argument     | Description                                |
|--------------|--------------------------------------------|
| `radau5_mem` | Opaque solver memory.                      |
| `jac`        | Jacobian function, or `NULL` for DQ.       |

### Radau5SetMassFn

```c
int Radau5SetMassFn(void* radau5_mem, Radau5MassFn mas, SUNMatrix M);
```

Sets the mass matrix function and provides a template mass matrix for DAE
problems of the form M(t) y' = f(t, y). The matrix `M` is **not** cloned — the
solver stores the pointer directly, and the user retains ownership. `M` is
**not** freed by `Radau5Free`.

| Argument     | Description                                          |
|--------------|------------------------------------------------------|
| `radau5_mem` | Opaque solver memory.                                |
| `mas`        | Mass matrix function.                                |
| `M`          | Template mass matrix (user-owned, not freed by solver). |

### Radau5SetSolOutFn

```c
int Radau5SetSolOutFn(void* radau5_mem, Radau5SolOutFn solout);
```

Sets the solution output callback invoked after each accepted step. Within the
callback, `Radau5Contr` may be called for continuous output.

| Argument     | Description                                |
|--------------|--------------------------------------------|
| `radau5_mem` | Opaque solver memory.                      |
| `solout`     | Output callback, or `NULL` to disable.     |

### Radau5SetUserData

```c
int Radau5SetUserData(void* radau5_mem, void* user_data);
```

Sets a pointer to user-defined data that is passed as the last argument to all
user-supplied callback functions (RHS, Jacobian, mass, solout).

| Argument     | Description                                |
|--------------|--------------------------------------------|
| `radau5_mem` | Opaque solver memory.                      |
| `user_data`  | Pointer to user data.                      |

### Radau5SetMaxNumSteps

```c
int Radau5SetMaxNumSteps(void* radau5_mem, long int mxsteps);
```

Sets the maximum number of step attempts allowed during a single call to
`Radau5Solve`. If this limit is reached, the solver returns
`RADAU5_TOO_MANY_STEPS`.

| Argument     | Description                                |
|--------------|--------------------------------------------|
| `radau5_mem` | Opaque solver memory.                      |
| `mxsteps`    | Maximum step attempts (default: 100000).   |

### Radau5SetMaxNewtonIter

```c
int Radau5SetMaxNewtonIter(void* radau5_mem, int maxnit);
```

Sets the maximum number of Newton iterations per step. If Newton fails to
converge within this limit, the step is rejected and retried with a smaller
step size.

| Argument     | Description                                |
|--------------|--------------------------------------------|
| `radau5_mem` | Opaque solver memory.                      |
| `maxnit`     | Maximum Newton iterations (default: 7).    |

### Radau5SetInitStep

```c
int Radau5SetInitStep(void* radau5_mem, sunrealtype h0);
```

Sets the initial step size. If not set, the default is 1e-6.

| Argument     | Description                                |
|--------------|--------------------------------------------|
| `radau5_mem` | Opaque solver memory.                      |
| `h0`         | Initial step size (> 0).                   |

### Radau5SetMaxStep

```c
int Radau5SetMaxStep(void* radau5_mem, sunrealtype hmax);
```

Sets the maximum allowed step size. A value of 0.0 (the default) means no
upper bound is imposed.

| Argument     | Description                                |
|--------------|--------------------------------------------|
| `radau5_mem` | Opaque solver memory.                      |
| `hmax`       | Maximum step size (>= 0; 0 = no limit).   |

### Radau5SetSafetyFactor

```c
int Radau5SetSafetyFactor(void* radau5_mem, sunrealtype safe);
```

Sets the safety factor used in step size prediction. The new step size is
computed as `h_new = h * safe * (tol / err)^(1/p)` where `p` depends on the
method order.

| Argument     | Description                                |
|--------------|--------------------------------------------|
| `radau5_mem` | Opaque solver memory.                      |
| `safe`       | Safety factor (default: 0.9). Typically in (0, 1). |

### Radau5SetStepSizeController

```c
int Radau5SetStepSizeController(void* radau5_mem, int pred);
```

Selects the step size control strategy.

| Value | Strategy    | Description                                              |
|-------|-------------|----------------------------------------------------------|
| 1     | Gustafsson  | Predictive controller using the ratio of consecutive errors (default). |
| 2     | Classical   | Standard controller based on the current error estimate only.          |

| Argument     | Description                                |
|--------------|--------------------------------------------|
| `radau5_mem` | Opaque solver memory.                      |
| `pred`       | Controller selection: 1 or 2.              |

### Radau5SetDAEIndex

```c
int Radau5SetDAEIndex(void* radau5_mem, sunindextype nind1,
                      sunindextype nind2, sunindextype nind3);
```

Specifies the DAE index-component partition. The first `nind1` components are
index-1 variables, the next `nind2` are index-2, and the next `nind3` are
index-3. The sum `nind1 + nind2 + nind3` must equal `n`. For pure ODE systems,
the default `(n, 0, 0)` is appropriate.

| Argument     | Description                                        |
|--------------|----------------------------------------------------|
| `radau5_mem` | Opaque solver memory.                              |
| `nind1`      | Number of index-1 components (default: n).         |
| `nind2`      | Number of index-2 components (default: 0).         |
| `nind3`      | Number of index-3 components (default: 0).         |

### Radau5SetStartNewton

```c
int Radau5SetStartNewton(void* radau5_mem, int startn);
```

Controls the starting values for the Newton iteration at each step.

| Value | Behavior                                                          |
|-------|-------------------------------------------------------------------|
| 0     | Use extrapolation from the previous step to predict starting values (default). |
| 1     | Use zero starting values for the Newton iteration.                |

| Argument     | Description                                |
|--------------|--------------------------------------------|
| `radau5_mem` | Opaque solver memory.                      |
| `startn`     | 0 = extrapolate, 1 = zero.                |

### Radau5SetSchurDecomp

```c
int Radau5SetSchurDecomp(void* radau5_mem, int use_schur);
```

Selects the transform used to decouple the 3-stage Newton system into
independent linear systems.

| Value | Transform              | Description                                    |
|-------|------------------------|------------------------------------------------|
| 0     | Eigenvalue (default)   | Classical Hairer-Wanner T/TI eigenvector transform. E1 and E2 are solved independently. |
| 1     | Schur                  | Orthogonal US/TS transform. Uses block back-substitution: solve E1 for the third stage, then substitute into E2 for the first two stages. |

The Schur decomposition can improve numerical stability for very stiff problems
at the cost of slightly more coupling between the linear solves.

| Argument     | Description                                |
|--------------|--------------------------------------------|
| `radau5_mem` | Opaque solver memory.                      |
| `use_schur`  | 0 = eigenvalue, 1 = Schur.                |

### Radau5SetSparsityPattern

```c
int Radau5SetSparsityPattern(void* radau5_mem, SUNMatrix S);
```

Sets the sparsity pattern for automatic sparse difference-quotient (DQ) Jacobian
computation using the Curtis-Powell-Reid column grouping technique. The solver
groups structurally independent columns so that multiple Jacobian columns can be
computed with a single RHS evaluation, reducing the cost from `n` evaluations to
`ngroups` (typically much less than `n` for sparse problems).

The column grouping is computed immediately when this function is called. Two
column orderings are tried (natural order and reverse-degree order) with
first-fit greedy graph coloring, and the ordering producing fewer groups is kept.

`S` must be a `SUNMATRIX_SPARSE` matrix in CSC format with dimensions n x n
matching the problem size. Its sparsity pattern (column pointers and row indices)
must be identical to the Jacobian template `J` passed to `Radau5SetLinearSolver`.
Typically, the same `SUNMatrix` is passed to both functions.

If both `Radau5SetJacFn` and `Radau5SetSparsityPattern` are called, the analytic
Jacobian takes priority (the sparsity pattern data is kept but unused).

Must be called after `Radau5Init` and before `Radau5Solve`. May be called
multiple times; each call frees previous grouping data and recomputes.

| Argument     | Description                                                  |
|--------------|--------------------------------------------------------------|
| `radau5_mem` | Opaque solver memory.                                        |
| `S`          | Sparse matrix (CSC) whose pattern defines the Jacobian sparsity. User retains ownership. |

| Return                | Description                                      |
|-----------------------|--------------------------------------------------|
| `RADAU5_SUCCESS`      | Sparsity pattern set and column grouping computed. |
| `RADAU5_MEM_NULL`     | `radau5_mem` is `NULL`.                          |
| `RADAU5_ILL_INPUT`    | `S` is `NULL`, not sparse, not CSC, or wrong dimensions. |
| `RADAU5_MEM_FAIL`     | Memory allocation failed during grouping computation. |

## 10. Main Solver Function

### Radau5Solve

```c
int Radau5Solve(void* radau5_mem, sunrealtype tout, N_Vector yout,
                sunrealtype* tret);
```

Integrates the ODE/DAE system from the current internal time toward `tout`.
On the first call, the solver applies the internal tolerance transformation
(see Section 8). The solver advances by adaptive steps until `tout` is reached
or an error/stop condition occurs.

On successful return, `yout` contains the solution at time `*tret`, and the
solver's internal state is updated so that subsequent calls continue from
`*tret`.

| Argument     | Description                                                    |
|--------------|----------------------------------------------------------------|
| `radau5_mem` | Opaque solver memory.                                          |
| `tout`       | Desired output time. Must satisfy `tout > t_current`.          |
| `yout`       | On return, contains the solution vector at time `*tret`.       |
| `tret`       | On return, the time actually reached by the solver.            |

| Return                     | Description                                          |
|----------------------------|------------------------------------------------------|
| `RADAU5_SUCCESS`           | Successfully reached `tout`.                         |
| `RADAU5_TSTOP_RETURN`      | Successfully reached a stop time.                    |
| `RADAU5_SOLOUT_RETURN`      | Integration halted by the `SolOut` callback.         |
| `RADAU5_MEM_NULL`          | `radau5_mem` is `NULL`.                              |
| `RADAU5_ILL_INPUT`         | Invalid input (e.g., solver not initialized).        |
| `RADAU5_TOO_MANY_STEPS`   | Maximum number of steps exceeded (`mxstep`).         |
| `RADAU5_STEP_TOO_SMALL`   | Step size became too small to make progress.         |
| `RADAU5_SINGULAR_MATRIX`  | Repeated singular matrix in LU decomposition.        |
| `RADAU5_CONV_FAILURE`     | Newton iteration failed to converge after step size reduction. |
| `RADAU5_LSETUP_FAIL`      | Linear solver setup (factorization) failed.          |
| `RADAU5_LSOLVE_FAIL`      | Linear solve failed.                                 |
| `RADAU5_RHSFUNC_FAIL`     | The user-supplied RHS function returned an error.    |
| `RADAU5_MEM_FAIL`         | Memory allocation failed during integration.         |

## 11. Optional Output Functions

The following functions retrieve solver statistics after a call to
`Radau5Solve`. Each returns `RADAU5_SUCCESS` on success or `RADAU5_MEM_NULL`
if the solver memory is `NULL`.

| Function                  | Output type      | Description                                          |
|---------------------------|------------------|------------------------------------------------------|
| `Radau5GetNumSteps`       | `long int*`      | Total number of step attempts (accepted + rejected). |
| `Radau5GetNumAccSteps`    | `long int*`      | Number of accepted steps.                            |
| `Radau5GetNumRejSteps`    | `long int*`      | Number of rejected steps.                            |
| `Radau5GetNumRhsEvals`    | `long int*`      | Total RHS evaluations (including Newton and error estimate). |
| `Radau5GetNumJacEvals`    | `long int*`      | Number of Jacobian evaluations (analytic or DQ).     |
| `Radau5GetNumDecomps`     | `long int*`      | Number of LU decompositions (E1 and E2 counted separately). |
| `Radau5GetNumLinSolves`   | `long int*`      | Number of linear solves.                             |
| `Radau5GetNumNewtonIters` | `long int*`      | Total Newton iterations across all steps.            |
| `Radau5GetCurrentStep`    | `sunrealtype*`   | Current (most recent) step size.                     |
| `Radau5GetCurrentTime`    | `sunrealtype*`   | Current internal time.                               |

### Detailed Signatures

```c
int Radau5GetNumSteps(void* radau5_mem, long int* nsteps);
int Radau5GetNumAccSteps(void* radau5_mem, long int* naccpt);
int Radau5GetNumRejSteps(void* radau5_mem, long int* nrejct);
int Radau5GetNumRhsEvals(void* radau5_mem, long int* nfcn);
int Radau5GetNumJacEvals(void* radau5_mem, long int* njac);
int Radau5GetNumDecomps(void* radau5_mem, long int* ndec);
int Radau5GetNumLinSolves(void* radau5_mem, long int* nsol);
int Radau5GetNumNewtonIters(void* radau5_mem, long int* nnewt);
int Radau5GetCurrentStep(void* radau5_mem, sunrealtype* hcur);
int Radau5GetCurrentTime(void* radau5_mem, sunrealtype* tcur);
```

All getters follow the same pattern:

| Argument     | Description                                    |
|--------------|------------------------------------------------|
| `radau5_mem` | Opaque solver memory.                          |
| (output ptr) | On return, contains the requested statistic.   |

| Return             | Description                              |
|--------------------|------------------------------------------|
| `RADAU5_SUCCESS`   | Value retrieved successfully.            |
| `RADAU5_MEM_NULL`  | `radau5_mem` is `NULL`.                  |

**Note on `ndec`:** Each call to `SUNLinSolSetup` increments the decomposition
counter once. Since E1 and E2 are factored separately, the C implementation
reports approximately twice the decomposition count compared to the original
Fortran code (which counts the pair as one decomposition).

**Note on `nfcn`:** The RHS evaluation count includes all evaluations: the
initial evaluation, 3 evaluations per Newton iteration (one at each collocation
point), and any additional evaluations for the error estimate refinement.

## 12. Continuous Output

### Radau5Contr

```c
sunrealtype Radau5Contr(void* radau5_mem, sunindextype i, sunrealtype t);
```

Returns the i-th component of the continuous output solution at time `t` using
polynomial interpolation based on the collocation stages. The interpolant is
a degree-3 polynomial that is accurate to the order of the method within the
current step interval.

This function is **only valid when called from within a `Radau5SolOutFn`
callback**, for times `t` in the interval `[t_old, t_new]` where `t_old` and
`t_new` are the step boundaries passed to the callback.

| Argument     | Description                                                  |
|--------------|--------------------------------------------------------------|
| `radau5_mem` | Opaque solver memory (the `radau5_mem` argument of `SolOut`).|
| `i`          | Component index (0-based: 0, 1, ..., n-1).                  |
| `t`          | Time at which to evaluate the interpolant, `t_old <= t <= t_new`. |

| Return         | Description                                            |
|----------------|--------------------------------------------------------|
| `sunrealtype`  | The interpolated value of y_i(t).                      |

Example usage within a `SolOut` callback:

```c
static int my_solout(long int nr, sunrealtype told, sunrealtype t,
                     N_Vector y, void* radau5_mem, void* user_data)
{
  /* Evaluate component 0 at 10 equally spaced points in [told, t] */
  for (int k = 0; k <= 10; k++) {
    sunrealtype s = told + (t - told) * k / 10.0;
    sunrealtype y0_s = Radau5Contr(radau5_mem, 0, s);
    printf("  t=%.6e  y[0]=%.14e\n", s, y0_s);
  }
  return 0;  /* continue integration */
}
```

## 13. Consistent Initial Conditions

### Radau5CalcIC

```c
int Radau5CalcIC(void* radau5_mem, N_Vector id);
```

Computes consistent initial conditions for index-1 DAE systems of the form
M y' = f(t, y) where M is singular. This function solves for the algebraic
variables so that the algebraic constraints are satisfied at the initial time.

The `id` vector identifies which components are differential and which are
algebraic:

| `id[i]` value | Component type |
|---------------|----------------|
| 1.0           | Differential   |
| 0.0           | Algebraic      |

The function modifies the internal solution vector (the initial conditions
passed to `Radau5Init`) in place. It must be called after `Radau5Init`,
`Radau5SetLinearSolver`, and `Radau5SetMassFn`, but before `Radau5Solve`.

| Argument     | Description                                                  |
|--------------|--------------------------------------------------------------|
| `radau5_mem` | Opaque solver memory.                                        |
| `id`         | Vector of length n with 1.0 for differential, 0.0 for algebraic components. |

| Return             | Description                                          |
|--------------------|------------------------------------------------------|
| `RADAU5_SUCCESS`   | Consistent initial conditions computed successfully. |
| `RADAU5_MEM_NULL`  | `radau5_mem` is `NULL`.                              |
| `RADAU5_ILL_INPUT` | `id` is `NULL` or solver not properly initialized.   |
| `RADAU5_IC_FAIL`   | The IC computation failed to converge.               |

## 14. Return Codes

The following table lists all return codes defined in `radau5.h`.

| Constant                  | Value | Description                                                        |
|---------------------------|-------|--------------------------------------------------------------------|
| `RADAU5_SUCCESS`          |   0   | Successful return.                                                 |
| `RADAU5_TSTOP_RETURN`     |   1   | Integration reached a prescribed stop time.                        |
| `RADAU5_SOLOUT_RETURN`     |   2   | Integration halted by the user's `SolOut` callback (negative return). |
| `RADAU5_MEM_NULL`         |  -1   | The `radau5_mem` argument is `NULL`.                               |
| `RADAU5_ILL_INPUT`        |  -2   | An illegal or invalid input was provided.                          |
| `RADAU5_TOO_MANY_STEPS`  |  -3   | The maximum number of step attempts (`mxstep`) was exceeded.       |
| `RADAU5_STEP_TOO_SMALL`  |  -4   | The step size fell below the machine-precision-limited minimum.     |
| `RADAU5_SINGULAR_MATRIX` |  -5   | Repeated singular matrix encountered during LU decomposition.      |
| `RADAU5_CONV_FAILURE`    |  -6   | Newton iteration failed to converge after step size reductions.    |
| `RADAU5_NEWT_PREDICT`    | -60   | Newton predicted slow convergence; step size already reduced internally. This is an internal code not normally seen by users. |
| `RADAU5_LSETUP_FAIL`     |  -7   | The linear solver setup (factorization) returned an error.         |
| `RADAU5_LSOLVE_FAIL`     |  -8   | The linear solve returned an error.                                |
| `RADAU5_RHSFUNC_FAIL`    |  -9   | The user-supplied RHS function returned a non-zero (error) value.  |
| `RADAU5_MEM_FAIL`        | -10   | A memory allocation failed.                                        |
| `RADAU5_IC_FAIL`         | -11   | Consistent initial condition computation failed.                   |
