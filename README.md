# RADAU5 — Variable-Order Implicit Runge-Kutta Solver (C / SUNDIALS)

A C implementation of the RADAU solver — a variable-order (3, 5, or 7 stages; orders 5, 9, 13) implicit Runge-Kutta method (Radau IIA family) for stiff ODEs and DAEs. Faithfully translated from the Fortran code by E. Hairer and G. Wanner (*Solving Ordinary Differential Equations II*, Springer 1996), using [SUNDIALS](https://github.com/LLNL/sundials) abstractions (N_Vector, SUNMatrix, SUNLinearSolver) for flexible linear algebra backends.

## Features

- **Variable-order**: ns=3 (order 5), ns=5 (order 9), ns=7 (order 13), with runtime order switching
- **Dense, band, and sparse (KLU)** Jacobian support
- **DQ (difference quotient) Jacobian** with Curtis-Powell-Reid column grouping for sparse systems
- **DAE support** — index-1, 2, and 3 with user-supplied mass matrix (dense or sparse)
- **Two Newton decoupling modes**: eigenvalue decomposition (classical) or Schur decomposition
- **Event detection (rootfinding)** via Illinois method with dense output interpolation
- **Discontinuity handling** via `Radau5ResetForDiscontinuity`
- **Continuous output** polynomial for dense output between steps
- **24 example problems** from the IVP Test Set and beyond

## Dependencies

- [SUNDIALS](https://github.com/LLNL/sundials) (v6.0+ recommended)
- CMake 3.18+
- C compiler (C99)
- Optional: SuiteSparse/KLU for sparse linear algebra

## Build

```bash
mkdir build && cd build
cmake .. -DSUNDIALS_DIR=/path/to/sundials/install
make -j$(nproc)
ctest   # runs all 201 tests (problems × modes × order configs)
```

## Quick Example

```c
#include "radau5.h"

/* RHS: y' = -1000*y */
int rhs(sunrealtype t, N_Vector y, N_Vector yd, void* ud) {
  N_VGetArrayPointer(yd)[0] = -1000.0 * N_VGetArrayPointer(y)[0];
  return 0;
}

int main() {
  SUNContext sunctx;
  SUNContext_Create(SUN_COMM_NULL, &sunctx);

  void* mem = Radau5Create(sunctx);
  Radau5SetOrderLimits(mem, 3, 7);  /* variable order: ns=3→7 */

  N_Vector y0 = N_VNew_Serial(1, sunctx);
  N_VGetArrayPointer(y0)[0] = 1.0;

  Radau5Init(mem, rhs, 0.0, y0);

  SUNMatrix J = SUNDenseMatrix(1, 1, sunctx);
  Radau5SetLinearSolver(mem, J, NULL);
  Radau5SStolerances(mem, 1e-8, 1e-10);

  N_Vector yout = N_VNew_Serial(1, sunctx);
  sunrealtype tret;
  Radau5Solve(mem, 1.0, yout, &tret);

  printf("y(1) = %e\n", N_VGetArrayPointer(yout)[0]);

  N_VDestroy(y0); N_VDestroy(yout);
  SUNMatDestroy(J); Radau5Free(&mem);
  SUNContext_Free(&sunctx);
  return 0;
}
```

## API Overview

```c
/* Create / Destroy */
void* Radau5Create(SUNContext sunctx);
void  Radau5Free(void** mem);

/* Initialize */
int Radau5Init(void* mem, Radau5RhsFn rhs, sunrealtype t0, N_Vector y0);

/* Linear solver setup (M is only needed when both J and M are sparse) */
int Radau5SetLinearSolver(void* mem, SUNMatrix J, SUNMatrix M);

/* Tolerances */
int Radau5SStolerances(void* mem, sunrealtype rtol, sunrealtype atol);

/* Optional configuration */
int Radau5SetNumStages(void* mem, int ns);         /* fixed order: ns ∈ {3,5,7} */
int Radau5SetOrderLimits(void* mem, int nsmin, int nsmax); /* variable order */
int Radau5SetJacFn(void* mem, Radau5JacFn jac);
int Radau5SetMassFn(void* mem, Radau5MassFn mas, SUNMatrix M);
int Radau5SetSchurDecomp(void* mem, int use_schur);
int Radau5RootInit(void* mem, int nrtfn, Radau5RootFn g);

/* Solve */
int Radau5Solve(void* mem, sunrealtype tout, N_Vector yout, sunrealtype* tret);

/* Continuous output (callable from SolOut callback) */
sunrealtype Radau5Contr(void* mem, sunindextype i, sunrealtype t);
```

## Examples

24 example programs covering stiff ODEs, banded/sparse systems, DAEs of index 1–3, and event detection:

| Problem | Type | n | Matrix |
|---------|------|---|--------|
| vdpol | ODE | 2 | dense |
| rober | ODE | 3 | dense |
| orego | ODE | 3 | dense |
| e5 | ODE | 4 | dense |
| hires | ODE | 8 | dense (DQ) |
| plei | ODE | 28 | dense |
| ringmod | ODE | 15 | dense |
| heat1d | ODE | 100 | band |
| medakzo | ODE | 400 | band |
| pollu | ODE | 20 | sparse/KLU |
| chemakzo | DAE-1 | 6 | dense + mass |
| transamp | DAE-1 | 8 | dense + mass |
| transamp_sparse | DAE-1 | 8 | sparse + sparse mass |
| andrews | DAE-3 | 27 | dense + mass |
| caraxis | DAE-3 | 10 | dense + mass |
| crank | DAE-2 | 24 | dense + mass |
| water | DAE-2 | 49 | dense + mass |
| fekete | DAE-2 | 160 | dense + mass |
| tba | DAE-1 | 350 | dense + mass |
| bounce | ODE | 2 | dense (rootfinding) |
| rocket | ODE | 2 | dense (rootfinding) |
| orbit | ODE | 4 | dense (rootfinding) |
| kepler | ODE | 4 | dense (rootfinding) |
| kneeode | ODE | 1 | dense (rootfinding) |

See [doc/examples.md](doc/examples.md) for full details.

## Documentation

- [doc/theory.md](doc/theory.md) — Mathematical theory
- [doc/user_guide.md](doc/user_guide.md) — API reference
- [doc/examples.md](doc/examples.md) — Example problems guide

## License

BSD-3-Clause. See [LICENSE](LICENSE).

## Acknowledgments

Based on the Fortran RADAU5/RADAU code by Ernst Hairer and Gerhard Wanner. Uses the [SUNDIALS](https://github.com/LLNL/sundials) library developed at Lawrence Livermore National Laboratory.
