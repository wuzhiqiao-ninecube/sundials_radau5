# RADAU5 — Implicit Runge-Kutta Solver (C / SUNDIALS)

A C implementation of the RADAU5 solver — a 3-stage, order-5 implicit Runge-Kutta method (Radau IIA family) for stiff ODEs and DAEs. Faithfully translated from the Fortran code by E. Hairer and G. Wanner (*Solving Ordinary Differential Equations II*, Springer 1996), using [SUNDIALS](https://github.com/LLNL/sundials) abstractions (N_Vector, SUNMatrix, SUNLinearSolver) for flexible linear algebra backends.

## Features

- **Dense, band, and sparse (KLU)** Jacobian support
- **DQ (difference quotient) Jacobian** with Curtis-Powell-Reid column grouping for sparse systems
- **DAE support** — index-1, 2, and 3 with user-supplied mass matrix (dense or sparse)
- **Two Newton decoupling modes**: eigenvalue decomposition (classical) or Schur decomposition
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
ctest   # runs all 48 tests (24 problems × 2 modes)
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
  N_Vector y0 = N_VNew_Serial(1, sunctx);
  N_VGetArrayPointer(y0)[0] = 1.0;

  Radau5Init(mem, rhs, 0.0, y0);

  SUNMatrix J = SUNDenseMatrix(1, 1, sunctx);
  Radau5SetLinearSolver(mem, J);
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

## Examples

24 example programs covering stiff ODEs, banded/sparse systems, and DAEs of index 1–3:

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
| andrews | DAE-3 | 27 | dense + mass |
| tba | DAE-1 | 350 | dense + mass |
| ... | | | |

See [doc/examples.md](doc/examples.md) for full details.

## Documentation

- [doc/theory.md](doc/theory.md) — Mathematical theory
- [doc/user_guide.md](doc/user_guide.md) — API reference
- [doc/examples.md](doc/examples.md) — Example problems guide

## License

BSD-3-Clause. See [LICENSE](LICENSE).

## Acknowledgments

Based on the Fortran RADAU5 code by Ernst Hairer and Gerhard Wanner. Uses the [SUNDIALS](https://github.com/LLNL/sundials) library developed at Lawrence Livermore National Laboratory.
