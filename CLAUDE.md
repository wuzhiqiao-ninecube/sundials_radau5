# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

C implementation of the RADAU5 solver — a 3-stage, order-5 implicit Runge-Kutta method (Radau IIA family) for stiff ODEs and DAEs. Faithfully translated from the Fortran code by Hairer & Wanner ("Solving ODEs II"), using SUNDIALS abstractions (N_Vector, SUNMatrix, SUNLinearSolver) to support dense, band, and sparse (KLU) linear algebra.

Supports two transform modes for the Newton system decoupling:
- **Eigenvalue decomposition** (default): classical Hairer-Wanner T/TI eigenvector transform
- **Schur decomposition** (optional): orthogonal US/TS transform via `Radau5SetSchurDecomp(mem, 1)`

## Build

```shell
# Quick build via bash script (clang, static linking)
cd /mnt/d/workon/sundials/thirdpart/radau5
bash build.sh radau5_vdpol radau5_rober   # builds into bin/

# CMake build (requires SUNDIALS installed)
mkdir build && cd build
cmake .. -DSUNDIALS_DIR=/path/to/sundials/install
make -j$(nproc)
ctest   # runs all 48 tests (24 problems × 2 modes)
```

The bash `build.sh` hardcodes paths to SUNDIALS build at `/mnt/d/workon/sundials/build` and uses clang. It compiles all radau5 source files together with each example into a standalone executable under `bin/`.

## Architecture

### Solver Flow (single step)

```
radau5_Step (radau5_step.c)
  ├─ Entry: skipdecomp? → label30 | caljac? → label20 | → label10
  ├─ label10: Compute Jacobian (analytic or DQ)
  ├─ label20: Build E1=fac1*M-J, E2=coupled 2n×2n; factor both
  ├─ label30: Newton iteration (radau5_newt.c)
  │    ├─ Evaluate RHS at 3 collocation points (t+c1*h, t+c2*h, t+h)
  │    ├─ Forward transform (TI or US^T) → form linear system RHS
  │    ├─ Eigenvalue mode: solve E1 + E2 independently
  │    ├─ Schur mode: block back-substitution (E1 for dF3, then E2 for [dF1;dF2])
  │    ├─ Back-transform (T or US) → update Z1,Z2,Z3
  │    └─ Convergence check (theta, faccon, fnewt)
  ├─ Error estimate (radau5_estrad.c)
  ├─ Step size control (Gustafsson or classical)
  └─ Accept/reject → update h, caljac, skipdecomp flags
```

### Key Design Decisions

- **E2 is 2n×2n real** (not n×n complex) because SUNDIALS has no complex matrix support. For sparse J, E2 is also sparse.
- **Eigenvalue mode**: E2 = `[[alphn*M-J, -betan*M],[betan*M, alphn*M-J]]` (antisymmetric off-diagonal). E1 and E2 are solved independently.
- **Schur mode**: E2 = `[[TS00/h*M-J, TS01/h*M],[TS10/h*M, TS11/h*M-J]]` (asymmetric off-diagonal). Uses block back-substitution: solve E1 for dF3 first, substitute into rhs1/rhs2, then solve E2 for [dF1;dF2]. Transform uses orthogonal US (T=US, TI=US^T).
- **Tolerance transformation**: Fortran RADAU5 internally transforms `rtol = 0.1 * rtol^(2/3)`. This is done once in `Radau5Solve` on first call (`tol_transformed` flag). The transformation modifies tolerances in-place (original values are not preserved). If the user re-sets tolerances via `Radau5SStolerances` or `Radau5SVtolerances` between `Radau5Solve` calls, `tol_transformed` is reset to 0 so the new tolerances are properly transformed on the next call.
- **Step reuse**: `skipdecomp` flag allows skipping both Jacobian and LU decomposition when `qt = hnew/h ∈ [quot1, quot2]` (default [1.0, 1.2]).
- **Newton slow-convergence**: Returns `RADAU5_NEWT_PREDICT` (not `CONV_FAILURE`) when Newton predicts slow convergence and already reduced h — caller skips the label78 h*=0.5 path.
- **DQ Jacobian for sparse**: Supported via column grouping (Curtis-Powell-Reid). User calls `Radau5SetSparsityPattern(mem, S)` to provide the CSC sparsity pattern; the solver computes a graph coloring to group structurally independent columns, then evaluates multiple Jacobian columns per RHS call. If no sparsity pattern is set and no analytic Jacobian is provided, returns `RADAU5_ILL_INPUT`.
- **RHS error handling**: RHS return value convention follows SUNDIALS: `0 = success, > 0 = recoverable, < 0 = unrecoverable`. Recoverable errors during Newton iteration, DQ Jacobian, and error estimation route through label78 (halve step, retry). Only singular matrix and RHS recoverable errors increment `nsing` (max 5 before fatal). Newton convergence failures do not count toward `nsing`. Unrecoverable errors and post-acceptance/initial RHS failures terminate immediately.
- **DQ Jacobian increment (dense)**: Uses Fortran radau5.f formula `DELT=sqrt(UROUND*max(1e-5,|y_j|))` rather than CVODE-style formula. This matches the original Fortran behavior for stiff DAE problems.
- **Discontinuity handling**: `Radau5ResetForDiscontinuity(mem, h0)` resets solver state at discontinuity boundaries (h, first, reject, nsing, faccon), matching Fortran radau5 entry behavior. The `first` flag is automatically cleared after the first accepted step, restoring extrapolation for subsequent steps.
- **Mass matrix M**: User-owned (not destroyed by `Radau5Free`).
- **Sparse M with sparse J**: When both M and J are sparse, E1's sparsity pattern is the union of J and M (computed lazily in `Radau5Solve` after mass matrix evaluation via `radau5_SparseUnion`). E2 NNZ = 2×nnz_union + 2×nnz_M. KLU solvers are recreated after the pattern merge.

### Source Files

| File | Purpose |
|------|---------|
| `radau5.c` | Create/Free/Init/Solve, all setters/getters, tolerance transformation |
| `radau5_impl.h` | `Radau5Mem_` struct, internal function prototypes |
| `radau5_linsys.c` | Method constants (eigen + Schur), DQ Jacobian (dense/band/sparse), BuildE1/E2, DecompE1/E2, ComputeScal, MassMult, SparseUnion/SparseLookup |
| `radau5_newt.c` | Simplified Newton iteration with eigenvalue/Schur branches |
| `radau5_step.c` | Single step attempt: Jac reuse logic, Newton, error estimate, accept/reject |
| `radau5_estrad.c` | Error estimation (ESTRAD from Fortran, transform-independent) |
| `radau5_contr.c` | Continuous output polynomial coefficients |
| `radau5_colgroup.c` | Column grouping (graph coloring) for sparse DQ Jacobian: first-fit greedy coloring with two orderings |
| `radau5_colgroup.h` | Internal header for column grouping functions |
| `radau5_ic.c` | Consistent initial conditions for index-1 DAEs |

### Schur Decomposition Details

The Schur option decomposes `A^{-1} = US * TS * US'` where US is orthogonal and TS is upper quasi-triangular:
```
TS = [ 2.6811  -8.4239  -4.0886 ]    US = [ 0.1387   0.0463   0.9893 ]
     [ 1.1046   2.6811   4.7002 ]         [-0.2296  -0.9702   0.0776 ]
     [ 0        0        3.6378 ]         [-0.9633   0.2379   0.1239 ]
```
- `use_schur` flag in `Radau5Mem_`, set via `Radau5SetSchurDecomp()` (re-runs `radau5_InitConstants`)
- `radau5_InitConstants()`: populates `US[3][3]`, `TS[3][3]`, overrides `u1=TS[2][2]`, sets `T=US`, `TI=US^T`
- `radau5_BuildE2()`: unified with generic `a00,a01,a10,a11` coefficients (eigenvalue: antisymmetric; Schur: from TS block)
- `radau5_newt.c`: Schur path uses full TS row scalings + block back-substitution; back-transform uses full 3×3 US
- `radau5_estrad.c`, `radau5_contr.c`: no changes needed (work in physical Z-space)

### Fortran Correspondence

The C code closely follows the Fortran structure. Key label mappings in `radau5_step.c`:
- **label10** = Fortran label 10 (compute Jacobian)
- **label20** = Fortran label 20 (build + factor E1, E2)
- **label30** = Fortran label 30 (Newton + error estimate + accept/reject)
- **label78** = Fortran label 78 (Newton failure / singular matrix recovery)

### Default Parameters (matching Fortran)

| Parameter | Default | Fortran variable |
|-----------|---------|-----------------|
| `safe` | 0.9 | SAFE |
| `facl` | 5.0 | FACL (upper bound on quot) |
| `facr` | 0.125 | FACR (lower bound on quot = 1/8) |
| `thet` | 0.001 | THET (Jacobian reuse threshold) |
| `quot1` | 1.0 | QUOT1 |
| `quot2` | 1.2 | QUOT2 |
| `nit` | 7 | NIT (max Newton iterations) |
| `pred` | 1 | PRED (1=Gustafsson, 2=classical) |
| `mxstep` | 100000 | NMAX |

## Testing

All 24 examples accept command-line arguments: `rtol atol h0 use_schur` (positional, optional). The `pollu` example additionally accepts a 5th argument `use_dq` (0 or 1) to test sparse DQ Jacobian via column grouping.

```shell
# Run with defaults
./bin/radau5_vdpol

# Run with custom tolerances and Schur mode
./bin/radau5_vdpol 1e-8 1e-8 1e-6 1

# CTest: all 48 tests (24 problems × eigenvalue + Schur)
cd build && ctest
ctest -R eigen           # eigenvalue mode only
ctest -R schur           # Schur mode only
ctest -R vdpol           # both modes for vdpol

# Compare with Fortran baseline
cd /mnt/d/workon/sundials/thirdpart/IVPtestset_2.4/build
echo "1e-6 1e-6 1e-6" | ./radau5_f_vdpol
```

Unit tests: `test_radau5_constants`, `test_radau5_api`, `test_radau5_build_e1`, `test_radau5_build_e2`, `test_radau5_dq_jac`, `test_radau5_colgroup`, `test_radau5_dq_jac_sparse`, `test_radau5_dahlquist`.

### Examples (24 problems)

#### Original 13 IVPtestset problems

| Example | Type | n | Matrix | Default rtol/atol/h0 |
|---------|------|---|--------|---------------------|
| vdpol | ODE | 2 | dense | 1e-6 / 1e-6 / 1e-6 |
| rober | ODE | 3 | dense | 1e-10 / 1e-14 / 1e-12 |
| orego | ODE | 3 | dense | 1e-6 / 1e-6 / 1e-6 |
| e5 | ODE | 4 | dense | 1e-4 / 1e-24 / 1e-6 |
| hires | ODE | 8 | dense (DQ) | 1e-6 / 1e-6 / 1e-6 |
| plei | ODE | 28 | dense | 1e-6 / 1e-6 / 1e-6 |
| ringmod | ODE | 15 | dense | 1e-7 / 1e-7 / 1e-6 |
| heat1d | ODE | 100 | band(1,1) DQ | 1e-10 / 1e-12 / 1e-6 |
| medakzo | ODE | 400 | band(2,2) DQ | 1e-7 / 1e-7 / 1e-6 |
| pollu | ODE | 20 | sparse/KLU (DQ or analytic) | 1e-6 / 1e-6 / 1e-6 |
| chemakzo | DAE-1 | 6 | dense + mass | 1e-7 / 1e-7 / 1e-10 |
| transamp | DAE-1 | 8 | dense + mass | 1e-6 / 1e-6 / 1e-6 |
| andrews | DAE-3 | 27 | dense + mass | 1e-6 / 1e-6 / 1e-10 |

#### Sparse mass matrix variant

| Example | Type | n | Matrix | Default rtol/atol/h0 |
|---------|------|---|--------|---------------------|
| transamp_sparse | DAE-1 | 8 | sparse J + sparse M / KLU | 1e-6 / 1e-6 / 1e-6 |

#### Additional IVPtestset problems

| Example | Type | n | Matrix | Default rtol/atol/h0 |
|---------|------|---|--------|---------------------|
| vdpolm | ODE | 2 | dense (mu=1000) | 1e-6 / 1e-6 / 1e-6 |
| caraxis | DAE-3 | 10 | dense + diag mass | 1e-6 / 1e-6 / 1e-10 |
| crank | DAE-2 | 24 | dense + diag mass | 1e-6 / 1e-6 / 1e-10 |
| water | DAE-2 | 49 | dense + diag mass | 1e-7 / 1e-7 / 1e-6 |
| beam | ODE | 80 | dense (DQ) | 1e-6 / 1e-6 / 1e-6 |
| fekete | DAE-2 | 160 | dense (DQ) + diag mass | 1e-6 / 1e-6 / 1e-6 |
| tba | DAE-1 | 350 | dense (DQ) + diag mass | 1e-5 / 1e-5 / 4e-5 |

#### ARKODE comparison

| Example | Type | n | Matrix | Default rtol/atol/h0 |
|---------|------|---|--------|---------------------|
| ark_analytic | ODE | 1 | dense (lambda=-1e6) | 1e-5 / 1e-10 / 1e-4 |


## Conventions

- Public API: `Radau5CamelCase`. Internal: `radau5_snake_case`. Macros: `RADAU5_UPPER`.
- All floating-point via `sunrealtype`, indices via `sunindextype`, constants via `SUN_RCONST()`.
- Fortran line references in comments (e.g., "Fortran radau5.f lines 793-1101").
- `ndec` counter increments once per `SUNLinSolSetup` call (E1 and E2 counted separately, so C ndec ≈ 2× Fortran ndec).
- `nfcn` counts all RHS evaluations including those inside Newton (3 per iteration) and error estimate refinement.

## Documentation

Detailed documentation in `doc/`:

| Document | Description |
|----------|-------------|
| [doc/theory.md](doc/theory.md) | Mathematical theory: Radau IIA method, Newton iteration (eigenvalue + Schur), error estimation, step size control, tolerance transformation, continuous output, DAE support |
| [doc/user_guide.md](doc/user_guide.md) | API reference: skeleton program, all function signatures, optional inputs/outputs, return codes, complete Van der Pol example |
| [doc/examples.md](doc/examples.md) | 13 example problems: equations, parameters, usage, reference solutions. Covers stiff/non-stiff ODEs, band/sparse, DAE index-1/2/3 |

## Sparse Mass Matrix Support

When both J and M are `SUNMATRIX_SPARSE` (CSC), the solver handles the union sparsity pattern automatically:

1. `Radau5SetLinearSolver(mem, J)` initially clones E1 from J's pattern (M not yet known)
2. `Radau5SetMassFn(mem, mas, M)` stores M
3. On first call to `Radau5Solve`, after mass matrix evaluation, if M is sparse:
   - E1 is rebuilt with `radau5_SparseUnion(J, M)` pattern
   - E2 NNZ = 2×nnz_union + 2×nnz_M
   - KLU solvers are recreated for the new matrices
   - Controlled by `sparse_ls_finalized` flag (one-time operation)

Helper functions in `radau5_linsys.c`:
- `radau5_SparseUnion(A, B, sunctx)` — merge two CSC sparsity patterns
- `radau5_SparseLookup(A, row, col)` — binary search for value in CSC column

## Sparse DQ Jacobian via Column Grouping

When the Jacobian matrix is sparse (CSC) and no analytic Jacobian is provided, the solver can compute a finite-difference Jacobian using the Curtis-Powell-Reid column grouping technique. This groups structurally independent columns so multiple columns are computed per RHS evaluation.

Usage:
1. `Radau5SetLinearSolver(mem, J)` — pass a CSC sparse template as usual
2. `Radau5SetSparsityPattern(mem, S)` — pass a CSC sparse matrix whose pattern defines which entries are nonzero. Typically the same matrix as `J`. This triggers column grouping computation immediately.
3. Do NOT call `Radau5SetJacFn` — omitting it activates the DQ path

The sparsity pattern `S` and the Jacobian template `J` must have identical CSC structure (same `colptrs`, same `rowinds`).

Column grouping algorithm (`radau5_colgroup.c`):
- Builds a row-to-column reverse mapping from the CSC pattern
- Runs first-fit greedy coloring in two column orderings (natural + reverse-degree)
- Keeps whichever ordering produces fewer groups
- All-zero columns are excluded (marked as group = -1)

`radau5_DQJacSparse` (`radau5_linsys.c`):
- Iterates over groups; for each group, perturbs all columns simultaneously
- One RHS evaluation per group (CVODE-style increment formula)
- Extracts Jacobian entries directly from the CSC pattern — O(nnz) total work
- For the pollu problem (n=20, nnz=86): ~7 groups instead of 20 column-by-column evaluations
