# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

C implementation of the RADAU solver — a variable-order (3, 5, or 7 stages; orders 5, 9, 13) implicit Runge-Kutta method (Radau IIA family) for stiff ODEs and DAEs. Faithfully translated from the Fortran code by Hairer & Wanner ("Solving ODEs II"), using SUNDIALS abstractions (N_Vector, SUNMatrix, SUNLinearSolver) to support dense, band, and sparse (KLU) linear algebra.

Supports variable stage count via `Radau5SetNumStages(mem, ns)` where ns ∈ {3, 5, 7}:
- **ns=3** (default): 3-stage, order 5 — equivalent to original RADAU5
- **ns=5**: 5-stage, order 9
- **ns=7**: 7-stage, order 13

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
ctest   # runs all 201 tests (problems × modes × order configs)
```

The bash `build.sh` hardcodes paths to SUNDIALS build at `/mnt/d/workon/sundials/build` and uses clang. It compiles all radau5 source files together with each example into a standalone executable under `bin/`.

## Architecture

### Solver Flow (single step)

```
radau5_Step (radau5_step.c)
  ├─ Entry: skipdecomp? → ikeep/label30 | caljac? → label20 | → label10
  ├─ label10: Compute Jacobian (analytic or DQ)
  ├─ label20: Variable-order selection (if variab), then Build E1/E2, factor all
  ├─ label30: Newton iteration (radau5_newt.c)
  │    ├─ Evaluate RHS at ns collocation points (t+c[k]*h, k=0..ns-1)
  │    ├─ Forward transform (TI or US^T) → form linear system RHS
  │    ├─ Eigenvalue mode: solve E1 + E2[0..npairs-1] independently
  │    ├─ Schur mode: block back-substitution (E1 for last stage, then E2[k] pairs)
  │    ├─ Back-transform (T or US) → update Z[0..ns-1]
  │    └─ Convergence check (theta, faccon, fnewt)
  ├─ Error estimate (radau5_estrad.c)
  ├─ Step size control (Gustafsson or classical, exponent 1/(ns+1))
  └─ Accept/reject → update h, caljac, skipdecomp flags
```

### Key Design Decisions

- **Variable-order support**: `Radau5SetNumStages(mem, ns)` selects a fixed ns=3 (order 5), ns=5 (order 9), or ns=7 (order 13). `Radau5SetOrderLimits(mem, nsmin, nsmax)` enables runtime variable-order selection. Both must be called before `Radau5Init`. Default is ns=3 for backward compatibility. The solver pre-allocates `RADAU5_NPAIRS_MAX=3` complex E2 systems and `RADAU5_NS_MAX=7` stage vectors regardless of the initial ns, so that `radau5_ChangeOrder` can switch orders without reallocation. Step-size exponent is `1/(ns+1)`.
- **Variable-order selection** (Fortran radau.f lines 900-935): When `variab=1` (nsmin < nsmax), the solver evaluates order-change criteria at label20 before each E1/E2 build. Order increases when `newt_prev > 1 && thetat <= vitu && hquot ∈ (hhod, hhou)` and `ichan > 10`. Order decreases when `thetat >= vitd` or unexpected rejection/Newton failure occurred. Parameters: `vitu=0.002`, `vitd=0.8`, `hhou=1.2`, `hhod=0.8`. After an order change, `radau5_ChangeOrder` reloads constants, updates `nit` and `fnewt`, and zeros starting values. The Gustafsson controller is skipped on the first step after an order change.
- **E2 is 2n×2n real** (not n×n complex) because SUNDIALS has no complex matrix support. For sparse J, E2 is also sparse. One E2 per complex eigenvalue pair: `E2[0..npairs-1]`.
- **Eigenvalue mode**: E2[k] = `[[alphn[k]*M-J, -betan[k]*M],[betan[k]*M, alphn[k]*M-J]]` (antisymmetric off-diagonal). E1 and all E2[k] are solved independently.
- **Schur mode**: E2[k] uses the k-th 2×2 diagonal block of TS: `[[TS[r0][r0]/h*M-J, TS[r0][r1]/h*M],[TS[r1][r0]/h*M, TS[r1][r1]/h*M-J]]` where r0=2k, r1=2k+1. Uses block back-substitution from bottom-right 1×1 block upward through each 2×2 block, with TS upper-triangular coupling terms.
- **Tolerance transformation**: Fortran RADAU5 internally transforms `rtol = 0.1 * rtol^(2/3)`. This is done once in `Radau5Solve` on first call (`tol_transformed` flag). The transformation modifies tolerances in-place (original values are not preserved). If the user re-sets tolerances via `Radau5SStolerances` or `Radau5SVtolerances` between `Radau5Solve` calls, `tol_transformed` is reset to 0 so the new tolerances are properly transformed on the next call.
- **Step reuse**: `skipdecomp` flag allows skipping both Jacobian and LU decomposition when `qt = hnew/h ∈ [quot1, quot2]` (default [1.0, 1.2]).
- **Newton slow-convergence**: Returns `RADAU5_NEWT_PREDICT` (not `CONV_FAILURE`) when Newton predicts slow convergence and already reduced h — caller skips the label78 h*=0.5 path.
- **DQ Jacobian for sparse**: Supported via column grouping (Curtis-Powell-Reid). User calls `Radau5SetSparsityPattern(mem, S)` to provide the CSC sparsity pattern; the solver computes a graph coloring to group structurally independent columns, then evaluates multiple Jacobian columns per RHS call. If no sparsity pattern is set and no analytic Jacobian is provided, returns `RADAU5_ILL_INPUT`.
- **RHS error handling**: RHS return value convention follows SUNDIALS: `0 = success, > 0 = recoverable, < 0 = unrecoverable`. Recoverable errors during Newton iteration, DQ Jacobian, and error estimation route through label78 (halve step, retry). Only singular matrix and RHS recoverable errors increment `nsing` (max 5 before fatal). Newton convergence failures do not count toward `nsing`. Unrecoverable errors and post-acceptance/initial RHS failures terminate immediately.
- **DQ Jacobian increment (dense)**: Uses Fortran radau5.f formula `DELT=sqrt(UROUND*max(1e-5,|y_j|))` rather than CVODE-style formula. This matches the original Fortran behavior for stiff DAE problems.
- **Discontinuity handling**: `Radau5ResetForDiscontinuity(mem, h0)` resets solver state at discontinuity boundaries (h, first, reject, nsing, faccon), matching Fortran radau5 entry behavior. The `first` flag is automatically cleared after the first accepted step, restoring extrapolation for subsequent steps.
- **Mass matrix M**: User-owned (not destroyed by `Radau5Free`).
- **Sparse M with sparse J**: When both M and J are sparse, E1's sparsity pattern is the union of J and M, computed eagerly in `Radau5SetLinearSolver(mem, J, M)` via `radau5_SparseUnion`. E2 NNZ = 2×nnz_union + 2×nnz_M. The M parameter is only needed for sparse pattern computation; pass NULL when J is dense/band or when there is no mass matrix.
- **Continuous output**: Uses Newton divided-difference algorithm (general for all ns). After each accepted step, computes `cont[0..ns]` coefficients where `cont[0]=ycur` and `cont[1..ns]` are divided differences on nodes `c[ns-1], c[ns-2], ..., c[0], 0`. Evaluation via `Radau5Contr` uses Horner form: `s=(t-xsol)/hsol+1`, then descends from `cont[ns]` through `cont[0]`. This matches the Fortran `CONTRA` function exactly.
- **Step extrapolation**: Newton starting values for the next step use the continuous output polynomial evaluated at new collocation points: `z[k] = c[k]*hquot * P(c[k]*hquot)` where P is the Horner evaluation of `cont[1..ns]`. The forward transform `f = TI * z` is then applied. This general algorithm (from Fortran radau.f ns=5/7 blocks) works for all ns values.
- **Event detection (rootfinding)**: User provides `Radau5RootFn g(t, y, gout, user_data)` via `Radau5RootInit(mem, nrtfn, g)`. After each accepted step, the solver evaluates g at the step endpoint and checks for sign changes vs the step start. If detected, the Illinois method (modified regula falsi) refines the root location using `Radau5Contr` dense output interpolation — no extra RHS evaluations needed. Solver stops at the root and returns `RADAU5_ROOT_RETURN = 3`. User queries `Radau5GetRootInfo` for which functions fired (+1 rising, -1 falling), then resumes with another `Radau5Solve` call. Direction filtering via `Radau5SetRootDirection`. Root tolerance: `100*uround*max(|tlo|,|thi|)` (matches SUNDIALS).

### Source Files

| File | Purpose |
|------|---------|
| `radau5.c` | Create/Free/Init/Solve, all setters/getters, tolerance transformation |
| `radau5_impl.h` | `Radau5Mem_` struct, internal function prototypes |
| `radau5_linsys.c` | Method constants (eigen + Schur), DQ Jacobian (dense/band/sparse), BuildE1/E2, DecompE1/E2, ChangeOrder, ComputeScal, MassMult, SparseUnion/SparseLookup |
| `radau5_newt.c` | Generalized Newton iteration for variable ns with eigenvalue/Schur branches |
| `radau5_step.c` | Single step attempt: Jac reuse logic, variable-order selection, Newton, error estimate, accept/reject, general extrapolation |
| `radau5_estrad.c` | Error estimation (ESTRAV from Fortran, general ns-term loop) |
| `radau5_contr.c` | Continuous output: Newton divided-difference coefficients (general for all ns) |
| `radau5_colgroup.c` | Column grouping (graph coloring) for sparse DQ Jacobian: first-fit greedy coloring with two orderings |
| `radau5_colgroup.h` | Internal header for column grouping functions |
| `radau5_ic.c` | Consistent initial conditions for index-1 DAEs |
| `radau5_root.c` | Rootfinding (event detection): sign-change detection, Illinois method refinement using continuous output |
| `radau5_root.h` | Internal header for rootfinding functions |

### Schur Decomposition Details

The Schur option decomposes `A^{-1} = US * TS * US'` where US is orthogonal and TS is upper quasi-triangular. For ns=3:
```
TS = [ 2.6811  -8.4239  -4.0886 ]    US = [ 0.1387   0.0463   0.9893 ]
     [ 1.1046   2.6811   4.7002 ]         [-0.2296  -0.9702   0.0776 ]
     [ 0        0        3.6378 ]         [-0.9633   0.2379   0.1239 ]
```
For ns=5: TS is 5×5 with two 2×2 blocks + one 1×1 block. For ns=7: TS is 7×7 with three 2×2 blocks + one 1×1 block.

- `use_schur` flag in `Radau5Mem_`, set via `Radau5SetSchurDecomp()` (re-runs `radau5_InitConstants`)
- `radau5_InitConstants()`: populates `US_mat[ns*ns]`, `TS_mat[ns*ns]`, overrides `u1=TS[ns-1][ns-1]`, sets `T=US`, `TI=US^T`
- `radau5_BuildE2(rmem, pair_idx, ...)`: uses the pair_idx-th 2×2 diagonal block of TS: `TS[r0:r1, r0:r1]` where `r0=2*pair_idx`
- `radau5_newt.c`: Schur path uses full TS row scalings + block back-substitution from bottom-right 1×1 block upward; back-transform uses full ns×ns US
- `radau5_estrad.c`, `radau5_contr.c`: no changes needed (work in physical Z-space)

### Fortran Correspondence

The C code closely follows the Fortran structure. Key label mappings in `radau5_step.c`:
- **label10** = Fortran label 10 (compute Jacobian)
- **label20** = Fortran label 20 (variable-order selection + build/factor E1, E2)
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
| `nit` | 7 | NIT (max Newton iterations, auto-adjusted for ns) |
| `pred` | 1 | PRED (1=Gustafsson, 2=classical) |
| `mxstep` | 100000 | NMAX |
| `vitu` | 0.002 | VITU (order increase threshold) |
| `vitd` | 0.8 | VITD (order decrease threshold) |
| `hhou` | 1.2 | HHOU (step ratio upper bound for order increase) |
| `hhod` | 0.8 | HHOD (step ratio lower bound for order increase) |

## Testing

All examples accept command-line arguments: `rtol atol h0 use_schur nsmin nsmax` (positional, optional). Default `nsmin=3 nsmax=7` enables variable-order mode. The `pollu` example has an extra `use_dq` argument before `nsmin nsmax`. Rootfinding examples (bounce, roberts_root, rocket, orbit, kepler, kneeode) take only `nsmin nsmax`.

```shell
# Run with defaults (variable order, nsmin=3, nsmax=7)
./bin/radau5_vdpol

# Run with fixed ns=5 and Schur mode
./bin/radau5_vdpol 1e-6 1e-6 1e-6 1 5 5

# Run with variable order
./bin/radau5_vdpol 1e-6 1e-6 1e-6 0 3 7

# CTest: 201 tests (problems × schur modes × order configs)
cd build && ctest
ctest -R eigen           # eigenvalue mode only
ctest -R schur           # Schur mode only
ctest -R ns5             # fixed ns=5 only
ctest -R ns7             # fixed ns=7 only
ctest -R nsvar           # variable order only
ctest -R vdpol           # all modes for vdpol

# Compare with Fortran baseline
cd /mnt/d/workon/sundials/thirdpart/IVPtestset_2.4/build
echo "1e-6 1e-6 1e-6" | ./radau5_f_vdpol
```

CTest configurations per problem:
- `radau5_<prob>_eigen_ns3` / `radau5_<prob>_schur_ns3` — fixed order 5 (classic RADAU5)
- `radau5_<prob>_eigen_ns5` / `radau5_<prob>_schur_ns5` — fixed order 9
- `radau5_<prob>_eigen_ns7` / `radau5_<prob>_schur_ns7` — fixed order 13
- `radau5_<prob>_eigen_nsvar` / `radau5_<prob>_schur_nsvar` — variable order 3→7

DAE-3 problems (andrews, caraxis) and tba skip ns=7 tests (known incompatibility with high-order methods for index-3 DAEs; tba is too slow at ns=7).

Unit tests: `test_radau5_constants`, `test_radau5_api`, `test_radau5_build_e1`, `test_radau5_build_e2`, `test_radau5_dq_jac`, `test_radau5_colgroup`, `test_radau5_dq_jac_sparse`.

### Examples (24 problems)

#### Original 13 IVPtestset problems

| Example | Type | n | Matrix | Default rtol/atol/h0 |
|---------|------|---|--------|---------------------|
| vdpol | ODE | 2 | dense | 1e-6 / 1e-6 / 1e-6 |
| rober | ODE | 3 | dense | 1e-10 / 1e-14 / 1e-12 |
| orego | ODE | 3 | dense | 1e-6 / 1e-6 / 1e-6 |
| e5 | ODE | 4 | dense | 1e-4 / 1e-24 / 1e-6 |
| hires | ODE | 8 | dense (DQ) | 1e-9 / 1e-9 / 1e-9 |
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

#### Event detection (rootfinding)

| Example | Type | n | Matrix | Default rtol/atol/h0 |
|---------|------|---|--------|---------------------|
| bounce | ODE | 2 | dense | 1e-10 / 1e-12 / 1e-4 |
| roberts_root | ODE | 3 | dense | 1e-4 / 1e-8 / 1e-4 |
| rocket | ODE | 2 | dense | 1e-5 / 1e-2 / 1e-3 |
| orbit | ODE | 4 | dense (DQ) | 1e-5 / 1e-4 / 1e-3 |
| kepler | ODE | 4 | dense (DQ) | 1e-8 / 1e-10 / 1e-3 |
| kneeode | ODE | 1 | dense | 1e-3 / 1e-6 / 1e-6 |


## Conventions

- Public API: `Radau5CamelCase`. Internal: `radau5_snake_case`. Macros: `RADAU5_UPPER`.
- All floating-point via `sunrealtype`, indices via `sunindextype`, constants via `SUN_RCONST()`.
- Fortran line references in comments (e.g., "Fortran radau5.f lines 793-1101").
- `ndec` counter increments once per `SUNLinSolSetup` call (E1 and E2 counted separately, so C ndec ≈ 2× Fortran ndec).
- `nfcn` counts all RHS evaluations including those inside Newton (ns per iteration) and error estimate refinement.

## Documentation

Detailed documentation in `doc/`:

| Document | Description |
|----------|-------------|
| [doc/theory.md](doc/theory.md) | Mathematical theory: Radau IIA method, Newton iteration (eigenvalue + Schur), error estimation, step size control, tolerance transformation, continuous output, DAE support |
| [doc/user_guide.md](doc/user_guide.md) | API reference: skeleton program, all function signatures, optional inputs/outputs, return codes, complete Van der Pol example |
| [doc/examples.md](doc/examples.md) | 13 example problems: equations, parameters, usage, reference solutions. Covers stiff/non-stiff ODEs, band/sparse, DAE index-1/2/3 |

## Sparse Mass Matrix Support

When both J and M are `SUNMATRIX_SPARSE` (CSC), the solver handles the union sparsity pattern automatically:

1. Allocate sparse Mt template with correct CSC sparsity pattern
2. `Radau5SetLinearSolver(mem, Jt, Mt)` — computes `union(J, M)` pattern for E1 and allocates E2 with `2×nnz_union + 2×nnz_M`
3. `Radau5SetMassFn(mem, mas, Mt)` — registers the mass function callback and M template

For dense/band J (even with a mass matrix), pass NULL as M to `Radau5SetLinearSolver` since the sparsity pattern is irrelevant for dense/band allocation.

Helper functions in `radau5_linsys.c`:
- `radau5_SparseUnion(A, B, sunctx)` — merge two CSC sparsity patterns
- `radau5_SparseLookup(A, row, col)` — binary search for value in CSC column

## Sparse DQ Jacobian via Column Grouping

When the Jacobian matrix is sparse (CSC) and no analytic Jacobian is provided, the solver can compute a finite-difference Jacobian using the Curtis-Powell-Reid column grouping technique. This groups structurally independent columns so multiple columns are computed per RHS evaluation.

Usage:
1. `Radau5SetLinearSolver(mem, J, NULL)` — pass a CSC sparse template as usual (NULL for M since no sparse mass matrix)
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
