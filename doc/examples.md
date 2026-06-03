# RADAU5 Examples Guide

## Introduction

The RADAU5 solver ships with 36 registered example programs covering stiff and non-stiff ODEs, banded and sparse systems, DAEs of index 1 through 3, event detection, and larger PDE/circuit benchmarks. Many problems come from the IVP Test Set; additional examples exercise band/sparse linear algebra, sparse mass matrices, rootfinding, discontinuity handling, and high-dimensional benchmark systems.

Most examples accept positional command-line arguments:

```bash
./radau5_<prob> [rtol] [atol] [h0] [use_schur] [nsmin] [nsmax]
```

- `rtol` — relative tolerance (default varies per problem)
- `atol` — absolute tolerance (default varies per problem)
- `h0` — initial step size (default varies per problem)
- `use_schur` — `0` for eigenvalue decomposition, `1` for Schur decomposition
- `nsmin`, `nsmax` — stage-count limits; each must be one of `3, 5, 7, 9, 11, 13`

Exceptions:

- `pollu` accepts an extra `use_dq` argument before the order limits:
  ```bash
  ./radau5_pollu [rtol] [atol] [h0] [use_schur] [use_dq] [nsmin] [nsmax]
  ```
- Rootfinding examples (`bounce`, `roberts_root`, `rocket`, `orbit`, `kepler`, `kneeode`) use hardcoded tolerances and accept only:
  ```bash
  ./radau5_<root_prob> [nsmin] [nsmax]
  ```

Build a single example with the shell script:

```bash
bash build.sh radau5_vdpol
./bin/radau5_vdpol
```

Or build all examples and run the regression matrix via CMake:

```bash
mkdir build && cd build
cmake .. -DSUNDIALS_DIR=/path/to/sundials/install
make -j$(nproc)
ctest          # problems x transform modes x order configurations
```

CTest registers the standard examples in eigenvalue and Schur modes across fixed-order configurations (`ns3`, `ns5`, `ns7`, `ns9`, `ns11`, `ns13`) and variable order (`nsvar`). Reduced-order tests (`ns3`, `ns5`, `nsvar` with `nsmin=3`, `nsmax=7`) are used for DAE-3 and large/expensive examples such as `andrews`, `caraxis`, `tba`, `tba_smooth`, and `ks`. Rootfinding examples are tested across order configurations without `rtol/atol/h0/use_schur` arguments.

## Summary Table

| Example | Type | n | Matrix | Jacobian / Notes | Default rtol / atol / h0 |
|---------|------|---|--------|------------------|--------------------------|
| vdpol | ODE | 2 | dense | analytic | 1e-6 / 1e-6 / 1e-6 |
| rober | ODE | 3 | dense | analytic | 1e-10 / 1e-14 / 1e-12 |
| orego | ODE | 3 | dense | analytic | 1e-6 / 1e-6 / 1e-6 |
| e5 | ODE | 4 | dense | analytic | 1e-4 / 1e-24 / 1e-6 |
| hires | ODE | 8 | dense | analytic | 1e-9 / 1e-9 / 1e-9 |
| plei | ODE | 28 | dense | analytic | 1e-6 / 1e-6 / 1e-6 |
| pollu | ODE | 20 | sparse/KLU | analytic or sparse DQ | 1e-6 / 1e-6 / 1e-6 |
| plate | ODE | 80 | sparse/KLU | sparse benchmark | 1e-10 / 1e-10 / 1e-2 |
| ringmod | ODE | 15 | dense | analytic | 1e-7 / 1e-7 / 1e-6 |
| heat1d | ODE | 100 | band(1,1) | DQ | 1e-10 / 1e-12 / 1e-6 |
| medakzo | ODE | 400 | band(2,2) | DQ | 1e-7 / 1e-7 / 1e-6 |
| chemakzo | DAE-1 | 6 | dense + mass | analytic | 1e-7 / 1e-7 / 1e-10 |
| transamp | DAE-1 | 8 | dense + mass | analytic | 1e-6 / 1e-6 / 1e-6 |
| transamp_sparse | DAE-1 | 8 | sparse J + sparse M/KLU | analytic | 1e-6 / 1e-6 / 1e-6 |
| andrews | DAE-3 | 27 | dense + mass | analytic; reduced-order tests | 1e-6 / 1e-6 / 1e-10 |
| ark_analytic | ODE | 1 | dense | analytic, ARKODE comparison | 1e-5 / 1e-10 / 1e-4 |
| vdpolm | ODE | 2 | dense | analytic, mu=1000 | 1e-6 / 1e-6 / 1e-6 |
| caraxis | DAE-3 | 10 | dense + mass | reduced-order tests | 1e-6 / 1e-6 / 1e-10 |
| beam | ODE | 80 | dense | DQ | 1e-6 / 1e-6 / 1e-6 |
| crank | DAE-2 | 24 | dense + mass | analytic | 1e-6 / 1e-6 / 1e-10 |
| water | DAE-2 | 49 | dense + mass | analytic | 1e-7 / 1e-7 / 1e-6 |
| fekete | DAE-2 | 160 | dense + mass | DQ | 1e-6 / 1e-6 / 1e-6 |
| pump | DAE-2 | 9 | dense + mass | analytic | 1e-6 / 1e-6 / 1e-12 |
| pump_smooth | DAE-2 | 9 | dense + mass | analytic, smoothed | 1e-7 / 1e-7 / 1e-3 |
| tba | DAE-1 | 350 | dense + mass | DQ, discontinuities; reduced-order tests | 1e-5 / 1e-5 / 4e-5 |
| tba_smooth | DAE-1 | 350 | dense + mass | DQ, smoothed; reduced-order tests | 1e-5 / 1e-5 / 4e-5 |
| tba2 | DAE-1 | 350 | dense + mass | DQ | problem default |
| bounce | ODE | 2 | dense | rootfinding | hardcoded |
| roberts_root | ODE | 3 | dense | rootfinding | hardcoded |
| rocket | ODE | 2 | dense | rootfinding | hardcoded |
| orbit | ODE | 4 | dense | rootfinding | hardcoded |
| kepler | ODE | 4 | dense | rootfinding | hardcoded |
| kneeode | ODE | 1 | dense | rootfinding | hardcoded |
| ks | ODE | 1022 | dense | DQ; reduced-order tests | 1e-8 / 1e-8 / 1e-6 |
| bruss | ODE | 2 | dense | analytic | 1e-6 / 1e-6 / 1e-6 |
| bruss_2d | ODE | 8192 | sparse/KLU | sparse PDE benchmark | 1e-6 / 1e-6 / 1e-4 |

## Problem Groups

### Original stiff ODE benchmarks

- `vdpol` — Van der Pol oscillator with strong stiffness.
- `rober` — Robertson chemical kinetics on a long integration interval.
- `orego` — Oregonator oscillatory chemical kinetics.
- `e5` — extremely stiff DETEST chemical kinetics problem.
- `hires` — HIRES chemical kinetics problem, using a dense analytic Jacobian.
- `ringmod` — ring modulator electronic circuit.

### Non-stiff and large ODE benchmarks

- `plei` — Pleiades 7-body gravitational problem.
- `beam` — 80-variable beam model using dense DQ Jacobians.
- `vdpolm` — Van der Pol variant with `mu=1000`.
- `ks` — 1022-variable Kuramoto-Sivashinsky benchmark; CTest uses reduced-order configurations for runtime.
- `bruss` — two-variable Brusselator reaction model.

### Banded and sparse ODE benchmarks

- `heat1d` — 1D heat equation using band storage.
- `medakzo` — 400-variable Akzo Nobel chemical problem using band storage.
- `pollu` — 20-variable pollution chemistry problem using sparse KLU. It can run with an analytic sparse Jacobian or sparse difference-quotient Jacobian by passing `use_dq=1`; sparse DQ uses `Radau5SetSparsityPattern` and Curtis-Powell-Reid column grouping.
- `plate` — sparse/KLU plate benchmark.
- `bruss_2d` — 8192-variable 2D Brusselator PDE benchmark using sparse/KLU storage.

### DAE benchmarks

- `chemakzo` — index-1 chemical DAE with dense mass matrix.
- `transamp` — index-1 transistor amplifier with dense mass matrix.
- `transamp_sparse` — sparse-J/sparse-M variant of `transamp` using KLU.
- `andrews` — index-3 Andrews squeezing mechanism; CTest uses reduced-order configurations.
- `caraxis` — index-3 car-axis problem; CTest uses reduced-order configurations.
- `crank` — index-2 slider-crank DAE.
- `water` — index-2 water-tube DAE.
- `fekete` — index-2 constrained dynamics problem.
- `pump` — index-2 pump DAE with analytic Jacobian.
- `pump_smooth` — smoothed pump DAE variant with analytic Jacobian.
- `tba` — 350-variable two-bit adding unit circuit DAE with discontinuities; CTest uses reduced-order configurations.
- `tba_smooth` — smoothed TBA variant; CTest uses reduced-order configurations.
- `tba2` — alternate TBA configuration.

### Rootfinding examples

The rootfinding examples demonstrate `Radau5RootInit`, `Radau5GetRootInfo`, and dense-output based event localization. They stop with `RADAU5_ROOT_RETURN` when a root is found and can be resumed by calling `Radau5Solve` again.

- `bounce` — bouncing ball impact events.
- `roberts_root` — Robertson kinetics with event functions.
- `rocket` — rocket trajectory event detection.
- `orbit` — orbital event detection.
- `kepler` — Kepler problem rootfinding.
- `kneeode` — scalar knee problem event detection.

Usage:

```bash
./bin/radau5_bounce
./bin/radau5_bounce 3 13
```

### ARKODE comparison

- `ark_analytic` — scalar stiff analytic test used for comparison with ARKODE-style examples.

## Notes on order configurations

`ns=3` gives the classic 3-stage order-5 RADAU5 method. Higher stage counts select higher-order Radau IIA formulas: `ns=5` (order 9), `ns=7` (order 13), `ns=9` (order 17), `ns=11` (order 21), and `ns=13` (order 25). Passing different `nsmin` and `nsmax` enables runtime variable-order selection.

Examples:

```bash
./bin/radau5_rober 1e-10 1e-14 1e-12 0 3 3    # fixed classic RADAU5
./bin/radau5_rober 1e-10 1e-14 1e-12 1 5 5    # fixed ns=5, Schur mode
./bin/radau5_rober 1e-10 1e-14 1e-12 0 3 13   # variable order
./bin/radau5_pollu 1e-6 1e-6 1e-6 0 1 3 13    # sparse DQ Jacobian
```
