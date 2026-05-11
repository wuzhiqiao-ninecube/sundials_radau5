# RADAU5 Examples Guide

## Introduction

The RADAU5 solver ships with 14 example programs covering stiff ODEs, large-scale banded/sparse
systems, and differential-algebraic equations (DAEs) of index 1 through 3, plus 5 event detection
(rootfinding) examples. Twelve of the fourteen ODE/DAE
problems come from the **IVP Test Set 2.4** (<http://www.dm.uniba.it/~testset/>); `heat1d` is a
standard 1D heat equation used to exercise the band matrix path, and `tba` (Two Bit Adding Unit)
is a large DAE-1 MOSFET circuit problem with 63 discontinuities.

All examples accept the same positional command-line arguments:

```
./radau5_<prob> [rtol] [atol] [h0] [use_schur]
```

- `rtol` — relative tolerance (default varies per problem)
- `atol` — absolute tolerance (default varies per problem)
- `h0` — initial step size (default varies per problem)
- `use_schur` — `0` for eigenvalue decomposition (default), `1` for Schur decomposition

Build a single example with the shell script:

```bash
bash build.sh radau5_vdpol
./bin/radau5_vdpol
```

Or build all examples and run the test suite via CMake:

```bash
mkdir build && cd build
cmake .. -DSUNDIALS_DIR=/path/to/sundials/install
make -j$(nproc)
ctest          # 28 tests: 14 problems x 2 modes (eigenvalue + Schur)
```

## Summary Table

| Example | Type | n | Matrix | Jacobian | Default rtol / atol / h0 |
|-----------|-------|-----|---------------|----------|--------------------------|
| vdpol | ODE | 2 | dense | analytic | 1e-6 / 1e-6 / 1e-6 |
| rober | ODE | 3 | dense | analytic | 1e-10 / 1e-14 / 1e-12 |
| orego | ODE | 3 | dense | analytic | 1e-6 / 1e-6 / 1e-6 |
| e5 | ODE | 4 | dense | analytic | 1e-4 / 1e-24 / 1e-6 |
| hires | ODE | 8 | dense | DQ | 1e-6 / 1e-6 / 1e-6 |
| ringmod | ODE | 15 | dense | analytic | 1e-7 / 1e-7 / 1e-6 |
| plei | ODE | 28 | dense | analytic | 1e-6 / 1e-6 / 1e-6 |
| heat1d | ODE | 100 | band(1,1) | DQ | 1e-10 / 1e-12 / 1e-6 |
| medakzo | ODE | 400 | band(2,2) | DQ | 1e-7 / 1e-7 / 1e-6 |
| pollu | ODE | 20 | sparse (KLU) | analytic | 1e-6 / 1e-6 / 1e-6 |
| chemakzo | DAE-1 | 6 | dense + mass | analytic | 1e-7 / 1e-7 / 1e-10 |
| transamp | DAE-1 | 8 | dense + mass | analytic | 1e-6 / 1e-6 / 1e-6 |
| andrews | DAE-3 | 27 | dense + mass | analytic | 1e-6 / 1e-6 / 1e-10 |
| tba | DAE-1 | 350 | dense + mass | DQ | 1e-5 / 1e-5 / 4e-5 |

---

## Stiff ODEs

### vdpol — Van der Pol Oscillator

**Source:** `examples/radau5_vdpol.c`
**Classification:** Stiff ODE, n = 2, very stiff ($\varepsilon = 10^{-6}$)

**Equations:**

$$y_1' = y_2$$
$$y_2' = \frac{(1 - y_1^2)\,y_2 - y_1}{\varepsilon}, \quad \varepsilon = 10^{-6}$$

**Initial conditions:** $y(0) = (2,\; 0)$, $\;t \in [0, 2]$

**Matrix:** Dense (2x2), analytic Jacobian.

**Usage:**

```bash
./bin/radau5_vdpol                    # defaults
./bin/radau5_vdpol 1e-8 1e-8 1e-6 1  # tight tol, Schur mode
```

**Reference solution at $t = 2$:**

| Component | Value |
|-----------|-------|
| $y_1$ | $1.7061677321\ldots$ |
| $y_2$ | $-0.8928097010\ldots$ |

---

### rober — Robertson Chemical Kinetics

**Source:** `examples/radau5_rober.c`
**Classification:** Stiff ODE, n = 3, very stiff, long integration interval

**Equations:**

$$y_1' = -0.04\,y_1 + 10^4\,y_2\,y_3$$
$$y_2' = 0.04\,y_1 - 10^4\,y_2\,y_3 - 3 \times 10^7\,y_2^2$$
$$y_3' = 3 \times 10^7\,y_2^2$$

**Initial conditions:** $y(0) = (1,\; 0,\; 0)$, $\;t \in [0, 10^{11}]$

**Matrix:** Dense (3x3), analytic Jacobian.

**Usage:**

```bash
./bin/radau5_rober
./bin/radau5_rober 1e-12 1e-16 1e-14 0
```

**Reference solution at $t = 10^{11}$:**

| Component | Value |
|-----------|-------|
| $y_1$ | $2.0833 \times 10^{-8}$ |
| $y_2$ | $8.3334 \times 10^{-13}$ |
| $y_3$ | $0.99999998\ldots$ |

---

### orego — Oregonator

**Source:** `examples/radau5_orego.c`
**Classification:** Stiff ODE, n = 3, oscillatory chemical kinetics

**Equations:**

$$y_1' = 77.27\bigl(y_2 + y_1(1 - 8.375 \times 10^{-6}\,y_1 - y_2)\bigr)$$
$$y_2' = \frac{y_3 - (1 + y_1)\,y_2}{77.27}$$
$$y_3' = 0.161\,(y_1 - y_3)$$

**Initial conditions:** $y(0) = (1,\; 2,\; 3)$, $\;t \in [0, 360]$

**Matrix:** Dense (3x3), analytic Jacobian.

**Usage:**

```bash
./bin/radau5_orego
./bin/radau5_orego 1e-9 1e-9 1e-6 1
```

**Reference solution at $t = 360$:**

| Component | Value |
|-----------|-------|
| $y_1$ | $1.000815\ldots$ |
| $y_2$ | $1228.178\ldots$ |
| $y_3$ | $132.0555\ldots$ |

---

### e5 — E5 Stiff Detest

**Source:** `examples/radau5_e5.c`
**Classification:** Extremely stiff ODE, n = 4, chemical kinetics, $t_{\text{end}} = 10^{13}$

**Equations:**

Define the production rates:

$$p_1 = A\,y_1, \quad p_2 = B\,y_1\,y_3, \quad p_3 = C_m\,y_2\,y_3, \quad p_4 = C\,y_4$$

where $A = 7.89 \times 10^{-10}$, $B = 1.1 \times 10^{7}$, $C_m = 1.13 \times 10^{9}$, $C = 1.13 \times 10^{3}$.

$$y_1' = -p_1 - p_2$$
$$y_2' = p_1 - p_3$$
$$y_4' = p_2 - p_4$$
$$y_3' = y_2' - y_4'$$

**Initial conditions:** $y(0) = (1.76 \times 10^{-3},\; 0,\; 0,\; 0)$, $\;t \in [0, 10^{13}]$

**Matrix:** Dense (4x4), analytic Jacobian.

**Usage:**

```bash
./bin/radau5_e5
./bin/radau5_e5 1e-6 1e-26 1e-8 1
```

**Reference solution at $t = 10^{13}$:**

All components decay to near zero. $y_2 \approx 8.868 \times 10^{-22}$,
$y_3 \approx 8.855 \times 10^{-22}$, $y_1 \approx 0$, $y_4 = 0$.

---

### hires — HIRES

**Source:** `examples/radau5_hires.c`
**Classification:** Stiff ODE, n = 8, DQ Jacobian

**Equations:**

$$y_1' = -1.71\,y_1 + 0.43\,y_2 + 8.32\,y_3 + 0.0007$$
$$y_2' = 1.71\,y_1 - 8.75\,y_2$$
$$y_3' = -10.03\,y_3 + 0.43\,y_4 + 0.035\,y_5$$
$$y_4' = 8.32\,y_2 + 1.71\,y_3 - 1.12\,y_4$$
$$y_5' = -1.745\,y_5 + 0.43\,(y_6 + y_7)$$
$$y_6' = -280\,y_6\,y_8 + 0.69\,y_4 + 1.71\,y_5 - 0.43\,y_6 + 0.69\,y_7$$
$$y_7' = 280\,y_6\,y_8 - 1.81\,y_7$$
$$y_8' = -y_7'$$

**Initial conditions:** $y(0) = (1,\; 0,\; 0,\; 0,\; 0,\; 0,\; 0,\; 0.0057)$, $\;t \in [0, 321.8122]$

**Matrix:** Dense (8x8), difference-quotient (DQ) Jacobian — no analytic Jacobian provided.

**Usage:**

```bash
./bin/radau5_hires
./bin/radau5_hires 1e-9 1e-9 1e-6 0
```

**Reference solution at $t = 321.8122$:**

| Component | Value |
|-----------|-------|
| $y_1$ | $7.3713 \times 10^{-4}$ |
| $y_2$ | $1.4425 \times 10^{-4}$ |
| $y_3$ | $5.8887 \times 10^{-5}$ |
| $y_4$ | $1.1757 \times 10^{-3}$ |
| $y_5$ | $2.3864 \times 10^{-3}$ |
| $y_6$ | $6.2390 \times 10^{-3}$ |
| $y_7$ | $2.8500 \times 10^{-3}$ |
| $y_8$ | $2.8500 \times 10^{-3}$ |

---

### ringmod — Ring Modulator

**Source:** `examples/radau5_ringmod.c`
**Classification:** Stiff ODE, n = 15, oscillatory electronic circuit

**Equations:**

The system models a ring modulator circuit with capacitors ($C$, $C_s$, $C_p$),
resistors ($R$, $R_p$, $R_i$, $R_c$, $R_{G1}$, $R_{G2}$, $R_{G3}$), and
inductors ($L_h$, $L_{s1}$, $L_{s2}$, $L_{s3}$). Four ideal diodes produce
nonlinear terms via $q(u) = \gamma\,(\exp(\delta\,u) - 1)$ where
$\gamma = 40.67286402 \times 10^{-9}$ and $\delta = 17.7493332$.

Circuit parameters: $C = 1.6 \times 10^{-8}$, $C_s = 2 \times 10^{-12}$,
$C_p = 10^{-8}$, $R = 25000$, $R_p = 50$, $L_h = 4.45$, $L_{s1} = 2 \times 10^{-3}$,
$L_{s2} = 5 \times 10^{-4}$, $L_{s3} = 5 \times 10^{-4}$, $R_{G1} = 36.3$,
$R_{G2} = 17.3$, $R_{G3} = 17.3$, $R_i = 50$, $R_c = 600$.

Input signals: $u_{\text{in1}}(t) = 0.5\sin(2000\pi t)$,
$u_{\text{in2}}(t) = 2\sin(20000\pi t)$.

Diode arguments:

$$u_{D1} = y_3 - y_5 - y_7 - u_{\text{in2}}, \quad u_{D2} = -y_4 + y_6 - y_7 - u_{\text{in2}}$$
$$u_{D3} = y_4 + y_5 + y_7 + u_{\text{in2}}, \quad u_{D4} = -y_3 - y_6 + y_7 + u_{\text{in2}}$$

Representative equations (capacitor voltages):

$$C\,y_1' = y_8 - \tfrac{1}{2}y_{10} + \tfrac{1}{2}y_{11} + y_{14} - y_1/R$$
$$C_s\,y_3' = y_{10} - q(u_{D1}) + q(u_{D4})$$
$$C_p\,y_7' = -y_7/R_p + q(u_{D1}) + q(u_{D2}) - q(u_{D3}) - q(u_{D4})$$

Inductor currents:

$$L_h\,y_8' = -y_1, \quad L_{s2}\,y_{10}' = \tfrac{1}{2}y_1 - y_3 - R_{G2}\,y_{10}$$
$$L_{s1}\,y_{14}' = -y_1 + u_{\text{in1}} - (R_i + R_{G1})\,y_{14}$$

**Initial conditions:** $y(0) = 0$ (all components), $\;t \in [0, 10^{-3}]$

**Matrix:** Dense (15x15), analytic Jacobian. Max steps raised to 500000.

**Usage:**

```bash
./bin/radau5_ringmod
./bin/radau5_ringmod 1e-9 1e-9 1e-8 1
```

**Reference solution at $t = 10^{-3}$:**

| Component | Value |
|-----------|-------|
| $y_1$ | $-0.02339\ldots$ |
| $y_3$ | $0.25830\ldots$ |
| $y_7$ | $0.11068\ldots$ |
| $y_8$ | $2.940 \times 10^{-7}$ |

---

## Non-stiff ODE

### plei — Pleiades 7-Body Problem

**Source:** `examples/radau5_plei.c`
**Classification:** Non-stiff ODE, n = 28, gravitational N-body problem

**Equations:**

Seven bodies with masses $m_j = j$ ($j = 1, \ldots, 7$) interact gravitationally in 2D.
State vector: 7 x-positions, 7 y-positions, 7 x-velocities, 7 y-velocities.

$$x_i' = \dot{x}_i, \quad y_i' = \dot{y}_i$$

$$\ddot{x}_i = \sum_{\substack{j=1 \\ j \neq i}}^{7} \frac{m_j\,(x_j - x_i)}{r_{ij}^3}, \quad
\ddot{y}_i = \sum_{\substack{j=1 \\ j \neq i}}^{7} \frac{m_j\,(y_j - y_i)}{r_{ij}^3}$$

where $r_{ij} = \sqrt{(x_j - x_i)^2 + (y_j - y_i)^2}$.

**Initial conditions:**

| Body | $x_i(0)$ | $y_i(0)$ | $\dot{x}_i(0)$ | $\dot{y}_i(0)$ |
|------|-----------|-----------|-----------------|-----------------|
| 1 | 3 | 3 | 0 | 0 |
| 2 | 3 | -3 | 0 | 0 |
| 3 | -1 | 2 | 0 | 0 |
| 4 | -3 | 0 | 0 | -1.25 |
| 5 | 2 | 0 | 0 | 1 |
| 6 | -2 | -4 | 1.75 | 0 |
| 7 | 2 | 4 | -1.5 | 0 |

$t \in [0, 3]$

**Matrix:** Dense (28x28), analytic Jacobian.

**Usage:**

```bash
./bin/radau5_plei
./bin/radau5_plei 1e-9 1e-9 1e-6 0
```

**Reference solution at $t = 3$ (positions only):**

| Component | Value |
|-----------|-------|
| $x_1$ | $0.3706\ldots$ |
| $x_2$ | $3.2373\ldots$ |
| $y_1$ | $-3.9434\ldots$ |
| $y_2$ | $-3.2714\ldots$ |

---

## Large-Scale ODEs (Band / Sparse)

### heat1d — 1D Heat Equation

**Source:** `examples/radau5_heat1d.c`
**Classification:** Mildly stiff ODE, n = 100, method-of-lines PDE discretization

**Equations:**

The 1D heat equation $u_t = u_{xx}$ on $[0, 1]$ with homogeneous Dirichlet boundary
conditions $u(0,t) = u(1,t) = 0$ is discretized with $n = 100$ interior points
($\Delta x = 1/101$):

$$y_i' = \frac{y_{i-1} - 2\,y_i + y_{i+1}}{\Delta x^2}, \quad i = 1, \ldots, 100$$

with $y_0 = y_{101} = 0$.

**Initial conditions:** $y_i(0) = \sin(\pi\,x_i)$ where $x_i = i\,\Delta x$, $\;t \in [0, 0.1]$

**Exact solution:** $u(x,t) = e^{-\pi^2 t}\sin(\pi x)$

**Matrix:** Band with $m_l = 1$, $m_u = 1$ (tridiagonal), DQ Jacobian.

**Usage:**

```bash
./bin/radau5_heat1d
./bin/radau5_heat1d 1e-12 1e-14 1e-8 1
```

**Reference:** Exact discrete solution. Error dominated by spatial discretization
$O(\Delta x^2) \approx 3 \times 10^{-5}$.

---

### medakzo — Medical Akzo Nobel

**Source:** `examples/radau5_medakzo.c`
**Classification:** Stiff ODE, n = 400, reaction-diffusion with discontinuity

**Equations:**

A reaction-diffusion system on $[0, 1]$ with $N = 200$ grid points, yielding 400 ODEs.
Each grid point $j$ has two unknowns $(y_{2j-1}, y_{2j})$. Let $\zeta_j = j/N$,
$\Delta\zeta = 1/N$:

$$\beta_j = \left(\frac{(\zeta_j - 1)^2}{c}\right)^2, \quad
\alpha_j = \frac{2(\zeta_j - 1)}{c}\,\beta_j^{1/2}$$

where $c = 4$, $k = 100$.

For interior points $j = 2, \ldots, N-1$:

$$y_{2j-1}' = \frac{\beta_j}{\Delta\zeta^2}(y_{2j-3} - 2\,y_{2j-1} + y_{2j+1})
+ \frac{\alpha_j}{2\,\Delta\zeta}(y_{2j+1} - y_{2j-3}) - k\,y_{2j-1}\,y_{2j}$$
$$y_{2j}' = -k\,y_{2j-1}\,y_{2j}$$

Boundary: $\phi(t) = 2$ for $t \le 5$, $\phi(t) = 0$ for $t > 5$ (discontinuity at $t = 5$).

**Initial conditions:** $y_{2j-1}(0) = 0$, $y_{2j}(0) = 1$ for $j = 1, \ldots, 200$, $\;t \in [0, 20]$

**Matrix:** Band with $m_l = 2$, $m_u = 2$ (pentadiagonal), DQ Jacobian.

The solver integrates in two phases (0 to 5, then 5 to 20) to handle the discontinuity.

**Usage:**

```bash
./bin/radau5_medakzo
./bin/radau5_medakzo 1e-9 1e-9 1e-8 0
```

**Reference solution at $t = 20$ (selected components):**

| Component | Value |
|-----------|-------|
| $y_{79}$ | $2.340 \times 10^{-4}$ |
| $y_{133}$ | $3.577 \times 10^{-4}$ |
| $y_{171}$ | $3.086 \times 10^{-4}$ |
| $y_{199}$ | $1.174 \times 10^{-4}$ |
| $y_{200}$ | $6.191 \times 10^{-6}$ |

---

### pollu — Pollution Problem

**Source:** `examples/radau5_pollu.c`
**Classification:** Stiff ODE, n = 20, atmospheric chemistry, sparse Jacobian

**Equations:**

A system of 20 species coupled by 25 reactions with rate constants $k_1, \ldots, k_{25}$.
Define the reaction rates:

$$r_1 = k_1\,y_1, \quad r_2 = k_2\,y_2\,y_4, \quad r_3 = k_3\,y_5\,y_2, \quad r_4 = k_4\,y_7, \quad r_5 = k_5\,y_7$$
$$r_6 = k_6\,y_7\,y_6, \quad r_7 = k_7\,y_9, \quad r_8 = k_8\,y_9\,y_6, \quad r_9 = k_9\,y_{11}\,y_2, \quad r_{10} = k_{10}\,y_{11}\,y_1$$
$$r_{11} = k_{11}\,y_{13}, \quad r_{12} = k_{12}\,y_{10}\,y_2, \quad r_{13} = k_{13}\,y_{14}, \quad r_{14} = k_{14}\,y_1\,y_6, \quad r_{15} = k_{15}\,y_3$$
$$r_{16} = k_{16}\,y_4, \quad r_{17} = k_{17}\,y_4, \quad r_{18} = k_{18}\,y_{16}, \quad r_{19} = k_{19}\,y_{16}, \quad r_{20} = k_{20}\,y_{17}\,y_6$$
$$r_{21} = k_{21}\,y_{19}, \quad r_{22} = k_{22}\,y_{19}, \quad r_{23} = k_{23}\,y_1\,y_4, \quad r_{24} = k_{24}\,y_{19}\,y_1, \quad r_{25} = k_{25}\,y_{20}$$

Rate constants:

| $k_1$ | $k_2$ | $k_3$ | $k_4$ | $k_5$ | $k_6$ | $k_7$ | $k_8$ | $k_9$ |
|--------|--------|--------|--------|--------|--------|--------|--------|--------|
| 0.35 | 26.6 | 1.23e4 | 8.6e-4 | 8.2e-4 | 1.5e4 | 1.3e-4 | 2.4e4 | 1.65e4 |

| $k_{10}$ | $k_{11}$ | $k_{12}$ | $k_{13}$ | $k_{14}$ | $k_{15}$ | $k_{16}$ | $k_{17}$ |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| 9.0e3 | 0.022 | 1.2e4 | 1.88 | 1.63e4 | 4.8e6 | 3.5e-4 | 0.0175 |

| $k_{18}$ | $k_{19}$ | $k_{20}$ | $k_{21}$ | $k_{22}$ | $k_{23}$ | $k_{24}$ | $k_{25}$ |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| 1.0e8 | 4.44e11 | 1.24e3 | 2.1 | 5.78 | 0.0474 | 1.78e3 | 3.12 |

Representative ODEs:

$$y_1' = -r_1 - r_{10} - r_{14} - r_{23} - r_{24} + r_2 + r_3 + r_9 + r_{11} + r_{12} + r_{22} + r_{25}$$
$$y_2' = -r_2 - r_3 - r_9 - r_{12} + r_1 + r_{21}$$

**Initial conditions:** $y_2(0) = 0.2$, $y_4(0) = 0.04$, $y_7(0) = 0.1$, $y_8(0) = 0.3$,
$y_9(0) = 0.01$, $y_{17}(0) = 0.007$, all others zero. $\;t \in [0, 60]$

**Matrix:** Sparse CSC (KLU), 86 nonzeros, analytic Jacobian. DQ Jacobian is not
supported for sparse matrices — an analytic Jacobian is required.

**Usage:**

```bash
./bin/radau5_pollu
./bin/radau5_pollu 1e-8 1e-8 1e-6 1
```

**Reference solution at $t = 60$ (selected components):**

| Component | Value |
|-----------|-------|
| $y_1$ | $5.646 \times 10^{-2}$ |
| $y_2$ | $1.342 \times 10^{-1}$ |
| $y_8$ | $3.245 \times 10^{-1}$ |
| $y_{15}$ | $8.965 \times 10^{-3}$ |

---

## DAEs with Mass Matrix

### chemakzo — Chemical Akzo Nobel

**Source:** `examples/radau5_chemakzo.c`
**Classification:** DAE index-1, n = 6, chemical kinetics with algebraic constraint

**Equations:**

Five differential equations plus one algebraic constraint. Define the reaction rates:

$$r_1 = k_1\,y_1^4\,\sqrt{|y_2|}, \quad r_2 = k_2\,y_3\,y_4, \quad r_3 = \frac{k_2}{K}\,y_1\,y_5$$
$$r_4 = k_3\,y_1\,y_4^2, \quad r_5 = k_4\,y_6^2\,\sqrt{|y_2|}$$
$$F_{\text{in}} = k_{La}\left(\frac{p_{O_2}}{H} - y_2\right)$$

where $k_1 = 18.7$, $k_2 = 0.58$, $k_3 = 0.09$, $k_4 = 0.42$, $K = 34.4$,
$k_{La} = 3.3$, $k_s = 115.83$, $p_{O_2} = 0.9$, $H = 737$.

$$y_1' = -2\,r_1 + r_2 - r_3 - r_4$$
$$y_2' = -\tfrac{1}{2}\,r_1 - r_4 - \tfrac{1}{2}\,r_5 + F_{\text{in}}$$
$$y_3' = r_1 - r_2 + r_3$$
$$y_4' = -r_2 + r_3 - 2\,r_4$$
$$y_5' = r_2 - r_3 + r_5$$
$$0 = k_s\,y_1\,y_4 - y_6 \quad \text{(algebraic)}$$

**Mass matrix:** $M = \operatorname{diag}(1, 1, 1, 1, 1, 0)$

**Initial conditions:** $y(0) = (0.444,\; 0.00123,\; 0,\; 0.007,\; 0,\; k_s \cdot 0.444 \cdot 0.007)$,
$\;t \in [0, 180]$

**Matrix:** Dense (6x6), analytic Jacobian. DAE index set via `Radau5SetDAEIndex(mem, 6, 0, 0)`.

**Usage:**

```bash
./bin/radau5_chemakzo
./bin/radau5_chemakzo 1e-9 1e-9 1e-12 1
```

**Reference solution at $t = 180$:**

| Component | Value |
|-----------|-------|
| $y_1$ | $0.11508\ldots$ |
| $y_2$ | $1.2038 \times 10^{-3}$ |
| $y_3$ | $0.16116\ldots$ |
| $y_4$ | $3.6562 \times 10^{-4}$ |
| $y_5$ | $1.7080 \times 10^{-2}$ |
| $y_6$ | $4.8735 \times 10^{-3}$ |

---

### transamp — Transistor Amplifier

**Source:** `examples/radau5_transamp.c`
**Classification:** DAE index-1, n = 8, oscillatory electronic circuit with banded mass matrix

**Equations:**

A 6-transistor amplifier circuit. The system is $M\,y' = f(t, y)$ where $M$ is a
banded mass matrix (bandwidth 1). Circuit parameters:

$$U_b = 6, \quad U_F = 0.026, \quad \alpha = 0.99, \quad \beta = 10^{-6}$$
$$R_0 = 1000, \quad R_k = 9000 \;(k = 1, \ldots, 9)$$
$$C_1 = 10^{-6}, \quad C_2 = 2 \times 10^{-6}, \quad C_3 = 3 \times 10^{-6}, \quad C_4 = 4 \times 10^{-6}, \quad C_5 = 5 \times 10^{-6}$$

Input signal: $U_e(t) = 0.1\sin(200\pi t)$.

Nonlinear diode terms:

$$g_1 = \beta\left(\exp\!\left(\frac{y_2 - y_3}{U_F}\right) - 1\right), \quad
g_2 = \beta\left(\exp\!\left(\frac{y_5 - y_6}{U_F}\right) - 1\right)$$

The RHS:

$$f_1 = \frac{y_1 - U_e(t)}{R_0}$$
$$f_2 = \frac{y_2}{R_1} + \frac{y_2 - U_b}{R_2} + (1 - \alpha)\,g_1$$
$$f_3 = \frac{y_3}{R_3} - g_1$$
$$f_4 = \frac{y_4 - U_b}{R_4} + \alpha\,g_1$$
$$f_5 = \frac{y_5}{R_5} + \frac{y_5 - U_b}{R_6} + (1 - \alpha)\,g_2$$
$$f_6 = \frac{y_6}{R_7} - g_2$$
$$f_7 = \frac{y_7 - U_b}{R_8} + \alpha\,g_2$$
$$f_8 = \frac{y_8}{R_9}$$

**Mass matrix:** Banded (tridiagonal), stored dense:

$$M = \begin{pmatrix}
-C_1 & C_1 & & & & & & \\
C_1 & -C_1 & & & & & & \\
& & -C_2 & & & & & \\
& & & -C_3 & C_3 & & & \\
& & & C_3 & -C_3 & & & \\
& & & & & -C_4 & & \\
& & & & & & -C_5 & C_5 \\
& & & & & & C_5 & -C_5
\end{pmatrix}$$

**Initial conditions:** $y(0) = (0,\; 3,\; 3,\; 6,\; 3,\; 3,\; 6,\; 0)$, $\;t \in [0, 0.2]$

(Derived from steady-state: $y_2 = y_3 = U_b/(R_2/R_1 + 1) = 3$, $y_4 = y_7 = U_b = 6$.)

**Matrix:** Dense (8x8), analytic Jacobian. DAE index: `Radau5SetDAEIndex(mem, 8, 0, 0)`.

**Usage:**

```bash
./bin/radau5_transamp
./bin/radau5_transamp 1e-8 1e-8 1e-8 1
```

**Reference solution at $t = 0.2$:**

| Component | Value |
|-----------|-------|
| $y_1$ | $-5.562 \times 10^{-3}$ |
| $y_2$ | $3.0065\ldots$ |
| $y_3$ | $2.8500\ldots$ |
| $y_4$ | $2.9264\ldots$ |
| $y_5$ | $2.7046\ldots$ |
| $y_6$ | $2.7618\ldots$ |
| $y_7$ | $4.7709\ldots$ |
| $y_8$ | $1.2370\ldots$ |

---

### andrews — Andrews' Squeezer Mechanism

**Source:** `examples/radau5_andrews.c`
**Classification:** DAE index-3, n = 27, constrained multibody dynamics

**Equations:**

A planar mechanism with 7 rigid bodies connected by joints, subject to 6 holonomic
constraints. The state vector contains 7 angles $(\beta, \theta, \gamma, \phi, \delta,
\omega, \epsilon)$, 7 angular velocities, 7 "accelerations" (from the second-order
formulation), and 6 Lagrange multipliers $\lambda_1, \ldots, \lambda_6$.

The system is $M\,y' = f(t, y)$ where:

- $y_{1\text{--}7}$: angles
- $y_{8\text{--}14}$: angular velocities ($y_{i+7} = y_i'$)
- $y_{15\text{--}21}$: from the inertia equation $\mathbf{M}_{\text{inertia}}\,\ddot{q} = F_v - G^T\lambda$
- $y_{22\text{--}27}$: Lagrange multipliers (algebraic)

The inertia matrix $\mathbf{M}_{\text{inertia}}$ is a 7x7 symmetric matrix depending on
the angles, with entries built from body masses ($m_1, \ldots, m_7$), moments of inertia
($I_1, \ldots, I_7$), and geometric parameters. A spring force with stiffness $c_0 = 4530$
and rest length $l_0 = 0.07785$ acts on body 3.

Physical parameters (selected):

| Parameter | Value | Parameter | Value |
|-----------|-------|-----------|-------|
| $m_1$ | 0.04325 | $I_1$ | $2.194 \times 10^{-6}$ |
| $m_2$ | 0.00365 | $I_2$ | $4.410 \times 10^{-7}$ |
| $m_7$ | 0.05498 | $I_7$ | $1.912 \times 10^{-5}$ |
| $c_0$ | 4530 | $l_0$ | 0.07785 |
| $d$ | 0.028 | $e$ | 0.02 |

The 6 position-level constraints $g(q) = 0$ enforce the closed kinematic loops:

$$g_1 = r_r\cos\beta - d\cos(\beta+\theta) - s_s\sin\gamma - x_B$$
$$g_2 = r_r\sin\beta - d\sin(\beta+\theta) + s_s\cos\gamma - y_B$$
$$g_3 = r_r\cos\beta - d\cos(\beta+\theta) - e\sin(\phi+\delta) - z_t\cos\delta - x_A$$
$$g_4 = r_r\sin\beta - d\sin(\beta+\theta) + e\cos(\phi+\delta) - z_t\sin\delta - y_A$$
$$g_5 = r_r\cos\beta - d\cos(\beta+\theta) - z_f\cos(\omega+\epsilon) - u\sin\epsilon - x_A$$
$$g_6 = r_r\sin\beta - d\sin(\beta+\theta) - z_f\sin(\omega+\epsilon) + u\cos\epsilon - y_A$$

The assembled RHS:

$$f_{1\text{--}14}: \quad y_{i+7}' = y_{i+14} \quad \text{(kinematic identity)}$$
$$f_{15\text{--}21}: \quad \mathbf{M}_{\text{inertia}}\,y_{15\text{--}21} + G^T\lambda = F_v$$
$$f_{22\text{--}27}: \quad g(q) = 0 \quad \text{(constraints)}$$

**Mass matrix:** $M = \operatorname{diag}(\underbrace{1, \ldots, 1}_{14},\; \underbrace{0, \ldots, 0}_{13})$

**DAE index:** `Radau5SetDAEIndex(mem, 7, 7, 13)` — 7 index-1, 7 index-2, 13 index-3 variables.

**Initial conditions:** Consistent initial state with nonzero angles, two nonzero
accelerations ($y_{15} = 14222.4\ldots$, $y_{16} = -10666.8\ldots$), and two nonzero
multipliers ($y_{22} = 98.567\ldots$, $y_{23} = -6.123\ldots$). $\;t \in [0, 0.03]$

**Matrix:** Dense (27x27), analytic Jacobian (approximate, from Hairer & Wanner p.540).

**Usage:**

```bash
./bin/radau5_andrews
./bin/radau5_andrews 1e-8 1e-8 1e-12 0
```

**Reference solution at $t = 0.03$ (angles only):**

| Component | Value |
|-----------|-------|
| $\beta$ ($y_1$) | $15.811\ldots$ |
| $\theta$ ($y_2$) | $-15.756\ldots$ |
| $\gamma$ ($y_3$) | $0.04082\ldots$ |
| $\phi$ ($y_4$) | $-0.53473\ldots$ |
| $\delta$ ($y_5$) | $0.52441\ldots$ |
| $\omega$ ($y_6$) | $0.53473\ldots$ |
| $\epsilon$ ($y_7$) | $1.04808\ldots$ |

---

## DAE with Discontinuities

### tba — Two Bit Adding Unit

**Source:** `examples/radau5_tba.c`
**Classification:** DAE index-1, n = 350, MOSFET circuit with 63 discontinuities

**Description:**

MOSFET circuit simulation of a two-bit adder computing:
$$A_1 \cdot 2 + A_0 + B_1 \cdot 2 + B_0 + C_{in} = C \cdot 4 + S_1 \cdot 2 + S_0$$

The system is formulated as a semi-explicit index-1 DAE:
$$M \cdot y' = f(t, y)$$

where $M = \text{diag}(1,\ldots,1,0,\ldots,0)$ with 175 ones and 175 zeros.

**State vector (0-based):**
- $y[0..174]$ = charges $g(U)$
- $y[175..349]$ = node voltages $U$

**Circuit structure:**
- 10 logic gates: 3 NOR, 4 ANDOI, 1 NAND, 1 ORANI, 1 ANDOIP
- 3 additional enhancement transistors in series
- Shichman-Hodges MOSFET model

**Integration domain:** $t \in [0, 320]$ with 63 discontinuities at $t = 5k$, $k = 1,\ldots,63$.
Five PULSE input signals have derivative jumps at multiples of 5. The solver integrates
segment-by-segment and calls `Radau5ResetForDiscontinuity` at each boundary.

**Matrix:** Dense (350×350), DQ Jacobian, diagonal mass matrix.

**Usage:**

```bash
./bin/radau5_tba
./bin/radau5_tba 1e-5 1e-5 4e-5 1
```

**Reference solution at $t = 315$ (from IVP Test Set, radau5 with rtol=atol=1e-5):**

| Component | Node | Signal | Value |
|-----------|------|--------|-------|
| $y[223]$ | 49 | $S_0$ | $0.2040419147\ldots$ |
| $y[304]$ | 130 | $S_1$ | $4.9972384557\ldots$ |
| $y[322]$ | 148 | $C$ | $0.2038985905\ldots$ |

---

## Event Detection (Rootfinding) Examples

The following examples demonstrate the rootfinding capability of RADAU5. They use
`Radau5RootInit` to register event functions and handle `RADAU5_ROOT_RETURN` in the
integration loop. These examples do NOT accept the standard `[rtol atol h0 use_schur]`
command-line arguments — they use hardcoded parameters appropriate for each problem.

### bounce — Bouncing Ball

**Source:** `examples/radau5_bounce.c`
**Ported from:** MATLAB `ballode.m`

**Equations:**

$$y_1' = y_2 \quad \text{(height)}$$
$$y_2' = -g \quad \text{(velocity, } g = 9.8\text{)}$$

**Initial conditions:** $y(0) = (0,\; 20)$, ball launched upward.

**Event:** $g_0 = y_1$ (height = 0), direction = $-1$ (falling only). Terminal.

**Behavior:** On each bounce, velocity is reversed with coefficient of restitution 0.9.
The solver is reinitialized after each event via `Radau5Init`. Runs for 10 bounces
or until $t = 30$.

**Usage:**

```bash
./bin/radau5_bounce
```

**Analytical first bounce time:** $t = 2v_0/g = 40/9.8 \approx 4.0816$

---

### roberts_root — Robertson Chemical Kinetics with Events

**Source:** `examples/radau5_roberts_root.c`
**Ported from:** SUNDIALS `cvRoberts_dns.c`

**Equations:**

$$y_1' = -0.04\,y_1 + 10^4\,y_2\,y_3$$
$$y_2' = 0.04\,y_1 - 10^4\,y_2\,y_3 - 3\times10^7\,y_2^2$$
$$y_3' = 3\times10^7\,y_2^2$$

**Initial conditions:** $y(0) = (1, 0, 0)$, $\;t \in [0, 4\times10^{10}]$

**Events (2 root functions, non-terminal):**
- $g_0 = y_1 - 10^{-4}$ (species 1 drops to $10^{-4}$)
- $g_1 = y_3 - 0.01$ (species 3 reaches 0.01)

**Behavior:** Since RADAU5 rootfinding is always terminal, the example resumes
integration after each root return to simulate non-terminal behavior. Two roots
are detected: $g_1$ fires near $t \approx 0.264$ and $g_0$ fires near $t \approx 2.08\times10^7$.

**Usage:**

```bash
./bin/radau5_roberts_root
```

---

### rocket — Rocket Ascent with Multi-Phase Events

**Source:** `examples/radau5_rocket.c`
**Ported from:** SUNDIALS `cvRocket_dns.c`

**Equations:**

$$H' = v \quad \text{(height)}$$
$$v' = a(t,v) \quad \text{(velocity)}$$

where the acceleration depends on engine state:
- Engine ON: $a = F/(M_r + M_{f0} - r\,t) - D\,v - g$
- Engine OFF: $a = -D\,v - g$

Parameters: $F = 2200$, $M_r = 10$, $M_{f0} = 1$, $r = 0.1$, $D = 0.3$, $g = 32$, $H_{cut} = 4000$.

**Initial conditions:** $y(0) = (0, 0)$

**Events (multi-phase):**
- Phase 1 (engine on, 2 root functions):
  - $g_0 = M_{f0} - r\,t$ (fuel exhaustion)
  - $g_1 = H - H_{cut}$ (height cutoff)
- Phase 2 (engine off, 1 root function):
  - $g_0 = v$ with direction $= -1$ (maximum height when velocity crosses zero falling)

**Behavior:** Demonstrates terminal events with solver restart and dynamic switching
of the root function set between phases. After engine cutoff, the solver is
reinitialized with a new root function. After max height is detected, rootfinding
is disabled for the descent phase.

**Usage:**

```bash
./bin/radau5_rocket
```

**Key events:** Engine cutoff at $t \approx 10.0$ (fuel exhaustion at $H \approx 4000$ ft),
maximum height $\approx 5262$ ft at $t \approx 16.2$.

---

### orbit — Restricted Three-Body Problem with Direction-Sensitive Events

**Source:** `examples/radau5_orbit.c`
**Ported from:** MATLAB `orbitode.m` (Shampine & Gordon, p.246)

**Equations (rotating frame):**

$$x' = v_x$$
$$y' = v_y$$
$$v_x' = 2\,v_y + x - \mu^*(x+\mu)/r_{13} - \mu(x-\mu^*)/r_{23}$$
$$v_y' = -2\,v_x + y - \mu^*\,y/r_{13} - \mu\,y/r_{23}$$

where $\mu = 1/82.45$, $\mu^* = 1 - \mu$, $r_{13} = ((x+\mu)^2 + y^2)^{3/2}$,
$r_{23} = ((x-\mu^*)^2 + y^2)^{3/2}$.

**Initial conditions:** $y(0) = (1.2,\; 0,\; 0,\; -1.04935750983031990726)$, $\;t \in [0, 7]$

**Events (2 root functions, same value, different directions):**

$$\frac{d}{dt}\|y - y_0\|^2 = 2\,((x - x_0)\,v_x + (y - y_0)\,v_y)$$

- $g_0$: direction $= +1$ → local minimum of distance (orbit returns to start). **Terminal.**
- $g_1$: direction $= -1$ → local maximum of distance (farthest point). Non-terminal.

**Behavior:** This is the canonical test for direction-sensitive event detection.
Both events share the same function value; only the crossing direction distinguishes
them. The orbit completes one period near $t \approx 6.19$ with the body returning
to within $\sim 2\times10^{-4}$ of the initial point.

**Usage:**

```bash
./bin/radau5_orbit
```

---

### kepler — Kepler Two-Body Orbit with Half-Orbit Counting

**Source:** `examples/radau5_kepler.c`
**Ported from:** SUNDIALS `ark_kepler.c`

**Equations (Hamiltonian):**

$$q_1' = p_1, \quad q_2' = p_2$$
$$p_1' = -q_1/(q_1^2 + q_2^2)^{3/2}, \quad p_2' = -q_2/(q_1^2 + q_2^2)^{3/2}$$

**Parameters:** Eccentricity $e = 0.6$, orbital period $T = 2\pi$.

**Initial conditions (perihelion):**
$q_1(0) = 1-e = 0.4$, $q_2(0) = 0$, $p_1(0) = 0$, $p_2(0) = \sqrt{(1+e)/(1-e)} = 2.0$

**Event:** $g_0 = q_2$ (y-coordinate crosses zero → half-orbit). No direction filter.
Non-terminal.

**Behavior:** Counts x-axis crossings over multiple orbits. Each crossing represents
a half-orbit. Also monitors energy conservation ($H = \frac{1}{2}(p_1^2+p_2^2) - 1/r$).

**Usage:**

```bash
./bin/radau5_kepler
```

**Expected:** For 16 orbits ($t_{final} = 32\pi$), approximately 31 half-orbit crossings
with energy drift depending on tolerance.

---

### kneeode — Knee Problem with Non-Negativity via Rootfinding

**Source:** `examples/radau5_kneeode.c`
**Ported from:** MATLAB `kneeode.m` (Dahlquist, Edsberg, Skollermo, Soderlind)

**Equations:**

$$\varepsilon\,y' = (1-x)\,y - y^2, \quad \varepsilon = 10^{-6}$$

**Initial conditions:** $y(0) = 1$, $\;x \in [0, 2]$

**Behavior:** The solution follows the isocline $y = 1-x$ for $x < 1$, then drops
to $y = 0$ for $x > 1$. Without non-negativity enforcement, numerical solvers
overshoot and produce incorrect negative solutions near the "knee" at $x = 1$.

**Event:** $g_0 = y$ (detect $y = 0$), direction $= -1$ (falling only). Terminal.

After the event fires (near $x \approx 1.006$), the solution is held at $y = 0$
for the remainder of the integration. This works because $f(t, 0) = 0$ — zero is
a stable equilibrium for $x > 1$.

**Key technique:** This example demonstrates using rootfinding as a workaround for
the lack of built-in `NonNegative` constraints in RADAU5. The approach is applicable
whenever the zero boundary is a stable equilibrium of the ODE.

**Usage:**

```bash
./bin/radau5_kneeode
```

The example runs twice (with and without the constraint) and prints a comparison table.
