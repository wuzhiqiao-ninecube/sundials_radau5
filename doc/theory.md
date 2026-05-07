# RADAU5: Mathematical Theory and Algorithmic Description

This document describes the mathematical foundations and algorithmic details of the
RADAU5 solver — a 3-stage, order-5 implicit Runge-Kutta method of the Radau IIA
family. The implementation follows the Fortran code of Hairer and Wanner
\[1\] and uses SUNDIALS abstractions for linear algebra.

---

## 1. Problem Formulation

RADAU5 solves initial value problems of the form

$$
M \, y'(t) = f(t, y), \qquad y(t_0) = y_0,
$$

where $y \in \mathbb{R}^n$, $f: \mathbb{R} \times \mathbb{R}^n \to \mathbb{R}^n$,
and $M \in \mathbb{R}^{n \times n}$ is a constant mass matrix provided by the user
via `Radau5SetMassFn`.

Three cases are supported:

- **Explicit ODE** ($M = I$): The identity mass matrix is the default. No mass
  function is needed, and the problem reduces to $y' = f(t, y)$.

- **Implicit ODE** ($M$ nonsingular): The mass matrix $M$ is nonsingular but not
  the identity. This arises in method-of-lines discretizations where the spatial
  mass matrix is non-diagonal (e.g., finite element methods).

- **Differential-algebraic equation** ($M$ singular): When $M$ is singular, the
  system is a DAE. RADAU5 supports DAEs of index up to 3, with the differential
  index communicated via `Radau5SetDAEIndex`. See [Section 8](#8-dae-support).

The mass matrix $M$ may be dense or banded. It is evaluated once (or whenever the
user function is called) and held constant throughout the integration.

---

## 2. The Radau IIA Method

RADAU5 implements the 3-stage Radau IIA collocation method, which is an implicit
Runge-Kutta method of order 5 with stage order 3. It is L-stable and therefore
well suited for stiff problems.

### 2.1 Collocation Nodes

The three collocation nodes are the zeros of

$$
\frac{d^2}{ds^2}\left[ s^2 (s - 1) \right] = 0 \quad \text{and} \quad s = 1,
$$

giving:

$$
c_1 = \frac{4 - \sqrt{6}}{10} \approx 0.15505, \qquad
c_2 = \frac{4 + \sqrt{6}}{10} \approx 0.64495, \qquad
c_3 = 1.
$$

The fact that $c_3 = 1$ is a defining property of Radau methods and means the
last stage coincides with the step endpoint.

### 2.2 Butcher Tableau

The method is defined by the Butcher tableau $(A, b, c)$ where
$A \in \mathbb{R}^{3 \times 3}$, $b \in \mathbb{R}^3$, and
$c = (c_1, c_2, c_3)^T$. The coefficients are:

$$
A = \begin{pmatrix}
\frac{88 - 7\sqrt{6}}{360} & \frac{296 - 169\sqrt{6}}{1800} & \frac{-2 + 3\sqrt{6}}{225} \\[6pt]
\frac{296 + 169\sqrt{6}}{1800} & \frac{88 + 7\sqrt{6}}{360} & \frac{-2 - 3\sqrt{6}}{225} \\[6pt]
\frac{16 - \sqrt{6}}{36} & \frac{16 + \sqrt{6}}{36} & \frac{1}{9}
\end{pmatrix}
$$

Since $c_3 = 1$, the weights satisfy $b^T = A_{3,:}$ (the last row of $A$), so:

$$
b = \left( \frac{16 - \sqrt{6}}{36}, \;\; \frac{16 + \sqrt{6}}{36}, \;\; \frac{1}{9} \right).
$$

### 2.3 Stage Equations

Given the current solution $y_n$ at time $t_n$ and step size $h$, the method
computes stage increments $Z_1, Z_2, Z_3 \in \mathbb{R}^n$ satisfying the
coupled nonlinear system:

$$
M \, Z_i = h \sum_{j=1}^{3} a_{ij} \, f\!\left(t_n + c_j h, \; y_n + Z_j\right),
\qquad i = 1, 2, 3.
$$

Once the stage increments are found, the new solution is:

$$
y_{n+1} = y_n + Z_3.
$$

This follows directly from $c_3 = 1$ and $b = A_{3,:}$, which implies
$y_{n+1} = y_n + h \sum_j b_j k_j = y_n + Z_3$ where $k_j$ are the
standard Runge-Kutta slopes and $Z_i = h \sum_j a_{ij} k_j$.

---

## 3. Simplified Newton Iteration

The $3n$-dimensional nonlinear system for $(Z_1, Z_2, Z_3)$ is solved by
simplified Newton iteration. The key idea is to transform the system so that
the $3n \times 3n$ linear system at each Newton step decouples into smaller
systems that can be solved independently.

### 3.1 The W-Transformation

Define the transformed variables $W_i$ via:

$$
W = A^{-1} Z, \qquad \text{i.e.,} \quad W_i = \sum_{j=1}^{3} (A^{-1})_{ij} \, Z_j.
$$

In terms of $W$, the stage equations become:

$$
\frac{1}{h} M \sum_{j=1}^{3} a_{ij} W_j = f\!\left(t_n + c_i h, \; y_n + Z_i\right),
\qquad i = 1, 2, 3,
$$

and the Newton correction $\Delta W$ satisfies the linear system:

$$
\left( I_3 \otimes \frac{M}{h} - A^{-1} \otimes J \right) \Delta W = \text{RHS},
$$

where $J \approx \partial f / \partial y$ is the Jacobian, evaluated once and
held fixed during the Newton iterations (simplified Newton). Here $I_3$ is the
$3 \times 3$ identity and $\otimes$ denotes the Kronecker product.

The matrix $A^{-1}$ has one real eigenvalue $\gamma$ and a complex conjugate
pair $\alpha \pm i\beta$. By diagonalizing or Schur-decomposing $A^{-1}$, the
$3n \times 3n$ system is reduced to smaller independent systems.

### 3.2 Eigenvalue Decomposition (Default Mode)

The inverse Butcher matrix $A^{-1}$ has eigenvalues:

$$
\gamma, \qquad \alpha + i\beta, \qquad \alpha - i\beta,
$$

and admits the diagonalization $T^{-1} A^{-1} T = \Lambda$ where
$\Lambda = \operatorname{diag}(\gamma, \alpha + i\beta, \alpha - i\beta)$.

The eigenvalues are computed from cube roots following the Fortran code:

$$
\gamma_{\text{raw}} = \frac{6 + 81^{1/3} - 9^{1/3}}{30},
$$

$$
\alpha_{\text{raw}} = \frac{12 - 81^{1/3} + 9^{1/3}}{60}, \qquad
\beta_{\text{raw}} = \frac{(81^{1/3} + 9^{1/3})\sqrt{3}}{60}.
$$

These are then inverted to give the eigenvalues of $A^{-1}$:

$$
\gamma = \frac{1}{\gamma_{\text{raw}}} \approx 3.6378, \qquad
\alpha = \frac{\alpha_{\text{raw}}}{\alpha_{\text{raw}}^2 + \beta_{\text{raw}}^2} \approx 2.6811, \qquad
\beta = \frac{\beta_{\text{raw}}}{\alpha_{\text{raw}}^2 + \beta_{\text{raw}}^2} \approx 2.7202.
$$

Applying the eigenvector transformation $\Delta F = T^{-1} \Delta W$ decouples
the $3n \times 3n$ system into:

1. An $n \times n$ real system:

$$
E_1 \, \delta_1 = r_1, \qquad E_1 = \frac{\gamma}{h} M - J.
$$

2. An $n \times n$ complex system:

$$
\left( \frac{\alpha + i\beta}{h} M - J \right) \delta_{2,3} = r_{2,3}.
$$

Since SUNDIALS has no complex matrix support, the complex system is realified
into a $2n \times 2n$ real system:

$$
E_2 = \begin{pmatrix}
\frac{\alpha}{h} M - J & -\frac{\beta}{h} M \\[4pt]
\frac{\beta}{h} M & \frac{\alpha}{h} M - J
\end{pmatrix}.
$$

Note the antisymmetric structure of the off-diagonal blocks:
the $(1,2)$ block is $-(\beta/h)M$ and the $(2,1)$ block is $+(\beta/h)M$.

The eigenvector matrices $T$ and $T^{-1}$ (denoted `TI` in the code) are
hardcoded constants transcribed from the Fortran source. The forward transform
(physical $\to$ eigenvalue space) uses $T^{-1}$, and the back-transform uses $T$.

### 3.3 Schur Decomposition (Optional Mode)

An alternative to the eigenvalue decomposition is the real Schur decomposition
of $A^{-1}$, enabled via `Radau5SetSchurDecomp(mem, 1)`. This computes an
orthogonal matrix $U_S$ and an upper quasi-triangular matrix $T_S$ such that:

$$
A^{-1} = U_S \, T_S \, U_S^T.
$$

The Schur form $T_S$ has the structure:

$$
T_S = \begin{pmatrix}
T_{S,11} & T_{S,12} & T_{S,13} \\
T_{S,10} & T_{S,11} & T_{S,23} \\
0 & 0 & T_{S,33}
\end{pmatrix}
$$

where the $2 \times 2$ upper-left block
$\bigl(\begin{smallmatrix} T_{S,11} & T_{S,12} \\ T_{S,10} & T_{S,11} \end{smallmatrix}\bigr)$
has eigenvalues $\alpha \pm i\beta$, and $T_{S,33} = \gamma$ is the real
eigenvalue. The precomputed numerical values are:

$$
T_S = \begin{pmatrix}
2.6811 & -8.4239 & -4.0886 \\
1.1046 & 2.6811 & 4.7002 \\
0 & 0 & 3.6378
\end{pmatrix},
\qquad
U_S = \begin{pmatrix}
0.1387 & 0.0463 & 0.9893 \\
-0.2296 & -0.9702 & 0.0776 \\
-0.9633 & 0.2379 & 0.1239
\end{pmatrix}.
$$

In Schur mode, the forward transform uses $U_S^T$ (since $U_S$ is orthogonal,
$U_S^{-1} = U_S^T$) and the back-transform uses $U_S$. The Newton system
becomes block upper triangular rather than diagonal, requiring block
back-substitution instead of independent solves.

The two linear systems in Schur mode are:

1. The $n \times n$ real system (from the $(3,3)$ block of $T_S$):

$$
E_1 = \frac{T_{S,33}}{h} M - J = \frac{\gamma}{h} M - J.
$$

This is identical to the eigenvalue mode $E_1$.

2. The $2n \times 2n$ coupled system (from the $2 \times 2$ Schur block):

$$
E_2 = \begin{pmatrix}
\frac{T_{S,11}}{h} M - J & \frac{T_{S,12}}{h} M \\[4pt]
\frac{T_{S,10}}{h} M & \frac{T_{S,11}}{h} M - J
\end{pmatrix}.
$$

Note that the off-diagonal blocks are **not** antisymmetric: $T_{S,12} \neq -T_{S,10}$
(numerically, $T_{S,12} \approx -8.42$ while $T_{S,10} \approx 1.10$). This is a
key difference from the eigenvalue mode.

The solve proceeds by block back-substitution:

1. Solve $E_1 \, \delta F_3 = r_3$ (the $n \times n$ real system).
2. Substitute $\delta F_3$ into the right-hand sides of the first two equations:
   $r_1 \leftarrow r_1 - (T_{S,13}/h) \, M \, \delta F_3$ and
   $r_2 \leftarrow r_2 - (T_{S,23}/h) \, M \, \delta F_3$.
3. Solve $E_2 \, (\delta F_1; \delta F_2) = (r_1; r_2)$ (the $2n \times 2n$ system).

The Schur mode has the advantage that $U_S$ is orthogonal, which can improve
numerical conditioning of the transforms compared to the potentially
ill-conditioned eigenvector matrix $T$.

### 3.4 Newton Iteration Procedure

Each Newton iteration proceeds as follows:

1. Evaluate the right-hand side $f$ at the three collocation points:
   $f_i = f(t_n + c_i h, \; y_n + Z_i)$ for $i = 1, 2, 3$.

2. Apply the forward transform $T^{-1}$ (eigenvalue mode) or $U_S^T$ (Schur mode)
   to map from physical space to the decoupled space.

3. Form the linear system right-hand sides incorporating the accumulated
   Newton corrections $F_1, F_2, F_3$ (the transformed-space iterates) and
   the mass matrix contributions.

4. Solve the decoupled systems:
   - Eigenvalue mode: solve $E_1$ and $E_2$ independently.
   - Schur mode: solve $E_1$ first, back-substitute, then solve $E_2$.

5. Accumulate the corrections: $F_i \leftarrow F_i + \delta F_i$.

6. Apply the back-transform $T$ (eigenvalue mode) or $U_S$ (Schur mode) to
   recover the physical-space stage increments $Z_1, Z_2, Z_3$.

7. Check convergence.

### 3.5 Convergence Monitoring

The convergence of the simplified Newton iteration is monitored using the
weighted RMS norm of the correction over all $3n$ components:

$$
\|\delta\| = \sqrt{\frac{1}{3n} \sum_{i=1}^{n}
  \left[ \left(\frac{\delta F_{1,i}}{w_i}\right)^2
       + \left(\frac{\delta F_{2,i}}{w_i}\right)^2
       + \left(\frac{\delta F_{3,i}}{w_i}\right)^2 \right]},
$$

where $w_i = \text{atol}_i + \text{rtol}_i \, |y_i|$ is the error weight
(see [Section 6](#6-tolerance-transformation)).

The contraction rate $\theta$ is estimated from successive correction norms.
After the second Newton iteration:

$$
\theta_{\text{new}} =
\begin{cases}
\|\delta^{(k)}\| / \|\delta^{(k-1)}\| & \text{if } k = 2, \\
\sqrt{(\|\delta^{(k)}\| / \|\delta^{(k-1)}\|) \cdot \theta_{\text{prev}}} & \text{if } k > 2.
\end{cases}
$$

The convergence factor is then:

$$
\text{faccon} = \frac{\theta}{1 - \theta}.
$$

The iteration converges when:

$$
\text{faccon} \cdot \|\delta\| \leq \text{fnewt},
$$

where $\text{fnewt}$ is the Newton convergence threshold derived from the
tolerances (see [Section 6](#6-tolerance-transformation)).

If the predicted number of remaining iterations would exceed `nit` (the maximum
allowed), the step size is reduced preemptively. Specifically, the predicted
convergence measure is:

$$
\text{dyth} = \text{faccon} \cdot \|\delta\| \cdot \theta^{(\text{nit} - 1 - k)} / \text{fnewt}.
$$

If $\text{dyth} \geq 1$, the step size is reduced by the factor:

$$
\text{hhfac} = 0.8 \cdot q_{\text{newt}}^{-1/(4 + \text{nit} - 1 - k)},
\qquad q_{\text{newt}} = \max(10^{-4}, \min(20, \text{dyth})),
$$

and the step is retried with the smaller $h$ (return code `RADAU5_NEWT_PREDICT`).

If $\theta \geq 0.99$, the iteration is considered divergent and the step is
rejected with `RADAU5_CONV_FAILURE`.

---

## 4. Error Estimation

After a successful Newton solve, the local error is estimated using the ESTRAD
procedure from the Fortran code. This provides an asymptotically correct
estimate of the local truncation error.

### 4.1 Error Coefficients

The error estimation uses three coefficients derived from the method:

$$
d_1 = -\frac{13 + 7\sqrt{6}}{3}, \qquad
d_2 = \frac{-13 + 7\sqrt{6}}{3}, \qquad
d_3 = -\frac{1}{3}.
$$

### 4.2 First Error Pass

The error estimate is computed in two stages. The first pass forms a linear
combination of the stage increments scaled by $1/h$:

$$
F_{2,i} = \frac{d_1}{h} Z_{1,i} + \frac{d_2}{h} Z_{2,i} + \frac{d_3}{h} Z_{3,i}.
$$

For problems with a general mass matrix $M$, the raw combination is first
computed and then multiplied by $M$:

$$
\text{raw}_i = \frac{d_1}{h} Z_{1,i} + \frac{d_2}{h} Z_{2,i} + \frac{d_3}{h} Z_{3,i},
\qquad F_2 = M \cdot \text{raw}.
$$

The error correction vector is then:

$$
\text{cont}_i = F_{2,i} + f(t_n, y_n)_i,
$$

where $f(t_n, y_n)$ is the right-hand side evaluated at the beginning of the
step (stored as `fn` in the code).

This vector is then solved against the already-factored $E_1$ matrix:

$$
E_1 \cdot \text{cont} = \text{cont} \qquad \text{(in-place solve)}.
$$

The weighted RMS error norm is:

$$
\text{err} = \max\!\left( \sqrt{\frac{1}{n} \sum_{i=1}^{n} \left(\frac{\text{cont}_i}{w_i}\right)^2}, \;\; 10^{-10} \right),
$$

where $w_i$ is the error weight vector.

### 4.3 Second Error Pass

If $\text{err} \geq 1$ and this is either the first step or a step following a
rejection, a second pass is performed to obtain a more accurate error estimate.
This guards against overestimation of the error on difficult initial transients.

The corrected solution is formed:

$$
\tilde{y}_i = y_{n,i} + \text{cont}_i.
$$

The right-hand side is re-evaluated at this corrected point:

$$
\tilde{f} = f(t_n, \tilde{y}),
$$

and the error correction is updated:

$$
\text{cont}_i = \tilde{f}_i + F_{2,i}.
$$

The $E_1$ solve and norm computation are repeated to produce the final error
estimate.

---

## 5. Step Size Control

The step size is adapted after each step attempt based on the error estimate.
RADAU5 uses a combination of a standard error-based controller and an optional
Gustafsson predictive controller.

### 5.1 Basic Step Size Formula

After computing the error estimate $\text{err}$, the proposed new step size is:

$$
h_{\text{new}} = \frac{h}{\text{quot}},
$$

where the quotient $\text{quot}$ is bounded:

$$
\text{quot} = \max\!\left(\text{facr}, \;\; \min\!\left(\text{facl}, \;\; \frac{\text{err}^{1/4}}{\text{fac}}\right)\right),
$$

and the safety factor $\text{fac}$ is:

$$
\text{fac} = \min\!\left(\text{safe}, \;\; \frac{\text{cfac}}{k + 2 \cdot \text{nit}}\right),
$$

with $\text{cfac} = \text{safe} \cdot (1 + 2 \cdot \text{nit})$ and $k$ the
number of Newton iterations actually performed. The exponent $1/4$ reflects the
order $p = 5$ of the method via $1/(p-1+1) = 1/(5-1+1)$... but in practice the
Fortran code uses $1/4$ which corresponds to the embedded error estimator's
effective order.

The bounds $\text{facl}$ and $\text{facr}$ limit how much the step size can
change in a single step:

- $\text{facl} = 5.0$: maximum value of $\text{quot}$, limiting step shrinkage
  to a factor of $1/5$.
- $\text{facr} = 0.125$: minimum value of $\text{quot}$, limiting step growth
  to a factor of $8$.

### 5.2 Gustafsson Predictive Controller

When `pred = 1` (the default), the Gustafsson predictive controller is applied
after step acceptance. This controller uses information from the previous
accepted step to improve the step size prediction:

$$
\text{facgus} = \frac{h_{\text{acc}}}{h} \cdot \frac{(\text{err}^2 / \text{erracc})^{1/4}}{\text{safe}},
$$

where $h_{\text{acc}}$ is the step size of the previous accepted step and
$\text{erracc}$ is the error estimate from the previous accepted step. The
quotient is then updated:

$$
\text{quot} = \max(\text{quot}, \; \text{facgus}),
$$

ensuring the Gustafsson controller can only make the step size more conservative
(smaller) than the basic controller would suggest.

### 5.3 Step Acceptance

When $\text{err} < 1$, the step is accepted. After acceptance:

- The solution is advanced: $y_{n+1} = y_n + Z_3$.
- Continuous output coefficients are updated (see [Section 7](#7-continuous-output)).
- Error weights are recomputed at the new solution.
- The right-hand side $f(t_{n+1}, y_{n+1})$ is evaluated for the next step.

The Jacobian and factorization reuse strategy is then determined:

- If $\theta \leq \text{thet}$ and $q_t = h_{\text{new}}/h \in [\text{quot1}, \text{quot2}]$:
  skip both Jacobian evaluation and matrix factorization on the next step
  (`skipdecomp = 1`).
- If $\theta \leq \text{thet}$ but $q_t$ is outside $[\text{quot1}, \text{quot2}]$:
  reuse the Jacobian but refactor $E_1$ and $E_2$ with the new $h$.
- If $\theta > \text{thet}$: recompute the Jacobian and refactor.

### 5.4 Step Rejection

When $\text{err} \geq 1$, the step is rejected:

- On the first step: $h \leftarrow 0.1 \cdot h$.
- Otherwise: $h \leftarrow h_{\text{new}}$ (from the error formula above).

The step is then retried with the reduced $h$.

### 5.5 Newton Failure Recovery

When the Newton iteration fails to converge or the matrix factorization
encounters a singularity (label 78 in the Fortran code):

$$
h \leftarrow 0.5 \cdot h,
$$

and the step is retried. Up to 5 consecutive failures are allowed before
returning `RADAU5_SINGULAR_MATRIX`.

When Newton predicts slow convergence (`RADAU5_NEWT_PREDICT`), the step size
has already been reduced inside the Newton routine, so the label-78 halving
is skipped.

### 5.6 Default Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `safe` | 0.9 | Safety factor for step size formula |
| `facl` | 5.0 | Upper bound on quot (limits step shrinkage to $1/5$) |
| `facr` | 0.125 | Lower bound on quot (limits step growth to $8\times$) |
| `thet` | 0.001 | Jacobian reuse threshold on $\theta$ |
| `quot1` | 1.0 | Lower bound on $h_{\text{new}}/h$ for skipping decomposition |
| `quot2` | 1.2 | Upper bound on $h_{\text{new}}/h$ for skipping decomposition |
| `nit` | 7 | Maximum Newton iterations per step |
| `pred` | 1 | Step size controller (1 = Gustafsson, 2 = classical) |
| `mxstep` | 100000 | Maximum number of steps |

---

## 6. Tolerance Transformation

RADAU5 follows the Fortran convention of internally transforming the
user-supplied tolerances before use. This transformation is performed once,
on the first call to `Radau5Solve`.

### 6.1 Tolerance Scaling

The user-supplied relative and absolute tolerances $(\text{rtol}, \text{atol})$
are transformed to internal values:

$$
\text{rtol}_{\text{new}} = 0.1 \cdot \text{rtol}^{2/3},
$$

$$
\text{atol}_{\text{new}} = \text{rtol}_{\text{new}} \cdot \frac{\text{atol}}{\text{rtol}}.
$$

This transformation has the effect of tightening the tolerances (since
$\text{rtol}^{2/3} < \text{rtol}$ for $\text{rtol} < 1$) while preserving
the ratio $\text{atol}/\text{rtol}$. It accounts for the fact that the
error estimator of the Radau IIA method tends to overestimate the true error,
so tighter internal tolerances yield solutions that meet the user's requested
accuracy.

For vector tolerances, the transformation is applied component-wise.

### 6.2 Newton Convergence Threshold

The Newton convergence threshold $\text{fnewt}$ is derived from the
(transformed) relative tolerance:

$$
\text{fnewt} = \max\!\left( 10 \cdot u_{\text{round}} / \text{rtol}, \;\;
\min(0.03, \; \sqrt{\text{rtol}}) \right),
$$

where $u_{\text{round}}$ is the unit roundoff. This ensures that Newton
converges to a precision commensurate with the requested accuracy, but never
demands convergence below machine precision.

### 6.3 Error Weight Vector

The error weight vector used in all norm computations is:

$$
w_i = \text{atol}_i + \text{rtol}_i \cdot |y_i|,
$$

where the tolerances are the transformed values. This vector is recomputed
whenever the solution $y$ changes (after each accepted step and at
initialization).

---

## 7. Continuous Output

After each accepted step, RADAU5 constructs a polynomial interpolant that
provides a continuous approximation to the solution between $t_n$ and
$t_{n+1} = t_n + h$. This is used for dense output within the `SolOut`
callback via the function `Radau5Contr(mem, i, t)`.

### 7.1 Interpolation Polynomial

The continuous output is a degree-4 polynomial in the normalized step
coordinate $s = (t - t_n) / h$:

$$
y_i(t_n + s \cdot h) = \text{cont1}_i + s \left( \text{cont2}_i
  + (s - c_{2,m1}) \left( \text{cont3}_i
  + (s - c_{1,m1}) \, \text{cont4}_i \right) \right),
$$

where $c_{1,m1} = c_1 - 1$ and $c_{2,m1} = c_2 - 1$.

### 7.2 Coefficient Computation

The four coefficient vectors are computed from the stage increments
$Z_1, Z_2, Z_3$ after step acceptance. Following the Fortran code
(lines 1017--1027 of `radau5.f`):

$$
\text{cont1}_i = y_{n+1,i},
$$

$$
\text{cont2}_i = \frac{Z_{2,i} - Z_{3,i}}{c_2 - 1},
$$

$$
ak = \frac{Z_{1,i} - Z_{2,i}}{c_1 - c_2}, \qquad
\text{acont3} = \frac{ak - Z_{1,i}/c_1}{c_2},
$$

$$
\text{cont3}_i = \frac{ak - \text{cont2}_i}{c_1 - 1},
$$

$$
\text{cont4}_i = \text{cont3}_i - \text{acont3}.
$$

The polynomial interpolates $y_n$ at $s = 0$ and $y_{n+1}$ at $s = 1$,
and also passes through the internal stage values at $s = c_1$ and $s = c_2$.

### 7.3 Extrapolation for Newton Starting Values

The continuous output polynomial is also used to provide starting values for
the Newton iteration on the next step. Given the new step size $h_{\text{new}}$
and the previous step size $h_{\text{old}}$, the ratio
$c_{3q} = h_{\text{new}} / h_{\text{old}}$ is used to extrapolate:

$$
Z_i^{(0)} = c_{iq} \left( \text{cont2} + (c_{iq} - c_{2,m1})
  \left( \text{cont3} + (c_{iq} - c_{1,m1}) \, \text{cont4} \right) \right),
$$

where $c_{1q} = c_1 \cdot c_{3q}$, $c_{2q} = c_2 \cdot c_{3q}$. On the first
step or when `startn = 1`, zero starting values are used instead.

---

## 8. DAE Support

RADAU5 can solve differential-algebraic equations of the form $M y' = f(t, y)$
where $M$ is singular. The differential index of the DAE is communicated to the
solver via `Radau5SetDAEIndex(mem, nind1, nind2, nind3)`.

### 8.1 Index Classification

The $n$ components of $y$ are partitioned into three groups according to their
differential index:

- **Index-1 variables** (first $n_1$ components): These are the differential
  variables. Standard error weighting is applied.

- **Index-2 variables** (next $n_2$ components): The error weights for these
  components are scaled by the step size ratio:

$$
w_i \leftarrow w_i / h_{\text{fac}}, \qquad i = n_1 + 1, \ldots, n_1 + n_2,
$$

where $h_{\text{fac}}$ tracks the current step size.

- **Index-3 variables** (next $n_3$ components): The error weights are scaled
  by the square of the step size ratio:

$$
w_i \leftarrow w_i / h_{\text{fac}}^2, \qquad i = n_1 + n_2 + 1, \ldots, n_1 + n_2 + n_3.
$$

This scaling reflects the fact that higher-index variables have larger local
errors relative to the step size, and the error weights must be adjusted
accordingly to maintain the correct convergence behavior.

### 8.2 Consistent Initial Conditions

For index-1 DAEs, consistent initial conditions can be computed via
`Radau5CalcIC(mem, id)`, where `id` is a vector identifying the differential
($\text{id}_i = 1$) and algebraic ($\text{id}_i = 0$) components. This
function adjusts the algebraic components so that the initial conditions
satisfy the algebraic constraints.

For index-2 and index-3 DAEs, the user is responsible for providing consistent
initial conditions.

---

## References

1. E. Hairer and G. Wanner, *Solving Ordinary Differential Equations II:
   Stiff and Differential-Algebraic Problems*, Springer Series in Computational
   Mathematics, Vol. 14, 2nd edition, 1996.

2. E. Hairer, S.P. Norsett, and G. Wanner, *Solving Ordinary Differential
   Equations I: Nonstiff Problems*, Springer Series in Computational
   Mathematics, Vol. 8, 2nd edition, 1993.

3. Ch. Lubich and M. Roche, "Rosenbrock methods for differential-algebraic
   systems with solution-dependent singular matrix multiplying the derivative",
   *Computing*, 43(4):325--342, 1990.

4. K. Gustafsson, "Control-theoretic techniques for stepsize selection in
   implicit Runge-Kutta methods", *ACM Trans. Math. Software*, 20(4):496--517,
   1994.

5. Test Set for IVP Solvers, University of Bari,
   [http://www.dm.uniba.it/~testset/](http://www.dm.uniba.it/~testset/).
