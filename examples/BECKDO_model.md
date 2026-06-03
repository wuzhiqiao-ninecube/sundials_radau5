# BECKDO 模型描述（Becker-Döring 方程）

## 英文原文

**BECKDO** — the Becker-Döring model describes the dynamics of a system with a large number of identical particles which can coagulate to form clusters. We let $y_k$ denote the expected number of $k$-particle clusters per unit volume. Assuming that

clusters can gain or loose only single particles, we are led to the system

$$
y_1' = -J_1 - \sum_{k=1}^{N-1} J_k, \quad y_N' = J_{N-1}
$$

$$
y_k' = J_{k-1} - J_k, \quad k = 2, 3, \ldots, N-1,
$$

where $J_k = y_1 y_k - b_{k+1} y_{k+1}$ and $b_{k+1} = \exp(k^2/3 - (k-1)^2/3)$. For a detailed description of this system we refer to the article by Carr, Duncan & Walshaw (1995). This equation is especially interesting because of its **metastability** (extremely slow variations in the solution over very long time intervals; see Fig. 10.6).

As initial condition we take

$$
y_1(0) = \varrho, \quad y_k(0) = 0 \quad \text{for} \quad k = 2, \ldots, N \tag{10.13}
$$

(no clusters at the beginning). It can be seen by differentiation that the density (total number of particles per unit volume)

$$
\sum_{k=1}^N k y_k \quad (= \varrho) \tag{10.14}
$$

is an invariant of the system (10.12). Most numerical schemes (in particular Runge-Kutta methods and multistep methods) preserve automatically such linear invariants in the absence of round-off errors. Whenever the relation (10.14) is not satisfactorily preserved, there is the possibility to re-establish it during the computations by projections (see “differential equations with invariants”, Sect. VII.2). This precautionary measure was not used in the subsequent numerical tests.

In order to be able to observe the metastable states of the system, the dimension $N$ has to be sufficiently large. Following the experiments of Carr, Duncan & Walshaw (1995) we take $N = 5000$ and $\rho = 7.5$, and consider the solution on the interval $0 \leq t \leq 10^{15}$. We compare the errors at $x_{\text{out}} = 1, 10, 10^2, 10^3, \ldots, 10^{15}$.

The Jacobian of this system is tri-diagonal with an additional non-zero first row and a non-zero first column. A Gershgorin test reveals that its eigenvalues can not go, except for the initial phase, beyond $-10$. Stiffness, in this example, is therefore not created by large eigenvalues of $J$, but by the extremely long integration interval.

**BECKDO** — for this problem, the stiff codes (the only ones which work) require the solution of linear systems of the form

$$
\begin{pmatrix}
u & v^T \\
w & T
\end{pmatrix}
\begin{pmatrix}
x \\
y
\end{pmatrix}
=
\begin{pmatrix}
a \\
b
\end{pmatrix},
$$

where $v, w, b$ are $(n-1)$-dimensional vectors and $T$ is a tri-diagonal matrix. Since the linear algebra routines are completely separated from the codes RADAU5, RODAS and SEULEX, it is easy to replace these routines by a special program which solves (10.19) efficiently as follows

$$
x = (a - v^T T^{-1} b) / (u - v^T T^{-1} w)
$$

$$
y = T^{-1} b - x T^{-1} w.
$$

(10.20)

It is not necessary to alter the stiff integrator itself.

Fig. 10.9 shows that, as usual, RODAS is best for low tolerances and RADAU5 is preferable for high precision. Not as usual is the fact that RODAS performs very badly for stringent tolerances. We explain this by the fact that the linear system (10.19) is sensitive to round-off errors, or, as Wilkinson would turn it, delivers a solution for a wrong Jacobian. Thus, the order of the Rosenbrock method drops to 1.

---

## 中文翻译

**BECKDO** — Becker-Döring 模型描述了一个包含大量相同粒子、可以聚集成团簇的系统的动力学。我们用 $y_k$ 表示单位体积内 $k$ 粒子团簇的期望数目。假设团簇只能获得或失去单个粒子，我们得到如下系统：

$$
y_1' = -J_1 - \sum_{k=1}^{N-1} J_k, \quad y_N' = J_{N-1}
$$

$$
y_k' = J_{k-1} - J_k, \quad k = 2, 3, \ldots, N-1,
$$

其中 $J_k = y_1 y_k - b_{k+1} y_{k+1}$，$b_{k+1} = \exp(k^2/3 - (k-1)^2/3)$。关于该系统的详细描述，请参阅 Carr、Duncan 和 Walshaw（1995）的文章。该方程因其**亚稳定性**（解在极长的时间区间内变化极其缓慢；见图 10.6）而特别有趣。

初始条件取为

$$
y_1(0) = \varrho, \quad y_k(0) = 0 \quad \text{对于} \quad k = 2, \ldots, N \tag{10.13}
$$

（开始时没有团簇）。通过微分可以看出，密度（单位体积内的总粒子数）

$$
\sum_{k=1}^N k y_k \quad (= \varrho) \tag{10.14}
$$

是系统 (10.12) 的一个不变量。大多数数值格式（特别是 Runge-Kutta 方法和多步方法）在没有舍入误差的情况下会自动保持这类线性不变量。当关系式 (10.14) 未能令人满意地保持时，可以通过投影在计算过程中重新建立它（参见“带不变量的微分方程”，第 VII.2 节）。在随后的数值测试中没有使用这一预防措施。

为了能够观察到系统的亚稳态，维数 $N$ 必须足够大。根据 Carr、Duncan 和 Walshaw（1995）的实验，我们取 $N = 5000$ 和 $\rho = 7.5$，并考虑区间 $0 \leq t \leq 10^{15}$ 上的解。我们在 $x_{\text{out}} = 1, 10, 10^2, 10^3, \ldots, 10^{15}$ 处比较误差。

该系统的 Jacobi 矩阵是三对角的，并多出一个非零的第一行和一个非零的第一列。Gershgorin 圆盘定理表明，除初始阶段外，其特征值不会小于 $-10$。因此，在这个例子中，刚性并非来自 $J$ 的大特征值，而是来自极长的积分区间。

**BECKDO** — 对于这个问题，刚性求解器（唯一能工作的方法）需要求解如下形式的线性系统：

$$
\begin{pmatrix}
u & v^T \\
w & T
\end{pmatrix}
\begin{pmatrix}
x \\
y
\end{pmatrix}
=
\begin{pmatrix}
a \\
b
\end{pmatrix},
$$

其中 $v, w, b$ 是 $(n-1)$ 维向量，$T$ 是三对角矩阵。由于线性代数例程与 RADAU5、RODAS 和 SEULEX 求解器完全分离，我们可以轻松地用专用程序替换这些例程，该程序高效地按如下方式求解 (10.19)：

$$
x = (a - v^T T^{-1} b) / (u - v^T T^{-1} w)
$$

$$
y = T^{-1} b - x T^{-1} w.
$$

(10.20)

无需修改刚性积分器本身。

图 10.9 显示，与通常情况一样，RODAS 在低精度要求下表现最佳，而 RADAU5 在高精度下更优。不同寻常的是，RODAS 在严格公差下表现非常差。我们将此归因于线性系统 (10.19) 对舍入误差敏感，或者用 Wilkinson 的话说，它给出了一个错误 Jacobi 矩阵的解。因此，Rosenbrock 方法的阶数降为 1。