# II.4 Practical Error Estimation and Step Size Selection

Ich glaube indessen, dass ein practischer Rechner sich meistens mit der geringeren Sicherheit begnügen wird, die er aus der Uebereinstimmung seiner Resultate für grössere und kleinere Schritte gewinnt.
(C. Runge 1895)

Even the simplified error estimates of Section II.3, which are content with the leading error term, are of little practical interest, because they require the computation and majorization of several partial derivatives of high orders. But the main advantage of Runge-Kutta methods, compared with Taylor series, is precisely that the computation of derivatives should be no longer necessary. However, since practical error estimates are necessary (on the one hand to ensure that the step sizes $h_i$ are chosen sufficiently small to yield the required precision of the computed results, and on the other hand to ensure that the step sizes are sufficiently large to avoid unnecessary computational work), we shall now discuss alternative methods for error estimates.

The oldest device, used by Runge in his numerical examples, is to repeat the computations with *halved* step sizes and to compare the results: those digits which haven't changed are assumed to be correct (“... woraus ich schliessen zu dürfen glaube ...”).

# II.4 实际误差估计与步长选择

然而，我相信一位实际的计算者通常会满足于从较大和较小步长的结果的一致性中获得的较低置信度。
（C. Runge 1895）

即使第II.3节中只保留主误差项的简化误差估计，也缺乏实际意义，因为它们需要计算并界定多个高阶偏导数。而龙格-库塔方法相较于泰勒级数的主要优势，恰恰在于不再需要计算这些导数。然而，由于实际误差估计是必要的（一方面要确保步长 $h_i$ 选得足够小以达到所需的计算精度，另一方面要确保步长足够大以避免不必要的计算工作量），我们现在将讨论替代的误差估计方法。

龙格在他数值例子中使用的最古老手段是：用*减半*的步长重复计算并比较结果——那些没有变化的数字被认为是正确的（“……由此我相信可以推断……”）。

---

# Richardson Extrapolation

... its usefulness for practical computations can hardly be overestimated.
(G. Birkhoff & G.C. Rota)

The idea of Richardson, announced in his classical paper Richardson (1910) which treats mainly partial differential equations, and explained in full detail in Richardson (1927), is to use more carefully the known behaviour of the error as a function of $h$.

Suppose that, with a given initial value $(x_0, y_0)$ and step size $h$, we compute two steps, using a fixed Runge-Kutta method of order $p$, and obtain the numerical results $y_1$ and $y_2$. We then compute, starting from $(x_0, y_0)$, one big step with step size $2h$ to obtain the solution $w$. The error of $y_1$ is known to be (Theorem 3.2)

$$
e_1 = y(x_0 + h) - y_1 = C \cdot h^{p+1} + \mathcal{O}(h^{p+2}) \tag{4.1}
$$

where $C$ contains the error coefficients of the method and the elementary differentials $F^J(t)(y_0)$ of order $p+1$. The error of $y_2$ is composed of two parts: the

# 理查德森外推

……它在实际计算中的用处怎么强调都不过分。
（G. Birkhoff & G.C. Rota）

理查德森的思想在其主要处理偏微分方程的经典论文 Richardson (1910) 中提出，并在 Richardson (1927) 中详细阐述，其核心是更精细地利用误差作为 $h$ 函数的已知性质。

假设给定初值 $(x_0, y_0)$ 和步长 $h$，使用一个固定的 $p$ 阶龙格-库塔方法计算两步，得到数值结果 $y_1$ 和 $y_2$。然后，从 $(x_0, y_0)$ 出发，用步长 $2h$ 计算一大步，得到解 $w$。已知 $y_1$ 的误差为（定理3.2）

$$
e_1 = y(x_0 + h) - y_1 = C \cdot h^{p+1} + \mathcal{O}(h^{p+2}) \tag{4.1}
$$

其中 $C$ 包含方法的误差系数和 $p+1$ 阶初等微分 $F^J(t)(y_0)$。$y_2$ 的误差由两部分组成：

---

transported error of the first step, which is

$$
\left( I + h \frac{\partial f}{\partial y} + \mathcal{O}(h^2) \right) e_1,
$$

and the local error of the second step, which is the same as (4.1), but with the elementary differentials evaluated at $y_1 = y_0 + \mathcal{O}(h)$. Thus we obtain

$$
e_2 = y(x_0 + 2h) - y_2 = \left( I + \mathcal{O}(h) \right) Ch^{p+1} + \left( C + \mathcal{O}(h) \right) h^{p+1} + \mathcal{O}(h^{p+2}) = 2Ch^{p+1} + \mathcal{O}(h^{p+2}). \tag{4.2}
$$

Similarly to (4.1), we have for the big step

$$
y(x_0 + 2h) - w = C(2h)^{p+1} + \mathcal{O}(h^{p+2}). \tag{4.3}
$$

Neglecting the terms $\mathcal{O}(h^{p+2})$, formulas (4.2) and (4.3) allow us to eliminate the unknown constant $C$ and to “extrapolate” a better value $\hat{y}_2$ for $y(x_0 + 2h)$, for which we obtain:

**Theorem 4.1.** Suppose that $y_2$ is the numerical result of two steps with step size $h$ of a Runge-Kutta method of order $p$, and $w$ is the result of one big step with step size $2h$. Then the error of $y_2$ can be extrapolated as

$$
y(x_0 + 2h) - y_2 = \frac{y_2 - w}{2^p - 1} + \mathcal{O}(h^{p+2}) \tag{4.4}
$$

and

$$
\hat{y}_2 = y_2 + \frac{y_2 - w}{2^p - 1} \tag{4.5}
$$

is an approximation of order $p + 1$ to $y(x_0 + 2h)$.

第一步的传递误差：

$$
\left( I + h \frac{\partial f}{\partial y} + \mathcal{O}(h^2) \right) e_1,
$$

以及第二步的局部误差，它与 (4.1) 形式相同，但初等微分在 $y_1 = y_0 + \mathcal{O}(h)$ 处取值。于是得到

$$
e_2 = y(x_0 + 2h) - y_2 = \left( I + \mathcal{O}(h) \right) Ch^{p+1} + \left( C + \mathcal{O}(h) \right) h^{p+1} + \mathcal{O}(h^{p+2}) = 2Ch^{p+1} + \mathcal{O}(h^{p+2}). \tag{4.2}
$$

类似于 (4.1)，对于大步长有

$$
y(x_0 + 2h) - w = C(2h)^{p+1} + \mathcal{O}(h^{p+2}). \tag{4.3}
$$

忽略 $\mathcal{O}(h^{p+2})$ 项，公式 (4.2) 和 (4.3) 允许我们消去未知常数 $C$，并“外推”出一个更好的 $y(x_0 + 2h)$ 的近似值 $\hat{y}_2$，从而得到：

**定理4.1.** 假设 $y_2$ 是使用 $p$ 阶龙格-库塔方法以步长 $h$ 计算两步得到的数值结果，$w$ 是以步长 $2h$ 计算一大步得到的结果。那么 $y_2$ 的误差可以外推为

$$
y(x_0 + 2h) - y_2 = \frac{y_2 - w}{2^p - 1} + \mathcal{O}(h^{p+2}) \tag{4.4}
$$

并且

$$
\hat{y}_2 = y_2 + \frac{y_2 - w}{2^p - 1} \tag{4.5}
$$

是 $y(x_0 + 2h)$ 的一个 $p+1$ 阶近似。

---

Formula (4.4) is a very simple device to estimate the error of $y_2$ and formula (4.5) allows one to increase the precision by one additional order (“... The better theory of the following sections is complicated, and tends thereby to suggest that the practice may also be complicated; whereas it is really simple.” Richardson).

公式 (4.4) 是一个估计 $y_2$ 误差的非常简单的技巧，而公式 (4.5) 则能将精度提高一阶（“……后面章节中更完善的理论是复杂的，从而容易让人误以为实践也很复杂；但实际上它非常简单。”——理查德森）。


# Embedded Runge-Kutta Formulas

Scraton is right in his criticism of Merson’s process, although Merson did not claim as much for his process as some people expect.
(R. England 1969)

The idea is, rather than using Richardson extrapolation, to construct Runge-Kutta formulas which themselves contain, besides the numerical approximation \( y_1 \), a second approximation \( \hat{y}_1 \). The difference then yields an estimate of the local error for the less precise result and can be used for step size control (see below). Since

# 嵌入式龙格-库塔公式

斯卡拉顿对默森方法的批评是正确的，尽管默森本人并未像某些人期望的那样为自己的方法宣称那么多。
（R. England 1969）

其思想是，不使用理查德森外推，而是构造龙格-库塔公式，使公式本身除了数值近似 $ y_1 $ 之外，还包含第二个近似 $ \hat{y}_1 $。两者的差值即可给出精度较低结果的局部误差估计，并可用于步长控制（见下文）。由于

it is at our disposal at every step, this gives more flexibility to the code and makes step rejections less expensive.

We consider two Runge-Kutta methods (one for $ y_1 $ and one for $ \hat{y}_1 $) such that both use the same function values. We thus have to find a scheme of coefficients (see (1.8')),

| 0          |               |               |            |                    |               |
|------------|---------------|---------------|------------|--------------------|---------------|
| $ c_2 $    | $ a_{21} $    |               |            |                    |               |
| $ c_3 $    | $ a_{31} $    | $ a_{32} $    |            |                    |               |
| $ \vdots $ | $ \vdots $    | $ \vdots $    |            |                    |               |
| $ c_s $    | $ a_{s1} $    | $ a_{s2} $    | $ \cdots $ | $ a_{s,s-1} $      |               |
|            | $ b_1 $       | $ b_2 $       | $ \cdots $ | $ b_{s-1} $        | $ b_s $       |
|            | $ \hat{b}_1 $ | $ \hat{b}_2 $ | $ \cdots $ | $ \hat{b}_{s-1} $  | $ \hat{b}_s $ |

(4.6)


在每个步骤中我们都可以自由使用，这给代码带来了更大的灵活性，并降低了步长被拒绝时的代价。

我们考虑两个龙格-库塔方法（一个用于 $ y_1 $，一个用于 $ \hat{y}_1 $），使得两者使用相同的函数值。因此我们需要找到一个系数表（参见 (1.8')）：

| 0          |               |               |            |                    |               |
|------------|---------------|---------------|------------|--------------------|---------------|
| $ c_2 $    | $ a_{21} $    |               |            |                    |               |
| $ c_3 $    | $ a_{31} $    | $ a_{32} $    |            |                    |               |
| $ \vdots $ | $ \vdots $    | $ \vdots $    |            |                    |               |
| $ c_s $    | $ a_{s1} $    | $ a_{s2} $    | $ \cdots $ | $ a_{s,s-1} $      |               |
|            | $ b_1 $       | $ b_2 $       | $ \cdots $ | $ b_{s-1} $        | $ b_s $       |
|            | $ \hat{b}_1 $ | $ \hat{b}_2 $ | $ \cdots $ | $ \hat{b}_{s-1} $  | $ \hat{b}_s $ |

(4.6)

such that

$$
y_1 = y_0 + h(b_1 k_1 + \ldots + b_s k_s) \tag{4.7}
$$

is of order $ p $, and

$$
\hat{y}_1 = y_0 + h(\hat{b}_1 k_1 + \ldots + \hat{b}_s k_s) \tag{4.7'}
$$

is of order $ \hat{p} $ (usually $ \hat{p} = p - 1 $ or $ \hat{p} = p + 1 $). The approximation $ y_1 $ is used to continue the integration.

From Theorem 2.13, we have to satisfy the conditions

$$
\sum_{j=1}^s b_j \Phi_j(t) = \frac{1}{\gamma(t)} \quad \text{for all trees of order } \leq p, \tag{4.8}
$$

---

使得

$$
y_1 = y_0 + h(b_1 k_1 + \ldots + b_s k_s) \tag{4.7}
$$

具有阶数 $ p $，并且

$$
\hat{y}_1 = y_0 + h(\hat{b}_1 k_1 + \ldots + \hat{b}_s k_s) \tag{4.7'}
$$

具有阶数 $ \hat{p} $（通常 $ \hat{p} = p - 1 $ 或 $ \hat{p} = p + 1 $）。近似解 $ y_1 $ 用于继续积分。

根据定理 2.13，我们需要满足条件

$$
\sum_{j=1}^s b_j \Phi_j(t) = \frac{1}{\gamma(t)} \quad \text{对于所有阶数 } \leq p \text{ 的树}, \tag{4.8}
$$

$$
\sum_{j=1}^s \hat{b}_j \Phi_j(t) = \frac{1}{\gamma(t)} \quad \text{for all trees of order} \leq \hat{p}. \tag{4.8'}
$$

The first methods of this type were proposed by Merson (1957), Ceschino (1962), and Zonneveld (1963). Those of Merson and Zonneveld are given in Tables 4.1 and 4.2. Here, “name $ p(\hat{p}) $” means that the order of $ y_1 $ is $ p $ and the order of the error estimator \( \hat{y}_1 \) is \( \hat{p} \). Merson’s \( \hat{y}_1 \) is of order 5 only for linear equations with constant coefficients; for nonlinear problems it is of order 3. This method works quite well and has been used very often, especially by NAG users. Further embedded methods were then derived by Sarafyan (1966), England (1969), and Fehlberg (1964, 1968, 1969). Let us start with the construction of some low order embedded methods.

**Methods of order 3(2).** It is a simple task to construct embedded formulas of order 3(2) with \( s = 3 \) stages. Just take a 3-stage method of order 3 (Exercise II.1.4) and put $ \hat{b}_3 = 0, \hat{b}_2 = 1/2c_2, \hat{b}_1 = 1 - 1/2c_2 $.

---

$$
\sum_{j=1}^s \hat{b}_j \Phi_j(t) = \frac{1}{\gamma(t)} \quad \text{对于所有阶数} \leq \hat{p} \text{的树}. \tag{4.8'}
$$

这类方法的第一个例子由 Merson (1957)、Ceschino (1962) 和 Zonneveld (1963) 提出。Merson 和 Zonneveld 的方法分别见表 4.1 和表 4.2。这里，“名称 $ p(\hat{p}) $” 表示 $ y_1 $ 的阶数为 $ p $，误差估计量 $ \hat{y}_1 $ 的阶数为 $ \hat{p} $. Merson 的 $ \hat{y}_1 $ 仅对常系数线性方程为 5 阶；对非线性问题则为 3 阶。该方法效果很好且经常被使用，尤其是被 NAG 库的用户。随后，Sarafyan (1966)、England (1969) 以及 Fehlberg (1964, 1968, 1969) 又推导出更多的嵌入式方法。我们先从构造一些低阶嵌入式方法开始。

**3(2) 阶方法。** 构造具有 $ s = 3 $ 级的 3(2) 阶嵌入式公式是一项简单的任务。只需采用一个 3 级 3 阶方法（习题 II.1.4），并令 $ \hat{b}_3 = 0, \hat{b}_2 = 1/(2c_2), \hat{b}_1 = 1 - 1/(2c_2) $.

---

$$
\hat{b}_1 = 2b_1 - 1/6, \quad \hat{b}_2 = 2(1 - c_2)b_2, \quad \hat{b}_3 = 2(1 - c_3)b_3, \quad \hat{b}_4 = 0, \quad \hat{b}_5 = 1/6. \tag{4.9}
$$

---

$$
\hat{b}_1 = 2b_1 - \frac{1}{6}, \quad \hat{b}_2 = 2(1 - c_2)b_2, \quad \hat{b}_3 = 2(1 - c_3)b_3, \quad \hat{b}_4 = 0, \quad \hat{b}_5 = \frac{1}{6}. \tag{4.9}
$$