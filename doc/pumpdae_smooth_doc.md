# 电荷泵 DAE 右端项光滑化方法

## 1. 问题背景

电荷泵（Charge Pump）电路的 DAE 模型包含 9 个方程，其中 MOSFET 的电荷函数 $Q_G$、$Q_S$、$Q_D$ 根据工作区域采用分段定义。这些函数在区域边界上**连续但不可微**（C0 而非 C1），导致：

- Jacobian 矩阵在边界处不连续
- Newton 迭代收敛困难
- 求解器被迫频繁减小步长甚至失败

## 2. 不可微点的定位

通过分析代码和文档，确认了以下不可微点：

| 位置 | 条件 | 性质 | 是否需要处理 |
|------|------|------|:---:|
| `max(ugd - vte, 0)` | `ugd = vte` | C0，导数从 0 跳到 1 | ✅ |
| Region B/C 分支 | `ugs = vte` | C0，导数跳变 | ✅ |
| Region A/B 分支 | `ugb = vfb` | **已经 C1** | ❌ |
| Swap 条件 | `vgs = vgd` | **对称性保证 C1** | ❌ |
| `Vin(t)` 梯形脉冲 | 斜率跳变点 | 已通过分段积分处理 | ❌ |

## 3. 光滑化工具

### 3.1 Softplus 函数（光滑 max）

将 `max(x, 0)` 替换为：

$$\text{softplus}(x, \varepsilon) = \varepsilon \cdot \ln(1 + e^{x/\varepsilon})$$

性质：
- $x \gg \varepsilon$ 时，$\approx x$
- $x \ll -\varepsilon$ 时，$\approx 0$
- 处处 $C^\infty$ 光滑
- 过渡宽度由 $\varepsilon$ 控制

导数为 sigmoid 函数：

$$\frac{d}{dx}\text{softplus}(x, \varepsilon) = \frac{1}{1 + e^{-x/\varepsilon}}$$

### 3.2 Sigmoid 函数（光滑阶跃）

用于混合两个区域的表达式：

$$\sigma(x, \varepsilon) = \frac{1}{1 + e^{-x/\varepsilon}}$$

性质：
- $x \gg \varepsilon$ 时，$\approx 1$
- $x \ll -\varepsilon$ 时，$\approx 0$
- 在 $x = 0$ 处值为 $0.5$，实现平滑过渡

导数：

$$\frac{d}{dx}\sigma(x, \varepsilon) = \frac{\sigma(1-\sigma)}{\varepsilon}$$

### 3.3 数值稳定性

实现中对 $z = x/\varepsilon$ 做了截断保护：

```matlab
function y = softplus(x, eps)
    z = x / eps;
    if z > 30
        y = x;              % 避免 exp 溢出
    elseif z < -30
        y = 0;              % exp(z) ≈ 0
    else
        y = eps * log(1 + exp(z));
    end
end
```

## 4. 具体光滑化操作

### 4.1 光滑化点 1：`max(ugd - vte, 0)`

**原始代码（三个电荷函数中都有）：**
```matlab
ugdt = max(ugd - vte, 0);
```

**光滑化后：**
```matlab
ugdt = softplus(ugd - vte, EPS);
```

同时对 `ugst = ugs - vte` 也做保护（因为 sigmoid 混合时 Region C 表达式在 B 区也会被求值）：
```matlab
ugst = softplus(ugs - vte, EPS);
```

### 4.2 光滑化点 2：Region B/C 边界（`ugs = vte`）

**原始代码（以 `fn_qsrc` 为例）：**
```matlab
if ugs <= vte
    q = 0;           % Region B
else
    q = -COX*(1/3)*H;  % Region C
end
```

**光滑化后：**
```matlab
sigma_BC = sigmoid(ugs - vte, EPS);
q = sigma_BC * q_C;   % 从 0 平滑过渡到 q_C
```

对于 `fn_qgate`，Region B 的值不为零，使用加权混合：
```matlab
sigma_BC = sigmoid(ugs - vte, EPS);
q = (1 - sigma_BC) * q_B + sigma_BC * q_C;
```

### 4.3 图示

```
原始 max(x,0):          光滑化 softplus(x,ε):

     /                        /
    /                        /
   /                      __/
--+-------  x=0      ---~´--------  x=0
  导数跳变 0→1          导数连续变化 0→1
```

```
原始 if/else:            光滑化 sigmoid 混合:

  0 |■■■■■|  q_C         0 |~~~~|  q_C
    |     |                 | /  |
    |     |                 |/   |
----+-----+-- ugs=vte   ---~----+-- ugs=vte
  q=0                     平滑过渡
```

## 5. Jacobian 的对应更新

光滑化后 Jacobian 必须与 RHS 严格一致。关键变化：

### 5.1 `ugdt` 的导数

原始：离散判断 `ugd > vte` 决定导数是 0 还是 1

光滑化后：
```matlab
dsp_d = dsoftplus(ugd - vte, EPS);   % sigmoid，连续值在 [0,1]
dugdt_dugd = dsp_d * 1;
dugdt_dugs = dsp_d * (-dvte_dubs);
dugdt_dugb = dsp_d * dvte_dubs;
```

### 5.2 区域混合的导数

对 $q = \sigma \cdot q_C$，用乘积法则：

$$\frac{dq}{dx} = \frac{d\sigma}{dx} \cdot q_C + \sigma \cdot \frac{dq_C}{dx}$$

对 $q = (1-\sigma) \cdot q_B + \sigma \cdot q_C$：

$$\frac{dq}{dx} = \frac{d\sigma}{dx} \cdot (q_C - q_B) + (1-\sigma) \cdot \frac{dq_B}{dx} + \sigma \cdot \frac{dq_C}{dx}$$

其中 $\sigma$ 对各变量的导数通过链式法则传递：
```matlab
dsigma_dugs = dsigmoid(ugs - vte, EPS) * (1 - dvte_dubs);
dsigma_dugb = dsigmoid(ugs - vte, EPS) * dvte_dubs;
dsigma_dugd = 0;   % vte 不依赖 ugd
```

## 6. 参数选择

光滑参数 $\varepsilon$ 控制过渡带宽度：

| $\varepsilon$ | 过渡带宽度 | 特点 |
|---|---|---|
| `1e-3` | ~几 mV | 求解最容易，但物理精度略有损失 |
| `1e-4` | ~0.1 mV | **推荐值**，平衡精度与光滑性 |
| `1e-5` | ~0.01 mV | 接近原始模型，光滑效果有限 |

选择依据：问题中电压变化量级为 0.01~20V，阈值电压 $V_{T0} = 0.2$V。$\varepsilon = 10^{-4}$ 相对于这些尺度足够小（不影响物理行为），但足以消除 Jacobian 的跳变。

## 7. 预期效果

| 指标 | 原始版本 | 光滑化版本（预期） |
|------|---------|-----------------|
| Newton 迭代收敛 | 边界附近可能失败 | 稳定收敛 |
| 步长拒绝率 | 较高 | 降低 |
| 总步数 | 较多 | 减少 |
| 精度 | 参考解 | 与参考解偏差 < $O(\varepsilon)$ |

## 8. 使用方法

```matlab
% 运行光滑化版本
pumpdae_smooth

% 修改光滑参数（文件第 29 行）
EPS_SMOOTH = 1e-4;  % 默认值，可调整
```

光滑化**不能替代**分段积分：`Vin(t)` 的梯形脉冲在时间上的不连续仍然需要在跳变点重启求解器。光滑化解决的是**空间变量**（电压）方向上电荷函数的不可微性。
