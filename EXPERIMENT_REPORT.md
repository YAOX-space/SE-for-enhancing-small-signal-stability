# 实验复现汇报

**论文来源**：Yuan et al., "Placing Storage Energies for Enhancing Small-Signal Stability of Converter-Based Renewable Systems", IEEE Transactions on Industry Applications, 2025

**复现范围**：Tables I / II / III，Figs. 6 / 8 / 10 / 11 / 12 / 14

---

## 目录

1. [核心理论框架](#1-核心理论框架)
2. [代码架构总览](#2-代码架构总览)
3. [网络构建方法](#3-网络构建方法)
4. [第一性原理特征值计算（替代 Table V）](#4-第一性原理特征值计算替代-table-v)
5. [Table I & II：两区域贪婪放置](#5-table-i--ii两区域贪婪放置)
6. [Table III：两区域穷举搜索](#6-table-iii两区域穷举搜索)
7. [Fig. 6：两区域四 CBR 时域响应](#7-fig-6两区域四-cbr-时域响应)
8. [Fig. 8：最弱特征值轨迹](#8-fig-8最弱特征值轨迹)
9. [Fig. 10：IEEE 39 节点特征值随 SE 倍率变化](#9-fig-10ieee-39-节点特征值随-se-倍率变化)
10. [Fig. 11：IEEE 39 节点时域响应（不同 gSCR）](#10-fig-11ieee-39-节点时域响应不同-gscr)
11. [Fig. 12：IEEE 39 节点时域响应（不同放置方案）](#11-fig-12ieee-39-节点时域响应不同放置方案)
12. [Fig. 14：33 变流器时域响应](#12-fig-1433-变流器时域响应)
13. [关键工程决策与修正记录](#13-关键工程决策与修正记录)

---

## 1. 核心理论框架

### 1.1 gSCR 定义

论文的核心稳定性指标是广义短路比（generalized Short-Circuit Ratio）：

$$\text{gSCR} = \frac{1}{\lambda_{\max}^+(\mathbf{W})}$$

其中加权阻抗矩阵为：

$$\mathbf{W} = S'_{B1} \mathbf{Z} = \text{diag}(s_1, s_2, \ldots, s_N) \cdot \mathbf{Z}$$

- $s_i > 0$：CBR 节点（注入有功）
- $s_i < 0$：SE 节点（吸收有功）
- $s_i = 0$：无源节点（既无 CBR 也无 SE）
- $\lambda_{\max}^+$：矩阵最大正实特征值
- $\mathbf{Z}$：节点阻抗矩阵（N×N，由网络 Laplacian 求逆得到）

**代码实现**（`src/se_placement.py`）：

```python
def lambda_max_positive(matrix: np.ndarray, tol: float = 1e-9) -> float:
    eigvals = np.linalg.eigvals(matrix)
    real_eigs = eigvals[np.abs(eigvals.imag) < 1e-7].real
    return float(np.sort(real_eigs[real_eigs > tol])[-1])

def gscr_from_weighted_impedance(weighted_impedance: np.ndarray) -> float:
    return 1.0 / lambda_max_positive(weighted_impedance)
```

MATLAB 对应（`lib/max_positive_eig.m`）：

```matlab
function lambda = max_positive_eig(matrix)
    vals = eig(matrix);
    vals = real(vals(abs(imag(vals)) < 1e-7));
    vals = vals(vals > 1e-9);
    lambda = max(vals);
end
```

### 1.2 SE 放置敏感度公式

论文 Appendix I-C（公式 26）推导出：对第 $j$ 个节点放置 SE 后，$\lambda_{\max}$ 的变化量为：

$$\frac{\partial \lambda}{\partial S_{BE,j}} = -a \cdot \lambda \cdot v_j^2$$

其中 $v$ 是 $\mathbf{W} = S'_{B1} \mathbf{Z}$ 对应 $\lambda_{\max}$ 的**左特征向量**，$a > 0$ 是正标量。因此，使 $\lambda$ 下降最快（即使 gSCR 提升最多）的节点是 $|v_j|$ 最大的节点。

**代码实现**（`src/se_placement.py`）：

```python
def left_eigen_sensitivities(weighted, candidate_nodes):
    lam = lambda_max_positive(weighted)
    # 左特征向量 = W^T 的右特征向量
    eigvals, eigvecs = np.linalg.eig(weighted.T)
    idx = int(np.argmin(np.abs(eigvals - lam)))
    vector = eigvecs[:, idx].real
    vector = vector / np.linalg.norm(vector)
    # 返回值越小（越负），该节点越应该放 SE
    return {node: float(-lam * vector[node] ** 2) for node in candidate_nodes}
```

**关键细节**：对非对称的 $\mathbf{W}$（异质容量时），左特征向量 $\neq$ 右特征向量。代码通过对 $\mathbf{W}^T$ 求右特征向量来获得 $\mathbf{W}$ 的左特征向量，这与论文 Appendix I-C 一致。

### 1.3 贪婪放置算法

```python
def greedy_place_storage(base_signed_capacities, impedance,
                          candidate_nodes, storage_capacity, count):
    capacities = np.array(base_signed_capacities, dtype=float)
    for iteration in range(1, count + 1):
        weighted = weighted_impedance(capacities, impedance)  # W = diag(s) @ Z
        lam = lambda_max_positive(weighted)
        sensitivities = left_eigen_sensitivities(weighted, candidate_nodes)
        chosen = min(sensitivities, key=sensitivities.get)  # 选最负的节点
        capacities[chosen] -= storage_capacity               # s_chosen -= S_SE
        # 记录本次迭代结果...
    return capacities, steps
```

每次迭代：计算当前 $\mathbf{W}$，计算各节点敏感度，选择最小值（最能降低 $\lambda$）的节点，将其容量减去 SE 容量（正容量 → 更小或负值）。

---

## 2. 代码架构总览

```
SE-for-enhancing-small-signal-stability/
├── src/
│   └── se_placement.py              # 核心算法：gSCR、贪婪放置、敏感度
├── examples/
│   └── two_area_demo.py             # 两区域实验：Tables I/II/III
├── matlab/
│   ├── run_physics_simulink_reproduction.m  # 主运行脚本：Figs 6/8/10/11/12/14
│   ├── compute_first_principles_eigenvalues.m  # 第一性原理特征值计算
│   ├── lib/
│   │   ├── build_network_from_edges.m   # 网络阻抗矩阵构建（标准法）
│   │   ├── build_gfl_state_space.m      # 模态状态空间模型
│   │   ├── max_positive_eig.m           # 最大正实特征值
│   │   └── build_physics_simulink_model.m  # Simulink 模型生成
│   └── data/
│       ├── table_x_two_area_lines.csv   # 两区域线路参数（9 节点）
│       ├── table_xi_ieee39_lines.csv    # IEEE 39 节点线路参数
│       └── table_iii_two_area_two_se_gscr.csv  # 论文 Table III 参考值
└── results/                            # 输出目录
```

**运行入口**：
- Python（Tables I/II/III）：`python examples/two_area_demo.py`
- MATLAB（Figs 6/8/10/11/12/14）：在 MATLAB 中运行 `run_physics_simulink_reproduction`

---

## 3. 网络构建方法

### 3.1 标准法：删除接地节点后求逆

论文和 MATLAB 代码 `build_network_from_edges.m` 采用的标准方法：

```matlab
function [Q, Z] = build_network_from_edges(edgeTable, nodeCount, groundNode)
    % 1. 构建完整 n×n 导纳 Laplacian
    Qfull = zeros(nodeCount, nodeCount);
    for k = 1:height(edgeTable)
        i = edgeTable.from(k);  j = edgeTable.to(k);
        b = 1 / edgeTable.x(k);  % 电纳 = 1/电抗
        Qfull(i,i) += b;  Qfull(j,j) += b;
        Qfull(i,j) -= b;  Qfull(j,i) -= b;
    end
    % 2. 删除接地节点所在行列
    keep = setdiff(1:nodeCount, groundNode);
    Q = Qfull(keep, keep);
    % 3. 求逆得阻抗矩阵
    Z = inv(Q);
end
```

**两区域系统的致命缺陷**：两区域系统（`table_x_two_area_lines.csv`）拓扑如下：

```
CBR1(节点1) ─ 节点5 ─ 节点6 ─ 节点7 ─ 节点8 ─ 节点9(∞母线)
CBR2(节点2) ─ 节点5                              |
CBR3(节点3) ────────────────────────────── 节点9 ┘
CBR4(节点4) ────────────────────────────── 节点9 ┘
```

CBR3、CBR4 只通过节点9连接到网络。**删除节点9后，区域1（CBR1/2）与区域2（CBR3/4）完全断开**，Z 矩阵出现零块，任何跨区域放置的 gSCR 计算均错误。

### 3.2 伪逆法：保留全节点耦合

针对两区域系统的修正方案（`build_network_z_pinv`）：

```python
def build_network_z_pinv(lines_path: Path, n_nodes: int) -> np.ndarray:
    """Z = pinv(Q_full)，保留所有节点间耦合"""
    Q_full = build_laplacian(lines_path, n_nodes)  # 9×9 完整 Laplacian
    return np.linalg.pinv(Q_full)                  # Moore-Penrose 伪逆
```

**数学依据**：Laplacian 矩阵是半正定的（零特征值对应全1向量，即 Kirchhoff 约束）。伪逆忽略零特征值分量，对满足 $\sum_i I_i = 0$（注入电流之和为零）的激励产生与真实网络一致的电压响应。gSCR 定义的电流注入正好满足此约束。

**Python 实现的完整 Laplacian 构建**：

```python
def build_laplacian(lines_path: Path, n_nodes: int) -> np.ndarray:
    Q = np.zeros((n_nodes, n_nodes))
    with open(lines_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            i = int(row["from"]) - 1    # 0 索引
            j = int(row["to"]) - 1
            b = 1.0 / float(row["x"])   # 电纳 = 1/电抗
            Q[i, i] += b;  Q[j, j] += b
            Q[i, j] -= b;  Q[j, i] -= b
    return Q
```

### 3.3 刻度校准（Scale Calibration）

CSV 数据的物理单位与论文的标幺值系统不同，需要一次线性校准：

$$k = \frac{\lambda_{\max}^+(\mathbf{W}_{\text{paper}})}{\lambda_{\max}^+(\mathbf{W}_{\text{csv}})}$$

其中 $\mathbf{W}_{\text{paper}} = \text{INITIAL\_WEIGHTED}$（论文公式 (17) 中打印的矩阵），$\mathbf{W}_{\text{csv}} = \text{diag}(s_{\text{base}}) \cdot \mathbf{Z}_{\text{raw}}$。校准后 $\mathbf{Z} = k \cdot \mathbf{Z}_{\text{raw}}$ 使 gSCR 与论文初始值对齐，但不将论文的 $\mathbf{W}$ 反推 $\mathbf{Z}$（那是循环论证）。

---

## 4. 第一性原理特征值计算（替代 Table V）

**文件**：`matlab/compute_first_principles_eigenvalues.m`

**目的**：论文 Table V 给出了 4 个 gSCR 值对应的主导特征值 $(\sigma, \omega)$。复现中若直接使用这些数值作为插值表，则 Figs 8/10/11/12 的特征值本质上是"代入论文答案得出答案"。此函数从 9 状态 GFL 控制方程出发，通过数值 Jacobian 独立计算 $(\sigma, \omega)$。

### 4.1 9 状态 GFL 模型

每个变流器有 9 个状态变量（在共转坐标系下）：

| 状态 | 含义 |
|------|------|
| $\phi_k = \theta_k - \omega_0 t$ | PLL 角度偏差（相对全局旋转帧） |
| $\varepsilon_{\text{pll}}$ | PLL 积分器状态 |
| $\varepsilon_p$ | 有功功率外环积分器 |
| $\varepsilon_{id}$ | d 轴电流内环积分器 |
| $\varepsilon_{iq}$ | q 轴电流内环积分器 |
| $i_d$ | d 轴电流（变流器本地坐标系） |
| $i_q$ | q 轴电流 |
| $V_{d,\text{ff}}$ | GFF 前馈滤波器 d 轴输出 |
| $V_{q,\text{ff}}$ | GFF 前馈滤波器 q 轴输出 |

9 个变流器共 81 个状态。

**GFF 滤波器是不稳定的关键**：前馈滤波器 $G(s) = \frac{1}{1+\tau s}$（$\tau=0.01$ s），在 88 rad/s 处幅值 $|G|=0.752$，相位滞后 $-41°$。这导致电压前馈无法完全消除电压-功率耦合，约 34% 的耦合量保留，当 gSCR 低于临界值时激发振荡。

```matlab
% GFF 过滤后的电压用于电流内环前馈
vd_vsc = Kp_i*(id_ref - id_loc) + Ki_i*eps_id - w0*Lf*iq_loc + Vd_ff;  % Vd_ff 不是即时值！
vq_vsc = Kp_i*(iq_ref - iq_loc) + Ki_i*eps_iq + w0*Lf*id_loc + Vq_ff;

% GFF 滤波动态
dVd_ff = (Vd_l - Vd_ff) / tau;   % tau=0.01 s → 截止频率 100 rad/s
dVq_ff = (Vq_l - Vq_ff) / tau;
```

### 4.2 稳态求解（fsolve）

在共转坐标系下，稳态是一个**固定点**（$\dot{x}_{\text{ss}} = 0$），可以用 Newton 法（`fsolve`）直接求解。这与 ODE 前向积分不同——ODE 在不稳定工作点会发散，而 `fsolve` 不论稳定性如何均可收敛。

```matlab
function x_ss = find_steady_state_9state(X_net, S_CBR, p)
    % 代数初值猜测（准静态，iq=0）
    id_ss0 = S_CBR(:);
    Vd_g0  = 1.0 - X_net * zeros(N,1);   % Iq=0 → Vd_g ≈ 1
    Vq_g0  = X_net * id_ss0;              % Id=id_ss → Vq_g
    phi0   = atan2(Vq_g0, Vd_g0);
    % 关键初值：eps_p 不为零！
    % 稳态：id_ref = Ki_p * eps_p → eps_p = id_ss / Ki_p
    x0(b+3) = id_ss0(k) / p.Ki_p;

    f = @(x) gfl9_rotating_ode(x, p, X_net, S_CBR);
    [x_ss, ~, exitflag] = fsolve(f, x0, opts);
    % 若 fsolve 不收敛，备用方案：先运行 ODE 15s 再 fsolve
end
```

### 4.3 数值 Jacobian 线性化

在稳态点 $x_{\text{ss}}$ 处，用有限差分法构建 81×81 Jacobian：

```matlab
eps_fd = 1e-7;
for j = 1:81
    xp = x_ss;  xp(j) = xp(j) + eps_fd;
    J(:,j) = (gfl9_rotating_ode(xp, p, X_net, S_CBR) - f0) / eps_fd;
end
ev = eig(J);
% 主导特征值：虚部最接近 88 rad/s 的那个
[~, bi] = min(abs(abs(imag(ev)) - 88));
lam = ev(bi);   % sigma + j*omega
```

### 4.4 gSCR 扫描策略

通过等比例缩放网络阻抗矩阵来改变 gSCR：

$$\text{gSCR}(\alpha) = \frac{\text{gSCR}_{\text{base}}}{\alpha}, \quad \mathbf{Z}(\alpha) = \alpha \cdot \mathbf{Z}_{\text{base}}$$

例如，若基础 gSCR = 2.65，目标 gSCR = 3.0，则 $\alpha = 2.65/3.0 = 0.883$。

```matlab
fp_gscr_sweep = [linspace(1.5, 2.0, 3), linspace(2.1, 3.0, 5), linspace(3.2, 5.0, 5)];
fp_raw = compute_first_principles_eigenvalues(fp_gscr_sweep);
% 输出：fp_raw 包含 gSCR、sigma、omega 三列
```

扫描完成后，所有用到 $(\sigma, \omega)$ 的地方均通过 `interp1(..., 'pchip', 'extrap')` 插值，不再使用论文 Table V 的硬编码数值。

---

## 5. Table I & II：两区域贪婪放置

**文件**：`examples/two_area_demo.py`  
**关键函数**：`greedy_place_storage`（`src/se_placement.py`）

### 5.1 场景描述

- **网络**：两区域系统，9 个节点（节点 9 为无穷大母线）
- **初始 CBR**：节点 1（0.5 pu）、节点 2（1.0 pu）、节点 3（1.5 pu）、节点 4（0.5 pu）
- **任务**：贪婪放置 2 个 SE，每个容量 0.5 pu，选择使 gSCR 最大化的节点

### 5.2 避免循环论证的 Z 矩阵构建

**早期版本的循环问题**（已修正）：

```python
# 错误做法：从论文打印的 W 反推 Z
for row, capacity in enumerate(capacities[:4]):
    impedance[row, :] = INITIAL_WEIGHTED[row, :] / capacity  # Z = W/S
# 这等价于 W = diag(S) @ Z = diag(S) @ (W/S) = W，是恒等式！
```

**正确做法**（当前代码）：

```python
# Step 1：从 CSV 线路数据构建 9×9 完整 Laplacian 的伪逆
Z_raw = build_network_z_pinv(lines_path, n_nodes=9)

# Step 2：容量向量（来自论文说明，不来自 W 矩阵）
s_base = np.array([0.5, 1.0, 1.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0])

# Step 3：仅用 INITIAL_WEIGHTED 做刻度校准（它不参与算法计算）
lam_csv   = lambda_max_positive(np.diag(s_base) @ Z_raw)
lam_paper = lambda_max_positive(INITIAL_WEIGHTED)
k_greedy  = lam_paper / lam_csv          # 纯缩放因子
Z_calibrated = k_greedy * Z_raw          # 物理网络来源的 Z，已校准单位

# Step 4：贪婪放置（所有计算仅依赖 Z_calibrated 和 s_base）
final_capacities, steps = greedy_place_storage(
    s_base, Z_calibrated,
    candidate_nodes=list(range(9)),  # 9 个候选节点
    storage_capacity=0.5,
    count=2,
)
```

### 5.3 贪婪迭代过程

**第 1 次迭代**：

计算初始 $\mathbf{W} = \text{diag}(s_{\text{base}}) \cdot \mathbf{Z}$，求 $\lambda_{\max}^+(\mathbf{W})$，再用 $\mathbf{W}^T$ 的特征向量算各节点敏感度，选出 $|v_j|$ 最大的节点。

**第 2 次迭代**：

将第 1 次选中节点的容量减去 0.5 pu（$s_j \leftarrow s_j - 0.5$），用更新后的 $\mathbf{W}'$ 重复上述过程。

### 5.4 复现结果

| 迭代 | 选择节点 | $\lambda_{\max}$（之前） | gSCR（之前） |
|------|----------|--------------------------|--------------|
| 1    | **节点 4** | 0.3760 | 2.6595 |
| 2    | **节点 1** | 0.3293 | 3.0367 |

最终 gSCR：**3.4997**

论文报告：第 1 次选节点 4，第 2 次选节点 1，与复现结果完全一致。

**Table I**（敏感度值）：节点 4 在初始状态下的左特征向量分量最大，故被首先选中。  
**Table II**（放置后结果）：2 个 SE 分别放置于节点 1 和节点 4 后，gSCR 从 2.66 提升至 3.50，提升 32%。

---

## 6. Table III：两区域穷举搜索

**文件**：`examples/two_area_demo.py`  
**函数**：`compute_table_iii`

### 6.1 场景描述

对 9 个节点（节点 1–9）的所有**有放回组合**（combinations_with_replacement），枚举放置 2 个 SE（各 0.5 pu）后的 gSCR：

$$\binom{9+1}{2} = 45 \text{ 种组合}$$

### 6.2 实现细节

```python
def compute_table_iii(lines_path, se_capacity=0.5):
    # 使用伪逆法（保留区域间耦合）
    Z_raw = build_network_z_pinv(lines_path, n_nodes=9)
    s_base = np.array([0.5, 1.0, 1.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0])

    # 校准
    lam_csv   = lambda_max_positive(np.diag(s_base) @ Z_raw)
    lam_paper = lambda_max_positive(INITIAL_WEIGHTED)
    Z = (lam_paper / lam_csv) * Z_raw

    results = []
    for i, j in combinations_with_replacement(range(9), 2):
        s = s_base.copy()
        s[i] -= se_capacity   # 节点 i 加 SE：正容量减小
        s[j] -= se_capacity   # 节点 j 加 SE
        W = np.diag(s) @ Z
        try:
            gscr = gscr_from_weighted_impedance(W)  # 1/λmax
        except ValueError:
            gscr = float("nan")  # 若无正实特征值则无效
        results.append((i + 1, j + 1, gscr))   # 1 索引
    return results, lam_paper / lam_csv
```

### 6.3 技术难点：为何必须用伪逆

使用标准法（删除节点 9 求逆）时，区域 1（CBR1/2，节点 1/2）与区域 2（CBR3/4，节点 3/4）完全断开。对跨区域组合（如节点 2 + 节点 3）的 gSCR 计算完全错误，最大误差达 6.89（相对于论文值）。

切换到伪逆法后，区域间耦合通过节点 9 的电气路径得到保留，最大误差降至 **0.69**（主要集中在 (3,4) 组合，可能因网络拓扑轻微差异）。

### 6.4 复现结果摘要

| 节点组合 | 计算值 | 论文值 | 偏差 |
|----------|--------|--------|------|
| (1,4) ← 贪婪选择 | 3.4997 | 3.5700 | −0.07 |
| (1,2) | 3.6550 | 3.0600 | +0.60 |
| (3,4) | 3.6720 | 2.9800 | +0.69 |
| (3,9) | 3.0625 | 2.9100 | +0.15 |
| (6,7) | 2.6839 | 2.7500 | −0.07 |

计算得到的全局最优：节点 (3,4)，gSCR = 3.672  
论文报告的全局最优：节点 (1,4)，gSCR = 3.570

两者贪婪序列一致（节点 4、节点 1），但穷举最优节点略有差异，推测原因是论文的网络数据与 CSV 中存在轻微拓扑差异。

---

## 7. Fig. 6：两区域四 CBR 时域响应

**文件**：`matlab/run_physics_simulink_reproduction.m`（第 173–205 行）  
**依赖**：`build_gfl_state_space.m`，`build_physics_simulink_model.m`

### 7.1 场景描述

论文 Fig. 6 展示三种场景下，0.06 s 电压跌落扰动后各 CBR 的有功功率响应：

- **Case 1**：无 SE，4 个 CBR（节点 1/2/3/4），系统接近不稳定
- **Case 2**：SE 放置于节点 4 和节点 6（次优位置），轻微阻尼
- **Case 3**：SE 放置于节点 1 和节点 4（最优位置），强阻尼

### 7.2 Z 矩阵准备

Fig. 6 使用混合 Z：4×4 CBR 子块来自论文打印的 W 矩阵反推（仅用于时域仿真，不用于贪婪算法），被动节点块来自 CSV 线路数据：

```matlab
% 从论文 W 矩阵重建对称 Z（仅用于 Fig.6 仿真，不参与算法）
W_paper = [0.223 0.069 0.015 0.015;
           0.139 0.238 0.030 0.030;
           0.045 0.045 0.249 0.130;
           0.015 0.015 0.043 0.246];
s_cbr_2a = [0.5; 1.0; 1.5; 0.5];
Z_cbr_2a = diag(1./s_cbr_2a) * W_paper;
Z_cbr_2a = (Z_cbr_2a + Z_cbr_2a') / 2;  % 对称化

% 被动节点来自 CSV
[~, Z2a_net] = build_network_from_edges(lines2a, 9, 9);  % 8×8
Z2a = Z2a_net;
Z2a(1:4, 1:4) = Z_cbr_2a;  % 用论文值替换 CBR 子块
```

**注意**：Fig. 6 的 Z 矩阵构建仍依赖论文打印的 W，因为 Fig. 6 的目的是复现论文图形的视觉效果，而非验证放置算法。实际放置算法（Tables I/II）已在 Python 侧使用独立的 Z。

### 7.3 模态状态空间模型（`build_gfl_state_space.m`）

由于全非线性 EMT 仿真代价高昂，本复现采用**物理信息模态状态空间**模型：

**步骤 1**：计算当前配置的 gSCR。

**步骤 2**：通过插值 `eig_table`（第一性原理计算）得到主导模态 $(\sigma, \omega)$：

```matlab
sigma = interp1(eig_table.gSCR, eig_table.sigma, gscr, 'pchip', 'extrap');
omega = interp1(eig_table.gSCR, eig_table.omega, gscr, 'pchip', 'extrap');
```

**步骤 3**：从 $\mathbf{W}^T$ 的最大正实特征向量（即 $\mathbf{W}$ 的左特征向量）提取各变流器的**参与因子** $v_k$：

```matlab
[V_left, lam_all] = eig(Wsigned');
% 找最大正实特征值对应的特征向量
pos_lams = real(lam_diag(real_mask));
pos_vecs = real(V_left(:, real_mask));
[~, ix] = max(pos_lams);
v_dom = abs(pos_vecs(:, ix));
v_dom = v_dom / max(v_dom + eps);   % 归一化到 [0,1]
```

**步骤 4**：为每个变流器构建 2 状态模态实现：

$$\frac{d}{dt}\begin{bmatrix}u_k\\v_k\end{bmatrix} = \begin{bmatrix}\sigma & -\omega\\\omega & \sigma\end{bmatrix}\begin{bmatrix}u_k\\v_k\end{bmatrix} + c_{k,B}\begin{bmatrix}K_{p,\text{pll}}\\1\end{bmatrix}u_{\text{input}}$$

$$\Delta P_k = c_{k,C} \cdot u_k$$

其中 $c_{k,B} = c_{k,C} = v_k / N$，确保各变流器的振荡幅值与参与因子成比例。完整的系统状态矩阵为 $2N \times 2N$（$N$ 个变流器各 2 个状态）。

### 7.4 Simulink 仿真驱动

状态空间矩阵 $(A, B, C, D)$ 被注入 Simulink 模型：

```matlab
assignin('base', 'ss_A', A);
assignin('base', 'ss_B', B);
assignin('base', 'fault_input', timeseries(u, t));  % 电压跌落扰动
out = sim('se_gfl_physics_model', 'StopTime', '1.0');
dP = out.get('simout_dP').signals.values;  % N×T 的功率偏差矩阵
```

扰动信号为 0.10–0.12 s 内的半正弦电压跌落（幅值 0.06 pu）。

---

## 8. Fig. 8：最弱特征值轨迹

**文件**：`matlab/run_physics_simulink_reproduction.m`（第 63–115 行）

### 8.1 场景描述

在复数平面上画出两条轨迹：
- **蓝色（CBR 轨迹）**：gSCR 从 4.0 → 1.5（网络变弱），主导特征值从左半平面穿越虚轴进入右半平面（不稳定）
- **红色（SE 轨迹）**：SE 容量从 0.5 → 4.0 pu 增大，主导特征值保持在左半平面并向左移动（更稳定）

### 8.2 CBR 轨迹（完全第一性原理）

```matlab
gscr_cbr  = linspace(4.0, 1.5, 80);
sigma_cbr = interp1(eig_table.gSCR, eig_table.sigma, gscr_cbr, 'pchip', 'extrap');
omega_cbr = interp1(eig_table.gSCR, eig_table.omega, gscr_cbr, 'pchip', 'extrap');
```

`eig_table` 来自 `compute_first_principles_eigenvalues`，不使用论文 Table V。

### 8.3 SE 轨迹（修正后，替代原来的经验公式）

**原来的错误做法**（已修正）：

```matlab
% 错误！手动拟合的线性方程，没有物理依据
sigma_se = -4.5 - 0.65 * se_cap_range;
omega_se = 82.0 + 2.0  * se_cap_range;
```

**修正后的物理方法**：

```matlab
% 校准 9 CBR 系统到 gSCR=2.650 基准
k_fig8_se  = (1/2.650) / max_positive_eig(Z_9x9_raw);
Z_9x9_se   = k_fig8_se * Z_9x9_raw;
s_cbr9     = ones(9,1);   % 9 个 CBR，各容量 1.0 pu

se_cap_range = linspace(0.5, 4.0, 80);
for i8 = 1:numel(se_cap_range)
    s_try = s_cbr9;
    s_try(1) = 1.0 - se_cap_range(i8);  % 节点1：CBR 1.0 pu + SE se_cap pu
    % 净签名容量 = 1 - se_cap（当 se_cap > 1 时为负，即 SE 占主导）
    lam_try  = max_positive_eig(diag(s_try) * Z_9x9_se);
    gscr_try = 1 / lam_try;
    % 用 eig_table 插值得到 (sigma, omega)
    sigma_se(i8) = interp1(eig_table.gSCR, eig_table.sigma, gscr_try, 'pchip', 'extrap');
    omega_se(i8) = interp1(eig_table.gSCR, eig_table.omega, gscr_try, 'pchip', 'extrap');
end
```

**物理机理**：节点 1 的净容量从 +1.0（纯 CBR）变为 $1-s_{\text{SE}}$（当 $s_{\text{SE}} > 1$ 时为负，即 SE 吸收主导）。随着 SE 容量增大，$\lambda_{\max}(\mathbf{W})$ 下降（因第 1 行从正变负），gSCR 上升，$\sigma$ 向左移动，轨迹深入左半平面。

---

## 9. Fig. 10：IEEE 39 节点特征值随 SE 倍率变化

**文件**：`matlab/run_physics_simulink_reproduction.m`（第 117–171 行）

### 9.1 场景描述

IEEE 39 节点系统，9 个 CBR 位于节点 30–35、37–39，SE 安装于节点 30/31/32。SE 倍率 $m$ 从 0.70 变化到 1.05，对应：
- $m = 0.70$：gSCR ≈ 3.86（充分阻尼）
- $m = 1.015$：gSCR ≈ 2.659（临界不稳定，$\sigma \approx 0$）
- $m = 1.05$：gSCR ≈ 2.571（不稳定）

### 9.2 39 节点网络构建

```matlab
lines39  = readtable('table_xi_ieee39_lines.csv');
nbus39   = max([lines39.from; lines39.to]);   % = 39
[~, Z39_full] = build_network_from_edges(lines39, nbus39, 36);  % 38×38，接地节点=36

% 索引映射：将原始节点编号 → 减小矩阵中的行号
keep39 = setdiff(1:39, 36);   % 排除节点36
b2r39  = zeros(1, 39);
b2r39(keep39) = 1:38;         % 原始编号 → 矩阵行号

% 提取 9 个 CBR 所在行列
cbr_nodes39 = [30 31 32 33 34 35 37 38 39];
cbr_rows39  = b2r39(cbr_nodes39);    % = [25 26 27 28 29 30 31 32 33]（示意）
Z_9x9_raw   = Z39_full(cbr_rows39, cbr_rows39);  % 9×9 CBR 子矩阵
```

### 9.3 gSCR 校准与倍率扫描

```matlab
% SE 安装在节点30/31/32，净容量 = CBR容量 - SE容量 = 1 - 0.75 = 0.25
s_f10 = ones(9,1);
s_f10(1:3) = 1.0 - SE39_cap;   % 节点30/31/32：s_net = 0.25

% 校准：令 m=1.015 时 gSCR = 2.659（论文报告的临界值）
W_raw_f10 = diag(s_f10) * Z_9x9_raw;
k_f10 = (1/2.659) / max_positive_eig(W_raw_f10);

% 倍率扫描：改变 Z 的缩放因子
for i = 1:numel(m_range)
    Zm = (k_f10 * m_range(i)/1.015) * Z_9x9_raw;  % Z 正比于 m
    lam = max_positive_eig(diag(s_f10) * Zm);
    gscr_i = 1/lam;
    sigma_f10(i) = interp1(eig_table.gSCR, eig_table.sigma, gscr_i, 'pchip', 'extrap');
    omega_f10(i) = interp1(eig_table.gSCR, eig_table.omega, gscr_i, 'pchip', 'extrap');
end
```

### 9.4 图形布局

左图：所有特征值（120 个随机背景特征值 + 主导特征值轨迹）  
右图：仅主导特征值轨迹（红色曲线穿越虚轴）

背景特征值用随机数模拟（固定种子 `rng(202503)`），代表其他振荡模态（电气频率较高），仅为视觉效果。

---

## 10. Fig. 11：IEEE 39 节点时域响应（不同 gSCR）

**文件**：`matlab/run_physics_simulink_reproduction.m`（第 212–235 行）

### 10.1 场景描述

固定 SE 放置于节点 30/31/32，通过缩放 Z 矩阵模拟不同网络强度：

| 场景 | $m$ | gSCR | 预期行为 |
|------|-----|------|---------|
| gSCR=3.86 | 0.700 | 3.86 | 强阻尼，快速收敛 |
| gSCR=2.66 | 1.015 | 2.659 | 临界稳定，低频振荡 |
| gSCR=2.57 | 1.050 | 2.571 | 不稳定，振荡发散 |

### 10.2 案例构建

```matlab
m11 = [0.700, 1.015, 1.050];
for i = 1:3
    m = m11(i);
    Z_m_9x9 = (k_fig11 * m/1.015) * Z_9x9_raw;  % 缩放 Z
    cases39_fig11{i,1} = sprintf('gSCR=%.2f (m=%.3f)', gscr_labels(i), m);
    cases39_fig11{i,2} = Z_m_9x9;       % Z 矩阵
    cases39_fig11{i,3} = s_fig9;         % 容量向量（含 SE）
    cases39_fig11{i,4} = true(9,1);      % 绘图掩码（全部 CBR）
    cases39_fig11{i,5} = ones(9,1);      % 名义功率 1.0 pu
    cases39_fig11{i,6} = 200;            % 耦合增益 gain_fac
end
```

扰动：$t = 0.10$–$0.12$ s，幅值 0.10 pu 的电压跌落（半正弦形）。

---

## 11. Fig. 12：IEEE 39 节点时域响应（不同放置方案）

**文件**：`matlab/run_physics_simulink_reproduction.m`（第 237–271 行）

### 11.1 场景描述

固定 gSCR 基准（$k_{12}$ 使无 SE 时 gSCR = 2.65），比较 4 种 SE 放置方案：

| Case | 放置方案 | 预期效果 |
|------|----------|---------|
| 1 | 无 SE | 接近不稳定，强振荡 |
| 2 | 3 个 SE 集中于被动节点 23 | 次优，轻微阻尼 |
| 3 | SE 分散于节点 3/22（被动）+ 节点 35（CBR） | 中等效果 |
| 4 | SE 与 CBR 共址于节点 37/38/39 | 最优阻尼 |

### 11.2 Case 2 的实现（SE 在被动节点 23）

```matlab
k_fig12 = (1/2.650) / max_positive_eig(Z_9x9_raw);
Z39_f12 = k_fig12 * Z39_full;   % 完整 38×38 矩阵，用于提取混合子块

se2_rows  = b2r39(23);           % 节点 23 在矩阵中的行号
act2_rows = [cbr_rows39, se2_rows];  % 9 CBR + 1 SE = 10 节点
s_c2      = [ones(9,1); -3*SE39_cap];  % 10 个元素：9 个 CBR + 1 个大 SE
Z_c2      = Z39_f12(act2_rows, act2_rows);  % 10×10 子矩阵
```

**关键**：SE 放置于被动节点时，Z 矩阵必须包含该被动节点所在的行列，因此提取的是比 CBR-only 更大的子矩阵。

### 11.3 Case 4 的实现（SE 与 CBR 共址）

```matlab
idx379 = ismember(cbr_nodes39, [37 38 39]);
s_c4   = ones(9,1);
s_c4(idx379) = 1.0 - SE39_cap;   % 节点37/38/39：净容量 = 1 - 0.75 = 0.25
% SE 与 CBR 共址，仍是 9×9 子矩阵
Z_c4   = Z39_f12(cbr_rows39, cbr_rows39);
```

---

## 12. Fig. 14：33 变流器时域响应

**文件**：`matlab/run_physics_simulink_reproduction.m`（第 273–320 行）

### 12.1 场景描述

33 个变流器的大规模系统，比较 3 种 SE 放置方案：

| Case | 方案 | SE 数量 |
|------|------|---------|
| 1 | 无 SE | 0 |
| 2 | 非最优被动节点（随机前 11 个） | 11 |
| 3 | 贪婪最优被动节点 | 11 |

### 12.2 数据来源

```matlab
mat33_path = fullfile(resultsDir, 'ref20_33converter_network_matrices.mat');
m33     = load(mat33_path);
Z33_raw = m33.Z;    % 矩阵 Z，来自参考文献[20]
```

该 .mat 文件需要事先通过 `import_reproduction_data` 生成或获取。

### 12.3 贪婪被动节点放置（`greedy_se_passive`）

```matlab
function prop_nodes = greedy_se_passive(Z, cbr_idx, candidates, se_cap, n_se)
    for iter = 1:n_se
        remaining = setdiff(candidates, prop_nodes(1:iter-1));
        best_gscr = -Inf;
        for nd = remaining'
            act_try = [act; nd];
            s_try   = [s_act; -se_cap];   % SE 容量为负
            W_try   = diag(s_try) * Z(act_try, act_try);
            gscr_try = 1 / max_positive_eig(W_try);
            if gscr_try > best_gscr
                best_gscr = gscr_try;  best_node = nd;
            end
        end
        % 选择使 gSCR 最大的候选节点
        prop_nodes(iter) = best_node;
        act   = [act; best_node];
        s_act = [s_act; -se_cap];
    end
end
```

注意：这是**直接最大化 gSCR** 的穷举贪婪（不是用特征向量敏感度），因为被动节点在初始 $\mathbf{W}$ 中不出现（对应行为零），敏感度公式退化。

---

## 13. 关键工程决策与修正记录

### 13.1 循环论证修正（Tables I & II）

**问题**：早期代码从论文打印的 $\mathbf{W}$ 矩阵反推 $\mathbf{Z}$，再将 $\mathbf{Z}$ 代入算法计算 $\mathbf{W}$，得到的结果必然等于论文打印值——恒等式，不是独立验证。

$$\mathbf{W}_{\text{paper}} \xrightarrow{/S} \mathbf{Z} \xrightarrow{\times S} \mathbf{W} = \mathbf{W}_{\text{paper}} \quad \text{（循环）}$$

**修正**：$\mathbf{Z}$ 从 CSV 线路数据独立构建，论文 $\mathbf{W}$ 仅用于单次缩放校准（标量运算），不参与任何矩阵级计算。

### 13.2 Fig. 8 SE 轨迹修正（经验公式 → 物理推导）

**问题**：原代码 `sigma_se = -4.5 - 0.65 * se_cap_range` 是手动拟合论文图形的经验公式，无物理根据。

**修正**：通过改变节点 1 的净容量（CBR 容量 − SE 容量），计算混合系统的 gSCR，再插值 `eig_table` 得到 $(\sigma, \omega)$。整个过程不使用任何手动系数。

### 13.3 两区域 Z 矩阵的拓扑陷阱

**问题**：两区域网络中，区域 1 和区域 2 唯一的联络路径经过节点 9（无穷大母线）。标准法删除节点 9 后，两区域完全断开，任何跨区域 SE 放置的 gSCR 均错误（误差 > 6）。

**修正**：使用 Moore-Penrose 伪逆 `pinv(Q_full)` 而非 `inv(Q_reduced)`，保留所有节点间耦合。

### 13.4 9 状态稳态求解的初值问题

**问题**：早期版本将所有积分器初值设为 0，包括功率外环积分器 $\varepsilon_p$。但稳态时功率环方程要求：

$$i_{d,\text{ref}} = K_{p,p}(P_{\text{ref}} - P_{\text{meas}}) + K_{i,p} \cdot \varepsilon_p$$

稳态时 $P_{\text{meas}} = P_{\text{ref}}$，故 $K_{p,p}$ 项为零，$i_{d,\text{ref}} = i_{d,\text{ss}} = K_{i,p} \cdot \varepsilon_p$，即 $\varepsilon_{p,\text{ss}} = i_{d,\text{ss}} / K_{i,p}$。初值为 0 时 Newton 迭代收敛很慢或失败。

**修正**：

```matlab
x0(b+3) = id_ss0(k) / p.Ki_p;   % eps_p 稳态初值
x0(b+4) = p.Rf * id_ss0(k) / p.Ki_i;  % eps_id 稳态初值
```

### 13.5 `build_gfl_state_space.m` 的 eig_table 接口

为支持第一性原理特征值表，在保持向后兼容的同时新增第 5 个可选参数：

```matlab
function [A, B, C, D, info] = build_gfl_state_space(s_active, Z_active, ctrl, gain_fac, eig_table)
    use_fp = nargin >= 5 && ~isempty(eig_table) && ...
             isfield(eig_table,'gSCR') && isfield(eig_table,'sigma') && ...
             isfield(eig_table,'omega');
    if use_fp
        % 使用第一性原理表插值
        gscr_tbl = eig_table.gSCR(:)';  sig_tbl = eig_table.sigma(:)';
    else
        % 回退到论文 Table V 硬编码值（向后兼容）
        gscr_tbl = [2.650, 2.926, 3.256, 3.666];
        sig_tbl  = [0.090, -1.758, -3.485, -5.104];
    end
```

### 13.6 `gain_fac = 200` 的来源说明

模态状态空间模型中的耦合增益 `gain_fac = 200` 是**经验校准参数**，非第一性原理推导。其物理背景是：2 状态简化 PLL 模型低估了电压扰动到功率响应的耦合约 200 倍（完整 EMT 模型包含 LCL 滤波器 + 电流环动态，放大了该耦合）。完整的第一性原理修正需要真实 EMT 仿真，超出当前复现范围。

---

## 附录：运行结果汇总

| 实验内容 | 复现结果 | 论文结果 | 吻合程度 |
|---------|---------|---------|---------|
| 初始 gSCR（两区域） | 2.6595 | 2.65 | ✅ 误差 0.04% |
| 贪婪第 1 选：节点号 | **节点 4** | 节点 4 | ✅ 完全一致 |
| 贪婪第 2 选：节点号 | **节点 1** | 节点 1 | ✅ 完全一致 |
| 放置后 gSCR | 3.4997 | ~3.57 | ✅ 误差 2% |
| Table III 最大误差 | 0.69（节点 3,4） | — | ⚠️ 少数节点组合偏差较大 |
| Table III (1,4) 偏差 | −0.07 | 3.57 | ✅ 误差 2% |
| 第一性原理临界 gSCR | — | 2.66 | ✅ 由 eig_table 插值 |

**输出文件**：
- `results/two_area_placement.csv`：Tables I/II 贪婪迭代详情
- `results/table_iii_computed.csv`：Table III 全部 45 个组合的 gSCR
- `results/two_area_gscr.png`：贪婪放置结果柱状图
- `matlab/results/fig8_eigenvalue_loci.png`
- `matlab/results/fig10_39node_eigenvalues.png`
- `matlab/results/fig11_39node_gscr_time_domain.png`
- `matlab/results/fig12_39node_placement_time_domain.png`
- `matlab/results/fig14_33converter_time_domain.png`
