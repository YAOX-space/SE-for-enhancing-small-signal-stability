# 论文复现：Storage Energies Placement for Small-Signal Stability

复现对象：

Yuan et al., **"Placing Storage Energies for Enhancing Small-Signal Stability of Converter-Based-Renewable Systems"**, IEEE Transactions on Industry Applications, 2025.

## 已复现内容

本仓库复现了论文第 IV-V 节的核心算法，并补充了 MATLAB 定性时域响应脚本：

1. 用加权网络阻抗矩阵计算广义短路比

   ```text
   gSCR = 1 / lambda_max_positive(S'_B1 Z)
   ```

2. 用最大正特征值的灵敏度进行储能选址：

   ```text
   d lambda_max / d S_BE,j = -a * lambda_max * v_j^2
   ```

3. 在论文两区四 CBR、9 节点示例上复现贪心选址结果：

   ```text
   第 1 个 SE: 节点 1
   第 2 个 SE: 节点 4
   ```

这对应论文 Table I-III 和 Fig. 4-6 附近的算法验证。

`matlab/` 目录下提供完整批处理脚本，生成对应 Fig. 8、Fig. 10、Fig. 6、Fig. 11、Fig. 12、Fig. 14 的定性插图。

## 复现范围与限制

### 算法与数值（精确复现）

| 内容 | 来源 | 验证结果 |
|---|---|---|
| gSCR 公式 (式19) | `src/se_placement.py` | λ_max = 0.376011，与论文 0.376 一致 |
| 灵敏度公式 (式16) | `src/se_placement.py` | 节点排名与 Table I/II 完全一致 |
| 两区贪心选址 | `examples/two_area_demo.py` | [1, 4]，与 Table III 全局最优一致 |
| IEEE 39 节点选址 | `matlab/run_physics_simulink_reproduction.m` | Table XII→节点38，Table XIII→节点37 |
| 33 变流器选址 | 同上（`greedy_se_passive`） | 独立计算所得11节点与论文一致 |

**33 变流器第一轮注意**：Table VIII 的灵敏度值截尾至 4 位小数，导致节点 45/56/67 三者并列（均为 -0.0122）。论文标注节点 67，与本仓库贪心算法从 Z 矩阵直接计算的结果一致，但截尾精度不足以从 Table VIII 数值单独验证此选择；原始未截尾数据方可区分。

### 时域仿真（定性趋势插图，非 EMT 复现）

Fig. 6、Fig. 11、Fig. 12、Fig. 14 的时域曲线**不是**论文 EMT 模型的逐块重建，而是由 `matlab/lib/build_gfl_state_space.m` 构造的模态降阶状态空间模型生成的**定性插图**。该模型以论文 Table V 报告的 (σ, ω) 四组数据点对 gSCR 做插值，再构建 2N×2N 状态空间，因此：

- 能定性复现"gSCR 越低阻尼越差、SE 放置越优收敛越快"的稳定性趋势；
- 不能用于独立验证论文的稳定性数值结论（模型本身已用论文数据标定）；
- 与作者 EMT 波形（含电流内环、LCL 滤波器、限幅器、dq 初始化）不等价。

从第一原理推导 GFL 变换器状态矩阵时，静态网络耦合项 O(0.02) 远不足以将 PLL 模式阻尼从 −13 移动到论文报告的 +0.09。差距来源于 Table IX 中未完整公开的电压前馈滤波器 1/(1+0.01s) 动态，该部分引入的相位滞后在 88 rad/s 处产生强电流-网络耦合，是复现真实极点轨迹的必要条件。

### 两区四 CBR 示例的范围限制

`examples/two_area_demo.py` 从论文式 (17) 打印的加权矩阵 W 反推 Z，只能恢复 Z 的前 4 行（CBR 行）及对应列；被动节点的 5×5 子块（行/列 4-8）保持为零。脚本中 `candidate_nodes=list(range(9))` 虽包含全部 9 个节点，但由于被动节点行全零，其灵敏度恒为零，算法不会选中它们——这与本示例的最优结果（CBR 节点 1、4）一致，因此结论正确。**若需超过两步的全节点贪心或将 SE 放在被动节点，应改用由 `table_x_two_area_lines.csv` 完整构建的 Z 矩阵。**

## 运行

**Python 两区四 CBR 示例（Windows）：**

```powershell
$env:PYTHONPATH="src"
py -3 examples/two_area_demo.py
```

**Python 两区四 CBR 示例（Linux/macOS）：**

```bash
PYTHONPATH=src python3 examples/two_area_demo.py
```

**MATLAB 全部复现实验：**

在 MATLAB 中进入仓库的 `matlab/` 子目录，然后运行：

```matlab
run_physics_simulink_reproduction
```

脚本通过 `fileparts(mfilename('fullpath'))` 自动定位自身路径，无需修改任何硬编码路径。

输出文件：

- `results/two_area_placement.csv`
- `results/two_area_gscr.png`
- `matlab/results/*.png`（Fig. 6/8/10/11/12/14）
- `matlab/results/ref20_33converter_network_matrices.mat`

## 文件结构

```text
src/se_placement.py                     # gSCR、特征值灵敏度、贪心选址（核心算法）
examples/two_area_demo.py               # 两区四 CBR 示例（仅复现前两步选址）
matlab/run_physics_simulink_reproduction.m  # 全部实验入口（单一脚本）
matlab/lib/build_gfl_state_space.m      # Table V 标定 → 模态状态空间（定性时域用）
matlab/lib/build_network_from_edges.m   # CSV 边表 → Laplacian → Z 矩阵
matlab/lib/max_positive_eig.m           # lambda_max+
matlab/data/                            # 从附录图片录入的参数表（CSV）
matlab/results/                         # 输出图和网络矩阵
results/                                # Python 脚本输出
```

## 数据来源说明

论文 PDF 只在正文式 (17) 给出两区示例的初始 `S'_B1 Z`（4 行×9 列）。仓库已根据附录图片录入 Table IX-XIII，包括两区线路、IEEE 39 节点线路、控制参数，以及第二/第三轮灵敏度表。CBR 容量比例 `S_B = [0.5, 1.0, 1.5, 0.5]` 由 Z 对称性反推得到。
