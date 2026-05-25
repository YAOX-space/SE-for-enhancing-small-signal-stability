# 论文复现：Storage Energies Placement for Small-Signal Stability

复现对象：

Yuan et al., **"Placing Storage Energies for Enhancing Small-Signal Stability of Converter-Based-Renewable Systems"**, IEEE Transactions on Industry Applications, 2025.

## 已复现内容

本仓库复现了论文第 IV-V 节的核心算法，并补充了 MATLAB 降阶时域响应脚本：

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

此外，`matlab/` 目录下提供完整批处理脚本，生成对应 Fig. 8、Fig. 10、Fig. 6、Fig. 11、Fig. 12、Fig. 14 的复现实验图。时域部分使用透明的 MATLAB 降阶响应模型，包含故障期间平滑有功下挫和故障清除后的 PLL 主导振荡，复现论文稳定性趋势和相对阻尼效果；由于论文没有提供原始 EMT `.slx`，该部分不是作者 EMT 模型的逐块重建。

## 运行

Python 两区四 CBR 示例：

```bash
PYTHONPATH=src python3 examples/two_area_demo.py
```

MATLAB 全部复现实验：

```bash
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('/Users/yao/Documents/zju/matlab'); run_all_reproduction"
```

单独构建并运行 Simulink 近似模型：

```bash
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('/Users/yao/Documents/zju/matlab'); run_simulink_reproduction"
```

输出文件：

- `results/two_area_placement.csv`
- `results/two_area_gscr.png`
- `matlab/results/*.png`
- `matlab/results/*.csv`
- `matlab/results/ref20_33converter_network_matrices.mat`
- `matlab/se_gfl_pll_small_signal_model.slx`

## 文件结构

```text
src/se_placement.py        # gSCR、特征值灵敏度、贪心选址
examples/two_area_demo.py  # 两区四 CBR 示例
matlab/                    # MATLAB 全部实验复现
matlab/data/               # 从附录图片录入的参数表
matlab/results/            # MATLAB 输出图、表和网络矩阵
results/                  # 运行后生成的结果
```

## 说明

论文 PDF 只在正文式 (17) 给出了两区示例的初始 `S'_B1 Z`。当前仓库已根据用户补充的附录图片录入 Table IX-XIII，包括两区线路、IEEE 39 节点线路、控制参数，以及第二/第三轮灵敏度表。

Python 示例脚本仍利用 `Z` 的对称性从打印矩阵反推 CBR 容量比例：

```text
S_B = [0.5, 1.0, 1.5, 0.5]
```

并据此重构与式 (17) 一致的阻抗矩阵部分。该实现可以复现论文给出的两步选址 `[1, 4]`。
