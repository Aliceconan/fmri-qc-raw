# fmri-qc-raw 指标定义与验证记录

每个指标的**精确定义、单位、怎么算出来的、验证过没有**。改代码或和别的
pipeline 对拍时看这里。

---

## 1. 头动参数

### 1.1 AFNI 原始输出

`3dvolreg -1Dfile` 输出 6 列，**旋转在前**：

```
roll  pitch  yaw  dS  dL  dP
```

| 列 | 含义 | 单位 |
|---|---|---|
| `roll` | 绕 I-S 轴（= z）旋转，CCW 为正 | degree |
| `pitch` | 绕 R-L 轴（= x）旋转 | degree |
| `yaw` | 绕 A-P 轴（= y）旋转 | degree |
| `dS` | 向 **S**uperior 方向位移 | mm |
| `dL` | 向 **L**eft 方向位移 | mm |
| `dP` | 向 **P**osterior 方向位移 | mm |

### 1.2 换算到 BIDS/RAS

有**两层**符号问题，必须同时处理，只处理一层会得到全反的结果：

1. **轴向**：AFNI 用 S/L/P，BIDS 用 RAS。`Left = -x`、`Posterior = -y`、`Superior = +z`
2. **方向**：3dvolreg 给的是「把该卷配回基准所需的校正量」，即**头动的逆**

两层叠加：

```
trans_x = +dL        trans_y = +dP        trans_z = -dS          (mm)
rot_x   = -pitch     rot_y   = -yaw       rot_z   = -roll        (deg → rad)
```

输出语义是**头相对基准卷动了多少**：`trans_x = +3` 表示头向右移动了 3mm。

### 1.3 验证记录

用合成数据验证（40×40×24，体素 3×3×3.5mm，affine 全正对角，
不对称脑 + 3 个偏心团块，高斯噪声 σ=6）：

| 真值 | 恢复值 |
|---|---|
| trans (3, 0, 0) mm | (3.00, −0.01, −0.01) |
| trans (0, 6, 0) mm | (0.00, 6.00, −0.01) |
| trans (0, 0, −3.5) mm | (0.00, −0.00, −3.50) |
| trans (−6, 3, 3.5) mm | (−5.99, 3.00, 3.50) |
| rot_z +3° | (0.01, −0.02, **+3.04**) |
| rot_z −4° | (0.01, −0.01, **−3.96**) |
| rot_z +6° | (0.00, −0.02, **+6.04**) |

**平移最大误差 0.009 mm，rot_z 最大误差 0.042°，轴间串扰 ≈ 0。**

> ⚠️ 合成测试数据必须**面内不对称**。用旋转对称的球/椭球会让刚体配准锁不住
> in-plane rotation，产生几度的假 `rot_z`，经 50mm 弧长折算后 FD 虚高 5–7 倍。
> 这个坑踩过一次。

---

## 2. 位移类指标

### 2.1 framewise_displacement（Power FD）

```
FD_t = Σ|Δtrans_i| + Σ|Δrot_i(rad)| × r,    r = 50 mm
```

Power et al. 2012 的定义，旋转按 50mm 等效脑半径折算成弧长。**单位 mm**。
第一帧为 `n/a`。基于一阶差分，与基准卷选择无关。

常用门限：任务态 0.5mm，静息态功能连接 0.2–0.3mm。

#### 与 fmriprep 完全一致（已实测验证）

fmriprep 的 FD 来自 `nipype.algorithms.confounds.FramewiseDisplacement`，
核心就三行（nipype 1.11.0，默认 `radius=50`、`normalize=False`）：

```python
diff = mpars[:-1, :6] - mpars[1:, :6]
diff[:, 3:6] *= self.inputs.radius        # 50
fd_res = np.abs(diff).sum(axis=1)
```

拿真实数据（任务A 的 aborted run-03，104 帧）把 3dvolreg 的原始 `.1D` 直接喂给
`FramewiseDisplacement(parameter_source='AFNI')`，和本工具的
`framewise_displacement` 列对比：

```
mean_FD : nipype 0.117835224   qc-raw 0.117835224
max_FD  : nipype 0.583230077   qc-raw 0.583230077
逐帧最大绝对差 = 2.08e-16      ← 浮点精度极限
```

**为什么符号约定不同却完全相等。** nipype 的 `normalize_mc_params(source='AFNI')`
只做重排 + 度转弧度，得到 `[dL, dP, dS, pitch, yaw, roll]`，**不做**本工具那两层
符号翻转（见 §1.2）。所以两边的 6 参数在 z 和三个旋转上差一个负号：

```
            trans_x   trans_y   trans_z    rot_x     rot_y     rot_z
nipype     -0.297300 -1.770400 -1.125900 -0.014860 -0.004372 -0.001883
qc-raw     -0.297300 -1.770400 +1.125900 +0.014860 +0.004372 +0.001883
```

但 FD 是 `|一阶差分|` 求和——**逐列翻符号、换轴序都不改变结果**。
nipype 不需要处理轴向和逆变换，因为它的输出只用来算 FD，不用来写 confounds
（fmriprep 的 `trans_*`/`rot_*` 列走的是 mcflirt→FSL 路径，不是 AFNI 路径）。

结论：**`framewise_displacement` 可以和 fmriprep 的同名列直接比**，
差异只会来自运动估计本身（3dvolreg vs mcflirt）和参考卷不同，不来自公式。

### 2.2 enorm_afni

```
enorm_t = ‖Δ(roll, pitch, yaw, dS, dL, dP)_t‖₂
```

afni_proc.py 的 `motion_enorm`。注意 AFNI **直接把 degree 当 mm 用**，不乘半径。
`afni_proc.py -regress_censor_motion 0.2/0.3` 的阈值就是按这个口径定的，
**不要拿 FD 的阈值套 enorm，也不要给 enorm 做单位换算**。

### 2.3 maxdisp / maxdisp_range / maxdisp_delt

来自 `3dvolreg -maxdisp1D`，是脑内体素（`3dAutomask -clfrac 0.33` 定义）
相对基准卷的**实际最大位移**，不依赖 50mm 这种经验半径，物理意义比 FD 更直接。

| 列 | 含义 | 依赖基准卷？ |
|---|---|---|
| `maxdisp` | 每个 TR 相对基准卷的最大体素位移 | **是** |
| `maxdisp_max` | 上者在 run 内的最大值 | **是**（含 run 间偏移） |
| `maxdisp_range` | `ptp(maxdisp)`，run 内极差 | 否 ← **判 run 内头动用这个** |
| `maxdisp_delt` | 相邻 TR 之间最大位移的变化（AFNI 写到 `*_delt` 文件） | 否 |

> ⚠️ `--motion-base global` 时基准卷可能在**别的 run** 里，`maxdisp_max`
> 因此含 run 间偏移。曾经拿它当「run 内最大位移」判定，导致真实数据上
> `maxdisp_max 12.54mm ≈ betrun_disp 12.72mm` —— 报的其实是同一件事。
> 已改用 `maxdisp_range`。

### 2.4 run 间头动

```
betrun_disp(A, B) = Σ|Δmedian_trans| + Σ|Δmedian_rot(rad)| × 50mm
```

取每个 run 的**中位姿态**（对 6 个参数逐列取中位数），两 run 中位姿态之差按
Power 公式折算成一个 mm 数。

| 列 | 参照点 | 用途 |
|---|---|---|
| `betrun_disp_vs_median` | 组内各 run 中位姿态的逐分量中位数 | **判定用这个** |
| `betrun_disp_vs_ref` | 配准基准卷所在的 run | 诊断/溯源 |
| `betrun_disp_vs_prev` | 采集顺序上的前一个 run | 定位「哪次挪了」 |

⚠️ **不能拿 `vs_ref` 判定**：基准卷按 min-outlier 选（信号最干净），
可能正好在头位分布边缘，会把结论判反。实测案例见
[`pitfalls.md`](pitfalls.md) §5a。

**只在 `--motion-base global` 下有效**，因为需要所有 run 处在同一坐标系。
`per-run` 模式下这两列是 `n/a`。

---

## 2.5 阈值出处（重要：分清哪些有依据、哪些是本工具自定的）

改阈值前先看这张表。**标 ❌ 的没有文献依据，请按自己项目重新标定。**

| 阈值 | 默认值 | 出处 |
|---|---|---|
| `fd_thresh` | 0.5 mm | ✅ Power et al. 2012 的 scrubbing 阈值。静息态功能连接常用更严的 0.2 |
| `outlier_thresh` | 0.05 | ✅ `afni_proc.py -regress_censor_outliers` 的默认值 |
| enorm 判读参考 | 0.2 / 0.3 | ✅ `afni_proc.py -regress_censor_motion` 的常用值 |
| `fd_mean_warn/fail` | 0.2 / 0.5 mm | ⚠️ 静息态文献常见的被试排除线（0.2 / 0.25 / 0.5 都有人用），无单一权威出处 |
| `pct_bad_warn/fail` | 10% / 25% | ⚠️ 「censored 帧超过 20–25% 就弃用」是常见做法，无权威出处 |
| `maxdisp_warn/fail_vox` | 1 / 2 个体素 | ❌ **本工具自定。**「头动不超过一个体素」是圈内经验法则，找不到可引的原始文献 |
| `betrun_warn` | 3 mm | ❌ **本工具自定**，无文献标准 |
| `betrun_fail` | 6 mm | ❌ 无文献标准，但有**实测锚点**：7T 1.2mm 数据上某 run 相对基准 9.05mm，该被试配准效果明显变差 |
| `brain_edge_warn/fail` | 0.005 / 0.02 | ❌ 自定，但**按实测标定**：一批覆盖完好的 7T 数据贴边占比 ≤0.0003 |
| `coverage_loss_warn/fail` | 5% / 15% | ❌ 自定，实测基线：同一 session 内组内覆盖差异 ≤2.3% |
| `tsnr_warn/fail` | 40 / 20 | ❌ **文献中不存在此类门限**（查证见 §4.1）。只作兜底，必须按自己协议标定。7T 1.2mm iso 实测标定值：20 / 12 |
| `tsnr_rel_warn/fail` | 70% / 50% | ⚠️ 相对队列中位数。做法本身有 MRIQC 规程背书，具体百分比是本工具自定 |
| `drift_warn_pct` | 5% | ❌ 自定 |
| `truncated_fail/warn_ratio` | 60% / 95% | ❌ 自定 |

**怎么给自己的项目定阈值**：拿一批你认可「质量没问题」的数据跑一遍，看各指标的
实际分布，把门限设在明显高于该基线的位置。上面 `brain_edge` 和 `coverage_loss`
就是这么定的。

---

## 2.6 run 间位移为什么判 fail

理论上 run 间位移会被预处理的配准消掉，看起来不该判 fail。**但实测不是这样**：
7T 1.2mm 数据上一个相对基准 9.05mm 的 run，配准效果明显变差。

原因是位移大到一定程度后，harm 不再只是"位置不同"：

1. **磁化率畸变随头位变化** —— EPI 的几何畸变取决于头在磁场里的位置，
   头位差太多，畸变模式就不一样，刚体配准对不齐
2. **fieldmap 失效** —— reverse-PE 采集时的头位只对邻近的 run 有效
3. **配准优化器落到局部极小** —— 初始偏移太大时刚体配准容易不收敛

所以本工具**保留 `betrun_fail`**，同时另外直接测一个具体后果：

| 指标 | 含义 |
|---|---|
| `brain_edge_frac` | 脑体素贴在 FOV 六个面上的最大占比 → 脑被切掉 |
| `coverage_rel_pct` | 该 run 的脑体素数 / 组内最多的 run → 看到的脑变少了 |

这两个来自 Phase A 已经算好的 automask，零额外成本。
位移大但覆盖没损失，也仍然 fail —— 因为畸变和配准的问题不体现在覆盖上。

---

## 2.7 func 与配对 fmap 的头位差（blip 矫正的前提）

blip / topup 假设 func 和 reverse-PE fmap 处在**同一头位**。头动过大则估出来的
场图不适用于该 run，矫正会失败。

### 配对规则

优先级：BIDS `IntendedFor` → `run` 实体号相同 → **采集时间最近**。
时间最近是最稳的兜底：fmap 少于 run 时（一个 fmap 服务多个 run），
只有时间距离才有物理意义。输出列 `fmap_pair_method` 记录用了哪条。

### ⚠ 刚体配准会把畸变差异算成头动

func 与 fmap 相位编码方向相反 → 几何畸变也相反。直接拿 6 dof 位移判定会
**系统性高估**。6 个自由度里只有 3 个不受污染：

| 分量 | 干净？ | 理由 |
|---|---|---|
| 垂直 PE 的两个平移 | ✅ | 沿 PE 的位移场造不出另外两轴的位移 |
| 绕 PE 轴的旋转 | ✅ | 该旋转只动另外两轴 |
| PE 轴平移 | ❌ | 直接叠加畸变 |
| 绕另外两轴的旋转 | ❌ | PE 位移随 z 变化 ≈ 绕 x 转；随 x 变化 ≈ 绕 z 转 |

所以 `fmap_disp_clean` 只累加干净的三个分量，**判定用它**；
`fmap_disp_pe` 单独报出（标注为头动+畸变混合），`fmap_disp_full` 仅供参考。

**局限（必须知道）**：func↔fmap 之间沿 PE 的**真实**平移、绕另两轴的**真实**
旋转抓不到 —— 单个时刻只有一种极性，物理上无法分离。这是个偏保守的检查：
报了一定有问题，没报不代表一定没问题。

### 阈值

`fmap_disp_warn_vox` / `fmap_disp_fail_vox`，单位是面内体素倍数，默认 1 / 2。
❌ 无文献依据，来自实践经验（超过一个体素的错位，blip 矫正的收益就被抵消）。

---

## 2.8 B0 漂移估计（7T 长时程扫描）

### 现象

同一段时间跨度，两个**同极性**比较给出符号相反的位移（实测，7T 1.2mm，2 小时）：

```
fmap→fmap（都是 PA）  ty = +2.98 mm
func→func（都是 AP）  ty = −5.54 mm
```

同一个头、同一段时间，真实头动不可能既 +3 又 −5.5。必然有一个随 PE 极性反号的
量在随时间变化 —— 那只能是几何畸变。

### 分解

```
真实头动(PE 轴) = (Δ_AP + Δ_PA)/2
畸变漂移        = (Δ_AP − Δ_PA)/2
```

实测：真实头动 ≈ **1.28 mm**，畸变漂移 ≈ **4.26 mm**。
也就是说该被试 PE 轴上被报出的 5–7mm「run 间位移」里，大部分是 B0 漂移不是头动。

### 独立佐证

`TotalReadoutTime = 56.6 ms`，畸变位移(体素) = ΔB0[Hz] × TRT：

```
4.3 mm = 3.55 体素  →  ΔB0 ≈ 63 Hz
```

7T 两小时扫描、梯度线圈发热导致的 B0 漂移，文献量级正是数十 Hz。
加上漂移是**单调**的（符合发热曲线，不像头动那样随机），解释自洽。

### ⚠ 这是估计，不是测量

- **刚体配准测不出畸变场**。它只给一个全局平移，是真实位移场在全脑上的
  某种加权平均。
- 上面的分解依赖「畸变严格反对称」的假设，且刚体拟合对畸变图像的响应
  并非严格线性。
- 结论「畸变随时间变了」是稳的（只需同极性比较符号相反）；
  **具体拆出的数值是量级估计**。

输出列：被试级的 `b0_drift_pe_mm`、`head_drift_pe_mm`。
超过 `b0_drift_warn`（默认 2mm）会提示：远端 run 用同一个 fmap 做 blip 可能对不齐。

---

## 3. 基准卷的选择

按 afni_proc.py 的 `align_opts` 惯例用 **min-outlier volume**：

1. `3dToutcount -automask -fraction -legendre` 得到每个 TR 的 outlier 比例
2. 跳过前 2 帧（躲开非稳态帧），取 outlier 最低的那一帧

`global` 模式下再进一步：同一几何组内，取**平均 outlier 最低的那个 run** 的
min-outlier 卷作为全组基准。

**几何分组键**：`(subject, session, dim_str, vox_str, pe_dir)`。
分辨率/矩阵/PE 不同的 run 之间做刚体配准没有意义，必须分开。

---

## 4. 信号质量指标

| 指标 | 怎么算 | 说明 |
|---|---|---|
| `outlier_fraction` | `3dToutcount -automask -fraction -legendre` | 每个 TR 中被判为时间序列离群的体素比例 |
| `dvars` | `3dTto1D -method dvars` | **在未做 MoCo 的原始数据上算**，因此比预处理后偏高，只做相对比较 |
| `tsnr_median` | `3dvolreg` 输出 → `3dTstat -cvarinv` → 脑掩膜内取中位数 | `-cvarinv` 会先去掉均值+线性趋势再算 stdev。**必须在 MoCo 之后算**，否则测的是头动 |
| `gs_drift_pct` | 全脑均值时间序列（跳过前 5 帧）线性拟合，`\|slope × N\| / mean × 100` | 扫描仪漂移 / 被试逐渐移出视野 |
| `n_nonsteady_est` | 全脑均值前 10 帧中，偏离稳态中位数 > 5×MAD 的**连续前缀**长度 | 未剔除的 dummy 或未达稳态 |

### 4.1 tSNR 没有可引的绝对门限（查过文献）

**结论：文献里不存在「tSNR 低于 X 就不合格」这种标准。** 查证过程：

- **[Triantafyllou et al. 2005, NeuroImage](https://pubmed.ncbi.nlm.nih.gov/15862224/)**
  —— 这是 tSNR 与场强/体素的经典参考。它给的是**物理关系**不是门限：
  tSNR 随体素体积线性上升，直到撞上生理噪声天花板；场强越高，天花板在越大的
  体素处就达到。所以 **7T 高分辨率数据处在热噪声主导区**，tSNR 天然低，
  而且大致正比于体素体积。
- **[Hutton et al. 2011, NeuroImage](https://pmc.ncbi.nlm.nih.gov/articles/PMC3115139/)**
  —— 7T 视觉皮层报到 tSNR 88.9（未做生理噪声校正）/ 147.6（校正后），
  但该论文用了 1.1×1.1×1.8、2×2×2、3×3×2 三套参数，公开摘要里没说清这个
  数对应哪一套，**不能拿来当 1.2mm 的参考**。
- **[Gorgolewski et al., 7T 测重测数据集](https://pmc.ncbi.nlm.nih.gov/articles/PMC4412153/)**
  —— 1.5mm iso，tSNR 只画在图里，正文/表里没有数值。
- **[Provins et al. 2023, MRIQC/fMRIPrep QC 规程](https://www.frontiersin.org/journals/neuroimaging/articles/10.3389/fnimg.2022.1073734/full)**
  —— 明确说排除标准**取决于具体项目**，做法是「筛出平均 tSNR 过低的那些被试」，
  即**相对本队列**判断，不给绝对数。

高分辨率 7T 的 tSNR 值在文献里多以 tSNR map 的形式出现在图里，很少表格化，
而且强烈依赖线圈、ROI、是否做生理噪声校正、是否 denoise。**没有可移植的绝对数。**

### 因此本工具用两道 tSNR 判定

| 判定 | 说明 |
|---|---|
| `tsnr_warn` / `tsnr_fail` | **绝对下限**，只当「明显坏掉」的兜底。⚠ 必须按自己的场强/分辨率/线圈重标定，默认值面向 3T 3mm |
| `tsnr_rel_warn` / `tsnr_rel_fail` | **占同协议队列中位数的百分比**（默认 70% / 50%）。这是 MRIQC 推荐的做法，换场强换序列都不用改，是更可靠的那一道 |

按 `(task, 矩阵, 体素)` 分组，组内 run 数 ≥ `tsnr_rel_min_n`（默认 3）才算相对值。

### 怎么给自己的协议定绝对下限

跑一批你认可没问题的数据，看实测分布，把 warn 设在中位数的 ~2/3、fail 设在 ~1/3。

实测记录（`3dTstat -cvarinv`，脑掩膜内中位数）：

| 协议 | 实测 tSNR | 建议 warn / fail |
|---|---|---|
| 7T 1.2mm iso GE-EPI, TR=2, 124 层, 32ch | 27.4 – 33.7（中位 29.9，8 个 run）| 20 / 12 |
| 3T ~3.4×3.4×4mm GE-EPI, TR=2 | 27.4 – 52.3（12 个 run，另一数据集）| — |

⚠ 上面第一行是**本项目实测标定**，不是文献值。换线圈、换 TE、换 flip angle
都要重测。

### 非稳态帧的判读

真实数据上第 1 帧常常是这样（3T，455 TR）：

```
global_signal   前10: [428.0 519.3 525.6 526.7 528.3 ...]   稳态中位数 531.2
outlier_fraction 前10: [0.685 0.046 0.021 0.017 0.018 ...]
```

第 1 帧全脑信号比稳态低 19%、68.5% 的体素是 outlier —— 这一帧确定不能用。
（信号偏**低**而不是偏高也很常见，取决于序列的预脉冲设置。）

---

## 5. 与 fmriprep confounds 对拍

`motion/*_desc-rawmotion_timeseries.tsv` 刻意用 BIDS confounds 列名，
可以直接和日后的 `desc-confounds_timeseries.tsv` 对拍：

| 对得上 | 说明 |
|---|---|
| `framewise_displacement` | 都是 Power FD、mm、r=50mm。相关应 > 0.95 |
| `trans_*` / `rot_*` 的**量级和差分** | 单位一致（mm / rad） |

| 可能对不上 | 原因 |
|---|---|
| `trans_*` / `rot_*` 的**绝对值** | 基准卷不同（本工具用 min-outlier，fmriprep 用自己的参考） |
| 整体**符号** | 各家对「报头动还是报校正量」的约定不同。本工具报头动（见 1.2） |
| `dvars` | 本工具在原始数据上算，fmriprep 在预处理后算，绝对值差很多 |

行数一定相同（= BOLD volume 数，含 dummy）。缺失值一律 `n/a`，不写空、不写
`NaN`、不写 `0`。
