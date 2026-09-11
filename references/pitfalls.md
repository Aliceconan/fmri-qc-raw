# fmri-qc-raw 已知陷阱

每条都是实际踩过或必然会踩的。按「会不会让你得出错误结论」排序。

---

## A. 会让你得出错误结论的

### 1. `maxdisp_max` 不是 run 内头动

`--motion-base global` 下所有 run 配到同一基准卷，基准卷可能在别的 run 里，
所以 `maxdisp_max` 含 run 间偏移。

**症状**：真实数据上 `maxdisp_max 12.54mm` 和 `betrun_disp_vs_ref 12.72mm`
几乎相等——它们报的是同一件事。

**用 `maxdisp_range`**（run 内极差，与基准卷无关）判 run 内头动。已改。

### 2. 3dvolreg 的 6 个参数是「校正量」，不是「头动」

参数描述的是把该卷配回基准所需的变换，即头动的**逆**。加上 AFNI 用 S/L/P
而 BIDS 用 RAS，一共两层符号翻转，只处理一层会得到全反的结果。
正确换算见 [`metrics.md`](metrics.md) §1.2。

**症状**：合成数据上「头向右移 3mm」，输出 `trans_x = -3`。

FD / enorm / maxdisp 基于差分或幅度，不受影响；只有逐点比对别家 confounds
时才会暴露。

### 3. 用众数推断期望 volume 数，在变长范式上全线误报

自定步调 / 被试按键结束的范式，run 长度本来就不等（实测 419–479 都正常）。
拿众数做严格比对会把大部分 run 判成 fail。

**现在的做法**：推断模式下用「同 task 最长 run」当参照，只在
`< truncated_fail_ratio`（默认 60%）时判 fail。显式配置 `expected_volumes`
才走严格比对。

### 4. tSNR 阈值必须跟着场强和分辨率走

默认 40/20 是按 3T 3mm EPI 定的。7T 1.2mm 各向同性的 tSNR 典型只有 10–25，
不改阈值就是全线误报。

| 数据 | tsnr_warn / fail 参考 |
|---|---|
| 3T 3.0mm | 40 / 20 |
| 3T 2.0mm | 25 / 12 |
| 7T 1.2mm iso | 15 / 8 |
| 7T 0.8mm | 8 / 4 |

### 5. 合成测试数据必须面内不对称

拿旋转对称的球/椭球做测试，刚体配准锁不住 in-plane rotation，会产生几度的
假 `rot_z`，经 50mm 弧长折算后 FD 虚高 5–7 倍，看起来像代码有 bug。

**症状**：完全没注入运动的合成 run，`fd_mean = 0.75mm`，且几乎全部来自
`rot_z`（`ptp(rot_z) ≈ 3.8°`），而三个平移的 `ptp` 都 < 0.12mm。

测试数据要做成椭球 + 几个偏心团块。修好后同样的数据 `fd_mean = 0.10mm`。

### 5a. run 间位移的参照点不能用配准基准卷

配准基准卷按 **min-outlier**（信号最干净）选，这对配准质量是对的，
**但它可能正好落在头位分布的边缘**，拿它当「run 间位移」的参照就会把结论判反。

**实测症状**（任务A 某被试，7 个 run）：被试在前两个 run 里往前挪了约 5mm 然后
安顿下来，run-03~07 头位彼此高度一致。基准卷选中了起始位置的 run-01，于是：

| run | vs run-01（错） | vs 组内中位姿态（对） |
|---|---|---|
| run-01 | 0.00 | **6.06** ← 真正的离群者 |
| run-03 | 5.88 | 0.81 |
| run-04 | **7.00 fail** | 1.02 |
| run-05 | **7.75 fail** | 1.88 |
| run-06 | **7.31 fail** | 1.36 |
| run-07 | **7.39 fail** | 1.40 |

4 个头位一致的 run 全判 fail，而唯一离群的 run-01 判 pass —— **结论完全反了**。

**修法**：两个角色分开。
- **配准基准卷**（喂给 `3dvolreg -base`）：仍用 min-outlier，配准质量优先
- **判定参照姿态**：用组内各 run 中位姿态的**逐分量中位数**，即「多数 run 共同的头位」

现在 tsv 里三列：
`betrun_disp_vs_median`（判定用）、`betrun_disp_vs_ref`（相对配准基准卷，诊断用）、
`betrun_disp_vs_prev`（相对采集顺序前一个，看什么时候挪的）。

**副产品观察**：同一项目的两个被试都是 run-01 离群 6–7mm、run-02 过渡、
之后稳定 —— 被试需要一两个 run 才安顿好。这类协议可以考虑让被试进磁体后多躺
几分钟再开始，或把第一个 run 当适应期。

### 5b. `base_volume` 必须记基准 run 的序号，不是自己的

`global` 模式下 8 个 run 共用 `run-06[303]` 这一个基准卷，但曾经把
`base_volume` 记成了**每个 run 自己的** min-outlier 序号（208/170/219/362…）。
计算用的是对的，错的只有落盘的溯源信息——**但这意味着别人拿
`base_run` + `base_volume` 复现不出同样的数**。

**症状**：按记录的基准卷独立重跑 3dvolreg，每个 TR 的 6 参数都差一个
近乎恒定的量（本例 Σ|Δ| ≈ 0.25，104 帧全都差同样多）。
**恒定偏移 = 基准卷不同**，随机抖动才是数值噪声，两者要分清。

修好后独立重跑与原输出的差异 < 5e-7 rad（只剩 tsv 写 6 位小数的舍入）。

现在 tsv 里有两列：`base_volume`（基准 run 的，复现用）和
`min_outlier_index`（该 run 自己的，诊断用）。

### 5c. 反向 PE 图和 func 的 affine 不一致时，`fslmerge` 会静默错位

topup 的标准流程是 `fslmerge -t AP_PA AP PA`，**但 fslmerge 只是按体素堆叠，
不看 affine，也不重采样**，输出沿用第一个输入的 header。所以只要 AP 和 PA
的体素网格朝向不同，合出来的 4D 文件就是错位的，而且没有任何报错。

`3dTcat` / `3dbucket` 同理。

**实测案例**（7T，见 `qc_shim.py`）：PA 的面内朝向比 AP 差 4.100°——
中心和法向都是精确复制的，只有 `dInPlaneRot` 没继承：AP 是技师定位时转出的
−4.100°，PA 是模板默认的整数 180.000°。

> **`dInPlaneRot` 与 PE 极性无关。** 相位编码反向由独立的极性标志
> （CSA image header 的 `PhaseEncodingDirectionPositive`）控制，不是靠把
> `dInPlaneRot` 加 180° 实现的。已在同机型另外两场数据上确认：func 与 reverse
> 的 `dInPlaneRot` 完全相同（都是 −0.300°），而 PE 极性确实相反（1 vs 0）。
> 所以**正确继承的 reverse，`dInPlaneRot` 应当和 func 一模一样**——本例中
> 应为 −4.100°。度量两者差异时对 180° 取模（矩形 FOV 转 180° 是同一个框）。

按体素强行对应的错位量：

| 距 FOV 中心 | 错位 | 体素数（1.2mm） |
|---|---|---|
| 25 mm | 1.79 mm | 1.49 |
| 50 mm（皮层典型） | 3.58 mm | 2.98 |
| 70 mm（大脑外缘） | 5.01 mm | 4.17 |
| 90 mm（FOV 边缘） | 6.44 mm | 5.37 |

中心几乎不动（0.06 mm），**错位随离中心距离线性增长**——这是这个 bug 的指纹：
中央区域看着还行，越往外缘越糟。topup 的模型里只有 PE 方向的位移、没有旋转项，
它会试图用 PE 位移去解释这个旋转，估出来的场因此被污染，后续所有
unwarp 和 BOLD→T1 配准跟着歪。

**怎么确认**：转完 BIDS 后直接比 affine，`fslhd` 的 `sto_xyz` 或

```bash
python3 -c "
import nibabel as nib, numpy as np
a=nib.load('..._dir-AP_..._bold.nii.gz').affine[:3,:3]
b=nib.load('..._dir-PA_..._epi.nii.gz').affine[:3,:3]
R=b@np.linalg.inv(a)
print('网格旋转 %.3f°'%np.degrees(np.arccos((np.trace(R)-1)/2)))"
```

不是 0 就有问题。**注意头动配准掩盖不了它**：这是采集几何差异，不是被试动了，
FD 再小也照样存在（案例里 Mean_FD 只有 0.08–0.13 mm）。

fMRIPrep / sdcflows 会先做刚体配准，比裸 topup 稳，但代价是多一次插值，
且不同版本处理程度不同。根治办法是改协议模板让 PA 正确继承 AP 的定位。

### 5d. 匀场框会冻结在协议模板值上，而 `lCoupleAdjVolTo` 不会告诉你

西门子协议里 adjustment volume 存的是一组绝对坐标。技师按被试头位重定 FOV 时，
**匀场框不一定跟着走**——跟不跟随取决于 UI 上的耦合状态，而那个状态
**没有忠实反映在 `sAdjData.lCoupleAdjVolTo` 里**。

**实测**（同一台 7T、同一份协议、6 个 session）：`lCoupleAdjVolTo` 全部为 1，
实际行为却分成两类。

**跟随组**（匀场框始终等于扫描框，偏移 0.00mm）：

- 另一批更早的 7T 数据：匀场框 `[6.69, −48.20, 15.62]`，等于扫描框。
- session A-1：中途技师重新定位，扫描框从 `[−5.44, 2.02, 15.76]` 移到
  `[−7.02, 15.82, 6.52]`（**16.68 mm、21.4°**），**匀场框同步跟到新位置**，
  全部 10 个 EPI series 偏移都是 0.00 mm。

**冻结组**（匀场框卡死在协议模板值 `[−5.44, 2.02, 15.76]`，法向 `[−0.046, 0.0546, 0.9975]`，
扫描框每场按被试头位重定）：

| session | 时间 | 扫描框中心 | 偏移 | 法向夹角 |
|---|---|---|---|---|
| A-2 | 第1天 场次2 | [−5.60, 4.25, 18.61] | 3.62 mm | 5.02° |
| B-1 | 第2天 场次1 | [0.42, −1.51, 16.41] | 6.87 mm | 7.66° |
| A-3 | 第2天 场次2 | [1.98, 7.27, 10.65] | 10.42 mm | 9.40° |
| B-2 | 第4天 场次1 | [3.59, 6.30, −0.16] | **18.80 mm** | **15.15°** |

**关键**：A-1 跟随，紧接着同一天的 A-2 就不跟随了，
之后再没跟随过。所以这不是机器或协议的固有行为，**是 A-1 与 A-2 之间某个设置被改动了**——
意味着它可以改回来。查那个时间点前后动过协议的什么。

**后果**：匀的是一块和成像区错开的脑组织。成像区内的 B0 均匀性没被优化，
眶额/颞极的 dropout 加重、EPI 畸变变大。7T 上磁化率效应正比于 B0，比 3T 敏感一倍多。
**这是采集期损失，后期补不回来**——fieldmap 能校几何畸变，但信号已经 dropout 的地方
没有信息可恢复。

**判据只能是直接比几何**（`qc_shim.py` 第一层），不能看 `lCoupleAdjVolTo`。
每场扫完立刻跑一遍，当场发现当场重定，比事后追悔便宜得多。

---

## B. 环境与运行

### 6. dcm2niix 不在 PATH，但 AFNI 自带一个

AFNI 装了就有 `~/abin/dcm2niix_afni`。dcm2bids 找的是 `dcm2niix`：

```bash
mkdir -p /tmp/bin && ln -sf ~/abin/dcm2niix_afni /tmp/bin/dcm2niix
export PATH=/tmp/bin:$PATH
```

比全局装一份干净，也不动用户环境。检查版本要 ≥ v1.0.20211006。

### 7. tSNR 是唯一会写大文件的步骤

`3dvolreg` 要真正输出配准后的数据才能算 tSNR。7T 470 vol × 150×150×124
单个 run 的中间文件约 2.6 GB，`--jobs 4` 峰值就是 10 GB 临时盘 + 相应内存。

大数据用 `--no-tsnr`，或把 `--jobs` 降到 2–3。脚本跑完会自动删中间文件
（`--keep-tmp` 可保留）。

### 8. `3dvolreg -prefix NULL` 才是真的不写盘

`-prefix NULL` 是 AFNI 的特殊值（大写），写成 `null` 或 `/dev/null` 都不行。

### 9. `-maxdisp1D foo.1D` 的 delta 文件叫 `foo.1D_delt`

不是 `_delta`。少一个 a。

### 10. `3dvolreg` 没有 `-quiet`

只有 `-verbose`。默认就不吵，不要加不存在的开关。

---

## C. 判读

### 11. 每个 run 都报「检出 N 个非稳态帧」通常是真的

实测某 3T 数据集 12/12 个 run 的第 1 帧都是坏的：

```
global_signal   前10: [428.0 519.3 525.6 ...]  稳态中位数 531.2   ← 低 19%
outlier_fraction 前10: [0.685 0.046 0.021 ...]                    ← 68.5% 体素离群
```

这是 dummy 没剔干净。确认后要么在预处理里丢掉，要么把
`nonsteady_expected` 设成实际数量止住告警——但**不要在没确认前就设**。

注意非稳态帧的信号可能偏**低**而不是偏高，取决于序列的预脉冲设置。

### 12. 「整体协议与队列多数派不同」是最该停下来看的告警

实测某数据集 sub-01 是 30 层 / 4.0mm，其余三人是 20 层 / 4.7mm——
被试间协议不一致，混进组分析会直接污染结果。

脚本区分两种情况：
- 该被试**所有** run 都与多数派不同 → 协议不同 → 直接 `fail`
- 只有部分 run 不同 → run 间参数漂移 → `warn`

### 12b. 不要把 pe_dir 放进跨被试一致性检查

不同 session 用相反的相位编码方向是**允许的**——只要该 session 内 func 与
fmap 相反就行。曾经把 `pe_dir` 放进 `consistency_fields`，结果两个方向不同
但各自都正确的 session，会有一个被判成「整体协议与队列多数派不同」→ fail。

真正该查的是配对关系，由 `check_pe_opposite` 负责。默认
`consistency_fields` 已移除 `pe_dir`。

### 13. FD 小但 maxdisp_range 大 = 慢漂移

不是急动。任务态 GLM 影响有限，静息态功能连接会受影响。看曲线判断，
不要只看单个数。

### 13b. 别人的脚本 mean_FD 比你小好几倍，先看是不是同一个被试

实测对比过一次：参考脚本 mean_FD 0.017–0.042，本工具 0.081–0.127，
差 3–5 倍，看起来像公式不同。实际是**两个不同被试**。

真要排除公式差异，直接把 3dvolreg 的原始 `.1D` 喂给 nipype 复算：

```python
from nipype.algorithms.confounds import FramewiseDisplacement
FramewiseDisplacement(in_file='mot.1D', parameter_source='AFNI').run()
```

本工具与它逐帧差 2e-16（见 [`metrics.md`](metrics.md) §2.1）。
所以 FD 对不上时，问题一定在运动估计或参考卷，不在公式。

### 14. DVARS 是在未做 MoCo 的原始数据上算的

因此绝对值比 fmriprep 的高很多，**只能用来在本批数据内部做相对比较**，
不要和文献阈值或 fmriprep 输出直接比。

### 15. fail 不等于排除

`fail` 的语义是「必须人看」，不是「已判死刑」。处置有好几种：整体排除 /
截掉尾部 / scrubbing / 保留但加协变量。决定要落到 `reason` 列或
`decision_log.md`，不许只记在脑子里。

---

## D. 与 bids-convert 的边界

### 16. qc_raw.py 不吃平铺 DICOM

必须先转 BIDS。bids-convert 的 `[7a] cleanup_aborted.py` 已经处理 1-vol 的
aborted run；qc_raw.py 抓的是它漏掉的「跑了一半才停」（例如 470 期望里只有
104 vol）以及所有质量指标。

### 17. 重复扫描的 run 编号要在转换阶段定死

实测数据里 DICOM 序列名 `run2` 实际是 run1（扫描时来不及改名），
`run2_2` 才是 run2。这种事 qc_raw.py 看不出来，只能在 bids-convert 的
`[3] 交互确认` 阶段问人，并写进 `decision_log.md`。

保留中断 run 时建议用 `acq-aborted` 标记而不是往后顺延 run 号，
这样 run 编号始终对应范式里的真实编号：

```
sub-01_task-A_acq-aborted_run-03_bold.nii.gz   ← 中断的那次
sub-01_task-A_run-03_bold.nii.gz               ← 重扫成功的
```
