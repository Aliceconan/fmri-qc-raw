---
name: fmri-qc-raw
description: 对 fMRI 数据做数据检查 / 质量检查 / 入库体检。输出每个 run 的 TR（volume）数、6 方向头动（run 内 + run 间）、采集是否齐全、参数是否一致、outlier/DVARS/tSNR，给出 pass/warn/fail 名单与带原因的 HTML 报告和文本表。只做 EPI-to-EPI 刚体估计，不做解剖或模板配准。当用户说「对这批数据做数据检查」「查一下数据质量」「数据能不能用」「头动大不大」「数据采集齐了没」「跑个 QC」，或需要在预处理前判断新到的数据是否合格时使用。主流程输入是 BIDS，拿到原始 DICOM 先用 bids-convert 转换。另含 qc_shim.py：查匀场框(shim box)与扫描框是否对齐、反向 PE 图是否继承了定位、跨 run 匀场结果是否一致、发射校准(reference amplitude)是否被手动覆盖——这一项只能读 DICOM 私有头，必须在转 BIDS 前跑，当用户问「匀场框和扫描框对齐没」「AP/PA 配不准」「topup/畸变校正有问题」「fieldmap 对不上」「reference amplitude 警告是什么」「翻转角不对」时也用本 skill。
---

# fmri-qc-raw · 原始数据入库体检

**SKILL_DIR（执行所有脚本时使用此路径）：**
```
SKILL_DIR=${CLAUDE_SKILL_DIR:-$HOME/.claude/skills/fmri-qc-raw}
```

数据刚到手时回答一个问题：**这批数据能不能用，不能用的是哪几个 run、为什么。**

## 先确认输入是什么

| 拿到的是 | 怎么做 |
|---|---|
| **BIDS 目录**（有 `sub-*/func/*_bold.nii.gz`）| 直接跑本 skill |
| **原始 DICOM**（平铺的 `.IMA` / `.dcm`）| 先跑一下 `qc_shim.py`（见 [0]，只有这一步能看匀场框，转完就没了），然后**用 `bids-convert` 转成 BIDS** 再跑主流程。`qc_raw.py` 不吃平铺 DICOM |
| 一堆散装 NIfTI，没按 BIDS 命名 | 用 `--pattern` 给 glob，但 run/task 解析会退化，完整性检查基本失效 |

在一个新文件夹里被要求「做数据检查」时，**第一步是 `ls` 看清楚这是哪种**，
不要假设。

不做配准（不配 T1、不配模板），只用 `3dvolreg` 做 EPI-to-EPI 刚体估计拿头动。
单被试 12 个 run（455 TR，3mm）在 6 并行下约 1 分钟。

## 这个 skill 检查什么

按「发现得越早越值钱」排序——**缺数据和协议漂移比头动更致命**，因为被试走了就补不回来：

| 层 | 查什么 | 典型问题 |
|---|---|---|
| **0 匀场**（DICOM 阶段）| 匀场框 vs 扫描框、跨 run 的定位与 shim 结果一致性 | 匀场框放歪、反向 PE 图没继承定位、AP/PA 不在同一 B0 场下 |
| **1 完整性** | run 齐不齐、volume 够不够、anat/fmap/json/events 在不在 | aborted run、扫到一半停机、sidecar 丢失 |
| **2 一致性** | TR / 层数 / 体素 / 矩阵 / TE / FA / PE 方向跨 run 跨被试是否一致 | 技师改了协议、某个被试用了另一套参数 |
| **3 头动** | 6 参数、FD、enorm、最大位移；**run 内**与**run 间**分开算 | 头动超阈、run 之间头位大幅偏移 |
| **4 信号** | outlier 比例、DVARS、tSNR、信号漂移、非稳态帧 | dummy 没剔干净、线圈异常、尖峰伪影 |

第 0 层单独由 `qc_shim.py` 做，**输入是 DICOM 不是 BIDS**，必须在转 BIDS 前跑
（信息转完就丢）；1–4 层由 `qc_raw.py` 做，输入是 BIDS。

## 工具依赖

| 工具 | 用途 |
|---|---|
| AFNI `3dvolreg` | 刚体 MoCo 估计 → 6 参数 + `-maxdisp1D` 最大体素位移 |
| AFNI `3dToutcount` | 逐 TR outlier 体素比例 + 选 min-outlier 基准卷 |
| AFNI `3dTstat` / `3dAutomask` / `3dmaskave` | 均值像、脑掩膜、全脑均值时间序列、tSNR |
| AFNI `3dTto1D` | DVARS |
| python3 | numpy、nibabel 必需；pyyaml（读配置）、matplotlib（画图）可选 |

## 目录结构

```
fmri-qc-raw/
  SKILL.md
  scripts/
    qc_raw.py            ← 主脚本：扫描 → 指标 → 判定 → 落盘
    qc_report.py         ← HTML 报告（被主脚本调用，也可单独重跑）
    qc_summary.py        ← 终端文本表（序列汇总 + 头动汇总），主脚本跑完自动写 summary.txt
    qc_shim.py           ← 匀场框 + 扫描框 + 匀场结果 + 发射校准检查。
                           【吃 DICOM，不吃 BIDS】独立于上面三个，
                           在转 BIDS 之前跑，见下面「[0] 匀场框检查」
  templates/
    qc_config.yaml       ← 阈值与期望值配置模板
  references/
    metrics.md           ← 每个指标的定义、单位、AFNI→BIDS 换算与验证记录
    pitfalls.md          ← 已知陷阱
```

## 与 bids-convert 联用

标准链路：

```
DICOM ──bids-convert──> bids/ ──qc_raw.py──> derivatives/qc-raw/
  │                       ↑                        │
  │                       └── 发现 aborted run ────┘
  │                           回到 bids-convert 的 [7a] 清理后重跑
  │
  └──qc_shim.py──> qc_shim.tsv        ← 只能在这一步做：匀场框不进 BIDS sidecar
```

- **`qc_raw.py` 必须先转 BIDS**，它吃 BIDS 布局，不吃平铺 DICOM。
- **`qc_shim.py` 反过来只吃 DICOM**，且必须在转换前跑——dcm2niix 不导出
  `sAdjData.*`，转完 BIDS 匀场框信息就没了。两条线互不依赖，可以并行。
- **撞名检测只在 DICOM 阶段有效**：`qc_shim.py` 第 [5] 层自动报同名序列，
  `qc_raw.py` 不查这个——转换后 BIDS 里若两个同名 series 已被合并或覆盖，
  就只剩一个文件，无从知道原本有两个。详见 pitfalls §17。
- bids-convert 的 `[6] 验证` / `[7a] cleanup_aborted.py` 已经查过文件数和 1-vol run；
  本 skill 是它的下一道门，查的是 **volume 数与多数派是否一致**（半截 run）
  以及 bids-convert 完全不看的头动与信号质量。
- 报出 `volume 数 X != 期望 Y` 时，回 bids-convert 的 [7a] 处理，
  清理完再跑一次。

## 工作流

```
[0] 匀场框检查 → qc_shim.py，**在转 BIDS 之前**、手里还是 DICOM 时跑（可选但便宜）
[1] 确认输入   → BIDS 根目录；确认是否已转完、是否需要先清 aborted run
[2] 配置阈值   → 复制 qc_config.yaml，按 3T/7T 调 tSNR 与体素相关阈值
[3] 跑         → qc_raw.py
[4] 读报告     → 只看 warn/fail，逐条判读
[5] 落决策     → 把处置写进 qc_summary 的 reason / decision_log
```

### [0] 匀场框检查（DICOM 阶段，转 BIDS 前）

**为什么必须在这一步**：匀场框（Siemens 的 adjustment volume）只存在于 DICOM
私有头 `(0029,1020)` 的 ASCCONV 文本段里。**dcm2niix 不会把它写进 BIDS sidecar**，
转完 BIDS 这个信息就永久丢了。所以它是本 skill 里唯一一个吃 DICOM 的脚本，
和 `qc_raw.py` 那条链完全独立。

```bash
# 最简：一场扫描的 DICOM 目录
python3 "$SKILL_DIR/scripts/qc_shim.py" --dicom /path/to/dicom_session

# 落盘 tsv + json
python3 "$SKILL_DIR/scripts/qc_shim.py" --dicom /path/to/dcm -o "$PROJECT_DIR/derivatives/qc-shim/sub-01"

# 核对本机型的坐标约定（换机器/换软件版本时值得跑一次）
python3 "$SKILL_DIR/scripts/qc_shim.py" --dicom /path/to/dcm --verify
```

只依赖 pydicom + numpy，不需要 AFNI。3982 个文件约 3.5 秒（每个 series 只读第一帧）。

查六层：

| 层 | 查什么 | 判定 |
|---|---|---|
| **series 内** | 匀场框 vs 扫描框：中心偏移、法向夹角、FOV 覆盖 | 偏移 > `--max-offset`(5mm) 或夹角 > `--max-angle`(3°) → `fail`；匀场框比扫描框小 → `warn` |
| **series 间·匀场框几何** | 匀场框几何分组，组内面内朝向是否统一 | 朝向差 > `--max-rot`(0.5°) → `warn` |
| **series 间·扫描框几何** | 同规格 EPI 的扫描框分组——func 与其反向 PE 图必须同组 | 分出多组 → `warn`（topup 的两个输入几何不同） |
| **series 间·匀场结果** | 实际 shim 电流 + 中心频率分组 | 同一定位下出现多套匀场结果 → `warn` |
| **series 间·发射校准** | reference amplitude 是否被手动钉死、同规格 EPI 内是否一致 | 带 `ucManualReferenceAmplitudeValid` → `warn`，并折算实际翻转角占标称的比例 |
| **series 名撞名** | 按 `SeriesDescription` 归组，同名即报 | 非定位像同名 → `warn`，并按帧数差异提示是中断+重扫还是两次完整采集 |

**别拿 `sAdjData.lCoupleAdjVolTo` 当判据。** 它看着像"匀场框是否耦合到 slice group"，
但实测 6 个 7T session **全部** `lCoupleAdjVolTo = 1`，实际行为却分成两类：2 场匀场框
严丝合缝跟着扫描框走（其中一场技师中途把 FOV 移了 16.68mm、转了 21.4°，匀场框同步跟到位，
偏移始终 0.00mm），另外 4 场匀场框冻结在协议模板值上纹丝不动，偏移 3.6 / 6.9 / 10.4 / 18.8 mm。
**这个字段的值和实际是否跟随无关**，只能当参考信息看（表里的 `耦合` 列照原样显示）。
唯一可靠的判据是直接比几何，也就是第一层在做的事——第一层永远要看，
不能因为 `耦合=1` 就跳过。

六层都会抓到东西，实测在 7T 数据上抓到过这些：

- **匀场框冻结在协议模板值，不跟 FOV 走**（最严重，见 pitfalls §5d）：连续 4 个
  session 的匀场框中心都是同一个 `[-5.44, 2.02, 15.76]`，而扫描框每次按被试头位
  重定，最差一场偏 18.8mm、法向差 15.15°。匀的是一块跟成像区错开的脑组织。
  同一台机器前一场还是跟随的，所以是设置被改动，不是固有行为。
- **func 与反向 PE 图定位不配对**（第 2b 层）：一场里所有 func 落一组、所有 reverse
  落另一组，差 0.87mm / 1.19°；另一场前 3 个 run 的 func 与 reverse 差 3.17mm，
  技师中途重新定位后才配上。这一层独立于匀场框分组——匀场框冻结时所有 EPI 会挤进
  同一个匀场几何组，只比匀场框看不出 func/reverse 之间的定位差异。
- **发射校准被手动钉死**（第 4 层，见 pitfalls §5f）：一个手动值 220.0 V 贯穿 8 场里的 7 场，
  但挂的位置会变——一场挂在 forward 的单个 run 上（造成**同一场内 run 之间**翻转角
  83.2% vs 100%），另外六场挂在所有 reverse 上。所以 func 和 reverse 两边都要查。
  判据是 `ucManualReferenceAmplitudeValid` 存不存在，
  **不能看 `bReferenceAmplitudeValid`**（两种情况都是 1）。
- **反向 PE 图的面内朝向没继承 func 的定位**：func 是 −4.100°（技师定位时转过），
  reverse 是整数 180.000°（模板默认值），对 180° 取模后净差 4.1°。
  中心和法向都是精确复制的，只有这一项没跟上。
  注意 `dInPlaneRot` **与 PE 极性无关**（极性另由 `PhaseEncodingDirectionPositive`
  区分），所以正确继承的 reverse 这一项应当与 func **完全相同**，不是差 180°——详见 pitfalls §5c。
  那个极性标志可以用来确认 func / reverse 确实反向（取值必须相反），
  但**不能拿它推绝对解剖方向**，会算反，见 pitfalls §5e。
- **AP 与 PA 用了两套匀场结果**：shim 电流差最多 189 DAC，f0 差 72 Hz。
  topup/SDC 的前提是 AP/PA 经历同一个 B0 场、只有 PE 极性相反，这个前提被破坏。

**别把这两条混为一谈**：4.1° 旋转对正方形 FOV 覆盖的组织影响很小，shim 差异更可能
只是「又独立 adjust 了一次」——间隔两小时的 7T，72 Hz 的 f0 漂移属正常范围。
从 DICOM 分辨不出到底是哪个原因，报告里也不下这个结论。

输出的 `零阶项参考量级: N voxel` 是拿 Δf0 除 `BandwidthPerPixelPhaseEncode` 折算的，
**不是实际位移**——扫描仪重设了 f0，常数项大部分会被 topup 吸收进场估计；
真正建模不了的是一/二阶 shim 的空间不均匀差异。别拿这个数去汇报。

产物 `qc_shim.tsv`（一行一个 series，24 列）+ `qc_shim.json`（含完整框几何与 shim 电流）。
有 `fail` 退出码为 1。

**厂商限制**：只有 Siemens 经典 DICOM/`.IMA` 能拿到完整几何。
Siemens XA 系列私有标签迁到 `(0021,xxxx)` 且导出时可能被剥掉；
GE `(0043,xxxx)` 一般没有独立 shim box（auto-shim 直接基于处方 FOV）；
Philips `(2005,xxxx)` 通常只能拿到模式不是几何。
读不到时脚本会明说原因并正常退出，不会假装通过。

### [1] 确认输入

启动时确认（不要默认假设）：

| 问题 | 影响 |
|---|---|
| 3T 还是 7T？体素多大？ | `tsnr_warn/fail` 必须跟着变，3T 的 40 拿到 7T 高分辨率上全会误报 |
| **有没有扫描协议表？** | 有就把每个 run 的期望 volume 数抄进 `expected_volumes`，走严格比对；没有就用比例推断 |
| 每个 task 几个 run、几个 fmap？ | `expected_runs` / `expected_fmap` |
| 采集时保留了几个 dummy？ | `nonsteady_expected`，不设的话每个 run 都会报非稳态帧 |

**要主动问有没有协议表。** 变长范式（run 长度本来就不等）下自动推断只能抓「明显中断」，
抓不到「比设计少 6 个 TR」这种——实测就有一个 run 期望 464 实际 470，
只有严格比对才看得见。

### [2] 配置阈值

```bash
cp "$SKILL_DIR/templates/qc_config.yaml" "$PROJECT_DIR/code/qc_config.yaml"
```

不给 `--config` 就用内置默认值（面向 3T 3mm EPI）。**7T 或高分辨率必须改 `tsnr_warn/fail`。**

### [3] 跑

```bash
# 最简：全部用默认阈值，输出到 bids/derivatives/qc-raw
python3 "$SKILL_DIR/scripts/qc_raw.py" --bids "$PROJECT_DIR/bids"

# 带配置 + 指定并行数
python3 "$SKILL_DIR/scripts/qc_raw.py" --bids "$PROJECT_DIR/bids" \
  --config "$PROJECT_DIR/code/qc_config.yaml" --jobs 8

# 只跑新到的被试
python3 "$SKILL_DIR/scripts/qc_raw.py" --bids "$PROJECT_DIR/bids" \
  --participant sub-05 sub-06

# 7T 大数据：关掉 tSNR（省一次 volreg 写盘）
python3 "$SKILL_DIR/scripts/qc_raw.py" --bids "$PROJECT_DIR/bids" --no-tsnr

# 非 BIDS 布局（自行给 glob）
python3 "$SKILL_DIR/scripts/qc_raw.py" --bids /data --pattern 'sub-*/*task*.nii.gz'

# 改完阈值重新判定，不重算指标（秒出；7T 全量重算要 20 分钟）
python3 "$SKILL_DIR/scripts/qc_raw.py" --rejudge "$PROJECT_DIR/bids/derivatives/qc-raw" \
  --config "$PROJECT_DIR/code/qc_config.yaml"
```

主要参数：

| 参数 | 说明 |
|---|---|
| `-o/--out` | 输出目录，默认 `<bids>/derivatives/qc-raw` |
| `--jobs` | 并行 run 数，默认 CPU 数的一半 |
| `--no-tsnr` | 跳过 tSNR（唯一需要写中间文件的步骤） |
| `--motion-base global\|per-run` | 默认 `global`：同几何组所有 run 配到同一基准卷，才能算 run 间头动 |
| `--keep-tmp` | 保留中间文件排查用 |
| `--rejudge OUTDIR` | 从已有 `qc_raw.json` 重新判定，不碰 AFNI。调阈值 / 补协议表用 |

退出码：有 fail 返回 1，可直接接 CI。

### [4] 产物

```
derivatives/qc-raw/
  qc_raw_runs.tsv        一行一个 run，53 列，含 status + reason
  qc_raw_subjects.tsv    一行一个 sub[_ses]，准入名单
  qc_raw.json            机器可读全量（含所用阈值，便于溯源）
  report.html            单文件报告，fail 的 run 默认展开曲线
  summary.txt            两张终端文本表：扫描序列汇总 + 头动汇总
  motion/
    sub-XX_task-YY_run-ZZ_desc-rawmotion_timeseries.tsv (+.json)
```

`motion/*.tsv` 用 **BIDS confounds 兼容列名**（`trans_x/y/z` mm、`rot_x/y/z` rad、
`framewise_displacement` mm），缺失一律 `n/a`。可以直接跟日后 fmriprep 的
`desc-confounds_timeseries.tsv` 对拍，验证预处理没把头动算错。

### [4b] 文本表

`summary.txt` 自动生成，也可单独重跑或换列名风格：

```bash
# 6 参数用 AFNI 原生列名 dS/dL/dP/Roll/Pitch/Yaw
python3 "$SKILL_DIR/scripts/qc_summary.py" --bids "$PROJECT_DIR/bids" \
  --qc "$PROJECT_DIR/bids/derivatives/qc-raw" --afni-names

# 只看序列表（不需要先跑 qc_raw.py，只读 header，秒出）
python3 "$SKILL_DIR/scripts/qc_summary.py" --bids "$PROJECT_DIR/bids"
```

序列表里的 `Phase_Dir` 是**从 affine + `PhaseEncodingDirection` 算出来的解剖方向**，
不是从序列名猜的。实测数据里名字叫 `AP_xxx_run` 的实际是 `P->A`——
Siemens 的 AP/PA 命名经常和实际方向相反，**以 Phase_Dir 为准**。
表里同时给出原始 `PE(BIDS)` 字段，方便自己核。

`summary.txt` 里还有一张 **func ↔ reverse-PE 配对检查表**。判定看
`头位差(净)` 那列（只含不受畸变污染的分量）；`PE轴(含畸变)` 那列**不参与判定**，
而且**不能当场漂移读**——func 与 reverse 的 PE 方向相反，同一个 B0 不均匀会把
两张图往相反方向推，所以它约等于 **2× 单向畸变量**。想看 run 内场漂移得去比
**相邻两个 reverse-PE**（PE 方向相同）。详见 `references/pitfalls.md` §5g，
那里有实测数字和命令。

### [5] 判读

报告只有三种状态，**fail 不等于排除**——它是「需要人看」的名单：

| 状态 | 含义 | 该做什么 |
|---|---|---|
| `pass` | 所有指标在门限内 | 放行 |
| `warn` | 有指标越线但不致命 | 看一眼曲线，多数可放行 |
| `fail` | volume 数不对、协议不同、或头动/信号严重超标 | 必须人工判读：整体排除？截尾？scrubbing？ |

处置策略写回 `reason` 列或项目的 `decision_log.md`，不要只记在脑子里。

## 关键实现事实（改代码前必读）

### FD 与 fmriprep 完全一致

`framewise_displacement` 用 Power FD（r=50mm），和
`nipype.algorithms.confounds.FramewiseDisplacement`（fmriprep 用的就是它）
**逐帧差 2e-16**，已在真实数据上实测。所以这一列可以直接和 fmriprep 的
同名列比；对不上就是运动估计或参考卷的问题，不是公式问题。

顺带：nipype 的 AFNI 参数归一化不做符号翻转，它的 6 参数和本工具在 z 和三个
旋转上差一个负号——但 FD 是 `|一阶差分|` 求和，不受影响。详见
[`references/metrics.md`](references/metrics.md) §2.1。

### 头动参数的坐标与符号

`3dvolreg -1Dfile` 输出 `roll pitch yaw dS dL dP`（旋转在前，degree；平移在后，mm），
并且给的是**把该卷配回基准所需的校正量，即头动的逆**。两层符号叠加后：

```
trans_x = +dL   trans_y = +dP   trans_z = -dS        (mm)
rot_x = -pitch  rot_y = -yaw    rot_z = -roll        (deg → rad)
```

已用已知位移/旋转的合成数据验证：平移误差 < 0.01 mm，旋转误差 < 0.05°，轴间串扰 ≈ 0。
输出的是**头动了多少**（`trans_x=+3` 表示头向右移动 3mm）。

### run 内 vs run 间必须分开

`--motion-base global` 下所有 run 配到同一个基准卷，基准卷可能在别的 run 里。因此：

- `maxdisp_max` = 相对**基准卷**的位移，**含 run 间偏移**，不能当 run 内头动
- `maxdisp_range` = 同 run 内极差，与基准卷无关，**这才是 run 内位移**（判定用它）
- `betrun_disp_vs_ref` / `_vs_prev` = 各 run 中位姿态之差，才是 run 间头动

FD 和 enorm 基于一阶差分，与基准卷无关，不受影响。

### 相位编码方向：查配对，不查全队列一致

**跨 session / 跨被试用相反的相位编码方向是允许的**，真正的硬性要求是
「同一 session 内，至少有一个 fmap 与 func 方向严格相反」——这是 blip/topup
矫正的前提。所以：

- `pe_dir` **不在**默认 `consistency_fields` 里（放进去会把正常的 session
  差异误报成协议不一致）
- 改由 `check_pe_opposite` 单独查，不满足直接 `fail`
- 同 session 内 func 各 run 方向不统一，也 `fail`

方向一律**从 affine + `PhaseEncodingDirection` 算**，不看序列名。
被试表输出 `func_pe` / `fmap_pe` / `pe_opposite_ok` 三列。

```
subject status func_pe fmap_pe  pe_opposite_ok  reason
 sub-01   pass    P->A    A->P            True
 sub-02   fail    P->A    P->A           False  没有与 func (P->A) 方向相反的 fmap，blip 矫正做不了
 sub-03   pass    A->P    P->A            True  ← 整体与 sub-01 相反，但自身配对正确，不报错
```

### 几何分组

分辨率/矩阵/PE 方向不同的 run 之间做刚体配准没有意义，脚本按
`(subject, session, dim, vox, pe_dir)` 分组，每组独立选基准卷，
run 间头动只在组内算。

### 期望值：配置优先，否则比例推断

`expected_volumes` 配了就**严格比对**（差 1 个 TR 都报）。没配则用「同 task 最长 run」
当参照做**比例判定**：< 60% 判 fail（疑似中断），< 95% 判 warn。

不能用众数做严格比对——自定步调 / 按键结束的范式 run 长度本来就不等
（实测 419–479 都是设计值），严格比对会全线误报。

`expected_volumes` 支持逐 run：

```yaml
expected_volumes:
  taskA: {"01": 464, "02": 470, "03": 419}   # 逐 run
  rest: 240                                  # 或该 task 统一
```

## 常见问题

### 每个 run 都报「检出 N 个非稳态帧」
多半是真的——采集时的 dummy 没剔。看 `motion/*.tsv` 的 `global_signal` 前几帧和
`outlier_fraction[0]`，如果第 1 帧全脑信号比稳态低 10%+ 且 outlier 占 50%+，
就是非稳态帧。确认后把 `nonsteady_expected` 设成实际数量，或在预处理里丢掉。

### tSNR 全线偏低
默认阈值按 3T 3mm EPI 定（典型 60–100）。7T 高分辨率、多回波、小体素的 tSNR
本来就低得多，必须改 `tsnr_warn/fail`，否则全是误报。

### 某个被试所有 run 都报「整体协议与队列多数派不同」
这是最值得停下来看的告警：该被试用了另一套采集参数（层数/体素/TE）。
混进组分析会直接污染结果。先确认是不是记录错误，再决定单独处理还是排除。

### FD 很小但 maxdisp_range 很大
慢漂移而非急动。任务态影响有限，静息态功能连接会受影响。看曲线判断。

### 合成数据上头动虚高
如果拿人造球体测试，面内旋转对称会让刚体配准锁不住 in-plane rotation，
造出几度的假旋转，经 50mm 弧长折算后 FD 虚高。测试数据必须做成不对称。

### 3dvolreg 报错找不到 base
`--motion-base global` 时基准卷写成 `path[idx]`。若某 run 的几何与组内其他 run
不同却被分到一组，说明分组键不够——检查 `dim_str` / `vox_str` / `pe_dir` 是否读到了。

## 按需参考

- [`scripts/qc_raw.py`](scripts/qc_raw.py) — 主脚本
- [`scripts/qc_report.py`](scripts/qc_report.py) — 报告生成
- [`scripts/qc_shim.py`](scripts/qc_shim.py) — 匀场框检查（吃 DICOM，转 BIDS 前跑）
- [`templates/qc_config.yaml`](templates/qc_config.yaml) — 阈值配置模板
- [`references/metrics.md`](references/metrics.md) — 指标定义、单位、换算与验证记录
- [`references/pitfalls.md`](references/pitfalls.md) — 已知陷阱
