# fmri-qc-raw

原始 fMRI 数据的入库体检工具。数据刚从扫描仪下来时回答一个问题：**这批数据能不能用，不能用的是哪几个 run、为什么。**

同时是一个 [Claude Code](https://claude.com/claude-code) skill——克隆到 `~/.claude/skills/fmri-qc-raw` 就能用自然语言驱动；不用 Claude Code 也可以直接当命令行脚本跑。

```bash
git clone https://github.com/Aliceconan/fmri-qc-raw ~/.claude/skills/fmri-qc-raw
```

> 命名说明：skill 和仓库叫 `fmri-qc-raw`，但脚本名（`qc_raw.py`）和 BIDS 输出路径（`derivatives/qc-raw/`）保持不变，这样已有的产物目录和 `--rejudge` 不会失效。

## 它查什么

按「发现得越早越值钱」排序。缺数据和协议漂移比头动更致命，因为被试走了就补不回来。

| 层 | 查什么 | 典型问题 |
|---|---|---|
| **0 匀场**（DICOM 阶段）| 匀场框 vs 扫描框、跨 run 的定位与 shim 结果一致性 | 匀场框放歪、反向 PE 图没继承定位、AP/PA 不在同一 B0 场下 |
| **1 完整性** | run 齐不齐、volume 够不够、anat/fmap/json/events 在不在 | aborted run、扫到一半停机、sidecar 丢失 |
| **2 一致性** | TR / 层数 / 体素 / 矩阵 / TE / FA / PE 方向跨 run 跨被试是否一致 | 技师改了协议、某个被试用了另一套参数 |
| **3 头动** | 6 参数、FD、enorm、最大位移；run 内与 run 间分开算 | 头动超阈、run 之间头位大幅偏移 |
| **4 信号** | outlier 比例、DVARS、tSNR、信号漂移、非稳态帧 | dummy 没剔干净、线圈异常、尖峰伪影 |

不做配准（不配 T1、不配模板），只用 `3dvolreg` 做 EPI-to-EPI 刚体估计拿头动。

## 两个入口，输入不一样

```
DICOM ──bids-convert──> bids/ ──qc_raw.py──> derivatives/qc-raw/
  │
  └──qc_shim.py──> qc_shim.tsv    ← 只能在这一步做：匀场框不进 BIDS sidecar
```

- **`qc_raw.py` 吃 BIDS**（第 1–4 层），不吃平铺 DICOM。
- **`qc_shim.py` 只吃 DICOM**（第 0 层），且必须在转换前跑。

### qc_raw.py

```bash
# 最简
python3 scripts/qc_raw.py --bids /path/to/bids

# 带配置 + 并行
python3 scripts/qc_raw.py --bids /path/to/bids --config qc_config.yaml --jobs 8

# 改完阈值重新判定，不重算指标（秒出）
python3 scripts/qc_raw.py --rejudge /path/to/bids/derivatives/qc-raw --config qc_config.yaml
```

产物：`qc_raw_runs.tsv`（一行一个 run，53 列，含 status + reason）、`qc_raw_subjects.tsv`（准入名单）、`report.html`、`summary.txt`，以及 BIDS confounds 兼容列名的逐 run 头动时间序列。有 fail 退出码为 1。

### qc_shim.py

匀场框（Siemens 的 adjustment volume）只存在于 DICOM 私有头 `(0029,1020)` 的 ASCCONV 文本段里，**dcm2niix 不会把它写进 BIDS sidecar**，转完 BIDS 这个信息就永久丢了。

```bash
python3 scripts/qc_shim.py --dicom /path/to/dicom_session
python3 scripts/qc_shim.py --dicom /path/to/dcm -o out/     # 落 tsv + json
python3 scripts/qc_shim.py --dicom /path/to/dcm --verify    # 核对本机型坐标约定
```

只依赖 pydicom + numpy，不需要 AFNI。~4000 个文件约 3.5 秒。

厂商限制：只有 Siemens 经典 DICOM/`.IMA` 能拿到完整几何。XA 系列私有标签迁到 `(0021,xxxx)` 且导出时可能被剥掉；GE 一般没有独立 shim box；Philips 通常只有匀场模式没有几何。读不到时会明说原因并正常退出，不会假装通过。

## 几个真实抓到过的问题

都在 [`references/pitfalls.md`](references/pitfalls.md) 里有完整记录和复现方法。

- **匀场框冻结在协议模板值上，不跟 FOV 走**。连续 4 个 session 的匀场框中心分毫不差，而扫描框每次按被试头位重定，最差一场偏 18.8 mm、法向差 15.15°——匀的是一块跟成像区错开的脑组织。而 `sAdjData.lCoupleAdjVolTo` 在跟随和不跟随的场次里都是 `1`，**这个字段不能当判据**。
- **反向 PE 图和 func 的 affine 不一致时，`fslmerge` 会静默错位**。fslmerge 只按体素堆叠，不看 affine 也不重采样。实测一例面内差 4.100°，中心几乎不动但错位随离中心距离线性增长，FOV 边缘 6.44 mm（5.4 voxel）。topup 的模型里没有旋转项，会拿 PE 位移去硬解释它。
- **run 间位移的参照点不能用配准基准卷**。基准卷按 min-outlier 选，可能正好落在头位分布边缘，于是 4 个头位一致的 run 全判 fail、唯一离群的那个判 pass——结论完全反了。

## 依赖

| | |
|---|---|
| `qc_raw.py` | AFNI（`3dvolreg` `3dToutcount` `3dTstat` `3dAutomask` `3dmaskave` `3dTto1D`）、numpy、nibabel；pyyaml 和 matplotlib 可选 |
| `qc_shim.py` | pydicom、numpy |

## 已验证的一致性

`framewise_displacement` 用 Power FD（r=50mm），与 `nipype.algorithms.confounds.FramewiseDisplacement`（fMRIPrep 用的就是它）在真实数据上**逐帧差 2e-16**。所以这一列可以直接和 fMRIPrep 的同名列对拍；对不上就是运动估计或参考卷的问题，不是公式问题。

换算细节、符号约定和全部验证记录见 [`references/metrics.md`](references/metrics.md)。

## 文档

- [`SKILL.md`](SKILL.md) — 完整用法、参数、判读指引
- [`references/metrics.md`](references/metrics.md) — 每个指标的定义、单位、AFNI→BIDS 换算与验证记录
- [`references/pitfalls.md`](references/pitfalls.md) — 已知陷阱，按「会不会让你得出错误结论」排序
- [`templates/qc_config.yaml`](templates/qc_config.yaml) — 阈值配置模板

文档是中文的。默认阈值面向 3T 3mm EPI，**7T 或高分辨率必须改 `tsnr_warn/fail`**，否则全是误报。

## License

MIT
