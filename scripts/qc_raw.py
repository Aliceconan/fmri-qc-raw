#!/usr/bin/env python3
"""
qc_raw.py — 原始 fMRI 数据入库体检

不做任何空间配准（不配到 T1、不配到模板），只做：
  1. 完整性 —— run 齐不齐、volume 够不够、anat/fmap/events 在不在
  2. 一致性 —— TR / 体素 / 矩阵 / TE / PE 方向跨 run 跨被试有没有漂
  3. 头动   —— 3dvolreg 提取 6 参数，run 内 + run 间
  4. 信号   —— outlier / DVARS / tSNR / 漂移 / 非稳态帧

依赖：AFNI (3dvolreg 3dToutcount 3dTstat 3dAutomask 3dmaskave 3dTto1D)
      python3 + numpy + nibabel (+ pyyaml 可选, matplotlib 可选)

用法：
    python3 qc_raw.py --bids /path/to/bids
    python3 qc_raw.py --bids BIDS -o BIDS/derivatives/qc-raw --jobs 8
    python3 qc_raw.py --bids BIDS --config code/qc_config.yaml
    python3 qc_raw.py --bids BIDS --participant sub-01 sub-02 --no-tsnr
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

try:
    import nibabel as nib
except ImportError:
    sys.exit("需要 nibabel：pip install nibabel")

# ---------------------------------------------------------------- 常量与默认阈值

# Power FD 把旋转折算成弧长时用的等效脑半径 (mm)
FD_HEAD_RADIUS = 50.0

DEFAULT_CFG = {
    # --- 逐 TR 判定阈值 ---
    "fd_thresh": 0.5,             # mm, Power FD 单帧阈值
    "outlier_thresh": 0.05,       # 单帧 outlier 体素比例阈值 (afni_proc 默认)
    # --- run 级判定 ---
    "pct_bad_warn": 10.0,         # 超阈值 TR 占比 (%)
    "pct_bad_fail": 25.0,
    "fd_mean_warn": 0.20,         # mm
    "fd_mean_fail": 0.50,
    # 单帧最大 FD（急动）。❌ 无文献标准，本工具自定。
    # 定在 1mm 是因为它约等于常见体素尺寸，且明显高于正常呼吸波动（~0.15mm）。
    "fd_max_warn": 1.0,
    "fd_max_fail": 3.0,
    "maxdisp_warn_vox": 1.0,      # run 内最大位移，单位=面内体素大小的倍数
    "maxdisp_fail_vox": 2.0,
    # tSNR 绝对下限。⚠ 文献中没有可引的门限（见 references/metrics.md §4.1），
    # 这两个数只当「明显坏掉」的兜底，必须按自己的场强/分辨率/线圈重标定。
    # 默认值面向 3T 3mm EPI。设为 None 可关掉绝对判定，只用相对判定。
    "tsnr_warn": 40.0,
    "tsnr_fail": 20.0,
    # tSNR 相对判定：占同协议队列中位数的百分比。这是 MRIQC 推荐的做法，
    # 不依赖任何绝对标准，换场强换序列都不用改。
    "tsnr_rel_warn": 70.0,
    "tsnr_rel_fail": 50.0,
    "tsnr_rel_min_n": 3,          # 同协议 run 数少于这个就没有可比基线，跳过
    "drift_warn_pct": 5.0,        # 全脑均值线性漂移占基线百分比
    # run 间位移 (mm)。没有文献标准。fail 线 6mm 的依据是实测：
    # 7T 1.2mm 数据上某 run 相对基准 9.05mm，该被试配准效果明显变差。
    "betrun_warn": 3.0,
    "betrun_fail": 6.0,
    # FOV 覆盖。阈值按实测标定：一批覆盖完好的 7T 数据贴边占比 ≤0.0003、
    # 组内覆盖差异 ≤2.3%，门限设在明显高于该基线的位置。
    "brain_edge_warn": 0.005,
    "brain_edge_fail": 0.02,
    "coverage_loss_warn": 5.0,    # %
    "coverage_loss_fail": 15.0,
    "nonsteady_expected": 0,      # 预期已剔除的 dummy 数；实测超出即告警
    # --- 完整性 ---
    "volume_tolerance": 0,        # 与期望 volume 数允许的偏差（仅 expected_volumes 显式配置时生效）
    # 推断模式下的截断判定：n_volumes / 同 task 最长 run 的比值
    "truncated_fail_ratio": 0.60,
    "truncated_warn_ratio": 0.95,
    "expected_volumes": {},       # {"task-xxx": 300} 显式期望；为空则用同 task 最长 run 推断
    "expected_runs": {},          # {"task-xxx": 6}
    "require_anat": True,
    "require_fmap": False,
    "expected_fmap": None,        # 协议里应有几个 fmap 文件；None = 不检查
    # --- 一致性（这些字段跨 run 必须一致，不一致即告警）---
    # 跨被试必须一致的字段。pe_dir 不在这里 —— 不同 session 用相反的相位编码
    # 方向是允许的，真正要保证的是「同 session 内 func 与 fmap 方向相反」，
    # 那个由 check_pe_opposite 单独查。
    "consistency_fields": ["tr", "n_slices", "vox_str", "dim_str", "te",
                           "flip_angle"],
    "check_pe_opposite": True,    # 查 func 与 fmap 的相位编码方向是否相反
    # func base 卷与配对 fmap 的头位差，单位=面内体素倍数。
    # 只统计不受畸变污染的分量（垂直 PE 的平移 + 绕 PE 轴的旋转）。
    "check_fmap_alignment": True,
    "fmap_disp_warn_vox": 1.0,
    "fmap_disp_fail_vox": 2.0,
    "presc_shift_warn": 2.0,      # mm，FOV 中心偏离组内中位规划
    "presc_angle_warn": 2.0,      # 度，层面朝向偏离
    "b0_drift_warn": 2.0,         # mm，整场 PE 轴 B0 漂移的告警线（自定）
    # --- 运行控制 ---
    "compute_tsnr": True,
    "motion_base": "global",      # global | per-run
}

JSON_REQUIRED_FIELDS = ["RepetitionTime", "EchoTime", "PhaseEncodingDirection",
                        "SliceTiming", "TotalReadoutTime", "FlipAngle"]

BIDS_ENTITY_RE = re.compile(r"(?:^|_)(sub|ses|task|acq|ce|dir|rec|run|echo|part)-([A-Za-z0-9]+)")


# ---------------------------------------------------------------- 小工具

def run_cmd(cmd, cwd=None, capture=True):
    """跑一条命令，失败抛异常，返回 stdout。"""
    p = subprocess.run(cmd, cwd=cwd, capture_output=capture, text=True)
    if p.returncode != 0:
        raise RuntimeError(
            f"命令失败 (exit {p.returncode}): {' '.join(map(str, cmd))}\n"
            f"stderr: {(p.stderr or '')[-2000:]}"
        )
    return p.stdout if capture else ""


def read_1d(path):
    """读 AFNI .1D 文件为 float 数组（跳过 # 注释行）。"""
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            rows.append([float(x) for x in line.split()])
    if not rows:
        return np.zeros((0, 0))
    return np.array(rows, dtype=float)


def parse_entities(name: str) -> dict:
    return {k: v for k, v in BIDS_ENTITY_RE.findall(name)}


def na(x, fmt="{:.4f}"):
    """输出到 tsv：缺失一律写 n/a（I2 契约要求）。"""
    if x is None:
        return "n/a"
    if isinstance(x, float) and (np.isnan(x) or np.isinf(x)):
        return "n/a"
    if isinstance(x, (bool, np.bool_)):
        return "True" if x else "False"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    if isinstance(x, float):
        return fmt.format(x)
    if isinstance(x, (list, tuple)):
        return ";".join(str(i) for i in x) if x else "n/a"
    s = str(x)
    return s if s else "n/a"      # 空字符串也必须写 n/a（I2 契约 6.3）


def pe_direction(affine, pe):
    """从 affine + BIDS PhaseEncodingDirection 算出解剖方向，如 'A->P'。

    序列名里的 AP/PA 经常和实际方向相反（Siemens 命名习惯），必须从 affine 算。
    """
    if not pe:
        return None
    axis = {"i": 0, "j": 1, "k": 2}.get(pe[0])
    if axis is None:
        return None
    v = affine[:3, axis] * (-1.0 if pe.endswith("-") else 1.0)
    n = np.linalg.norm(v)
    if n == 0:
        return None
    v = v / n
    k = int(np.argmax(np.abs(v)))
    pos, neg = ["R", "A", "S"][k], ["L", "P", "I"][k]
    return f"{neg}->{pos}" if v[k] > 0 else f"{pos}->{neg}"


def opposite_pe(a, b):
    """两个 'X->Y' 方向是否严格相反。"""
    if not a or not b:
        return None
    return a == "->".join(reversed(b.split("->")))


def check_afni():
    missing = [t for t in ("3dvolreg", "3dToutcount", "3dTstat", "3dAutomask",
                           "3dmaskave", "3dTto1D") if shutil.which(t) is None]
    if missing:
        sys.exit(f"以下 AFNI 程序不在 PATH 中：{', '.join(missing)}\n"
                 f"检查 ~/abin 是否在 PATH，或 source AFNI 的环境。")


# ---------------------------------------------------------------- 数据扫描

def scan_runs(bids_dir: Path, participants=None, tasks=None, pattern=None):
    """扫描 BIDS 目录下的 func BOLD 文件。pattern 给定时走通配符模式。"""
    if pattern:
        files = sorted(bids_dir.glob(pattern))
    else:
        files = sorted(bids_dir.glob("sub-*/**/func/*_bold.nii*"))
    runs = []
    for f in files:
        if f.name.startswith("."):
            continue
        ents = parse_entities(f.name)
        sub = ents.get("sub")
        if sub is None:
            continue
        if participants and f"sub-{sub}" not in participants and sub not in participants:
            continue
        if tasks and ents.get("task") not in tasks:
            continue
        runs.append({
            "path": str(f),
            "subject": f"sub-{sub}",
            "session": f"ses-{ents['ses']}" if "ses" in ents else "",
            "task": ents.get("task", "unknown"),
            "acq": ents.get("acq", ""),
            "dir": ents.get("dir", ""),
            "echo": ents.get("echo", ""),
            "run": ents.get("run", ""),
            "name": f.name.split(".nii")[0],
        })
    return runs


def sidecar_json(nii_path: Path):
    """就近读 sidecar；不做完整的 BIDS inheritance，只回溯到 bids 根。"""
    stem = nii_path.name.split(".nii")[0]
    cand = nii_path.parent / f"{stem}.json"
    if cand.exists():
        try:
            return json.loads(cand.read_text()), str(cand)
        except json.JSONDecodeError:
            return None, str(cand)
    return None, None


def read_header(run):
    """只读 header，不载入体数据。"""
    img = nib.load(run["path"])
    hdr = img.header
    shape = img.shape
    zooms = [float(z) for z in hdr.get_zooms()]
    nvol = int(shape[3]) if len(shape) > 3 else 1
    tr_hdr = float(zooms[3]) if len(zooms) > 3 else None

    js, js_path = sidecar_json(Path(run["path"]))
    js = js or {}
    tr = float(js["RepetitionTime"]) if "RepetitionTime" in js else tr_hdr
    st = js.get("SliceTiming")

    # 斜切角度：affine 的旋转部分相对于坐标轴的最大偏角
    aff = img.affine[:3, :3]
    norm = aff / np.linalg.norm(aff, axis=0, keepdims=True)
    obliq = float(np.degrees(np.arccos(np.clip(np.abs(norm).max(axis=0), -1, 1))).max())

    info = {
        "n_volumes": nvol,
        "tr": tr,
        "tr_header": tr_hdr,
        "dim_x": int(shape[0]), "dim_y": int(shape[1]), "dim_z": int(shape[2]),
        "n_slices": int(shape[2]),
        "vox_x": zooms[0], "vox_y": zooms[1], "vox_z": zooms[2],
        "vox_str": "x".join(f"{z:.2f}" for z in zooms[:3]),
        "dim_str": "x".join(str(s) for s in shape[:3]),
        "te": float(js["EchoTime"]) if "EchoTime" in js else None,
        "flip_angle": float(js["FlipAngle"]) if "FlipAngle" in js else None,
        "pe_dir": js.get("PhaseEncodingDirection"),
        "pe_anat": pe_direction(img.affine, js.get("PhaseEncodingDirection")),
        "readout": js.get("TotalReadoutTime"),
        "multiband": js.get("MultibandAccelerationFactor"),
        "slice_timing_n": len(st) if isinstance(st, list) else None,
        "obliquity_deg": obliq,
        # FOV 中心（世界坐标）与层面朝向矩阵：用来检测「层面重新规划」——
        # 技师中途重新定位会改变这两个，那是操作动作不是头动，必须分开报
        "fov_centre": (aff @ ((np.array(shape[:3]) - 1) / 2.0)
                       + img.affine[:3, 3]).tolist(),
        "orient_cols": (aff / np.linalg.norm(aff, axis=0)).T.tolist(),
        "duration_sec": (nvol * tr) if tr else None,
        "has_json": js_path is not None and bool(js),
        "json_missing_fields": [k for k in JSON_REQUIRED_FIELDS if k not in js],
        "voxel_inplane": float(np.mean(zooms[:2])),
        "scan_date": js.get("AcquisitionDateTime") or js.get("AcquisitionTime"),
        "series_desc": js.get("SeriesDescription"),
    }
    # TR 一致性自查：header 与 sidecar 打架是 dcm2niix/3drefit 的经典坑
    info["tr_mismatch"] = bool(
        tr and tr_hdr and abs(tr - tr_hdr) > 1e-3
    )
    return info


def scan_time_key(v):
    """把 AcquisitionDateTime / AcquisitionTime 转成可排序的数值。

    不能直接拿字符串比：dcm2niix 写出来的 AcquisitionTime 秒数**不补零**
    （见过 '15:18:3.447500'），字符串比较会把 :3 排到 :25 后面。
    解析失败时退回原字符串，至少不会崩。
    """
    if not v:
        return (1, "")
    t = str(v)
    if "T" in t:                      # ISO datetime，日期部分定长可直接比
        date, _, t = t.partition("T")
    else:
        date = ""
    try:
        parts = t.split(":")
        sec = float(parts[2]) if len(parts) > 2 else 0.0
        return (0, date, int(parts[0]) * 3600 + int(parts[1]) * 60 + sec)
    except (ValueError, IndexError):
        return (1, str(v))


# 解剖轴字母 -> BIDS/RAS 索引（0=x 左右, 1=y 前后, 2=z 上下）
ANAT_AXIS = {"R": 0, "L": 0, "A": 1, "P": 1, "S": 2, "I": 2}


def scan_fmaps(bids_dir: Path, subject, session):
    """收集该 session 的 reverse-PE fmap（*_epi.nii*）。"""
    sdir = bids_dir / subject / session if session else bids_dir / subject
    fdir = sdir / "fmap"
    if not fdir.exists():
        return []
    out = []
    for f in sorted(fdir.glob("*_epi.nii*")):
        js, _ = sidecar_json(f)
        js = js or {}
        try:
            aff = nib.load(str(f)).affine
        except Exception:                                     # noqa: BLE001
            continue
        out.append({
            "path": str(f), "name": f.name.split(".nii")[0],
            "run": parse_entities(f.name).get("run"),
            "pe_dir": js.get("PhaseEncodingDirection"),
            "pe_anat": pe_direction(aff, js.get("PhaseEncodingDirection")),
            "scan_date": js.get("AcquisitionDateTime") or js.get("AcquisitionTime"),
            "intended_for": js.get("IntendedFor"),
        })
    return out


def pair_fmap(run, fmaps):
    """给一个 func run 找配对的 fmap。

    优先级：BIDS 的 IntendedFor → run 实体号相同 → 采集时间最近。
    时间最近是最稳的兜底：fmap 少于 run 时（一个 fmap 服务多个 run），
    只有时间距离才有物理意义。
    """
    if not fmaps:
        return None, None
    for fm in fmaps:
        it = fm.get("intended_for")
        if not it:
            continue
        it = [it] if isinstance(it, str) else list(it)
        if any(Path(run["path"]).name in str(x) for x in it):
            return fm, "IntendedFor"
    if run.get("run"):
        same = [fm for fm in fmaps if fm.get("run") == run["run"]]
        if len(same) == 1:
            return same[0], "run 实体"
    tk = scan_time_key(run.get("scan_date"))
    if tk[0] == 0:
        cand = [fm for fm in fmaps if scan_time_key(fm.get("scan_date"))[0] == 0]
        if cand:
            return min(cand, key=lambda fm:
                       abs(scan_time_key(fm["scan_date"])[2] - tk[2])), "采集时间最近"
    return fmaps[0], "回退第一个"


def find_events(run):
    p = Path(run["path"])
    ev = p.parent / f"{run['name'].replace('_bold','')}_events.tsv"
    if not ev.exists():
        return None, None
    try:
        n = sum(1 for i, _ in enumerate(ev.open()) if i > 0)
    except OSError:
        n = None
    return str(ev), n


# ---------------------------------------------------------------- Phase A：逐 run 独立指标

def phase_a(run, tmp_root):
    """mask + outcount + 全脑均值 + 非稳态估计。不涉及 run 之间的关系。"""
    out = {"name": run["name"], "error": None}
    try:
        tmp = Path(tmp_root) / run["name"]
        tmp.mkdir(parents=True, exist_ok=True)
        src = run["path"]

        mean_f = tmp / "mean.nii.gz"
        mask_f = tmp / "mask.nii.gz"
        run_cmd(["3dTstat", "-mean", "-prefix", str(mean_f), src])
        run_cmd(["3dAutomask", "-q", "-clfrac", "0.5", "-prefix", str(mask_f), str(mean_f)])

        # 逐 TR outlier 体素比例
        oc = run_cmd(["3dToutcount", "-automask", "-fraction", "-legendre", src])
        outfrac = np.array([float(x) for x in oc.split()], dtype=float)

        # 逐 TR 全脑均值（稳态 / 漂移用）
        gs_txt = run_cmd(["3dmaskave", "-quiet", "-mask", str(mask_f), src])
        gs = np.array([float(x) for x in gs_txt.split()], dtype=float)

        # 掩膜统计：脑体素数 + 六个面的贴边占比。
        # 用来判断 run 间位移有没有把脑切出 FOV（位移大不一定切，切了一定糟）。
        mk = nib.load(str(mask_f)).get_fdata() > 0
        tot = int(mk.sum())
        faces = [mk[0], mk[-1], mk[:, 0], mk[:, -1], mk[:, :, 0], mk[:, :, -1]]
        edge = max(float(x.sum()) / tot for x in faces) if tot else None

        out.update({
            "mask": str(mask_f),
            "brain_voxels": tot,
            "brain_edge_frac": edge,
            "outlier_fraction": outfrac.tolist(),
            "global_signal": gs.tolist(),
            "min_outlier_index": int(np.argmin(outfrac[2:]) + 2) if len(outfrac) > 3
                                  else int(np.argmin(outfrac)) if len(outfrac) else 0,
            "outlier_frac_mean": float(np.mean(outfrac)) if len(outfrac) else None,
        })
    except Exception as e:                                    # noqa: BLE001
        out["error"] = f"{type(e).__name__}: {e}"
    return out


def estimate_nonsteady(gs, max_check=10):
    """前几帧显著偏离稳态基线 → 未剔除的 dummy / 未达稳态。"""
    gs = np.asarray(gs, dtype=float)
    if gs.size < 15:
        return 0
    base = gs[max_check:]
    med = float(np.median(base))
    mad = float(np.median(np.abs(base - med))) * 1.4826
    if mad <= 0:
        mad = float(np.std(base)) or 1e-9
    n = 0
    for i in range(min(max_check, gs.size)):
        if abs(gs[i] - med) > 5 * mad:
            n = i + 1
        else:
            break
    return n


def estimate_drift_pct(gs, skip=5):
    """线性漂移幅度占基线的百分比。"""
    gs = np.asarray(gs, dtype=float)[skip:]
    if gs.size < 10:
        return None
    x = np.arange(gs.size, dtype=float)
    slope, intercept = np.polyfit(x, gs, 1)
    base = float(np.mean(gs))
    if base == 0:
        return None
    return float(abs(slope * gs.size) / abs(base) * 100.0)


# ---------------------------------------------------------------- Phase C：3dvolreg

def afni_to_bids_motion(mp):
    """
    3dvolreg -1Dfile 列序：roll pitch yaw dS dL dP
        roll  = 绕 I-S 轴 (z) 旋转, 单位 degree (CCW)
        pitch = 绕 R-L 轴 (x) 旋转
        yaw   = 绕 A-P 轴 (y) 旋转
        dS/dL/dP = 向 Superior / Left / Posterior 的位移, 单位 mm

    两个必须同时处理的符号问题（已用已知位移的合成数据实测验证）：
      1) 轴向：AFNI 用 S/L/P，BIDS/RAS 用 R/A/S，L=-x、P=-y、S=+z
      2) 方向：3dvolreg 给的是「把该卷配回基准所需的校正量」，即头动的**逆**

    两者叠加后，头相对基准卷的实际位移为：
        trans_x = +dL,  trans_y = +dP,  trans_z = -dS
        rot_x   = -pitch, rot_y = -yaw, rot_z = -roll   (deg -> rad)

    即：本函数输出的是「头动了多少」（trans_x=+3 表示头向右移动 3mm），
    而不是 3dvolreg 原始数值。FD / enorm / maxdisp 都基于差分或幅度，
    与符号约定无关；只有逐点比对别家 confounds 时才需要留意这层换算。
    """
    roll, pitch, yaw, dS, dL, dP = (mp[:, i] for i in range(6))
    return np.column_stack([
        dL, dP, -dS,
        np.radians(-pitch), np.radians(-yaw), np.radians(-roll),
    ])


def fd_power(bids_mp, radius=FD_HEAD_RADIUS):
    """Power et al. 2012 FD：|Δtrans| 之和 + |Δrot(rad)|*r 之和，单位 mm。"""
    if bids_mp.shape[0] < 2:
        return np.zeros(bids_mp.shape[0])
    d = np.abs(np.diff(bids_mp, axis=0))
    d[:, 3:] *= radius
    return np.concatenate([[np.nan], d.sum(axis=1)])


def enorm_afni(mp):
    """afni_proc.py 的 motion_enorm：6 个原始参数一阶差分的欧氏范数。
    注意 AFNI 在这里直接把 degree 当 mm 用（不乘半径），
    标准 censor 阈值 0.2/0.3 就是按这个口径定的，不要换算。"""
    if mp.shape[0] < 2:
        return np.zeros(mp.shape[0])
    d = np.diff(mp, axis=0)
    return np.concatenate([[np.nan], np.sqrt((d ** 2).sum(axis=1))])


def phase_c(run, tmp_root, base_spec, mask, compute_tsnr):
    """3dvolreg（对 base_spec）+ DVARS + 可选 tSNR。"""
    out = {"name": run["name"], "error": None}
    try:
        tmp = Path(tmp_root) / run["name"]
        tmp.mkdir(parents=True, exist_ok=True)
        src = run["path"]
        mot_f = tmp / "motion.1D"
        maxd_f = tmp / "maxdisp.1D"
        volreg_prefix = str(tmp / "volreg.nii.gz") if compute_tsnr else "NULL"

        cmd = ["3dvolreg", "-zpad", "4",
               "-base", base_spec,
               "-prefix", volreg_prefix,
               "-1Dfile", str(mot_f),
               "-maxdisp1D", str(maxd_f),
               "-overwrite", src]
        run_cmd(cmd)

        mp = read_1d(mot_f)                       # (T, 6) AFNI 口径
        maxdisp = read_1d(maxd_f).ravel()         # 每个 TR 相对 base 的最大体素位移
        delt_p = Path(str(maxd_f) + "_delt")
        maxdisp_delt = read_1d(delt_p).ravel() if delt_p.exists() else np.array([])

        out["motion_afni"] = mp.tolist()
        out["maxdisp"] = maxdisp.tolist()
        out["maxdisp_delt"] = maxdisp_delt.tolist()

        # DVARS（在原始数据上算，未做 MoCo，故偏高，仅作相对比较）
        dv_f = tmp / "dvars.1D"
        try:
            run_cmd(["3dTto1D", "-input", src, "-mask", mask,
                     "-method", "dvars", "-prefix", str(dv_f)])
            out["dvars"] = read_1d(dv_f).ravel().tolist()
        except RuntimeError as e:
            out["dvars"] = []
            out["dvars_error"] = str(e)[:300]

        # tSNR：必须在 MoCo 之后算才有意义
        if compute_tsnr and Path(volreg_prefix).exists():
            tsnr_f = tmp / "tsnr.nii.gz"
            run_cmd(["3dTstat", "-cvarinv", "-prefix", str(tsnr_f), volreg_prefix])
            t = nib.load(str(tsnr_f)).get_fdata()
            m = nib.load(mask).get_fdata() > 0
            vals = t[m]
            vals = vals[np.isfinite(vals) & (vals > 0)]
            out["tsnr_median"] = float(np.median(vals)) if vals.size else None
            os.remove(volreg_prefix)
            os.remove(tsnr_f)
    except Exception as e:                                    # noqa: BLE001
        out["error"] = f"{type(e).__name__}: {e}"
    return out


# ---------------------------------------------------------------- Phase F：func vs fmap

def phase_f(run, fmap, base_idx, tmp_root):
    """把配对 fmap 的均值刚体配到 func run 的基准卷，拆出不受畸变污染的分量。

    ⚠ func 与 fmap 相位编码方向相反 → 几何畸变相反。刚体配准会把「畸变差异」
    算进位移，**PE 轴上的平移不能当头动看**。6 个自由度里只有 3 个干净：
      - 垂直于 PE 的两个平移（y 方向的位移场造不出 x/z 位移）
      - 绕 PE 轴的旋转（该旋转只动另外两轴）
    受污染的：PE 轴平移，以及绕另外两轴的旋转
      （y 位移随 z 变化 ≈ 绕 x 转；随 x 变化 ≈ 绕 z 转）
    代价：func↔fmap 之间沿 PE 的真实平移、绕另两轴的真实旋转抓不到 ——
    单个时刻只有一种极性，物理上无法分离。
    """
    out = {"name": run["name"], "error": None}
    try:
        tmp = Path(tmp_root) / (run["name"] + "_fmap")
        tmp.mkdir(parents=True, exist_ok=True)
        fm_mean = tmp / "fmap_mean.nii.gz"
        fbase = tmp / "func_base.nii.gz"
        mot = tmp / "f2b.1D"
        run_cmd(["3dTstat", "-mean", "-prefix", str(fm_mean),
                 "-overwrite", fmap["path"]])
        run_cmd(["3dTcat", "-prefix", str(fbase), "-overwrite",
                 f"{run['path']}[{base_idx}]"])
        run_cmd(["3dvolreg", "-zpad", "4", "-base", str(fbase), "-prefix", "NULL",
                 "-1Dfile", str(mot), "-overwrite", str(fm_mean)])
        mp = afni_to_bids_motion(read_1d(mot))
        if mp.shape[0] == 0:
            raise RuntimeError("3dvolreg 未产出参数")
        d = np.abs(mp[0])

        pe = run.get("pe_anat") or fmap.get("pe_anat")
        ax = ANAT_AXIS.get(pe.split("->")[0]) if pe else 1     # 默认按 y
        perp = [i for i in range(3) if i != ax]

        out.update({
            "fmap_name": fmap["name"],
            "fmap_pe_dir": fmap.get("pe_dir"),
            "fmap_pe_anat": fmap.get("pe_anat"),
            "pe_opposite_pair": opposite_pe(run.get("pe_anat"),
                                            fmap.get("pe_anat")),
            "fmap_disp_clean": float(d[perp[0]] + d[perp[1]]
                                     + d[3 + ax] * FD_HEAD_RADIUS),
            "fmap_disp_pe": float(d[ax]),
            "fmap_disp_full": float(d[:3].sum() + d[3:].sum() * FD_HEAD_RADIUS),
            "fmap_params": [float(x) for x in mp[0]],
            "fmap_pe_axis": "xyz"[ax],
        })
    except Exception as e:                                    # noqa: BLE001
        out["error"] = f"{type(e).__name__}: {e}"
    return out


def phase_fmap_traj(fmaps, tmp_root):
    """把同 session 的 fmap 互相配准（同极性，畸变相同 → 纯头动轨迹）。

    用来估计 B0 漂移：func 轨迹含 +畸变漂移，fmap 轨迹含 −畸变漂移，
    两者之差的一半就是畸变漂移的量级估计。
    """
    if len(fmaps) < 2:
        return {}
    tmp = Path(tmp_root) / "fmap_traj"
    tmp.mkdir(parents=True, exist_ok=True)
    means = {}
    for fm in fmaps:
        m = tmp / f"{fm['name']}_mean.nii.gz"
        run_cmd(["3dTstat", "-mean", "-prefix", str(m), "-overwrite", fm["path"]])
        means[fm["name"]] = str(m)
    ref = fmaps[0]["name"]
    traj = {}
    for fm in fmaps:
        if fm["name"] == ref:
            traj[fm["name"]] = np.zeros(6)
            continue
        f1 = tmp / f"{fm['name']}_to_ref.1D"
        run_cmd(["3dvolreg", "-zpad", "4", "-base", means[ref], "-prefix", "NULL",
                 "-1Dfile", str(f1), "-overwrite", means[fm["name"]]])
        mp = afni_to_bids_motion(read_1d(f1))
        traj[fm["name"]] = mp[0] if mp.shape[0] else np.zeros(6)
    return traj


# ---------------------------------------------------------------- 判定

def judge_run(r, cfg):
    """返回 (status, [reasons])。fail 优先于 warn。"""
    fails, warns = [], []

    if r.get("qc_error"):
        return "fail", [f"处理失败: {r['qc_error'][:120]}"]

    # --- 完整性 ---
    # 显式配置了期望值就严格比对；靠众数推断时不能严格比对——很多范式本来就是
    # 变长的（自定步调 / 被试按键结束），严格比对会全线误报。推断模式下只认
    # 「明显短于同 task 最长 run」= 中断，其余只提示。
    exp, src = r.get("expected_volumes"), r.get("expected_volumes_source")
    if exp is not None:
        diff = r["n_volumes"] - exp
        if src == "config":
            if abs(diff) > cfg["volume_tolerance"]:
                fails.append(f"volume 数 {r['n_volumes']} != 期望 {exp} ({diff:+d})")
        else:
            ratio = r["n_volumes"] / exp if exp else 1.0
            if ratio < cfg["truncated_fail_ratio"]:
                fails.append(f"volume 数 {r['n_volumes']} 仅为同 task 最长 run "
                             f"({exp}) 的 {ratio:.0%}，疑似中断")
            elif ratio < cfg["truncated_warn_ratio"]:
                warns.append(f"volume 数 {r['n_volumes']} 短于同 task 最长 run "
                             f"({exp}) 的 {ratio:.0%}")
    if not r.get("has_json"):
        warns.append("缺 sidecar json")
    elif r.get("json_missing_fields"):
        warns.append("json 缺字段: " + ",".join(r["json_missing_fields"]))
    if r.get("tr_mismatch"):
        warns.append(f"TR 不一致: header {r['tr_header']} vs json {r['tr']}")
    if r.get("n_volumes", 0) < 10:
        fails.append(f"volume 数过少 ({r['n_volumes']})，疑似中断 run")

    # --- 一致性（跨 run 众数比对，在上层填好 inconsistent_fields）---
    if r.get("inconsistent_fields"):
        warns.append("参数与多数 run 不一致: " + ",".join(r["inconsistent_fields"]))

    # --- 层面规划（操作动作，非头动）---
    ps, pa = r.get("presc_shift_mm"), r.get("presc_angle_deg")
    if ps is not None and pa is not None:
        if ps > cfg["presc_shift_warn"] or pa > cfg["presc_angle_warn"]:
            warns.append(f"层面规划与组内多数不同（中心偏 {ps:.1f}mm，朝向偏 "
                         f"{pa:.2f}°）—— 中途重新定位过，不是头动")

    # --- 头动 ---
    pct = r.get("pct_fd_gt_thresh")
    if pct is not None:
        if pct > cfg["pct_bad_fail"]:
            fails.append(f"FD>{cfg['fd_thresh']}mm 的 TR 占 {pct:.1f}%")
        elif pct > cfg["pct_bad_warn"]:
            warns.append(f"FD>{cfg['fd_thresh']}mm 的 TR 占 {pct:.1f}%")
    # 单帧急动：均值和占比都抓不到它 —— 一次 2mm 的猛动只贡献 1/450 的占比，
    # 但会污染后续多个 TR（spin-history），必须单独判。实测漏过一次：
    # run-07 TR 212 FD=1.97mm、outlier 从 2.7% 跳到 20% 且不回落，
    # 而 fd_mean 0.140、超阈占比 0.93%，两道全过。
    fdx = r.get("fd_max")
    if fdx is not None:
        n_sp = r.get("n_fd_spikes") or 0
        where = f"（TR {r.get('fd_max_tr')}）" if r.get("fd_max_tr") else ""
        if fdx > cfg["fd_max_fail"]:
            fails.append(f"单帧最大 FD {fdx:.2f}mm{where}")
        elif fdx > cfg["fd_max_warn"]:
            warns.append(f"单帧最大 FD {fdx:.2f}mm{where}"
                         + (f"，共 {n_sp} 次急动" if n_sp > 1 else ""))

    fdm = r.get("fd_mean")
    if fdm is not None:
        if fdm > cfg["fd_mean_fail"]:
            fails.append(f"平均 FD {fdm:.3f}mm")
        elif fdm > cfg["fd_mean_warn"]:
            warns.append(f"平均 FD {fdm:.3f}mm")
    md, vox = r.get("maxdisp_range"), r.get("voxel_inplane")
    if md is not None and vox:
        if md > cfg["maxdisp_fail_vox"] * vox:
            fails.append(f"run 内位移范围 {md:.2f}mm (>{cfg['maxdisp_fail_vox']} 体素)")
        elif md > cfg["maxdisp_warn_vox"] * vox:
            warns.append(f"run 内位移范围 {md:.2f}mm (>{cfg['maxdisp_warn_vox']} 体素)")
    # 判定用「相对组内中位姿态」的位移；vs_ref（相对配准基准卷）只作诊断输出，
    # 因为基准卷是按信号质量选的，可能本身就在头位分布边缘。
    br = r.get("betrun_disp_vs_median")
    if br is None:
        br = r.get("betrun_disp_vs_ref")
    if br is not None:
        if br > cfg["betrun_fail"]:
            fails.append(f"头位偏离本 session 中位姿态 {br:.2f}mm，配准可能失败")
        elif br > cfg["betrun_warn"]:
            warns.append(f"头位偏离本 session 中位姿态 {br:.2f}mm")

    # --- FOV 覆盖：run 间位移的一种具体后果，独立于位移量本身 ---
    ef = r.get("brain_edge_frac")
    if ef is not None:
        if ef > cfg["brain_edge_fail"]:
            fails.append(f"{ef:.2%} 的脑体素贴在 FOV 边界上，脑被切了")
        elif ef > cfg["brain_edge_warn"]:
            warns.append(f"{ef:.2%} 的脑体素贴在 FOV 边界上")
    cov = r.get("coverage_rel_pct")
    if cov is not None:
        loss = 100.0 - cov
        if loss > cfg["coverage_loss_fail"]:
            fails.append(f"脑覆盖比组内最好的 run 少 {loss:.1f}%")
        elif loss > cfg["coverage_loss_warn"]:
            warns.append(f"脑覆盖比组内最好的 run 少 {loss:.1f}%")

    # --- func 与配对 fmap 的头位差（blip 矫正的前提）---
    fdc, vox = r.get("fmap_disp_clean"), r.get("voxel_inplane")
    if fdc is not None and vox:
        pe = r.get("fmap_disp_pe")
        extra = f"；PE 轴另有 {pe:.2f}mm（头动与畸变混合，不计入判定）" if pe else ""
        if fdc > cfg["fmap_disp_fail_vox"] * vox:
            fails.append(f"func 与 fmap 头位差 {fdc:.2f}mm "
                         f"(>{cfg['fmap_disp_fail_vox']} 体素)，blip 矫正会失败{extra}")
        elif fdc > cfg["fmap_disp_warn_vox"] * vox:
            warns.append(f"func 与 fmap 头位差 {fdc:.2f}mm "
                         f"(>{cfg['fmap_disp_warn_vox']} 体素){extra}")
    if r.get("pe_opposite_pair") is False:
        fails.append(f"该 run ({r.get('pe_anat')}) 与配对 fmap "
                     f"({r.get('fmap_pe_anat')}) 方向不相反，blip 矫正做不了")
    if r.get("fmap_error"):
        warns.append(f"fmap 比较失败: {r['fmap_error'][:80]}")

    # --- 信号 ---
    po = r.get("pct_out_gt_thresh")
    if po is not None:
        if po > cfg["pct_bad_fail"]:
            fails.append(f"outlier>{cfg['outlier_thresh']} 的 TR 占 {po:.1f}%")
        elif po > cfg["pct_bad_warn"]:
            warns.append(f"outlier>{cfg['outlier_thresh']} 的 TR 占 {po:.1f}%")
    # tSNR 两道：绝对下限（明显坏掉）+ 队列相对（真正可靠的那道）
    ts = r.get("tsnr_median")
    if ts is not None and cfg.get("tsnr_fail") and ts < cfg["tsnr_fail"]:
        fails.append(f"tSNR 中位数 {ts:.1f}，低于绝对下限 {cfg['tsnr_fail']}")
    elif ts is not None and cfg.get("tsnr_warn") and ts < cfg["tsnr_warn"]:
        warns.append(f"tSNR 中位数 {ts:.1f}，低于绝对下限 {cfg['tsnr_warn']}")
    tr_ = r.get("tsnr_rel_pct")
    if tr_ is not None:
        if tr_ < cfg["tsnr_rel_fail"]:
            fails.append(f"tSNR 只有同协议队列中位数的 {tr_:.0f}%")
        elif tr_ < cfg["tsnr_rel_warn"]:
            warns.append(f"tSNR 只有同协议队列中位数的 {tr_:.0f}%")
    dr = r.get("gs_drift_pct")
    if dr is not None and dr > cfg["drift_warn_pct"]:
        warns.append(f"信号漂移 {dr:.1f}%")
    ns = r.get("n_nonsteady_est")
    if ns is not None and ns > cfg["nonsteady_expected"]:
        warns.append(f"检出 {ns} 个非稳态帧（预期 {cfg['nonsteady_expected']}）")

    if fails:
        return "fail", fails + warns
    if warns:
        return "warn", warns
    return "pass", []


# ---------------------------------------------------------------- 主流程

def main():
    ap = argparse.ArgumentParser(
        description="原始 fMRI 数据入库体检（完整性 / 一致性 / 头动 / 信号）")
    ap.add_argument("--bids", help="BIDS 根目录（--rejudge 时不需要）")
    ap.add_argument("--rejudge", metavar="OUTDIR",
                    help="从已有输出目录的 qc_raw.json 重新判定，不重算指标。"
                         "改完阈值用这个，秒出结果")
    ap.add_argument("-o", "--out", default=None,
                    help="输出目录，默认 <bids>/derivatives/qc-raw")
    ap.add_argument("--config", default=None, help="qc_config.yaml（覆盖默认阈值）")
    ap.add_argument("--participant", nargs="+", default=None, help="只跑指定被试")
    ap.add_argument("--task", nargs="+", default=None, help="只跑指定 task")
    ap.add_argument("--pattern", default=None,
                    help="非 BIDS 布局时的 glob，如 'sub-*/*bold.nii.gz'")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 4) // 2))
    ap.add_argument("--tsnr", dest="tsnr", action="store_true", default=None)
    ap.add_argument("--no-tsnr", dest="tsnr", action="store_false",
                    help="跳过 tSNR（省一次 volreg 写盘，7T 大数据建议关）")
    ap.add_argument("--motion-base", choices=["global", "per-run"], default=None,
                    help="global=同几何组所有 run 配到同一基准卷（默认，可算 run 间头动）")
    ap.add_argument("--keep-tmp", action="store_true")
    ap.add_argument("--no-report", action="store_true", help="不生成 HTML")
    args = ap.parse_args()

    if args.rejudge:
        return rejudge(Path(args.rejudge), args.config, not args.no_report)
    if not args.bids:
        ap.error("需要 --bids（或用 --rejudge 重判已有结果）")

    check_afni()

    cfg = dict(DEFAULT_CFG)
    if args.config:
        import yaml
        user = yaml.safe_load(Path(args.config).read_text()) or {}
        cfg.update(user.get("qc_raw", user))
    if args.tsnr is not None:
        cfg["compute_tsnr"] = args.tsnr
    if args.motion_base:
        cfg["motion_base"] = args.motion_base

    bids = Path(args.bids).resolve()
    out_dir = Path(args.out).resolve() if args.out else bids / "derivatives" / "qc-raw"
    (out_dir / "motion").mkdir(parents=True, exist_ok=True)

    runs = scan_runs(bids, args.participant, args.task, args.pattern)
    if not runs:
        sys.exit(f"在 {bids} 下没找到 BOLD 文件（试试 --pattern）")
    print(f"[qc-raw] 找到 {len(runs)} 个 BOLD run，{args.jobs} 并行", flush=True)

    # ---- header ----
    for r in runs:
        try:
            r.update(read_header(r))
            ev, nev = find_events(r)
            r["events_path"], r["n_events"] = ev, nev
        except Exception as e:                                # noqa: BLE001
            r["qc_error"] = f"读 header 失败: {e}"

    tmp_root = tempfile.mkdtemp(prefix="qcraw_", dir=str(out_dir))
    t0 = time.time()

    # ---- Phase A ----
    print("[qc-raw] Phase A: mask / outlier / 全脑均值 ...", flush=True)
    a_res = {}
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(phase_a, r, tmp_root): r["name"] for r in runs
                if not r.get("qc_error")}
        for i, fut in enumerate(as_completed(futs), 1):
            res = fut.result()
            a_res[res["name"]] = res
            print(f"  [{i}/{len(futs)}] {res['name']}"
                  f"{'  ERROR' if res['error'] else ''}", flush=True)
    for r in runs:
        a = a_res.get(r["name"])
        if a is None:
            continue
        if a["error"]:
            r["qc_error"] = a["error"]
            continue
        r["_mask"] = a["mask"]
        r["_outfrac"] = np.array(a["outlier_fraction"])
        r["_gs"] = np.array(a["global_signal"])
        r["_min_out_idx"] = a["min_outlier_index"]
        r["outlier_frac_mean"] = a["outlier_frac_mean"]
        r["brain_voxels"] = a.get("brain_voxels")
        r["brain_edge_frac"] = a.get("brain_edge_frac")

    # ---- 几何分组 + 选基准卷 ----
    # 分辨率/矩阵/PE 不同的 run 之间做刚体配准没有意义，必须分组
    groups = defaultdict(list)
    for r in runs:
        if r.get("qc_error"):
            continue
        key = (r["subject"], r["session"], r.get("dim_str"), r.get("vox_str"),
               r.get("pe_dir"))
        groups[key].append(r)

    for key, grp in groups.items():
        if cfg["motion_base"] == "per-run":
            for r in grp:
                r["_base_spec"] = f"{r['path']}[{r['_min_out_idx']}]"
                r["_base_run"] = r["name"]
                r["_base_vol"] = r["_min_out_idx"]
            continue
        # 全局基准：组内平均 outlier 最低的 run 的 min-outlier 卷
        ref = min(grp, key=lambda r: r.get("outlier_frac_mean") or 1e9)
        spec = f"{ref['path']}[{ref['_min_out_idx']}]"
        for r in grp:
            r["_base_spec"] = spec
            r["_base_run"] = ref["name"]
            # 注意：是 ref 的 min-outlier 序号，不是 r 自己的。
            # 记错了会导致别人拿 base_run+base_volume 复现不出同样的结果
            r["_base_vol"] = ref["_min_out_idx"]

    # ---- Phase C ----
    print("[qc-raw] Phase C: 3dvolreg 提取头动"
          f"{' + tSNR' if cfg['compute_tsnr'] else ''} ...", flush=True)
    c_res = {}
    todo = [r for r in runs if not r.get("qc_error")]
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(phase_c, r, tmp_root, r["_base_spec"], r["_mask"],
                          cfg["compute_tsnr"]): r["name"] for r in todo}
        for i, fut in enumerate(as_completed(futs), 1):
            res = fut.result()
            c_res[res["name"]] = res
            print(f"  [{i}/{len(futs)}] {res['name']}"
                  f"{'  ERROR' if res['error'] else ''}", flush=True)

    # ---- Phase F：func base 卷 vs 配对 fmap ----
    if cfg.get("check_fmap_alignment"):
        print("[qc-raw] Phase F: func 与配对 fmap 的头位比较 ...", flush=True)
        fm_cache = {}
        jobs = []
        for r in todo:
            key = (r["subject"], r["session"])
            if key not in fm_cache:
                fm_cache[key] = scan_fmaps(bids, r["subject"], r["session"])
            fm, how = pair_fmap(r, fm_cache[key])
            if fm is None:
                continue
            r["fmap_pair_method"] = how
            jobs.append((r, fm))
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futs = {ex.submit(phase_f, r, fm, r["_min_out_idx"], tmp_root): r["name"]
                    for r, fm in jobs}
            for i, fut in enumerate(as_completed(futs), 1):
                res = fut.result()
                tgt = next(r for r in todo if r["name"] == res["name"])
                if res["error"]:
                    tgt["fmap_error"] = res["error"]
                else:
                    tgt.update({k: v for k, v in res.items()
                                if k not in ("name", "error")})
                print(f"  [{i}/{len(futs)}] {res['name']}"
                      f"{'  ERROR' if res['error'] else ''}", flush=True)
        # fmap 之间的轨迹（同极性 → 纯头动），用来估 B0 漂移
        for key, fms in fm_cache.items():
            try:
                traj = phase_fmap_traj(fms, tmp_root)
            except Exception as e:                            # noqa: BLE001
                print(f"[qc-raw] fmap 轨迹失败 {key}: {e}")
                continue
            for r in todo:
                if (r["subject"], r["session"]) == key and r.get("fmap_name"):
                    t = traj.get(r["fmap_name"])
                    if t is not None:
                        r["fmap_pose"] = [float(x) for x in t]

    # ---- 汇总逐 run 指标 ----
    for r in runs:
        c = c_res.get(r["name"])
        if c is None:
            continue
        if c["error"]:
            r["qc_error"] = c["error"]
            continue
        mp = np.array(c["motion_afni"])
        if mp.ndim != 2 or mp.shape[0] == 0:
            r["qc_error"] = "3dvolreg 未产出运动参数"
            continue
        bids_mp = afni_to_bids_motion(mp)
        fd = fd_power(bids_mp)
        en = enorm_afni(mp)
        maxdisp = np.array(c["maxdisp"], dtype=float)
        maxdelt = np.array(c["maxdisp_delt"], dtype=float)
        dvars = np.array(c["dvars"], dtype=float)
        outfrac = r["_outfrac"]

        fd_v = fd[1:] if fd.size > 1 else fd
        r["_bids_mp"], r["_fd"], r["_enorm"] = bids_mp, fd, en
        r["_maxdisp"], r["_maxdisp_delt"], r["_dvars"] = maxdisp, maxdelt, dvars

        n_bad_fd = int(np.sum(fd_v > cfg["fd_thresh"]))
        n_bad_out = int(np.sum(outfrac > cfg["outlier_thresh"]))
        r.update({
            "base_run": r["_base_run"],
            "base_volume": r["_base_vol"],
            "min_outlier_index": r["_min_out_idx"],
            "fd_mean": float(np.nanmean(fd_v)) if fd_v.size else None,
            "fd_max": float(np.nanmax(fd_v)) if fd_v.size else None,
            "fd_p95": float(np.nanpercentile(fd_v, 95)) if fd_v.size else None,
            "n_fd_gt_thresh": n_bad_fd,
            "pct_fd_gt_thresh": 100.0 * n_bad_fd / fd_v.size if fd_v.size else None,
            "enorm_mean": float(np.nanmean(en[1:])) if en.size > 1 else None,
            "enorm_max": float(np.nanmax(en[1:])) if en.size > 1 else None,
            # maxdisp_max 是相对「基准卷」的绝对位移；global base 模式下基准卷可能
            # 在别的 run 里，所以它含 run 间偏移，不能当作 run 内头动。
            # maxdisp_range = 同一 run 内的极差，与基准卷无关，才是 run 内位移。
            "maxdisp_max": float(np.nanmax(maxdisp)) if maxdisp.size else None,
            "maxdisp_range": float(np.ptp(maxdisp)) if maxdisp.size else None,
            "maxdisp_delt_max": float(np.nanmax(maxdelt)) if maxdelt.size else None,
            "trans_range_max": float(np.ptp(bids_mp[:, :3], axis=0).max()),
            "rot_range_max_deg": float(np.degrees(np.ptp(bids_mp[:, 3:], axis=0).max())),
            "outlier_frac_max": float(outfrac.max()) if outfrac.size else None,
            "n_out_gt_thresh": n_bad_out,
            "pct_out_gt_thresh": 100.0 * n_bad_out / outfrac.size if outfrac.size else None,
            "dvars_mean": float(np.mean(dvars[1:])) if dvars.size > 1 else None,
            "dvars_p95": float(np.percentile(dvars[1:], 95)) if dvars.size > 1 else None,
            "tsnr_median": c.get("tsnr_median"),
            "n_nonsteady_est": estimate_nonsteady(r["_gs"]),
            "gs_drift_pct": estimate_drift_pct(r["_gs"]),
        })

    # ---- run 间头动 ----
    # global base 模式下所有 run 在同一坐标系，run 的中位姿态差即 run 间位移
    for key, grp in groups.items():
        # 组内相对覆盖率：以覆盖最多的 run 为 100%。
        # 掉得多 = 该 run 看到的脑比别人少，通常是挪出 FOV。
        vols = [r["brain_voxels"] for r in grp if r.get("brain_voxels")]
        if vols:
            mx = max(vols)
            for r in grp:
                if r.get("brain_voxels"):
                    r["coverage_rel_pct"] = 100.0 * r["brain_voxels"] / mx

        ok = [r for r in grp if r.get("_bids_mp") is not None and not r.get("qc_error")]
        if not ok:
            continue
        # 按**采集时间**排，不是按 run 编号 —— vs_prev 的语义是「采集顺序上的
        # 前一个 run」，重扫/补扫时编号顺序和采集顺序会不一致。
        # 没有 AcquisitionTime 时退回编号排序。
        ok.sort(key=lambda r: (scan_time_key(r.get("scan_date")),
                               r["task"], r["run"], r["name"]))
        if cfg["motion_base"] != "global":
            for r in ok:
                r["betrun_disp_vs_ref"] = None
                r["betrun_disp_vs_prev"] = None
            continue
        med = {r["name"]: np.median(r["_bids_mp"], axis=0) for r in ok}
        ref_name = ok[0]["_base_run"]
        ref_med = med.get(ref_name, med[ok[0]["name"]])

        def disp(a, b):
            d = np.abs(a - b)
            return float(d[:3].sum() + (d[3:] * FD_HEAD_RADIUS).sum())

        prev = None
        for r in ok:
            r["pose_median"] = [float(x) for x in med[r["name"]]]
            r["betrun_ref_run"] = ref_name
            r["betrun_disp_vs_ref"] = disp(med[r["name"]], ref_med)
            # 组内第一个 run 没有「前一个」，必须是 n/a 而不是 0.0 ——
            # 写 0 会被读成「和前一个 run 完全没动」，语义相反
            r["betrun_disp_vs_prev"] = (disp(med[r["name"]], med[prev])
                                        if prev else None)
            prev = r["name"]

    rc = finalize(runs, cfg, bids, out_dir, write_motion=True,
                  want_report=not args.no_report, elapsed=time.time() - t0)

    if not args.keep_tmp:
        shutil.rmtree(tmp_root, ignore_errors=True)
    return rc


def lookup_expected_volumes(cfg, task, run_id):
    """expected_volumes 支持两种写法：
         taskA: 470           → 该 task 所有 run 都是 470
         taskA: {"01": 464, "02": 470}  → 逐 run 指定
    找不到返回 None（走比例推断）。
    """
    cfgv = cfg.get("expected_volumes") or {}
    entry = cfgv.get(task, cfgv.get(f"task-{task}"))
    if entry is None:
        return None
    if isinstance(entry, dict):
        # run 号可能写成 "01" / "1" / 1，都认
        for k in (run_id, str(run_id).lstrip("0"), f"{str(run_id).lstrip('0'):0>2}"):
            if k in entry:
                return entry[k]
            if str(k) in {str(x) for x in entry}:
                return {str(x): v for x, v in entry.items()}[str(k)]
        return None
    return entry


def finalize(runs, cfg, bids, out_dir, write_motion, want_report, elapsed=None):
    """一致性比对 → 期望值 → 判定 → 落盘。指标已算好，本函数不碰 AFNI。"""
    out_dir = Path(out_dir)

    # ---- 一致性：与「多数派」比对 ----
    fields = cfg["consistency_fields"]
    by_task = defaultdict(list)
    for r in runs:
        by_task[r["task"]].append(r)
    for task, grp in by_task.items():
        modes = {}
        for f in fields:
            vals = [r.get(f) for r in grp if r.get(f) is not None]
            if vals:
                modes[f] = Counter(map(str, vals)).most_common(1)[0][0]
        for r in grp:
            bad = [f for f in fields
                   if r.get(f) is not None and f in modes and str(r[f]) != modes[f]]
            r["inconsistent_fields"] = bad

    # ---- 期望 volume 数：显式配置优先，否则用同 task 最长 run 做比例判定 ----
    for task, grp in by_task.items():
        # 用「同 task 最长 run」而非众数：变长范式下众数没有意义，
        # 而「最长的那个」是这个 task 实际能跑满的长度，截断判定用它才稳
        vols = [r["n_volumes"] for r in grp if r.get("n_volumes")]
        fallback = max(vols) if vols else None
        for r in grp:
            exp = lookup_expected_volumes(cfg, task, r.get("run"))
            if exp is not None:
                r["expected_volumes"] = exp
                r["expected_volumes_source"] = "config"
            else:
                r["expected_volumes"] = fallback
                r["expected_volumes_source"] = "inferred_max"

    # ---- 层面重新规划检测 ----
    # 技师中途重新定位（新 localizer + 重新规划层面）会让 FOV 中心和层面朝向
    # 跳变。这不是头动，是操作动作，但后果一样严重：不同规划下的 run 处在
    # 不同体素网格、磁化率畸变模式也不同。以组内**中位**规划为参照。
    presc = defaultdict(list)
    for r in runs:
        if r.get("fov_centre") and r.get("orient_cols"):
            presc[(r.get("subject"), r.get("session"), r.get("task"),
                   r.get("dim_str"), r.get("vox_str"))].append(r)
    for _, grp in presc.items():
        cen = np.median(np.array([r["fov_centre"] for r in grp], dtype=float), axis=0)
        mats = [np.array(r["orient_cols"], dtype=float).T for r in grp]
        ref = mats[int(np.argmin([sum(np.linalg.norm(m - n) for n in mats)
                                  for m in mats]))]
        for r, m in zip(grp, mats):
            r["presc_shift_mm"] = float(np.linalg.norm(
                np.array(r["fov_centre"], dtype=float) - cen))
            M = m @ np.linalg.inv(ref)
            r["presc_angle_deg"] = float(np.degrees(np.arccos(
                np.clip((np.trace(M) - 1) / 2, -1, 1))))

    # ---- 急动定位：峰值 TR + 尖峰个数 ----
    # 放在 finalize 而不是 main，因为 --rejudge 会把 motion tsv 读回来，
    # 这样调完 fd_max_warn 阈值不用重算 20 分钟就能重新定位尖峰。
    for r in runs:
        fd = r.get("_fd")
        if fd is None or not np.any(np.isfinite(fd)):
            continue
        fd = np.asarray(fd, dtype=float)
        r["fd_max_tr"] = int(np.nanargmax(fd))
        sp = np.where(fd > cfg["fd_max_warn"])[0]
        r["n_fd_spikes"] = int(sp.size)
        r["spike_trs"] = ";".join(str(int(i)) for i in sp[:20]) or None

    # ---- run 间位移的参照姿态：用组内**中位姿态**，不是配准基准卷 ----
    # 配准基准卷按 min-outlier 选（为了配准质量），它可能正好落在头位分布的边缘。
    # 实测踩过：被试在前两个 run 安顿下来，基准卷选中起始位置的 run-01，
    # 结果 4 个头位彼此一致的 run 全被判 fail，而真正离群的 run-01 判 pass。
    # 判定必须相对「多数 run 共同的头位」，即逐分量中位数。
    pose_groups = defaultdict(list)
    for r in runs:
        if r.get("pose_median"):
            pose_groups[(r.get("subject"), r.get("session"),
                         r.get("dim_str"), r.get("vox_str"),
                         r.get("pe_dir"))].append(r)
    for _, grp in pose_groups.items():
        poses = np.array([r["pose_median"] for r in grp], dtype=float)
        centre = np.median(poses, axis=0)
        for r in grp:
            d = np.abs(np.array(r["pose_median"], dtype=float) - centre)
            r["betrun_disp_vs_median"] = float(
                d[:3].sum() + (d[3:] * FD_HEAD_RADIUS).sum())

    # ---- 相对 tSNR：同协议队列内部比 ----
    # 文献里没有可引的 tSNR 绝对门限（见 references/metrics.md §4.1），
    # MRIQC 的建议做法就是「相对本队列找离群」。按 (task, 矩阵, 体素) 分组，
    # 组内 run 数 >= tsnr_rel_min_n 才有可比基线。
    by_proto = defaultdict(list)
    for r in runs:
        by_proto[(r.get("task"), r.get("dim_str"), r.get("vox_str"))].append(r)
    for _, grp in by_proto.items():
        vals = [r["tsnr_median"] for r in grp if r.get("tsnr_median")]
        if len(vals) < cfg["tsnr_rel_min_n"]:
            continue
        med = float(np.median(vals))
        if med <= 0:
            continue
        for r in grp:
            if r.get("tsnr_median"):
                r["tsnr_rel_pct"] = 100.0 * r["tsnr_median"] / med

    # ---- 判定 ----
    for r in runs:
        r["status"], reasons = judge_run(r, cfg)
        r["reason"] = "; ".join(reasons) if reasons else ""

    # ---- 写逐 run 运动时间序列（BIDS confounds 兼容列名）----
    if write_motion:
        for r in runs:
            if r.get("_bids_mp") is None:
                continue
            write_motion_tsv(out_dir / "motion", r, cfg)

    subj_rows = summarize_subjects(runs, bids, cfg)

    write_runs_tsv(out_dir / "qc_raw_runs.tsv", runs)
    write_subjects_tsv(out_dir / "qc_raw_subjects.tsv", subj_rows)
    write_json(out_dir / "qc_raw.json", runs, subj_rows, cfg, bids, out_dir)

    if want_report:
        try:
            sys.path.insert(0, str(Path(__file__).parent))
            from qc_report import build_report
            if not any("_bids_mp" in r for r in runs):
                from qc_report import attach_timeseries
                attach_timeseries(out_dir, runs)
            build_report(out_dir, runs, subj_rows, cfg)
            print(f"[qc-raw] 报告: {out_dir/'report.html'}")
        except Exception as e:                                # noqa: BLE001
            print(f"[qc-raw] 报告生成失败（指标已落盘）: {e}")

    try:
        sys.path.insert(0, str(Path(__file__).parent))
        from qc_summary import build_text
        (out_dir / "summary.txt").write_text(
            build_text(bids, out_dir, afni_names=True), encoding="utf-8")
    except Exception as e:                                    # noqa: BLE001
        print(f"[qc-raw] summary.txt 生成失败（其余已落盘）: {e}")

    n_fail = sum(1 for r in runs if r["status"] == "fail")
    n_warn = sum(1 for r in runs if r["status"] == "warn")
    tail = f"，用时 {elapsed:.0f}s" if elapsed is not None else ""
    print(f"\n[qc-raw] 完成 {len(runs)} run{tail}")
    print(f"[qc-raw] pass {len(runs)-n_fail-n_warn} / warn {n_warn} / fail {n_fail}")
    print(f"[qc-raw] 输出: {out_dir}")
    for r in sorted(runs, key=lambda x: x["name"]):
        if r["status"] != "pass":
            print(f"  [{r['status'].upper()}] {r['name']}: {r['reason']}")
    return 1 if n_fail else 0


def rejudge(out_dir: Path, cfg_path, want_report=True):
    """从已有的 qc_raw.json 重新判定，不重算任何指标。改阈值后用。"""
    out_dir = Path(out_dir).resolve()
    payload = json.loads((out_dir / "qc_raw.json").read_text())
    runs = payload["Runs"]
    cfg = dict(DEFAULT_CFG)
    cfg.update(payload.get("Thresholds") or {})
    if cfg_path:
        import yaml
        user = yaml.safe_load(Path(cfg_path).read_text()) or {}
        cfg.update(user.get("qc_raw", user))
    # 重读 header：只读文件头，不碰体数据，几乎瞬时。
    # 这样后来新增的 header 派生字段（如 pe_anat）在旧结果上也能补齐，
    # 不必为了加一列就重算 20 分钟。
    n_refresh = 0
    for r in runs:
        if not os.path.exists(r.get("path", "")):
            continue
        try:
            r.update(read_header(r))
            n_refresh += 1
        except Exception:                                     # noqa: BLE001
            pass
    # 旧结果里没有 pose_median（新增字段），从磁盘上的 motion tsv 补算 ——
    # 只是读几个小文本文件，不碰 AFNI
    # 无条件把 motion tsv 读回来：pose_median、急动定位都要用，
    # 只是读几个小文本文件，不碰 AFNI
    n_pose = 0
    try:
        from qc_report import attach_timeseries
        attach_timeseries(out_dir, runs)
        for r in runs:
            if r.get("_bids_mp") is not None and not r.get("pose_median"):
                r["pose_median"] = [float(x) for x in
                                    np.median(r["_bids_mp"], axis=0)]
                n_pose += 1
    except Exception as e:                                    # noqa: BLE001
        print(f"[qc-raw] 时间序列读回失败: {e}")
    # 回填 fmap 的方向信息：只读 sidecar，不重跑配准
    n_fm = 0
    need = [r for r in runs if r.get("fmap_name") and not r.get("fmap_pe_anat")]
    if need:
        cache = {}
        for r in need:
            key = (r.get("subject"), r.get("session"))
            if key not in cache:
                cache[key] = {f["name"]: f for f in
                              scan_fmaps(Path(payload["SourceBIDS"]), *key)}
            fm = cache[key].get(r["fmap_name"])
            if fm:
                r["fmap_pe_dir"] = fm.get("pe_dir")
                r["fmap_pe_anat"] = fm.get("pe_anat")
                r["pe_opposite_pair"] = opposite_pe(r.get("pe_anat"),
                                                    fm.get("pe_anat"))
                n_fm += 1
    print(f"[qc-raw] 重判 {len(runs)} run（不重算指标；{n_refresh} 个已刷新 header"
          f"{f'，{n_pose} 个补算 pose_median' if n_pose else ''}"
          f"{f'，{n_fm} 个回填 fmap 方向' if n_fm else ''}）")
    return finalize(runs, cfg, Path(payload["SourceBIDS"]), out_dir,
                    write_motion=False, want_report=want_report)


# ---------------------------------------------------------------- 输出

MOTION_COLS = ["trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z",
               "framewise_displacement", "enorm_afni", "maxdisp", "maxdisp_delt",
               "outlier_fraction", "dvars", "global_signal"]


def write_motion_tsv(mdir: Path, r, cfg):
    mdir.mkdir(parents=True, exist_ok=True)
    n = r["_bids_mp"].shape[0]

    def pad(a, first_na=False):
        a = np.asarray(a, dtype=float).ravel()
        if a.size == 0:
            return np.full(n, np.nan)
        if a.size < n:
            a = np.concatenate([np.full(n - a.size, np.nan), a]) if first_na else \
                np.concatenate([a, np.full(n - a.size, np.nan)])
        return a[:n]

    cols = np.column_stack([
        r["_bids_mp"],
        pad(r["_fd"]), pad(r["_enorm"]),
        pad(r["_maxdisp"]), pad(r.get("_maxdisp_delt", []), first_na=True),
        pad(r["_outfrac"]), pad(r["_dvars"], first_na=True), pad(r["_gs"]),
    ])
    base = f"{r['name'].replace('_bold','')}_desc-rawmotion_timeseries"
    with (mdir / f"{base}.tsv").open("w") as f:
        f.write("\t".join(MOTION_COLS) + "\n")
        for row in cols:
            f.write("\t".join("n/a" if not np.isfinite(v) else f"{v:.6f}"
                              for v in row) + "\n")
    meta = {
        "Sources": [r["path"]],
        "MotionEstimation": {
            "Software": "AFNI 3dvolreg",
            "BaseRun": r.get("base_run"),
            "BaseVolume": r.get("base_volume"),
            "BaseSelection": "min-outlier volume of BaseRun "
                             "(3dToutcount -automask -fraction, 跳过前 2 帧)",
            "BaseSpec": f"{r.get('base_run')}[{r.get('base_volume')}]",
            "ThisRunMinOutlierIndex": r.get("min_outlier_index"),
            "Registration": "rigid-body 6dof, EPI-to-EPI only (no anatomical/template alignment)",
        },
        "trans_x": {"Description": "位移，+x=Right (= -dL)", "Units": "mm"},
        "trans_y": {"Description": "位移，+y=Anterior (= -dP)", "Units": "mm"},
        "trans_z": {"Description": "位移，+z=Superior (= dS)", "Units": "mm"},
        "rot_x": {"Description": "绕 R-L 轴旋转 (AFNI pitch)", "Units": "rad"},
        "rot_y": {"Description": "绕 A-P 轴旋转 (AFNI yaw)", "Units": "rad"},
        "rot_z": {"Description": "绕 I-S 轴旋转 (AFNI roll)", "Units": "rad"},
        "framewise_displacement": {
            "Description": f"Power FD, 旋转按 r={FD_HEAD_RADIUS}mm 折算弧长",
            "Units": "mm"},
        "enorm_afni": {
            "Description": "AFNI motion_enorm：6 个原始参数(deg/mm混用)一阶差分的欧氏范数，"
                           "对应 afni_proc.py 的 0.2/0.3 censor 阈值",
            "Units": "arbitrary"},
        "maxdisp": {"Description": "3dvolreg -maxdisp1D：脑内体素相对基准卷的最大位移",
                    "Units": "mm"},
        "maxdisp_delt": {"Description": "相邻 TR 之间的最大体素位移变化", "Units": "mm"},
        "outlier_fraction": {"Description": "3dToutcount -automask -fraction",
                             "Units": "fraction"},
        "dvars": {"Description": "3dTto1D -method dvars，在未做 MoCo 的原始数据上计算",
                  "Units": "arbitrary"},
        "global_signal": {"Description": "automask 内全脑均值", "Units": "arbitrary"},
    }
    (mdir / f"{base}.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2))


RUN_COLS = [
    ("subject", "s"), ("session", "s"), ("task", "s"), ("acq", "s"), ("dir", "s"),
    ("run", "s"), ("status", "s"), ("reason", "s"),
    ("n_volumes", "i"), ("expected_volumes", "i"), ("expected_volumes_source", "s"),
    ("tr", "f3"), ("duration_sec", "f1"), ("n_slices", "i"),
    ("dim_str", "s"), ("vox_str", "s"), ("te", "f4"), ("flip_angle", "f1"),
    ("pe_dir", "s"), ("pe_anat", "s"), ("multiband", "s"),
    ("slice_timing_n", "i"),
    ("obliquity_deg", "f2"), ("presc_shift_mm", "f2"),
    ("presc_angle_deg", "f2"), ("scan_date", "s"), ("has_json", "s"), ("json_missing_fields", "l"),
    ("inconsistent_fields", "l"), ("n_events", "i"),
    ("fd_mean", "f4"), ("fd_max", "f4"), ("fd_max_tr", "i"),
    ("n_fd_spikes", "i"), ("spike_trs", "s"), ("fd_p95", "f4"),
    ("n_fd_gt_thresh", "i"), ("pct_fd_gt_thresh", "f2"),
    ("enorm_mean", "f4"), ("enorm_max", "f4"),
    ("maxdisp_range", "f3"), ("maxdisp_max", "f3"), ("maxdisp_delt_max", "f3"),
    ("trans_range_max", "f3"), ("rot_range_max_deg", "f3"),
    ("betrun_disp_vs_median", "f3"), ("betrun_ref_run", "s"),
    ("betrun_disp_vs_ref", "f3"), ("betrun_disp_vs_prev", "f3"),
    ("brain_voxels", "i"), ("brain_edge_frac", "f5"),
    ("coverage_rel_pct", "f2"),
    ("fmap_name", "s"), ("fmap_pe_dir", "s"), ("fmap_pe_anat", "s"),
    ("pe_opposite_pair", "s"), ("fmap_pair_method", "s"),
    ("fmap_disp_clean", "f3"), ("fmap_disp_pe", "f3"), ("fmap_disp_full", "f3"),
    ("outlier_frac_mean", "f4"), ("outlier_frac_max", "f4"),
    ("n_out_gt_thresh", "i"), ("pct_out_gt_thresh", "f2"),
    ("dvars_mean", "f2"), ("dvars_p95", "f2"),
    ("tsnr_median", "f2"), ("tsnr_rel_pct", "f1"),
    ("gs_drift_pct", "f2"), ("n_nonsteady_est", "i"),
    ("base_run", "s"), ("base_volume", "i"), ("min_outlier_index", "i"),
    ("path", "s"),
]
FMT = {"f1": "{:.1f}", "f2": "{:.2f}", "f3": "{:.3f}", "f4": "{:.4f}",
       "f5": "{:.5f}"}


def write_runs_tsv(path: Path, runs):
    with path.open("w") as f:
        f.write("\t".join(c for c, _ in RUN_COLS) + "\n")
        for r in sorted(runs, key=lambda x: (x["subject"], x["session"],
                                             x["task"], x["run"], x["name"])):
            f.write("\t".join(na(r.get(c), FMT.get(t, "{:.4f}"))
                              for c, t in RUN_COLS) + "\n")


SUBJ_COLS = ["subject", "session", "status", "reason", "n_runs", "n_runs_expected",
             "missing_runs", "protocol_mismatch", "n_fail", "n_warn",
             "n_anat", "n_fmap", "func_pe", "fmap_pe", "pe_opposite_ok",
             "fmap_disp_clean_max", "b0_drift_pe_mm", "head_drift_pe_mm",
             "fd_mean_median", "fd_mean_worst", "pct_bad_worst",
             "tsnr_median_min", "betrun_disp_max", "tasks"]


def summarize_subjects(runs, bids: Path, cfg):
    by_sub = defaultdict(list)
    for r in runs:
        by_sub[(r["subject"], r["session"])].append(r)

    # 每个 task 的期望 run 数：配置优先，否则用被试间众数
    exp_runs_cfg = cfg.get("expected_runs") or {}
    task_counts = defaultdict(list)
    for (sub, ses), grp in by_sub.items():
        c = Counter(r["task"] for r in grp)
        for t, n in c.items():
            task_counts[t].append(n)
    exp_runs = {}
    for t, ns in task_counts.items():
        v = exp_runs_cfg.get(t, exp_runs_cfg.get(f"task-{t}"))
        exp_runs[t] = v if v is not None else Counter(ns).most_common(1)[0][0]

    rows = []
    for (sub, ses), grp in sorted(by_sub.items()):
        sdir = bids / sub / ses if ses else bids / sub
        anat_files = sorted(sdir.glob("anat/*.nii*")) if sdir.exists() else []
        fmap_files = sorted(sdir.glob("fmap/*.nii*")) if sdir.exists() else []
        n_anat, n_fmap = len(anat_files), len(fmap_files)

        # func 的相位编码方向（本 session 内）
        func_pes = {r["pe_anat"] for r in grp if r.get("pe_anat")}
        func_pe = ",".join(sorted(func_pes)) if func_pes else None
        # fmap 的方向：直接读 sidecar + affine
        fmap_pes = set()
        for f in fmap_files:
            js, _ = sidecar_json(f)
            if not js:
                continue
            try:
                d = pe_direction(nib.load(str(f)).affine,
                                 js.get("PhaseEncodingDirection"))
            except Exception:                                 # noqa: BLE001
                d = None
            if d:
                fmap_pes.add(d)
        fmap_pe = ",".join(sorted(fmap_pes)) if fmap_pes else None

        # blip 矫正的硬性前提：至少有一个 fmap 与 func 方向严格相反。
        # 跨 session 用不同方向是允许的，所以只在 session 内部查。
        pe_opposite = None
        if cfg.get("check_pe_opposite") and func_pes and fmap_pes:
            pe_opposite = any(opposite_pe(a, b)
                              for a in func_pes for b in fmap_pes)
        c = Counter(r["task"] for r in grp)
        missing = []
        for t, exp in exp_runs.items():
            if t in c and c[t] < exp:
                missing.append(f"{t}:{c[t]}/{exp}")
        fails = [r for r in grp if r["status"] == "fail"]
        warns = [r for r in grp if r["status"] == "warn"]

        # 整个被试的所有 run 都和队列多数派不一致 → 是协议不同，不是 run 内漂移。
        # 这比单个 run 参数漂移严重得多：混进组分析会直接污染结果。
        inc_sets = [set(r.get("inconsistent_fields") or []) for r in grp]
        shared_inc = set.intersection(*inc_sets) if inc_sets else set()
        any_inc = set().union(*inc_sets) if inc_sets else set()
        protocol_mismatch = ",".join(sorted(shared_inc))

        reasons = []
        if shared_inc:
            reasons.append("整体协议与队列多数派不同: " + ",".join(sorted(shared_inc)))
        elif any_inc:
            reasons.append("部分 run 参数漂移: " + ",".join(sorted(any_inc)))
        if missing:
            reasons.append("run 不齐: " + ",".join(missing))
        if cfg["require_anat"] and n_anat == 0:
            reasons.append("缺 anat")
        if cfg["require_fmap"] and n_fmap == 0:
            reasons.append("缺 fmap")
        exp_fmap = cfg.get("expected_fmap")
        if exp_fmap and n_fmap < exp_fmap:
            reasons.append(f"fmap 只有 {n_fmap} 个，协议要 {exp_fmap} 个")
        if len(func_pes) > 1:
            reasons.append(f"session 内 func 的相位编码方向不统一: {func_pe}")
        if pe_opposite is False:
            reasons.append(f"没有与 func ({func_pe}) 方向相反的 fmap"
                           f"（现有 fmap: {fmap_pe}），blip 矫正做不了")
        if fails:
            reasons.append(f"{len(fails)} 个 run 不合格: " +
                           ",".join(r["name"] for r in fails))
        fmap_short = bool(cfg.get("expected_fmap")) and n_fmap < cfg["expected_fmap"]
        status = "fail" if (missing or fails or shared_inc or fmap_short or
                            pe_opposite is False or len(func_pes) > 1 or
                            (cfg["require_anat"] and n_anat == 0)) else \
                 ("warn" if warns or reasons else "pass")
        if warns and not fails:
            reasons.append(f"{len(warns)} 个 run 有告警")

        # --- B0 漂移估计 ---
        # func 轨迹 = 真实头动 + 该极性的畸变漂移
        # fmap 轨迹 = 真实头动 - 该极性的畸变漂移（PE 相反）
        # 两者之差的一半 ≈ 畸变漂移；之和的一半 ≈ 真实头动。
        # ⚠ 这依赖「畸变严格反对称」的假设，而且刚体拟合只给全局平均位移，
        #   不是畸变场本身 —— 是量级估计，不是测量。
        b0_drift = head_drift_pe = None
        paired = [r for r in grp if r.get("pose_median") and r.get("fmap_pose")]
        if len(paired) >= 2:
            pe = next((r.get("pe_anat") for r in paired if r.get("pe_anat")), None)
            ax = ANAT_AXIS.get(pe.split("->")[0]) if pe else 1
            paired.sort(key=lambda r: scan_time_key(r.get("scan_date")))
            d_func = paired[-1]["pose_median"][ax] - paired[0]["pose_median"][ax]
            d_fmap = paired[-1]["fmap_pose"][ax] - paired[0]["fmap_pose"][ax]
            b0_drift = abs(d_func - d_fmap) / 2.0
            head_drift_pe = abs(d_func + d_fmap) / 2.0
        if b0_drift is not None and b0_drift > cfg["b0_drift_warn"]:
            reasons.append(f"整场 B0 漂移约 {b0_drift:.1f}mm（PE 轴），"
                           f"远端 run 用同一个 fmap 做 blip 可能对不齐")

        vals = lambda k: [r[k] for r in grp if r.get(k) is not None]  # noqa: E731
        fdm, pb = vals("fd_mean"), vals("pct_fd_gt_thresh")
        ts, br = vals("tsnr_median"), vals("betrun_disp_vs_ref")
        rows.append({
            "subject": sub, "session": ses, "status": status,
            "reason": "; ".join(reasons),
            "n_runs": len(grp),
            "n_runs_expected": sum(exp_runs.get(t, 0) for t in c),
            "missing_runs": ",".join(missing),
            "protocol_mismatch": protocol_mismatch,
            "func_pe": func_pe, "fmap_pe": fmap_pe,
            "b0_drift_pe_mm": b0_drift, "head_drift_pe_mm": head_drift_pe,
            "fmap_disp_clean_max": max(
                [r["fmap_disp_clean"] for r in grp
                 if r.get("fmap_disp_clean") is not None] or [0]) or None,
            "pe_opposite_ok": pe_opposite,
            "n_fail": len(fails), "n_warn": len(warns),
            "n_anat": n_anat, "n_fmap": n_fmap,
            "fd_mean_median": float(np.median(fdm)) if fdm else None,
            "fd_mean_worst": max(fdm) if fdm else None,
            "pct_bad_worst": max(pb) if pb else None,
            "tsnr_median_min": min(ts) if ts else None,
            "betrun_disp_max": max(br) if br else None,
            "tasks": ",".join(f"{t}x{n}" for t, n in sorted(c.items())),
        })
    return rows


def write_subjects_tsv(path: Path, rows):
    with path.open("w") as f:
        f.write("\t".join(SUBJ_COLS) + "\n")
        for r in rows:
            f.write("\t".join(na(r.get(c), "{:.3f}") for c in SUBJ_COLS) + "\n")


def write_json(path: Path, runs, subj_rows, cfg, bids, out_dir):
    clean = []
    for r in runs:
        clean.append({k: (v if not isinstance(v, np.generic) else v.item())
                      for k, v in r.items() if not k.startswith("_")})
    payload = {
        "GeneratedAt": datetime.now(timezone.utc).isoformat(),
        "Tool": "qc-raw",
        "SourceBIDS": str(bids),
        "OutputDir": str(out_dir),
        "Thresholds": cfg,
        "Subjects": subj_rows,
        "Runs": clean,
    }
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=str))


if __name__ == "__main__":
    sys.exit(main())
