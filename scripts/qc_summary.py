#!/usr/bin/env python3
"""
qc_summary.py — 把 qc-raw 的结果打成终端文本表

两张表：
  === 扫描序列信息汇总 ===   每个序列的 volume 数 / 体素 / FOV / 相位编码方向
  === 头动情况汇总 ===       每个 run 的 Mean_FD / Max_FD 与 6 参数变化范围

用法：
    python3 qc_summary.py --bids <bids>                     # 只出序列表
    python3 qc_summary.py --bids <bids> --qc <qc-raw 输出>   # 两张表都出
    python3 qc_summary.py --bids <bids> --qc <dir> --afni-names   # 6 参数用
                                                            # AFNI 原生列名
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from pathlib import Path

import numpy as np

try:
    import nibabel as nib
except ImportError:
    raise SystemExit("需要 nibabel：pip install nibabel")


sys.path.insert(0, str(Path(__file__).parent))
from qc_raw import pe_direction          # noqa: E402  单一实现，避免两份逻辑漂移


def series_rows(bids: Path):
    rows = []
    for f in sorted(glob.glob(str(bids / "sub-*/**/*.nii*"), recursive=True)):
        if "/derivatives/" in f:
            continue
        img = nib.load(f)
        shape = img.shape
        zooms = [float(z) for z in img.header.get_zooms()]
        nvol = int(shape[3]) if len(shape) > 3 else 1
        js_path = f.split(".nii")[0] + ".json"
        js = {}
        if os.path.exists(js_path):
            try:
                js = json.loads(Path(js_path).read_text())
            except json.JSONDecodeError:
                pass
        name = js.get("SeriesDescription") or Path(f).name.split(".nii")[0]
        vox = " x ".join(f"{z:.2f}" for z in zooms[:3])
        fov = " x ".join(f"{shape[i]*zooms[i]:.0f}" for i in range(3))
        rows.append({
            "sequence": name,
            "trs": nvol,
            "resolution": f"{vox} mm",
            "fov": f"{fov} mm",
            "matrix": " x ".join(str(s) for s in shape[:3]),
            "pe": pe_direction(img.affine,
                               js.get("PhaseEncodingDirection")) or "N/A",
            "pe_bids": js.get("PhaseEncodingDirection") or "N/A",
            "path": f,
        })
    return rows


AFNI_COLS = ["dS(mm)", "dL(mm)", "dP(mm)", "Roll(°)", "Pitch(°)", "Yaw(°)"]
BIDS_COLS = ["trans_x", "trans_y", "trans_z", "rot_x(°)", "rot_y(°)", "rot_z(°)"]
# BIDS 列 → AFNI 列的位置映射。取的是极差，符号无关，所以只换名字和顺序：
#   dS ↔ trans_z, dL ↔ trans_x, dP ↔ trans_y
#   Roll ↔ rot_z,  Pitch ↔ rot_x, Yaw ↔ rot_y
AFNI_ORDER = [2, 0, 1, 5, 3, 4]


def motion_rows(qc_dir: Path, afni_names=False):
    rows = []
    for f in sorted(glob.glob(str(qc_dir / "motion" / "*_timeseries.tsv"))):
        d = np.genfromtxt(f, delimiter="\t", names=True,
                          missing_values="n/a", filling_values=np.nan)
        mp = np.column_stack([d[c] for c in d.dtype.names[:6]])
        rng = np.ptp(mp, axis=0)
        rng[3:] = np.degrees(rng[3:])          # rad → deg
        fd = d["framewise_displacement"]
        name = Path(f).name.replace("_desc-rawmotion_timeseries.tsv", "")
        vals = rng[AFNI_ORDER] if afni_names else rng
        rows.append({
            "run": name,
            "mean_fd": float(np.nanmean(fd[1:])),
            "max_fd": float(np.nanmax(fd[1:])),
            "vals": vals,
            "nvol": mp.shape[0],
        })
    return rows


FMAP_COLS = ["Run", "func_PE", "func_Phase_Dir", "配对fmap", "fmap_PE",
             "fmap_Phase_Dir", "方向相反?", "配对依据", "头位差(净)", "PE轴(含畸变)"]


def fmap_rows(qc_dir: Path):
    """从 qc_raw.json 读 func↔fmap 配对信息。"""
    f = Path(qc_dir) / "qc_raw.json"
    if not f.exists():
        return []
    runs = json.loads(f.read_text()).get("Runs", [])
    out = []
    for r in sorted(runs, key=lambda x: str(x.get("name"))):
        if not r.get("fmap_name"):
            continue
        fmt = lambda v: "n/a" if v is None else (f"{v:.3f}" if isinstance(v, float)
                                                 else str(v))   # noqa: E731
        out.append([r["name"], fmt(r.get("pe_dir")), fmt(r.get("pe_anat")),
                    fmt(r.get("fmap_name")), fmt(r.get("fmap_pe_dir")),
                    fmt(r.get("fmap_pe_anat")),
                    {True: "是", False: "否"}.get(r.get("pe_opposite_pair"), "n/a"),
                    fmt(r.get("fmap_pair_method")),
                    fmt(r.get("fmap_disp_clean")), fmt(r.get("fmap_disp_pe"))])
    return out


def render(headers, rows, aligns=None):
    """定宽文本表。"""
    ncol = len(headers)
    aligns = aligns or [">"] * ncol
    w = [len(h) for h in headers]
    for r in rows:
        for i, c in enumerate(r):
            w[i] = max(w[i], len(str(c)))
    out = ["  ".join(f"{h:{aligns[i]}{w[i]}s}" for i, h in enumerate(headers))]
    for r in rows:
        out.append("  ".join(f"{str(c):{aligns[i]}{w[i]}s}"
                             for i, c in enumerate(r)))
    return "\n".join(out)


def build_text(bids: Path, qc_dir=None, afni_names=True):
    """生成两张表的完整文本。qc_raw.py 跑完会调它写 summary.txt。"""
    L = ["=== 扫描序列信息汇总 ===", ""]
    srows = series_rows(Path(bids))
    L.append(render(["Sequence", "TRs", "Resolution", "FOV", "Matrix",
                     "PE(BIDS)", "Phase_Dir"],
                    [[r["sequence"], str(r["trs"]), r["resolution"], r["fov"],
                      r["matrix"], r["pe_bids"], r["pe"]] for r in srows],
                    aligns=[">"] * 7))
    L += ["",
          "Phase_Dir 由 affine + PhaseEncodingDirection 算出"
          "（BIDS 规范：字母是图像轴，无 '-' 则沿该轴递增方向）。",
          "⚠ 序列名里的 AP/PA 经常和实际方向相反，以 Phase_Dir 为准，不要信名字。"]
    if qc_dir is None:
        return "\n".join(L)

    mrows = motion_rows(Path(qc_dir), afni_names)
    cols = AFNI_COLS if afni_names else BIDS_COLS
    L += ["", "", "=== 头动情况汇总 ===",
          "Mean_FD / Max_FD = Power FD（Σ|Δ平移| + Σ|Δ旋转(rad)|×50mm），单位 mm",
          "                   与 fmriprep/nipype FramewiseDisplacement 逐帧一致（差 2e-16）",
          f"{'/'.join(c.split('(')[0] for c in cols)} = 该 run 内的变化范围"
          "（peak-to-peak），平移 mm / 旋转 deg", ""]
    L.append(render(["Run", "TRs", "Mean_FD", "Max_FD"] + cols,
                    [[r["run"], str(r["nvol"]), f"{r['mean_fd']:.3f}",
                      f"{r['max_fd']:.3f}"] + [f"{v:.3f}" for v in r["vals"]]
                     for r in mrows]))

    frows = fmap_rows(qc_dir)
    if frows:
        L += ["", "", "=== func ↔ reverse-PE 配对检查 ===",
              "方向必须相反，否则 blip / topup 矫正做不了。",
              "头位差(净) 只含不受畸变污染的分量（垂直 PE 的平移 + 绕 PE 轴的旋转）；",
              "PE 轴那列混了头动与畸变差异，**不参与判定**，仅供参考。", ""]
        L.append(render(FMAP_COLS, frows))
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description="qc-raw 结果的文本表呈现")
    ap.add_argument("--bids", required=True)
    ap.add_argument("--qc", help="qc-raw 输出目录（给了才出头动表）")
    ap.add_argument("--afni-names", action="store_true",
                    help="头动 6 参数用 AFNI 原生列名 dS/dL/dP/Roll/Pitch/Yaw")
    ap.add_argument("--out", help="写到文件而不是只打印")
    args = ap.parse_args()

    txt = build_text(Path(args.bids).resolve(),
                     Path(args.qc).resolve() if args.qc else None,
                     args.afni_names)
    print(txt)
    if args.out:
        Path(args.out).write_text(txt, encoding="utf-8")
        print(f"\n已写入 {args.out}")
    return 0


def _unused_old_main(args):
    bids = Path(args.bids).resolve()

    print("=== 扫描序列信息汇总 ===\n")
    srows = series_rows(bids)
    print(render(["Sequence", "TRs", "Resolution", "FOV", "Matrix",
                  "PE(BIDS)", "Phase_Dir"],
                 [[r["sequence"], str(r["trs"]), r["resolution"], r["fov"],
                   r["matrix"], r["pe_bids"], r["pe"]] for r in srows],
                 aligns=[">"] * 7))
    print("\nPhase_Dir 由 affine + PhaseEncodingDirection 算出（BIDS 规范：字母是图像轴，"
          "无 '-' 则沿该轴递增方向）。\n"
          "⚠ 序列名里的 AP/PA 经常和实际方向相反，以 Phase_Dir 为准，不要信名字。")

    if not args.qc:
        print("\n（未给 --qc，跳过头动表）")
        return 0

    qc = Path(args.qc).resolve()
    mrows = motion_rows(qc, args.afni_names)
    cols = AFNI_COLS if args.afni_names else BIDS_COLS
    print("\n\n=== 头动情况汇总 ===")
    print("Mean_FD / Max_FD = Power FD（Σ|Δ平移| + Σ|Δ旋转(rad)|×50mm），单位 mm")
    print(f"{'/'.join(c.split('(')[0] for c in cols)} = 该 run 内的变化范围"
          f"（peak-to-peak），平移 mm / 旋转 deg\n")
    print(render(["Run", "TRs", "Mean_FD", "Max_FD"] + cols,
                 [[r["run"], str(r["nvol"]), f"{r['mean_fd']:.3f}",
                   f"{r['max_fd']:.3f}"] + [f"{v:.3f}" for v in r["vals"]]
                  for r in mrows]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
