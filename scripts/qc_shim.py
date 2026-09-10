#!/usr/bin/env python3
"""
qc_shim.py — 匀场框 / 扫描框几何与匀场结果一致性检查（**吃 DICOM，不吃 BIDS**）

匀场框（Siemens: adjustment volume）只存在于 DICOM 私有头里，dcm2niix 不会把它
写进 BIDS sidecar。所以这一项必须在**转 BIDS 之前**、拿到原始 DICOM 时查，
查完了这批 DICOM 可以照常走 bids-convert。

三层检查：
  [1] series 内   匀场框 vs 扫描框是否重合（中心偏移 / 法向夹角 / FOV 覆盖）
  [2] series 间   匀场框几何分组——同一次定位的 run 应当落在同一组
  [3] series 间   实际匀场结果分组（shim 电流 + 中心频率）——框一样 ≠ 匀场结果一样

用法：
    python3 qc_shim.py --dicom <DICOM目录>
    python3 qc_shim.py --dicom <dir> -o <输出目录>      # 落 tsv + json
    python3 qc_shim.py --dicom <dir> --verify           # 打印坐标约定核对
    python3 qc_shim.py --dicom <dir> --max-offset 3     # 收紧 series 内对齐阈值

退出码：有 fail 返回 1（可接 CI）。

数据来源（Siemens 经典 DICOM / .IMA）：私有标签 (0029,1020) CSA Series Header
里的 ASCCONV 文本段。
    sAdjData.sAdjVolume.sPosition.dSag/dCor/dTra   匀场框中心
    sAdjData.sAdjVolume.sNormal.dSag/dCor/dTra     匀场框法向
    sAdjData.sAdjVolume.dThickness/dReadoutFOV/dPhaseFOV/dInPlaneRot
    sAdjData.lCoupleAdjVolTo                       匀场框是否耦合到 slice group
    sSliceArray.asSlice[i].*                       扫描框（结构同上）
    sGRADSPEC.alShimCurrent[0..14]                 实际 shim 电流
    sTXSPEC.asNucleusInfo[0].lFrequency            实际中心频率
坐标约定已实测核对：ASCCONV 的 (dSag,dCor,dTra) 即 DICOM LPS 的 (x,y,z)，
法向与 ImageOrientationPatient 叉乘逐位吻合，无需转换（--verify 可复现）。
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import re
import sys
from pathlib import Path

import numpy as np

try:
    import pydicom
except ImportError:
    raise SystemExit("需要 pydicom：pip install pydicom")


CSA_SERIES_TAG = (0x0029, 0x1020)
ACQ_MATRIX_TAG = (0x0018, 0x1310)
BWPPPE_TAG = (0x0019, 0x1028)          # BandwidthPerPixelPhaseEncode


# ---------------------------------------------------------------- 解析
def read_ascconv(path):
    """读一个 DICOM，返回 (dataset, ASCCONV 键值 dict)。非西门子 / 无 CSA 返回 {}。"""
    ds = pydicom.dcmread(path, stop_before_pixels=True)
    if CSA_SERIES_TAG not in ds:
        return ds, {}
    txt = ds[CSA_SERIES_TAG].value.decode("latin-1", errors="replace")
    m = re.search(r"### ASCCONV BEGIN.*?###(.*?)### ASCCONV END ###", txt, re.S)
    if not m:
        return ds, {}
    out = {}
    for line in m.group(1).splitlines():
        if "=" not in line:
            continue
        k, _, v = line.partition("=")
        k, v = k.strip(), v.strip()
        if not k or "__attribute__" in k:
            continue
        try:
            out[k] = int(v, 0) if re.fullmatch(r"[+-]?(0x[0-9a-fA-F]+|\d+)", v) else float(v)
        except ValueError:
            out[k] = v.strip('"')
    return ds, out


def _v(d, prefix, default=(0.0, 0.0, 0.0)):
    """ASCCONV 省略等于默认值的键，取不到就回退（几乎总是 0）。"""
    return np.array([d.get(f"{prefix}.dSag", default[0]),
                     d.get(f"{prefix}.dCor", default[1]),
                     d.get(f"{prefix}.dTra", default[2])], float)


def _unit(v):
    n = float(np.linalg.norm(v))
    return v / n if n else v


def shim_box(d):
    return dict(center=_v(d, "sAdjData.sAdjVolume.sPosition"),
                normal=_unit(_v(d, "sAdjData.sAdjVolume.sNormal", (0, 0, 1))),
                thick=float(d.get("sAdjData.sAdjVolume.dThickness", 0.0)),
                ro=float(d.get("sAdjData.sAdjVolume.dReadoutFOV", 0.0)),
                pe=float(d.get("sAdjData.sAdjVolume.dPhaseFOV", 0.0)),
                rot=math.degrees(d.get("sAdjData.sAdjVolume.dInPlaneRot", 0.0)))


def scan_box(d):
    """扫描框。2D 多层：slab 中心取首末层中点，slab 厚取覆盖范围。"""
    n = int(d.get("sSliceArray.lSize", 1) or 1)
    p0 = _v(d, "sSliceArray.asSlice[0].sPosition")
    pN = _v(d, f"sSliceArray.asSlice[{n-1}].sPosition") if n > 1 else p0
    th = float(d.get("sSliceArray.asSlice[0].dThickness", 0.0))
    return dict(center=(p0 + pN) / 2.0,
                normal=_unit(_v(d, "sSliceArray.asSlice[0].sNormal", (0, 0, 1))),
                thick=float(np.linalg.norm(pN - p0)) + th if n > 1 else th,
                ro=float(d.get("sSliceArray.asSlice[0].dReadoutFOV", 0.0)),
                pe=float(d.get("sSliceArray.asSlice[0].dPhaseFOV", 0.0)),
                rot=math.degrees(d.get("sSliceArray.asSlice[0].dInPlaneRot", 0.0)),
                n_slices=n, slice_thick=th)


def compare(sh, sc):
    """匀场框相对扫描框的偏差，拆成沿法向与面内两部分。"""
    delta = sh["center"] - sc["center"]
    cosang = abs(float(np.dot(sh["normal"], sc["normal"])))
    ang = math.degrees(math.acos(max(-1.0, min(1.0, cosang))))
    along = float(np.dot(delta, sc["normal"]))
    return dict(dist=float(np.linalg.norm(delta)), along=along,
                inplane=float(np.linalg.norm(delta - along * sc["normal"])),
                angle=ang,
                d_ro=sh["ro"] - sc["ro"], d_pe=sh["pe"] - sc["pe"],
                d_th=sh["thick"] - sc["thick"])


def rot_delta(a, b):
    """面内朝向差。对 180° 取模——PE 极性反转会把 rot 加 180°，那不是真差异。"""
    d = abs(a - b) % 180.0
    return min(d, 180.0 - d)


# ---------------------------------------------------------------- 扫描目录
def index_series(root: Path, exts=("*.IMA", "*.dcm", "*.DCM")):
    """按 SeriesNumber 归组，每组只留第一个文件（协议在 series 级别，逐帧相同）。"""
    files = []
    for e in exts:
        files += glob.glob(str(root / "**" / e), recursive=True)
        files += glob.glob(str(root / e))
    files = sorted(set(files))
    if not files:
        raise SystemExit(f"目录里没找到 DICOM（{'/'.join(exts)}）：{root}")

    first, count = {}, {}
    for f in files:
        try:
            hdr = pydicom.dcmread(f, stop_before_pixels=True,
                                  specific_tags=["SeriesNumber"])
            sn = int(hdr.SeriesNumber)
        except Exception:
            continue
        count[sn] = count.get(sn, 0) + 1
        if sn not in first:
            first[sn] = f
    return first, count


def slice_center_from_dicom(ds):
    """用标准 DICOM 字段算首层中心。mosaic 下 (Cols-1)/2 恰等于
    (mosaic-single)/2 + (single-1)/2，故对 mosaic 与非 mosaic 同样成立。"""
    iop = np.array(ds.ImageOrientationPatient, float)
    ipp = np.array(ds.ImagePositionPatient, float)
    ps = [float(x) for x in ds.PixelSpacing]
    return ipp + iop[:3] * ps[1] * (ds.Columns - 1) / 2 + iop[3:] * ps[0] * (ds.Rows - 1) / 2


# ---------------------------------------------------------------- 主流程
def collect(root: Path):
    first, count = index_series(root)
    rows = []
    for sn in sorted(first):
        ds, d = read_ascconv(first[sn])
        desc = str(getattr(ds, "SeriesDescription", "?"))
        rec = dict(series=sn, desc=desc, nfiles=count.get(sn, 0), file=first[sn],
                   has_shim="sAdjData.sAdjVolume.dThickness" in d, ds=ds)
        if not d:
            rec["note"] = "无 CSA/ASCCONV（非西门子或已去私有头）"
            rows.append(rec)
            continue
        sc = scan_box(d)
        rec.update(n_slices=sc["n_slices"], scan=sc,
                   coupled=int(d.get("sAdjData.lCoupleAdjVolTo", -1)),
                   shim_mode=int(d.get("sAdjData.uiAdjShimMode", 0)),
                   adj_prot_id=int(d.get("sAdjData.uiAdjProtID", -1)),
                   f0=int(d.get("sTXSPEC.asNucleusInfo[0].lFrequency", 0)),
                   shim=[int(d.get(f"sGRADSPEC.alShimCurrent[{i}]", 0)) for i in range(9)],
                   bwpppe=float(ds[BWPPPE_TAG].value) if BWPPPE_TAG in ds else None,
                   ascconv=d, ds=ds)
        if rec["has_shim"]:
            sh = shim_box(d)
            rec.update(shimbox=sh, cmp=compare(sh, sc))
        else:
            rec["note"] = "该 series 无匀场框记录（scout / 派生图常见）"
        rows.append(rec)
    return rows


def judge(rows, max_offset, max_angle, max_rot):
    """写入 status/reason。返回是否有 fail。"""
    have = [r for r in rows if r.get("has_shim")]

    # [1] series 内：匀场框 vs 扫描框
    for r in rows:
        r["status"], r["reason"] = "pass", []
        if not r.get("has_shim"):
            r["status"] = "skip"
            continue
        c = r["cmp"]
        if c["dist"] > max_offset:
            r["status"] = "fail"
            r["reason"].append(f"匀场框中心偏离扫描框 {c['dist']:.1f}mm "
                               f"(沿法向 {c['along']:+.1f} / 面内 {c['inplane']:.1f})")
        if c["angle"] > max_angle:
            r["status"] = "fail"
            r["reason"].append(f"匀场框法向与扫描框夹角 {c['angle']:.2f}°")
        if min(r["shimbox"]["ro"], r["shimbox"]["pe"]) > 0 and (
                c["d_ro"] < -1 or c["d_pe"] < -1 or c["d_th"] < -1):
            if r["status"] == "pass":
                r["status"] = "warn"
            r["reason"].append(f"匀场框小于扫描框 (ΔRO={c['d_ro']:+.0f} "
                               f"ΔPE={c['d_pe']:+.0f} ΔTH={c['d_th']:+.0f} mm)")

    # [2] 几何分组：同法向 + 同中心 = 同一次定位；组内比面内朝向
    geo = {}
    for r in have:
        sh = r["shimbox"]
        key = (tuple(sh["center"].round(1)), tuple(np.abs(sh["normal"]).round(3)),
               round(sh["thick"]), round(sh["ro"]), round(sh["pe"]))
        geo.setdefault(key, []).append(r)
    for members in geo.values():
        base = members[0]["shimbox"]["rot"]
        for r in members:
            dr = rot_delta(r["shimbox"]["rot"], base)
            r["rot_vs_group"] = dr
            if dr > max_rot:
                if r["status"] == "pass":
                    r["status"] = "warn"
                r["reason"].append(
                    f"面内朝向与同定位组差 {dr:.2f}°（组基准 {base:.3f}°，本 series "
                    f"{r['shimbox']['rot']:.3f}°）——多半是协议模板没继承定位")

    # [2b] 扫描框分组：同规格的 EPI（func 与其反向 PE 配对）应当落在同一组。
    # 独立于匀场框分组——匀场框冻结时所有 EPI 会挤进同一个匀场几何组，
    # 只比匀场框看不出 func/reverse 之间的定位差异。
    spec = {}
    for r in have:
        sc = r["scan"]
        if sc["n_slices"] <= 1:            # 只管 2D 多层（EPI），3D 解剖不参与
            continue
        spec.setdefault((round(sc["ro"]), round(sc["pe"]), round(sc["thick"]),
                         sc["n_slices"]), []).append(r)
    scangrp = {}
    for members in spec.values():
        sub = {}
        for r in members:
            sc = r["scan"]
            key = (tuple(sc["center"].round(2)), tuple(np.abs(sc["normal"]).round(4)),
                   round(sc["rot"] % 180.0, 2))
            sub.setdefault(key, []).append(r)
        for i, (key, ms) in enumerate(sub.items(), 1):
            for r in ms:
                r["scan_group"] = i
        if len(sub) > 1:
            for r in members:
                if r["status"] == "pass":
                    r["status"] = "warn"
                r["reason"].append(
                    f"同规格 EPI 的扫描框分成 {len(sub)} 组（本 series 属第 "
                    f"{r['scan_group']} 组）——func 与反向 PE 图定位可能不配对")
        scangrp.update({(k, id(members)): v for k, v in sub.items()})

    # [3] 匀场结果分组：框一样不代表 adjust 结果一样
    shimgrp = {}
    for r in have:
        shimgrp.setdefault((tuple(r["shim"]), r["f0"]), []).append(r)
    for key, members in shimgrp.items():
        for r in members:
            r["shim_group"] = list(shimgrp).index(key) + 1

    # 同一几何组里出现多套 shim → 反向 PE 配对做畸变校正的前提被破坏
    for members in geo.values():
        groups = {r["shim_group"] for r in members}
        if len(groups) > 1:
            for r in members:
                if r["status"] == "pass":
                    r["status"] = "warn"
                r["reason"].append(
                    f"同一定位下存在 {len(groups)} 套匀场结果（本 series 属第 "
                    f"{r['shim_group']} 套）——反向 PE 配对不在同一 B0 场下")

    return geo, shimgrp, spec, any(r["status"] == "fail" for r in rows)


# ---------------------------------------------------------------- 输出
def print_table(rows):
    print(f"{'Ser':>4} {'序列名':<40} {'层':>4} {'偏移mm':>7} {'法向°':>6} "
          f"{'朝向°':>6} {'shim组':>6} {'耦合':>4} {'状态':>5}  原因")
    print("-" * 132)
    for r in rows:
        if not r.get("has_shim"):
            print(f"{r['series']:>4} {r['desc'][:40]:<40} {'—':>4} {'—':>7} {'—':>6} "
                  f"{'—':>6} {'—':>6} {'—':>4} {'skip':>5}  {r.get('note','')}")
            continue
        c = r["cmp"]
        print(f"{r['series']:>4} {r['desc'][:40]:<40} {r['n_slices']:>4} "
              f"{c['dist']:>7.1f} {c['angle']:>6.2f} {r.get('rot_vs_group',0):>6.2f} "
              f"{r['shim_group']:>6} {r['coupled']:>4} {r['status']:>5}  "
              f"{'; '.join(r['reason'])}")


def print_groups(geo, shimgrp):
    """匀场结果的差异只在**同一几何组内**比较——跨定位、跨序列类型比没有意义
    （拿 MP2RAGE 的 shim 去减 EPI 的，差值和折算出的体素位移都是假的）。"""
    print("\n=== 匀场框几何分组（同一次定位应落在同一组）===")
    for i, (key, members) in enumerate(geo.items(), 1):
        c = np.array(key[0])
        print(f"  几何组#{i}  center=({c[0]:7.2f},{c[1]:7.2f},{c[2]:7.2f}) "
              f"FOV={key[3]:.0f}×{key[4]:.0f}×{key[2]:.0f}mm")
        print(f"            ← series {','.join(str(m['series']) for m in members)}")

        rots = sorted({round(m["shimbox"]["rot"], 3) for m in members})
        if len(rots) > 1:
            print(f"            ⚠ 组内面内朝向不统一: {rots}°")

        # 组内的匀场结果细分
        sub = {}
        for m in members:
            sub.setdefault(m["shim_group"], []).append(m)
        if len(sub) == 1:
            print(f"            匀场结果统一（shim组#{list(sub)[0]}）")
            continue
        print(f"            ⚠ 组内有 {len(sub)} 套匀场结果：")
        base_g = None
        for gid, ms in sorted(sub.items()):
            shim, f0 = tuple(ms[0]["shim"]), ms[0]["f0"]
            tag = f"shim组#{gid}  f0={f0}"
            sers = ",".join(str(m["series"]) for m in ms)
            if base_g is None:
                base_g, base_f0 = shim, f0
                print(f"              {tag}  (组内基准)  ← series {sers}")
            else:
                dv = [a - b for a, b in zip(shim, base_g)]
                df = f0 - base_f0
                print(f"              {tag}  ← series {sers}")
                print(f"                Δvs组内基准: shim={dv}  Δf0={df:+d}Hz")
                bw = next((m["bwpppe"] for m in ms if m.get("bwpppe")), None)
                if bw and df:
                    print(f"                零阶项参考量级: {abs(df)/bw:.1f} voxel "
                          f"(BWPPPE={bw} Hz/px)。注意这不是实际位移——扫描仪重设了 f0，"
                          f"常数项大部分被 topup 吸收；真正建模不了的是一/二阶差异。")

    print("\n=== 实际匀场结果一览（shim 电流 + 中心频率）===")
    _print_shim_list(shimgrp)


def print_scan_groups(spec):
    """同规格 EPI 的扫描框分组。func 和它的反向 PE 图必须落在同一组，
    否则 topup/SDC 拿到的是两个几何不同的输入。"""
    print("\n=== EPI 扫描框分组（func 与反向 PE 图应落在同一组）===")
    for (ro, pe, th, ns), members in spec.items():
        sub = {}
        for r in members:
            sub.setdefault(r["scan_group"], []).append(r)
        print(f"  规格 {ro:.0f}×{pe:.0f}×{th:.0f}mm / {ns} 层：")
        if len(sub) == 1:
            print(f"    ✓ 全部 {len(members)} 个 series 定位一致")
            continue
        print(f"    ⚠ 分成 {len(sub)} 组")
        base = None
        for gid, ms in sorted(sub.items()):
            sc = ms[0]["scan"]
            names = ", ".join(f"{m['series']}:{m['desc'][:26]}" for m in ms)
            print(f"      扫描框组#{gid}  center={sc['center'].round(2)} "
                  f"rot={sc['rot']:.3f}°")
            print(f"        {names}")
            if base is None:
                base = sc
            else:
                d = float(np.linalg.norm(sc["center"] - base["center"]))
                a = math.degrees(math.acos(max(-1, min(1, abs(float(
                    np.dot(sc["normal"], base["normal"])))))))
                print(f"        Δvs组1: 中心 {d:.2f}mm  法向 {a:.2f}°  "
                      f"面内 {rot_delta(sc['rot'], base['rot']):.2f}°")


def _print_shim_list(shimgrp):
    for (shim, f0), members in shimgrp.items():
        gid = members[0]["shim_group"]
        print(f"  shim组#{gid}  f0={f0}  shim[0:9]={list(shim)}")
        print(f"            ← series {','.join(str(m['series']) for m in members)}")


def write_out(rows, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    cols = ["series", "desc", "nfiles", "n_slices", "status", "coupled", "shim_mode",
            "adj_prot_id", "shim_group", "f0", "offset_mm", "along_mm", "inplane_mm",
            "normal_deg", "rot_vs_group_deg", "shimbox_ro", "shimbox_pe", "shimbox_th",
            "scanbox_ro", "scanbox_pe", "scanbox_th", "reason"]
    tsv = outdir / "qc_shim.tsv"
    with tsv.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            if not r.get("has_shim"):
                fh.write(f"{r['series']}\t{r['desc']}\t{r['nfiles']}\tn/a\tskip"
                         + "\tn/a" * 16 + f"\t{r.get('note','')}\n")
                continue
            c, sh, sc = r["cmp"], r["shimbox"], r["scan"]
            fh.write("\t".join(str(x) for x in [
                r["series"], r["desc"], r["nfiles"], r["n_slices"], r["status"],
                r["coupled"], r["shim_mode"], r["adj_prot_id"], r["shim_group"], r["f0"],
                f"{c['dist']:.3f}", f"{c['along']:.3f}", f"{c['inplane']:.3f}",
                f"{c['angle']:.4f}", f"{r.get('rot_vs_group',0):.4f}",
                f"{sh['ro']:.1f}", f"{sh['pe']:.1f}", f"{sh['thick']:.1f}",
                f"{sc['ro']:.1f}", f"{sc['pe']:.1f}", f"{sc['thick']:.1f}",
                "; ".join(r["reason"]) or "n/a"]) + "\n")

    js = outdir / "qc_shim.json"
    payload = []
    for r in rows:
        e = {k: r[k] for k in ("series", "desc", "nfiles", "status") if k in r}
        e["reason"] = r.get("reason", [])
        if r.get("has_shim"):
            e.update(coupled=r["coupled"], shim_mode=r["shim_mode"],
                     adj_prot_id=r["adj_prot_id"], shim_group=r["shim_group"],
                     f0=r["f0"], shim_current=r["shim"],
                     shimbox={k: (v.tolist() if isinstance(v, np.ndarray) else v)
                              for k, v in r["shimbox"].items()},
                     scanbox={k: (v.tolist() if isinstance(v, np.ndarray) else v)
                              for k, v in r["scan"].items()},
                     compare=r["cmp"], rot_vs_group=r.get("rot_vs_group"))
        else:
            e["note"] = r.get("note", "")
        payload.append(e)
    js.write_text(json.dumps(payload, indent=2, ensure_ascii=False))
    print(f"\n已写出：\n  {tsv}\n  {js}")


def do_verify(rows):
    """用标准 DICOM 字段独立复算，核对 ASCCONV 的坐标约定。

    只挑 2D 多层序列：3D 序列（如 MP2RAGE）的 asSlice[0].sPosition 存的是整个
    slab 的中心，而 DICOM 首帧是第一个 partition，两者本来就不该相等，拿来核对
    只会得出假的偏差。
    """
    print("=== 坐标约定核对（ASCCONV vs 标准 DICOM 字段）===")
    cands = [r for r in rows
             if r.get("has_shim") and "ds" in r and r["scan"]["n_slices"] > 1
             and all(hasattr(r["ds"], a) for a in
                     ("ImageOrientationPatient", "ImagePositionPatient", "PixelSpacing"))]
    if not cands:
        print("  没有可核对的 2D 多层序列（3D 序列的 asSlice[0] 是 slab 中心，不可比）。")
        return

    d_pos, d_ang = [], []
    for r in cands[:3]:
        ds, d = r["ds"], r["ascconv"]
        cen = slice_center_from_dicom(ds)
        asc = _v(d, "sSliceArray.asSlice[0].sPosition")
        iop = np.array(ds.ImageOrientationPatient, float)
        nrm_dcm = np.cross(iop[:3], iop[3:])
        nrm_asc = r["scan"]["normal"]
        dp = float(np.linalg.norm(asc - cen))
        da = math.degrees(math.acos(
            max(-1, min(1, abs(float(np.dot(nrm_asc, nrm_dcm)))))))
        d_pos.append(dp)
        d_ang.append(da)
        print(f"\n  Series {r['series']}  {r['desc'][:46]}  ({r['scan']['n_slices']} 层 2D)")
        print(f"    首层中心  ASCCONV {asc.round(3)}   DICOM {cen.round(3)}   差 {dp:.3f} mm")
        print(f"    层法向    ASCCONV {nrm_asc.round(5)}   DICOM {nrm_dcm.round(5)}   "
              f"夹角 {da:.4f}°")
        if float(np.dot(nrm_asc, nrm_dcm)) < 0:
            print(f"    （法向整体反号：同一条轴、层序方向相反，不影响几何比较）")

    mp, ma = max(d_pos), max(d_ang)
    print(f"\n  最大偏差：位置 {mp:.3f} mm，法向 {ma:.4f}°")
    if mp < 1.0 and ma < 0.1:
        print("  → 一致。(dSag,dCor,dTra) 即 DICOM LPS 的 (x,y,z)，无需坐标转换。")
        print("    亚毫米残差属正常：ASCCONV 存 FOV 中心，DICOM 算像素矩阵中心，"
              "有面内旋转时两者差半个体素量级。")
    else:
        print("  → 偏差偏大，本机型的坐标约定可能与预期不同，用本工具的结果前请先核实。")


def main():
    ap = argparse.ArgumentParser(
        description="匀场框/扫描框几何与匀场结果一致性检查（输入 DICOM，非 BIDS）")
    ap.add_argument("--dicom", required=True, help="DICOM 目录（一场扫描）")
    ap.add_argument("-o", "--out", help="输出目录，写 qc_shim.tsv / qc_shim.json")
    ap.add_argument("--max-offset", type=float, default=5.0,
                    help="series 内匀场框中心偏移上限 mm，超过判 fail（默认 5）")
    ap.add_argument("--max-angle", type=float, default=3.0,
                    help="匀场框法向夹角上限 度，超过判 fail（默认 3）")
    ap.add_argument("--max-rot", type=float, default=0.5,
                    help="同定位组内面内朝向差上限 度，超过判 warn（默认 0.5）")
    ap.add_argument("--verify", action="store_true", help="打印坐标约定核对后退出")
    args = ap.parse_args()

    root = Path(args.dicom).expanduser()
    if not root.is_dir():
        raise SystemExit(f"不是目录：{root}")

    rows = collect(root)
    if not any(r.get("has_shim") for r in rows):
        vendors = sorted({str(getattr(r["ds"], "Manufacturer", "?")).strip()
                          for r in rows if "ds" in r} or {"?"})
        models = sorted({str(getattr(r["ds"], "ManufacturerModelName", "?")).strip()
                         for r in rows if "ds" in r})
        print(f"没有任何 series 带匀场框记录。")
        print(f"  检测到的厂商/机型：{', '.join(vendors)} / {', '.join(models)}")
        if not any("SIEMENS" in v.upper() for v in vendors):
            print("  → 非西门子设备。本工具依赖 Siemens CSA 私有头 (0029,1020) 的 "
                  "ASCCONV，其他厂商没有等价的公开结构：")
            print("     GE (0043,xxxx) 一般不存独立 shim box（auto-shim 直接基于处方 FOV）；")
            print("     Philips (2005,xxxx) 通常只有匀场模式没有几何；")
            print("     UIH / 联影未见公开的匀场框字段。")
        else:
            print("  → 是西门子但读不到。多半是 XA 系列（私有标签迁到 (0021,xxxx)，"
                  "且导出时可能被剥掉），或匿名化去掉了私有组。")
        return 0

    if args.verify:
        do_verify(rows)
        return 0

    geo, shimgrp, spec, has_fail = judge(rows, args.max_offset, args.max_angle,
                                         args.max_rot)
    print_table(rows)
    print_groups(geo, shimgrp)
    if spec:
        print_scan_groups(spec)

    n_f = sum(1 for r in rows if r["status"] == "fail")
    n_w = sum(1 for r in rows if r["status"] == "warn")
    n_p = sum(1 for r in rows if r["status"] == "pass")
    print(f"\n汇总：pass {n_p} / warn {n_w} / fail {n_f} "
          f"（另有 {sum(1 for r in rows if r['status']=='skip')} 个 series 无匀场框记录）")

    if args.out:
        write_out(rows, Path(args.out).expanduser())
    return 1 if has_fail else 0


if __name__ == "__main__":
    sys.exit(main())
