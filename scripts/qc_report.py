#!/usr/bin/env python3
"""
qc_report.py — 把 qc_raw.py 的结果渲染成单文件 HTML 报告

被 qc_raw.py 自动调用；也可以单独跑（从 qc_raw.json + motion/*.tsv 重建）：
    python3 qc_report.py /path/to/derivatives/qc-raw
"""
from __future__ import annotations

import base64
import io
import json
import sys
from pathlib import Path

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False

STATUS_COLOR = {"pass": "#1f9d55", "warn": "#c98a00", "fail": "#d33"}


def _fig_to_data_uri(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=96, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def plot_run(r, cfg):
    """三联图：6 参数 / FD+maxdisp / outlier+DVARS。"""
    if not HAVE_MPL or r.get("_bids_mp") is None:
        return None
    mp = np.asarray(r["_bids_mp"])
    fd = np.asarray(r.get("_fd", []), dtype=float)
    md = np.asarray(r.get("_maxdisp", []), dtype=float)
    of = np.asarray(r.get("_outfrac", []), dtype=float)
    dv = np.asarray(r.get("_dvars", []), dtype=float)
    t = np.arange(mp.shape[0])

    fig, axes = plt.subplots(3, 1, figsize=(11, 6.2), sharex=True,
                             gridspec_kw={"height_ratios": [2, 1.3, 1.3]})
    ax = axes[0]
    for i, lab in enumerate(["trans_x", "trans_y", "trans_z"]):
        ax.plot(t, mp[:, i], lw=0.9, label=lab)
    ax.set_ylabel("translation (mm)")
    ax.legend(ncol=3, fontsize=7, loc="upper right", framealpha=0.6)
    ax2 = ax.twinx()
    for i, lab in enumerate(["rot_x", "rot_y", "rot_z"]):
        ax2.plot(t, np.degrees(mp[:, 3 + i]), lw=0.7, ls="--", alpha=0.75, label=lab)
    ax2.set_ylabel("rotation (deg)")
    ax2.legend(ncol=3, fontsize=7, loc="lower right", framealpha=0.6)

    ax = axes[1]
    if fd.size:
        ax.plot(t[:fd.size], fd, lw=0.9, color="#333", label="FD (Power, mm)")
        ax.axhline(cfg["fd_thresh"], color="#d33", lw=0.8, ls=":",
                   label=f"thr {cfg['fd_thresh']}")
        bad = np.where(fd > cfg["fd_thresh"])[0]
        if bad.size:
            ax.scatter(bad, fd[bad], s=6, color="#d33", zorder=3)
    if md.size:
        ax.plot(t[:md.size], md, lw=0.7, color="#0a7", alpha=0.8, label="maxdisp (mm)")
    ax.set_ylabel("displacement (mm)")
    ax.legend(ncol=3, fontsize=7, loc="upper right", framealpha=0.6)

    ax = axes[2]
    if of.size:
        ax.plot(t[:of.size], of, lw=0.9, color="#36c", label="outlier fraction")
        ax.axhline(cfg["outlier_thresh"], color="#d33", lw=0.8, ls=":")
    ax.set_ylabel("outlier")
    ax.set_xlabel("TR")
    if dv.size > 1:
        ax3 = ax.twinx()
        ax3.plot(t[:dv.size], dv, lw=0.7, color="#a5a", alpha=0.8, label="DVARS")
        ax3.set_ylabel("DVARS")
        ax3.legend(fontsize=7, loc="lower right", framealpha=0.6)
    ax.legend(fontsize=7, loc="upper right", framealpha=0.6)

    fig.suptitle(r["name"], fontsize=9)
    fig.tight_layout()
    return _fig_to_data_uri(fig)


CSS = """
:root{--bg:#fff;--fg:#1b1b1b;--mut:#666;--line:#e3e3e3;--card:#fafafa}
@media (prefers-color-scheme:dark){
 :root{--bg:#15171a;--fg:#e6e6e6;--mut:#9aa;--line:#2c3036;--card:#1c1f24}}
*{box-sizing:border-box}
body{margin:0;padding:24px;background:var(--bg);color:var(--fg);
 font:14px/1.55 -apple-system,BlinkMacSystemFont,"PingFang SC","Helvetica Neue",sans-serif}
h1{font-size:20px;margin:0 0 4px} h2{font-size:16px;margin:28px 0 10px}
.meta{color:var(--mut);font-size:12px;margin-bottom:18px}
.cards{display:flex;gap:10px;flex-wrap:wrap;margin:12px 0 22px}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;
 padding:10px 16px;min-width:110px}
.card .n{font-size:22px;font-weight:600} .card .l{font-size:11px;color:var(--mut)}
.tw{overflow-x:auto;border:1px solid var(--line);border-radius:8px}
table{border-collapse:collapse;width:100%;font-size:12px;white-space:nowrap}
th,td{padding:5px 9px;border-bottom:1px solid var(--line);text-align:left}
th{background:var(--card);position:sticky;top:0;font-weight:600}
tr.fail td{background:rgba(221,51,51,.10)} tr.warn td{background:rgba(201,138,0,.12)}
.badge{display:inline-block;padding:1px 7px;border-radius:9px;color:#fff;font-size:11px}
details{border:1px solid var(--line);border-radius:8px;margin:8px 0;background:var(--card)}
summary{cursor:pointer;padding:8px 12px;font-size:13px}
details img{width:100%;max-width:1100px;display:block;padding:0 12px 12px}
.reason{color:#b35;font-size:11px;white-space:normal;max-width:520px}
code{background:var(--card);padding:1px 5px;border-radius:4px;font-size:12px}
"""

RUN_TABLE_COLS = [
    ("name", "run"), ("status", "状态"), ("n_volumes", "TR数"),
    ("expected_volumes", "期望"), ("tr", "TR(s)"), ("vox_str", "体素"),
    ("fd_mean", "FD均值"), ("fd_max", "FD最大"), ("fd_max_tr", "峰值TR"),
    ("n_fd_spikes", "急动次数"), ("pct_fd_gt_thresh", "FD超阈%"),
    ("maxdisp_range", "run内位移"), ("betrun_disp_vs_ref", "run间位移"),
    ("pct_out_gt_thresh", "outlier超阈%"), ("tsnr_median", "tSNR"),
    ("n_nonsteady_est", "非稳态帧"), ("reason", "原因"),
]


def _cell(r, k):
    v = r.get(k)
    if v is None or v == "":
        return "n/a" if k != "reason" else ""
    if isinstance(v, float):
        return f"{v:.3f}" if abs(v) < 100 else f"{v:.1f}"
    return str(v)


def build_report(out_dir: Path, runs, subj_rows, cfg):
    out_dir = Path(out_dir)
    n = len(runs)
    nf = sum(1 for r in runs if r["status"] == "fail")
    nw = sum(1 for r in runs if r["status"] == "warn")

    h = ["<style>", CSS, "</style>",
         "<h1>原始数据 QC 报告 · qc-raw</h1>",
         f"<div class=meta>数据源 <code>{out_dir}</code> ｜ "
         f"仅刚体 EPI-to-EPI 头动估计，未做任何解剖/模板配准</div>",
         "<div class=cards>",
         f"<div class=card><div class=n>{n}</div><div class=l>run 总数</div></div>",
         f"<div class=card><div class=n style='color:{STATUS_COLOR['pass']}'>"
         f"{n-nf-nw}</div><div class=l>pass</div></div>",
         f"<div class=card><div class=n style='color:{STATUS_COLOR['warn']}'>"
         f"{nw}</div><div class=l>warn</div></div>",
         f"<div class=card><div class=n style='color:{STATUS_COLOR['fail']}'>"
         f"{nf}</div><div class=l>fail</div></div>",
         f"<div class=card><div class=n>{len(subj_rows)}</div>"
         f"<div class=l>被试(×session)</div></div>",
         "</div>"]

    # 被试表
    h.append("<h2>被试级</h2><div class=tw><table><tr>")
    scols = ["subject", "session", "status", "n_runs", "n_runs_expected",
             "missing_runs", "protocol_mismatch", "n_anat", "n_fmap", "fd_mean_worst",
             "tsnr_median_min", "betrun_disp_max", "tasks", "reason"]
    h += [f"<th>{c}</th>" for c in scols] + ["</tr>"]
    for s in subj_rows:
        h.append(f"<tr class={s['status']}>")
        for c in scols:
            cls = " class=reason" if c == "reason" else ""
            if c == "status":
                h.append(f"<td><span class=badge style='background:"
                         f"{STATUS_COLOR[s['status']]}'>{s['status']}</span></td>")
            else:
                h.append(f"<td{cls}>{_cell(s, c)}</td>")
        h.append("</tr>")
    h.append("</table></div>")

    # run 表
    h.append("<h2>run 级</h2><div class=tw><table><tr>")
    h += [f"<th>{lab}</th>" for _, lab in RUN_TABLE_COLS] + ["</tr>"]
    for r in sorted(runs, key=lambda x: (x["subject"], x["session"],
                                         x["task"], x["run"], x["name"])):
        h.append(f"<tr class={r['status']}>")
        for k, _ in RUN_TABLE_COLS:
            if k == "status":
                h.append(f"<td><span class=badge style='background:"
                         f"{STATUS_COLOR[r['status']]}'>{r['status']}</span></td>")
            else:
                cls = " class=reason" if k == "reason" else ""
                h.append(f"<td{cls}>{_cell(r, k)}</td>")
        h.append("</tr>")
    h.append("</table></div>")

    # 曲线
    h.append("<h2>逐 run 运动曲线</h2>")
    if not HAVE_MPL:
        h.append("<p class=meta>未安装 matplotlib，跳过绘图。</p>")
    for r in sorted(runs, key=lambda x: ({"fail": 0, "warn": 1, "pass": 2}[x["status"]],
                                         x["name"])):
        uri = plot_run(r, cfg)
        openattr = " open" if r["status"] in ("fail", "warn") else ""
        badge = (f"<span class=badge style='background:{STATUS_COLOR[r['status']]}'>"
                 f"{r['status']}</span>")
        h.append(f"<details{openattr}><summary>{badge} {r['name']} "
                 f"— {r.get('reason','') or '正常'}</summary>")
        h.append(f"<img src='{uri}'>" if uri else "<p class=meta>无曲线</p>")
        h.append("</details>")

    # 阈值
    h.append("<h2>判定阈值</h2><div class=tw><table><tr><th>参数</th><th>值</th></tr>")
    for k, v in cfg.items():
        h.append(f"<tr><td>{k}</td><td>{v}</td></tr>")
    h.append("</table></div>")

    (out_dir / "report.html").write_text(
        f"<!doctype html><meta charset=utf-8><title>qc-raw 报告</title>"
        + "".join(h), encoding="utf-8")


def attach_timeseries(out_dir: Path, runs):
    """把 motion/*.tsv 的曲线读回 run dict，供画图用（--rejudge 与独立重跑都要）。"""
    mdir = Path(out_dir) / "motion"
    for r in runs:
        f = mdir / f"{r['name'].replace('_bold','')}_desc-rawmotion_timeseries.tsv"
        if not f.exists():
            continue
        raw = np.genfromtxt(f, delimiter="\t", names=True, missing_values="n/a",
                            filling_values=np.nan)
        cols = raw.dtype.names
        r["_bids_mp"] = np.column_stack([raw[c] for c in cols[:6]])
        r["_fd"] = raw["framewise_displacement"]
        r["_maxdisp"] = raw["maxdisp"]
        r["_outfrac"] = raw["outlier_fraction"]
        r["_dvars"] = raw["dvars"]
    return runs


def _rebuild_from_disk(out_dir: Path):
    payload = json.loads((out_dir / "qc_raw.json").read_text())
    runs, cfg = payload["Runs"], payload["Thresholds"]
    attach_timeseries(out_dir, runs)
    return runs, payload["Subjects"], cfg


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("用法: qc_report.py <qc-raw 输出目录>")
    d = Path(sys.argv[1]).resolve()
    runs, subj, cfg = _rebuild_from_disk(d)
    build_report(d, runs, subj, cfg)
    print(f"报告: {d/'report.html'}")
