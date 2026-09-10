"""Text/CSV report and matplotlib figures for GNSS QC."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .constants import IOD_JUMP_THRES, SYS_NAME
from .quality import QCResult, SystemQC
from .rinex_reader import RinexObs, summarize_header


def format_report(obs: RinexObs, qc: QCResult) -> str:
    lines: List[str] = []
    lines.append("=" * 72)
    lines.append("GNSS Observation Data Quality Assessment")
    lines.append("Reference: BD 420022-2019 北斗/GNSS测量型接收机观测数据质量评估方法")
    lines.append("=" * 72)
    lines.append("")
    lines.append(summarize_header(obs.header))
    lines.append(f"Epochs loaded: {qc.n_epochs_file}")
    lines.append("")

    for sys, sq in qc.systems.items():
        pair = qc.pairs.get(sys)
        lines.append("-" * 72)
        lines.append(f"System: {sys} ({sq.name})")
        if pair:
            lines.append(
                f"  Dual-freq pair: {pair.c1}/{pair.l1} + {pair.c2}/{pair.l2} "
                f"(bands {pair.band1}/{pair.band2})"
            )
        slip_txt = f"{sq.slip_ratio:.1f}" if np.isfinite(sq.slip_ratio) else "Inf (no slips)"
        lines.append(f"  Data completeness DIs (dual-freq): {sq.completeness:.2f} %")
        lines.append(f"  Observed epochs (sum over sats): {sq.n_epochs_obs}")
        lines.append(f"  Cycle-slip events: {sq.n_slip_epochs}")
        lines.append(f"  Observations per slip (o/slps): {slip_txt}")
        lines.append(f"  Multipath RMS MP1: {sq.mp1_rms:.4f} m" if np.isfinite(sq.mp1_rms) else "  Multipath RMS MP1: N/A")
        lines.append(f"  Multipath RMS MP2: {sq.mp2_rms:.4f} m" if np.isfinite(sq.mp2_rms) else "  Multipath RMS MP2: N/A")
        lines.append(f"  IOD RMS: {sq.iod_rms:.5f} m/s" if np.isfinite(sq.iod_rms) else "  IOD RMS: N/A")
        lines.append(f"  IOD jumps (|IOD|>{IOD_JUMP_THRES} m/s): {sq.n_iod_jumps}")
        if sq.pr_noise:
            pr = ", ".join(f"{k}={v:.4f} m" for k, v in sq.pr_noise.items())
            lines.append(f"  Pseudorange noise: {pr}")
        if sq.ph_noise:
            ph = ", ".join(f"{k}={v:.5f} cyc" for k, v in sq.ph_noise.items())
            lines.append(f"  Carrier-phase noise: {ph}")
        if sq.cnr_mean:
            cn = ", ".join(f"{k}={v:.2f} dBHz" for k, v in sq.cnr_mean.items())
            lines.append(f"  Mean CNR: {cn}")
        lines.append("  Per-satellite summary:")
        lines.append(
            f"    {'SAT':<5} {'N':>6} {'Slip':>5} {'MP1':>8} {'MP2':>8} {'IOD':>9} {'IOD#':>5}"
        )
        for sat, s in sorted(sq.sats.items()):
            lines.append(
                f"    {sat:<5} {s.n_epochs:6d} {s.n_slips:5d} "
                f"{_f(s.mp1_rms,8,3)} {_f(s.mp2_rms,8,3)} {_f(s.iod_rms,9,5)} {s.n_iod_jumps:5d}"
            )
        lines.append("")

    lines.append("=" * 72)
    lines.append("Notes:")
    lines.append("  - Completeness uses dual-frequency valid epochs / theoretical epochs.")
    lines.append("  - Cycle slips: MW recursive test + GF consecutive/poly residual (BD 6.2).")
    lines.append("  - Multipath: dual-freq MP combination with sliding-window detrending (BD 6.3).")
    lines.append("  - IOD: dI/dt from dual-freq phase ionosphere proxy (BD 6.4).")
    lines.append("  - Noise: triple-difference RMS (BD 6.5 / 6.6).")
    lines.append("=" * 72)
    return "\n".join(lines)


def _f(v: float, width: int, prec: int) -> str:
    if not np.isfinite(v):
        return f"{'N/A':>{width}}"
    return f"{v:{width}.{prec}f}"


def write_csv_summary(qc: QCResult, path: Path) -> None:
    rows = ["sys,sat,n_epochs,n_complete,n_slips,n_outliers,mp1_rms_m,mp2_rms_m,iod_rms_mps,n_iod_jumps"]
    for sys, sq in qc.systems.items():
        for sat, s in sorted(sq.sats.items()):
            rows.append(
                f"{sys},{sat},{s.n_epochs},{s.n_complete},{s.n_slips},{s.n_outliers},"
                f"{s.mp1_rms},{s.mp2_rms},{s.iod_rms},{s.n_iod_jumps}"
            )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def _pick_example_sats(sq: SystemQC, k: int = 4) -> List[str]:
    ranked = sorted(sq.sats.items(), key=lambda kv: kv[1].n_complete, reverse=True)
    return [s for s, _ in ranked[:k]]


def make_plots(obs: RinexObs, qc: QCResult, outdir: Path, max_sats_per_sys: int = 4) -> List[Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    paths: List[Path] = []

    # Overview bar charts
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    sys_list = list(qc.systems.keys())
    names = [SYS_NAME.get(s, s) for s in sys_list]
    comp = [qc.systems[s].completeness for s in sys_list]
    slips = []
    for s in sys_list:
        sr = qc.systems[s].slip_ratio
        if np.isfinite(sr):
            slips.append(sr)
        else:
            # No slips detected: show observed epochs as a high "good" ratio
            slips.append(float(max(qc.systems[s].n_epochs_obs, 1)))
    mp1 = [qc.systems[s].mp1_rms for s in sys_list]
    iod = [qc.systems[s].iod_rms for s in sys_list]

    axes[0, 0].bar(names, comp, color="#2a6f97")
    axes[0, 0].set_ylabel("Completeness (%)")
    axes[0, 0].set_title("Data Completeness")
    axes[0, 0].tick_params(axis="x", rotation=20)

    axes[0, 1].bar(names, slips, color="#bc4749")
    axes[0, 1].set_ylabel("Obs / Slip")
    axes[0, 1].set_title("Cycle-Slip Ratio (higher is better)")
    axes[0, 1].tick_params(axis="x", rotation=20)

    axes[1, 0].bar(names, [m if np.isfinite(m) else 0 for m in mp1], color="#6a994e")
    axes[1, 0].set_ylabel("MP1 RMS (m)")
    axes[1, 0].set_title("Multipath MP1 RMS")
    axes[1, 0].tick_params(axis="x", rotation=20)

    axes[1, 1].bar(names, [m if np.isfinite(m) else 0 for m in iod], color="#e09f3e")
    axes[1, 1].set_ylabel("IOD RMS (m/s)")
    axes[1, 1].set_title("Ionospheric Delay Rate RMS")
    axes[1, 1].tick_params(axis="x", rotation=20)

    fig.suptitle("GNSS QC Overview (BD 420022-2019)", fontsize=13)
    fig.tight_layout()
    p = outdir / "qc_overview.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    paths.append(p)

    # Per-system time series for top satellites
    for sys, sq in qc.systems.items():
        sats = _pick_example_sats(sq, max_sats_per_sys)
        if not sats:
            continue
        fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
        for sat in sats:
            s = sq.sats[sat]
            if s.t is None:
                continue
            tmin = s.t / 60.0
            if s.mp1 is not None:
                axes[0].plot(tmin, s.mp1, lw=0.8, label=sat)
            if s.iod is not None:
                axes[1].plot(tmin, s.iod, lw=0.8, label=sat)
            if s.gf is not None:
                axes[2].plot(tmin, s.gf, lw=0.8, label=sat)
                if s.slip_flags is not None and np.any(s.slip_flags):
                    axes[2].scatter(
                        tmin[s.slip_flags],
                        s.gf[s.slip_flags],
                        c="red",
                        s=18,
                        zorder=5,
                        marker="x",
                    )
        axes[0].set_ylabel("MP1 (m)")
        axes[0].set_title(f"{SYS_NAME.get(sys, sys)} multipath / IOD / GF")
        axes[0].legend(loc="upper right", fontsize=8, ncol=4)
        axes[0].grid(True, alpha=0.3)

        axes[1].axhline(IOD_JUMP_THRES, color="r", ls="--", lw=0.8)
        axes[1].axhline(-IOD_JUMP_THRES, color="r", ls="--", lw=0.8)
        axes[1].set_ylabel("IOD (m/s)")
        axes[1].grid(True, alpha=0.3)

        axes[2].set_ylabel("GF (m)")
        axes[2].set_xlabel("Time (min)")
        axes[2].grid(True, alpha=0.3)

        fig.tight_layout()
        p = outdir / f"qc_series_{sys}.png"
        fig.savefig(p, dpi=140)
        plt.close(fig)
        paths.append(p)

    # CNR distribution if available
    fig, ax = plt.subplots(figsize=(9, 4.5))
    plotted = False
    for sys, sq in qc.systems.items():
        vals = []
        for s in sq.sats.values():
            for v in s.cnr_mean.values():
                if np.isfinite(v):
                    vals.append(v)
        if vals:
            ax.hist(vals, bins=15, alpha=0.55, label=SYS_NAME.get(sys, sys))
            plotted = True
    if plotted:
        ax.set_xlabel("Mean CNR (dBHz)")
        ax.set_ylabel("Satellite count")
        ax.set_title("Carrier-to-Noise Ratio Distribution")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        p = outdir / "qc_cnr.png"
        fig.savefig(p, dpi=140)
        paths.append(p)
    plt.close(fig)

    return paths
