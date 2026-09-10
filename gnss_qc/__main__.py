#!/usr/bin/env python3
"""CLI: python -m gnss_qc <rinex.26o> [--out Result_QC] [--max-epochs N]"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .generate_sample import generate_sample_rinex
from .quality import analyze
from .report import format_report, make_plots, write_csv_summary
from .rinex_reader import read_rinex_obs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="GNSS RINEX observation quality analysis (BD 420022-2019)"
    )
    parser.add_argument(
        "rinex",
        nargs="?",
        default=None,
        help="RINEX 3 observation file (.26o / .rnx). If omitted, a sample is generated.",
    )
    parser.add_argument(
        "--out",
        default="Result_QC",
        help="Output directory for report, CSV and plots (default: Result_QC)",
    )
    parser.add_argument(
        "--max-epochs",
        type=int,
        default=None,
        help="Optional limit on number of epochs to read (for quick tests)",
    )
    parser.add_argument(
        "--systems",
        default=None,
        help="Comma-separated systems to analyze, e.g. G,C,E (default: all in file)",
    )
    parser.add_argument(
        "--sample-epochs",
        type=int,
        default=180,
        help="When generating sample file, number of epochs (default: 180)",
    )
    args = parser.parse_args(argv)

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.rinex is None:
        sample = Path("data") / "sample_yygc1.26o"
        print(f"No RINEX given; generating sample -> {sample}")
        generate_sample_rinex(sample, n_epochs=args.sample_epochs)
        rinex_path = sample
    else:
        rinex_path = Path(args.rinex)
        if not rinex_path.exists():
            print(f"ERROR: file not found: {rinex_path}", file=sys.stderr)
            return 1

    print(f"Reading {rinex_path} ...")
    obs = read_rinex_obs(rinex_path, max_epochs=args.max_epochs)
    print(f"  header lines={obs.header.header_lines}, epochs={len(obs.epochs)}")
    print(f"  systems={list(obs.header.obs_types.keys())}, sats={len(obs.sat_list())}")

    systems = None
    if args.systems:
        systems = [s.strip().upper() for s in args.systems.split(",") if s.strip()]

    print("Running quality assessment (BD 420022-2019) ...")
    qc = analyze(obs, systems=systems)

    report = format_report(obs, qc)
    report_path = outdir / "qc_report.txt"
    report_path.write_text(report, encoding="utf-8")
    print(report)
    print(f"\nReport saved: {report_path}")

    csv_path = outdir / "qc_summary.csv"
    write_csv_summary(qc, csv_path)
    print(f"CSV saved:    {csv_path}")

    print("Generating plots ...")
    plots = make_plots(obs, qc, outdir)
    for p in plots:
        print(f"  plot: {p}")

    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
