"""Generate a synthetic RINEX 3.04 observation file for demos/tests."""

from __future__ import annotations

import math
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

from .constants import C_LIGHT, FREQ, wavelength


OBS_TYPES = {
    "G": ["C1C", "C2W", "C5X", "L1C", "L2W", "L5X", "D1C", "D2W", "D5X", "S1C", "S2W", "S5X"],
    "R": ["C1C", "C2C", "L1C", "L2C", "S1C", "S2C", "D1C", "D2C"],
    "C": ["C2I", "C7I", "C6I", "L2I", "L7I", "L6I", "S2I", "S7I", "S6I", "D2I", "D7I", "D6I"],
    "E": ["C1X", "C5X", "C7X", "L1X", "L5X", "L7X", "S1X", "S5X", "S7X", "D1X", "D5X", "D7X"],
    "J": ["C1C", "C2X", "C5X", "L1C", "L2X", "L5X", "D1C", "D2X", "D5X", "S1C", "S2X"],
}

# sats used in sample
SATS = {
    "G": ["G05", "G15", "G24", "G29"],
    "R": ["R09", "R10", "R22"],
    "C": ["C01", "C06", "C11", "C23"],
    "E": ["E07", "E13", "E26"],
    "J": ["J03"],
}

GLO_CH = {9: -2, 10: -7, 22: 0}


def _fmt_obs_types(sys: str, types: list) -> list:
    lines = []
    n = len(types)
    # first line up to 13 types
    first = types[:13]
    body = f"{sys}  {n:3d}" + "".join(f" {t:<3s}" for t in first)
    body = body.ljust(60) + "SYS / # / OBS TYPES"
    lines.append(body)
    rest = types[13:]
    while rest:
        chunk = rest[:13]
        rest = rest[13:]
        body = "      " + "".join(f" {t:<3s}" for t in chunk)
        body = body.ljust(60) + "SYS / # / OBS TYPES"
        lines.append(body)
    return lines


def _range_geom(t_sec: float, sat: str) -> float:
    """Simple smooth geometric range model (~2.2e7 m) with sat-dependent offset."""
    seed = sum(ord(c) for c in sat)
    r0 = 22.0e6 + (seed % 500) * 1e3
    # slow variation
    return r0 + 800.0 * math.sin(2 * math.pi * t_sec / 3600.0 + seed) + 50.0 * t_sec * 0.01


def _iono(t_sec: float, sat: str) -> float:
    seed = sum(ord(c) for c in sat)
    return 3.0 + 0.4 * math.sin(2 * math.pi * t_sec / 1800.0 + seed / 10.0)


def generate_sample_rinex(
    path: str | Path,
    n_epochs: int = 180,
    interval: float = 1.0,
    start: datetime | None = None,
    inject_slips: bool = True,
) -> Path:
    """Write a RINEX 3.04 mixed observation file resembling the user's .26o layout."""
    path = Path(path)
    start = start or datetime(2026, 7, 3, 7, 4, 4)
    end = start + timedelta(seconds=(n_epochs - 1) * interval)

    lines: list[str] = []
    lines.append(f"{3.04:9.2f}           OBSERVATION DATA    M (MIXED)           RINEX VERSION / TYPE")
    lines.append(f"{'GNSS QC SAMPLE':<40s}{'PC':<20s}PGM / RUN BY / DATE")
    lines.append(f"{'Synthetic data for BD420022 QC demo':<60s}COMMENT")
    lines.append(f"{'yygc1':<60s}MARKER NAME")
    lines.append(f"{'':<60s}MARKER NUMBER")
    lines.append(f"{'CURSOR':<20s}{'AGENT':<20s}{'':<20s}OBSERVER / AGENCY")
    lines.append(f"{'SIMREC':<20s}{'GENERIC':<20s}{'':<20s}REC # / TYPE / VERS")
    lines.append(f"{'ANT000':<20s}{'NONE':<20s}{'':<20s}ANT # / TYPE")
    x, y, z = -2245742.9269, 5064323.2954, 3150178.4754
    lines.append(f"{x:14.4f}{y:14.4f}{z:14.4f}                  APPROX POSITION XYZ")
    lines.append(f"{0.0:14.4f}{0.0:14.4f}{0.0:14.4f}                  ANTENNA: DELTA H/E/N")
    for sys, types in OBS_TYPES.items():
        lines.extend(_fmt_obs_types(sys, types))
    lines.append(f"{interval:10.3f}                                                  INTERVAL")
    lines.append(
        f"{start.year:6d}{start.month:6d}{start.day:6d}{start.hour:6d}{start.minute:6d}"
        f"{start.second + start.microsecond/1e6:13.7f}     GPS         TIME OF FIRST OBS"
    )
    lines.append(
        f"{end.year:6d}{end.month:6d}{end.day:6d}{end.hour:6d}{end.minute:6d}"
        f"{end.second + end.microsecond/1e6:13.7f}     GPS         TIME OF LAST OBS"
    )
    # GLONASS channels
    slots = sorted(GLO_CH.keys())
    glo = f"{len(slots):3d}"
    for s in slots:
        glo += f" {s:2d} {GLO_CH[s]:2d}"
    lines.append(glo.ljust(60) + "GLONASS SLOT / FRQ #")
    lines.append(f"{18:6d}                                                      LEAP SECONDS")
    lines.append(f"{'':<60s}END OF HEADER")

    rng = np.random.default_rng(42)

    # Maintain phase state per sat/band
    phase_state: dict[tuple[str, str], float] = {}

    all_sats = [s for sys in SATS for s in SATS[sys]]

    for k in range(n_epochs):
        t = start + timedelta(seconds=k * interval)
        t_sec = k * interval
        # epoch header
        lines.append(
            f"> {t.year:4d} {t.month:2d} {t.day:2d} {t.hour:2d} {t.minute:2d} "
            f"{t.second + t.microsecond/1e6:10.7f}  0 {len(all_sats):2d}"
        )
        for sat in all_sats:
            sys = sat[0]
            types = OBS_TYPES[sys]
            geom = _range_geom(t_sec, sat)
            iono = _iono(t_sec, sat)
            # mild multipath
            mp = 0.15 * math.sin(2 * math.pi * t_sec / 120.0 + sum(ord(c) for c in sat))
            elev_proxy = 30 + 20 * math.sin(t_sec / 600.0 + ord(sat[-1]))

            fields = []
            # Build dual-freq values depending on system
            if sys == "G":
                bands = [("1", "C1C", "L1C", "D1C", "S1C"),
                         ("2", "C2W", "L2W", "D2W", "S2W"),
                         ("5", "C5X", "L5X", "D5X", "S5X")]
            elif sys == "R":
                bands = [("1", "C1C", "L1C", "D1C", "S1C"),
                         ("2", "C2C", "L2C", "D2C", "S2C")]
            elif sys == "C":
                bands = [("2", "C2I", "L2I", "D2I", "S2I"),
                         ("7", "C7I", "L7I", "D7I", "S7I"),
                         ("6", "C6I", "L6I", "D6I", "S6I")]
            elif sys == "E":
                bands = [("1", "C1X", "L1X", "D1X", "S1X"),
                         ("5", "C5X", "L5X", "D5X", "S5X"),
                         ("7", "C7X", "L7X", "D7X", "S7X")]
            else:  # J
                bands = [("1", "C1C", "L1C", "D1C", "S1C"),
                         ("2", "C2X", "L2X", "D2X", "S2X"),
                         ("5", "C5X", "L5X", "D5X", "S5X")]

            values: dict[str, float] = {}
            for bi, (band, c_code, l_code, d_code, s_code) in enumerate(bands):
                f = FREQ[sys][band]
                # iono scale 1/f^2
                iono_m = iono * (FREQ[sys][bands[0][0]] ** 2) / (f ** 2)
                pr_noise = rng.normal(0, 0.3)
                ph_noise = rng.normal(0, 0.003)  # cycles
                rho = geom + iono_m + (mp if bi == 0 else -mp * 0.3) + pr_noise
                key = (sat, band)
                if key not in phase_state:
                    phase_state[key] = (geom - iono_m) / wavelength(f)
                # advance phase approximately with geometry change
                phase_state[key] = (geom - iono_m) / wavelength(f) + ph_noise
                # inject cycle slip
                if inject_slips and sat in ("G05", "C06") and k in (60, 120):
                    if bi == 0:
                        phase_state[key] += 5.0  # 5-cycle slip on L1
                L = phase_state[key]
                dop = rng.normal(1000 if bi == 0 else -500, 50)
                snr = 35 + elev_proxy / 5.0 + rng.normal(0, 0.8) - bi * 1.5
                values[c_code] = rho
                values[l_code] = L
                values[d_code] = dop
                values[s_code] = snr

            # format line
            line = f"{sat:<3s}"
            for code in types:
                if code in values and math.isfinite(values[code]):
                    v = values[code]
                    if code.startswith("L"):
                        line += f"{v:14.3f}  "
                    elif code.startswith("S"):
                        line += f"{v:14.3f}  "
                    elif code.startswith("D"):
                        line += f"{v:14.3f}  "
                    else:
                        line += f"{v:14.3f}  "
                else:
                    line += " " * 16
            lines.append(line.rstrip())

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path


if __name__ == "__main__":
    out = Path(__file__).resolve().parents[1] / "data" / "sample_yygc1.26o"
    p = generate_sample_rinex(out, n_epochs=180)
    print(f"Wrote {p} ({p.stat().st_size} bytes)")
