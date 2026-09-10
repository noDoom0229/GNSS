"""RINEX 3.x observation file reader (.YYo / .rnx)."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

from .constants import SYS_NAME


@dataclass
class RinexHeader:
    version: float = 3.0
    file_type: str = "O"
    marker_name: str = ""
    approx_pos: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    interval: float = 1.0
    time_first: Optional[datetime] = None
    time_last: Optional[datetime] = None
    leap_seconds: int = 0
    # sys -> ordered list of obs types, e.g. {"G": ["C1C", "L1C", ...]}
    obs_types: Dict[str, List[str]] = field(default_factory=dict)
    # GLONASS slot -> frequency channel number
    glo_channels: Dict[int, int] = field(default_factory=dict)
    comments: List[str] = field(default_factory=list)
    header_lines: int = 0


@dataclass
class Epoch:
    time: datetime
    flag: int
    sats: List[str]
    # sat -> obs_code -> value (NaN if missing)
    data: Dict[str, Dict[str, float]]
    # sat -> obs_code -> LLI (optional)
    lli: Dict[str, Dict[str, int]] = field(default_factory=dict)


@dataclass
class RinexObs:
    header: RinexHeader
    epochs: List[Epoch]

    @property
    def systems(self) -> List[str]:
        return sorted(self.header.obs_types.keys())

    def sat_list(self, sys: Optional[str] = None) -> List[str]:
        sats = set()
        for ep in self.epochs:
            for s in ep.sats:
                if sys is None or s[0] == sys:
                    sats.add(s)
        return sorted(sats)

    def times(self) -> np.ndarray:
        return np.array([ep.time for ep in self.epochs], dtype=object)


def _parse_float(token: str) -> float:
    token = token.strip()
    if not token or token == "0.000":
        # RINEX often uses blank for missing; "0.000" may be real zero (Doppler)
        # Caller decides; here parse literally.
        try:
            return float(token) if token else float("nan")
        except ValueError:
            return float("nan")
    try:
        return float(token)
    except ValueError:
        return float("nan")


def _label(line: str) -> str:
    return line[60:].rstrip("\r\n") if len(line) >= 60 else ""


def _parse_obs_types_block(lines: List[str], start: int) -> Tuple[str, List[str], int]:
    """Parse SYS / # / OBS TYPES starting at `start`. Returns (sys, types, next_index)."""
    line = lines[start]
    sys = line[0].strip()
    n = int(line[3:6])
    types: List[str] = []
    # First line: up to 13 types after columns 7-58
    chunk = line[7:60]
    types.extend([chunk[i : i + 4].strip() for i in range(0, len(chunk), 4) if chunk[i : i + 4].strip()])
    idx = start + 1
    while len(types) < n and idx < len(lines):
        cont = lines[idx]
        if "SYS / # / OBS TYPES" not in _label(cont) and cont[0].strip():
            break
        chunk = cont[7:60]
        types.extend([chunk[i : i + 4].strip() for i in range(0, len(chunk), 4) if chunk[i : i + 4].strip()])
        idx += 1
    return sys, types[:n], idx


def parse_header(lines: List[str]) -> Tuple[RinexHeader, int]:
    hdr = RinexHeader()
    i = 0
    while i < len(lines):
        line = lines[i].rstrip("\n")
        lab = _label(line)
        if lab == "END OF HEADER":
            hdr.header_lines = i + 1
            return hdr, i + 1
        if lab == "RINEX VERSION / TYPE":
            hdr.version = float(line[0:9].strip() or 3.0)
            hdr.file_type = line[20:21].strip() or "O"
        elif lab == "MARKER NAME":
            hdr.marker_name = line[0:60].strip()
        elif lab == "APPROX POSITION XYZ":
            try:
                hdr.approx_pos = (
                    float(line[0:14]),
                    float(line[14:28]),
                    float(line[28:42]),
                )
            except ValueError:
                pass
        elif lab == "INTERVAL":
            try:
                hdr.interval = float(line[0:10].strip())
            except ValueError:
                pass
        elif lab == "TIME OF FIRST OBS":
            try:
                y, m, d, hh, mm = [int(line[j : j + 6]) for j in range(0, 30, 6)]
                sec = float(line[30:43])
                hdr.time_first = datetime(y, m, d, hh, mm) + timedelta(seconds=sec)
            except Exception:
                pass
        elif lab == "TIME OF LAST OBS":
            try:
                y, m, d, hh, mm = [int(line[j : j + 6]) for j in range(0, 30, 6)]
                sec = float(line[30:43])
                hdr.time_last = datetime(y, m, d, hh, mm) + timedelta(seconds=sec)
            except Exception:
                pass
        elif lab == "LEAP SECONDS":
            try:
                hdr.leap_seconds = int(line[0:6].strip())
            except ValueError:
                pass
        elif lab == "SYS / # / OBS TYPES":
            sys, types, nxt = _parse_obs_types_block(lines, i)
            hdr.obs_types[sys] = types
            i = nxt
            continue
        elif lab == "GLONASS SLOT / FRQ #":
            # Format: nsat then pairs of slot/channel; may span multiple lines
            try:
                nsat = int(line[0:3])
            except ValueError:
                nsat = 0
            tokens = line[4:60].split()
            pairs: List[str] = list(tokens)
            j = i + 1
            while len(pairs) < 2 * nsat and j < len(lines) and "GLONASS SLOT / FRQ #" in _label(lines[j]):
                pairs.extend(lines[j][4:60].split())
                j += 1
            for k in range(0, min(len(pairs) - 1, 2 * nsat), 2):
                try:
                    slot = int(pairs[k])
                    ch = int(pairs[k + 1])
                    hdr.glo_channels[slot] = ch
                except ValueError:
                    continue
            i = j
            continue
        elif lab == "COMMENT":
            hdr.comments.append(line[0:60].rstrip())
        i += 1
    hdr.header_lines = len(lines)
    return hdr, len(lines)


def _parse_epoch_header(line: str) -> Tuple[datetime, int, int]:
    # > YYYY MM DD HH MM SS.sssssss  flag  nsat
    parts = line[1:].split()
    y, m, d, hh, mm = map(int, parts[0:5])
    sec = float(parts[5])
    flag = int(parts[6])
    nsat = int(parts[7])
    t = datetime(y, m, d, hh, mm) + timedelta(seconds=sec)
    return t, flag, nsat


def _parse_obs_line(line: str, obs_types: List[str]) -> Tuple[str, Dict[str, float], Dict[str, int]]:
    sat = line[0:3].strip()
    values: Dict[str, float] = {}
    lli_map: Dict[str, int] = {}
    # Each observation field is 16 chars: 14 value + 1 LLI + 1 SSI
    for i, code in enumerate(obs_types):
        start = 3 + i * 16
        field = line[start : start + 16] if start < len(line) else ""
        if len(field) < 14 or not field[:14].strip():
            values[code] = float("nan")
            continue
        values[code] = _parse_float(field[:14])
        if len(field) >= 15 and field[14].strip():
            try:
                lli_map[code] = int(field[14])
            except ValueError:
                pass
    return sat, values, lli_map


def read_rinex_obs(path: str | Path, max_epochs: Optional[int] = None) -> RinexObs:
    """Read a RINEX 3 observation file into memory."""
    path = Path(path)
    text = path.read_text(encoding="utf-8", errors="replace")
    # Normalize CRLF
    lines = text.splitlines()
    header, idx = parse_header(lines)
    if not header.obs_types:
        raise ValueError(f"No SYS / # / OBS TYPES found in {path}")

    epochs: List[Epoch] = []
    i = idx
    nlines = len(lines)
    while i < nlines:
        line = lines[i]
        if not line.startswith(">"):
            i += 1
            continue
        try:
            t, flag, nsat = _parse_epoch_header(line)
        except Exception:
            i += 1
            continue
        i += 1
        sats: List[str] = []
        data: Dict[str, Dict[str, float]] = {}
        lli: Dict[str, Dict[str, int]] = {}
        for _ in range(nsat):
            if i >= nlines:
                break
            sat_line = lines[i]
            i += 1
            if not sat_line or sat_line.startswith(">"):
                i -= 1
                break
            sys = sat_line[0]
            types = header.obs_types.get(sys, [])
            sat, vals, lli_s = _parse_obs_line(sat_line, types)
            sats.append(sat)
            data[sat] = vals
            if lli_s:
                lli[sat] = lli_s
        # Special event epochs (flag>1) may have fewer sat lines; still store
        epochs.append(Epoch(time=t, flag=flag, sats=sats, data=data, lli=lli))
        if max_epochs is not None and len(epochs) >= max_epochs:
            break

    if header.time_first is None and epochs:
        header.time_first = epochs[0].time
    if header.time_last is None and epochs:
        header.time_last = epochs[-1].time
    return RinexObs(header=header, epochs=epochs)


def extract_series(
    obs: RinexObs,
    sat: str,
    codes: Iterable[str],
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """
    Extract time series for one satellite and observation codes.
    Returns (times as float seconds from first epoch, {code: values}).
    """
    t0 = obs.epochs[0].time if obs.epochs else None
    times: List[float] = []
    series: Dict[str, List[float]] = {c: [] for c in codes}
    for ep in obs.epochs:
        if sat not in ep.data:
            continue
        dt = (ep.time - t0).total_seconds() if t0 else 0.0
        times.append(dt)
        row = ep.data[sat]
        for c in codes:
            series[c].append(row.get(c, float("nan")))
    return np.asarray(times, dtype=float), {c: np.asarray(v, dtype=float) for c, v in series.items()}


def summarize_header(header: RinexHeader) -> str:
    lines = [
        f"RINEX {header.version:.2f}  marker={header.marker_name or '-'}",
        f"Approx XYZ: {header.approx_pos}",
        f"Interval: {header.interval} s",
        f"First: {header.time_first}  Last: {header.time_last}",
        f"Leap seconds: {header.leap_seconds}",
        "Observation types:",
    ]
    for sys, types in header.obs_types.items():
        name = SYS_NAME.get(sys, sys)
        lines.append(f"  {sys} ({name}): {len(types)} -> {', '.join(types)}")
    if header.glo_channels:
        lines.append(f"GLONASS channels: {len(header.glo_channels)} slots")
    return "\n".join(lines)
