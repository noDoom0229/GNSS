"""Quality metrics per BD 420022-2019."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from .constants import (
    C_LIGHT,
    C_MS,
    CLOCK_JUMP_XI,
    DUAL_PAIRS,
    FREQ,
    GF_LAMBDA_FACTOR,
    IOD_JUMP_THRES,
    MP_WINDOW,
    MW_SIGMA_FACTOR,
    SYS_NAME,
    band_from_obs,
    glonass_freq,
    wavelength,
)
from .rinex_reader import RinexObs


@dataclass
class PairConfig:
    sys: str
    band1: str
    band2: str
    c1: str
    l1: str
    c2: str
    l2: str
    s1: Optional[str] = None
    s2: Optional[str] = None
    f1: float = 0.0
    f2: float = 0.0


@dataclass
class SatQC:
    sat: str
    n_epochs: int = 0
    n_complete: int = 0
    n_slips: int = 0
    n_outliers: int = 0
    n_iod_jumps: int = 0
    mp1_rms: float = float("nan")
    mp2_rms: float = float("nan")
    iod_rms: float = float("nan")
    pr_noise: Dict[str, float] = field(default_factory=dict)
    ph_noise: Dict[str, float] = field(default_factory=dict)
    cnr_mean: Dict[str, float] = field(default_factory=dict)
    # time series for plotting
    t: Optional[np.ndarray] = None
    mp1: Optional[np.ndarray] = None
    mp2: Optional[np.ndarray] = None
    iod: Optional[np.ndarray] = None
    gf: Optional[np.ndarray] = None
    mw: Optional[np.ndarray] = None
    slip_flags: Optional[np.ndarray] = None


@dataclass
class SystemQC:
    sys: str
    name: str
    completeness: float = float("nan")  # %
    slip_ratio: float = float("nan")    # obs/slip (o/slps)
    n_epochs_obs: int = 0
    n_slip_epochs: int = 0
    n_iod_jumps: int = 0
    mp1_rms: float = float("nan")
    mp2_rms: float = float("nan")
    iod_rms: float = float("nan")
    pr_noise: Dict[str, float] = field(default_factory=dict)
    ph_noise: Dict[str, float] = field(default_factory=dict)
    cnr_mean: Dict[str, float] = field(default_factory=dict)
    sats: Dict[str, SatQC] = field(default_factory=dict)


@dataclass
class QCResult:
    systems: Dict[str, SystemQC]
    interval: float
    n_epochs_file: int
    clock_jumps: List[float] = field(default_factory=list)
    pairs: Dict[str, PairConfig] = field(default_factory=dict)


def _find_code(types: List[str], kind: str, band: str) -> Optional[str]:
    """Find first observation code of kind C/L/S/D for given band digit."""
    for t in types:
        if t.startswith(kind) and len(t) >= 2 and t[1] == band:
            return t
    return None


def select_pair(obs: RinexObs, sys: str) -> Optional[PairConfig]:
    types = obs.header.obs_types.get(sys, [])
    if not types:
        return None
    for b1, b2 in DUAL_PAIRS.get(sys, []):
        c1 = _find_code(types, "C", b1)
        l1 = _find_code(types, "L", b1)
        c2 = _find_code(types, "C", b2)
        l2 = _find_code(types, "L", b2)
        if c1 and l1 and c2 and l2:
            s1 = _find_code(types, "S", b1)
            s2 = _find_code(types, "S", b2)
            f1 = FREQ.get(sys, {}).get(b1, 0.0)
            f2 = FREQ.get(sys, {}).get(b2, 0.0)
            return PairConfig(sys, b1, b2, c1, l1, c2, l2, s1, s2, f1, f2)
    return None


def _freq_for_sat(obs: RinexObs, sat: str, band: str, default: float) -> float:
    if sat[0] != "R":
        return default
    try:
        slot = int(sat[1:])
    except ValueError:
        return default
    ch = obs.header.glo_channels.get(slot, 0)
    return glonass_freq(band, ch)


def _phase_m(L_cycles: np.ndarray, freq: float) -> np.ndarray:
    return L_cycles * wavelength(freq)


def mw_combination(phi1: np.ndarray, phi2: np.ndarray, rho1: np.ndarray, rho2: np.ndarray,
                   f1: float, f2: float) -> np.ndarray:
    """Melbourne-Wübbena combination in meters (BD formula 3)."""
    return (f1 * phi1 - f2 * phi2) / (f1 - f2) - (f1 * rho1 + f2 * rho2) / (f1 + f2)


def gf_combination(phi1: np.ndarray, phi2: np.ndarray) -> np.ndarray:
    """Geometry-free phase combination (m): L_GF = φ2 - φ1."""
    return phi2 - phi1


def multipath(rho1, rho2, phi1, phi2, f1, f2) -> Tuple[np.ndarray, np.ndarray]:
    """Pseudorange multipath combinations MP1, MP2 (BD formula 16)."""
    a = (f1 ** 2 + f2 ** 2) / (f1 ** 2 - f2 ** 2)
    b = 2 * f2 ** 2 / (f1 ** 2 - f2 ** 2)
    # MP1 = ρ1 - ((f1^2+f2^2)/(f1^2-f2^2)) φ1 + (2 f2^2/(f1^2-f2^2)) φ2
    mp1 = rho1 - a * phi1 + b * phi2
    # MP2 = ρ2 - ((f1^2+f2^2)/(f2^2-f1^2)) φ2 + (2 f1^2/(f2^2-f1^2)) φ1
    a2 = (f1 ** 2 + f2 ** 2) / (f2 ** 2 - f1 ** 2)
    b2 = 2 * f1 ** 2 / (f2 ** 2 - f1 ** 2)
    mp2 = rho2 - a2 * phi2 + b2 * phi1
    return mp1, mp2


def ionosphere_delay(phi1, phi2, f1, f2) -> Tuple[np.ndarray, np.ndarray]:
    """Ionospheric delay proxies I1, I2 from dual-frequency phase (BD formula 18)."""
    i1 = (f2 ** 2) / (f1 ** 2 - f2 ** 2) * (phi1 - phi2)
    i2 = (f1 ** 2) / (f1 ** 2 - f2 ** 2) * (phi1 - phi2)
    return i1, i2


def _detrend_mp(mp: np.ndarray, slip: np.ndarray, window: int = MP_WINDOW) -> np.ndarray:
    """Remove ambiguity bias per arc via sliding-window mean (BD formula 17)."""
    out = np.full_like(mp, np.nan, dtype=float)
    n = len(mp)
    # Split into arcs at slips / gaps
    valid = np.isfinite(mp)
    arc_id = np.zeros(n, dtype=int)
    cur = 0
    for i in range(n):
        if i > 0 and (slip[i] or not valid[i - 1]):
            cur += 1
        arc_id[i] = cur
    for a in np.unique(arc_id):
        idx = np.where(arc_id == a)[0]
        if len(idx) < 3:
            continue
        vals = mp[idx].copy()
        # centered moving average
        half = max(1, window // 2)
        smooth = np.full_like(vals, np.nan)
        for k in range(len(vals)):
            lo = max(0, k - half)
            hi = min(len(vals), k + half + 1)
            seg = vals[lo:hi]
            seg = seg[np.isfinite(seg)]
            if len(seg) >= 3:
                smooth[k] = vals[k] - np.mean(seg)
        out[idx] = smooth
    return out


def detect_mw_slips(mw: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    MW recursive mean/variance cycle-slip & outlier detection (BD 6.2.2).
    Returns (slip_flags, outlier_flags) boolean arrays aligned with mw.
    """
    n = len(mw)
    slips = np.zeros(n, dtype=bool)
    outliers = np.zeros(n, dtype=bool)
    if n == 0:
        return slips, outliers

    mean = 0.0
    var = 0.0
    count = 0
    i = 0
    while i < n:
        if not np.isfinite(mw[i]):
            # gap resets arc
            if count > 0:
                mean = 0.0
                var = 0.0
                count = 0
            i += 1
            continue
        if count == 0:
            mean = mw[i]
            var = 0.0
            count = 1
            i += 1
            continue
        sigma = np.sqrt(var) if var > 0 else 0.0
        # bootstrap sigma from early samples
        thr = MW_SIGMA_FACTOR * max(sigma, 0.5)  # floor 0.5 m ≈ 0.6 cycle WL
        if abs(mw[i] - mean) >= thr and count >= 5:
            # look ahead one epoch to distinguish outlier vs slip
            j = i + 1
            while j < n and not np.isfinite(mw[j]):
                j += 1
            if j >= n:
                slips[i] = True
                mean = mw[i]
                var = 0.0
                count = 1
                i += 1
                continue
            # predicted mean/var including next without i
            mean_next = ((count - 0) * mean + mw[j]) / (count + 1) if False else None
            # BD: compute mean/var from previous arc for epoch i+1
            mean_ip1 = (count * mean + mw[j]) / (count + 1)
            var_ip1 = ((count - 1) / count) * var + (mw[j] - mean) ** 2 / (count + 1) if count > 0 else 0.0
            sigma_ip1 = np.sqrt(max(var_ip1, 0.0))
            thr_ip1 = MW_SIGMA_FACTOR * max(sigma_ip1, 0.5)
            i_over = True
            ip1_over = abs(mw[j] - mean) >= thr  # vs previous mean
            # simpler practical rule (BD intent):
            # if next also over and |mw[i]-mw[j]| small -> slip; else if next ok or jump back -> outlier
            if abs(mw[j] - mean) < thr:
                outliers[i] = True
                # skip outlier, keep arc
                i = j  # will process j normally next? better skip i
                i = i  # noqa
                i += 1
                continue
            if abs(mw[i] - mw[j]) <= thr and abs(mw[j] - mean) >= thr:
                slips[i] = True
                mean = mw[i]
                var = 0.0
                count = 1
                i += 1
                continue
            # default: treat as slip start of new arc
            slips[i] = True
            mean = mw[i]
            var = 0.0
            count = 1
            i += 1
            continue
        # update recursive mean / variance
        count += 1
        delta = mw[i] - mean
        mean += delta / count
        if count > 1:
            var = ((count - 2) / (count - 1)) * var + (delta ** 2) / count if count > 2 else (delta ** 2) / 2
            # stable Welford-like:
        # recompute with Welford for robustness
        # (overwrite with proper Welford on the fly)
        i += 1

    # Second pass with clean Welford + slip flags (more reliable)
    slips[:] = False
    outliers[:] = False
    mean = 0.0
    m2 = 0.0
    count = 0
    for i in range(n):
        if not np.isfinite(mw[i]):
            mean = 0.0
            m2 = 0.0
            count = 0
            continue
        if count == 0:
            mean = mw[i]
            m2 = 0.0
            count = 1
            continue
        sigma = np.sqrt(m2 / (count - 1)) if count > 1 else 0.0
        thr = MW_SIGMA_FACTOR * max(sigma, 0.5)
        if abs(mw[i] - mean) >= thr and count >= 5:
            # peek next valid
            j = i + 1
            while j < n and not np.isfinite(mw[j]):
                j += 1
            if j < n and abs(mw[j] - mean) < thr:
                outliers[i] = True
                continue
            if j < n and abs(mw[i] - mw[j]) <= max(thr, 1.0) and abs(mw[j] - mean) >= thr:
                slips[i] = True
            else:
                slips[i] = True
            mean = mw[i]
            m2 = 0.0
            count = 1
            continue
        count += 1
        delta = mw[i] - mean
        mean += delta / count
        m2 += delta * (mw[i] - mean)
    return slips, outliers


def detect_gf_slips(gf: np.ndarray, lam1: float, lam2: float, already: np.ndarray) -> np.ndarray:
    """
    Supplemental GF cycle-slip detection (BD 6.2.3 simplified).
    Uses consecutive GF difference vs 6*(λ1-λ2) when poly fit is short;
    for longer arcs uses polynomial residual jumps.
    """
    n = len(gf)
    slips = already.copy()
    thr = GF_LAMBDA_FACTOR * abs(lam1 - lam2)
    # consecutive difference on arcs without existing slips
    prev = np.nan
    for i in range(n):
        if not np.isfinite(gf[i]):
            prev = np.nan
            continue
        if slips[i]:
            prev = gf[i]
            continue
        if np.isfinite(prev) and abs(gf[i] - prev) > thr:
            slips[i] = True
        prev = gf[i]

    # Polynomial residual check on longer continuous arcs
    arc_start = 0
    while arc_start < n:
        while arc_start < n and (not np.isfinite(gf[arc_start]) or slips[arc_start]):
            arc_start += 1
        arc_end = arc_start
        while arc_end < n and np.isfinite(gf[arc_end]) and not slips[arc_end]:
            arc_end += 1
        length = arc_end - arc_start
        if length >= 30:
            x = np.arange(length, dtype=float)
            y = gf[arc_start:arc_end]
            q = 6 if length / 100.0 >= 6 else max(1, int(length / 100.0) + 1)
            q = min(q, length - 1)
            try:
                coef = np.polyfit(x, y, q)
                qfit = np.polyval(coef, x)
                resid = y - qfit
                # jump in residual consecutive difference
                for k in range(1, length - 1):
                    d0 = abs(resid[k] - resid[k - 1])
                    d1 = abs(resid[k + 1] - resid[k])
                    if d0 > thr and d1 > thr * 0.5:
                        slips[arc_start + k] = True
            except Exception:
                pass
        arc_start = max(arc_end, arc_start + 1)
    return slips


def triple_diff_noise(x: np.ndarray, slip: Optional[np.ndarray] = None) -> float:
    """RMS of triple differences / sqrt(8*(N-1)/N) ≈ formula (21)/(23)."""
    if slip is not None:
        x = x.copy()
        x[slip] = np.nan
    # process continuous finite segments
    vals = []
    seg: List[float] = []
    for v in x:
        if np.isfinite(v):
            seg.append(float(v))
        else:
            if len(seg) >= 4:
                a = np.asarray(seg)
                d1 = np.diff(a)
                d2 = np.diff(d1)
                d3 = np.diff(d2)
                vals.extend(d3.tolist())
            seg = []
    if len(seg) >= 4:
        a = np.asarray(seg)
        d1 = np.diff(a)
        d2 = np.diff(d1)
        d3 = np.diff(d2)
        vals.extend(d3.tolist())
    if len(vals) < 2:
        return float("nan")
    arr = np.asarray(vals, dtype=float)
    n = len(arr)
    # BD 420022-2019 eq.(21)/(23): σ = sqrt( Σ(Δ³)² / (8*(N-1)) )
    return float(np.sqrt(np.sum(arr ** 2) / (8.0 * (n - 1))))

def _rms(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return float("nan")
    return float(np.sqrt(np.mean(x ** 2)))


def analyze_satellite(obs: RinexObs, sat: str, pair: PairConfig) -> SatQC:
    qc = SatQC(sat=sat)
    epochs = obs.epochs
    n = len(epochs)
    t0 = epochs[0].time
    t = np.full(n, np.nan)
    rho1 = np.full(n, np.nan)
    rho2 = np.full(n, np.nan)
    L1c = np.full(n, np.nan)
    L2c = np.full(n, np.nan)
    s1 = np.full(n, np.nan)
    s2 = np.full(n, np.nan)

    f1 = _freq_for_sat(obs, sat, pair.band1, pair.f1)
    f2 = _freq_for_sat(obs, sat, pair.band2, pair.f2)
    if f1 <= 0 or f2 <= 0:
        return qc

    for i, ep in enumerate(epochs):
        t[i] = (ep.time - t0).total_seconds()
        if sat not in ep.data:
            continue
        row = ep.data[sat]
        rho1[i] = row.get(pair.c1, np.nan)
        rho2[i] = row.get(pair.c2, np.nan)
        L1c[i] = row.get(pair.l1, np.nan)
        L2c[i] = row.get(pair.l2, np.nan)
        if pair.s1:
            s1[i] = row.get(pair.s1, np.nan)
        if pair.s2:
            s2[i] = row.get(pair.s2, np.nan)

    dual = np.isfinite(rho1) & np.isfinite(rho2) & np.isfinite(L1c) & np.isfinite(L2c)
    qc.n_epochs = int(np.sum(np.isfinite(rho1) | np.isfinite(L1c)))
    qc.n_complete = int(np.sum(dual))
    if qc.n_complete < 5:
        return qc

    phi1 = _phase_m(L1c, f1)
    phi2 = _phase_m(L2c, f2)
    mw = mw_combination(phi1, phi2, rho1, rho2, f1, f2)
    gf = gf_combination(phi1, phi2)
    mw[~dual] = np.nan
    gf[~dual] = np.nan

    slips, outliers = detect_mw_slips(mw)
    slips = detect_gf_slips(gf, wavelength(f1), wavelength(f2), slips)
    qc.n_slips = int(np.sum(slips))
    qc.n_outliers = int(np.sum(outliers))
    qc.slip_flags = slips

    # Multipath
    mp1_raw, mp2_raw = multipath(rho1, rho2, phi1, phi2, f1, f2)
    mp1_raw[~dual] = np.nan
    mp2_raw[~dual] = np.nan
    mp1 = _detrend_mp(mp1_raw, slips)
    mp2 = _detrend_mp(mp2_raw, slips)
    qc.mp1 = mp1
    qc.mp2 = mp2
    qc.mp1_rms = _rms(mp1)
    qc.mp2_rms = _rms(mp2)

    # Ionosphere / IOD
    i1, _i2 = ionosphere_delay(phi1, phi2, f1, f2)
    i1[~dual] = np.nan
    iod = np.full(n, np.nan)
    for i in range(1, n):
        if slips[i]:
            continue
        if np.isfinite(i1[i]) and np.isfinite(i1[i - 1]) and (t[i] - t[i - 1]) > 0:
            # skip across slip at i-1 boundary
            if slips[i - 1]:
                continue
            iod[i] = (i1[i] - i1[i - 1]) / (t[i] - t[i - 1])
    qc.iod = iod
    qc.iod_rms = _rms(iod)
    qc.n_iod_jumps = int(np.sum(np.abs(iod[np.isfinite(iod)]) > IOD_JUMP_THRES))

    # Noise
    qc.pr_noise[pair.c1] = triple_diff_noise(rho1)
    qc.pr_noise[pair.c2] = triple_diff_noise(rho2)
    qc.ph_noise[pair.l1] = triple_diff_noise(L1c, slips)  # cycles
    qc.ph_noise[pair.l2] = triple_diff_noise(L2c, slips)

    # CNR
    if pair.s1:
        qc.cnr_mean[pair.s1] = float(np.nanmean(s1)) if np.any(np.isfinite(s1)) else float("nan")
    if pair.s2:
        qc.cnr_mean[pair.s2] = float(np.nanmean(s2)) if np.any(np.isfinite(s2)) else float("nan")

    qc.t = t
    qc.gf = gf
    qc.mw = mw
    return qc


def theoretical_epochs(obs: RinexObs) -> int:
    hdr = obs.header
    if hdr.time_first and hdr.time_last and hdr.interval > 0:
        span = (hdr.time_last - hdr.time_first).total_seconds()
        return int(round(span / hdr.interval)) + 1
    return len(obs.epochs)


def analyze(obs: RinexObs, systems: Optional[List[str]] = None) -> QCResult:
    """Run BD 420022-2019 quality assessment on a RINEX observation dataset."""
    systems = systems or list(obs.header.obs_types.keys())
    pairs: Dict[str, PairConfig] = {}
    result_systems: Dict[str, SystemQC] = {}
    n_theo = theoretical_epochs(obs)

    for sys in systems:
        pair = select_pair(obs, sys)
        if pair is None:
            continue
        pairs[sys] = pair
        sq = SystemQC(sys=sys, name=SYS_NAME.get(sys, sys))
        sats = obs.sat_list(sys)
        mp1s, mp2s, iods = [], [], []
        pr_acc: Dict[str, List[float]] = defaultdict(list)
        ph_acc: Dict[str, List[float]] = defaultdict(list)
        cnr_acc: Dict[str, List[float]] = defaultdict(list)
        total_obs = 0
        total_slips = 0
        total_complete = 0
        total_iod_j = 0

        for sat in sats:
            sat_qc = analyze_satellite(obs, sat, pair)
            if sat_qc.n_epochs == 0:
                continue
            sq.sats[sat] = sat_qc
            total_obs += sat_qc.n_epochs
            total_slips += sat_qc.n_slips
            total_complete += sat_qc.n_complete
            total_iod_j += sat_qc.n_iod_jumps
            if np.isfinite(sat_qc.mp1_rms):
                mp1s.append(sat_qc.mp1_rms)
            if np.isfinite(sat_qc.mp2_rms):
                mp2s.append(sat_qc.mp2_rms)
            if np.isfinite(sat_qc.iod_rms):
                iods.append(sat_qc.iod_rms)
            for k, v in sat_qc.pr_noise.items():
                if np.isfinite(v):
                    pr_acc[k].append(v)
            for k, v in sat_qc.ph_noise.items():
                if np.isfinite(v):
                    ph_acc[k].append(v)
            for k, v in sat_qc.cnr_mean.items():
                if np.isfinite(v):
                    cnr_acc[k].append(v)

        # Completeness: use dual-freq complete epochs vs theoretical * n_sats approx
        # Per standard: sum actual / sum theoretical over sats
        n_sats = max(len(sq.sats), 1)
        theo_sum = n_theo * n_sats
        sq.completeness = 100.0 * total_complete / theo_sum if theo_sum > 0 else float("nan")
        sq.n_epochs_obs = total_obs
        sq.n_slip_epochs = total_slips
        sq.slip_ratio = (total_obs / total_slips) if total_slips > 0 else float("inf")
        sq.n_iod_jumps = total_iod_j
        sq.mp1_rms = float(np.mean(mp1s)) if mp1s else float("nan")
        sq.mp2_rms = float(np.mean(mp2s)) if mp2s else float("nan")
        sq.iod_rms = float(np.mean(iods)) if iods else float("nan")
        sq.pr_noise = {k: float(np.mean(v)) for k, v in pr_acc.items()}
        sq.ph_noise = {k: float(np.mean(v)) for k, v in ph_acc.items()}
        sq.cnr_mean = {k: float(np.mean(v)) for k, v in cnr_acc.items()}
        result_systems[sys] = sq

    return QCResult(
        systems=result_systems,
        interval=obs.header.interval,
        n_epochs_file=len(obs.epochs),
        pairs=pairs,
    )
