#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
GNSS 观测数据质量分析（可读教学版）
依据：BD 420022-2019《北斗/GNSS测量型接收机观测数据质量评估方法》
--------------------------------------------------------------------------------
本文件把主要算法集中写在一起，方便阅读与对照标准公式。
工程化分包版本见目录 gnss_qc/。

用法：
    python GNSS_QC_可读代码.py data/sample_yygc1.26o
    python GNSS_QC_可读代码.py your.26o --out Result_QC
================================================================================
"""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# 可选绘图（没有 matplotlib 也能出文字报告）
# ---------------------------------------------------------------------------
try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False


# =============================================================================
# 1. 常数与频率（附录 A）
# =============================================================================

C_LIGHT = 299792458.0  # 光速 m/s

# 各系统载波频率 Hz（RINEX 频点数字 -> 频率）
FREQ = {
    "G": {"1": 1575.42e6, "2": 1227.60e6, "5": 1176.45e6},          # GPS L1/L2/L5
    "R": {"1": 1602.0e6, "2": 1246.0e6, "3": 1202.025e6},           # GLONASS（标称）
    "E": {"1": 1575.42e6, "5": 1176.45e6, "7": 1207.14e6, "6": 1278.75e6},
    "C": {"2": 1561.098e6, "7": 1207.14e6, "6": 1268.52e6, "1": 1575.42e6, "5": 1176.45e6},
    "J": {"1": 1575.42e6, "2": 1227.60e6, "5": 1176.45e6},
}

SYS_NAME = {"G": "GPS", "R": "GLONASS", "E": "Galileo", "C": "BDS", "J": "QZSS"}

# 附录 B：优先选用频率差较大的双频组合
DUAL_PAIRS = {
    "G": [("1", "2"), ("1", "5")],       # L1/L2, L1/L5
    "R": [("1", "2")],
    "E": [("1", "5"), ("1", "7")],
    "C": [("2", "6"), ("2", "7")],       # B1I/B3I, B1I/B2
    "J": [("1", "2"), ("1", "5")],
}

MW_K = 4.0          # MW 超限倍数（标准 6.2.2）
GF_K = 6.0          # GF 周跳阈值倍数（标准 6.2.3）
IOD_JUMP = 0.07     # 电离层延迟变化率跳变阈值 m/s（标准 6.4）
MP_WIN = 50         # 多路径滑动窗口历元数（标准 6.3）


def wavelength(f_hz: float) -> float:
    """波长 λ = c / f，单位 m。"""
    return C_LIGHT / f_hz


def glonass_freq(band: str, channel: int) -> float:
    """GLONASS FDMA：G1=1602+k*9/16 MHz，G2=1246+k*7/16 MHz。"""
    if band == "1":
        return (1602.0 + channel * 9.0 / 16.0) * 1e6
    if band == "2":
        return (1246.0 + channel * 7.0 / 16.0) * 1e6
    return FREQ["R"][band]


# =============================================================================
# 2. RINEX 3.x 读取（.26o）
# =============================================================================

@dataclass
class Header:
    version: float = 3.04
    marker: str = ""
    xyz: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    interval: float = 1.0
    t_first: Optional[datetime] = None
    t_last: Optional[datetime] = None
    obs_types: Dict[str, List[str]] = field(default_factory=dict)  # 系统 -> 观测类型列表
    glo_ch: Dict[int, int] = field(default_factory=dict)           # GLONASS 频道号


@dataclass
class Epoch:
    time: datetime
    flag: int
    data: Dict[str, Dict[str, float]]  # 卫星号 -> {观测类型: 数值}


def _label(line: str) -> str:
    """RINEX 标签在第 61 列起。"""
    return line[60:].rstrip() if len(line) >= 60 else ""


def read_rinex(path: str | Path, max_epochs: Optional[int] = None) -> Tuple[Header, List[Epoch]]:
    """
    读取 RINEX 3 观测文件。
    - 文件头：SYS / # / OBS TYPES 定义每系统列顺序
    - 数据区：'>' 开头为历元行，其后为各卫星观测行（每字段宽 16）
    """
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    hdr = Header()
    i = 0

    # ---------- 解析文件头 ----------
    while i < len(lines):
        line = lines[i]
        lab = _label(line)
        if lab == "END OF HEADER":
            i += 1
            break
        if lab == "RINEX VERSION / TYPE":
            hdr.version = float(line[0:9].strip() or 3.04)
        elif lab == "MARKER NAME":
            hdr.marker = line[:60].strip()
        elif lab == "APPROX POSITION XYZ":
            hdr.xyz = (float(line[0:14]), float(line[14:28]), float(line[28:42]))
        elif lab == "INTERVAL":
            hdr.interval = float(line[0:10].strip())
        elif lab in ("TIME OF FIRST OBS", "TIME OF LAST OBS"):
            y, m, d, hh, mm = [int(line[j:j + 6]) for j in range(0, 30, 6)]
            sec = float(line[30:43])
            t = datetime(y, m, d, hh, mm) + timedelta(seconds=sec)
            if lab.startswith("TIME OF FIRST"):
                hdr.t_first = t
            else:
                hdr.t_last = t
        elif lab == "SYS / # / OBS TYPES":
            sys = line[0].strip()
            n = int(line[3:6])
            types: List[str] = []
            chunk = line[7:60]
            types += [chunk[k:k + 4].strip() for k in range(0, len(chunk), 4) if chunk[k:k + 4].strip()]
            i += 1
            while len(types) < n and i < len(lines):
                chunk = lines[i][7:60]
                types += [chunk[k:k + 4].strip() for k in range(0, len(chunk), 4) if chunk[k:k + 4].strip()]
                i += 1
            hdr.obs_types[sys] = types[:n]
            continue
        elif lab == "GLONASS SLOT / FRQ #":
            # 形如：  3  9 -2 10 -7 22  0
            toks = line[0:60].split()
            try:
                nsat = int(toks[0])
                for k in range(1, 1 + 2 * nsat, 2):
                    hdr.glo_ch[int(toks[k])] = int(toks[k + 1])
            except (ValueError, IndexError):
                pass
        i += 1

    # ---------- 解析历元数据 ----------
    epochs: List[Epoch] = []
    while i < len(lines):
        line = lines[i]
        if not line.startswith(">"):
            i += 1
            continue
        # > YYYY MM DD HH MM SS.sssssss  flag  nsat
        p = line[1:].split()
        y, m, d, hh, mm = map(int, p[0:5])
        sec = float(p[5])
        flag, nsat = int(p[6]), int(p[7])
        t = datetime(y, m, d, hh, mm) + timedelta(seconds=sec)
        i += 1
        data: Dict[str, Dict[str, float]] = {}
        for _ in range(nsat):
            if i >= len(lines) or lines[i].startswith(">"):
                break
            sline = lines[i]
            i += 1
            sat = sline[0:3].strip()
            types = hdr.obs_types.get(sat[0], [])
            vals: Dict[str, float] = {}
            for j, code in enumerate(types):
                # 每个观测值占 16 列：14 数值 + 1 LLI + 1 SSI
                field = sline[3 + j * 16: 3 + j * 16 + 14]
                if not field.strip():
                    vals[code] = float("nan")
                else:
                    try:
                        vals[code] = float(field)
                    except ValueError:
                        vals[code] = float("nan")
            data[sat] = vals
        epochs.append(Epoch(time=t, flag=flag, data=data))
        if max_epochs and len(epochs) >= max_epochs:
            break

    if hdr.t_first is None and epochs:
        hdr.t_first = epochs[0].time
    if hdr.t_last is None and epochs:
        hdr.t_last = epochs[-1].time
    return hdr, epochs


# =============================================================================
# 3. 线性组合（标准公式）
# =============================================================================

def mw_combo(phi1, phi2, rho1, rho2, f1, f2):
    """
    公式(3) Melbourne-Wübbena 组合（单位 m）：
        L_MW = (f1·φ1 - f2·φ2)/(f1-f2) - (f1·ρ1 + f2·ρ2)/(f1+f2)
    其中 φ 为以米为单位的载波相位。
    """
    return (f1 * phi1 - f2 * phi2) / (f1 - f2) - (f1 * rho1 + f2 * rho2) / (f1 + f2)


def gf_combo(phi1, phi2):
    """公式(8) 无几何距离组合：L_GF = φ2 - φ1（m）。"""
    return phi2 - phi1


def multipath_combo(rho1, rho2, phi1, phi2, f1, f2):
    """
    公式(16) 伪距多路径组合：
        MP1 = ρ1 - (f1²+f2²)/(f1²-f2²)·φ1 + 2·f2²/(f1²-f2²)·φ2
        MP2 = ρ2 - (f1²+f2²)/(f2²-f1²)·φ2 + 2·f1²/(f2²-f1²)·φ1
    """
    a = (f1 ** 2 + f2 ** 2) / (f1 ** 2 - f2 ** 2)
    b = 2 * f2 ** 2 / (f1 ** 2 - f2 ** 2)
    mp1 = rho1 - a * phi1 + b * phi2
    a2 = (f1 ** 2 + f2 ** 2) / (f2 ** 2 - f1 ** 2)
    b2 = 2 * f1 ** 2 / (f2 ** 2 - f1 ** 2)
    mp2 = rho2 - a2 * phi2 + b2 * phi1
    return mp1, mp2


def iono_delay(phi1, phi2, f1, f2):
    """
    公式(18) 电离层延迟（含模糊度常数项）：
        I1 = f2²/(f1²-f2²)·(φ1-φ2)
    """
    return (f2 ** 2) / (f1 ** 2 - f2 ** 2) * (phi1 - phi2)


# =============================================================================
# 4. 周跳探测（MW + GF）
# =============================================================================

def detect_mw_slips(mw: np.ndarray) -> np.ndarray:
    """
    标准 6.2.2：对 MW 做递推均值/方差，
    若 |L_MW(i) - 均值| ≥ 4σ，再结合下一历元判断粗差或周跳。
    返回与 mw 同长度的布尔数组（True = 该历元发生周跳）。
    """
    n = len(mw)
    slips = np.zeros(n, dtype=bool)
    mean = m2 = 0.0
    count = 0

    for i in range(n):
        if not np.isfinite(mw[i]):
            mean = m2 = 0.0
            count = 0
            continue
        if count == 0:
            mean, m2, count = mw[i], 0.0, 1
            continue

        sigma = math.sqrt(m2 / (count - 1)) if count > 1 else 0.0
        thr = MW_K * max(sigma, 0.5)  # 0.5 m 作为起步门限

        if abs(mw[i] - mean) >= thr and count >= 5:
            # 看下一有效历元：若回到原弧段 -> 粗差；否则 -> 周跳，新开弧段
            j = i + 1
            while j < n and not np.isfinite(mw[j]):
                j += 1
            if j < n and abs(mw[j] - mean) < thr:
                continue  # 粗差，不重置弧段
            slips[i] = True
            mean, m2, count = mw[i], 0.0, 1
            continue

        # Welford 在线更新均值/方差
        count += 1
        delta = mw[i] - mean
        mean += delta / count
        m2 += delta * (mw[i] - mean)
    return slips


def detect_gf_slips(gf: np.ndarray, lam1: float, lam2: float, slips: np.ndarray) -> np.ndarray:
    """
    标准 6.2.3 补充：相邻历元 GF 跳变超过 6·|λ1-λ2| 则判为周跳。
    """
    thr = GF_K * abs(lam1 - lam2)
    out = slips.copy()
    prev = np.nan
    for i in range(len(gf)):
        if not np.isfinite(gf[i]):
            prev = np.nan
            continue
        if out[i]:
            prev = gf[i]
            continue
        if np.isfinite(prev) and abs(gf[i] - prev) > thr:
            out[i] = True
        prev = gf[i]
    return out


# =============================================================================
# 5. 多路径去模糊度 / 噪声 / IOD
# =============================================================================

def detrend_mp(mp: np.ndarray, slips: np.ndarray, window: int = MP_WIN) -> np.ndarray:
    """
    公式(17)：无周跳弧段内，用滑动窗口均值去掉整周模糊度常数，
    得到多路径误差序列。
    """
    out = np.full_like(mp, np.nan)
    n = len(mp)
    arc = np.zeros(n, dtype=int)
    cur = 0
    for i in range(n):
        if i > 0 and (slips[i] or not np.isfinite(mp[i - 1])):
            cur += 1
        arc[i] = cur

    half = max(1, window // 2)
    for a in np.unique(arc):
        idx = np.where(arc == a)[0]
        if len(idx) < 3:
            continue
        vals = mp[idx]
        smooth = np.full_like(vals, np.nan)
        for k in range(len(vals)):
            seg = vals[max(0, k - half): min(len(vals), k + half + 1)]
            seg = seg[np.isfinite(seg)]
            if len(seg) >= 3:
                smooth[k] = vals[k] - np.mean(seg)
        out[idx] = smooth
    return out


def triple_diff_noise(x: np.ndarray, slips: Optional[np.ndarray] = None) -> float:
    """
    公式(21)/(23)：三次差均方根
        σ = sqrt( Σ(Δ³)² / (8·(N-1)) )
    相位噪声计算前应剔除周跳历元。
    """
    y = x.copy()
    if slips is not None:
        y[slips] = np.nan

    d3_list: List[float] = []
    seg: List[float] = []
    for v in y:
        if np.isfinite(v):
            seg.append(float(v))
        else:
            if len(seg) >= 4:
                a = np.asarray(seg)
                d3_list.extend(np.diff(np.diff(np.diff(a))).tolist())
            seg = []
    if len(seg) >= 4:
        a = np.asarray(seg)
        d3_list.extend(np.diff(np.diff(np.diff(a))).tolist())

    if len(d3_list) < 2:
        return float("nan")
    arr = np.asarray(d3_list)
    return float(np.sqrt(np.sum(arr ** 2) / (8.0 * (len(arr) - 1))))


def rms(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    return float(np.sqrt(np.mean(x ** 2))) if len(x) else float("nan")


# =============================================================================
# 6. 单星 / 单系统质量分析
# =============================================================================

def find_code(types: List[str], kind: str, band: str) -> Optional[str]:
    """在观测类型列表中找某频点的 C/L/S 码，例如 band='1', kind='L' -> L1C。"""
    for t in types:
        if t.startswith(kind) and len(t) >= 2 and t[1] == band:
            return t
    return None


def select_pair(hdr: Header, sys: str):
    """按附录 B 选择该系统的双频观测类型。"""
    types = hdr.obs_types.get(sys, [])
    for b1, b2 in DUAL_PAIRS.get(sys, []):
        c1, l1 = find_code(types, "C", b1), find_code(types, "L", b1)
        c2, l2 = find_code(types, "C", b2), find_code(types, "L", b2)
        if c1 and l1 and c2 and l2:
            s1, s2 = find_code(types, "S", b1), find_code(types, "S", b2)
            f1, f2 = FREQ[sys][b1], FREQ[sys][b2]
            return dict(b1=b1, b2=b2, c1=c1, l1=l1, c2=c2, l2=l2, s1=s1, s2=s2, f1=f1, f2=f2)
    return None


def analyze_sat(hdr: Header, epochs: List[Epoch], sat: str, pair: dict) -> dict:
    """对一颗卫星做完整质量分析，返回指标与绘图序列。"""
    n = len(epochs)
    t0 = epochs[0].time
    t = np.full(n, np.nan)
    rho1 = np.full(n, np.nan)
    rho2 = np.full(n, np.nan)
    L1 = np.full(n, np.nan)
    L2 = np.full(n, np.nan)
    S1 = np.full(n, np.nan)
    S2 = np.full(n, np.nan)

    # GLONASS 用频道号修正频率
    f1, f2 = pair["f1"], pair["f2"]
    if sat[0] == "R":
        ch = hdr.glo_ch.get(int(sat[1:]), 0)
        f1 = glonass_freq(pair["b1"], ch)
        f2 = glonass_freq(pair["b2"], ch)

    for i, ep in enumerate(epochs):
        t[i] = (ep.time - t0).total_seconds()
        if sat not in ep.data:
            continue
        row = ep.data[sat]
        rho1[i] = row.get(pair["c1"], np.nan)
        rho2[i] = row.get(pair["c2"], np.nan)
        L1[i] = row.get(pair["l1"], np.nan)
        L2[i] = row.get(pair["l2"], np.nan)
        if pair["s1"]:
            S1[i] = row.get(pair["s1"], np.nan)
        if pair["s2"]:
            S2[i] = row.get(pair["s2"], np.nan)

    dual = np.isfinite(rho1) & np.isfinite(rho2) & np.isfinite(L1) & np.isfinite(L2)
    n_ep = int(np.sum(np.isfinite(rho1) | np.isfinite(L1)))
    n_ok = int(np.sum(dual))
    if n_ok < 5:
        return {"sat": sat, "n": n_ep, "n_ok": n_ok, "slips": 0}

    # 相位周 -> 米
    phi1 = L1 * wavelength(f1)
    phi2 = L2 * wavelength(f2)

    mw = mw_combo(phi1, phi2, rho1, rho2, f1, f2)
    gf = gf_combo(phi1, phi2)
    mw[~dual] = np.nan
    gf[~dual] = np.nan

    slips = detect_mw_slips(mw)
    slips = detect_gf_slips(gf, wavelength(f1), wavelength(f2), slips)

    mp1_raw, mp2_raw = multipath_combo(rho1, rho2, phi1, phi2, f1, f2)
    mp1_raw[~dual] = np.nan
    mp2_raw[~dual] = np.nan
    mp1 = detrend_mp(mp1_raw, slips)
    mp2 = detrend_mp(mp2_raw, slips)

    # 公式(19) IOD = ΔI / Δt
    I1 = iono_delay(phi1, phi2, f1, f2)
    I1[~dual] = np.nan
    iod = np.full(n, np.nan)
    for i in range(1, n):
        if slips[i] or slips[i - 1]:
            continue
        dt = t[i] - t[i - 1]
        if dt > 0 and np.isfinite(I1[i]) and np.isfinite(I1[i - 1]):
            iod[i] = (I1[i] - I1[i - 1]) / dt

    return {
        "sat": sat,
        "n": n_ep,
        "n_ok": n_ok,
        "slips": int(np.sum(slips)),
        "mp1_rms": rms(mp1),
        "mp2_rms": rms(mp2),
        "iod_rms": rms(iod),
        "iod_jumps": int(np.sum(np.abs(iod[np.isfinite(iod)]) > IOD_JUMP)),
        "pr_noise": {
            pair["c1"]: triple_diff_noise(rho1),
            pair["c2"]: triple_diff_noise(rho2),
        },
        "ph_noise": {
            pair["l1"]: triple_diff_noise(L1, slips),
            pair["l2"]: triple_diff_noise(L2, slips),
        },
        "cnr": {
            k: float(np.nanmean(v))
            for k, v in ((pair["s1"], S1), (pair["s2"], S2))
            if k and np.any(np.isfinite(v))
        },
        "t": t,
        "mp1": mp1,
        "iod": iod,
        "gf": gf,
        "slips_flag": slips,
    }


def analyze_all(hdr: Header, epochs: List[Epoch], systems: Optional[List[str]] = None) -> dict:
    """按系统汇总质量指标。"""
    systems = systems or list(hdr.obs_types.keys())
    # 理论历元数
    if hdr.t_first and hdr.t_last and hdr.interval > 0:
        n_theo = int(round((hdr.t_last - hdr.t_first).total_seconds() / hdr.interval)) + 1
    else:
        n_theo = len(epochs)

    result = {}
    for sys in systems:
        pair = select_pair(hdr, sys)
        if not pair:
            continue
        sats = sorted({s for ep in epochs for s in ep.data if s[0] == sys})
        sat_res = {}
        total_n = total_ok = total_slip = 0
        mp1s, mp2s, iods = [], [], []
        pr_acc, ph_acc, cnr_acc = defaultdict(list), defaultdict(list), defaultdict(list)

        for sat in sats:
            r = analyze_sat(hdr, epochs, sat, pair)
            if r.get("n", 0) == 0:
                continue
            sat_res[sat] = r
            total_n += r["n"]
            total_ok += r["n_ok"]
            total_slip += r["slips"]
            for key, bucket in (("mp1_rms", mp1s), ("mp2_rms", mp2s), ("iod_rms", iods)):
                if np.isfinite(r.get(key, np.nan)):
                    bucket.append(r[key])
            for k, v in r.get("pr_noise", {}).items():
                if np.isfinite(v):
                    pr_acc[k].append(v)
            for k, v in r.get("ph_noise", {}).items():
                if np.isfinite(v):
                    ph_acc[k].append(v)
            for k, v in r.get("cnr", {}).items():
                cnr_acc[k].append(v)

        n_sats = max(len(sat_res), 1)
        completeness = 100.0 * total_ok / (n_theo * n_sats)
        slip_ratio = (total_n / total_slip) if total_slip > 0 else float("inf")

        result[sys] = {
            "name": SYS_NAME.get(sys, sys),
            "pair": pair,
            "completeness": completeness,
            "slip_ratio": slip_ratio,
            "n_obs": total_n,
            "n_slips": total_slip,
            "mp1_rms": float(np.mean(mp1s)) if mp1s else float("nan"),
            "mp2_rms": float(np.mean(mp2s)) if mp2s else float("nan"),
            "iod_rms": float(np.mean(iods)) if iods else float("nan"),
            "pr_noise": {k: float(np.mean(v)) for k, v in pr_acc.items()},
            "ph_noise": {k: float(np.mean(v)) for k, v in ph_acc.items()},
            "cnr": {k: float(np.mean(v)) for k, v in cnr_acc.items()},
            "sats": sat_res,
        }
    return result


# =============================================================================
# 7. 报告与绘图
# =============================================================================

def format_report(hdr: Header, epochs: List[Epoch], qc: dict) -> str:
    lines = [
        "=" * 70,
        "GNSS 观测数据质量评估报告（BD 420022-2019）",
        "=" * 70,
        f"测站: {hdr.marker or '-'}    RINEX {hdr.version:.2f}",
        f"近似坐标 XYZ: {hdr.xyz}",
        f"采样间隔: {hdr.interval} s    历元数: {len(epochs)}",
        f"起止时间: {hdr.t_first}  ~  {hdr.t_last}",
        "",
    ]
    for sys, sq in qc.items():
        p = sq["pair"]
        sr = f"{sq['slip_ratio']:.1f}" if np.isfinite(sq["slip_ratio"]) else "Inf(无周跳)"
        lines += [
            "-" * 70,
            f"系统 {sys} ({sq['name']})  双频: {p['c1']}/{p['l1']} + {p['c2']}/{p['l2']}",
            f"  完整率 DIs      : {sq['completeness']:.2f} %",
            f"  周跳比 o/slps   : {sr}   (周跳数={sq['n_slips']})",
            f"  多路径 MP1 RMS  : {sq['mp1_rms']:.4f} m",
            f"  多路径 MP2 RMS  : {sq['mp2_rms']:.4f} m",
            f"  电离层 IOD RMS  : {sq['iod_rms']:.5f} m/s",
        ]
        if sq["pr_noise"]:
            lines.append("  伪距噪声        : " + ", ".join(f"{k}={v:.4f}m" for k, v in sq["pr_noise"].items()))
        if sq["ph_noise"]:
            lines.append("  相位噪声        : " + ", ".join(f"{k}={v:.5f}周" for k, v in sq["ph_noise"].items()))
        if sq["cnr"]:
            lines.append("  平均载噪比      : " + ", ".join(f"{k}={v:.1f}dBHz" for k, v in sq["cnr"].items()))
        lines.append(f"  {'SAT':<5} {'N':>6} {'Slip':>5} {'MP1':>8} {'MP2':>8} {'IOD':>9}")
        for sat, r in sorted(sq["sats"].items()):
            lines.append(
                f"  {sat:<5} {r['n']:6d} {r['slips']:5d} "
                f"{r.get('mp1_rms', float('nan')):8.3f} "
                f"{r.get('mp2_rms', float('nan')):8.3f} "
                f"{r.get('iod_rms', float('nan')):9.5f}"
            )
        lines.append("")
    lines.append("=" * 70)
    return "\n".join(lines)


def make_plots(qc: dict, outdir: Path) -> List[Path]:
    if not HAS_PLOT:
        return []
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []

    # 总览柱状图
    sys_list = list(qc.keys())
    names = [qc[s]["name"] for s in sys_list]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    axes[0, 0].bar(names, [qc[s]["completeness"] for s in sys_list], color="#2a6f97")
    axes[0, 0].set_title("数据完整率 (%)")
    axes[0, 1].bar(
        names,
        [qc[s]["slip_ratio"] if np.isfinite(qc[s]["slip_ratio"]) else qc[s]["n_obs"] for s in sys_list],
        color="#bc4749",
    )
    axes[0, 1].set_title("周跳比 o/slps（越大越好）")
    axes[1, 0].bar(names, [qc[s]["mp1_rms"] for s in sys_list], color="#6a994e")
    axes[1, 0].set_title("多路径 MP1 RMS (m)")
    axes[1, 1].bar(names, [qc[s]["iod_rms"] for s in sys_list], color="#e09f3e")
    axes[1, 1].set_title("电离层变化率 IOD RMS (m/s)")
    for ax in axes.ravel():
        ax.tick_params(axis="x", rotation=20)
        ax.grid(True, alpha=0.3)
    fig.suptitle("GNSS 质量评估总览 (BD 420022-2019)")
    fig.tight_layout()
    p = outdir / "qc_overview.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    paths.append(p)

    # 各系统时序（取观测最多的几颗星）
    for sys, sq in qc.items():
        ranked = sorted(sq["sats"].items(), key=lambda kv: kv[1].get("n_ok", 0), reverse=True)[:4]
        if not ranked:
            continue
        fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
        for sat, r in ranked:
            if "t" not in r:
                continue
            tm = r["t"] / 60.0
            axes[0].plot(tm, r["mp1"], lw=0.8, label=sat)
            axes[1].plot(tm, r["iod"], lw=0.8, label=sat)
            axes[2].plot(tm, r["gf"], lw=0.8, label=sat)
            fl = r["slips_flag"]
            if np.any(fl):
                axes[2].scatter(tm[fl], r["gf"][fl], c="red", s=18, marker="x", zorder=5)
        axes[0].set_ylabel("MP1 (m)")
        axes[0].legend(fontsize=8, ncol=4)
        axes[0].set_title(f"{sq['name']} 多路径 / IOD / GF")
        axes[1].axhline(IOD_JUMP, color="r", ls="--", lw=0.8)
        axes[1].axhline(-IOD_JUMP, color="r", ls="--", lw=0.8)
        axes[1].set_ylabel("IOD (m/s)")
        axes[2].set_ylabel("GF (m)")
        axes[2].set_xlabel("时间 (min)")
        for ax in axes:
            ax.grid(True, alpha=0.3)
        fig.tight_layout()
        p = outdir / f"qc_series_{sys}.png"
        fig.savefig(p, dpi=140)
        plt.close(fig)
        paths.append(p)
    return paths


# =============================================================================
# 8. 主程序
# =============================================================================

def main():
    ap = argparse.ArgumentParser(description="GNSS RINEX 质量分析（可读教学版）")
    ap.add_argument("rinex", nargs="?", default="data/sample_yygc1.26o", help=".26o 观测文件路径")
    ap.add_argument("--out", default="Result_QC", help="输出目录")
    ap.add_argument("--max-epochs", type=int, default=None, help="最多读取历元数（调试用）")
    ap.add_argument("--systems", default=None, help="例如 G,C,E")
    args = ap.parse_args()

    path = Path(args.rinex)
    if not path.exists():
        print(f"找不到文件: {path}")
        print("可先运行: python -m gnss_qc   （自动生成示例 data/sample_yygc1.26o）")
        return 1

    print(f"读取 {path} ...")
    hdr, epochs = read_rinex(path, max_epochs=args.max_epochs)
    print(f"  系统: {list(hdr.obs_types.keys())}, 历元: {len(epochs)}")

    systems = [s.strip().upper() for s in args.systems.split(",")] if args.systems else None
    print("计算质量指标 ...")
    qc = analyze_all(hdr, epochs, systems)

    report = format_report(hdr, epochs, qc)
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "qc_report.txt").write_text(report, encoding="utf-8")
    print(report)
    print(f"\n报告已保存: {outdir / 'qc_report.txt'}")

    for p in make_plots(qc, outdir):
        print(f"  图: {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
