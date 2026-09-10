"""GNSS frequencies and quality-assessment defaults (BD 420022-2019)."""

from __future__ import annotations

C_LIGHT = 299792458.0  # m/s

# Carrier frequencies (Hz)
FREQ = {
    "G": {  # GPS
        "1": 1575.42e6,   # L1
        "2": 1227.60e6,   # L2
        "5": 1176.45e6,   # L5
    },
    "R": {  # GLONASS nominal (channel k applied separately)
        "1": 1602.0e6,
        "2": 1246.0e6,
        "3": 1202.025e6,
    },
    "E": {  # Galileo
        "1": 1575.42e6,   # E1
        "5": 1176.45e6,   # E5a
        "7": 1207.14e6,   # E5b
        "6": 1278.75e6,   # E6
        "8": 1191.795e6,  # E5
    },
    "C": {  # BDS
        "2": 1561.098e6,  # B1I
        "1": 1575.42e6,   # B1C
        "7": 1207.14e6,   # B2I / B2b
        "5": 1176.45e6,   # B2a
        "6": 1268.52e6,   # B3I
    },
    "J": {  # QZSS (same as GPS L-band)
        "1": 1575.42e6,
        "2": 1227.60e6,
        "5": 1176.45e6,
        "6": 1278.75e6,
    },
}

SYS_NAME = {
    "G": "GPS",
    "R": "GLONASS",
    "E": "Galileo",
    "C": "BDS",
    "J": "QZSS",
    "I": "IRNSS",
    "S": "SBAS",
}

# Preferred dual-frequency pairs (band1, band2), largest Δf first (Appendix B)
DUAL_PAIRS = {
    "G": [("1", "2"), ("1", "5")],
    "R": [("1", "2"), ("1", "3")],
    "E": [("1", "5"), ("1", "7"), ("1", "6")],
    "C": [("2", "6"), ("2", "7"), ("1", "5")],
    "J": [("1", "2"), ("1", "5")],
}

# Thresholds from BD 420022-2019
MW_SIGMA_FACTOR = 4.0          # MW outlier / slip test
GF_LAMBDA_FACTOR = 6.0         # GF slip test multiplier
IOD_JUMP_THRES = 0.07          # m/s ionospheric jump
MP_WINDOW = 50                 # sliding window epochs
CLOCK_JUMP_XI = 4.0            # m, noise experience value
C_MS = C_LIGHT * 1e-3          # meters in 1 ms


def wavelength(freq_hz: float) -> float:
    return C_LIGHT / freq_hz


def glonass_freq(band: str, channel: int) -> float:
    """GLONASS FDMA frequency for slot channel k."""
    if band == "1":
        return (1602.0 + channel * 9.0 / 16.0) * 1e6
    if band == "2":
        return (1246.0 + channel * 7.0 / 16.0) * 1e6
    return FREQ["R"].get(band, FREQ["R"]["1"])


def band_from_obs(obs_code: str) -> str:
    """Extract frequency band digit from RINEX 3 observation code (e.g. L1C -> 1)."""
    if len(obs_code) < 2:
        return ""
    return obs_code[1]
