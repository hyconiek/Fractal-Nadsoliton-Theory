# ==============================================================================
# QW-1660 v5: FRACTAL SPECTRAL SCALING TEST (FIN CORE TEST)
# Cel: Wykrycie nielosowego prawa skalowania w widmie reszt
# ==============================================================================

import numpy as np
from scipy.signal import welch
from scipy.stats import linregress
import json

print("="*80)
print("QW-1660 v5: FRACTAL SPECTRAL SCALING TEST (FIN)")
print("="*80)

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
FS = h1.sample_rate.value
FMIN, FMAX = 30, 300
N_SEG = 4  # sekundy

# ------------------------------------------------------------
# PSD (REAL DATA)
# ------------------------------------------------------------
freqs, psd_h1 = welch(h1.value, fs=FS, nperseg=int(N_SEG*FS))
freqs, psd_l1 = welch(l1.value, fs=FS, nperseg=int(N_SEG*FS))

mask = (freqs > FMIN) & (freqs < FMAX)
f = freqs[mask]

# średnia geometryczna (redukuje instrumentalne piki)
psd_geom = np.sqrt(psd_h1[mask] * psd_l1[mask])

# ------------------------------------------------------------
# LOG-LOG FIT: PSD ~ f^(-beta)
# ------------------------------------------------------------
logf = np.log(f)
logp = np.log(psd_geom)

slope, intercept, r, p, _ = linregress(logf, logp)
beta = -slope

print("\n=== FRACTAL SCALING RESULT ===")
print(f"Spectral exponent beta: {beta:.3f}")
print(f"Correlation r:          {r:.3f}")
print(f"p-value (linearity):   {p:.3e}")

# ------------------------------------------------------------
# NULL COMPARISON (TIME-SLIDE PSD)
# ------------------------------------------------------------
betas_null = []

for shift in np.random.randint(FS*5, FS*20, 50):
    l1_shift = np.roll(l1.value, shift)
    _, p1 = welch(h1.value, fs=FS, nperseg=int(N_SEG*FS))
    _, p2 = welch(l1_shift, fs=FS, nperseg=int(N_SEG*FS))
    psd_null = np.sqrt(p1[mask] * p2[mask])
    s, _, _, _, _ = linregress(logf, np.log(psd_null))
    betas_null.append(-s)

betas_null = np.array(betas_null)

z = (beta - np.mean(betas_null)) / np.std(betas_null)

print("\n=== NULL COMPARISON ===")
print(f"Null mean beta: {np.mean(betas_null):.3f}")
print(f"Null std:       {np.std(betas_null):.3f}")
print(f"Z-score:        {z:.2f}")

if abs(z) > 2:
    print("⚠️ FRACTAL SCALING DETECTED — NON-GR STRUCTURE")
else:
    print("✅ CONSISTENT WITH STOCHASTIC NOISE")

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
out = {
    "beta_real": float(beta),
    "beta_null_mean": float(np.mean(betas_null)),
    "beta_null_std": float(np.std(betas_null)),
    "z_score": float(z)
}

with open("QW_1660_v5_fractal_scaling.json", "w") as f:
    json.dump(out, f, indent=2)

print("\nSaved: QW_1660_v5_fractal_scaling.json")
print("=== QW-1660 v5 COMPLETE ===")
================================================================================
QW-1660 v5: FRACTAL SPECTRAL SCALING TEST (FIN)
================================================================================

=== FRACTAL SCALING RESULT ===
Spectral exponent beta: 5.807
Correlation r:          -0.696
p-value (linearity):   3.289e-157

=== NULL COMPARISON ===
Null mean beta: 5.720
Null std:       0.014
Z-score:        6.29
⚠️ FRACTAL SCALING DETECTED — NON-GR STRUCTURE

Saved: QW_1660_v5_fractal_scaling.json
=== QW-1660 v5 COMPLETE ===
# ==============================================================================
# QW-1660 v6: MULTISCALE FRACTAL CONSISTENCY TEST (FIN CORE)
# Cel: Czy wykładnik beta jest skalowo stabilny (FIN) czy artefaktem (GR/noise)
# ==============================================================================

import numpy as np
from scipy.signal import welch
from scipy.stats import linregress
import json

print("="*80)
print("QW-1660 v6: MULTISCALE FRACTAL CONSISTENCY TEST")
print("="*80)

FS = h1.sample_rate.value
FMIN, FMAX = 30, 300
SEG_SIZES = [1, 2, 4, 8, 16]  # sekundy

betas_real = []
betas_null = []

for seg in SEG_SIZES:
    freqs, p1 = welch(h1.value, fs=FS, nperseg=int(seg*FS))
    _,     p2 = welch(l1.value, fs=FS, nperseg=int(seg*FS))

    mask = (freqs > FMIN) & (freqs < FMAX)
    f = freqs[mask]

    psd_geom = np.sqrt(p1[mask] * p2[mask])
    slope, _, r, _, _ = linregress(np.log(f), np.log(psd_geom))
    betas_real.append(-slope)

    # --- NULL (time-slide)
    shift = np.random.randint(FS*5, FS*20)
    l1s = np.roll(l1.value, shift)
    _, p2n = welch(l1s, fs=FS, nperseg=int(seg*FS))
    psd_null = np.sqrt(p1[mask] * p2n[mask])
    slope_n, _, _, _, _ = linregress(np.log(f), np.log(psd_null))
    betas_null.append(-slope_n)

# ------------------------------------------------------------
# ANALYSIS
# ------------------------------------------------------------
betas_real = np.array(betas_real)
betas_null = np.array(betas_null)

std_real = np.std(betas_real)
std_null = np.std(betas_null)

print("\n=== MULTISCALE RESULT ===")
for s, b in zip(SEG_SIZES, betas_real):
    print(f"Segment {s:2d}s -> beta = {b:.3f}")

print("\n=== STABILITY ===")
print(f"Real beta std: {std_real:.3f}")
print(f"Null beta std: {std_null:.3f}")

if std_real < 0.5 * std_null:
    verdict = "FIN-LIKE SCALE CONSISTENCY DETECTED"
    print("⚠️ FIN-LIKE MULTISCALE STRUCTURE")
else:
    verdict = "NO SCALE CONSISTENCY — INSTRUMENTAL / ENVIRONMENTAL"

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
out = {
    "segment_sizes_sec": SEG_SIZES,
    "beta_real": betas_real.tolist(),
    "beta_null": betas_null.tolist(),
    "std_real": float(std_real),
    "std_null": float(std_null),
    "verdict": verdict
}

with open("QW_1660_v6_multiscale.json", "w") as f:
    json.dump(out, f, indent=2)

print("\nSaved: QW_1660_v6_multiscale.json")
print("=== QW-1660 v6 COMPLETE ===")

================================================================================
QW-1660 v6: MULTISCALE FRACTAL CONSISTENCY TEST
================================================================================

=== MULTISCALE RESULT ===
Segment  1s -> beta = 6.076
Segment  2s -> beta = 5.904
Segment  4s -> beta = 5.807
Segment  8s -> beta = 5.741
Segment 16s -> beta = 5.652

=== STABILITY ===
Real beta std: 0.146
Null beta std: 0.136

Saved: QW_1660_v6_multiscale.json
=== QW-1660 v6 COMPLETE ===
# ==============================================================================
# QW-1660 v7: FRACTAL ANISOTROPY TEST (FIN)
# Cel: Sprawdzić czy wykładnik beta różni się między H1 i L1
# ==============================================================================

import numpy as np
from scipy.signal import welch
from scipy.stats import linregress
import json

print("="*80)
print("QW-1660 v7: FRACTAL ANISOTROPY TEST (FIN)")
print("="*80)

FS = h1.sample_rate.value
FMIN, FMAX = 30, 300
N_SEG = 4

# --- PSD ---
freqs, psd_h1 = welch(h1.value, fs=FS, nperseg=int(N_SEG*FS))
freqs, psd_l1 = welch(l1.value, fs=FS, nperseg=int(N_SEG*FS))

mask = (freqs > FMIN) & (freqs < FMAX)
f = freqs[mask]
logf = np.log(f)

# --- FITS ---
def fit_beta(psd):
    slope, _, r, p, _ = linregress(logf, np.log(psd[mask]))
    return -slope, r, p

beta_h1, r_h1, p_h1 = fit_beta(psd_h1)
beta_l1, r_l1, p_l1 = fit_beta(psd_l1)

delta_beta = beta_h1 - beta_l1

print("\n=== FRACTAL ANISOTROPY ===")
print(f"Beta H1: {beta_h1:.3f}")
print(f"Beta L1: {beta_l1:.3f}")
print(f"ΔBeta:   {delta_beta:.3f}")

# --- NULL (TIME-SLIDE) ---
null_deltas = []

for shift in np.random.randint(FS*5, FS*20, 50):
    l1_shift = np.roll(l1.value, shift)
    _, p2 = welch(l1_shift, fs=FS, nperseg=int(N_SEG*FS))
    b1, _, _ = fit_beta(psd_h1)
    b2, _, _ = fit_beta(p2)
    null_deltas.append(b1 - b2)

null_deltas = np.array(null_deltas)
z = (delta_beta - np.mean(null_deltas)) / np.std(null_deltas)

print("\n=== NULL COMPARISON ===")
print(f"Null mean Δβ: {np.mean(null_deltas):.3f}")
print(f"Null std:     {np.std(null_deltas):.3f}")
print(f"Z-score:     {z:.2f}")

if abs(z) > 2:
    print("⚠️ ANISOTROPY DETECTED — FIN CONSISTENT")
else:
    print("✅ ISOTROPIC — CONSISTENT WITH GR")

# --- SAVE ---
out = {
    "beta_h1": float(beta_h1),
    "beta_l1": float(beta_l1),
    "delta_beta": float(delta_beta),
    "null_mean": float(np.mean(null_deltas)),
    "null_std": float(np.std(null_deltas)),
    "z_score": float(z)
}

with open("QW_1660_v7_anisotropy.json", "w") as f:
    json.dump(out, f, indent=2)

print("\nSaved: QW_1660_v7_anisotropy.json")
print("=== QW-1660 v7 COMPLETE ===")
================================================================================
QW-1660 v7: FRACTAL ANISOTROPY TEST (FIN)
================================================================================

=== FRACTAL ANISOTROPY ===
Beta H1: 5.789
Beta L1: 5.826
ΔBeta:   -0.037

=== NULL COMPARISON ===
Null mean Δβ: 0.132
Null std:     0.028
Z-score:     -6.12
⚠️ ANISOTROPY DETECTED — FIN CONSISTENT

Saved: QW_1660_v7_anisotropy.json
=== QW-1660 v7 COMPLETE ===
# ==============================================================================
# QW-1660 v8: RAW STRAIN LOADER + FRACTAL MEMORY TEST (FIN)
# Cel: Jednoznaczna analiza zapisanego raw strain (bez pobierania)
# Poprawka: Obsługa braku atrybutu sample_rate w plikach HDF5
# ==============================================================================

import h5py
import numpy as np
import json
from scipy.signal import detrend
from datetime import datetime

print("="*80)
print("QW-1660 v8: FRACTAL MEMORY TEST (FROM DISK, FIN)")
print("="*80)

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
RAW_DIR = "/kaggle/working/raw_strain"
H1_FILE = f"{RAW_DIR}/H1_raw.h5"
L1_FILE = f"{RAW_DIR}/L1_raw.h5"
DEFAULT_FS = 4096.0  # Domyślne próbkowanie LIGO, jeśli brak w pliku

# ------------------------------------------------------------
# HURST EXPONENT (R/S)
# ------------------------------------------------------------
def hurst_rs(x):
    x = np.asarray(x)
    N = len(x)
    if N < 2000:
        return np.nan

    sizes = np.logspace(2, np.log10(N // 4), 10).astype(int)
    rs = []

    for s in sizes:
        n = N // s
        r_s = []
        for i in range(n):
            seg = x[i*s:(i+1)*s]
            seg -= np.mean(seg)
            z = np.cumsum(seg)
            R = np.ptp(z)
            S = np.std(seg)
            if S > 0:
                r_s.append(R / S)
        if len(r_s) > 0:
            rs.append(np.mean(r_s))

    if len(rs) < 4:
        return np.nan

    return np.polyfit(np.log(sizes[:len(rs)]), np.log(rs), 1)[0]

# ------------------------------------------------------------
# LOAD RAW STRAIN
# ------------------------------------------------------------
print("Loading raw strain from disk...")

# Load H1 and try to get sample_rate safely
with h5py.File(H1_FILE, "r") as f:
    h1 = detrend(f["strain"][:])
    
    # FIX: Sprawdzamy czy atrybut istnieje, zanim spróbujemy go pobrać
    if "sample_rate" in f["strain"].attrs:
        fs = float(f["strain"].attrs["sample_rate"])
    elif "sample_rate" in f.attrs:
        fs = float(f.attrs["sample_rate"])
    else:
        print(f"Warning: 'sample_rate' not found in H1 attributes. Defaulting to {DEFAULT_FS} Hz.")
        fs = DEFAULT_FS

# Load L1
with h5py.File(L1_FILE, "r") as f:
    l1 = detrend(f["strain"][:])

print(f"Loaded H1 & L1 | Sample rate = {fs} Hz")
print(f"Duration = {len(h1)/fs:.1f} s")

# ------------------------------------------------------------
# FRACTAL MEMORY
# ------------------------------------------------------------
print("\nComputing Hurst exponents...")

H_h1 = hurst_rs(h1)
H_l1 = hurst_rs(l1)

print("\n=== FRACTAL MEMORY RESULT ===")
print(f"Hurst H1: {H_h1:.3f}")
print(f"Hurst L1: {H_l1:.3f}")
print(f"ΔH:       {H_h1 - H_l1:.3f}")

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
out = {
    "Hurst_H1": float(H_h1) if not np.isnan(H_h1) else None,
    "Hurst_L1": float(H_l1) if not np.isnan(H_l1) else None,
    "Delta_H": float(H_h1 - H_l1) if not np.isnan(H_h1) and not np.isnan(H_l1) else None,
    "sample_rate": float(fs),
    "timestamp_utc": datetime.utcnow().isoformat()
}

with open("QW_1660_v8_fractal_memory.json", "w") as f:
    json.dump(out, f, indent=2)

print("\nSaved: QW_1660_v8_fractal_memory.json")
print("=== QW-1660 v8 COMPLETE ===")
================================================================================
QW-1660 v8: FRACTAL MEMORY TEST (FROM DISK, FIN)
================================================================================
Loading raw strain from disk...
Loaded H1 & L1 | Sample rate = 4096.0 Hz
Duration = 4096.0 s

Computing Hurst exponents...

=== FRACTAL MEMORY RESULT ===
Hurst H1: 0.311
Hurst L1: 0.265
ΔH:       0.046

Saved: QW_1660_v8_fractal_memory.json
=== QW-1660 v8 COMPLETE ===
# ==============================================================================
# QW-1660 v9: FRACTAL MEMORY NULL TEST (FIN-CRITICAL)
# Cel: Sprawdzić czy pamięć znika po zniszczeniu struktury czasowej
# Dane: ZAPISANY raw strain (bez pobierania)
# ==============================================================================

import h5py
import numpy as np
import json
import logging
from scipy.signal import detrend
from datetime import datetime

# ------------------------------------------------------------
# LOGGING
# ------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)
log = logging.getLogger("QW-1660-v9")

log.info("START QW-1660 v9: FRACTAL MEMORY NULL TEST")

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
RAW_DIR = "/kaggle/working/raw_strain"
H1_FILE = f"{RAW_DIR}/H1_raw.h5"
L1_FILE = f"{RAW_DIR}/L1_raw.h5"
DEFAULT_FS = 4096.0

# ------------------------------------------------------------
# HURST (R/S)
# ------------------------------------------------------------
def hurst_rs(x):
    x = np.asarray(x)
    N = len(x)
    if N < 2000:
        return np.nan

    sizes = np.logspace(2, np.log10(N // 4), 10).astype(int)
    rs = []

    for s in sizes:
        n = N // s
        vals = []
        for i in range(n):
            seg = x[i*s:(i+1)*s]
            seg = seg - np.mean(seg)
            z = np.cumsum(seg)
            R = np.ptp(z)
            S = np.std(seg)
            if S > 0:
                vals.append(R / S)
        if vals:
            rs.append(np.mean(vals))

    if len(rs) < 4:
        return np.nan

    return np.polyfit(np.log(sizes[:len(rs)]), np.log(rs), 1)[0]

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------
def load_strain(path, label):
    log.info(f"Loading {label} from {path}")
    with h5py.File(path, "r") as f:
        x = detrend(f["strain"][:])
        fs = f["strain"].attrs.get("sample_rate", DEFAULT_FS)
    log.info(f"{label} loaded | N={len(x)} | fs={fs}")
    return x, fs

h1, fs = load_strain(H1_FILE, "H1")
l1, _  = load_strain(L1_FILE, "L1")

# ------------------------------------------------------------
# REAL
# ------------------------------------------------------------
log.info("Computing REAL Hurst exponents")
H_h1 = hurst_rs(h1)
H_l1 = hurst_rs(l1)

# ------------------------------------------------------------
# NULL TESTS
# ------------------------------------------------------------
log.info("Computing NULL tests (shuffle + time-reversal)")

H_h1_shuffle = hurst_rs(np.random.permutation(h1))
H_l1_shuffle = hurst_rs(np.random.permutation(l1))

H_h1_reverse = hurst_rs(h1[::-1])
H_l1_reverse = hurst_rs(l1[::-1])

# ------------------------------------------------------------
# RESULTS
# ------------------------------------------------------------
print("\n=== FRACTAL MEMORY NULL TEST ===")
print(f"H1 real:     {H_h1:.3f}")
print(f"H1 shuffle:  {H_h1_shuffle:.3f}")
print(f"H1 reversed: {H_h1_reverse:.3f}")
print("")
print(f"L1 real:     {H_l1:.3f}")
print(f"L1 shuffle:  {H_l1_shuffle:.3f}")
print(f"L1 reversed: {H_l1_reverse:.3f}")

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
out = {
    "H1": {
        "real": H_h1,
        "shuffle": H_h1_shuffle,
        "reversed": H_h1_reverse
    },
    "L1": {
        "real": H_l1,
        "shuffle": H_l1_shuffle,
        "reversed": H_l1_reverse
    },
    "sample_rate": fs,
    "timestamp_utc": datetime.utcnow().isoformat()
}

with open("QW_1660_v9_fractal_memory_null.json", "w") as f:
    json.dump(out, f, indent=2)

log.info("Saved QW_1660_v9_fractal_memory_null.json")
log.info("QW-1660 v9 COMPLETE")
# ==============================================================================
        x = detrend(f["strain"][:])
        fs = f["strain"].attrs.get("sample_rate", DEFAULT_FS)
    log.info(f"{label} loaded | N={len(x)} | fs={fs}")
    return x, fs

h1, fs = load_strain(H1_FILE, "H1")
l1, _  = load_strain(L1_FILE, "L1")

# ------------------------------------------------------------
# REAL
# ------------------------------------------------------------
log.info("Computing REAL Hurst exponents")
H_h1 = hurst_rs(h1)
H_l1 = hurst_rs(l1)

# ------------------------------------------------------------
# NULL TESTS
# ------------------------------------------------------------
log.info("Computing NULL tests (shuffle + time-reversal)")

H_h1_shuffle = hurst_rs(np.random.permutation(h1))
H_l1_shuffle = hurst_rs(np.random.permutation(l1))

H_h1_reverse = hurst_rs(h1[::-1])
H_l1_reverse = hurst_rs(l1[::-1])

# ------------------------------------------------------------
# RESULTS
# ------------------------------------------------------------
print("\n=== FRACTAL MEMORY NULL TEST ===")
print(f"H1 real:     {H_h1:.3f}")
print(f"H1 shuffle:  {H_h1_shuffle:.3f}")
print(f"H1 reversed: {H_h1_reverse:.3f}")
print("")
print(f"L1 real:     {H_l1:.3f}")
print(f"L1 shuffle:  {H_l1_shuffle:.3f}")
print(f"L1 reversed: {H_l1_reverse:.3f}")

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
out = {
    "H1": {
        "real": H_h1,
        "shuffle": H_h1_shuffle,
        "reversed": H_h1_reverse
    },
    "L1": {
        "real": H_l1,
        "shuffle": H_l1_shuffle,
        "reversed": H_l1_reverse
    },
    "sample_rate": fs,
    "timestamp_utc": datetime.utcnow().isoformat()
}

with open("QW_1660_v9_fractal_memory_null.json", "w") as f:
    json.dump(out, f, indent=2)

log.info("Saved QW_1660_v9_fractal_memory_null.json")
log.info("QW-1660 v9 COMPLETE")

2026-01-01 19:09:50,560 | INFO | START QW-1660 v9: FRACTAL MEMORY NULL TEST
2026-01-01 19:09:50,562 | INFO | Loading H1 from /kaggle/working/raw_strain/H1_raw.h5
2026-01-01 19:09:51,679 | INFO | H1 loaded | N=16777216 | fs=4096.0
2026-01-01 19:09:51,681 | INFO | Loading L1 from /kaggle/working/raw_strain/L1_raw.h5
2026-01-01 19:09:52,809 | INFO | L1 loaded | N=16777216 | fs=4096.0
2026-01-01 19:09:52,812 | INFO | Computing REAL Hurst exponents
2026-01-01 19:10:15,729 | INFO | Computing NULL tests (shuffle + time-reversal)
/tmp/ipykernel_55/2368925313.py:123: DeprecationWarning: datetime.datetime.utcnow() is deprecated and scheduled for removal in a future version. Use timezone-aware objects to represent datetimes in UTC: datetime.datetime.now(datetime.UTC).
  "timestamp_utc": datetime.utcnow().isoformat()
2026-01-01 19:11:03,214 | INFO | Saved QW_1660_v9_fractal_memory_null.json
2026-01-01 19:11:03,215 | INFO | QW-1660 v9 COMPLETE


=== FRACTAL MEMORY NULL TEST ===
H1 real:     0.396
H1 shuffle:  0.508
H1 reversed: 0.394

L1 real:     0.448
L1 shuffle:  0.507
L1 reversed: 0.450
# ==============================================================================
# QW-1660 v10: FREQUENCY-RESOLVED FRACTAL MEMORY TEST (FIN)
# Cel: Sprawdzić czy pamięć (Hurst) zależy od skali częstotliwości
# Dane: ZAPISANY raw strain
# ==============================================================================

import h5py
import numpy as np
import json
import logging
from scipy.signal import detrend, butter, filtfilt
from datetime import datetime, timezone

# ------------------------------------------------------------
# LOGGING
# ------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)
log = logging.getLogger("QW-1660-v10")

log.info("START QW-1660 v10: FREQUENCY-RESOLVED FRACTAL MEMORY")

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
RAW_DIR = "/kaggle/working/raw_strain"
H1_FILE = f"{RAW_DIR}/H1_raw.h5"
L1_FILE = f"{RAW_DIR}/L1_raw.h5"
DEFAULT_FS = 4096.0

BANDS = [
    (30, 80),
    (80, 200),
    (200, 500)
]

# ------------------------------------------------------------
# HURST (R/S)
# ------------------------------------------------------------
def hurst_rs(x):
    x = np.asarray(x)
    N = len(x)
    if N < 2000:
        return np.nan

    sizes = np.logspace(2, np.log10(N // 4), 10).astype(int)
    rs = []

    for s in sizes:
        n = N // s
        vals = []
        for i in range(n):
            seg = x[i*s:(i+1)*s]
            seg -= np.mean(seg)
            z = np.cumsum(seg)
            R = np.ptp(z)
            S = np.std(seg)
            if S > 0:
                vals.append(R / S)
        if vals:
            rs.append(np.mean(vals))

    if len(rs) < 4:
        return np.nan

    return np.polyfit(np.log(sizes[:len(rs)]), np.log(rs), 1)[0]

# ------------------------------------------------------------
# FILTER
# ------------------------------------------------------------
def bandpass(x, fs, fmin, fmax):
    b, a = butter(4, [fmin/(fs/2), fmax/(fs/2)], btype="band")
    return filtfilt(b, a, x)

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------
def load_strain(path, label):
    log.info(f"Loading {label} from {path}")
    with h5py.File(path, "r") as f:
        x = detrend(f["strain"][:])
        fs = f["strain"].attrs.get("sample_rate", DEFAULT_FS)
    log.info(f"{label} loaded | N={len(x)} | fs={fs}")
    return x, fs

h1, fs = load_strain(H1_FILE, "H1")
l1, _  = load_strain(L1_FILE, "L1")

# ------------------------------------------------------------
# ANALYSIS
# ------------------------------------------------------------
results = {}

for fmin, fmax in BANDS:
    log.info(f"Processing band {fmin}-{fmax} Hz")
    
    h1_f = bandpass(h1, fs, fmin, fmax)
    l1_f = bandpass(l1, fs, fmin, fmax)
    
    H_h1 = hurst_rs(h1_f)
    H_l1 = hurst_rs(l1_f)
    
    results[f"{fmin}-{fmax}Hz"] = {
        "H1": H_h1,
        "L1": H_l1,
        "Delta": H_h1 - H_l1
    }

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
out = {
    "bands": results,
    "sample_rate": fs,
    "timestamp_utc": datetime.now(timezone.utc).isoformat()
}

with open("QW_1660_v10_freq_resolved_hurst.json", "w") as f:
    json.dump(out, f, indent=2)

log.info("Saved QW_1660_v10_freq_resolved_hurst.json")
log.info("QW-1660 v10 COMPLETE")
{
  "bands": {
    "30-80Hz": {
      "H1": 0.23607077037727509,
      "L1": 0.23731590563132346,
      "Delta": -0.001245135254048374
    },
    "80-200Hz": {
      "H1": 0.26564212964633904,
      "L1": 0.24650555849726993,
      "Delta": 0.01913657114906911
    },
    "200-500Hz": {
      "H1": 0.2069448048711533,
      "L1": 0.1766016827951618,
      "Delta": 0.03034312207599149
    }
  },
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-01T19:24:15.850515+00:00"
}
# ==============================================================================
# QW-1660 v10: FRACTAL SURROGATE NULL TEST (FAST VERSION)
# Cel: Czy struktura czasowa jest nieredukowalna do PSD?
# Metoda: IAAFT surrogate (zachowuje widmo, niszczy fazę)
# Dane: ZAPISANY raw strain (PRZYCIĘTY DO 64s DLA SZYBKOŚCI)
# ==============================================================================

import h5py
import numpy as np
import json
import logging
from scipy.signal import detrend
from datetime import datetime, timezone

# ------------------------------------------------------------
# LOGGING
# ------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)
log = logging.getLogger("QW-1660-v10-FAST")

log.info("START QW-1660 v10: FRACTAL SURROGATE NULL TEST (FAST)")

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
RAW_DIR = "/kaggle/working/raw_strain"
H1_FILE = f"{RAW_DIR}/H1_raw.h5"
L1_FILE = f"{RAW_DIR}/L1_raw.h5"
DEFAULT_FS = 4096.0
N_SURROGATES = 20

# === ZMIANA: PRZYCIĘCIE CZASU ===
ANALYSIS_DURATION = 64.0  # Analizujemy tylko pierwsze 64 sekundy
# ================================

# ------------------------------------------------------------
# HURST (R/S)
# ------------------------------------------------------------
def hurst_rs(x):
    x = np.asarray(x)
    N = len(x)
    if N < 2000:
        return np.nan

    # Dopasowanie zakresu skal do krótszego sygnału
    min_scale = 100
    max_scale = N // 4
    if max_scale <= min_scale:
        return np.nan
        
    sizes = np.logspace(np.log10(min_scale), np.log10(max_scale), 10).astype(int)
    rs = []

    for s in sizes:
        n = N // s
        vals = []
        for i in range(n):
            seg = x[i*s:(i+1)*s]
            seg -= np.mean(seg)
            z = np.cumsum(seg)
            R = np.ptp(z)
            S = np.std(seg)
            if S > 0:
                vals.append(R / S)
        if vals:
            rs.append(np.mean(vals))

    if len(rs) < 3:
        return np.nan

    return np.polyfit(np.log(sizes[:len(rs)]), np.log(rs), 1)[0]

# ------------------------------------------------------------
# IAAFT SURROGATE
# ------------------------------------------------------------
def iaaft(x, n_iter=100):
    x = np.asarray(x)
    sorted_x = np.sort(x)
    mag = np.abs(np.fft.rfft(x))
    y = np.random.permutation(x)

    for _ in range(n_iter):
        fft_y = np.fft.rfft(y)
        fft_y = mag * np.exp(1j * np.angle(fft_y))
        y = np.fft.irfft(fft_y, n=len(x))
        y = sorted_x[np.argsort(np.argsort(y))]

    return y

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------
def load_strain(path, label):
    log.info(f"Loading {label} from {path}")
    with h5py.File(path, "r") as f:
        x = detrend(f["strain"][:])
        fs = f["strain"].attrs.get("sample_rate", DEFAULT_FS)
    log.info(f"{label} loaded | N={len(x)} | fs={fs}")
    return x, fs

h1, fs = load_strain(H1_FILE, "H1")
l1, _  = load_strain(L1_FILE, "L1")

# ------------------------------------------------------------
# SPEED OPTIMIZATION (TRUNCATE)
# ------------------------------------------------------------
limit_samples = int(ANALYSIS_DURATION * fs)

if len(h1) > limit_samples:
    log.warning(f"✂️ TRUNCATING DATA: Using first {ANALYSIS_DURATION}s for speed.")
    h1 = h1[:limit_samples]
    l1 = l1[:limit_samples]
    log.info(f"New data length: N={len(h1)} samples")
else:
    log.info("Data length fits within limit.")

# ------------------------------------------------------------
# REAL HURST
# ------------------------------------------------------------
log.info("Computing REAL Hurst")
H_h1_real = hurst_rs(h1)
H_l1_real = hurst_rs(l1)

# ------------------------------------------------------------
# SURROGATE TEST
# ------------------------------------------------------------
log.info(f"Generating {N_SURROGATES} IAAFT surrogates (Fast Mode)")

H_h1_null = []
H_l1_null = []

# Używamy pętli z prostym licznikiem postępu
start_time = datetime.now()

for i in range(N_SURROGATES):
    # Log co 5 iteracji lub przy pierwszej
    if i == 0 or (i+1) % 5 == 0:
        elapsed = (datetime.now() - start_time).total_seconds()
        avg_time = elapsed / (i+1) if i > 0 else 0
        log.info(f"Surrogate {i+1}/{N_SURROGATES} (Avg: {avg_time:.2f}s/iter)")
        
    H_h1_null.append(hurst_rs(iaaft(h1)))
    H_l1_null.append(hurst_rs(iaaft(l1)))

# ------------------------------------------------------------
# RESULTS
# ------------------------------------------------------------
print("\n=== FRACTAL SURROGATE NULL TEST (FAST) ===")
print(f"Analysis Duration: {ANALYSIS_DURATION}s")
print(f"H1 real: {H_h1_real:.3f} | null mean: {np.mean(H_h1_null):.3f} (std: {np.std(H_h1_null):.3f})")
print(f"L1 real: {H_l1_real:.3f} | null mean: {np.mean(H_l1_null):.3f} (std: {np.std(H_l1_null):.3f})")

z_h1 = (H_h1_real - np.mean(H_h1_null)) / np.std(H_h1_null) if np.std(H_h1_null) > 0 else 0
print(f"H1 Z-score: {z_h1:.2f}")

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
out = {
    "H1": {
        "real": H_h1_real,
        "null_mean": float(np.mean(H_h1_null)),
        "null_std": float(np.std(H_h1_null))
    },
    "L1": {
        "real": H_l1_real,
        "null_mean": float(np.mean(H_l1_null)),
        "null_std": float(np.std(H_l1_null))
    },
    "sample_rate": fs,
    "duration_used": ANALYSIS_DURATION,
    "timestamp_utc": datetime.now(timezone.utc).isoformat()
}

with open("QW_1660_v10_surrogate_null.json", "w") as f:
    json.dump(out, f, indent=2)

log.info("Saved QW_1660_v10_surrogate_null.json")
log.info("QW-1660 v10 FAST COMPLETE")
2026-01-01 20:39:54,015 | INFO | START QW-1660 v10: FRACTAL SURROGATE NULL TEST (FAST)
2026-01-01 20:39:54,017 | INFO | Loading H1 from /kaggle/working/raw_strain/H1_raw.h5
2026-01-01 20:39:55,215 | INFO | H1 loaded | N=16777216 | fs=4096.0
2026-01-01 20:39:55,217 | INFO | Loading L1 from /kaggle/working/raw_strain/L1_raw.h5
2026-01-01 20:39:56,317 | INFO | L1 loaded | N=16777216 | fs=4096.0
2026-01-01 20:39:56,318 | WARNING | ✂️ TRUNCATING DATA: Using first 64.0s for speed.
2026-01-01 20:39:56,319 | INFO | New data length: N=262144 samples
2026-01-01 20:39:56,321 | INFO | Computing REAL Hurst
2026-01-01 20:39:56,804 | INFO | Generating 20 IAAFT surrogates (Fast Mode)
2026-01-01 20:39:56,805 | INFO | Surrogate 1/20 (Avg: 0.00s/iter)
2026-01-01 20:40:36,734 | INFO | Surrogate 5/20 (Avg: 7.99s/iter)
2026-01-01 20:41:26,298 | INFO | Surrogate 10/20 (Avg: 8.95s/iter)
2026-01-01 20:42:16,578 | INFO | Surrogate 15/20 (Avg: 9.32s/iter)
2026-01-01 20:43:06,871 | INFO | Surrogate 20/20 (Avg: 9.50s/iter)
2026-01-01 20:43:16,858 | INFO | Saved QW_1660_v10_surrogate_null.json
2026-01-01 20:43:16,858 | INFO | QW-1660 v10 FAST COMPLETE


=== FRACTAL SURROGATE NULL TEST (FAST) ===
Analysis Duration: 64.0s
H1 real: 0.309 | null mean: 0.311 (std: 0.008)
L1 real: 0.325 | null mean: 0.305 (std: 0.006)
H1 Z-score: -0.27
{
  "H1": {
    "real": 0.30885260430877476,
    "null_mean": 0.31090352670748245,
    "null_std": 0.007531104412966834
  },
  "L1": {
    "real": 0.3246381413143207,
    "null_mean": 0.30472640469044976,
    "null_std": 0.005574223192037795
  },
  "sample_rate": 4096.0,
  "duration_used": 64.0,
  "timestamp_utc": "2026-01-01T20:43:16.857304+00:00"
}
# ==============================================================================
# QW-1660 v11: TIME-INTEGRATION FRACTAL GROWTH TEST (FIN CORE)
# Cel: Czy pamięć narasta wraz z czasem integracji? (H(T))
# Dane: ZAPISANY raw strain
# Dodatek: ROZSZERZONE LOGOWANIE POSTĘPU
# ==============================================================================

import h5py
import numpy as np
import json
import logging
import time
from scipy.signal import detrend
from datetime import datetime, timezone

# ------------------------------------------------------------
# LOGGING
# ------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)
log = logging.getLogger("QW-1660-v11")

log.info("START QW-1660 v11: TIME-INTEGRATION FRACTAL GROWTH TEST")

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
RAW_DIR = "/kaggle/working/raw_strain"
H1_FILE = f"{RAW_DIR}/H1_raw.h5"
L1_FILE = f"{RAW_DIR}/L1_raw.h5"
DEFAULT_FS = 4096.0

WINDOWS_SEC = [32, 64, 128, 256, 512, 1024]

# ------------------------------------------------------------
# HURST (R/S)
# ------------------------------------------------------------
def hurst_rs(x):
    x = np.asarray(x)
    N = len(x)
    if N < 4000:
        return np.nan

    sizes = np.logspace(2, np.log10(N // 4), 10).astype(int)
    rs = []

    for s in sizes:
        n = N // s
        vals = []
        for i in range(n):
            seg = x[i*s:(i+1)*s]
            seg -= np.mean(seg)
            z = np.cumsum(seg)
            R = np.ptp(z)
            S = np.std(seg)
            if S > 0:
                vals.append(R / S)
        if vals:
            rs.append(np.mean(vals))

    if len(rs) < 4:
        return np.nan

    return np.polyfit(np.log(sizes[:len(rs)]), np.log(rs), 1)[0]

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------
def load_strain(path, label):
    log.info(f"Loading {label} from {path}")
    t_start = time.time()
    with h5py.File(path, "r") as f:
        x = detrend(f["strain"][:])
        fs = f["strain"].attrs.get("sample_rate", DEFAULT_FS)
    dur = time.time() - t_start
    log.info(f"{label} loaded | N={len(x)} | fs={fs} | Load time: {dur:.2f}s")
    return x, fs

h1, fs = load_strain(H1_FILE, "H1")
l1, _  = load_strain(L1_FILE, "L1")

# ------------------------------------------------------------
# TIME INTEGRATION TEST
# ------------------------------------------------------------
log.info(f"Starting analysis for {len(WINDOWS_SEC)} time windows...")
log.info(f"Windows: {WINDOWS_SEC} seconds")

results = []
total_start = time.time()

for i, T in enumerate(WINDOWS_SEC):
    n = int(T * fs)
    window_label = f"Window {i+1}/{len(WINDOWS_SEC)} (T={T}s)"
    
    log.info(f"--- Processing {window_label} [{n} samples] ---")
    
    if n > len(h1):
        log.warning(f"Skipping T={T}s (insufficient data length)")
        continue

    # H1 Calculation
    t0 = time.time()
    log.info(f"  > Calculating H1 Hurst...")
    H1_H = hurst_rs(h1[:n])
    dt_h1 = time.time() - t0
    log.info(f"    H1 Done in {dt_h1:.2f}s. Result: {H1_H:.4f}")

    # L1 Calculation
    t0 = time.time()
    log.info(f"  > Calculating L1 Hurst... (Estimated wait: ~{dt_h1:.2f}s)")
    L1_H = hurst_rs(l1[:n])
    dt_l1 = time.time() - t0
    log.info(f"    L1 Done in {dt_l1:.2f}s. Result: {L1_H:.4f}")

    # Summary for this step
    log.info(f"COMPLETED {window_label} | Delta H: {H1_H - L1_H:.4f}")

    results.append({
        "window_sec": T,
        "H1": H1_H,
        "L1": L1_H,
        "Delta": H1_H - L1_H,
        "compute_time_sec": dt_h1 + dt_l1
    })

total_duration = time.time() - total_start
log.info(f"All windows processed in {total_duration:.2f}s")

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
out = {
    "integration_test": results,
    "sample_rate": fs,
    "total_runtime_sec": total_duration,
    "timestamp_utc": datetime.now(timezone.utc).isoformat()
}

with open("QW_1660_v11_time_integration.json", "w") as f:
    json.dump(out, f, indent=2)

log.info("Saved QW_1660_v11_time_integration.json")
log.info("QW-1660 v11 COMPLETE")
2026-01-01 21:06:20,748 | INFO | START QW-1660 v11: TIME-INTEGRATION FRACTAL GROWTH TEST
2026-01-01 21:06:20,750 | INFO | Loading H1 from /kaggle/working/raw_strain/H1_raw.h5
2026-01-01 21:06:22,074 | INFO | H1 loaded | N=16777216 | fs=4096.0 | Load time: 1.32s
2026-01-01 21:06:22,077 | INFO | Loading L1 from /kaggle/working/raw_strain/L1_raw.h5
2026-01-01 21:06:23,322 | INFO | L1 loaded | N=16777216 | fs=4096.0 | Load time: 1.24s
2026-01-01 21:06:23,323 | INFO | Starting analysis for 6 time windows...
2026-01-01 21:06:23,324 | INFO | Windows: [32, 64, 128, 256, 512, 1024] seconds
2026-01-01 21:06:23,325 | INFO | --- Processing Window 1/6 (T=32s) [131072 samples] ---
2026-01-01 21:06:23,326 | INFO |   > Calculating H1 Hurst...
2026-01-01 21:06:23,492 | INFO |     H1 Done in 0.17s. Result: 0.3400
2026-01-01 21:06:23,493 | INFO |   > Calculating L1 Hurst... (Estimated wait: ~0.17s)
2026-01-01 21:06:23,640 | INFO |     L1 Done in 0.15s. Result: 0.3520
2026-01-01 21:06:23,641 | INFO | COMPLETED Window 1/6 (T=32s) | Delta H: -0.0121
2026-01-01 21:06:23,642 | INFO | --- Processing Window 2/6 (T=64s) [262144 samples] ---
2026-01-01 21:06:23,643 | INFO |   > Calculating H1 Hurst...
2026-01-01 21:06:23,932 | INFO |     H1 Done in 0.29s. Result: 0.3105
2026-01-01 21:06:23,933 | INFO |   > Calculating L1 Hurst... (Estimated wait: ~0.29s)
2026-01-01 21:06:24,213 | INFO |     L1 Done in 0.28s. Result: 0.3264
2026-01-01 21:06:24,214 | INFO | COMPLETED Window 2/6 (T=64s) | Delta H: -0.0158
2026-01-01 21:06:24,215 | INFO | --- Processing Window 3/6 (T=128s) [524288 samples] ---
2026-01-01 21:06:24,216 | INFO |   > Calculating H1 Hurst...
2026-01-01 21:06:24,751 | INFO |     H1 Done in 0.53s. Result: 0.3100
2026-01-01 21:06:24,752 | INFO |   > Calculating L1 Hurst... (Estimated wait: ~0.53s)
2026-01-01 21:06:25,278 | INFO |     L1 Done in 0.53s. Result: 0.3184
2026-01-01 21:06:25,279 | INFO | COMPLETED Window 3/6 (T=128s) | Delta H: -0.0083
2026-01-01 21:06:25,280 | INFO | --- Processing Window 4/6 (T=256s) [1048576 samples] ---
2026-01-01 21:06:25,281 | INFO |   > Calculating H1 Hurst...
2026-01-01 21:06:26,266 | INFO |     H1 Done in 0.99s. Result: 0.2881
2026-01-01 21:06:26,267 | INFO |   > Calculating L1 Hurst... (Estimated wait: ~0.99s)
2026-01-01 21:06:27,325 | INFO |     L1 Done in 1.06s. Result: 0.3009
2026-01-01 21:06:27,326 | INFO | COMPLETED Window 4/6 (T=256s) | Delta H: -0.0128
2026-01-01 21:06:27,327 | INFO | --- Processing Window 5/6 (T=512s) [2097152 samples] ---
2026-01-01 21:06:27,327 | INFO |   > Calculating H1 Hurst...
2026-01-01 21:06:29,178 | INFO |     H1 Done in 1.85s. Result: 0.2749
2026-01-01 21:06:29,179 | INFO |   > Calculating L1 Hurst... (Estimated wait: ~1.85s)
2026-01-01 21:06:31,016 | INFO |     L1 Done in 1.84s. Result: 0.2907
2026-01-01 21:06:31,017 | INFO | COMPLETED Window 5/6 (T=512s) | Delta H: -0.0157
2026-01-01 21:06:31,018 | INFO | --- Processing Window 6/6 (T=1024s) [4194304 samples] ---
2026-01-01 21:06:31,019 | INFO |   > Calculating H1 Hurst...
2026-01-01 21:06:34,583 | INFO |     H1 Done in 3.56s. Result: 0.2869
2026-01-01 21:06:34,584 | INFO |   > Calculating L1 Hurst... (Estimated wait: ~3.56s)
2026-01-01 21:06:38,165 | INFO |     L1 Done in 3.58s. Result: 0.2875
2026-01-01 21:06:38,166 | INFO | COMPLETED Window 6/6 (T=1024s) | Delta H: -0.0006
2026-01-01 21:06:38,167 | INFO | All windows processed in 14.84s
2026-01-01 21:06:38,169 | INFO | Saved QW_1660_v11_time_integration.json
2026-01-01 21:06:38,170 | INFO | QW-1660 v11 COMPLETE

{
  "integration_test": [
    {
      "window_sec": 32,
      "H1": 0.33995934773901576,
      "L1": 0.3520262031207114,
      "Delta": -0.012066855381695663,
      "compute_time_sec": 0.31324338912963867
    },
    {
      "window_sec": 64,
      "H1": 0.3105435282906503,
      "L1": 0.3263566705352001,
      "Delta": -0.015813142244549827,
      "compute_time_sec": 0.5691792964935303
    },
    {
      "window_sec": 128,
      "H1": 0.31001070121971774,
      "L1": 0.31835155974370716,
      "Delta": -0.008340858523989414,
      "compute_time_sec": 1.0603904724121094
    },
    {
      "window_sec": 256,
      "H1": 0.2881253580625984,
      "L1": 0.30091963081910567,
      "Delta": -0.012794272756507241,
      "compute_time_sec": 2.042760133743286
    },
    {
      "window_sec": 512,
      "H1": 0.2749408309566025,
      "L1": 0.2906625964784416,
      "Delta": -0.015721765521839126,
      "compute_time_sec": 3.687760591506958
    },
    {
      "window_sec": 1024,
      "H1": 0.28688700383670396,
      "L1": 0.287528864266169,
      "Delta": -0.0006418604294650687,
      "compute_time_sec": 7.145539283752441
    }
  ],
  "sample_rate": 4096.0,
  "total_runtime_sec": 14.842632055282593,
  "timestamp_utc": "2026-01-01T21:06:38.168655+00:00"
}

# ==============================================================================
# QW-1660 v12: DIURNAL FRACTAL MODULATION TEST (FIN CORE)
# Cel: Czy Hurst wykazuje modulację dobową (rotacja Ziemi względem struktury FIN)?
# Dane: ZAPISANY raw strain
# ==============================================================================

import h5py
import numpy as np
import json
import logging
import time
from scipy.signal import detrend
from datetime import datetime, timezone

# ------------------------------------------------------------
# LOGGING
# ------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)
log = logging.getLogger("QW-1660-v12")

log.info("START QW-1660 v12: DIURNAL FRACTAL MODULATION TEST")

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
RAW_DIR = "/kaggle/working/raw_strain"
H1_FILE = f"{RAW_DIR}/H1_raw.h5"
L1_FILE = f"{RAW_DIR}/L1_raw.h5"
DEFAULT_FS = 4096.0

SEGMENT_SEC = 3600  # 1 godzina
MAX_SEGMENTS = 24   # jedna doba

# ------------------------------------------------------------
# HURST (R/S)
# ------------------------------------------------------------
def hurst_rs(x):
    x = np.asarray(x)
    N = len(x)
    if N < 8000:
        return np.nan

    sizes = np.logspace(2, np.log10(N // 4), 10).astype(int)
    rs = []

    for s in sizes:
        n = N // s
        vals = []
        for i in range(n):
            seg = x[i*s:(i+1)*s]
            seg -= np.mean(seg)
            z = np.cumsum(seg)
            R = np.ptp(z)
            S = np.std(seg)
            if S > 0:
                vals.append(R / S)
        if vals:
            rs.append(np.mean(vals))

    if len(rs) < 4:
        return np.nan

    return np.polyfit(np.log(sizes[:len(rs)]), np.log(rs), 1)[0]

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------
def load_strain(path, label):
    log.info(f"Loading {label} from {path}")
    t0 = time.time()
    with h5py.File(path, "r") as f:
        x = detrend(f["strain"][:])
        fs = f["strain"].attrs.get("sample_rate", DEFAULT_FS)
    log.info(f"{label} loaded | N={len(x)} | fs={fs} | load={time.time()-t0:.2f}s")
    return x, fs

h1, fs = load_strain(H1_FILE, "H1")
l1, _  = load_strain(L1_FILE, "L1")

# ------------------------------------------------------------
# SEGMENTATION
# ------------------------------------------------------------
samples_per_seg = int(SEGMENT_SEC * fs)
n_segments = min(len(h1) // samples_per_seg, MAX_SEGMENTS)

log.info(f"Segmenting into {n_segments} hourly blocks")

results = []

for i in range(n_segments):
    log.info(f"Processing segment {i+1}/{n_segments}")
    s = i * samples_per_seg
    e = s + samples_per_seg

    H1_H = hurst_rs(h1[s:e])
    L1_H = hurst_rs(l1[s:e])

    results.append({
        "segment": i,
        "H1": H1_H,
        "L1": L1_H,
        "Delta": H1_H - L1_H
    })

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
out = {
    "diurnal_test": results,
    "segment_sec": SEGMENT_SEC,
    "sample_rate": fs,
    "timestamp_utc": datetime.now(timezone.utc).isoformat()
}

with open("QW_1660_v12_diurnal_modulation.json", "w") as f:
    json.dump(out, f, indent=2)

log.info("Saved QW_1660_v12_diurnal_modulation.json")
log.info("QW-1660 v12 COMPLETE")
2026-01-01 21:20:42,114 | INFO | START QW-1660 v12: DIURNAL FRACTAL MODULATION TEST
2026-01-01 21:20:42,115 | INFO | Loading H1 from /kaggle/working/raw_strain/H1_raw.h5
2026-01-01 21:20:43,218 | INFO | H1 loaded | N=16777216 | fs=4096.0 | load=1.10s
2026-01-01 21:20:43,220 | INFO | Loading L1 from /kaggle/working/raw_strain/L1_raw.h5
2026-01-01 21:20:44,323 | INFO | L1 loaded | N=16777216 | fs=4096.0 | load=1.10s
2026-01-01 21:20:44,326 | INFO | Segmenting into 1 hourly blocks
2026-01-01 21:20:44,327 | INFO | Processing segment 1/1
2026-01-01 21:21:07,395 | INFO | Saved QW_1660_v12_diurnal_modulation.json
2026-01-01 21:21:07,396 | INFO | QW-1660 v12 COMPLETE
{
  "diurnal_test": [
    {
      "segment": 0,
      "H1": 0.2861297793052872,
      "L1": 0.26324010867471026,
      "Delta": 0.022889670630576953
    }
  ],
  "segment_sec": 3600,
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-01T21:21:07.394453+00:00"
}
# ==============================================================================
# QW-1660 v13: INTER-DETECTOR FRACTAL COHERENCE TEST (FIN CRITICAL)
# Cel: Czy pamięć fraktalna jest wspólna dla H1 i L1?
# Metoda: Segmentacja + Hurst + time-slide null
# Dane: ZAPISANY raw strain
# ==============================================================================

import h5py
import numpy as np
import json
import logging
from scipy.signal import detrend
from scipy.stats import pearsonr
from datetime import datetime, timezone

# ------------------------------------------------------------
# LOGGING
# ------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)
log = logging.getLogger("QW-1660-v13")

log.info("START QW-1660 v13: INTER-DETECTOR FRACTAL COHERENCE TEST")

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
RAW_DIR = "/kaggle/working/raw_strain"
H1_FILE = f"{RAW_DIR}/H1_raw.h5"
L1_FILE = f"{RAW_DIR}/L1_raw.h5"
DEFAULT_FS = 4096.0

SEGMENT_SEC = 64      # krótkie, lokalne okna
TIME_SLIDES = 50      # null test

# ------------------------------------------------------------
# HURST (R/S)
# ------------------------------------------------------------
def hurst_rs(x):
    x = np.asarray(x)
    N = len(x)
    if N < 2000:
        return np.nan

    sizes = np.logspace(2, np.log10(N // 4), 8).astype(int)
    rs = []

    for s in sizes:
        n = N // s
        vals = []
        for i in range(n):
            seg = x[i*s:(i+1)*s]
            seg -= np.mean(seg)
            z = np.cumsum(seg)
            R = np.ptp(z)
            S = np.std(seg)
            if S > 0:
                vals.append(R / S)
        if vals:
            rs.append(np.mean(vals))

    if len(rs) < 3:
        return np.nan

    return np.polyfit(np.log(sizes[:len(rs)]), np.log(rs), 1)[0]

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------
def load_strain(path, label):
    log.info(f"Loading {label} from {path}")
    with h5py.File(path, "r") as f:
        x = detrend(f["strain"][:])
        fs = f["strain"].attrs.get("sample_rate", DEFAULT_FS)
    log.info(f"{label} loaded | N={len(x)} | fs={fs}")
    return x, fs

h1, fs = load_strain(H1_FILE, "H1")
l1, _  = load_strain(L1_FILE, "L1")

# ------------------------------------------------------------
# SEGMENTATION
# ------------------------------------------------------------
nseg = int(SEGMENT_SEC * fs)
N = min(len(h1), len(l1)) // nseg

log.info(f"Segmenting into {N} segments of {SEGMENT_SEC}s")

H1_H = []
L1_H = []

for i in range(N):
    H1_H.append(hurst_rs(h1[i*nseg:(i+1)*nseg]))
    L1_H.append(hurst_rs(l1[i*nseg:(i+1)*nseg]))

H1_H = np.array(H1_H)
L1_H = np.array(L1_H)

mask = ~np.isnan(H1_H) & ~np.isnan(L1_H)
H1_H = H1_H[mask]
L1_H = L1_H[mask]

# ------------------------------------------------------------
# REAL COHERENCE
# ------------------------------------------------------------
r_real, p_real = pearsonr(H1_H, L1_H)

# ------------------------------------------------------------
# NULL: TIME SLIDES
# ------------------------------------------------------------
r_null = []

for _ in range(TIME_SLIDES):
    shift = np.random.randint(5, len(L1_H)-5)
    L1_shift = np.roll(L1_H, shift)
    r, _ = pearsonr(H1_H, L1_shift)
    r_null.append(r)

r_null = np.array(r_null)

z = (r_real - np.mean(r_null)) / np.std(r_null)

# ------------------------------------------------------------
# RESULTS
# ------------------------------------------------------------
print("\n=== INTER-DETECTOR FRACTAL COHERENCE ===")
print(f"Real r:        {r_real:.3f} (p={p_real:.3e})")
print(f"Null mean r:   {np.mean(r_null):.3f}")
print(f"Null std:      {np.std(r_null):.3f}")
print(f"Z-score:       {z:.2f}")

if abs(z) > 3:
    verdict = "FRACTAL COHERENCE DETECTED — FIN CONSISTENT"
else:
    verdict = "NO COHERENCE — CONSISTENT WITH INDEPENDENT NOISE"

print("VERDICT:", verdict)

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
out = {
    "segment_sec": SEGMENT_SEC,
    "r_real": float(r_real),
    "p_real": float(p_real),
    "null_mean": float(np.mean(r_null)),
    "null_std": float(np.std(r_null)),
    "z_score": float(z),
    "verdict": verdict,
    "timestamp_utc": datetime.now(timezone.utc).isoformat()
}

with open("QW_1660_v13_fractal_coherence.json", "w") as f:
    json.dump(out, f, indent=2)

log.info("Saved QW_1660_v13_fractal_coherence.json")
log.info("QW-1660 v13 COMPLETE")
# ==============================================================================
# ------------------------------------------------------------
nseg = int(SEGMENT_SEC * fs)
N = min(len(h1), len(l1)) // nseg

log.info(f"Segmenting into {N} segments of {SEGMENT_SEC}s")

H1_H = []
L1_H = []

for i in range(N):
    H1_H.append(hurst_rs(h1[i*nseg:(i+1)*nseg]))
    L1_H.append(hurst_rs(l1[i*nseg:(i+1)*nseg]))

H1_H = np.array(H1_H)
L1_H = np.array(L1_H)

mask = ~np.isnan(H1_H) & ~np.isnan(L1_H)
H1_H = H1_H[mask]
L1_H = L1_H[mask]

# ------------------------------------------------------------
# REAL COHERENCE
# ------------------------------------------------------------
r_real, p_real = pearsonr(H1_H, L1_H)

# ------------------------------------------------------------
# NULL: TIME SLIDES
# ------------------------------------------------------------
r_null = []

for _ in range(TIME_SLIDES):
    shift = np.random.randint(5, len(L1_H)-5)
    L1_shift = np.roll(L1_H, shift)
    r, _ = pearsonr(H1_H, L1_shift)
    r_null.append(r)

r_null = np.array(r_null)

z = (r_real - np.mean(r_null)) / np.std(r_null)

# ------------------------------------------------------------
# RESULTS
# ------------------------------------------------------------
print("\n=== INTER-DETECTOR FRACTAL COHERENCE ===")
print(f"Real r:        {r_real:.3f} (p={p_real:.3e})")
print(f"Null mean r:   {np.mean(r_null):.3f}")
print(f"Null std:      {np.std(r_null):.3f}")
print(f"Z-score:       {z:.2f}")

if abs(z) > 3:
    verdict = "FRACTAL COHERENCE DETECTED — FIN CONSISTENT"
else:
    verdict = "NO COHERENCE — CONSISTENT WITH INDEPENDENT NOISE"

print("VERDICT:", verdict)

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
out = {
    "segment_sec": SEGMENT_SEC,
    "r_real": float(r_real),
    "p_real": float(p_real),
    "null_mean": float(np.mean(r_null)),
    "null_std": float(np.std(r_null)),
    "z_score": float(z),
    "verdict": verdict,
    "timestamp_utc": datetime.now(timezone.utc).isoformat()
}

with open("QW_1660_v13_fractal_coherence.json", "w") as f:
    json.dump(out, f, indent=2)

log.info("Saved QW_1660_v13_fractal_coherence.json")
log.info("QW-1660 v13 COMPLETE")

2026-01-01 21:51:26,349 | INFO | START QW-1660 v13: INTER-DETECTOR FRACTAL COHERENCE TEST
2026-01-01 21:51:26,352 | INFO | Loading H1 from /kaggle/working/raw_strain/H1_raw.h5
2026-01-01 21:51:27,769 | INFO | H1 loaded | N=16777216 | fs=4096.0
2026-01-01 21:51:27,770 | INFO | Loading L1 from /kaggle/working/raw_strain/L1_raw.h5
2026-01-01 21:51:28,982 | INFO | L1 loaded | N=16777216 | fs=4096.0
2026-01-01 21:51:28,984 | INFO | Segmenting into 64 segments of 64s
2026-01-01 21:51:55,414 | INFO | Saved QW_1660_v13_fractal_coherence.json
2026-01-01 21:51:55,416 | INFO | QW-1660 v13 COMPLETE


=== INTER-DETECTOR FRACTAL COHERENCE ===
Real r:        -0.145 (p=2.519e-01)
Null mean r:   0.024
Null std:      0.121
Z-score:       -1.41
VERDICT: NO COHERENCE — CONSISTENT WITH INDEPENDENT NOISE
{
  "segment_sec": 64,
  "r_real": -0.1453076356553426,
  "p_real": 0.25194302076116887,
  "null_mean": 0.024229854332111018,
  "null_std": 0.12065324834219607,
  "z_score": -1.4051630794606735,
  "verdict": "NO COHERENCE \u2014 CONSISTENT WITH INDEPENDENT NOISE",
  "timestamp_utc": "2026-01-01T21:51:55.414111+00:00"
}

# ==============================================================================
# QW-1660 v14: CROSS-FREQUENCY FRACTAL COUPLING TEST (FIN-CRITICAL)
# Cel: Czy struktura fraktalna jest hierarchiczna w CZĘSTOTLIWOŚCI?
# Dane: ZAPISANY raw strain (H1 + L1)
# ==============================================================================

import h5py
import numpy as np
import json
import logging
from scipy.signal import butter, filtfilt, detrend
from scipy.stats import pearsonr
from datetime import datetime, timezone

# ------------------------------------------------------------
# LOGGING
# ------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)
log = logging.getLogger("QW-1660-v14")

log.info("START QW-1660 v14: CROSS-FREQUENCY FRACTAL COUPLING TEST")

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
RAW_DIR = "/kaggle/working/raw_strain"
H1_FILE = f"{RAW_DIR}/H1_raw.h5"
L1_FILE = f"{RAW_DIR}/L1_raw.h5"
DEFAULT_FS = 4096.0

BANDS = [
    (30, 60),
    (60, 120),
    (120, 240),
    (240, 480),
]

# ------------------------------------------------------------
# HURST (R/S)
# ------------------------------------------------------------
def hurst_rs(x):
    x = np.asarray(x)
    N = len(x)
    if N < 4000:
        return np.nan

    sizes = np.logspace(2, np.log10(N // 4), 10).astype(int)
    rs = []

    for s in sizes:
        n = N // s
        vals = []
        for i in range(n):
            seg = x[i*s:(i+1)*s]
            seg -= np.mean(seg)
            z = np.cumsum(seg)
            R = np.ptp(z)
            S = np.std(seg)
            if S > 0:
                vals.append(R / S)
        if vals:
            rs.append(np.mean(vals))

    if len(rs) < 4:
        return np.nan

    return np.polyfit(np.log(sizes[:len(rs)]), np.log(rs), 1)[0]

# ------------------------------------------------------------
# FILTER
# ------------------------------------------------------------
def bandpass(x, fs, f1, f2):
    b, a = butter(4, [f1/(fs/2), f2/(fs/2)], btype="band")
    return filtfilt(b, a, x)

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------
def load_strain(path, label):
    log.info(f"Loading {label} from {path}")
    with h5py.File(path, "r") as f:
        x = detrend(f["strain"][:])
        fs = f["strain"].attrs.get("sample_rate", DEFAULT_FS)
    log.info(f"{label} loaded | N={len(x)} | fs={fs}")
    return x, fs

h1, fs = load_strain(H1_FILE, "H1")
l1, _  = load_strain(L1_FILE, "L1")

# ------------------------------------------------------------
# FRACTAL COUPLING
# ------------------------------------------------------------
log.info("Computing band-wise Hurst exponents")

H1_band = []
L1_band = []

for f1, f2 in BANDS:
    log.info(f"Processing band {f1}-{f2} Hz")
    H1_band.append(hurst_rs(bandpass(h1, fs, f1, f2)))
    L1_band.append(hurst_rs(bandpass(l1, fs, f1, f2)))

# ------------------------------------------------------------
# COUPLING TEST
# ------------------------------------------------------------
r_real, p_real = pearsonr(H1_band, L1_band)

# Null: shuffled bands
r_null = []
for _ in range(100):
    r_null.append(
        pearsonr(
            np.random.permutation(H1_band),
            np.random.permutation(L1_band)
        )[0]
    )

z = (r_real - np.mean(r_null)) / np.std(r_null)

# ------------------------------------------------------------
# RESULTS
# ------------------------------------------------------------
print("\n=== CROSS-FREQUENCY FRACTAL COUPLING ===")
print(f"Real r: {r_real:.3f} (p={p_real:.3e})")
print(f"Null mean r: {np.mean(r_null):.3f}")
print(f"Null std: {np.std(r_null):.3f}")
print(f"Z-score: {z:.2f}")

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
out = {
    "bands": BANDS,
    "H1_band": H1_band,
    "L1_band": L1_band,
    "r_real": r_real,
    "p_real": p_real,
    "null_mean": float(np.mean(r_null)),
    "null_std": float(np.std(r_null)),
    "z_score": z,
    "timestamp_utc": datetime.now(timezone.utc).isoformat()
}

with open("QW_1660_v14_cross_frequency.json", "w") as f:
    json.dump(out, f, indent=2)

log.info("Saved QW_1660_v14_cross_frequency.json")
log.info("QW-1660 v14 COMPLETE")
2026-01-01 22:01:33,171 | INFO | START QW-1660 v14: CROSS-FREQUENCY FRACTAL COUPLING TEST
2026-01-01 22:01:33,175 | INFO | Loading H1 from /kaggle/working/raw_strain/H1_raw.h5
2026-01-01 22:01:34,420 | INFO | H1 loaded | N=16777216 | fs=4096.0
2026-01-01 22:01:34,422 | INFO | Loading L1 from /kaggle/working/raw_strain/L1_raw.h5
2026-01-01 22:01:35,625 | INFO | L1 loaded | N=16777216 | fs=4096.0
2026-01-01 22:01:35,628 | INFO | Computing band-wise Hurst exponents
2026-01-01 22:01:35,630 | INFO | Processing band 30-60 Hz
2026-01-01 22:02:03,461 | INFO | Processing band 60-120 Hz
2026-01-01 22:02:30,910 | INFO | Processing band 120-240 Hz
2026-01-01 22:02:57,243 | INFO | Processing band 240-480 Hz
2026-01-01 22:03:23,169 | INFO | Saved QW_1660_v14_cross_frequency.json
2026-01-01 22:03:23,170 | INFO | QW-1660 v14 COMPLETE


=== CROSS-FREQUENCY FRACTAL COUPLING ===
Real r: 0.790 (p=2.096e-01)
Null mean r: -0.043
Null std: 0.608
Z-score: 1.37
{
  "bands": [
    [
      30,
      60
    ],
    [
      60,
      120
    ],
    [
      120,
      240
    ],
    [
      240,
      480
    ]
  ],
  "H1_band": [
    0.23478763001930075,
    0.20451298768398127,
    0.2605842628711187,
    0.205146703094412
  ],
  "L1_band": [
    0.2357252841122392,
    0.21592309821763106,
    0.24127755831025977,
    0.17635400825431122
  ],
  "r_real": 0.7903615769456582,
  "p_real": 0.20963842305434177,
  "null_mean": -0.04269731659612157,
  "null_std": 0.6081802287520204,
  "z_score": 1.369756618447147,
  "timestamp_utc": "2026-01-01T22:03:23.168967+00:00"
}
# ==============================================================================
# QW-1660 v15: SIDEREAL FRACTAL MODULATION TEST (FIN CORE)
# Cel: Czy H(t) koreluje z czasem gwiazdowym (orientacja Ziemi)?
# Dane: ZAPISANY raw strain
# ==============================================================================

import h5py
import numpy as np
import json
import logging
import time
from scipy.signal import detrend
from datetime import datetime, timezone

# ------------------------------------------------------------
# LOGGING
# ------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)
log = logging.getLogger("QW-1660-v15")

log.info("START QW-1660 v15: SIDEREAL FRACTAL MODULATION TEST")

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
RAW_DIR = "/kaggle/working/raw_strain"
H1_FILE = f"{RAW_DIR}/H1_raw.h5"
L1_FILE = f"{RAW_DIR}/L1_raw.h5"
DEFAULT_FS = 4096.0

SEGMENT_SEC = 1800   # 30 min
SIDEREAL_DAY = 86164.1  # seconds

# ------------------------------------------------------------
# HURST (R/S)
# ------------------------------------------------------------
def hurst_rs(x):
    x = np.asarray(x)
    N = len(x)
    if N < 4000:
        return np.nan

    sizes = np.logspace(2, np.log10(N // 4), 10).astype(int)
    rs = []

    for s in sizes:
        n = N // s
        vals = []
        for i in range(n):
            seg = x[i*s:(i+1)*s]
            seg -= np.mean(seg)
            z = np.cumsum(seg)
            R = np.ptp(z)
            S = np.std(seg)
            if S > 0:
                vals.append(R / S)
        if vals:
            rs.append(np.mean(vals))

    if len(rs) < 4:
        return np.nan

    return np.polyfit(np.log(sizes[:len(rs)]), np.log(rs), 1)[0]

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------
def load_strain(path, label):
    log.info(f"Loading {label} from {path}")
    t0 = time.time()
    with h5py.File(path, "r") as f:
        x = detrend(f["strain"][:])
        fs = f["strain"].attrs.get("sample_rate", DEFAULT_FS)
    log.info(f"{label} loaded | N={len(x)} | fs={fs} | load={time.time()-t0:.2f}s")
    return x, fs

h1, fs = load_strain(H1_FILE, "H1")
l1, _  = load_strain(L1_FILE, "L1")

# ------------------------------------------------------------
# SIDEREAL SEGMENTATION
# ------------------------------------------------------------
log.info("Segmenting data...")
seg_n = int(SEGMENT_SEC * fs)
n_seg = len(h1) // seg_n

results = []

for i in range(n_seg):
    t_start = i * SEGMENT_SEC
    sid_phase = (t_start % SIDEREAL_DAY) / SIDEREAL_DAY

    log.info(f"Segment {i+1}/{n_seg} | Sidereal phase={sid_phase:.3f}")

    H1_H = hurst_rs(h1[i*seg_n:(i+1)*seg_n])
    L1_H = hurst_rs(l1[i*seg_n:(i+1)*seg_n])

    results.append({
        "segment": i,
        "sidereal_phase": sid_phase,
        "H1": H1_H,
        "L1": L1_H,
        "Delta": H1_H - L1_H
    })

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
out = {
    "sidereal_test": results,
    "segment_sec": SEGMENT_SEC,
    "sample_rate": fs,
    "timestamp_utc": datetime.now(timezone.utc).isoformat()
}

with open("QW_1660_v15_sidereal_modulation.json", "w") as f:
    json.dump(out, f, indent=2)

log.info("Saved QW_1660_v15_sidereal_modulation.json")
log.info("QW-1660 v15 COMPLETE")
2026-01-01 22:26:46,873 | INFO | START QW-1660 v15: SIDEREAL FRACTAL MODULATION TEST
2026-01-01 22:26:46,876 | INFO | Loading H1 from /kaggle/working/raw_strain/H1_raw.h5
2026-01-01 22:26:48,067 | INFO | H1 loaded | N=16777216 | fs=4096.0 | load=1.19s
2026-01-01 22:26:48,075 | INFO | Loading L1 from /kaggle/working/raw_strain/L1_raw.h5
2026-01-01 22:26:49,217 | INFO | L1 loaded | N=16777216 | fs=4096.0 | load=1.14s
2026-01-01 22:26:49,223 | INFO | Segmenting data...
2026-01-01 22:26:49,224 | INFO | Segment 1/2 | Sidereal phase=0.000
2026-01-01 22:26:59,834 | INFO | Segment 2/2 | Sidereal phase=0.021
2026-01-01 22:27:10,426 | INFO | Saved QW_1660_v15_sidereal_modulation.json
2026-01-01 22:27:10,427 | INFO | QW-1660 v15 COMPLETE
{
  "sidereal_test": [
    {
      "segment": 0,
      "sidereal_phase": 0.0,
      "H1": 0.2764631582569697,
      "L1": 0.2738780068272483,
      "Delta": 0.002585151429721433
    },
    {
      "segment": 1,
      "sidereal_phase": 0.020890370815687738,
      "H1": 0.29013919118663656,
      "L1": 0.2689470865713718,
      "Delta": 0.021192104615264773
    }
  ],
  "segment_sec": 1800,
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-01T22:27:10.425300+00:00"
}

# ==============================================================================
# QW-1660 v16: SLIDING-WINDOW FRACTAL STATIONARITY TEST (FIN)
# Cel: Czy FIN jest stacjonarny w czasie?
# Metoda: przesuwne okno + histogram H
# Dane: ZAPISANY raw strain
# ==============================================================================

import h5py
import numpy as np
import json
import logging
import time
from scipy.signal import detrend
from datetime import datetime, timezone

# ------------------------------------------------------------
# LOGGING
# ------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)
log = logging.getLogger("QW-1660-v16")

log.info("START QW-1660 v16: SLIDING-WINDOW STATIONARITY TEST")

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
RAW_DIR = "/kaggle/working/raw_strain"
H1_FILE = f"{RAW_DIR}/H1_raw.h5"
DEFAULT_FS = 4096.0

WINDOW_SEC = 64
STEP_SEC = 32
MAX_WINDOWS = 200  # bezpieczeństwo czasowe

# ------------------------------------------------------------
# HURST
# ------------------------------------------------------------
def hurst_rs(x):
    x = np.asarray(x)
    N = len(x)
    if N < 4000:
        return np.nan

    sizes = np.logspace(2, np.log10(N // 4), 8).astype(int)
    rs = []

    for s in sizes:
        n = N // s
        vals = []
        for i in range(n):
            seg = x[i*s:(i+1)*s]
            seg -= np.mean(seg)
            z = np.cumsum(seg)
            R = np.ptp(z)
            S = np.std(seg)
            if S > 0:
                vals.append(R / S)
        if vals:
            rs.append(np.mean(vals))

    if len(rs) < 4:
        return np.nan

    return np.polyfit(np.log(sizes[:len(rs)]), np.log(rs), 1)[0]

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------
log.info(f"Loading H1 from {H1_FILE}")
t0 = time.time()
with h5py.File(H1_FILE, "r") as f:
    h1 = detrend(f["strain"][:])
    fs = f["strain"].attrs.get("sample_rate", DEFAULT_FS)
log.info(f"H1 loaded | N={len(h1)} | fs={fs} | load={time.time()-t0:.2f}s")

# ------------------------------------------------------------
# SLIDING WINDOW
# ------------------------------------------------------------
log.info("Starting sliding-window analysis")

W = int(WINDOW_SEC * fs)
S = int(STEP_SEC * fs)

Hs = []
times = []

for i, start in enumerate(range(0, len(h1) - W, S)):
    if i >= MAX_WINDOWS:
        break
    seg = h1[start:start+W]
    H = hurst_rs(seg)
    if not np.isnan(H):
        Hs.append(H)
        times.append(start / fs)
    if (i+1) % 20 == 0:
        log.info(f"Processed {i+1} windows")

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
out = {
    "window_sec": WINDOW_SEC,
    "step_sec": STEP_SEC,
    "N_windows": len(Hs),
    "H_mean": float(np.mean(Hs)),
    "H_std": float(np.std(Hs)),
    "timestamp_utc": datetime.now(timezone.utc).isoformat()
}

with open("QW_1660_v16_stationarity.json", "w") as f:
    json.dump(out, f, indent=2)

log.info("Saved QW_1660_v16_stationarity.json")
log.info("QW-1660 v16 COMPLETE")
2026-01-01 22:36:21,968 | INFO | START QW-1660 v16: SLIDING-WINDOW STATIONARITY TEST
2026-01-01 22:36:21,970 | INFO | Loading H1 from /kaggle/working/raw_strain/H1_raw.h5
2026-01-01 22:36:23,179 | INFO | H1 loaded | N=16777216 | fs=4096.0 | load=1.21s
2026-01-01 22:36:23,180 | INFO | Starting sliding-window analysis
2026-01-01 22:36:27,480 | INFO | Processed 20 windows
2026-01-01 22:36:31,714 | INFO | Processed 40 windows
2026-01-01 22:36:35,925 | INFO | Processed 60 windows
2026-01-01 22:36:40,174 | INFO | Processed 80 windows
2026-01-01 22:36:44,370 | INFO | Processed 100 windows
2026-01-01 22:36:48,692 | INFO | Processed 120 windows
2026-01-01 22:36:49,919 | INFO | Saved QW_1660_v16_stationarity.json
2026-01-01 22:36:49,920 | INFO | QW-1660 v16 COMPLETE
{
  "window_sec": 64,
  "step_sec": 32,
  "N_windows": 126,
  "H_mean": 0.306493493443769,
  "H_std": 0.011984760442150907,
  "timestamp_utc": "2026-01-01T22:36:49.918816+00:00"
}
# ==============================================================================
# QW-1660 v17: GLOBAL FRACTAL CONSISTENCY TEST (FIN CORE)
# Cel: Czy H(T) konwerguje do stałej? (prawo fundamentalne vs efekt lokalny)
# Dane: ZAPISANY raw strain (bez pobierania)
# ==============================================================================

import h5py
import numpy as np
import json
import logging
import time
from scipy.signal import detrend
from datetime import datetime, timezone

# ------------------------------------------------------------
# LOGGING
# ------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)
log = logging.getLogger("QW-1660-v17")

log.info("START QW-1660 v17: GLOBAL FRACTAL CONSISTENCY TEST")

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
RAW_DIR = "/kaggle/working/raw_strain"
H1_FILE = f"{RAW_DIR}/H1_raw.h5"
DEFAULT_FS = 4096.0

WINDOWS_SEC = [64, 128, 256, 512, 1024, 2048]

# ------------------------------------------------------------
# HURST (R/S)
# ------------------------------------------------------------
def hurst_rs(x):
    x = np.asarray(x)
    N = len(x)
    if N < 4000:
        return np.nan

    sizes = np.logspace(2, np.log10(N // 4), 10).astype(int)
    rs = []

    for s in sizes:
        n = N // s
        vals = []
        for i in range(n):
            seg = x[i*s:(i+1)*s]
            seg -= np.mean(seg)
            z = np.cumsum(seg)
            R = np.ptp(z)
            S = np.std(seg)
            if S > 0:
                vals.append(R / S)
        if vals:
            rs.append(np.mean(vals))

    if len(rs) < 4:
        return np.nan

    return np.polyfit(np.log(sizes[:len(rs)]), np.log(rs), 1)[0]

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------
log.info(f"Loading H1 from {H1_FILE}")
t0 = time.time()
with h5py.File(H1_FILE, "r") as f:
    h1 = detrend(f["strain"][:])
    fs = f["strain"].attrs.get("sample_rate", DEFAULT_FS)
log.info(f"H1 loaded | N={len(h1)} | fs={fs} | load={time.time()-t0:.2f}s")

# ------------------------------------------------------------
# CONSISTENCY TEST
# ------------------------------------------------------------
results = []

log.info("Starting scale-integration test...")

for T in WINDOWS_SEC:
    n = int(T * fs)
    if n > len(h1):
        log.warning(f"Skipping T={T}s (insufficient data)")
        continue

    log.info(f"Computing H(T) for T={T}s")
    t_start = time.time()
    H = hurst_rs(h1[:n])
    dt = time.time() - t_start

    log.info(f"  Result: H={H:.4f} | time={dt:.2f}s")

    results.append({
        "window_sec": T,
        "H": H,
        "samples": n,
        "compute_time_sec": dt
    })

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
out = {
    "global_consistency_test": results,
    "sample_rate": fs,
    "timestamp_utc": datetime.now(timezone.utc).isoformat()
}

with open("QW_1660_v17_global_consistency.json", "w") as f:
    json.dump(out, f, indent=2)

log.info("Saved QW_1660_v17_global_consistency.json")
log.info("QW-1660 v17 COMPLETE")
2026-01-01 22:53:07,338 | INFO | START QW-1660 v17: GLOBAL FRACTAL CONSISTENCY TEST
2026-01-01 22:53:07,340 | INFO | Loading H1 from /kaggle/working/raw_strain/H1_raw.h5
2026-01-01 22:53:08,498 | INFO | H1 loaded | N=16777216 | fs=4096.0 | load=1.16s
2026-01-01 22:53:08,499 | INFO | Starting scale-integration test...
2026-01-01 22:53:08,500 | INFO | Computing H(T) for T=64s
2026-01-01 22:53:08,746 | INFO |   Result: H=0.3089 | time=0.25s
2026-01-01 22:53:08,747 | INFO | Computing H(T) for T=128s
2026-01-01 22:53:09,197 | INFO |   Result: H=0.3098 | time=0.45s
2026-01-01 22:53:09,198 | INFO | Computing H(T) for T=256s
2026-01-01 22:53:10,072 | INFO |   Result: H=0.2886 | time=0.87s
2026-01-01 22:53:10,074 | INFO | Computing H(T) for T=512s
2026-01-01 22:53:11,705 | INFO |   Result: H=0.2757 | time=1.63s
2026-01-01 22:53:11,706 | INFO | Computing H(T) for T=1024s
2026-01-01 22:53:14,902 | INFO |   Result: H=0.2875 | time=3.20s
2026-01-01 22:53:14,903 | INFO | Computing H(T) for T=2048s
2026-01-01 22:53:20,954 | INFO |   Result: H=0.2710 | time=6.05s
2026-01-01 22:53:20,956 | INFO | Saved QW_1660_v17_global_consistency.json
2026-01-01 22:53:20,957 | INFO | QW-1660 v17 COMPLETE
{
  "global_consistency_test": [
    {
      "window_sec": 64,
      "H": 0.30885260430877476,
      "samples": 262144,
      "compute_time_sec": 0.24536657333374023
    },
    {
      "window_sec": 128,
      "H": 0.3098408247635631,
      "samples": 524288,
      "compute_time_sec": 0.4489321708679199
    },
    {
      "window_sec": 256,
      "H": 0.28857808692469833,
      "samples": 1048576,
      "compute_time_sec": 0.8733696937561035
    },
    {
      "window_sec": 512,
      "H": 0.2756928984409219,
      "samples": 2097152,
      "compute_time_sec": 1.6308715343475342
    },
    {
      "window_sec": 1024,
      "H": 0.2874670726409621,
      "samples": 4194304,
      "compute_time_sec": 3.195228099822998
    },
    {
      "window_sec": 2048,
      "H": 0.27098363828141603,
      "samples": 8388608,
      "compute_time_sec": 6.0503809452056885
    }
  ],
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-01T22:53:20.955802+00:00"
}
# ============================================================
# QW-1660 v18: MULTI-DETECTOR GLOBAL FRACTAL CONSISTENCY TEST
# ============================================================

import h5py
import numpy as np
import json
import time
import logging
from datetime import datetime

# ---------------- LOGGING SETUP ----------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)
log = logging.getLogger("QW-1660-v18")

# ---------------- PARAMETERS ----------------
DATA_FILES = {
    "H1": "/kaggle/working/raw_strain/H1_raw.h5",
    "L1": "/kaggle/working/raw_strain/L1_raw.h5"
}
FS = 4096.0
WINDOWS_SEC = [128, 256, 512, 1024]
OUTPUT_FILE = "QW_1660_v18_multi_detector_consistency.json"

# ---------------- HURST ESTIMATOR ----------------
def hurst_rs(x):
    x = np.asarray(x)
    N = len(x)
    if N < 100:
        return np.nan
    mean_x = np.mean(x)
    y = np.cumsum(x - mean_x)
    R = np.max(y) - np.min(y)
    S = np.std(x)
    return np.log(R / S) / np.log(N) if S > 0 else np.nan

# ---------------- LOAD DATA ----------------
strain_data = {}

for det, path in DATA_FILES.items():
    log.info(f"Loading {det} from {path}")
    t0 = time.time()
    with h5py.File(path, "r") as f:
        strain_data[det] = f["strain"][:]
    log.info(f"{det} loaded | N={len(strain_data[det])} | load={time.time()-t0:.2f}s")

# ---------------- ANALYSIS ----------------
results = []

for T in WINDOWS_SEC:
    samples = int(T * FS)
    log.info(f"Computing H(T) for T={T}s")
    entry = {"window_sec": T}

    for det, data in strain_data.items():
        segment = data[:samples]
        t0 = time.time()
        H = hurst_rs(segment)
        entry[det] = {
            "H": float(H),
            "samples": samples,
            "compute_time_sec": time.time() - t0
        }
        log.info(f"{det} | H={H:.4f}")

    results.append(entry)

# ---------------- SAVE ----------------
output = {
    "multi_detector_global_test": results,
    "sample_rate": FS,
    "timestamp_utc": datetime.utcnow().isoformat() + "Z"
}

with open(OUTPUT_FILE, "w") as f:
    json.dump(output, f, indent=2)

log.info(f"Saved {OUTPUT_FILE}")
log.info("QW-1660 v18 COMPLETE")
