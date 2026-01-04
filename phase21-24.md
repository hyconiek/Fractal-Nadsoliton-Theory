# ==============================================================================
# QW-1660 v21: FRACTAL STABILITY ACROSS EPOCHS (ROBUST FETCH - FIXED)
# Cel: Stabilność fraktalna między epokami (O3a / O3b)
# Dane: AUTOMATYCZNE POBIERANIE raw strain (Open Data)
# Estymator: MULTI-SCALE R/S (spójny z v10, v18, v19)
# ==============================================================================

!pip install --quiet gwpy h5py scipy

import os
import time
import json
import h5py
import numpy as np
import logging
from gwpy.timeseries import TimeSeries
from gwpy.segments import DataQualityFlag
from scipy.signal import detrend
from datetime import datetime, timezone

# ------------------------------------------------------------
# LOGGING
# ------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)
log = logging.getLogger("QW-1660-v21")
log.info("START QW-1660 v21: FRACTAL STABILITY ACROSS EPOCHS (FIXED)")

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
RAW_DIR = "/kaggle/working/raw_strain"
os.makedirs(RAW_DIR, exist_ok=True)

FS = 4096
WINDOW_SEC = 512

EPOCHS = {
    "O3a_H1": {
        "det": "H1",
        "gps_start": 1238166018,
        "gps_end":   1245946818
    },
    "O3b_H1": {
        "det": "H1",
        "gps_start": 1245946818,
        "gps_end":   1253946818  # Skrócony zakres, aby przyspieszyć pobieranie listy segmentów
    }
}

# ------------------------------------------------------------
# MULTI-SCALE R/S HURST
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
# FETCH & CACHE RAW STRAIN
# ------------------------------------------------------------
def fetch_and_save(epoch, cfg):
    path = f"{RAW_DIR}/{epoch}.h5"

    if os.path.exists(path):
        log.info(f"{epoch} already cached → {path}")
        return path

    log.info(f"{epoch} missing → FETCHING raw strain")

    try:
        # FIX: Używamy fetch_open_data zamiast query, aby uniknąć błędów auth i argumentów.
        # Flaga to np. 'H1_DATA' dla danych publicznych.
        dq_flag = f"{cfg['det']}_DATA"
        
        log.info(f"Checking segments availability for {dq_flag}...")
        
        # Pobieramy flagi dostępności danych publicznych
        segs = DataQualityFlag.fetch_open_data(
            dq_flag,
            cfg["gps_start"],
            cfg["gps_end"]
        )

        # Znajdujemy najdłuższy ciągły segment danych
        longest = max(segs.active, key=lambda s: s[1] - s[0])
        
        # Margines bezpieczeństwa +100s od początku segmentu
        gps0 = longest[0] + 100
        
        # Sprawdzenie czy segment jest wystarczająco długi
        if (longest[1] - gps0) < WINDOW_SEC:
            log.warning(f"Segment too short for {epoch}")
            return None

        log.info(f"{epoch} | Found segment. Using GPS {gps0} → {gps0 + WINDOW_SEC}")

        # Pobranie właściwego sygnału
        ts = TimeSeries.fetch_open_data(
            cfg["det"],
            gps0,
            gps0 + WINDOW_SEC,
            verbose=True
        )

        if ts.sample_rate.value > FS:
            ts = ts.resample(FS)

        ts = ts.notch(60).notch(120).notch(180).bandpass(20, 1000)

        with h5py.File(path, "w") as f:
            d = f.create_dataset("strain", data=ts.value)
            d.attrs["sample_rate"] = FS

        log.info(f"{epoch} saved → {path}")
        return path

    except Exception as e:
        log.error(f"{epoch} FETCH FAILED: {e}")
        return None

# ------------------------------------------------------------
# ANALYSIS
# ------------------------------------------------------------
results = []

for epoch, cfg in EPOCHS.items():
    path = fetch_and_save(epoch, cfg)
    if path is None:
        continue

    log.info(f"Loading {epoch} from {path}")
    t0 = time.time()

    with h5py.File(path, "r") as f:
        x = detrend(f["strain"][:])

    log.info(f"{epoch} loaded | N={len(x)} | load={time.time()-t0:.2f}s")

    t1 = time.time()
    H = hurst_rs(x)

    log.info(f"{epoch} | H={H:.4f} | compute={time.time()-t1:.2f}s")

    results.append({
        "epoch": epoch,
        "window_sec": WINDOW_SEC,
        "samples": len(x),
        "H": float(H)
    })

# ------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------
Hs = [r["H"] for r in results if not np.isnan(r["H"])]

summary = {
    "H_mean": float(np.mean(Hs)) if Hs else None,
    "H_std": float(np.std(Hs)) if len(Hs) > 1 else None,
    "N_epochs": len(Hs)
}

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
out = {
    "cross_epoch_fractal_stability": results,
    "summary": summary,
    "sample_rate": FS,
    "timestamp_utc": datetime.now(timezone.utc).isoformat()
}

with open("QW_1660_v21_cross_epoch_stability.json", "w") as f:
    json.dump(out, f, indent=2)

log.info("Saved QW_1660_v21_cross_epoch_stability.json")
log.info("QW-1660 v21 COMPLETE")

2026-01-04 00:26:02,344 | INFO | START QW-1660 v21: FRACTAL STABILITY ACROSS EPOCHS (FIXED)
2026-01-04 00:26:02,346 | INFO | O3a_H1 missing → FETCHING raw strain
2026-01-04 00:26:02,347 | INFO | Checking segments availability for H1_DATA...
2026-01-04 00:26:17,507 | INFO | O3a_H1 | Found segment. Using GPS 1244083082 → 1244083594

Fetched 1 URLs from gwosc.org for [1244083082 .. 1244083594))
Reading data... 

2026-01-04 00:26:26,159 | INFO | O3a_H1 saved → /kaggle/working/raw_strain/O3a_H1.h5
2026-01-04 00:26:26,159 | INFO | Loading O3a_H1 from /kaggle/working/raw_strain/O3a_H1.h5

[Done]

2026-01-04 00:26:26,236 | INFO | O3a_H1 loaded | N=2097152 | load=0.08s
2026-01-04 00:26:27,433 | INFO | O3a_H1 | H=0.2900 | compute=1.20s
2026-01-04 00:26:27,433 | INFO | O3b_H1 missing → FETCHING raw strain
2026-01-04 00:26:27,434 | INFO | Checking segments availability for H1_DATA...
2026-01-04 00:26:35,685 | INFO | O3b_H1 | Found segment. Using GPS 1249546813 → 1249547325

Fetched 1 URLs from gwosc.org for [1249546813 .. 1249547325))
Reading data... 

2026-01-04 00:26:46,558 | INFO | O3b_H1 saved → /kaggle/working/raw_strain/O3b_H1.h5
2026-01-04 00:26:46,559 | INFO | Loading O3b_H1 from /kaggle/working/raw_strain/O3b_H1.h5

[Done]

2026-01-04 00:26:46,642 | INFO | O3b_H1 loaded | N=2097152 | load=0.08s
2026-01-04 00:26:47,893 | INFO | O3b_H1 | H=0.3607 | compute=1.25s
2026-01-04 00:26:47,895 | INFO | Saved QW_1660_v21_cross_epoch_stability.json
2026-01-04 00:26:47,895 | INFO | QW-1660 v21 COMPLETE
{
  "cross_epoch_fractal_stability": [
    {
      "epoch": "O3a_H1",
      "window_sec": 512,
      "samples": 2097152,
      "H": 0.2900189379914974
    },
    {
      "epoch": "O3b_H1",
      "window_sec": 512,
      "samples": 2097152,
      "H": 0.36065896220129573
    }
  ],
  "summary": {
    "H_mean": 0.3253389500963966,
    "H_std": 0.03532001210489916,
    "N_epochs": 2
  },
  "sample_rate": 4096,
  "timestamp_utc": "2026-01-04T00:26:47.894554+00:00"
}

# ==============================================================================
# QW-1660 v22: FRACTAL STABILITY ACROSS EPOCHS (ROBUST FETCH - O3a/O3b/O4)
# Cel: Stabilność fraktalna między epokami (O3a / O3b / O4)
# Dane: AUTOMATYCZNE POBIERANIE raw strain (jak v2)
# Estymator: MULTI-SCALE R/S (spójny z v10, v18, v19)
# ==============================================================================

!pip install --quiet gwpy h5py scipy

import os
import time
import json
import h5py
import numpy as np
import logging
from gwpy.timeseries import TimeSeries
from gwpy.segments import DataQualityFlag
from scipy.signal import detrend
from datetime import datetime, timezone

# ------------------------------------------------------------
# LOGGING
# ------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)
log = logging.getLogger("QW-1660-v22")
log.info("START QW-1660 v22: FRACTAL STABILITY ACROSS EPOCHS")

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
RAW_DIR = "/kaggle/working/raw_strain"
os.makedirs(RAW_DIR, exist_ok=True)

FS = 4096
WINDOW_SEC = 512

# Definicja epok i zakresów GPS
EPOCHS = {
    "O3a_H1": {
        "det": "H1",
        "gps_start": 1238166018,
        "gps_end":   1245946818
    },
    "O3b_H1": {
        "det": "H1",
        "gps_start": 1245946818,
        "gps_end":   1269363618
    },
    "O4_H1": {
        "det": "H1",
        "gps_start": 1368979218,
        "gps_end":   1388534418 
    }
}

# ------------------------------------------------------------
# MULTI-SCALE R/S HURST
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
# FETCH & CACHE RAW STRAIN
# ------------------------------------------------------------
def fetch_and_save(epoch, cfg):
    path = f"{RAW_DIR}/{epoch}.h5"

    if os.path.exists(path):
        log.info(f"{epoch} already cached → {path}")
        return path

    log.info(f"{epoch} missing → FETCHING raw strain")

    try:
        # FIX: Usunięto argument 'verbose', który powodował błąd w get_segments()
        log.info(f"Querying segments for {epoch} ({cfg['gps_start']} - {cfg['gps_end']})...")
        segs = DataQualityFlag.fetch_open_data(
            f"{cfg['det']}_DATA",
            cfg["gps_start"],
            cfg["gps_end"]
        )

        # Sprawdzenie czy znaleziono jakiekolwiek segmenty
        if len(segs.active) == 0:
            log.warning(f"{epoch} | No active segments found (data might not be public yet).")
            return None

        # Wybór najdłuższego segmentu
        longest = max(segs.active, key=lambda s: s[1] - s[0])
        gps0 = longest[0] + 100

        log.info(f"{epoch} | Found segment, downloading GPS {gps0} → {gps0 + WINDOW_SEC}")

        ts = TimeSeries.fetch_open_data(
            cfg["det"],
            gps0,
            gps0 + WINDOW_SEC,
            verbose=True
        )

        if ts.sample_rate.value > FS:
            ts = ts.resample(FS)

        ts = ts.notch(60).notch(120).notch(180).bandpass(20, 1000)

        with h5py.File(path, "w") as f:
            d = f.create_dataset("strain", data=ts.value)
            d.attrs["sample_rate"] = FS
            d.attrs["gps_start"] = gps0

        log.info(f"{epoch} saved → {path}")
        return path

    except Exception as e:
        log.error(f"{epoch} FETCH FAILED: {e}")
        return None

# ------------------------------------------------------------
# ANALYSIS
# ------------------------------------------------------------
results = []

for epoch, cfg in EPOCHS.items():
    path = fetch_and_save(epoch, cfg)
    if path is None:
        continue

    log.info(f"Loading {epoch} from {path}")
    t0 = time.time()

    with h5py.File(path, "r") as f:
        x = detrend(f["strain"][:])

    log.info(f"{epoch} loaded | N={len(x)} | load={time.time()-t0:.2f}s")

    t1 = time.time()
    H = hurst_rs(x)

    log.info(f"{epoch} | H={H:.4f} | compute={time.time()-t1:.2f}s")

    results.append({
        "epoch": epoch,
        "window_sec": WINDOW_SEC,
        "samples": len(x),
        "H": float(H)
    })

# ------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------
Hs = [r["H"] for r in results if not np.isnan(r["H"])]

summary = {
    "H_mean": float(np.mean(Hs)) if Hs else None,
    "H_std": float(np.std(Hs)) if len(Hs) > 1 else None,
    "N_epochs": len(Hs),
    "epochs_processed": [r["epoch"] for r in results]
}

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
out = {
    "cross_epoch_fractal_stability": results,
    "summary": summary,
    "sample_rate": FS,
    "timestamp_utc": datetime.now(timezone.utc).isoformat()
}

with open("QW_1660_v22_cross_epoch_stability.json", "w") as f:
    json.dump(out, f, indent=2)

log.info("Saved QW_1660_v22_cross_epoch_stability.json")
log.info("QW-1660 v22 COMPLETE")

2026-01-04 00:32:38,711 | INFO | START QW-1660 v22: FRACTAL STABILITY ACROSS EPOCHS
2026-01-04 00:32:38,714 | INFO | O3a_H1 already cached → /kaggle/working/raw_strain/O3a_H1.h5
2026-01-04 00:32:38,714 | INFO | Loading O3a_H1 from /kaggle/working/raw_strain/O3a_H1.h5
2026-01-04 00:32:38,814 | INFO | O3a_H1 loaded | N=2097152 | load=0.10s
2026-01-04 00:32:40,035 | INFO | O3a_H1 | H=0.2900 | compute=1.22s
2026-01-04 00:32:40,035 | INFO | O3b_H1 already cached → /kaggle/working/raw_strain/O3b_H1.h5
2026-01-04 00:32:40,036 | INFO | Loading O3b_H1 from /kaggle/working/raw_strain/O3b_H1.h5
2026-01-04 00:32:40,112 | INFO | O3b_H1 loaded | N=2097152 | load=0.08s
2026-01-04 00:32:41,322 | INFO | O3b_H1 | H=0.3607 | compute=1.21s
2026-01-04 00:32:41,323 | INFO | O4_H1 missing → FETCHING raw strain
2026-01-04 00:32:41,324 | INFO | Querying segments for O4_H1 (1368979218 - 1388534418)...
2026-01-04 00:33:02,806 | INFO | O4_H1 | Found segment, downloading GPS 1376866665 → 1376867177

Fetched 1 URLs from gwosc.org for [1376866665 .. 1376867177))
Reading data... 

2026-01-04 00:33:12,124 | INFO | O4_H1 saved → /kaggle/working/raw_strain/O4_H1.h5
2026-01-04 00:33:12,126 | INFO | Loading O4_H1 from /kaggle/working/raw_strain/O4_H1.h5

[Done]

2026-01-04 00:33:12,213 | INFO | O4_H1 loaded | N=2097152 | load=0.09s
2026-01-04 00:33:13,408 | INFO | O4_H1 | H=0.3301 | compute=1.19s
2026-01-04 00:33:13,410 | INFO | Saved QW_1660_v22_cross_epoch_stability.json
2026-01-04 00:33:13,411 | INFO | QW-1660 v22 COMPLETE
{
  "cross_epoch_fractal_stability": [
    {
      "epoch": "O3a_H1",
      "window_sec": 512,
      "samples": 2097152,
      "H": 0.2900189379914974
    },
    {
      "epoch": "O3b_H1",
      "window_sec": 512,
      "samples": 2097152,
      "H": 0.36065896220129573
    },
    {
      "epoch": "O4_H1",
      "window_sec": 512,
      "samples": 2097152,
      "H": 0.33009336350893775
    }
  ],
  "summary": {
    "H_mean": 0.32692375456724365,
    "H_std": 0.0289256295894912,
    "N_epochs": 3,
    "epochs_processed": [
      "O3a_H1",
      "O3b_H1",
      "O4_H1"
    ]
  },
  "sample_rate": 4096,
  "timestamp_utc": "2026-01-04T00:33:13.409766+00:00"
}
# ==============================================================================
# QW-1660 v23: CROSS-HURST FRACTAL CONSISTENCY (H1–L1–V1)
# Cel: Sprawdzenie fraktalnej korelacji MIĘDZY detektorami
# Metoda: CROSS-R/S (wspólne okna, wspólny GPS)
# Dane: AUTOMATYCZNE POBIERANIE raw strain (jak v22)
# ==============================================================================

!pip install --quiet gwpy h5py scipy

import os
import time
import json
import h5py
import numpy as np
import logging
from gwpy.timeseries import TimeSeries
from gwpy.segments import DataQualityFlag
from scipy.signal import detrend
from datetime import datetime, timezone

# ------------------------------------------------------------
# LOGGING
# ------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)
log = logging.getLogger("QW-1660-v23")
log.info("START QW-1660 v23: CROSS-HURST FRACTAL CONSISTENCY")

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
RAW_DIR = "/kaggle/working/raw_strain"
os.makedirs(RAW_DIR, exist_ok=True)

FS = 4096
WINDOW_SEC = 512

# Analiza na O3 (najpewniejsze wspólne dane H1–L1–V1)
EPOCH = {
    "gps_start": 1238166018,
    "gps_end":   1269363618
}

DETECTORS = ["H1", "L1", "V1"]

# ------------------------------------------------------------
# CROSS-R/S HURST
# ------------------------------------------------------------
def cross_hurst_rs(x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    N = min(len(x), len(y))
    if N < 4000:
        return np.nan

    sizes = np.logspace(2, np.log10(N // 4), 10).astype(int)
    rs = []

    for s in sizes:
        n = N // s
        vals = []
        for i in range(n):
            xs = x[i*s:(i+1)*s] - np.mean(x[i*s:(i+1)*s])
            ys = y[i*s:(i+1)*s] - np.mean(y[i*s:(i+1)*s])

            z = np.cumsum(xs + ys)  # wspólna trajektoria
            R = np.ptp(z)
            S = np.std(xs + ys)

            if S > 0:
                vals.append(R / S)

        if vals:
            rs.append(np.mean(vals))

    if len(rs) < 4:
        return np.nan

    return np.polyfit(np.log(sizes[:len(rs)]), np.log(rs), 1)[0]

# ------------------------------------------------------------
# FETCH & CACHE (jak v22)
# ------------------------------------------------------------
def fetch_detector(det):
    path = f"{RAW_DIR}/{det}_shared.h5"

    if os.path.exists(path):
        log.info(f"{det} cached → {path}")
        return path

    log.info(f"{det} missing → FETCHING raw strain")

    segs = DataQualityFlag.fetch_open_data(
        f"{det}_DATA",
        EPOCH["gps_start"],
        EPOCH["gps_end"]
    )

    if len(segs.active) == 0:
        log.warning(f"{det}: no public segments")
        return None

    longest = max(segs.active, key=lambda s: s[1] - s[0])
    gps0 = longest[0] + 100

    ts = TimeSeries.fetch_open_data(
        det,
        gps0,
        gps0 + WINDOW_SEC,
        verbose=True
    )

    if ts.sample_rate.value > FS:
        ts = ts.resample(FS)

    ts = ts.notch(60).notch(120).notch(180).bandpass(20, 1000)

    with h5py.File(path, "w") as f:
        d = f.create_dataset("strain", data=ts.value)
        d.attrs["sample_rate"] = FS
        d.attrs["gps_start"] = gps0

    log.info(f"{det} saved → {path}")
    return path

# ------------------------------------------------------------
# LOAD ALL DETECTORS
# ------------------------------------------------------------
data = {}

for det in DETECTORS:
    path = fetch_detector(det)
    if path is None:
        continue

    with h5py.File(path, "r") as f:
        data[det] = detrend(f["strain"][:])

# ------------------------------------------------------------
# CROSS-HURST ANALYSIS
# ------------------------------------------------------------
pairs = [("H1", "L1"), ("H1", "V1"), ("L1", "V1")]
results = []

for a, b in pairs:
    if a not in data or b not in data:
        continue

    log.info(f"Computing cross-Hurst {a}–{b}")
    t0 = time.time()
    Hxy = cross_hurst_rs(data[a], data[b])

    log.info(f"{a}-{b} | Hxy={Hxy:.4f} | t={time.time()-t0:.2f}s")

    results.append({
        "pair": f"{a}-{b}",
        "H_cross": float(Hxy),
        "samples": min(len(data[a]), len(data[b]))
    })

# ------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------
Hs = [r["H_cross"] for r in results if not np.isnan(r["H_cross"])]

summary = {
    "H_cross_mean": float(np.mean(Hs)) if Hs else None,
    "H_cross_std": float(np.std(Hs)) if len(Hs) > 1 else None,
    "pairs": [r["pair"] for r in results]
}

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
out = {
    "cross_hurst_results": results,
    "summary": summary,
    "sample_rate": FS,
    "window_sec": WINDOW_SEC,
    "timestamp_utc": datetime.now(timezone.utc).isoformat()
}

with open("QW_1660_v23_cross_hurst.json", "w") as f:
    json.dump(out, f, indent=2)

log.info("Saved QW_1660_v23_cross_hurst.json")
log.info("QW-1660 v23 COMPLETE")
2026-01-04 00:59:33,394 | INFO | START QW-1660 v23: CROSS-HURST FRACTAL CONSISTENCY
2026-01-04 00:59:33,396 | INFO | H1 missing → FETCHING raw strain

Fetched 1 URLs from gwosc.org for [1266965017 .. 1266965529))
Reading data... 

2026-01-04 01:00:01,992 | INFO | H1 saved → /kaggle/working/raw_strain/H1_shared.h5

[Done]

2026-01-04 01:00:02,076 | INFO | L1 missing → FETCHING raw strain

Fetched 2 URLs from gwosc.org for [1266884513 .. 1266885025))
Reading data... [Done]
Reading data... 

2026-01-04 01:00:36,301 | INFO | L1 saved → /kaggle/working/raw_strain/L1_shared.h5

[Done]

2026-01-04 01:00:36,393 | INFO | V1 missing → FETCHING raw strain
/usr/local/lib/python3.12/dist-packages/gwosc/api/v2.py:83: UserWarning: Your request will need to fetch for 128 pages. Try to constraint your request if possible.
  warnings.warn(warn_message)
2026-01-04 01:00:36,871 | WARNING | Your request will need to fetch for 128 pages. Try to constraint your request if possible.

Fetched 1 URLs from gwosc.org for [1262676141 .. 1262676653))
Reading data... 

2026-01-04 01:01:47,333 | INFO | V1 saved → /kaggle/working/raw_strain/V1_shared.h5

[Done]

2026-01-04 01:01:47,421 | INFO | Computing cross-Hurst H1–L1
2026-01-04 01:01:48,932 | INFO | H1-L1 | Hxy=0.1155 | t=1.51s
2026-01-04 01:01:48,933 | INFO | Computing cross-Hurst H1–V1
2026-01-04 01:01:50,412 | INFO | H1-V1 | Hxy=0.0780 | t=1.48s
2026-01-04 01:01:50,413 | INFO | Computing cross-Hurst L1–V1
2026-01-04 01:01:51,855 | INFO | L1-V1 | Hxy=0.1907 | t=1.44s
2026-01-04 01:01:51,856 | INFO | Saved QW_1660_v23_cross_hurst.json
2026-01-04 01:01:51,857 | INFO | QW-1660 v23 COMPLETE
{
  "cross_hurst_results": [
    {
      "pair": "H1-L1",
      "H_cross": 0.11552275914205934,
      "samples": 2097152
    },
    {
      "pair": "H1-V1",
      "H_cross": 0.07803517778157153,
      "samples": 2097152
    },
    {
      "pair": "L1-V1",
      "H_cross": 0.19066426135733577,
      "samples": 2097152
    }
  ],
  "summary": {
    "H_cross_mean": 0.12807406609365554,
    "H_cross_std": 0.04682932910348679,
    "pairs": [
      "H1-L1",
      "H1-V1",
      "L1-V1"
    ]
  },
  "sample_rate": 4096,
  "window_sec": 512,
  "timestamp_utc": "2026-01-04T01:01:51.856503+00:00"
}

# ==============================================================================
# QW-1660 v24: CROSS-HURST NULL-MODEL VALIDATION
# Cel: Walidacja istotności cross-Hurst (H1-L1 / H1-V1 / L1-V1)
# Metoda: REAL vs SHUFFLE vs PHASE-RANDOMIZED
# Estymator: MULTI-SCALE CROSS R/S (spójny z v23)
# ==============================================================================

!pip install --quiet gwpy h5py scipy

import os
import time
import json
import h5py
import numpy as np
import logging
from scipy.signal import detrend
from scipy.fft import rfft, irfft
from datetime import datetime, timezone

# ------------------------------------------------------------
# LOGGING
# ------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s"
)
log = logging.getLogger("QW-1660-v24")
log.info("START QW-1660 v24: CROSS-HURST NULL-MODEL VALIDATION")

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
RAW_DIR = "/kaggle/working/raw_strain"
FS = 4096
WINDOW_SEC = 512
N = FS * WINDOW_SEC

PAIRS = [
    ("H1", "L1"),
    ("H1", "V1"),
    ("L1", "V1")
]

FILES = {
    "H1": f"{RAW_DIR}/H1_shared.h5",
    "L1": f"{RAW_DIR}/L1_shared.h5",
    "V1": f"{RAW_DIR}/V1_shared.h5"
}

# ------------------------------------------------------------
# MULTI-SCALE CROSS R/S HURST
# ------------------------------------------------------------
def cross_hurst_rs(x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    N = min(len(x), len(y))

    sizes = np.logspace(2, np.log10(N // 4), 10).astype(int)
    rs = []

    for s in sizes:
        n = N // s
        vals = []
        for i in range(n):
            xs = x[i*s:(i+1)*s] - np.mean(x[i*s:(i+1)*s])
            ys = y[i*s:(i+1)*s] - np.mean(y[i*s:(i+1)*s])

            z = np.cumsum(xs * ys)
            R = np.ptp(z)
            S = np.std(xs) * np.std(ys)

            if S > 0:
                vals.append(R / S)

        if vals:
            rs.append(np.mean(vals))

    if len(rs) < 4:
        return np.nan

    return np.polyfit(np.log(sizes[:len(rs)]), np.log(rs), 1)[0]

# ------------------------------------------------------------
# NULL MODELS
# ------------------------------------------------------------
def shuffle(x):
    return np.random.permutation(x)

def phase_randomize(x):
    X = rfft(x)
    phases = np.exp(1j * np.random.uniform(0, 2*np.pi, len(X)))
    return irfft(np.abs(X) * phases, n=len(x))

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------
data = {}
for det, path in FILES.items():
    log.info(f"Loading {det} from {path}")
    with h5py.File(path, "r") as f:
        data[det] = detrend(f["strain"][:N])

# ------------------------------------------------------------
# ANALYSIS
# ------------------------------------------------------------
results = []

for d1, d2 in PAIRS:
    log.info(f"Processing pair {d1}-{d2}")

    x = data[d1]
    y = data[d2]

    t0 = time.time()
    H_real = cross_hurst_rs(x, y)

    H_shuffle = cross_hurst_rs(shuffle(x), shuffle(y))
    H_phase = cross_hurst_rs(phase_randomize(x), phase_randomize(y))

    dt = time.time() - t0

    log.info(
        f"{d1}-{d2} | H_real={H_real:.4f} | "
        f"H_shuffle={H_shuffle:.4f} | H_phase={H_phase:.4f}"
    )

    results.append({
        "pair": f"{d1}-{d2}",
        "H_real": float(H_real),
        "H_shuffle": float(H_shuffle),
        "H_phase_randomized": float(H_phase),
        "delta_real_shuffle": float(H_real - H_shuffle),
        "delta_real_phase": float(H_real - H_phase),
        "compute_time_sec": dt
    })

# ------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------
deltas = [r["delta_real_phase"] for r in results]

summary = {
    "mean_delta_real_phase": float(np.mean(deltas)),
    "std_delta_real_phase": float(np.std(deltas)),
    "N_pairs": len(results)
}

# ------------------------------------------------------------
# SAVE
# ------------------------------------------------------------
out = {
    "cross_hurst_null_validation": results,
    "summary": summary,
    "sample_rate": FS,
    "window_sec": WINDOW_SEC,
    "timestamp_utc": datetime.now(timezone.utc).isoformat()
}

with open("QW_1660_v24_cross_hurst_null_validation.json", "w") as f:
    json.dump(out, f, indent=2)

log.info("Saved QW_1660_v24_cross_hurst_null_validation.json")
log.info("QW-1660 v24 COMPLETE")

2026-01-04 01:25:05,964 | INFO | START QW-1660 v24: CROSS-HURST NULL-MODEL VALIDATION
2026-01-04 01:25:05,967 | INFO | Loading H1 from /kaggle/working/raw_strain/H1_shared.h5
2026-01-04 01:25:06,067 | INFO | Loading L1 from /kaggle/working/raw_strain/L1_shared.h5
2026-01-04 01:25:06,206 | INFO | Loading V1 from /kaggle/working/raw_strain/V1_shared.h5
2026-01-04 01:25:06,335 | INFO | Processing pair H1-L1
2026-01-04 01:25:12,583 | INFO | H1-L1 | H_real=0.3611 | H_shuffle=0.5310 | H_phase=0.3968
2026-01-04 01:25:12,583 | INFO | Processing pair H1-V1
2026-01-04 01:25:18,841 | INFO | H1-V1 | H_real=0.2151 | H_shuffle=0.5220 | H_phase=0.2362
2026-01-04 01:25:18,841 | INFO | Processing pair L1-V1
2026-01-04 01:25:25,094 | INFO | L1-V1 | H_real=0.2700 | H_shuffle=0.5040 | H_phase=0.2834
2026-01-04 01:25:25,096 | INFO | Saved QW_1660_v24_cross_hurst_null_validation.json
2026-01-04 01:25:25,097 | INFO | QW-1660 v24 COMPLETE
{
  "cross_hurst_null_validation": [
    {
      "pair": "H1-L1",
      "H_real": 0.3610605334380803,
      "H_shuffle": 0.5310469043003473,
      "H_phase_randomized": 0.39679047633700987,
      "delta_real_shuffle": -0.16998637086226698,
      "delta_real_phase": -0.035729942898929556,
      "compute_time_sec": 6.245482921600342
    },
    {
      "pair": "H1-V1",
      "H_real": 0.2150698637573661,
      "H_shuffle": 0.5219537647569182,
      "H_phase_randomized": 0.23618825276589883,
      "delta_real_shuffle": -0.30688390099955215,
      "delta_real_phase": -0.021118389008532745,
      "compute_time_sec": 6.256558418273926
    },
    {
      "pair": "L1-V1",
      "H_real": 0.27003525002927997,
      "H_shuffle": 0.5039965798656929,
      "H_phase_randomized": 0.28341320075874976,
      "delta_real_shuffle": -0.23396132983641293,
      "delta_real_phase": -0.01337795072946979,
      "compute_time_sec": 6.2524378299713135
    }
  ],
  "summary": {
    "mean_delta_real_phase": -0.023408760878977364,
    "std_delta_real_phase": 0.00926776639630034,
    "N_pairs": 3
  },
  "sample_rate": 4096,
  "window_sec": 512,
  "timestamp_utc": "2026-01-04T01:25:25.095853+00:00"
}
