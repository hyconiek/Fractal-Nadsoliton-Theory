import os, json, logging
import numpy as np
import h5py
from gwpy.timeseries import TimeSeries
from scipy.signal import detrend

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

RAW_DIR = "./raw_strain_unfiltered"
FS = 4096
WINDOW_SEC = 512
GPS = 1266965117 # Baseline quiet GPS

def fetch_pure_strain(det, gps):
    path = f"{RAW_DIR}/{det}_unfiltered_{gps}.h5"
    if os.path.exists(path):
        with h5py.File(path, "r") as f:
            return f["strain"][:]
            
    ts = TimeSeries.fetch_open_data(det, gps, gps + WINDOW_SEC, verbose=False, cache=True)
    if ts.sample_rate.value > FS:
        ts = ts.resample(FS)
    ts = ts.notch(60).notch(120).notch(180) # Generic line removal
    val = ts.value
    with h5py.File(path, "w") as f:
        f.create_dataset("strain", data=val)
    return val

def cross_mfdfa_q0(x, y, scales):
    z = np.cumsum((x - np.mean(x)) * (y - np.mean(y)))
    F, valid_scales = [], []
    N = len(z)
    for s in scales:
        n = N // s
        if n == 0: continue
        rms = []
        for i in range(n):
            seg = z[i*s:(i+1)*s]
            p = np.polyfit(np.arange(s), seg, 1)
            var = np.mean((seg - np.polyval(p, np.arange(s)))**2)
            if var > 0: rms.append(var)
        rms = np.array(rms)
        if len(rms) > 0:
            F.append(np.exp(0.5 * np.mean(np.log(rms + 1e-300))))
            valid_scales.append(s)
    if len(F) < 3: return np.nan
    return float(np.polyfit(np.log(valid_scales), np.log(F), 1)[0])

def main():
    log.info("START Phase 62: Inter-Observatory Test (H1-V1)")
    
    # Try fetching V1 data. Note that V1 might not be online exactly at this GPS.
    # We will wrap it in try-except
    try:
        x_h1 = detrend(fetch_pure_strain("H1", GPS))
        x_v1 = detrend(fetch_pure_strain("V1", GPS))
    except Exception as e:
        log.error(f"Failed to fetch V1 data: {e}. Defaulting to a different GPS if needed.")
        return

    N = min(len(x_h1), len(x_v1))
    x_h1 = x_h1[:N]
    x_v1 = x_v1[:N]
    
    scales = np.unique(np.logspace(3, np.log10(N//4), 15).astype(int))
    
    h_cross = cross_mfdfa_q0(x_h1, x_v1, scales)
    
    log.info(f"H1-V1 Cross-H = {h_cross:.4f}")
    
    out = {
        "H1_V1_Cross_H": h_cross,
        "Verdict": "If ~0.50, the anomaly is an H1-L1 paired architectural artifact. If ~0.31, it's truly global across distinct architectures."
    }
    
    with open("QW_1660_v62_InterObservatory.json", "w") as f:
        json.dump(out, f, indent=2)

    print("\n--- INTER-OBSERVATORY TEST ---")
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
