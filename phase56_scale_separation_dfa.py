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
    with h5py.File(path, "r") as f:
        return f["strain"][:]

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
    log.info("START Phase 56: Exact Scale Separation")
    x_h1 = detrend(fetch_pure_strain("H1", GPS))
    x_l1 = detrend(fetch_pure_strain("L1", GPS))
    
    N_samples = len(x_h1)
    
    # Scale division: 10^5 samples is approx 24.4 seconds
    split_scale = 100000 
    
    # Short scales: e.g., 4096 (1s) to 100000 (~24s), avoiding really short scales to bypass high-freq noise
    scales_short = np.logspace(np.log10(4096), np.log10(split_scale), 15).astype(int)
    scales_short = np.unique(scales_short)
    
    # Long scales: > 100000 up to N/4 (~128s)
    scales_long = np.logspace(np.log10(split_scale), np.log10(N_samples//4), 15).astype(int)
    scales_long = np.unique(scales_long)
    
    # Full baseline scales for reference
    scales_full = np.logspace(3, np.log10(N_samples//4), 15).astype(int)
    scales_full = np.unique(scales_full)
    
    log.info(f"Computing Cross-H for SHORT scales ({scales_short[0]} to {scales_short[-1]} samples)...")
    h_cross_short = cross_mfdfa_q0(x_h1, x_l1, scales_short)
    
    log.info(f"Computing Cross-H for LONG scales ({scales_long[0]} to {scales_long[-1]} samples)...")
    h_cross_long = cross_mfdfa_q0(x_h1, x_l1, scales_long)
    
    log.info(f"Computing Cross-H for FULL scales ({scales_full[0]} to {scales_full[-1]} samples)...")
    h_cross_full = cross_mfdfa_q0(x_h1, x_l1, scales_full)

    out = {
        "H_cross_short_scales_1s_to_25s": h_cross_short,
        "H_cross_long_scales_25s_to_128s": h_cross_long,
        "H_cross_full_baseline": h_cross_full,
        "Verdict": ""
    }
    
    if h_cross_long < 0.35 and h_cross_short > 0.45:
        out["Verdict"] = "The anti-persistence is strictly a LONG-SCALE phenomenon (>25s), bypassing rapid feedback loops. Supports global/structural origin."
    elif h_cross_short < 0.35 and h_cross_long > 0.45:
        out["Verdict"] = "The anti-persistence is strictly a SHORT-SCALE phenomenon (<25s), indicative of fast feedback loop synchrony."
    else:
        out["Verdict"] = "Anti-persistence is scale-invariant across both short and long bins."
        
    with open("QW_1660_v56_ScaleSeparation.json", "w") as f:
        json.dump(out, f, indent=2)
        
    print("\n--- SCALE SEPARATION RESULTS ---")
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
