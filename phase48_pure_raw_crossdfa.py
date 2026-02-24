import os, json, logging, time
import numpy as np
import h5py
from gwpy.timeseries import TimeSeries
from scipy.signal import detrend

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
log = logging.getLogger("QW-1660-v48")
log.info("START QW-1660 v48: PURE RAW CROSS-MF-DFA (NO BANDPASS)")

RAW_DIR = "./raw_strain_unfiltered"
os.makedirs(RAW_DIR, exist_ok=True)
FS = 4096
WINDOW_SEC = 512
GPS = 1266965117 # known good GPS from previous scans

def cross_mfdfa_q0(x, y, scales):
    # Cross-covariance profile
    z = np.cumsum((x - np.mean(x)) * (y - np.mean(y)))
    F = []
    valid_scales = []
    N = len(z)
    for s in scales:
        n = N // s
        if n == 0: continue
        rms = []
        for i in range(n):
            seg = z[i*s:(i+1)*s]
            p = np.polyfit(np.arange(s), seg, 1)
            var = np.mean((seg - np.polyval(p, np.arange(s)))**2)
            if var > 0:
                rms.append(var)
        rms = np.array(rms)
        if len(rms) > 0:
            F.append(np.exp(0.5 * np.mean(np.log(rms + 1e-300))))
            valid_scales.append(s)
    
    if len(F) < 3: return np.nan
    H = np.polyfit(np.log(valid_scales), np.log(F), 1)[0]
    return float(H)

def fetch_pure_strain(det, gps):
    path = f"{RAW_DIR}/{det}_unfiltered.h5"
    if os.path.exists(path):
        log.info(f"Loading cached pure {det} strain from {path}")
        with h5py.File(path, "r") as f:
            return f["strain"][:]
            
    log.info(f"Fetching pure {det} strain...")
    try:
        ts = TimeSeries.fetch_open_data(det, gps, gps + WINDOW_SEC, verbose=False)
        if ts.sample_rate.value > FS:
            ts = ts.resample(FS)
        
        # APPLY ONLY NOTCH FILTERS (NO BANDPASS)
        ts = ts.notch(60).notch(120).notch(180)
        
        val = ts.value
        with h5py.File(path, "w") as f:
            f.create_dataset("strain", data=val)
        return val
    except Exception as e:
        log.error(f"Failed to fetch {det} data: {e}")
        return None

def main():
    # 1. Fetch pure unfiltered data for both detectors
    x_h1 = fetch_pure_strain("H1", GPS)
    x_l1 = fetch_pure_strain("L1", GPS)
    
    if x_h1 is None or x_l1 is None:
        log.error("Aborting, missing data.")
        return
        
    x_h1 = detrend(x_h1)
    x_l1 = detrend(x_l1)

    # We use the original long scales (10^3 to N/4 samples)
    # This probes the long-range behavior (low frequencies)
    scales_long = np.logspace(3, np.log10(len(x_h1)//4), 15).astype(int)
    scales_long = np.unique(scales_long)
    
    log.info("Computing Cross-MF-DFA H(q=0) on PURE RAW data...")
    H_cross = cross_mfdfa_q0(x_h1, x_l1, scales_long)
    log.info(f"Pure Cross-H = {H_cross:.4f}")

    # Compute individual pure H for reference
    def mfdfa_q0(x, scales):
        y = np.cumsum(x - np.mean(x))
        F = []
        valid_scales = []
        N = len(y)
        for s in scales:
            n = N // s
            if n == 0: continue
            rms = []
            for i in range(n):
                seg = y[i*s:(i+1)*s]
                p = np.polyfit(np.arange(s), seg, 1)
                var = np.mean((seg - np.polyval(p, np.arange(s)))**2)
                if var > 0:
                    rms.append(var)
            rms = np.array(rms)
            if len(rms) > 0:
                F.append(np.exp(0.5 * np.mean(np.log(rms + 1e-300))))
                valid_scales.append(s)
        if len(F) < 3: return np.nan
        return float(np.polyfit(np.log(valid_scales), np.log(F), 1)[0])

    H1_pure = mfdfa_q0(x_h1, scales_long)
    L1_pure = mfdfa_q0(x_l1, scales_long)
    log.info(f"Pure H1_H = {H1_pure:.4f}")
    log.info(f"Pure L1_H = {L1_pure:.4f}")

    # Interpretation
    verdict = ""
    if H_cross > 0.1:
        verdict = "SURPRISING: The control loops / seismic background of H1 and L1 are strongly correlated at long distances."
    elif abs(H_cross) < 0.05:
        verdict = "EXPECTED: H1 and L1 have independent feedback loops / seismic backgrounds with mean-reverting (H~0.04) behavior. No common fractal structure."
    else:
        verdict = "INCONCLUSIVE: Weak but non-zero correlation."

    out = {
        "H1_Pure_H": H1_pure,
        "L1_Pure_H": L1_pure,
        "Cross_H1_L1_Pure_H": H_cross,
        "Interpretation": verdict
    }

    with open("QW_1660_v48_Pure_Raw_CrossDFA.json", "w") as f:
        json.dump(out, f, indent=2)

    log.info("QW-1660 v48 COMPLETE")
    print("\n--- RESULTS ---")
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
