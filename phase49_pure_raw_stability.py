import numpy as np
import h5py
import os
import logging
import json
from gwpy.timeseries import TimeSeries
from scipy.signal import detrend

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

RAW_DIR = "./raw_strain_unfiltered"
os.makedirs(RAW_DIR, exist_ok=True)
FS = 4096
WINDOW_SEC = 512

def fetch_pure_strain(det, gps):
    path = f"{RAW_DIR}/{det}_unfiltered_{gps}.h5"
    if os.path.exists(path):
        log.info(f"Loading cached pure {det} strain for {gps} from {path}")
        with h5py.File(path, "r") as f:
            return f["strain"][:]
            
    log.info(f"Fetching pure {det} strain for {gps}...")
    try:
        ts = TimeSeries.fetch_open_data(det, gps, gps + WINDOW_SEC, verbose=False, cache=True)
        if ts.sample_rate.value > FS:
            ts = ts.resample(FS)
        
        # APPLY ONLY NOTCH FILTERS (NO BANDPASS)
        ts = ts.notch(60).notch(120).notch(180)
        
        val = ts.value
        with h5py.File(path, "w") as f:
            f.create_dataset("strain", data=val)
        return val
    except Exception as e:
        log.error(f"Failed to fetch {det} data for GPS {gps}: {e}")
        return None

def cross_mfdfa_q0(x, y, scales):
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

def main():
    base_gps = 1266965117 # Baseline
    
    gps_times = [
        base_gps,
        base_gps + 86400,          # +1 dzien (Ziemia przekrecona, inne wplywy, ksiezyc)
        base_gps + 86400 * 2,      # +2 dni
        base_gps + 86400 * 7,      # +1 tydzien
        1253326744                 # Random inny czas O3 (~160 dni wczesniej)
    ]
    
    results = {}
    
    for gps in gps_times:
        log.info(f"--- Processing GPS {gps} ---")
        x_h1 = fetch_pure_strain("H1", gps)
        x_l1 = fetch_pure_strain("L1", gps)
        
        if x_h1 is None or x_l1 is None:
            log.warning(f"Skipping GPS {gps} due to missing data.")
            continue
            
        x_h1 = detrend(x_h1)
        x_l1 = detrend(x_l1)
        
        # Test tylko na dugich skalach, by zlapac przestrzen pomiedzy urzadzeniami.
        scales_long = np.logspace(3, np.log10(len(x_h1)//4), 15).astype(int)
        scales_long = np.unique(scales_long)
        
        h_cross = cross_mfdfa_q0(x_h1, x_l1, scales_long)
        h1_pure = mfdfa_q0(x_h1, scales_long)
        l1_pure = mfdfa_q0(x_l1, scales_long)
        
        results[str(gps)] = {
            "H1_Pure": h1_pure,
            "L1_Pure": l1_pure,
            "Cross_H1_L1": h_cross
        }
        
        log.info(f"GPS {gps} -> Cross-H: {h_cross:.4f}, H1: {h1_pure:.4f}, L1: {l1_pure:.4f}")
        
    with open("QW_1660_v49_CrossHurst_Stability.json", "w") as f:
        json.dump(results, f, indent=2)
        
    print("\n--- STABILITY RESULTS ---")
    print(json.dumps(results, indent=2))

if __name__ == "__main__":
    main()
