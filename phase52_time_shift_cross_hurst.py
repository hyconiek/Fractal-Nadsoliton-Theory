import os, json, logging
import numpy as np
import h5py
from gwpy.timeseries import TimeSeries
from scipy.signal import detrend

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

RAW_DIR = "./raw_strain_unfiltered"
os.makedirs(RAW_DIR, exist_ok=True)
FS = 4096
WINDOW_SEC = 512
GPS = 1266965117 # Baseline quiet GPS
# Use a bit longer window to allow for time shifting
FETCH_WINDOW = WINDOW_SEC + 120

def fetch_pure_strain(det, gps):
    path = f"{RAW_DIR}/{det}_unfiltered_{gps}_long.h5"
    if os.path.exists(path):
        log.info(f"Loading cached pure {det} strain for {gps} from {path}")
        with h5py.File(path, "r") as f:
            return f["strain"][:]
            
    log.info(f"Fetching pure {det} strain for {gps} (window {FETCH_WINDOW}s)...")
    try:
        ts = TimeSeries.fetch_open_data(det, gps, gps + FETCH_WINDOW, verbose=False, cache=True)
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

def main():
    x_h1_full = fetch_pure_strain("H1", GPS)
    x_l1_full = fetch_pure_strain("L1", GPS)
    
    if x_h1_full is None or x_l1_full is None:
        log.error("Missing data, aborting.")
        return
        
    shifts_sec = [0.0, 0.1, 1.0, 5.0, 10.0, 50.0, 100.0]
    
    base_length = WINDOW_SEC * FS
    scales = np.logspace(3, np.log10(base_length//4), 15).astype(int)
    scales = np.unique(scales)
    
    results = {}
    
    log.info("Starting Time-Shift Decorrelation Test...")
    
    for shift in shifts_sec:
        log.info(f"--- Shift = {shift} s ---")
        shift_samples = int(shift * FS)
        
        # x is stationary at the beginning
        x = x_h1_full[:base_length]
        
        # y is shifted by `shift_samples` forward
        y = x_l1_full[shift_samples:shift_samples + base_length]
        
        x = detrend(x)
        y = detrend(y)
        
        h_cross = cross_mfdfa_q0(x, y, scales)
        log.info(f"Cross-H (tau={shift}s): {h_cross:.4f}")
        
        results[str(shift)] = h_cross
        
    with open("QW_1660_v52_TimeShift.json", "w") as f:
        json.dump(results, f, indent=2)
        
    print("\n--- TIME SHIFT DECORRELATION RESULTS ---")
    print(f"{'Shift (s)':>10} | {'Cross-H':>10}")
    print("-" * 25)
    for shift in shifts_sec:
        print(f"{shift:10.1f} | {results[str(shift)]:10.4f}")

if __name__ == "__main__":
    main()
