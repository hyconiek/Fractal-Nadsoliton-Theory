import os, json, logging
import numpy as np
import h5py
from scipy.signal import butter, filtfilt, detrend

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

RAW_DIR = "./raw_strain_unfiltered"
GPS = 1266965117 
N_SAMPLES = 524288 # 128s
FS = 4096

def fetch_empirical_data():
    path_h1 = f"{RAW_DIR}/H1_unfiltered_{GPS}.h5"
    path_l1 = f"{RAW_DIR}/L1_unfiltered_{GPS}.h5"
    with h5py.File(path_h1, "r") as f: x_h1 = f["strain"][:N_SAMPLES]
    with h5py.File(path_l1, "r") as f: x_l1 = f["strain"][:N_SAMPLES]
    return x_h1, x_l1

def bandpass_filter(data, lowcut, highcut, fs, order=4):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return filtfilt(b, a, data)

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
    log.info("START Phase 64: Band-Pass Scale Test")
    x_h1, x_l1 = fetch_empirical_data()
    
    bands = [
        (40, 80),
        (100, 200),
        (300, 500)
    ]
    
    scales = np.unique(np.logspace(3, np.log10(N_SAMPLES//4), 15).astype(int))
    
    results = {}
    
    for (low, high) in bands:
        log.info(f"Filtering {low}-{high} Hz...")
        h1_bp = bandpass_filter(x_h1, low, high, FS)
        l1_bp = bandpass_filter(x_l1, low, high, FS)
        
        hc = cross_mfdfa_q0(h1_bp, l1_bp, scales)
        log.info(f"Band {low}-{high} Hz Cross-H = {hc:.4f}")
        
        results[f"{low}-{high}Hz"] = hc
        
    with open("QW_1660_v64_Bandpass.json", "w") as f:
        json.dump(results, f, indent=2)

    print("\n--- BANDPASS TEST RESULTS ---")
    print(json.dumps(results, indent=2))

if __name__ == "__main__":
    main()
