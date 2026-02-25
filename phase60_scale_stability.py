import os, json, logging
import numpy as np
import h5py

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

RAW_DIR = "./raw_strain_unfiltered"
GPS = 1266965117
N_SAMPLES = 524288 # 128s

def fetch_empirical_data():
    path_h1 = f"{RAW_DIR}/H1_unfiltered_{GPS}.h5"
    path_l1 = f"{RAW_DIR}/L1_unfiltered_{GPS}.h5"
    with h5py.File(path_h1, "r") as f: x_h1 = f["strain"][:N_SAMPLES]
    with h5py.File(path_l1, "r") as f: x_l1 = f["strain"][:N_SAMPLES]
    return x_h1, x_l1

def cross_mfdfa_q0_local(x, y, scales):
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
            
    if len(valid_scales) < 3: return np.nan
    return float(np.polyfit(np.log(valid_scales), np.log(F), 1)[0])

def main():
    log.info("START Phase 60: Scale Stability Analysis")
    x_h1, x_l1 = fetch_empirical_data()
    
    # We want 5 sequential scale ranges spanning from 1000 to N/4
    # Let's say: [1000, 3000], [3000, 10000], [10000, 30000], [30000, 100000], [100000, N/4]
    
    scale_boundaries = [1000, 3000, 10000, 30000, 100000, N_SAMPLES//4]
    
    results = {}
    
    for i in range(5):
        s_min = scale_boundaries[i]
        s_max = scale_boundaries[i+1]
        
        scales = np.unique(np.logspace(np.log10(s_min), np.log10(s_max), 15).astype(int))
        
        h_cross = cross_mfdfa_q0_local(x_h1, x_l1, scales)
        
        key = f"Scale_{s_min}_to_{s_max}"
        results[key] = {
            "s_min": s_min,
            "s_max": s_max,
            "H_cross": h_cross
        }
        log.info(f"{key}: H_cross = {h_cross:.4f}")
        
    with open("QW_1660_v60_ScaleStability.json", "w") as f:
        json.dump(results, f, indent=2)

    print("\n--- SCALE STABILITY ANALYSIS ---")
    print(json.dumps(results, indent=2))

if __name__ == "__main__":
    main()
