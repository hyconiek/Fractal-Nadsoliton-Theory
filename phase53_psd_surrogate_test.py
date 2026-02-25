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

def phase_randomize(x, seed=None):
    if seed is not None:
        np.random.seed(seed)
    X = np.fft.rfft(x)
    phases = np.random.uniform(0, 2 * np.pi, len(X))
    # Keep DC component and Nyquist (if even length) phase to 0 or pi for purely real signal, 
    # but strictly speaking rfft handles the conjugate symmetry.
    # We can just randomize all phases since irfft only cares about the positive frequencies.
    # To keep the mean the same, set phase of DC to 0
    phases[0] = 0.0
    if len(x) % 2 == 0:
        phases[-1] = 0.0 # Nyquist
        
    X_surr = np.abs(X) * np.exp(1j * phases)
    x_surr = np.fft.irfft(X_surr, n=len(x))
    return x_surr

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
    x_h1 = fetch_pure_strain("H1", GPS)
    x_l1 = fetch_pure_strain("L1", GPS)
    
    if x_h1 is None or x_l1 is None:
        log.error("Missing data.")
        return
        
    x_h1 = detrend(x_h1)
    x_l1 = detrend(x_l1)
    
    scales_long = np.logspace(3, np.log10(len(x_h1)//4), 15).astype(int)
    scales_long = np.unique(scales_long)
    
    log.info("Computing Original Cross-H...")
    h_cross_orig = cross_mfdfa_q0(x_h1, x_l1, scales_long)
    log.info(f"Original Cross-H: {h_cross_orig:.4f}")
    
    N_SURR = 20
    surr_cross_h = []
    
    log.info("Starting Phase Randomization Surrogates...")
    for i in range(N_SURR):
        x_h1_surr = phase_randomize(x_h1, seed=i*100)
        x_l1_surr = phase_randomize(x_l1, seed=i*100+1)
        
        hc = cross_mfdfa_q0(x_h1_surr, x_l1_surr, scales_long)
        surr_cross_h.append(hc)
        log.info(f"Surrogate {i+1}/{N_SURR} -> Cross-H: {hc:.4f}")
        
    mean_surr = float(np.mean(surr_cross_h))
    std_surr = float(np.std(surr_cross_h))
    
    log.info(f"Surrogate Cross-H = {mean_surr:.4f} \u00b1 {std_surr:.4f}")
    
    out = {
        "original_cross_H": h_cross_orig,
        "surrogate_cross_H_mean": mean_surr,
        "surrogate_cross_H_std": std_surr,
        "N_surrogates": N_SURR
    }
    
    with open("QW_1660_v53_PSDSurrogate.json", "w") as f:
        json.dump(out, f, indent=2)
        
    print("\n--- PSD-MATCHED SURROGATE RESULTS ---")
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
