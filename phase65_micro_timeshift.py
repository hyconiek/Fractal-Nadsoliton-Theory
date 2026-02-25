import os, json, logging
import numpy as np
import h5py
from scipy.signal import correlate, detrend

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

RAW_DIR = "./raw_strain_unfiltered"
GPS = 1266965117
FS = 4096
N_SAMPLES = 524288 # 128s

def fetch_empirical_data():
    path_h1 = f"{RAW_DIR}/H1_unfiltered_{GPS}.h5"
    path_l1 = f"{RAW_DIR}/L1_unfiltered_{GPS}.h5"
    with h5py.File(path_h1, "r") as f: x_h1 = f["strain"][:N_SAMPLES]
    with h5py.File(path_l1, "r") as f: x_l1 = f["strain"][:N_SAMPLES]
    return detrend(x_h1), detrend(x_l1)

def exact_spectral_whitening(x):
    X_f = np.fft.rfft(x)
    phases = np.angle(X_f)
    X_white = np.exp(1j * phases)
    X_white[0] = 0.0
    if len(x) % 2 == 0: X_white[-1] = np.real(X_white[-1])
    x_white = np.fft.irfft(X_white, n=len(x))
    return detrend(x_white)

def main():
    log.info("START Phase 65: Precision Time-Shift (Light Travel Time)")
    x_h1, x_l1 = fetch_empirical_data()
    
    # We will compute the cross-correlation function. Since the raw data is heavily
    # dominated by the low-frequency PSD envelope, standard cross-correlation will
    # just show a giant broad peak at 0 or drift. To find a true 10ms (10 millisecond)
    # physical travel time, we MUST whiten the data first, otherwise the envelope obscures it.
    
    log.info("Whitening data for high-resolution timing...")
    hw = exact_spectral_whitening(x_h1)
    lw = exact_spectral_whitening(x_l1)
    
    # We want to check delays from -100 ms to +100 ms
    # 100 ms at 4096 Hz = 0.1 * 4096 = 409.6 samples
    max_lag_samples = 410
    
    log.info("Computing cross-correlation...")
    # Compute correlation only for the small lag window to save memory
    # correlate(mode='valid') is tricky with same length, so we will manually slice
    
    # Center chunk of L1
    center_start = max_lag_samples
    center_end = len(lw) - max_lag_samples
    lw_center = lw[center_start:center_end]
    
    lags = np.arange(-max_lag_samples, max_lag_samples + 1)
    corr_values = []
    
    for lag in lags:
        # lag > 0 means H1 is shifted forward relative to L1
        h1_shifted = hw[center_start + lag : center_end + lag]
        # Pearson correlation
        c = np.corrcoef(h1_shifted, lw_center)[0, 1]
        corr_values.append(c)
        
    corr_values = np.array(corr_values)
    
    # Find the peak
    best_lag_idx = np.argmax(np.abs(corr_values))
    best_lag_samples = lags[best_lag_idx]
    best_lag_ms = (best_lag_samples / FS) * 1000.0
    best_corr = corr_values[best_lag_idx]
    
    # Light travel distance Hanford-Livingston is ~3000 km -> ~10 ms
    log.info(f"Max correlation found at {best_lag_ms:.2f} ms with value {best_corr:.5f}")
    
    # Save a window around [-20ms, +20ms]
    results = {}
    for lag, c in zip(lags, corr_values):
        t_ms = (lag / FS) * 1000.0
        if -25.0 <= t_ms <= 25.0:
            results[f"{t_ms:.2f} ms"] = float(c)
            
    out = {
        "Best_Lag_ms": best_lag_ms,
        "Best_Correlation": best_corr,
        "Verdict": "If peak is at ~10ms, it is a true gravitational/light-speed wave. If peak is exactly 0ms, it's global synchronous noise. If no clear peak, there is no broadband phase correlation.",
        "Correlation_Window_ms": results
    }

    with open("QW_1660_v65_MicroTimeShift.json", "w") as f:
        json.dump(out, f, indent=2)

    print("\n--- PRECISION TIME-SHIFT RESULTS ---")
    print(f"Max correlation found at {best_lag_ms:.2f} ms with value {best_corr:.5f}")
    if abs(best_lag_ms) > 1.0:
        print(f"Verdict: {out['Verdict']}")

if __name__ == "__main__":
    main()
