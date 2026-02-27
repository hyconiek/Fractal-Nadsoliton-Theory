import os, json, logging
import numpy as np
import h5py
from scipy.signal import detrend
from scipy.stats import linregress

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
    log.info("START Phase 66: Scale-Dependent Correlation Delay")
    x_h1, x_l1 = fetch_empirical_data()
    
    log.info("Whitening empirical data...")
    hw = exact_spectral_whitening(x_h1)
    lw = exact_spectral_whitening(x_l1)
    
    # We define scales (window sizes in samples)
    scales = np.unique(np.logspace(np.log10(128), np.log10(N_SAMPLES//8), 20).astype(int))
    
    max_lag_ms = 1000.0 # Search window +/- 1 second
    max_lag_samples = int((max_lag_ms / 1000.0) * FS)
    
    tau_max_values = []
    valid_scales = []
    peaks_corr = []
    
    log.info(f"Computing scale-dependent cross-correlation over {len(scales)} scales...")
    
    # We will compute correlation for each window scale by averaging the correlation 
    # of many non-overlapping windows of that size.
    for s in scales:
        n_windows = len(hw) // s
        if n_windows == 0: continue
        
        # We need windows that allow shifts up to max_lag_samples
        # So effective window size to compare is slightly smaller, or we just pad
        
        # Let's use a simpler approach: for a given scale S, we slide a window of size S across L1 
        # against a static window of size S in H1, within local bounds of +/- max_lag_samples
        
        # To avoid edge effects, we take pairs of segments from the middle
        # Actually average over several segments
        avg_correlogram = np.zeros(2 * max_lag_samples + 1)
        k_count = 0
        
        lags = np.arange(-max_lag_samples, max_lag_samples + 1)
        
        for i in range(n_windows):
            start = i * s
            end = start + s
            
            # Ensure we can shift L1 by max_lag_samples in both directions
            if start - max_lag_samples < 0 or end + max_lag_samples >= len(lw):
                continue
                
            h_seg = hw[start:end]
            h_seg_norm = (h_seg - np.mean(h_seg)) / (np.std(h_seg) + 1e-12)
            
            # Vectorized fast correlation over all lags simultaneously
            l_segs = np.array([lw[start + lag : end + lag] for lag in lags])
            l_segs_mean = np.mean(l_segs, axis=1, keepdims=True)
            l_segs_std = np.std(l_segs, axis=1, keepdims=True) + 1e-12
            l_segs_norm = (l_segs - l_segs_mean) / l_segs_std
            
            local_corr = np.dot(l_segs_norm, h_seg_norm) / s
            
            avg_correlogram += local_corr
            k_count += 1
            
        if k_count == 0: continue
            
        avg_correlogram /= k_count
        
        # Find tau_max (absolute delay magnitude since we just want tau(s) relation)
        # We might have negative or positive delay. We are interested in the absolute delay magnitude
        # or we just find the peak of the absolute correlation.
        best_idx = np.argmax(np.abs(avg_correlogram))
        best_lag_samples = lags[best_idx]
        best_tau_ms = abs(best_lag_samples) / FS * 1000.0
        
        valid_scales.append(s)
        tau_max_values.append(best_tau_ms)
        peaks_corr.append(float(avg_correlogram[best_idx]))
        
        log.info(f"Scale {s} samples (~{s/FS:.2f}s): tau_max = {best_tau_ms:.2f} ms (corr={avg_correlogram[best_idx]:.5f})")

    valid_scales = np.array(valid_scales)
    tau_max_values = np.array(tau_max_values)
    
    # We need strictly positive tau_max to take log.
    # If tau_max == 0, we can drop it or bound it to a small value e.g. 0.01 ms
    non_zero_idx = tau_max_values > 0
    
    log_s = np.log10(valid_scales[non_zero_idx])
    log_tau = np.log10(tau_max_values[non_zero_idx])
    
    if len(log_s) > 2:
        slope, intercept, r_value, p_value, std_err = linregress(log_s, log_tau)
        log.info(f"Linear Fit log(tau_max) vs log(s): Slope = {slope:.4f} \u00b1 {std_err:.4f}")
    else:
        slope, std_err, r_value, p_value = 0, 0, 0, 0
        log.info("Not enough non-zero tau_max values to perform regression.")

    results_data = {
        "Scales_samples": [int(x) for x in valid_scales],
        "Tau_max_ms": [float(x) for x in tau_max_values],
        "Peak_Correlations": [float(x) for x in peaks_corr],
        "Regression": {
            "Slope_H_obs": float(slope),
            "Std_Err": float(std_err),
            "R_squared": float(r_value**2),
            "P_value": float(p_value)
        },
        "Verdict": ""
    }
    
    if abs(slope - 0.23) < 0.05 and p_value < 0.05:
        results_data["Verdict"] = "FIN CONFIRMED: Deterministic scale-dependent delay strongly matches H=0.23!"
    elif slope < 0.05:
        results_data["Verdict"] = "FIN FALSIFIED: Delay \u03c4_max is effectively independent of scale (slope ~ 0)."
    else:
        results_data["Verdict"] = f"Unclear: Slope is {slope:.4f}, doesn't exactly match FIN but implies some non-trivial scaling."
        
    with open("QW_1660_v66_ScaleDelayFunction.json", "w") as f:
        json.dump(results_data, f, indent=2)

    print("\n--- SCALE-DEPENDENT DELAY FUNCTION RESULTS ---")
    print(json.dumps(results_data["Regression"], indent=2))
    print(f"Verdict: {results_data['Verdict']}")

if __name__ == "__main__":
    main()
