"""
Phase 76: True FIN Functional on Whitened LIGO Strain
Implements the mathematically exact formulation of FIN on H1 and L1 data:
F_FIN[x,y] = H( MF-DFA_q=0 [ C(x,y) ] )
Author: Antigravity & User | Date: Feb 2026
"""

import os
import json
import logging
import numpy as np
import matplotlib.pyplot as plt
from gwpy.timeseries import TimeSeries
from scipy.signal import detrend
from scipy.stats import linregress
from scipy.optimize import curve_fit
import warnings

warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
log = logging.getLogger("Phase76")

GPS = 1266965117
WINDOW_SEC = 512
FS = 4096

def mf_dfa_q0(Z, min_scale=15):
    """
    Computes the log-Poisson q=0 Limit of MF-DFA.
    """
    N = len(Z)
    if N < 4 * min_scale:
        return np.nan, [], []

    # 15 scales logarithmically spaced
    scales = np.unique(np.logspace(np.log10(min_scale), np.log10(N//4), 30).astype(int))
    F_q0 = []
    valid_scales = []
    
    for s in scales:
        k = N // s
        if k < 4: continue
        var_list = []
        
        for i in range(k):
             seg = Z[i*s:(i+1)*s]
             t = np.arange(s)
             pz = np.polyfit(t, seg, 1)
             var = np.mean((seg - np.polyval(pz, t))**2)
             if var > 1e-15:  
                 var_list.append(var)
        
        if len(var_list) > 2:
            F_q0.append(np.exp(0.5 * np.mean(np.log(var_list))))
            valid_scales.append(s)

    valid_scales = np.array(valid_scales)
    F_q0 = np.array(F_q0)
    
    if len(valid_scales) < 3:
        return np.nan, [], []
        
    slope, _, _, _, _ = linregress(np.log(valid_scales), np.log(F_q0))
    return slope, valid_scales, F_q0

def main():
    log.info(f"Fetching H1 and L1 raw strain for GPS {GPS} ({WINDOW_SEC}s)...")
    
    try:
        ts_h1 = TimeSeries.fetch_open_data("H1", GPS, GPS + WINDOW_SEC, verbose=False)
        ts_l1 = TimeSeries.fetch_open_data("L1", GPS, GPS + WINDOW_SEC, verbose=False)
        
        if ts_h1.sample_rate.value > FS: ts_h1 = ts_h1.resample(FS)
        if ts_l1.sample_rate.value > FS: ts_l1 = ts_l1.resample(FS)
    except Exception as e:
        log.error(f"Failed to fetch data: {e}. Check internet/GWOSC.")
        return

    log.info("Whitening data to utterly destroy PSD/amplitude artifacts...")
    # Whitening flattens the power spectral density (PSD), extracting the pure phase!
    ts_h1_white = ts_h1.whiten(fduration=4)
    ts_l1_white = ts_l1.whiten(fduration=4)
    
    # Crop edges corrupted by whitening filter
    crop_sec = 4
    x1 = ts_h1_white.value[int(crop_sec*FS) : -int(crop_sec*FS)]
    x2 = ts_l1_white.value[int(crop_sec*FS) : -int(crop_sec*FS)]
    
    # Bandpass standard LIGO band to clean out DC/Extreme HF noise
    x1 = detrend(x1)
    x2 = detrend(x2)
    
    log.info("Standardizing Whitened Streams (Non-Energetic Axiom A2)...")
    x1 = (x1 - np.mean(x1)) / np.std(x1)
    x2 = (x2 - np.mean(x2)) / np.std(x2)
    
    log.info("Applying Relational Operator (Axiom A1): z = x1 * x2")
    z = x1 * x2
    
    log.info("Generating cumulative profile for MF-DFA(q=0)...")
    Z = np.cumsum(z - np.mean(z))
    
    log.info("Running Log-Poisson True Estimator: H( MF-DFA_{q=0} [ Z ] )")
    H_q0, scales, F_vals = mf_dfa_q0(Z, min_scale=30)
    
    log.info(f">>> TRUE FIN FUNCTIONAL H_{{q=0}} = {H_q0:.4f} <<<")
    
    # Save the output plot
    plt.figure(figsize=(10, 6))
    plt.loglog(scales, F_vals, 'o', label=f'Log-Poisson Fluctuation F_0(s)\nTrue $H_{{q=0}} = {H_q0:.4f}$', color='purple')
    
    # Fit line
    p = np.polyfit(np.log(scales), np.log(F_vals), 1)
    y_fit = np.exp(np.polyval(p, np.log(scales)))
    plt.loglog(scales, y_fit, 'r-', linewidth=2, label=f'Fit (Slope: {p[0]:.4f})')
    
    plt.xlabel('Scale s (samples)')
    plt.ylabel('Fluctuation $F_0(s)$')
    plt.title(f'Phase 76: True FIN Functional\nWhitened Raw Strain (GPS {GPS})')
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.savefig("phase76_true_fin_functional.png", dpi=200, bbox_inches='tight')
    plt.close()

    out = {
        "gps": GPS,
        "duration": WINDOW_SEC,
        "fs": FS,
        "H_q0_whitened": H_q0,
        "interpretation": "If H_q0 > 0.15, FIN survives Whitening and exhibits true phase-space log-Poisson structure. If H_q0 ~ 0, flat noise."
    }
    
    with open("phase76_true_fin_results.json", "w") as f:
        json.dump(out, f, indent=4)
        
    log.info("Phase 76 Completed successfully.")

if __name__ == "__main__":
    main()
