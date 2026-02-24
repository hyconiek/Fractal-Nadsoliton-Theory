"""
Phase 50: Monte Carlo Test — Does Cross-MF-DFA inflate the true Hurst exponent?

Hypothesis: If the true background has H=0.23, and we measure it through
two locally-damped channels (H_local ~ 0.04), does Cross-MF-DFA output
~0.31 (observed) or ~0.23 (true)?

Method:
1. Generate fractional Gaussian noise (fGn) with exact H=0.23 (the "shared background")
2. Generate two independent fGn signals with H=0.04 (local damping noise)
3. Mix: x = background + amplitude * local_noise_1
         y = background + amplitude * local_noise_2
4. Measure H(x), H(y), and Cross-H(x,y) via MF-DFA q=0
5. Repeat N_trials times and report statistics

If Cross-H ~ 0.23: the 0.31 we measured in real data is a DIFFERENT phenomenon
If Cross-H ~ 0.31: there exists a systematic inflation from the measurement method
"""

import numpy as np
import json
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

def generate_fgn(n, H, seed=None):
    """Generate fractional Gaussian noise with Hurst exponent H using Hosking method."""
    if seed is not None:
        np.random.seed(seed)
    
    # Use spectral method (Davies-Harte) for speed
    # Power spectrum of fGn: S(f) ~ |f|^(-(2H-1)) for fBm derivative
    # We generate fBm first, then differentiate
    
    # Simple spectral synthesis
    N = n
    f = np.fft.fftfreq(N)
    f[0] = 1e-10  # avoid division by zero
    
    # Power spectrum of fGn: S(f) ~ |f|^(1-2H) ... but we want fGn directly
    # fGn spectrum: S(f) ~ |f|^(-(2H-1)) for 0 < H < 1
    # Actually for fGn: S(f) ~ |2*sin(pi*f)|^2 * |f|^(-(2H+1)) but simplified:
    power = np.abs(f) ** (-(2*H - 1))
    power[0] = 0
    
    # Random phases
    phases = np.random.uniform(0, 2*np.pi, N)
    spectrum = np.sqrt(power) * np.exp(1j * phases)
    
    # Enforce conjugate symmetry for real output
    if N % 2 == 0:
        spectrum[N//2] = np.abs(spectrum[N//2])
    for i in range(1, N//2):
        spectrum[N-i] = np.conj(spectrum[i])
    spectrum[0] = 0
    
    x = np.real(np.fft.ifft(spectrum))
    x = (x - np.mean(x)) / (np.std(x) + 1e-30)
    return x

def mfdfa_q0(x, scales):
    """Standard MF-DFA at q=0."""
    y = np.cumsum(x - np.mean(x))
    F = []
    valid_scales = []
    N = len(y)
    for s in scales:
        n = N // s
        if n == 0:
            continue
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
    if len(F) < 3:
        return np.nan
    return float(np.polyfit(np.log(valid_scales), np.log(F), 1)[0])

def cross_mfdfa_q0(x, y, scales):
    """Cross-MF-DFA at q=0."""
    z = np.cumsum((x - np.mean(x)) * (y - np.mean(y)))
    F = []
    valid_scales = []
    N = len(z)
    for s in scales:
        n = N // s
        if n == 0:
            continue
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
    if len(F) < 3:
        return np.nan
    return float(np.polyfit(np.log(valid_scales), np.log(F), 1)[0])

def main():
    N = 2**20  # ~1M samples (comparable to 512s * 4096 Hz in real data)
    N_TRIALS = 20
    H_TRUE = 0.23       # Hypothesized true background Hurst
    H_LOCAL = 0.04       # Local damping Hurst (measured in real data)
    SNR_LOCAL = 20.0     # Local noise dominates background by factor 20

    scales = np.logspace(3, np.log10(N//4), 15).astype(int)
    scales = np.unique(scales)

    log.info(f"Monte Carlo Cross-Hurst Inflation Test")
    log.info(f"N={N}, trials={N_TRIALS}, H_true={H_TRUE}, H_local={H_LOCAL}, SNR_local={SNR_LOCAL}")

    results_Hx = []
    results_Hy = []
    results_Hcross = []

    # Also test: what does Cross-H give for two INDEPENDENT H=0.04 signals? (null model)
    results_null_cross = []

    for trial in range(N_TRIALS):
        log.info(f"Trial {trial+1}/{N_TRIALS}")
        
        # Generate shared background with H=0.23
        bg = generate_fgn(N, H_TRUE, seed=trial*100)
        
        # Generate independent local noise with H=0.04
        noise1 = generate_fgn(N, H_LOCAL, seed=trial*100 + 1)
        noise2 = generate_fgn(N, H_LOCAL, seed=trial*100 + 2)
        
        # Mix: local noise dominates
        x = bg + SNR_LOCAL * noise1
        y = bg + SNR_LOCAL * noise2
        
        Hx = mfdfa_q0(x, scales)
        Hy = mfdfa_q0(y, scales)
        Hcross = cross_mfdfa_q0(x, y, scales)
        
        results_Hx.append(Hx)
        results_Hy.append(Hy)
        results_Hcross.append(Hcross)
        
        log.info(f"  H(x)={Hx:.4f}, H(y)={Hy:.4f}, Cross-H={Hcross:.4f}")

        # Null model: two purely independent H=0.04 signals (no shared background)
        null1 = generate_fgn(N, H_LOCAL, seed=trial*100 + 50)
        null2 = generate_fgn(N, H_LOCAL, seed=trial*100 + 51)
        Hcross_null = cross_mfdfa_q0(null1, null2, scales)
        results_null_cross.append(Hcross_null)
        log.info(f"  Null Cross-H={Hcross_null:.4f}")

    out = {
        "config": {
            "N": N,
            "N_trials": N_TRIALS,
            "H_true_background": H_TRUE,
            "H_local_damping": H_LOCAL,
            "SNR_local_over_background": SNR_LOCAL,
        },
        "mixed_signal": {
            "H_x_mean": float(np.mean(results_Hx)),
            "H_x_std": float(np.std(results_Hx)),
            "H_y_mean": float(np.mean(results_Hy)),
            "H_y_std": float(np.std(results_Hy)),
            "Cross_H_mean": float(np.mean(results_Hcross)),
            "Cross_H_std": float(np.std(results_Hcross)),
        },
        "null_model": {
            "Cross_H_null_mean": float(np.mean(results_null_cross)),
            "Cross_H_null_std": float(np.std(results_null_cross)),
        },
        "verdict": ""
    }

    cross_mean = out["mixed_signal"]["Cross_H_mean"]
    null_mean = out["null_model"]["Cross_H_null_mean"]

    if abs(cross_mean - 0.23) < 0.03:
        out["verdict"] = "Cross-MF-DFA recovers true H=0.23. The measured 0.31 in real data is a DIFFERENT, independent phenomenon."
    elif abs(cross_mean - 0.31) < 0.03:
        out["verdict"] = "Cross-MF-DFA inflates H=0.23 to ~0.31 through damping. The 0.31 measured in real data MAY be the Weinberg angle seen through instrument response."
    else:
        out["verdict"] = f"Cross-MF-DFA gives {cross_mean:.3f}, which matches neither 0.23 nor 0.31. Relationship is complex."
    
    with open("QW_1660_v50_MonteCarlo_CrossHurst.json", "w") as f:
        json.dump(out, f, indent=2)

    print("\n--- MONTE CARLO RESULTS ---")
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
