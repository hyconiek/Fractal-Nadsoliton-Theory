import numpy as np
import json
import logging
import multiprocessing

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

def generate_fgn(n, H, seed=None):
    if seed is not None: np.random.seed(seed)
    f = np.fft.fftfreq(n)
    f[0] = 1e-10
    power = np.abs(f) ** (-(2*H - 1))
    power[0] = 0
    phases = np.random.uniform(0, 2*np.pi, n)
    spectrum = np.sqrt(power) * np.exp(1j * phases)
    if n % 2 == 0: spectrum[n//2] = np.abs(spectrum[n//2])
    for i in range(1, n//2): spectrum[n-i] = np.conj(spectrum[i])
    spectrum[0] = 0
    x = np.real(np.fft.ifft(spectrum))
    x = (x - np.mean(x)) / (np.std(x) + 1e-30)
    return x

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

def run_trial(trial_idx):
    # N=131072 (128K) is enough for stable MF-DFA but fast enough for 500 trials
    N = 131072
    H_LOCAL = 0.04
    scales = np.logspace(3, np.log10(N//4), 15).astype(int)
    scales = np.unique(scales)
    
    n1 = generate_fgn(N, H_LOCAL, seed=trial_idx*100 + 50)
    n2 = generate_fgn(N, H_LOCAL, seed=trial_idx*100 + 51)
    
    hc = cross_mfdfa_q0(n1, n2, scales)
    return hc

def main():
    N_TRIALS = 500
    log.info(f"START Phase 57: High-Statistics Monte Carlo Null Test ({N_TRIALS} trials)")
    
    # Run in parallel to save time
    num_cores = multiprocessing.cpu_count()
    log.info(f"Using {num_cores} cores...")
    
    with multiprocessing.Pool(num_cores) as p:
        results = p.map(run_trial, range(N_TRIALS))
    
    results = np.array(results)
    mean_val = np.mean(results)
    std_val = np.std(results)
    min_val = np.min(results)
    percentile_1 = np.percentile(results, 1)
    
    log.info(f"Null Model Mean: {mean_val:.4f} \u00b1 {std_val:.4f}")
    log.info(f"Minimum observed Cross-H: {min_val:.4f}")
    
    # Bootstrap CI for mean
    boot_means = []
    for _ in range(1000):
        sample = np.random.choice(results, size=N_TRIALS, replace=True)
        boot_means.append(np.mean(sample))
    ci_lower = float(np.percentile(boot_means, 2.5))
    ci_upper = float(np.percentile(boot_means, 97.5))
    
    # Calculate sigma exactly if observed is 0.311
    obs_h = 0.311
    sigma_deviation = (mean_val - obs_h) / std_val
    
    out = {
        "Config": {
            "N_trials": N_TRIALS,
            "N_samples_per_trial": 131072,
            "H_local": 0.04
        },
        "Results": {
            "Null_Cross_H_Mean": float(mean_val),
            "Null_Cross_H_Std": float(std_val),
            "Null_Cross_H_Min": float(min_val),
            "Null_1st_Percentile": float(percentile_1),
        },
        "Bootstrap_95_CI": [ci_lower, ci_upper],
        "Real_LIGO_Observation": obs_h,
        "Statistical_Significance": f"{sigma_deviation:.2f} sigma",
        "Verdict": ""
    }
    
    if min_val > obs_h:
        out["Verdict"] = f"Strong anomaly confirmed at {sigma_deviation:.2f} sigma. The observed {obs_h} NEVER occurs in {N_TRIALS} null trials. Real cross-detector anti-persistence exists."
    else:
        out["Verdict"] = f"Anomaly is weaker than thought. {obs_h} occurred in the null model tails."
        
    with open("QW_1660_v57_HighStatsMC.json", "w") as f:
        json.dump(out, f, indent=2)
        
    print("\n--- HIGH STATS MONTE CARLO RESULTS ---")
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
