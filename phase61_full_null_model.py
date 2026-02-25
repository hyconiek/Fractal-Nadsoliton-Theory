import os, json, logging
import numpy as np
import h5py
import multiprocessing

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

RAW_DIR = "./raw_strain_unfiltered"
GPS = 1266965117
N_SAMPLES = 131072

def fetch_empirical_spectra():
    path_h1 = f"{RAW_DIR}/H1_unfiltered_{GPS}.h5"
    path_l1 = f"{RAW_DIR}/L1_unfiltered_{GPS}.h5"
    with h5py.File(path_h1, "r") as f: x_h1 = f["strain"][:N_SAMPLES]
    with h5py.File(path_l1, "r") as f: x_l1 = f["strain"][:N_SAMPLES]
    R_h1 = np.fft.rfft(x_h1)
    R_l1 = np.fft.rfft(x_l1)
    return np.abs(R_h1), np.abs(R_l1)

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

G_H1, G_L1 = fetch_empirical_spectra()
SCALES = np.unique(np.logspace(3, np.log10(N_SAMPLES//4), 15).astype(int))

def run_trial(seed):
    np.random.seed(seed)
    n = len(G_H1)
    
    phase_h1 = np.random.uniform(0, 2*np.pi, n)
    phase_l1 = np.random.uniform(0, 2*np.pi, n)
    phase_h1[0], phase_l1[0] = 0, 0
    if (N_SAMPLES % 2) == 0: phase_h1[-1], phase_l1[-1] = 0, 0
    
    X_h1 = G_H1 * np.exp(1j * phase_h1)
    X_l1 = G_L1 * np.exp(1j * phase_l1)
    
    x_h1 = np.fft.irfft(X_h1, n=N_SAMPLES)
    x_l1 = np.fft.irfft(X_l1, n=N_SAMPLES)
    
    return cross_mfdfa_q0(x_h1, x_l1, SCALES)

def main():
    N_TRIALS = 1000
    log.info(f"START Phase 61: Full Realistic Null Model ({N_TRIALS} trials)")
    
    num_cores = multiprocessing.cpu_count()
    with multiprocessing.Pool(num_cores) as p:
        results = p.map(run_trial, range(N_TRIALS))
        
    results = np.array(results)
    
    m_null = np.mean(results)
    s_null = np.std(results)
    min_null = np.min(results)
    p1 = np.percentile(results, 1)
    
    # The real actual observation from LIGO data baseline length
    obs_h = 0.311
    
    z_score = (obs_h - m_null) / s_null
    
    log.info(f"Realistic Null Mean: {m_null:.4f} \u00b1 {s_null:.4f}")
    log.info(f"Minimum observed H_cross in 1000 trials: {min_null:.4f}")
    log.info(f"Z-Score against realistic null: {z_score:.2f} sigma")
    
    out = {
        "Config": {
            "N_trials": N_TRIALS,
            "N_samples": N_SAMPLES,
            "Description": "Independent random phases applied strictly to empirical H1 and L1 PSDs"
        },
        "Results": {
            "Mean_H_cross": float(m_null),
            "Std_H_cross": float(s_null),
            "Min_H_cross": float(min_null),
            "P1_H_cross": float(p1)
        },
        "Observation": obs_h,
        "Significance": f"{z_score:.2f} sigma",
        "Verdict": ""
    }
    
    if abs(z_score) < 3.0:
        out["Verdict"] = f"The {obs_h} anomaly is completely destroyed! It falls well within the {z_score:.2f} sigma bound of the purely realistic instrumental envelope model. The effect IS emergent from the PSD architecture."
    else:
        out["Verdict"] = f"Anomaly survives at {z_score:.2f} sigma even against realistic instrumental envelopes! The FIN correlation is still required."
        
    with open("QW_1660_v61_FullNullModel.json", "w") as f:
        json.dump(out, f, indent=2)

    print("\n--- FULL REALISTIC NULL MODEL ---")
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
