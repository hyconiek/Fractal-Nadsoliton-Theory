import os, json, logging
import numpy as np
import h5py
from gwpy.timeseries import TimeSeries
import multiprocessing

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

RAW_DIR = "./raw_strain_unfiltered"
FS = 4096
GPS = 1266965117 # Baseline
N_SAMPLES = 131072 # 128K samples

def fetch_empirical_spectra():
    # Fetch 128K samples of real data to get exact empirical PSD envelope
    path_h1 = f"{RAW_DIR}/H1_unfiltered_{GPS}.h5"
    path_l1 = f"{RAW_DIR}/L1_unfiltered_{GPS}.h5"
    
    with h5py.File(path_h1, "r") as f:
        x_h1 = f["strain"][:N_SAMPLES]
    with h5py.File(path_l1, "r") as f:
        x_l1 = f["strain"][:N_SAMPLES]
        
    R_h1 = np.fft.rfft(x_h1)
    R_l1 = np.fft.rfft(x_l1)
    
    return np.abs(R_h1), np.abs(R_l1)

def generate_fgn_spectrum(n, H, seed=None):
    if seed is not None: np.random.seed(seed)
    f = np.fft.rfftfreq(n)
    f[0] = 1e-10
    power = np.abs(f) ** (-(2*H - 1))
    power[0] = 0
    phases = np.random.uniform(0, 2*np.pi, len(f))
    phases[0] = 0.0
    if n % 2 == 0: phases[-1] = 0.0
    spectrum = np.sqrt(power) * np.exp(1j * phases)
    return spectrum

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

# Global empirical spectra for workers
G_H1, G_L1 = fetch_empirical_spectra()
SCALES = np.unique(np.logspace(3, np.log10(N_SAMPLES//4), 15).astype(int))

def run_trial(args):
    h_true, seed = args
    # Generate true global background spectrum
    X_in = generate_fgn_spectrum(N_SAMPLES, h_true, seed=seed)
    
    # Apply empirical transfer functions
    X_out_h1 = G_H1 * np.exp(1j * np.angle(X_in))
    X_out_l1 = G_L1 * np.exp(1j * np.angle(X_in))
    
    x_h1 = np.fft.irfft(X_out_h1, n=N_SAMPLES)
    x_l1 = np.fft.irfft(X_out_l1, n=N_SAMPLES)
    
    return cross_mfdfa_q0(x_h1, x_l1, SCALES)

def main():
    h_true_values = np.arange(0.05, 0.46, 0.01)
    h_true_values = np.round(h_true_values, 2)
    N_TRIALS = 100
    
    log.info(f"START Phase 58: Transformation Mapping. {len(h_true_values)} H_true steps, {N_TRIALS} trials each.")
    
    results = {}
    
    num_cores = multiprocessing.cpu_count()
    log.info(f"Using {num_cores} cores...")
    
    with multiprocessing.Pool(num_cores) as p:
        for h_true in h_true_values:
            log.info(f"Processing H_true = {h_true:.2f}...")
            # Prepare arguments for trials
            args = [(h_true, int(h_true*1000) * N_TRIALS + i) for i in range(N_TRIALS)]
            h_obs_trials = p.map(run_trial, args)
            
            mean_obs = np.mean(h_obs_trials)
            std_obs = np.std(h_obs_trials)
            
            results[f"{h_true:.2f}"] = {
                "mean": float(mean_obs),
                "std": float(std_obs)
            }
            log.info(f"  -> H_obs = {mean_obs:.4f} \u00b1 {std_obs:.4f}")

    with open("QW_1660_v58_HTrue_Mapping.json", "w") as f:
        json.dump(results, f, indent=2)
        
    print("\n--- TRANSFORMATION MAPPING F(H_true) -> H_obs ---")
    for h in h_true_values:
        key = f"{h:.2f}"
        marker = " <=== FIN" if abs(h - 0.23) < 0.005 else ""
        if abs(results[key]["mean"] - 0.31) < 0.015:
            marker += "  [MATCHES 0.31]"
        print(f"H_true = {h:.2f}  =>  H_obs = {results[key]['mean']:.4f} \u00b1 {results[key]['std']:.4f}{marker}")

if __name__ == "__main__":
    main()
