import os, json, logging
import numpy as np
import h5py
from scipy.signal import detrend

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

def exact_spectral_whitening(x):
    """
    Flattens the PSD perfectly to white noise (Gaussian amplitude = 1 for all f),
    but preserves the original phase angles identically.
    """
    X_f = np.fft.rfft(x)
    phases = np.angle(X_f)
    
    # Set amplitude to a constant flat 1.0 (excluding DC and Nyquist which we keep 0/real)
    X_white = np.exp(1j * phases)
    X_white[0] = 0.0 # Remove DC bias
    if len(x) % 2 == 0: X_white[-1] = np.real(X_white[-1])
    
    x_white = np.fft.irfft(X_white, n=len(x))
    return detrend(x_white)

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
    log.info("START Phase 63: Envelope-Invariant Whitening Test")
    x_h1, x_l1 = fetch_empirical_data()
    
    log.info("Applying absolute spectral whitening (preserving only true phase)...")
    hw1 = exact_spectral_whitening(x_h1)
    lw1 = exact_spectral_whitening(x_l1)
    
    scales = np.unique(np.logspace(3, np.log10(N_SAMPLES//4), 15).astype(int))
    
    log.info("Computing Cross-H on purely phase-correlated whitened signal...")
    h_cross = cross_mfdfa_q0(hw1, lw1, scales)
    log.info(f"Whitened Cross-H = {h_cross:.4f}")
    
    out = {
        "Whitened_Cross_H": h_cross,
        "Verdict": "If ~0.50, there is NO phase structure (FIN is falsified as a phase effect). If significantly diff from 0.50, FIN exists as pure phase correlation."
    }
    
    with open("QW_1660_v63_Whitening.json", "w") as f:
        json.dump(out, f, indent=2)

    print("\n--- WHITENING TEST RESULTS ---")
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
