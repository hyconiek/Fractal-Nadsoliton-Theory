import os, json, logging
import numpy as np
import h5py
import multiprocessing

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

RAW_DIR = "./raw_strain_unfiltered"
GPS = 1266965117 # Baseline
N_SAMPLES = 524288 # Use 128 seconds to capture standard scales

def fetch_empirical_data():
    path_h1 = f"{RAW_DIR}/H1_unfiltered_{GPS}.h5"
    path_l1 = f"{RAW_DIR}/L1_unfiltered_{GPS}.h5"
    
    with h5py.File(path_h1, "r") as f:
        x_h1 = f["strain"][:N_SAMPLES]
    with h5py.File(path_l1, "r") as f:
        x_l1 = f["strain"][:N_SAMPLES]
        
    return x_h1, x_l1

def power_law_spectrum(n, H):
    f = np.fft.rfftfreq(n)
    f[0] = 1e-10
    power = np.abs(f) ** (-(2*H - 1))
    power[0] = 0
    return np.sqrt(power)

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

SCALES = np.unique(np.logspace(3, np.log10(N_SAMPLES//4), 15).astype(int))

def run_test_1(args):
    """Test 1: Independent Envelope Control
    Random independent phases, but preserve shared empirical PSD amplitudes."""
    G_H1, G_L1, seed = args
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

def run_test_2(args):
    """Test 2: Envelope Swap Test
    Keep REAL phases (true cross-correlation remains), but swap L1 amplitude for a generic 1/f noise amplitude."""
    x_h1, x_l1, seed = args
    np.random.seed(seed)
    
    R_h1 = np.fft.rfft(x_h1)
    R_l1 = np.fft.rfft(x_l1)
    
    phase_h1 = np.angle(R_h1)
    phase_l1 = np.angle(R_l1)
    
    G_H1 = np.abs(R_h1)
    
    # Generic 1/f envelope (pink noise, H=0.5 -> beta=1 => power prop to 1/f => amplitude prop to 1/sqrt(f))
    G_generic = power_law_spectrum(N_SAMPLES, H=0.5)
    
    # Scale G_generic to somewhat match variance
    G_generic = G_generic * (np.sum(np.abs(R_l1)) / np.sum(G_generic))
    
    X_h1 = G_H1 * np.exp(1j * phase_h1)       # Real H1
    X_l1 = G_generic * np.exp(1j * phase_l1)  # Real phase, generic amplitude
    
    x_h1_out = np.fft.irfft(X_h1, n=N_SAMPLES)
    x_l1_out = np.fft.irfft(X_l1, n=N_SAMPLES)
    
    return cross_mfdfa_q0(x_h1_out, x_l1_out, SCALES)

def main():
    log.info("START Phase 59: Envelope-Swap Tests")
    x_h1, x_l1 = fetch_empirical_data()
    G_H1 = np.abs(np.fft.rfft(x_h1))
    G_L1 = np.abs(np.fft.rfft(x_l1))
    
    N_TRIALS = 100
    num_cores = multiprocessing.cpu_count()
    
    # TEST 1
    log.info("Running Test 1: Independent Phases + Real Empirical Envelopes...")
    with multiprocessing.Pool(num_cores) as p:
        h_cross_t1 = p.map(run_test_1, [(G_H1, G_L1, i) for i in range(N_TRIALS)])
    
    m_t1 = np.mean(h_cross_t1)
    s_t1 = np.std(h_cross_t1)
    log.info(f"Test 1 Cross-H = {m_t1:.4f} \u00b1 {s_t1:.4f}")
    
    # TEST 2
    log.info("Running Test 2: Real Phases + Swapped Generic L1 Envelope...")
    with multiprocessing.Pool(num_cores) as p:
        # Note Phase 59 actually only requires computing it once since phases are fixed?
        # NO, we can just compute it once because we are using the REAL phases and swapping the generic envelope.
        # Wait, if we just swap the generic envelope, there is no stochasticity unless we jitter something, 
        # but the request said "swap L1 PSD with random 1/f PSD". Let's run it 1 time since it's deterministic.
        hc_t2 = run_test_2((x_h1, x_l1, 42))
        
    log.info(f"Test 2 Cross-H = {hc_t2:.4f}")
    
    out = {
        "Test1_Independent_Phases_Real_Envelope": {
            "mean": m_t1,
            "std": s_t1,
            "verdict": "If ~0.50, then the envelope ALONE does not cause the anomaly. If ~0.31, then the envelope alone is responsible."
        },
        "Test2_Real_Phases_Generic_Envelope": {
            "value": hc_t2,
            "verdict": "If ~0.50, then removing the shared envelope destroys the anomaly. If ~0.31, the anomaly is robustly in the phase."
        }
    }

    with open("QW_1660_v59_EnvelopeSwap.json", "w") as f:
        json.dump(out, f, indent=2)
        
    print("\n--- ENVELOPE SWAP TEST ---")
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
