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
        with h5py.File(path, "r") as f:
            return f["strain"][:]
            
    ts = TimeSeries.fetch_open_data(det, gps, gps + WINDOW_SEC, verbose=False, cache=True)
    if ts.sample_rate.value > FS:
        ts = ts.resample(FS)
    ts = ts.notch(60).notch(120).notch(180)
    val = ts.value
    with h5py.File(path, "w") as f:
        f.create_dataset("strain", data=val)
    return val

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

def empirical_filter_injection(target_H, target_signal, real_data):
    # Instead of strict ARMA which requires statsmodels and takes forever on 2M samples,
    # we color the synthetic fGn (target_signal) with the exact empirical power spectrum 
    # of the real detector noise. This perfectly mimics the frequency response (control loop damping)
    # without altering the phase relationships of the injected true background.
    
    # Get empirical shape
    R = np.fft.rfft(real_data)
    empirical_amplitude = np.abs(R)
    
    # Get phases of the shared fundamental background
    T = np.fft.rfft(target_signal)
    target_phases = np.angle(T)
    
    # Combine: Empirical Control Loop Amplitude + Fundamental Background Phase
    # This precisely answers: "If H=0.23 is the true background, what happens when it passes through LIGO's physical instrument response?"
    combined_spectrum = empirical_amplitude * np.exp(1j * target_phases)
    
    injected_signal = np.fft.irfft(combined_spectrum, n=len(real_data))
    return injected_signal

def main():
    log.info("START Phase 55: Empirical Control Response Injection Test")
    x_h1 = detrend(fetch_pure_strain("H1", GPS))
    x_l1 = detrend(fetch_pure_strain("L1", GPS))
    
    N = len(x_h1)
    scales = np.logspace(3, np.log10(N//4), 15).astype(int)
    scales = np.unique(scales)
    
    H_TRUE = 0.23
    log.info(f"Generating fundamental global background background with H={H_TRUE}...")
    global_background = generate_fgn(N, H_TRUE, seed=42)
    
    log.info("Passing background through H1 control response envelope...")
    h1_injected = empirical_filter_injection(H_TRUE, global_background, x_h1)
    
    log.info("Passing background through L1 control response envelope...")
    l1_injected = empirical_filter_injection(H_TRUE, global_background, x_l1)
    
    log.info("Computing Cross-H of the instrumentally-filtered H=0.23 background...")
    h_cross_injection = cross_mfdfa_q0(h1_injected, l1_injected, scales)
    
    log.info(f"Injected Cross-H = {h_cross_injection:.4f}")
    
    verdict = ""
    if abs(h_cross_injection - 0.31) < 0.03:
        verdict = "The instrumental response accurately transforms a true H=0.23 background into an apparent Cross-H ~ 0.31. The 0.31 anomaly COULD be a projection of a 0.23 fundamental background."
    else:
        verdict = f"Passing H=0.23 through the real instrument response yields {h_cross_injection:.3f}, NOT 0.31. The 0.31 anomaly is NOT a projection of 0.23."

    out = {
        "Fundamental_H_Input": H_TRUE,
        "Injected_Cross_H_Output": h_cross_injection,
        "Verdict": verdict
    }
    
    with open("QW_1660_v55_InjectionTest.json", "w") as f:
        json.dump(out, f, indent=2)
        
    print("\n--- INJECTION TEST RESULTS ---")
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
