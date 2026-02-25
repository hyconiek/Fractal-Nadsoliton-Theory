import os, json, logging
import numpy as np
import h5py
from gwpy.timeseries import TimeSeries
from scipy.signal import detrend, csd, coherence

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

RAW_DIR = "./raw_strain_unfiltered"
os.makedirs(RAW_DIR, exist_ok=True)
FS_TARGET = 1 # 1 Hz
WINDOW_SEC = 4096 # 4096 seconds to allow up to 1000s scales
GPS = 1266965117 # Baseline

def fetch_and_downsample(det, gps):
    path = f"{RAW_DIR}/{det}_downsampled_1Hz_{gps}.h5"
    if os.path.exists(path):
        log.info(f"Loading cached {det} 1Hz strain for {gps} from {path}")
        with h5py.File(path, "r") as f:
            return f["strain"][:]
            
    log.info(f"Fetching {det} strain for {WINDOW_SEC}s (may take time)...")
    try:
        ts = TimeSeries.fetch_open_data(det, gps, gps + WINDOW_SEC, verbose=False, cache=True)
        
        # Antialiasing: The resample method in gwpy uses scipy resample, 
        # actually scipy.signal.decimate is better, but let's downsample
        log.info(f"Resampling {det} to {FS_TARGET} Hz...")
        ts_1hz = ts.resample(FS_TARGET)
        
        log.info(f"Lowpass filtering {det} below 0.1 Hz...")
        # A simple lowpass
        ts_01hz = ts_1hz.lowpass(0.1)
        
        val = ts_01hz.value
        with h5py.File(path, "w") as f:
            f.create_dataset("strain", data=val)
        return val
    except Exception as e:
        log.error(f"Failed to fetch {det} data: {e}")
        return None

def cross_mfdfa_q0(x, y, scales):
    z = np.cumsum((x - np.mean(x)) * (y - np.mean(y)))
    F = []
    valid_scales = []
    N = len(z)
    for s in scales:
        n = N // s
        if n == 0: continue
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
    
    if len(F) < 3: return np.nan
    H = np.polyfit(np.log(valid_scales), np.log(F), 1)[0]
    return float(H)

def main():
    x_h1 = fetch_and_downsample("H1", GPS)
    x_l1 = fetch_and_downsample("L1", GPS)
    
    if x_h1 is None or x_l1 is None:
        log.error("Missing data.")
        return
        
    x_h1 = detrend(x_h1)
    x_l1 = detrend(x_l1)
    
    # 1. MF-DFA on scales 10 to 1000s
    # target fs is 1 Hz, so 10 to 1000 samples
    scales = np.logspace(1, 3, 15).astype(int)
    scales = np.unique(scales)
    
    log.info("Computing Ultra-Low-Freq Cross-H...")
    h_cross = cross_mfdfa_q0(x_h1, x_l1, scales)
    log.info(f"Ultra-Low-Freq Cross-H (10-1000s): {h_cross:.4f}")
    
    # 2. CSD and Coherence
    log.info("Computing Cross Spectral Density (CSD) and Coherence...")
    # Using Welch's method with segments of e.g. 512 samples (~512 seconds)
    nperseg = 512
    f_csd, Pxy = csd(x_h1, x_l1, fs=FS_TARGET, nperseg=nperseg)
    f_coh, Cxy = coherence(x_h1, x_l1, fs=FS_TARGET, nperseg=nperseg)
    
    # Integrate coherence in very low frequencies (0.01 to 0.1 Hz)
    idx_low = (f_coh >= 0.01) & (f_coh <= 0.1)
    mean_coh_low = np.mean(Cxy[idx_low])
    
    # Also evaluate phase relation: if Re(CSD) is systematically negative, we have energetic anti-correlation
    Pxy_low = Pxy[idx_low]
    mean_real_csd = np.mean(np.real(Pxy_low))
    mean_imag_csd = np.mean(np.imag(Pxy_low))
    
    out = {
        "Config": {
            "GPS": GPS,
            "Fs_target": FS_TARGET,
            "Cutoff_Hz": 0.1,
            "Window_sec": WINDOW_SEC
        },
        "Cross_H_LF": h_cross,
        "Coherence_0_01_to_0_1_Hz": float(mean_coh_low),
        "Mean_Real_CSD_LF": float(mean_real_csd),
        "Mean_Imag_CSD_LF": float(mean_imag_csd),
        "Verdict": ""
    }
    
    if h_cross < 0.35 and mean_coh_low < 0.1:
        out["Verdict"] = "Cross-H shows anti-persistence but Coherence is very low. This points to a non-linear or phase-only statistical effect without true energetic correlation."
    elif h_cross < 0.35 and mean_coh_low > 0.3:
        out["Verdict"] = "Cross-H shows anti-persistence AND there is significant energetic coherence. Points to physical/instrumental correlation."
    else:
        out["Verdict"] = "Effect disappears at ultra-low frequencies or interpretation is mixed."

    with open("QW_1660_v54_LowFreq_CSD.json", "w") as f:
        json.dump(out, f, indent=2)

    print("\n--- LOW FREQ (<0.1 Hz) ISOLATION & CSD ---")
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
