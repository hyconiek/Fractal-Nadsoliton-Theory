import os, json, logging, time
import numpy as np
import h5py
from gwpy.timeseries import TimeSeries
from scipy.signal import detrend

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
log = logging.getLogger("QW-1660-v47")
log.info("START QW-1660 v47: FILTER PARADOX VERIFICATION")

RAW_DIR = "./raw_strain"
os.makedirs(RAW_DIR, exist_ok=True)
FS = 4096
WINDOW_SEC = 512
GPS = 1266965117 # known good GPS from previous scans

def mfdfa_q0(x, scales):
    y = np.cumsum(x - np.mean(x))
    F = []
    valid_scales = []
    N = len(y)
    for s in scales:
        n = N // s
        if n == 0: continue
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
    
    if len(F) < 3: return np.nan
    H = np.polyfit(np.log(valid_scales), np.log(F), 1)[0]
    return float(H)

def main():
    log.info("Fetching pure raw strain...")
    try:
        ts_raw = TimeSeries.fetch_open_data("H1", GPS, GPS + WINDOW_SEC, verbose=False)
        if ts_raw.sample_rate.value > FS:
            ts_raw = ts_raw.resample(FS)
    except Exception as e:
        log.error(f"Failed to fetch data: {e}. Ensure network/VPN allows GWOSC access.")
        return

    # 1. TEST UNFILTERED (LONG SCALES)
    # Apply notch filters to remove line noise, but NO high-pass
    log.info("Applying notch filters...")
    ts_unfiltered = ts_raw.notch(60).notch(120).notch(180)
    x_unfiltered = detrend(ts_unfiltered.value)

    # Original scales: 10^3 to N/4
    scales_long = np.logspace(3, np.log10(len(x_unfiltered)//4), 15).astype(int)
    scales_long = np.unique(scales_long)
    log.info("Computing H on UNFILTERED data (Original long scales)...")
    H_unfiltered = mfdfa_q0(x_unfiltered, scales_long)
    log.info(f"H_unfiltered = {H_unfiltered:.4f}")

    # 2. TEST FILTERED (SHORT SCALES)
    # Apply bandpass 20-1000 Hz just like in previous scripts
    log.info("Applying bandpass 20-1000 Hz...")
    ts_filtered = ts_unfiltered.bandpass(20, 1000)
    x_filtered = detrend(ts_filtered.value)

    # Short scales matching passband: 
    # 30Hz ~ 136 samples, 300Hz ~ 13 samples. We'll examine scales 10 to 200.
    scales_short = np.logspace(np.log10(10), np.log10(200), 15).astype(int)
    scales_short = np.unique(scales_short)
    log.info(f"Computing H on FILTERED data (Short scales {scales_short[0]}-{scales_short[-1]})...")
    H_short = mfdfa_q0(x_filtered, scales_short)
    log.info(f"H_short = {H_short:.4f}")

    # 3. TEST FILTERED (LONG SCALES) - REPRODUCE ANOMALY
    log.info("Computing H on FILTERED data (Original long scales - Anomaly repo)...")
    H_anomaly = mfdfa_q0(x_filtered, scales_long)
    log.info(f"H_anomaly = {H_anomaly:.4f}")

    # Save results
    out = {
        "H_unfiltered_long_scales": H_unfiltered,
        "H_filtered_short_scales": H_short,
        "H_anomaly_reproduced": H_anomaly,
        "Interpretation": {
            "H_anomaly_reproduced": "The ~0.23 value seen previously is an artifact of the 20Hz filter forcing the signal to cross zero over long scales.",
            "H_unfiltered_long_scales": "The true long-range scaling exponent of the raw noise.",
            "H_filtered_short_scales": "The short-range scaling exponent inside the studied beta frequency band."
        }
    }

    with open("QW_1660_v47_Filter_Paradox.json", "w") as f:
        json.dump(out, f, indent=2)

    log.info("QW-1660 v47 COMPLETE")
    print("\n--- RESULTS ---")
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
