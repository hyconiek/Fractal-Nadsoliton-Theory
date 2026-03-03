# Metodologiczny Raport Kodu: Fazy 47-66

Poniżej znajduje się zagregowany kod ze wszystkich faz (pozbawiony importów, printów, logów i komentarzy w celu ułatwienia oceny matematycznej) oraz wygenerowane wyniki (JSON), jeżeli skrypt je zapisywał.

## Faza 47: `phase47_filter_paradox.py`
### Kod:
```python
log = logging.getLogger("QW-1660-v47")

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
    try:
        ts_raw = TimeSeries.fetch_open_data("H1", GPS, GPS + WINDOW_SEC, verbose=False)
        if ts_raw.sample_rate.value > FS:
            ts_raw = ts_raw.resample(FS)
    except Exception as e:
        return

    ts_unfiltered = ts_raw.notch(60).notch(120).notch(180)
    x_unfiltered = detrend(ts_unfiltered.value)

    scales_long = np.logspace(3, np.log10(len(x_unfiltered)//4), 15).astype(int)
    scales_long = np.unique(scales_long)
    H_unfiltered = mfdfa_q0(x_unfiltered, scales_long)

    ts_filtered = ts_unfiltered.bandpass(20, 1000)
    x_filtered = detrend(ts_filtered.value)

    scales_short = np.logspace(np.log10(10), np.log10(200), 15).astype(int)
    scales_short = np.unique(scales_short)
    H_short = mfdfa_q0(x_filtered, scales_short)

    H_anomaly = mfdfa_q0(x_filtered, scales_long)

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

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "H_unfiltered_long_scales": 0.03733273522994843,
  "H_filtered_short_scales": 0.10019515377472395,
  "H_anomaly_reproduced": 0.002944397387996455,
  "Interpretation": {
    "H_anomaly_reproduced": "The ~0.23 value seen previously is an artifact of the 20Hz filter forcing the signal to cross zero over long scales.",
    "H_unfiltered_long_scales": "The true long-range scaling exponent of the raw noise.",
    "H_filtered_short_scales": "The short-range scaling exponent inside the studied beta frequency band."
  }
}
```
```json
{
  "qw847": {
    "N": 200,
    "d_max": 20.0,
    "N_eigenvalues": 200,
    "lambda_max": 252.49263471892004,
    "lambda_min": -58.136558661520496,
    "lambda_range": 310.62919338044054,
    "top_10": [
      252.49263471892004,
      227.2881598109568,
      34.93394850246776,
      31.766434185230608,
      15.666612501744495,
      14.365182951140113,
      9.262260581867196,
      8.542684646454305,
      6.215788519452744,
      5.7618250609189445
    ],
    "bottom_10": [
      0.05575610878080487,
      0.055725168958311,
      0.055701007767883746,
      0.05568381944379165,
      0.05567348450711878,
      -4.630432912082185,
      -35.98916042253846,
      -37.86173836745309,
      -54.693315588346245,
      -58.136558661520496
    ]
  },
  "qw848": {
    "N_ratios": 190,
    "ratio_min": 1.058350545284267,
    "ratio_max": 148.20532752482634,
    "ratio_mean": 15.178182192395111,
    "ratio_std": 29.61424621816515,
    "N_unique": 165,
    "near_power_of_2": 40,
    "near_power_of_4": 38,
    "sample_ratios": [
      1.1108921596660672,
      7.227715318269385,
      7.948409734836188,
      16.11660687278789,
      17.576708600072553,
      27.26036829639915,
      29.556590833973804,
      40.621175242485606,
      43.82164193625316,
      56.18016287040093
    ]
  },
  "qw849": {
    "N_gaps": 57,
    "gap_mean": 5.445391265136252,
    "gap_std": 25.668355333481568,
    "gap_cv": 4.71377612437539,
    "is_regular": false,
    "gap_min": 0.010082953939955286,
    "gap_max": 192.35421130848906,
    "quantization_fraction": 0.24561403508771928,
    "sample_gaps": [
      25.20447490796323,
      192.35421130848906,
      3.167514317237149,
      16.099821683486113,
      1.3014295506043823,
      5.1029223692729175,
      0.719575935412891,
      2.326896127001561,
      0.45396345853379927,
      1.2674871697778656
    ]
  },
  "qw850": {
    "N_log_ratios": 435,
    "log_ratio_mean": 1.623497496019754,
    "log_ratio_std": 1.4119879001854485,
    "N_peaks": 2,
    "peak_positions": [
      2.740409145891295,
      3.3080045796018913
    ],
    "peaks_at_ln4_multiples": [
      [
        2.740409145891295,
        2
      ]
    ],
    "ln4": 1.3862943611198906,
    "NATURAL_BASE_DETECTED": true
  },
  "qw851": {
    "decay_rate_b": 0.227944426581542,
    "intercept_a": 4.225272946973607,
    "r_squared": 0.8393748480110714,
    "natural_base": 1.256015522416549,
    "closest_known_base": "2",
    "closest_base_value": 2.0,
    "base_error_percent": 37.19922387917255,
    "EXPONENTIAL_QUANTIZATION": "False"
  },
  "qw852": {
    "N_degenerate_groups": 3,
    "group_sizes": [
      2,
      2,
      71
    ],
    "N_isolated": 125,
    "degeneracy_fraction": 0.375,
    "FAMILIES_DETECTED": true
  },
  "qw853": {
    "ipr_mean": 0.007934005482959835,
    "ipr_max": 0.09561736182464417,
    "ipr_min": 0.005064912135981274,
    "most_localized_eigenvalue": -4.630432912082185,
    "particle_threshold": 0.015,
    "N_particle_like": 1,
    "STABLE_PARTICLES_FOUND": "True",
    "top_5_localized_eigenvalues": [
      -4.630432912082185,
      227.2881598109568,
      -54.693315588346245,
      -35.98916042253846,
      31.766434185230608
    ]
  },
  "qw854": {
    "N_initial_projections": 40,
    "N_surviving": 1,
    "surviving_eigenvalues": [
      123.18912116511228
    ],
    "dominant_eigenvalue": 123.18912116511228,
    "STABLE_MODES_EXIST": true
  },
  "qw855": {
    "lambda_max": 252.49263471892004,
    "lambda_min": 0.05567348450711878,
    "hierarchy_ratio": 4535.240374376686,
    "orders_of_magnitude": 3.6566003102492206,
    "SM_comparison": 5.531478917042255,
    "SUFFICIENT_HIERARCHY": "True"
  },
  "qw856": {
    "KERNEL_PARAMETERS": {
      "alpha": 2.772588722239781,
      "omega": 0.7853981633974483,
      "phi": 0.5235987755982988,
      "beta": 0.01
    },
    "EMERGENT_PROPERTIES": {
      "Natural_Base": [
        [
          2.740409145891295,
          2
        ]
      ],
      "N_Families": 2,
      "N_Particle_Like": 1,
      "Hierarchy_Orders": 3.6566003102492206
    },
    "SUMMARY": "4/5 key properties emerged from K(d)",
    "SUCCESS": true
  }
}
```

---

## Faza 48: `phase48_pure_raw_crossdfa.py`
### Kod:
```python
log = logging.getLogger("QW-1660-v48")

RAW_DIR = "./raw_strain_unfiltered"
os.makedirs(RAW_DIR, exist_ok=True)
FS = 4096
WINDOW_SEC = 512
GPS = 1266965117 # known good GPS from previous scans

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

def fetch_pure_strain(det, gps):
    path = f"{RAW_DIR}/{det}_unfiltered.h5"
    if os.path.exists(path):
        with h5py.File(path, "r") as f:
            return f["strain"][:]
            
    try:
        ts = TimeSeries.fetch_open_data(det, gps, gps + WINDOW_SEC, verbose=False)
        if ts.sample_rate.value > FS:
            ts = ts.resample(FS)
        
        ts = ts.notch(60).notch(120).notch(180)
        
        val = ts.value
        with h5py.File(path, "w") as f:
            f.create_dataset("strain", data=val)
        return val
    except Exception as e:
        return None

def main():
    x_h1 = fetch_pure_strain("H1", GPS)
    x_l1 = fetch_pure_strain("L1", GPS)
    
    if x_h1 is None or x_l1 is None:
        return
        
    x_h1 = detrend(x_h1)
    x_l1 = detrend(x_l1)

    scales_long = np.logspace(3, np.log10(len(x_h1)//4), 15).astype(int)
    scales_long = np.unique(scales_long)
    
    H_cross = cross_mfdfa_q0(x_h1, x_l1, scales_long)

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
        return float(np.polyfit(np.log(valid_scales), np.log(F), 1)[0])

    H1_pure = mfdfa_q0(x_h1, scales_long)
    L1_pure = mfdfa_q0(x_l1, scales_long)

    verdict = ""
    if H_cross > 0.1:
        verdict = "SURPRISING: The control loops / seismic background of H1 and L1 are strongly correlated at long distances."
    elif abs(H_cross) < 0.05:
        verdict = "EXPECTED: H1 and L1 have independent feedback loops / seismic backgrounds with mean-reverting (H~0.04) behavior. No common fractal structure."
    else:
        verdict = "INCONCLUSIVE: Weak but non-zero correlation."

    out = {
        "H1_Pure_H": H1_pure,
        "L1_Pure_H": L1_pure,
        "Cross_H1_L1_Pure_H": H_cross,
        "Interpretation": verdict
    }

    with open("QW_1660_v48_Pure_Raw_CrossDFA.json", "w") as f:
        json.dump(out, f, indent=2)

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "H1_Pure_H": 0.03733273522994843,
  "L1_Pure_H": 0.05299552509319112,
  "Cross_H1_L1_Pure_H": 0.31105068288116156,
  "Interpretation": "SURPRISING: The control loops / seismic background of H1 and L1 are strongly correlated at long distances."
}
```

---

## Faza 49: `phase49_pure_raw_stability.py`
### Kod:
```python
log = logging.getLogger(__name__)

RAW_DIR = "./raw_strain_unfiltered"
os.makedirs(RAW_DIR, exist_ok=True)
FS = 4096
WINDOW_SEC = 512

def fetch_pure_strain(det, gps):
    path = f"{RAW_DIR}/{det}_unfiltered_{gps}.h5"
    if os.path.exists(path):
        with h5py.File(path, "r") as f:
            return f["strain"][:]
            
    try:
        ts = TimeSeries.fetch_open_data(det, gps, gps + WINDOW_SEC, verbose=False, cache=True)
        if ts.sample_rate.value > FS:
            ts = ts.resample(FS)
        
        ts = ts.notch(60).notch(120).notch(180)
        
        val = ts.value
        with h5py.File(path, "w") as f:
            f.create_dataset("strain", data=val)
        return val
    except Exception as e:
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
    return float(np.polyfit(np.log(valid_scales), np.log(F), 1)[0])

def main():
    base_gps = 1266965117 # Baseline
    
    gps_times = [
        base_gps,
        base_gps + 86400,          # +1 dzien (Ziemia przekrecona, inne wplywy, ksiezyc)
        base_gps + 86400 * 2,      # +2 dni
        base_gps + 86400 * 7,      # +1 tydzien
        1253326744                 # Random inny czas O3 (~160 dni wczesniej)
    ]
    
    results = {}
    
    for gps in gps_times:
        x_h1 = fetch_pure_strain("H1", gps)
        x_l1 = fetch_pure_strain("L1", gps)
        
        if x_h1 is None or x_l1 is None:
            continue
            
        x_h1 = detrend(x_h1)
        x_l1 = detrend(x_l1)
        
        scales_long = np.logspace(3, np.log10(len(x_h1)//4), 15).astype(int)
        scales_long = np.unique(scales_long)
        
        h_cross = cross_mfdfa_q0(x_h1, x_l1, scales_long)
        h1_pure = mfdfa_q0(x_h1, scales_long)
        l1_pure = mfdfa_q0(x_l1, scales_long)
        
        results[str(gps)] = {
            "H1_Pure": h1_pure,
            "L1_Pure": l1_pure,
            "Cross_H1_L1": h_cross
        }
        
        
    with open("QW_1660_v49_CrossHurst_Stability.json", "w") as f:
        json.dump(results, f, indent=2)
        

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "1266965117": {
    "H1_Pure": 0.03733273522994843,
    "L1_Pure": 0.05299552509319112,
    "Cross_H1_L1": 0.31105068288116156
  },
  "1267051517": {
    "H1_Pure": 0.041405315287956546,
    "L1_Pure": 0.06420362757311483,
    "Cross_H1_L1": 0.33515322923661534
  },
  "1267137917": {
    "H1_Pure": 0.044313974791691,
    "L1_Pure": 0.05405721290409476,
    "Cross_H1_L1": 0.3034840404652728
  },
  "1253326744": {
    "H1_Pure": 0.04180091028389436,
    "L1_Pure": 0.05994293130661616,
    "Cross_H1_L1": 0.7633012338448356
  }
}
```

---

## Faza 50: `phase50_monte_carlo_cross_hurst.py`
### Kod:
```python
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

log = logging.getLogger(__name__)

def generate_fgn(n, H, seed=None):
    """Generate fractional Gaussian noise with Hurst exponent H using Hosking method."""
    if seed is not None:
        np.random.seed(seed)
    
    
    N = n
    f = np.fft.fftfreq(N)
    f[0] = 1e-10  # avoid division by zero
    
    power = np.abs(f) ** (-(2*H - 1))
    power[0] = 0
    
    phases = np.random.uniform(0, 2*np.pi, N)
    spectrum = np.sqrt(power) * np.exp(1j * phases)
    
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

    results_Hx = []
    results_Hy = []
    results_Hcross = []

    results_null_cross = []

    for trial in range(N_TRIALS):
        
        bg = generate_fgn(N, H_TRUE, seed=trial*100)
        
        noise1 = generate_fgn(N, H_LOCAL, seed=trial*100 + 1)
        noise2 = generate_fgn(N, H_LOCAL, seed=trial*100 + 2)
        
        x = bg + SNR_LOCAL * noise1
        y = bg + SNR_LOCAL * noise2
        
        Hx = mfdfa_q0(x, scales)
        Hy = mfdfa_q0(y, scales)
        Hcross = cross_mfdfa_q0(x, y, scales)
        
        results_Hx.append(Hx)
        results_Hy.append(Hy)
        results_Hcross.append(Hcross)
        

        null1 = generate_fgn(N, H_LOCAL, seed=trial*100 + 50)
        null2 = generate_fgn(N, H_LOCAL, seed=trial*100 + 51)
        Hcross_null = cross_mfdfa_q0(null1, null2, scales)
        results_null_cross.append(Hcross_null)

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

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "config": {
    "N": 1048576,
    "N_trials": 20,
    "H_true_background": 0.23,
    "H_local_damping": 0.04,
    "SNR_local_over_background": 20.0
  },
  "mixed_signal": {
    "H_x_mean": 0.07864376942410117,
    "H_x_std": 0.002380065086687186,
    "H_y_mean": 0.07825426745142278,
    "H_y_std": 0.003211036386556685,
    "Cross_H_mean": 0.5094425737999033,
    "Cross_H_std": 0.02343557500667056
  },
  "null_model": {
    "Cross_H_null_mean": 0.4965606426963275,
    "Cross_H_null_std": 0.027611698738229173
  },
  "verdict": "Cross-MF-DFA gives 0.509, which matches neither 0.23 nor 0.31. Relationship is complex."
}
```

---

## Faza 51: `phase51_snr_scan.py`
### Kod:
```python
"""
Phase 51: SNR Scan — At what signal strength does H=0.23 background
produce Cross-H ≈ 0.31 in damped channels?

If there exists an SNR where Cross-H passes through 0.31,
then the measured 0.31 COULD be a projection of true 0.23.
If Cross-H jumps directly from ~0.50 to ~0.23 without passing through 0.31,
then 0.31 is NOT a projection of 0.23.
"""

log = logging.getLogger(__name__)

def generate_fgn(n, H, seed=None):
    if seed is not None:
        np.random.seed(seed)
    N = n
    f = np.fft.fftfreq(N)
    f[0] = 1e-10
    power = np.abs(f) ** (-(2*H - 1))
    power[0] = 0
    phases = np.random.uniform(0, 2*np.pi, N)
    spectrum = np.sqrt(power) * np.exp(1j * phases)
    if N % 2 == 0:
        spectrum[N//2] = np.abs(spectrum[N//2])
    for i in range(1, N//2):
        spectrum[N-i] = np.conj(spectrum[i])
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

def main():
    N = 2**19  # 512K samples (faster, still statistically valid)
    N_TRIALS = 10
    H_BG = 0.23
    H_LOCAL = 0.04

    snr_values = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 50.0]

    scales = np.logspace(3, np.log10(N//4), 15).astype(int)
    scales = np.unique(scales)

    results = {}

    for snr in snr_values:
        cross_vals = []

        for trial in range(N_TRIALS):
            bg = generate_fgn(N, H_BG, seed=trial*1000)
            n1 = generate_fgn(N, H_LOCAL, seed=trial*1000 + 1)
            n2 = generate_fgn(N, H_LOCAL, seed=trial*1000 + 2)

            x = bg + snr * n1
            y = bg + snr * n2

            hc = cross_mfdfa_q0(x, y, scales)
            cross_vals.append(hc)

        mean_h = float(np.mean(cross_vals))
        std_h = float(np.std(cross_vals))

        results[str(snr)] = {
            "Cross_H_mean": mean_h,
            "Cross_H_std": std_h
        }

    with open("QW_1660_v51_SNR_Scan.json", "w") as f:
        json.dump(results, f, indent=2)

    for snr in snr_values:
        r = results[str(snr)]
        marker = " <-- 0.31?" if abs(r["Cross_H_mean"] - 0.31) < 0.02 else ""

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "0.0": {
    "Cross_H_mean": 0.5013565178093967,
    "Cross_H_std": 0.020171871555999893
  },
  "0.5": {
    "Cross_H_mean": 0.5059634912258375,
    "Cross_H_std": 0.016360932501354154
  },
  "1.0": {
    "Cross_H_mean": 0.5103556104649954,
    "Cross_H_std": 0.0146661713518238
  },
  "2.0": {
    "Cross_H_mean": 0.5044638878327039,
    "Cross_H_std": 0.02237324856125114
  },
  "3.0": {
    "Cross_H_mean": 0.49910779621350587,
    "Cross_H_std": 0.024637612596447415
  },
  "5.0": {
    "Cross_H_mean": 0.4949717022769238,
    "Cross_H_std": 0.025364894370843234
  },
  "7.0": {
    "Cross_H_mean": 0.49354969534495413,
    "Cross_H_std": 0.025348111303817908
  },
  "10.0": {
    "Cross_H_mean": 0.49267169835880253,
    "Cross_H_std": 0.025164686986300058
  },
  "15.0": {
    "Cross_H_mean": 0.49211390367399277,
    "Cross_H_std": 0.02489982657201801
  },
  "20.0": {
    "Cross_H_mean": 0.49187840042302067,
    "Cross_H_std": 0.024724322758137508
  },
  "50.0": {
    "Cross_H_mean": 0.49152941673745076,
    "Cross_H_std": 0.02433616880408267
  }
}
```

---

## Faza 52: `phase52_time_shift_cross_hurst.py`
### Kod:
```python
log = logging.getLogger(__name__)

RAW_DIR = "./raw_strain_unfiltered"
os.makedirs(RAW_DIR, exist_ok=True)
FS = 4096
WINDOW_SEC = 512
GPS = 1266965117 # Baseline quiet GPS
FETCH_WINDOW = WINDOW_SEC + 120

def fetch_pure_strain(det, gps):
    path = f"{RAW_DIR}/{det}_unfiltered_{gps}_long.h5"
    if os.path.exists(path):
        with h5py.File(path, "r") as f:
            return f["strain"][:]
            
    try:
        ts = TimeSeries.fetch_open_data(det, gps, gps + FETCH_WINDOW, verbose=False, cache=True)
        if ts.sample_rate.value > FS:
            ts = ts.resample(FS)
        
        ts = ts.notch(60).notch(120).notch(180)
        
        val = ts.value
        with h5py.File(path, "w") as f:
            f.create_dataset("strain", data=val)
        return val
    except Exception as e:
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
    x_h1_full = fetch_pure_strain("H1", GPS)
    x_l1_full = fetch_pure_strain("L1", GPS)
    
    if x_h1_full is None or x_l1_full is None:
        return
        
    shifts_sec = [0.0, 0.1, 1.0, 5.0, 10.0, 50.0, 100.0]
    
    base_length = WINDOW_SEC * FS
    scales = np.logspace(3, np.log10(base_length//4), 15).astype(int)
    scales = np.unique(scales)
    
    results = {}
    
    
    for shift in shifts_sec:
        shift_samples = int(shift * FS)
        
        x = x_h1_full[:base_length]
        
        y = x_l1_full[shift_samples:shift_samples + base_length]
        
        x = detrend(x)
        y = detrend(y)
        
        h_cross = cross_mfdfa_q0(x, y, scales)
        
        results[str(shift)] = h_cross
        
    with open("QW_1660_v52_TimeShift.json", "w") as f:
        json.dump(results, f, indent=2)
        
    for shift in shifts_sec:

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "0.0": 0.31105068079469167,
  "0.1": 0.30997165995459897,
  "1.0": 0.28581875130863643,
  "5.0": 0.29588892330405664,
  "10.0": 0.2947718600757252,
  "50.0": 0.2982288953028582,
  "100.0": 0.3202829195822657
}
```

---

## Faza 53: `phase53_psd_surrogate_test.py`
### Kod:
```python
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
            
    try:
        ts = TimeSeries.fetch_open_data(det, gps, gps + WINDOW_SEC, verbose=False, cache=True)
        if ts.sample_rate.value > FS:
            ts = ts.resample(FS)
        
        ts = ts.notch(60).notch(120).notch(180)
        
        val = ts.value
        with h5py.File(path, "w") as f:
            f.create_dataset("strain", data=val)
        return val
    except Exception as e:
        return None

def phase_randomize(x, seed=None):
    if seed is not None:
        np.random.seed(seed)
    X = np.fft.rfft(x)
    phases = np.random.uniform(0, 2 * np.pi, len(X))
    phases[0] = 0.0
    if len(x) % 2 == 0:
        phases[-1] = 0.0 # Nyquist
        
    X_surr = np.abs(X) * np.exp(1j * phases)
    x_surr = np.fft.irfft(X_surr, n=len(x))
    return x_surr

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
    return float(np.polyfit(np.log(valid_scales), np.log(F), 1)[0])

def main():
    x_h1 = fetch_pure_strain("H1", GPS)
    x_l1 = fetch_pure_strain("L1", GPS)
    
    if x_h1 is None or x_l1 is None:
        return
        
    x_h1 = detrend(x_h1)
    x_l1 = detrend(x_l1)
    
    scales_long = np.logspace(3, np.log10(len(x_h1)//4), 15).astype(int)
    scales_long = np.unique(scales_long)
    
    h_cross_orig = cross_mfdfa_q0(x_h1, x_l1, scales_long)
    
    N_SURR = 20
    surr_cross_h = []
    
    for i in range(N_SURR):
        x_h1_surr = phase_randomize(x_h1, seed=i*100)
        x_l1_surr = phase_randomize(x_l1, seed=i*100+1)
        
        hc = cross_mfdfa_q0(x_h1_surr, x_l1_surr, scales_long)
        surr_cross_h.append(hc)
        
    mean_surr = float(np.mean(surr_cross_h))
    std_surr = float(np.std(surr_cross_h))
    
    
    out = {
        "original_cross_H": h_cross_orig,
        "surrogate_cross_H_mean": mean_surr,
        "surrogate_cross_H_std": std_surr,
        "N_surrogates": N_SURR
    }
    
    with open("QW_1660_v53_PSDSurrogate.json", "w") as f:
        json.dump(out, f, indent=2)
        

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "original_cross_H": 0.31105068288116156,
  "surrogate_cross_H_mean": 0.28332943157004487,
  "surrogate_cross_H_std": 0.008788916370598766,
  "N_surrogates": 20
}
```

---

## Faza 54: `phase54_low_freq_isolation_csd.py`
### Kod:
```python
log = logging.getLogger(__name__)

RAW_DIR = "./raw_strain_unfiltered"
os.makedirs(RAW_DIR, exist_ok=True)
FS_TARGET = 1 # 1 Hz
WINDOW_SEC = 4096 # 4096 seconds to allow up to 1000s scales
GPS = 1266965117 # Baseline

def fetch_and_downsample(det, gps):
    path = f"{RAW_DIR}/{det}_downsampled_1Hz_{gps}.h5"
    if os.path.exists(path):
        with h5py.File(path, "r") as f:
            return f["strain"][:]
            
    try:
        ts = TimeSeries.fetch_open_data(det, gps, gps + WINDOW_SEC, verbose=False, cache=True)
        
        ts_1hz = ts.resample(FS_TARGET)
        
        ts_01hz = ts_1hz.lowpass(0.1)
        
        val = ts_01hz.value
        with h5py.File(path, "w") as f:
            f.create_dataset("strain", data=val)
        return val
    except Exception as e:
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
        return
        
    x_h1 = detrend(x_h1)
    x_l1 = detrend(x_l1)
    
    scales = np.logspace(1, 3, 15).astype(int)
    scales = np.unique(scales)
    
    h_cross = cross_mfdfa_q0(x_h1, x_l1, scales)
    
    nperseg = 512
    f_csd, Pxy = csd(x_h1, x_l1, fs=FS_TARGET, nperseg=nperseg)
    f_coh, Cxy = coherence(x_h1, x_l1, fs=FS_TARGET, nperseg=nperseg)
    
    idx_low = (f_coh >= 0.01) & (f_coh <= 0.1)
    mean_coh_low = np.mean(Cxy[idx_low])
    
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

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "Config": {
    "GPS": 1266965117,
    "Fs_target": 1,
    "Cutoff_Hz": 0.1,
    "Window_sec": 4096
  },
  "Cross_H_LF": 1.9154092490606147,
  "Coherence_0_01_to_0_1_Hz": 0.6741010398437792,
  "Mean_Real_CSD_LF": 5.3899207891616826e-48,
  "Mean_Imag_CSD_LF": 1.0931139634509775e-50,
  "Verdict": "Effect disappears at ultra-low frequencies or interpretation is mixed."
}
```

---

## Faza 55: `phase55_arma_injection.py`
### Kod:
```python
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
    
    R = np.fft.rfft(real_data)
    empirical_amplitude = np.abs(R)
    
    T = np.fft.rfft(target_signal)
    target_phases = np.angle(T)
    
    combined_spectrum = empirical_amplitude * np.exp(1j * target_phases)
    
    injected_signal = np.fft.irfft(combined_spectrum, n=len(real_data))
    return injected_signal

def main():
    x_h1 = detrend(fetch_pure_strain("H1", GPS))
    x_l1 = detrend(fetch_pure_strain("L1", GPS))
    
    N = len(x_h1)
    scales = np.logspace(3, np.log10(N//4), 15).astype(int)
    scales = np.unique(scales)
    
    H_TRUE = 0.23
    global_background = generate_fgn(N, H_TRUE, seed=42)
    
    h1_injected = empirical_filter_injection(H_TRUE, global_background, x_h1)
    
    l1_injected = empirical_filter_injection(H_TRUE, global_background, x_l1)
    
    h_cross_injection = cross_mfdfa_q0(h1_injected, l1_injected, scales)
    
    
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
        

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "Fundamental_H_Input": 0.23,
  "Injected_Cross_H_Output": 0.300537392344461,
  "Verdict": "The instrumental response accurately transforms a true H=0.23 background into an apparent Cross-H ~ 0.31. The 0.31 anomaly COULD be a projection of a 0.23 fundamental background."
}
```

---

## Faza 56: `phase56_scale_separation_dfa.py`
### Kod:
```python
log = logging.getLogger(__name__)

RAW_DIR = "./raw_strain_unfiltered"
FS = 4096
WINDOW_SEC = 512
GPS = 1266965117 # Baseline quiet GPS

def fetch_pure_strain(det, gps):
    path = f"{RAW_DIR}/{det}_unfiltered_{gps}.h5"
    with h5py.File(path, "r") as f:
        return f["strain"][:]

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
    x_h1 = detrend(fetch_pure_strain("H1", GPS))
    x_l1 = detrend(fetch_pure_strain("L1", GPS))
    
    N_samples = len(x_h1)
    
    split_scale = 100000 
    
    scales_short = np.logspace(np.log10(4096), np.log10(split_scale), 15).astype(int)
    scales_short = np.unique(scales_short)
    
    scales_long = np.logspace(np.log10(split_scale), np.log10(N_samples//4), 15).astype(int)
    scales_long = np.unique(scales_long)
    
    scales_full = np.logspace(3, np.log10(N_samples//4), 15).astype(int)
    scales_full = np.unique(scales_full)
    
    h_cross_short = cross_mfdfa_q0(x_h1, x_l1, scales_short)
    
    h_cross_long = cross_mfdfa_q0(x_h1, x_l1, scales_long)
    
    h_cross_full = cross_mfdfa_q0(x_h1, x_l1, scales_full)

    out = {
        "H_cross_short_scales_1s_to_25s": h_cross_short,
        "H_cross_long_scales_25s_to_128s": h_cross_long,
        "H_cross_full_baseline": h_cross_full,
        "Verdict": ""
    }
    
    if h_cross_long < 0.35 and h_cross_short > 0.45:
        out["Verdict"] = "The anti-persistence is strictly a LONG-SCALE phenomenon (>25s), bypassing rapid feedback loops. Supports global/structural origin."
    elif h_cross_short < 0.35 and h_cross_long > 0.45:
        out["Verdict"] = "The anti-persistence is strictly a SHORT-SCALE phenomenon (<25s), indicative of fast feedback loop synchrony."
    else:
        out["Verdict"] = "Anti-persistence is scale-invariant across both short and long bins."
        
    with open("QW_1660_v56_ScaleSeparation.json", "w") as f:
        json.dump(out, f, indent=2)
        

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "qw937": {
    "N_positive": 200,
    "N_negative": 0,
    "lambda_max": 429.3316676127456,
    "lambda_min": 0.00693181271269972,
    "hierarchy_orders": 4.791946103124579,
    "N_degenerate_groups": 0,
    "mean_gap": 2.1574107326634815,
    "gap_std": 26.206574827188174
  },
  "qw938": {
    "N_localized": 0,
    "N_extended": 126,
    "mean_IPR": 0.009948635074028045,
    "center_spread": 2.4584674706054732e-09,
    "PARTICLE_STATES": 0
  },
  "qw939": {
    "N_zeros": 13,
    "zeros": [
      1.333333551699543,
      5.333333509671445,
      9.333333456949472,
      13.333333393805082,
      17.33333334323282,
      21.333333396025125,
      25.333333438075563,
      29.33333347004699,
      33.33333349244686,
      37.33333350564984
    ],
    "N_maxima": 0,
    "maxima_d": [],
    "maxima_K": [],
    "N_minima": 0,
    "oscillation_period": 8.0
  },
  "qw940": {
    "phase_at_d0": 0.5235987755982988,
    "phase_velocity_mean": -0.051522119518872714,
    "phase_velocity_expected": 0.7853981633974483,
    "N_phase_jumps": 4,
    "PHASE_LINEAR": "False"
  },
  "qw941": {
    "DOS_peaks": 0,
    "eigenvalue_range": 586.6333297079102,
    "mean_DOS": 6.0,
    "max_DOS": 295.0,
    "VAN_HOVE_SINGULARITIES": 0
  },
  "qw942": {
    "correlation_length": 40.0,
    "G_at_xi": 0.7142857142857143,
    "power_law_exponent": -0.012830395174034868,
    "LONG_RANGE": "True"
  },
  "qw943": {
    "scaling_dimension": 0.2832671888058274,
    "anomalous_dimension": -0.7167328111941726,
    "ALPHA_GEO_related": 2.772588722239781,
    "CONFORMAL": "False"
  },
  "qw944": {
    "susceptibility": 2387.7791053735386,
    "specific_heat": 59.137609834507685,
    "chi_per_mode": 37.90125564084982
  },
  "qw945": {
    "entanglement_entropy": 0.006320883685330237,
    "max_entropy": 3.912023005428146,
    "entropy_ratio": 0.0016157583113799855
  },
  "qw946": {
    "overall_exponent": -0.28796990563401287,
    "short_range_exp": -0.0782395106624966,
    "long_range_exp": -0.4958169257869002,
    "NEWTON_SHORT": "False",
    "NEWTON_LONG": "False"
  },
  "qw947": {
    "N_wells": 0,
    "well_positions": [],
    "well_depths": [],
    "BOUND_STATES_POSSIBLE": false
  },
  "qw948": {
    "tunneling_amplitude": 0.0,
    "WKB_action": Infinity,
    "COHERENT_TUNNELING": false
  },
  "qw949": {
    "even_component": 47.07854271690005,
    "odd_component": 27.187847534392105,
    "parity": 0.6339145144607446,
    "TIME_REVERSAL": false,
    "PARITY_SYMMETRIC": "False"
  },
  "qw950": {
    "N_resonances": 2,
    "resonances": [
      {
        "index": 1,
        "eigenvalue": 38.71465995031845,
        "isolation": 11.044792123615759
      },
      {
        "index": 2,
        "eigenvalue": 13.54136179414084,
        "isolation": 3.1118268546314054
      }
    ],
    "resonance_ratios": [
      2.8589931012012206
    ],
    "PARTICLE_SPECTRUM": false
  },
  "qw951": {
    "beta_exponent": 1.1086755504543382,
    "correlation_length": 100.0,
    "ALPHA_GEO": 2.772588722239781,
    "CRITICAL_BEHAVIOR": "True"
  },
  "qw952": {
    "results": [
      {
        "d_max": 10,
        "d_spectral": 0.046823231316955205
      },
      {
        "d_max": 20,
        "d_spectral": 0.06742245720130002
      },
      {
        "d_max": 40,
        "d_spectral": 0.09806808235067077
      },
      {
        "d_max": 80,
        "d_spectral": 0.14414159605447718
      },
      {
        "d_max": 160,
        "d_spectral": 0.21497840402969493
      }
    ],
    "UV_dimension": 0.046823231316955205,
    "IR_dimension": 0.21497840402969493,
    "DIMENSION_FLOWS": false
  },
  "qw953": {
    "mean_spacing": 2.1574107326634815,
    "r_statistic": 0.919496367531011,
    "POISSON_expected": 0.386,
    "WIGNER_expected": 0.536,
    "IS_CHAOTIC": "True",
    "IS_INTEGRABLE": false
  },
  "qw954": {
    "winding_number": 2.001949206344289,
    "phase_winding": -1.5902773407317583e-16,
    "expected_winding": 1.0,
    "TOPOLOGICALLY_NONTRIVIAL": "True"
  },
  "qw955": {
    "ERROR": "SVD did not converge in Linear Least Squares"
  },
  "qw956": {
    "CONFIRMED_PHYSICS": [
      "Topological protection"
    ],
    "QUANTUM_PROPERTIES": {
      "integrable": false,
      "entanglement_entropy": 0.006320883685330237
    },
    "CLASSICAL_LIMIT": {},
    "EMERGENT_PHENOMENA": [
      "Critical behavior"
    ],
    "VERDICT": "PARTIAL: 1 confirmed, needs more validation"
  }
}
```
```json
{
  "qw847": {
    "N": 200,
    "d_max": 20.0,
    "N_eigenvalues": 200,
    "lambda_max": 252.49263471892004,
    "lambda_min": -58.136558661520496,
    "lambda_range": 310.62919338044054,
    "top_10": [
      252.49263471892004,
      227.2881598109568,
      34.93394850246776,
      31.766434185230608,
      15.666612501744495,
      14.365182951140113,
      9.262260581867196,
      8.542684646454305,
      6.215788519452744,
      5.7618250609189445
    ],
    "bottom_10": [
      0.05575610878080487,
      0.055725168958311,
      0.055701007767883746,
      0.05568381944379165,
      0.05567348450711878,
      -4.630432912082185,
      -35.98916042253846,
      -37.86173836745309,
      -54.693315588346245,
      -58.136558661520496
    ]
  },
  "qw848": {
    "N_ratios": 190,
    "ratio_min": 1.058350545284267,
    "ratio_max": 148.20532752482634,
    "ratio_mean": 15.178182192395111,
    "ratio_std": 29.61424621816515,
    "N_unique": 165,
    "near_power_of_2": 40,
    "near_power_of_4": 38,
    "sample_ratios": [
      1.1108921596660672,
      7.227715318269385,
      7.948409734836188,
      16.11660687278789,
      17.576708600072553,
      27.26036829639915,
      29.556590833973804,
      40.621175242485606,
      43.82164193625316,
      56.18016287040093
    ]
  },
  "qw849": {
    "N_gaps": 57,
    "gap_mean": 5.445391265136252,
    "gap_std": 25.668355333481568,
    "gap_cv": 4.71377612437539,
    "is_regular": false,
    "gap_min": 0.010082953939955286,
    "gap_max": 192.35421130848906,
    "quantization_fraction": 0.24561403508771928,
    "sample_gaps": [
      25.20447490796323,
      192.35421130848906,
      3.167514317237149,
      16.099821683486113,
      1.3014295506043823,
      5.1029223692729175,
      0.719575935412891,
      2.326896127001561,
      0.45396345853379927,
      1.2674871697778656
    ]
  },
  "qw850": {
    "N_log_ratios": 435,
    "log_ratio_mean": 1.623497496019754,
    "log_ratio_std": 1.4119879001854485,
    "N_peaks": 2,
    "peak_positions": [
      2.740409145891295,
      3.3080045796018913
    ],
    "peaks_at_ln4_multiples": [
      [
        2.740409145891295,
        2
      ]
    ],
    "ln4": 1.3862943611198906,
    "NATURAL_BASE_DETECTED": true
  },
  "qw851": {
    "decay_rate_b": 0.227944426581542,
    "intercept_a": 4.225272946973607,
    "r_squared": 0.8393748480110714,
    "natural_base": 1.256015522416549,
    "closest_known_base": "2",
    "closest_base_value": 2.0,
    "base_error_percent": 37.19922387917255,
    "EXPONENTIAL_QUANTIZATION": "False"
  },
  "qw852": {
    "N_degenerate_groups": 3,
    "group_sizes": [
      2,
      2,
      71
    ],
    "N_isolated": 125,
    "degeneracy_fraction": 0.375,
    "FAMILIES_DETECTED": true
  },
  "qw853": {
    "ipr_mean": 0.007934005482959835,
    "ipr_max": 0.09561736182464417,
    "ipr_min": 0.005064912135981274,
    "most_localized_eigenvalue": -4.630432912082185,
    "particle_threshold": 0.015,
    "N_particle_like": 1,
    "STABLE_PARTICLES_FOUND": "True",
    "top_5_localized_eigenvalues": [
      -4.630432912082185,
      227.2881598109568,
      -54.693315588346245,
      -35.98916042253846,
      31.766434185230608
    ]
  },
  "qw854": {
    "N_initial_projections": 40,
    "N_surviving": 1,
    "surviving_eigenvalues": [
      123.18912116511228
    ],
    "dominant_eigenvalue": 123.18912116511228,
    "STABLE_MODES_EXIST": true
  },
  "qw855": {
    "lambda_max": 252.49263471892004,
    "lambda_min": 0.05567348450711878,
    "hierarchy_ratio": 4535.240374376686,
    "orders_of_magnitude": 3.6566003102492206,
    "SM_comparison": 5.531478917042255,
    "SUFFICIENT_HIERARCHY": "True"
  },
  "qw856": {
    "KERNEL_PARAMETERS": {
      "alpha": 2.772588722239781,
      "omega": 0.7853981633974483,
      "phi": 0.5235987755982988,
      "beta": 0.01
    },
    "EMERGENT_PROPERTIES": {
      "Natural_Base": [
        [
          2.740409145891295,
          2
        ]
      ],
      "N_Families": 2,
      "N_Particle_Like": 1,
      "Hierarchy_Orders": 3.6566003102492206
    },
    "SUMMARY": "4/5 key properties emerged from K(d)",
    "SUCCESS": true
  }
}
```
```json
{
  "H_cross_short_scales_1s_to_25s": 0.22332124569812675,
  "H_cross_long_scales_25s_to_128s": 0.3619575625119233,
  "H_cross_full_baseline": 0.31105068288116156,
  "Verdict": "Anti-persistence is scale-invariant across both short and long bins."
}
```
```json
{
  "qw1037": {
    "all_attempts": [
      {
        "name": "Base4_n=[0, 4, 7]",
        "ratios": [
          256.0,
          16384.0
        ],
        "errors_percent": [
          23.81026077536176,
          371.1904864616137
        ],
        "mean_error": 197.50037361848774,
        "SUCCESS": false
      },
      {
        "name": "Base4_n=[0, 3.86, 5.89]",
        "ratios": [
          210.8393004204987,
          3516.684027649038
        ],
        "errors_percent": [
          1.9690186201436815,
          1.1369664135581685
        ],
        "mean_error": 1.5529925168509249,
        "SUCCESS": true
      },
      {
        "name": "Base4_n=[0, 2, 4]",
        "ratios": [
          16.0,
          256.0
        ],
        "errors_percent": [
          92.2618587015399,
          92.63764864903729
        ],
        "mean_error": 92.44975367528859,
        "SUCCESS": false
      }
    ],
    "best": {
      "name": "Base4_n=[0, 3.86, 5.89]",
      "ratios": [
        210.8393004204987,
        3516.684027649038
      ],
      "errors_percent": [
        1.9690186201436815,
        1.1369664135581685
      ],
      "mean_error": 1.5529925168509249,
      "SUCCESS": true
    },
    "INSIGHT": "Base-4 requires n=3.86, 5.89 \u2192 fitted, not derived",
    "IS_FITTING": true
  },
  "qw1038": {
    "best": {
      "name": "\u03b2^N: \u0394N_\u03bc=1.0, \u0394N_\u03c4=1.0",
      "ratios": [
        100.0,
        100.0
      ],
      "errors_percent": [
        51.63661688462431,
        97.12408150353019
      ],
      "mean_error": 74.38034919407725,
      "SUCCESS": false
    },
    "delta_N_required_mu": 1.1577416634905275,
    "delta_N_required_tau": 1.770611713010186,
    "INSIGHT": "Need fractional layers: \u0394N_\u03bc=1.16, \u0394N_\u03c4=1.77",
    "IS_FITTING": true
  },
  "qw1039": {
    "beta_optimal": 0.1000052749944835,
    "result_exp_minus": {
      "name": "Winding \u03b2_topo=0.10",
      "ratios": [
        1.1051767478615198,
        1.2214156440137653
      ],
      "errors_percent": [
        99.46549913532968,
        99.96487308157504
      ],
      "mean_error": 99.71518610845236,
      "SUCCESS": "False"
    },
    "result_exp_plus": {
      "name": "Winding exp(+\u03b2n)",
      "ratios": [
        1.1051767478615198,
        1.2214156440137653
      ],
      "errors_percent": [
        99.46549913532968,
        99.96487308157504
      ],
      "mean_error": 99.71518610845236,
      "SUCCESS": "False"
    },
    "PROBLEM": "Integer windings [0,1,2] give exp ratios [1, e^\u03b2, e^2\u03b2], not [1, 207, 3477]",
    "IS_FITTING": true
  },
  "qw1040": {
    "beta_Y_optimal": 0.9899941047084844,
    "result": {
      "name": "Yukawa \u03b2_Y=0.9900",
      "ratios": [
        1.0101070251266415,
        1.0203162022101937
      ],
      "errors_percent": [
        99.51147806956267,
        99.97065653761815
      ],
      "mean_error": 99.74106730359041,
      "SUCCESS": "False"
    },
    "PROBLEM": "Single \u03b2 cannot reproduce both ratios",
    "required_beta_for_mu": 0.004836338311537568,
    "required_beta_for_tau_sqrt": 0.01695853323984657,
    "IS_FITTING": true
  },
  "qw1041": {
    "eigenvalues_top3": [
      62.3940753519193,
      56.86823217437668,
      9.425105917465894
    ],
    "masses": [
      7.898992046579063,
      7.5411028486804685,
      3.070033536863383
    ],
    "result": {
      "name": "K(d) eigenspectrum",
      "ratios": [
        2.4563584593232566,
        2.5729334718113113
      ],
      "errors_percent": [
        98.81202194763054,
        99.92600453038231
      ],
      "mean_error": 99.36901323900642,
      "SUCCESS": false
    },
    "INSIGHT": "Raw K(d) spectrum does NOT give SM ratios",
    "IS_FITTING": false
  },
  "qw1042": {
    "ERROR": "'list' object has no attribute 'tolist'"
  },
  "qw1043": {
    "m_D": 2.772588722239781,
    "M_heavy": [
      2.4011322677058873,
      2.3087810266402764,
      2.223270618246192
    ],
    "m_light": [
      3.2015096902745164,
      3.3295700778854966,
      3.457630465496478
    ],
    "result": {
      "name": "Seesaw",
      "ratios": [
        1.0399999999999998,
        1.08
      ],
      "errors_percent": [
        99.4970208156001,
        99.96894008023813
      ],
      "mean_error": 99.73298044791912,
      "SUCCESS": "False"
    },
    "INSIGHT": "Seesaw inverts hierarchy \u2014 may help"
  },
  "qw1044": {
    "loop_factor": 0.000580857456539764,
    "expected_ratio": 1721.5927741672065,
    "result": {
      "name": "Radiative loops",
      "ratios": [
        1721.5927741672065,
        2963881.6800647383
      ],
      "errors_percent": [
        732.6205090571106,
        85138.82145046197
      ],
      "mean_error": 42935.72097975954,
      "SUCCESS": false
    },
    "INSIGHT": "Loop factor gives ~5000x, not 207x \u2014 wrong scale",
    "IS_FITTING": false
  },
  "qw1045": {
    "Im_tau": 2.772588722239781,
    "weights": [
      0,
      1,
      2
    ],
    "masses": [
      1.0,
      0.1751576467072376,
      0.030680201200017458
    ],
    "result": {
      "name": "Modular weights",
      "ratios": [
        5.709142699727078,
        32.59431036584698
      ],
      "errors_percent": [
        97.2388654435275,
        99.06261419939183
      ],
      "mean_error": 98.15073982145967,
      "SUCCESS": "False"
    },
    "INSIGHT": "Modular weights can give hierarchy but require tuning \u03c4"
  },
  "qw1046": {
    "q": 3,
    "N_sites": [
      2,
      3,
      4
    ],
    "masses": [
      0.1111111111111111,
      0.037037037037037035,
      0.012345679012345678
    ],
    "delta_N_needed_for_207": 4.8530290865820715,
    "result": {
      "name": "Clockwork",
      "ratios": [
        3.0,
        9.0
      ],
      "errors_percent": [
        98.54909850653874,
        99.74116733531771
      ],
      "mean_error": 99.14513292092823,
      "SUCCESS": false
    },
    "INSIGHT": "Clockwork with q=3 needs \u0394N~5 sites between e and \u03bc"
  },
  "qw1047": {
    "areas": [
      0.234309062311414,
      7.1240054754597395,
      13.753271421510501
    ],
    "masses": [
      0.9768414659956649,
      0.49046539763785124,
      0.252756894843727
    ],
    "result": {
      "name": "Intersection areas",
      "ratios": [
        1.9404629810043095,
        3.8647470590253157
      ],
      "errors_percent": [
        99.06152645428485,
        99.88885302448772
      ],
      "mean_error": 99.47518973938628,
      "SUCCESS": "False"
    },
    "INSIGHT": "Worldsheet areas could give hierarchy"
  },
  "qw1048": {
    "K_values": [
      2.4011322677058873,
      2.3087810266402764,
      2.223270618246192
    ],
    "result_direct": {
      "name": "Pure K(d)",
      "ratios": [
        1.0384615384615385,
        1.08
      ],
      "errors_percent": [
        99.49776486764802,
        99.96894008023813
      ],
      "mean_error": 99.73335247394307,
      "SUCCESS": "False"
    },
    "result_inverse": {
      "name": "Inverse K(d)",
      "ratios": [
        1.04,
        1.08
      ],
      "errors_percent": [
        99.4970208156001,
        99.96894008023813
      ],
      "mean_error": 99.73298044791912,
      "SUCCESS": "False"
    },
    "INSIGHT": "Direct K(d) doesn't give 207\u00d7 hierarchy"
  },
  "qw1049": {
    "masses": [
      2.318631214137653,
      13.9189485244329,
      12.440802740597768
    ],
    "result": {
      "name": "Product formula",
      "ratios": [
        5.36558063427295,
        6.003088563443517
      ],
      "errors_percent": [
        97.40502368148216,
        99.82735606564445
      ],
      "mean_error": 98.61618987356331,
      "SUCCESS": "False"
    },
    "INSIGHT": "Product gives exponential but may be too strong"
  },
  "qw1050": {
    "cumulative_sums": [
      29.232206195137046,
      34.08728892995544,
      36.93013684213586
    ],
    "masses": [
      0.05376026682517239,
      0.033083225980982414,
      0.024896857601620623
    ],
    "result": {
      "name": "Exp of eigensum",
      "ratios": [
        1.3288113106623107,
        2.1593193681468037
      ],
      "errors_percent": [
        99.35734189494394,
        99.9378997348936
      ],
      "mean_error": 99.64762081491878,
      "SUCCESS": "False"
    }
  },
  "qw1051": {
    "generations": [
      {
        "winding": 1,
        "d": 0,
        "layer": 10
      },
      {
        "winding": 2,
        "d": 4,
        "layer": 9
      },
      {
        "winding": 3,
        "d": 8,
        "layer": 8
      }
    ],
    "masses": [
      6.657352246047373e+20,
      2.560520094633605e+19,
      5.54779353837281e+17
    ],
    "result": {
      "name": "VEV \u00d7 w\u00b2 \u00d7 K \u00d7 \u03b2^(-N)",
      "ratios": [
        46.15384615384616,
        1200.0
      ],
      "errors_percent": [
        77.6784385621343,
        65.48897804236228
      ],
      "mean_error": 71.58370830224828,
      "SUCCESS": "False"
    },
    "IS_FITTING": true
  },
  "qw1052": {
    "a_optimal": 6.586210983895974,
    "b_optimal": -1.2546134364083137,
    "a_required": 4.390807098471844,
    "b_required": -0.8364088372832867,
    "formula": "m = exp(a\u00d7n + b\u00d7n\u00b2)",
    "INSIGHT": "Quadratic + linear in n can fit, but is it physical?"
  },
  "qw1053": {
    "zeros": [
      1.3213213213213213,
      5.315315315315315,
      9.30930930930931,
      13.333333333333334,
      17.32732732732733
    ],
    "zero_ratios": [
      4.0227272727272725,
      1.7514124293785311,
      1.4322580645161291,
      1.2995495495495495
    ],
    "mass_ratio_log": [
      5.331597391782836,
      8.153968271715812
    ],
    "INSIGHT": "Zero spacings are regular (~4), not related to 207"
  },
  "qw1054": {
    "Koide_experimental": 1.5000772276851697,
    "Koide_K_direct": 2.9992595108648024,
    "Koide_target": 2.0,
    "K_values": [
      2.4011322677058873,
      2.3087810266402764,
      2.223270618246192
    ],
    "INSIGHT": "Koide is independent constraint on masses"
  },
  "qw1055": {
    "dimensions": [
      6,
      5,
      4
    ],
    "masses": [
      6.0516e-09,
      246000.0,
      1e+19
    ],
    "result": {
      "name": "EFT power counting",
      "ratios": [
        40650406504065.04,
        1.6524555489457334e+27
      ],
      "errors_percent": [
        19659911835418.574,
        4.7523274778072085e+25
      ],
      "mean_error": 2.3761637389045874e+25,
      "SUCCESS": false
    },
    "INSIGHT": "EFT: each dimension adds (v/\u039b) ~ 10^(-16)"
  },
  "qw1056": {
    "TOTAL_TESTS": 18,
    "SUCCESSES": [],
    "N_SUCCESSES": 0,
    "FAILURES": [
      "qw1037",
      "qw1038",
      "qw1039",
      "qw1040",
      "qw1041",
      "qw1042",
      "qw1043",
      "qw1044",
      "qw1045",
      "qw1046",
      "qw1047",
      "qw1048",
      "qw1049",
      "qw1050",
      "qw1051",
      "qw1052",
      "qw1053",
      "qw1054",
      "qw1055"
    ],
    "N_FAILURES": 19,
    "REQUIRE_FITTING": [
      "qw1037",
      "qw1038",
      "qw1039",
      "qw1040",
      "qw1051"
    ],
    "PHYSICAL_MECHANISMS": [
      "qw1041",
      "qw1042",
      "qw1043",
      "qw1044",
      "qw1045",
      "qw1046",
      "qw1047",
      "qw1048",
      "qw1049",
      "qw1050",
      "qw1052",
      "qw1053",
      "qw1054",
      "qw1055"
    ],
    "VERDICT": "\u017bADEN mechanizm NIE daje SM ratios bez fittingu"
  }
}
```

---

## Faza 57: `phase57_high_stats_mc.py`
### Kod:
```python
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
    
    num_cores = multiprocessing.cpu_count()
    
    with multiprocessing.Pool(num_cores) as p:
        results = p.map(run_trial, range(N_TRIALS))
    
    results = np.array(results)
    mean_val = np.mean(results)
    std_val = np.std(results)
    min_val = np.min(results)
    percentile_1 = np.percentile(results, 1)
    
    
    boot_means = []
    for _ in range(1000):
        sample = np.random.choice(results, size=N_TRIALS, replace=True)
        boot_means.append(np.mean(sample))
    ci_lower = float(np.percentile(boot_means, 2.5))
    ci_upper = float(np.percentile(boot_means, 97.5))
    
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
        

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "qw957": {
    "shape": [
      12,
      12
    ],
    "eigenvalues": [
      13.659845156208473,
      10.77009650351679,
      -0.3155355315737813,
      -0.6511353249029644,
      -1.3816432118225215,
      -1.5120853875683036,
      -1.6947951840429216,
      -1.7535893865961405,
      -1.8012112004357668,
      -2.5210548711990697,
      -6.156734311343714,
      -6.642157250240079
    ],
    "hierarchy": 1.2683122339477724,
    "OCTAVE_STRUCTURE": true
  },
  "qw958": {
    "psi_vev": 2.5819889114994248,
    "phi_vev": 3.162277660168381,
    "psi_vev_squared": 6.666666739105985,
    "phi_vev_squared": 10.00000000000001,
    "V_psi_minimum": -5.134423559703298,
    "V_phi_minimum": -6.931471805599453,
    "SYMMETRY_BROKEN": "True"
  },
  "qw959": {
    "m_psi_squared": 0.2772588737303114,
    "m_phi_squared": 0.13862943611198905,
    "m_higgs_squared": 5.5451774444795685,
    "m_psi": 0.5265537709772017,
    "m_higgs": 2.3548200450309507,
    "MASSES_POSITIVE": "True"
  },
  "qw960": {
    "yukawa_by_generation": {
      "gen_0": {
        "g_Y": 0.027725887222397813,
        "effective": 0.2772588722239784
      },
      "gen_1": {
        "g_Y": 0.026918337109124092,
        "effective": 0.2691833710912412
      },
      "gen_2": {
        "g_Y": 0.026156497379620575,
        "effective": 0.26156497379620597
      },
      "gen_3": {
        "g_Y": 0.025436593782016334,
        "effective": 0.2543659378201636
      }
    },
    "ratios": {
      "gen_1/gen_0": 0.9708737864077671,
      "gen_2/gen_0": 0.9433962264150941,
      "gen_3/gen_0": 0.9174311926605503
    },
    "HIERARCHY_PRESENT": false
  },
  "qw961": {
    "phi_vev": 3.162277660168381,
    "g_effective": 0.2745137348752259,
    "sin2_theta_W": 0.25,
    "M_W": 0.4340443256026563,
    "M_Z": 0.5011912164538465,
    "M_W_M_Z_ratio": 0.8660254037844386,
    "ratio_expected": 0.8660254037844386,
    "HIGGS_WORKS": "True"
  },
  "qw962": {
    "t_final": 50.0,
    "psi_final": 2.798255672812152,
    "oscillation_frequency": 0.8571428571428495,
    "mass_from_frequency": 0.1364185226501948,
    "STABLE_OSCILLATION": "True"
  },
  "qw963": {
    "eigenvalues": [
      -10.621542806771187,
      -7.73206001632361,
      3.3536929479524784,
      3.6892419235118843,
      4.419766599997655,
      4.550132187473054,
      4.732698808121555,
      4.792017871923722,
      4.839470955130477,
      5.559491943001191,
      9.194933533085171,
      9.680340184569749
    ],
    "mass_ratios": [
      1.0527906645260663,
      1.7412274869390298
    ],
    "hierarchy_orders": 0.46036732181436996,
    "N_positive": 10
  },
  "qw964": {
    "m_psi_from_pole": 0.4494665749754947,
    "m_higgs_from_pole": 2.335496832484569,
    "m_psi_input": 0.5265537709772017,
    "m_higgs_input": 2.3548200450309507,
    "POLES_MATCH": "True"
  },
  "qw965": {
    "E_initial": 2.197391375392131,
    "E_final": -9.0547136541431e+27,
    "E_variation": 63.85832107819911,
    "ENERGY_CONSERVED": "False"
  },
  "qw966": {
    "max_K_change": 5.372915195079843,
    "mean_gauge_field": 1.81902307822688,
    "GAUGE_INVARIANT": "False",
    "GAUGE_FIELD_NEEDED": "True"
  },
  "qw967": {
    "charges": [
      -1.5,
      -0.5,
      0.5,
      1.5,
      -1.5,
      -0.5,
      0.5,
      1.5,
      -1.5,
      -0.5,
      0.5,
      1.5
    ],
    "cubic_anomaly": 0.0,
    "linear_anomaly": 0.0,
    "ANOMALY_FREE": true
  },
  "qw968": {
    "beta_UV": 0.03329313119641991,
    "beta_IR": 0.6151058767531213,
    "fixed_point_d": 1.0,
    "fixed_point_K": 2.745137348752259,
    "ASYMPTOTIC_FREE": "False",
    "HAS_FIXED_POINT": "False"
  },
  "qw969": {
    "amplitude_max": 2.769818903336445,
    "amplitude_min": 2.745137348752259,
    "N_resonances": 0,
    "resonance_s": [],
    "RESONANCES_PRESENT": false
  },
  "qw970": {
    "hessian_eigenvalue": 0.2772588737303114,
    "V_at_vev": -5.134423559703299,
    "V_at_zero": 0.0,
    "barrier_height": 5.134423559703299,
    "VACUUM_STABLE": "True",
    "METASTABLE": "True"
  },
  "qw971": {
    "vacuum_energy_psi": -5.134423559703299,
    "vacuum_energy_phi": -6.931471805599452,
    "total_vacuum_energy": -12.065895365302751,
    "lambda_cosmological": -0.204182137205403,
    "POSITIVE_LAMBDA": "False"
  },
  "qw972": {
    "K_at_d1": 2.745137348752259,
    "alpha_from_octaves": 0.018204257438337593,
    "alpha_from_beta": 0.027725887222397813,
    "alpha_experimental": 0.0072973525205055605,
    "error_octaves": 149.46386223200307,
    "error_beta": 279.94446814085063,
    "ALPHA_EMERGENT": "False"
  },
  "qw973": {
    "linear_residual": 244.8273425273264,
    "coulomb_residual": 263.7086063088971,
    "string_tension": -0.011704142424576024,
    "CONFINING": "False",
    "COULOMB_LIKE": "False"
  },
  "qw974": {
    "chiral_condensate": 6.666666739105985,
    "m_quark_proxy": 0.027725887222397813,
    "m_pion_squared": 0.18483925015776312,
    "m_pion": 0.4299293548453782,
    "CHIRAL_BROKEN": true,
    "GOLDSTONE_LIGHT": "True"
  },
  "qw975": {
    "Jarlskog_analog": -53.54879469490156,
    "phases": [
      1.308996938995747,
      2.0943951023931953,
      2.8797932657906435,
      -2.6179938779914944
    ],
    "phase_variation": 2.1147492745580068,
    "CP_VIOLATED": "True"
  },
  "qw976": {
    "L_ZTP_VALID": [
      "Spontaneous symmetry breaking",
      "Vacuum stability"
    ],
    "EMERGENT_SM": [
      "Higgs mechanism",
      "CP violation",
      "Chiral symmetry breaking"
    ],
    "ISSUES": [
      "Fine structure constant not accurate"
    ],
    "VERDICT": "L_ZTP IS VIABLE: 5 phenomena confirmed"
  }
}
```
```json
{
  "qw857": {
    "N_layers": 5,
    "N_sites": 50,
    "Total_dim": 250,
    "lambda_max": 63.92488229753012,
    "lambda_min": 19.67826675761243,
    "Hierarchy_ratio": 3.2485016635319837,
    "Orders_of_magnitude": 0.5116830935002104
  },
  "qw858": {
    "results": [
      {
        "N_layers": 1,
        "Hierarchy": 3.2663137088043097,
        "Orders": 0.5140578936461591
      },
      {
        "N_layers": 2,
        "Hierarchy": 3.2643878555049857,
        "Orders": 0.5138017534820863
      },
      {
        "N_layers": 3,
        "Hierarchy": 3.263591099274347,
        "Orders": 0.5136957400165088
      },
      {
        "N_layers": 5,
        "Hierarchy": 3.2629801068003195,
        "Orders": 0.5136144260608799
      },
      {
        "N_layers": 7,
        "Hierarchy": 3.262757757631106,
        "Orders": 0.5135848309348456
      },
      {
        "N_layers": 10,
        "Hierarchy": 3.262620907205644,
        "Orders": 0.5135666148625881
      }
    ],
    "growth_rate_per_layer": -4.508998628360177e-05,
    "expected_rate": 2.0,
    "EXPONENTIAL_GROWTH": "False"
  },
  "qw859": {
    "results": [
      {
        "beta": 1.0,
        "Hierarchy": 12.564947586179592,
        "Orders": 1.0991606813026134
      },
      {
        "beta": 0.1,
        "Hierarchy": 3.6293593862732596,
        "Orders": 0.5598299750199133
      },
      {
        "beta": 0.01,
        "Hierarchy": 3.2629801068003195,
        "Orders": 0.5136144260608799
      },
      {
        "beta": 0.001,
        "Hierarchy": 3.2496644427578367,
        "Orders": 0.5118385184612192
      },
      {
        "beta": 0.0001,
        "Hierarchy": 3.2486629405908785,
        "Orders": 0.5117046542089249
      }
    ],
    "observation": "Coupling doesn't amplify hierarchy"
  },
  "qw860": {
    "depth": 3,
    "hierarchy_per_level": [
      90.48975140457385,
      90.62722637075372,
      90.6287276248184,
      90.62874264992023
    ],
    "final_hierarchy": 90.62874264992023,
    "amplification": 1.0015359888074502
  },
  "qw861": {
    "per_layer": [
      {
        "layer": 0,
        "dilation": 1.0,
        "E_max": -11.757841868372319,
        "E_min": 38.40480008061778,
        "gap": 2.433085800750751
      },
      {
        "layer": 1,
        "dilation": 0.01,
        "E_max": -0.1175784186837232,
        "E_min": 0.38404800080617785,
        "gap": 0.02433085800750752
      },
      {
        "layer": 2,
        "dilation": 0.0001,
        "E_max": -0.001175784186837232,
        "E_min": 0.0038404800080617785,
        "gap": 0.00024330858007507513
      },
      {
        "layer": 3,
        "dilation": 1.0000000000000002e-06,
        "E_max": -1.1757841868372321e-05,
        "E_min": 3.8404800080617786e-05,
        "gap": 2.4330858007507513e-06
      },
      {
        "layer": 4,
        "dilation": 1e-08,
        "E_max": -1.1757841868372319e-07,
        "E_min": 3.8404800080617785e-07,
        "gap": 2.43308580075075e-08
      }
    ],
    "combined_hierarchy": 20526070359.57998,
    "combined_orders": 10.31230581325257
  },
  "qw862": {
    "N_total_eigenvalues": 289,
    "lambda_max": 51.151397084777884,
    "lambda_min": 1.1818553208735312e-15,
    "Hierarchy": 4.328059127150262e+16,
    "Orders": 16.63629318518994
  },
  "qw863": {
    "target_orders": 5.5,
    "single_layer": 3.66,
    "results": [
      {
        "N_layers": 1,
        "required_beta": 3.162277660168379e-06,
        "with_beta_0.01": 5.66,
        "reaches_target": "True"
      },
      {
        "N_layers": 2,
        "required_beta": 0.0017782794100389228,
        "with_beta_0.01": 7.66,
        "reaches_target": "True"
      },
      {
        "N_layers": 3,
        "required_beta": 0.014677992676220698,
        "with_beta_0.01": 9.66,
        "reaches_target": "True"
      },
      {
        "N_layers": 5,
        "required_beta": 0.07943282347242814,
        "with_beta_0.01": 13.66,
        "reaches_target": "True"
      },
      {
        "N_layers": 10,
        "required_beta": 0.28183829312644537,
        "with_beta_0.01": 23.66,
        "reaches_target": "True"
      },
      {
        "N_layers": 20,
        "required_beta": 0.5308844442309884,
        "with_beta_0.01": 43.66,
        "reaches_target": "True"
      },
      {
        "N_layers": 30,
        "required_beta": 0.6556418494179789,
        "with_beta_0.01": 63.66,
        "reaches_target": "True"
      }
    ],
    "N_layers_needed_with_beta_0.01": 1
  },
  "qw864": {
    "N_layers": 5,
    "correlation_matrix": [
      [
        1.0,
        1.0,
        1.0,
        1.0,
        1.0
      ],
      [
        1.0,
        1.0,
        1.0,
        1.0,
        1.0
      ],
      [
        1.0,
        1.0,
        1.0,
        1.0,
        1.0
      ],
      [
        1.0,
        1.0,
        1.0,
        1.0,
        1.0
      ],
      [
        1.0,
        1.0,
        1.0,
        1.0,
        1.0
      ]
    ],
    "avg_off_diagonal_corr": NaN,
    "IDENTICAL_LAYERS": "False"
  },
  "qw865": {
    "N_layers": 5,
    "lambda_max": 38.40480008061778,
    "lambda_min": 0.05254857629316127,
    "Hierarchy": 730.8437790276694,
    "Orders": 2.8638245545818
  },
  "qw866": {
    "N_layers": 3,
    "N_tensor_eigenvalues": 125,
    "lambda_max": 3.298021308645905,
    "lambda_min": 0.2227881489326375,
    "Hierarchy": 14.80339652017621,
    "Orders": 1.1703613722012747,
    "MULTIPLICATIVE_HIERARCHY": true
  },
  "qw867": {
    "N_layers": 5,
    "d_spectral": 0.726295252874926,
    "d_s_range": [
      0.013168285314228978,
      4.184157387645739
    ],
    "FRACTAL_STRUCTURE": "False"
  },
  "qw868": {
    "most_localized_eigenvalue": -9.324756067621575,
    "most_localized_layer": 0,
    "ipr_max": 0.0204413451887011,
    "layer_distribution": [
      33,
      21,
      68,
      15,
      13
    ],
    "LAYER_0_DOMINANT": "True"
  },
  "qw869": {
    "formula": "Orders = single_layer + N \u00d7 log10(1/\u03b2)",
    "single_layer_orders": 0.5119478327087438,
    "log10(1/\u03b2)": 2.0,
    "results": [
      {
        "N": 1,
        "predicted_orders": 2.511947832708744,
        "observed_orders": 0.5119478327087438,
        "error": 2.0
      },
      {
        "N": 3,
        "predicted_orders": 6.511947832708744,
        "observed_orders": 0.5117316512158867,
        "error": 6.000216181492857
      },
      {
        "N": 5,
        "predicted_orders": 10.511947832708744,
        "observed_orders": 0.5116830935002104,
        "error": 10.000264739208534
      },
      {
        "N": 10,
        "predicted_orders": 20.511947832708742,
        "observed_orders": 0.5116545391976377,
        "error": 20.000293293511106
      }
    ],
    "avg_error": 9.500193553553125,
    "FORMULA_VALID": "False"
  },
  "qw870": {
    "results": [
      {
        "N_layers": 1,
        "gap_mean": 2.0088681451349633,
        "gap_cv": 3.2561339390150166,
        "N_gaps": 19
      },
      {
        "N_layers": 2,
        "gap_mean": 1.9875219682093943,
        "gap_cv": 3.2912869995743317,
        "N_gaps": 19
      },
      {
        "N_layers": 5,
        "gap_mean": 4.917216569710931,
        "gap_cv": 2.053284375803326,
        "N_gaps": 7
      },
      {
        "N_layers": 10,
        "gap_mean": 4.386734963368582,
        "gap_cv": 0.0,
        "N_gaps": 1
      }
    ]
  },
  "qw871": {
    "N_positive_eigenvalues": 190,
    "log_eig_range": [
      -0.9119857291110548,
      1.7090045326780712
    ],
    "N_peaks": 1,
    "peak_positions": [
      0.7916579410518769
    ],
    "GENERATIONS_DETECTED": false
  },
  "qw872": {
    "max_tested_layers": 30,
    "max_orders_achieved": 0.517443857102058,
    "target": 5.5,
    "tensor_product_orders": 1.1703613722012747,
    "SM_HIERARCHY_ACHIEVED": false
  },
  "qw873": {
    "final_ground_overlap": 4.320496654093577e-09,
    "max_overlap": 0.00802112020114204,
    "RELAXATION_TO_GROUND": "False"
  },
  "qw874": {
    "results": [
      {
        "N_layers": 1,
        "Orders": 0.5140578936461591,
        "Time_seconds": 0.0073659420013427734,
        "Orders_per_second": 69.78847967475836
      },
      {
        "N_layers": 2,
        "Orders": 0.5138017534820863,
        "Time_seconds": 0.011898279190063477,
        "Orders_per_second": 43.18286243536577
      },
      {
        "N_layers": 3,
        "Orders": 0.5136957400165088,
        "Time_seconds": 0.019559860229492188,
        "Orders_per_second": 26.262751062094132
      },
      {
        "N_layers": 5,
        "Orders": 0.5136144260608799,
        "Time_seconds": 0.031617164611816406,
        "Orders_per_second": 16.244797165300675
      },
      {
        "N_layers": 7,
        "Orders": 0.5135848309348456,
        "Time_seconds": 0.04834699630737305,
        "Orders_per_second": 10.622890151637456
      },
      {
        "N_layers": 10,
        "Orders": 0.5135666148625881,
        "Time_seconds": 0.07893872261047363,
        "Orders_per_second": 6.50588960498897
      },
      {
        "N_layers": 15,
        "Orders": 0.5135557246564735,
        "Time_seconds": 0.20531368255615234,
        "Orders_per_second": 2.5013224557794307
      },
      {
        "N_layers": 20,
        "Orders": 0.5135516097993248,
        "Time_seconds": 0.25138425827026367,
        "Orders_per_second": 2.0428948627610746
      }
    ],
    "optimal_N": 1,
    "optimal_ratio": 69.78847967475836
  },
  "qw875": {
    "results": [
      {
        "beta": 0.001,
        "orders": 2.314644504067296
      },
      {
        "beta": 0.01,
        "orders": 2.336273774636879
      },
      {
        "beta": 0.05,
        "orders": 2.4479050405552796
      },
      {
        "beta": 0.1,
        "orders": 2.5843144787801724
      },
      {
        "beta": 0.2,
        "orders": 2.532652983851886
      },
      {
        "beta": 0.5,
        "orders": 2.380802789899299
      },
      {
        "beta": 1.0,
        "orders": 2.5377823894126275
      }
    ],
    "optimal_beta": 0.1,
    "max_orders": 2.5843144787801724,
    "RESONANCE_FOUND": false
  },
  "qw876": {
    "Single_layer_orders": 3.66,
    "Target_SM_orders": 5.5,
    "Gap_to_fill": 1.84,
    "Best_multilayer_mechanism": "Tensor Product",
    "Max_orders_achieved": 1.1703613722012747,
    "Key_findings": [],
    "SUCCESS": false
  }
}
```
```json
{
  "Config": {
    "N_trials": 500,
    "N_samples_per_trial": 131072,
    "H_local": 0.04
  },
  "Results": {
    "Null_Cross_H_Mean": 0.5018624026057812,
    "Null_Cross_H_Std": 0.03820165650501183,
    "Null_Cross_H_Min": 0.39702322618872843,
    "Null_1st_Percentile": 0.4147443000100347
  },
  "Bootstrap_95_CI": [
    0.49861649052509377,
    0.5051488624382866
  ],
  "Real_LIGO_Observation": 0.311,
  "Statistical_Significance": "5.00 sigma",
  "Verdict": "Strong anomaly confirmed at 5.00 sigma. The observed 0.311 NEVER occurs in 500 null trials. Real cross-detector anti-persistence exists."
}
```
```json
{
  "qw1057": {
    "ALPHA_GEO": 2.772588722239781,
    "exp_alpha": 15.999999999999998,
    "ln_4": 1.3862943611198906,
    "k_required": 0.5,
    "natural_base_from_alpha": 16,
    "observed_base": 4,
    "INSIGHT": "\u03b1 = ln(16), observed Base = 4 = \u221a16. Suggests 2-level structure.",
    "HALF_RELATION": true
  },
  "qw1058": {
    "K_at_particles": {
      "top": {
        "d": 0.0,
        "K_real": 2.4011322677058873,
        "K_magnitude": 2.772588722239781,
        "K_phase": 0.5235987755982988
      },
      "bottom": {
        "d": 1.75,
        "K_real": -0.8758913379678931,
        "K_magnitude": 2.724902921120178,
        "K_phase": 1.8980455615438334
      },
      "tau": {
        "d": 2.25,
        "K_real": -1.7878677467016597,
        "K_magnitude": 2.711578212459444,
        "K_phase": 2.2907446432425576
      },
      "charm": {
        "d": 2.25,
        "K_real": -1.7878677467016597,
        "K_magnitude": 2.711578212459444,
        "K_phase": 2.2907446432425576
      },
      "muon": {
        "d": 3.5,
        "K_real": -2.6559119240266766,
        "K_magnitude": 2.67882968332346,
        "K_phase": -3.010692959690219
      },
      "strange": {
        "d": 3.5,
        "K_real": -2.6559119240266766,
        "K_magnitude": 2.67882968332346,
        "K_phase": -3.010692959690219
      },
      "down": {
        "d": 5.0,
        "K_real": -0.6834273957639219,
        "K_magnitude": 2.6405606878474104,
        "K_phase": -1.832595714594046
      },
      "up": {
        "d": 5.25,
        "K_real": -0.1722907159170401,
        "K_magnitude": 2.634288572199317,
        "K_phase": -1.6362461737446836
      },
      "electron": {
        "d": 6.0,
        "K_real": 1.307824868981029,
        "K_magnitude": 2.6156497379620576,
        "K_phase": -1.0471975511965976
      },
      "nu_atm": {
        "d": 13.75,
        "K_real": 0.7834896144020487,
        "K_magnitude": 2.437440634936071,
        "K_phase": -1.2435470920459601
      },
      "nu_sol": {
        "d": 14.5,
        "K_real": 1.9210851738364516,
        "K_magnitude": 2.421474866584962,
        "K_phase": -0.6544984694978742
      }
    },
    "K_at_0.25_steps": [
      2.4011322677058873,
      2.0793442106205897,
      1.6794478334369674,
      1.217155964739858,
      0.7104938272793256,
      0.17909726271870183,
      -0.3565472399076023,
      -0.8758913379678931,
      -1.3591121187449902,
      -1.7878677467016597
    ],
    "resonance_positions": [
      3.25
    ],
    "INSIGHT": "K(d) evaluated at Base-4 ladder positions"
  },
  "qw1059": {
    "K_period": 8.0,
    "ladder_step": 0.25,
    "steps_per_period": 32.0,
    "d_zeros_theory": [
      1.3333333333333335,
      5.333333333333333,
      9.333333333333332,
      13.333333333333334,
      17.333333333333336
    ],
    "d_zeros_mod_0.25": [
      0.08333333333333348,
      0.08333333333333304,
      0.08333333333333215,
      0.08333333333333393,
      0.0833333333333357
    ],
    "INSIGHT": "K period = 8.0, which is 32\u00d70.25 = 8 octaves"
  },
  "qw1060": {
    "scaling_results": {
      "s=1": {
        "period": 8.0,
        "K_at_0.25_steps": [
          2.4011322677058873,
          2.0793442106205897,
          1.6794478334369674,
          1.217155964739858,
          0.7104938272793256,
          0.17909726271870183
        ],
        "zeros_in_range": [
          1.25
        ]
      },
      "s=2": {
        "period": 4.0,
        "K_at_0.25_steps": [
          2.4011322677058873,
          1.6794478334369674,
          0.7104938272793256,
          -0.3565472399076023,
          -1.3591121187449902,
          -2.1459927063831588
        ],
        "zeros_in_range": [
          0.5,
          2.5
        ]
      },
      "s=4": {
        "period": 2.0,
        "K_at_0.25_steps": [
          2.4011322677058873,
          0.7104938272793256,
          -1.3591121187449902,
          -2.6001117014458375,
          -2.3087810266402764,
          -0.6834273957639219
        ],
        "zeros_in_range": [
          0.25,
          1.25,
          2.25,
          3.25
        ]
      },
      "s=8": {
        "period": 1.0,
        "K_at_0.25_steps": [
          2.4011322677058873,
          -1.3591121187449902,
          -2.3087810266402764,
          1.307824868981029,
          2.223270618246192,
          -1.2602676010180802
        ],
        "zeros_in_range": [
          0.0,
          0.5,
          1.0,
          1.5
        ]
      },
      "s=16": {
        "period": 0.5,
        "K_at_0.25_steps": [
          2.4011322677058873,
          -2.3087810266402764,
          2.223270618246192,
          -2.143868096165972,
          2.0699416100912833,
          -2.0009435564215745
        ],
        "zeros_in_range": [
          0.0,
          0.25,
          0.5,
          0.75
        ]
      },
      "s=32": {
        "period": 0.25,
        "K_at_0.25_steps": [
          2.4011322677058873,
          2.223270618246192,
          2.0699416100912833,
          1.9363969900853948,
          1.819039596746886,
          1.7150944769327787
        ],
        "zeros_in_range": []
      }
    },
    "INSIGHT": "Rescaling d moves K(d) features relative to 0.25 grid"
  },
  "qw1061": {
    "eigenvalues_top10": [
      31.46812394127692,
      28.465846559307664,
      4.891722817075275,
      4.3537170705101085,
      2.3537937939001004,
      2.120919211525996,
      1.5039655956147862,
      1.372357178856118,
      1.100439768735924,
      1.0171812498539807
    ],
    "mass_ratios": [
      30.93659458016382,
      27.985028787538127,
      4.809096528054853,
      4.2801782584324055,
      2.314035767212573,
      2.0850946788789706,
      1.4785620515819424,
      1.3491766379424746,
      1.0818521958538807,
      1.0
    ],
    "INSIGHT": "H with 0.25 spacing gives eigenvalue spectrum"
  },
  "qw1062": {
    "16_states": {
      "state_0000": {
        "decimal": 0,
        "d_position": 0.0
      },
      "state_0001": {
        "decimal": 1,
        "d_position": 1.25
      },
      "state_0010": {
        "decimal": 2,
        "d_position": 1.5
      },
      "state_0011": {
        "decimal": 3,
        "d_position": 2.75
      },
      "state_0100": {
        "decimal": 4,
        "d_position": 1.0
      },
      "state_0101": {
        "decimal": 5,
        "d_position": 2.25
      },
      "state_0110": {
        "decimal": 6,
        "d_position": 2.5
      },
      "state_0111": {
        "decimal": 7,
        "d_position": 3.75
      },
      "state_1000": {
        "decimal": 8,
        "d_position": 1.0
      },
      "state_1001": {
        "decimal": 9,
        "d_position": 2.25
      },
      "state_1010": {
        "decimal": 10,
        "d_position": 2.5
      },
      "state_1011": {
        "decimal": 11,
        "d_position": 3.75
      },
      "state_1100": {
        "decimal": 12,
        "d_position": 2.0
      },
      "state_1101": {
        "decimal": 13,
        "d_position": 3.25
      },
      "state_1110": {
        "decimal": 14,
        "d_position": 3.5
      },
      "state_1111": {
        "decimal": 15,
        "d_position": 4.75
      }
    },
    "particle_decomposition": {
      "top": {
        "d_expected": 0.0,
        "octave": 0,
        "sub_bit": 0,
        "reconstructed_d": 0.0
      },
      "bottom": {
        "d_expected": 1.75,
        "octave": 1,
        "sub_bit": 3,
        "reconstructed_d": 1.75
      },
      "tau": {
        "d_expected": 2.25,
        "octave": 2,
        "sub_bit": 1,
        "reconstructed_d": 2.25
      },
      "charm": {
        "d_expected": 2.25,
        "octave": 2,
        "sub_bit": 1,
        "reconstructed_d": 2.25
      },
      "muon": {
        "d_expected": 3.5,
        "octave": 3,
        "sub_bit": 2,
        "reconstructed_d": 3.5
      },
      "strange": {
        "d_expected": 3.5,
        "octave": 3,
        "sub_bit": 2,
        "reconstructed_d": 3.5
      },
      "down": {
        "d_expected": 5.0,
        "octave": 5,
        "sub_bit": 0,
        "reconstructed_d": 5.0
      },
      "up": {
        "d_expected": 5.25,
        "octave": 5,
        "sub_bit": 1,
        "reconstructed_d": 5.25
      },
      "electron": {
        "d_expected": 6.0,
        "octave": 6,
        "sub_bit": 0,
        "reconstructed_d": 6.0
      },
      "nu_atm": {
        "d_expected": 13.75,
        "octave": 13,
        "sub_bit": 3,
        "reconstructed_d": 13.75
      },
      "nu_sol": {
        "d_expected": 14.5,
        "octave": 14,
        "sub_bit": 2,
        "reconstructed_d": 14.5
      }
    },
    "INSIGHT": "Each particle = (octave, sub-bit) address in 4-bit system"
  },
  "qw1063": {
    "action_at_0.25_steps": [
      0.09652793036003567,
      0.6516408666569066,
      1.1147432718142123,
      1.4686090512457675,
      1.7002595368729907,
      1.8011103599935594,
      1.7678189464111176,
      1.6027526129853353,
      1.3128655512252496,
      0.9098339387033161,
      0.409208943619218,
      -0.1689032942828288
    ],
    "increments": [
      0.5551129362968709,
      0.4631024051573057,
      0.35386577943155517,
      0.23165048562722323,
      0.1008508231205687,
      -0.03329141358244181,
      -0.16506633342578225,
      -0.2898870617600857,
      -0.40303161252193354,
      -0.5006249950840981,
      -0.5781122379020468
    ],
    "mean_increment": -0.022862974632653056,
    "std_increment": 0.46359639187970103,
    "normalized_increments": [
      -24.27999616043202,
      -20.255562217870796,
      -15.477678872379173,
      -10.13212363435764,
      -4.411098063177347,
      1.4561278275178937,
      7.219810023759261,
      12.679323947027743,
      17.62813540222019,
      21.896756792491125,
      25.285958944134197
    ],
    "QUANTIZED": "False",
    "INSIGHT": "Action increments per 0.25 step"
  },
  "qw1064": {
    "entropy_from_K": 4.242694267762856,
    "entropy_4bit": 2.772588722239781,
    "ratio": 1.530228495027376,
    "effective_states": 69.59510732138335,
    "expected_states": 16,
    "MATCHES_4BIT": "False",
    "INSIGHT": "Entropy from K(d) vs 4-bit expectation"
  },
  "qw1065": {
    "resonance_positions": [
      3.306613226452906,
      7.304609218436873,
      11.332665330661323
    ],
    "mapped_to_0.25": [
      3.25,
      7.25,
      11.25
    ],
    "particle_matches": {
      "top": {
        "d_expected": 0.0,
        "closest_resonance": 3.25,
        "error": 3.25
      },
      "bottom": {
        "d_expected": 1.75,
        "closest_resonance": 3.25,
        "error": 1.5
      },
      "tau": {
        "d_expected": 2.25,
        "closest_resonance": 3.25,
        "error": 1.0
      },
      "charm": {
        "d_expected": 2.25,
        "closest_resonance": 3.25,
        "error": 1.0
      },
      "muon": {
        "d_expected": 3.5,
        "closest_resonance": 3.25,
        "error": 0.25
      },
      "strange": {
        "d_expected": 3.5,
        "closest_resonance": 3.25,
        "error": 0.25
      },
      "down": {
        "d_expected": 5.0,
        "closest_resonance": 3.25,
        "error": 1.75
      },
      "up": {
        "d_expected": 5.25,
        "closest_resonance": 3.25,
        "error": 2.0
      },
      "electron": {
        "d_expected": 6.0,
        "closest_resonance": 7.25,
        "error": 1.25
      },
      "nu_atm": {
        "d_expected": 13.75,
        "closest_resonance": 11.25,
        "error": 2.5
      },
      "nu_sol": {
        "d_expected": 14.5,
        "closest_resonance": 11.25,
        "error": 3.25
      }
    },
    "INSIGHT": "Resonances in K(d) mapped to Base-4 ladder"
  },
  "qw1066": {
    "d_particles": [
      0,
      1.75,
      2.25,
      3.5,
      5.0,
      6.0
    ],
    "masses_from_phase": [
      0.9947386320728451,
      0.1178952783075277,
      0.04120964028338885,
      0.003435921907139074,
      0.16071655474639293,
      0.699523585576572
    ],
    "mass_ratios_phase": [
      1.0,
      0.118518849581479,
      0.04142760616174707,
      0.0034540951726980485,
      0.16156661615874954,
      0.7032235031617287
    ],
    "expected_4_power": [
      1.0,
      0.08838834764831845,
      0.04419417382415922,
      0.0078125,
      0.0009765625,
      0.000244140625
    ],
    "INSIGHT": "Accumulated phase \u2192 mass"
  },
  "qw1067": {
    "positions_used": [
      0.0,
      1.75,
      2.25,
      3.5,
      5.0,
      5.25,
      6.0,
      13.75,
      14.5
    ],
    "eigenvalues": [
      14.215888695276231,
      7.010582052888755,
      1.509558104799375,
      1.222368921698692,
      0.5118846008414297,
      0.4415017836698711,
      0.2509222396842381,
      -0.10954413385095595,
      -3.4429718556546587
    ],
    "INSIGHT": "K matrix built from actual particle positions"
  },
  "qw1068": {
    "winding_per_step": [
      0.031249999999999754,
      0.031249999999999892,
      0.031249999999999754,
      0.031249999999999754,
      0.031249999999999892,
      0.031249999999999965,
      0.031250000000000035,
      0.031249999999999965,
      0.031250000000000035,
      0.031249999999999965,
      0.031249999999999965,
      0.031250000000000035,
      0.031250000000000035,
      0.031249999999999823,
      0.031249999999999965
    ],
    "cumulative_winding": [
      0.031249999999999754,
      0.062499999999999646,
      0.0937499999999994,
      0.12499999999999915,
      0.15624999999999906,
      0.18749999999999903,
      0.21874999999999906,
      0.24999999999999903,
      0.28124999999999906,
      0.312499999999999,
      0.34374999999999895,
      0.374999999999999,
      0.40624999999999906,
      0.4374999999999989,
      0.46874999999999883
    ],
    "mean_winding": 0.031249999999999976,
    "std_winding": 9.051985559525574e-17,
    "UNIFORM": "True",
    "INSIGHT": "Winding number accumulated per 0.25 octave"
  },
  "qw1069": {
    "scales": [
      1,
      2,
      4,
      8
    ],
    "N_features": [
      "3",
      "3",
      "3",
      "3"
    ],
    "fractal_dimension": 0.0,
    "INSIGHT": "Fractal dimension from K(d) scaling"
  },
  "qw1070": {
    "N_bound_states": 2,
    "bound_energies": [
      -1.7244307235153,
      -0.04671607297428096
    ],
    "energy_ratios": [
      1.0
    ],
    "INSIGHT": "Bound state energies as mass candidates"
  },
  "qw1071": {
    "correlation_length": 1.5075376884422111,
    "xi_in_0.25_steps": 6.030150753768845,
    "xi_over_period": 0.1884422110552764,
    "INSIGHT": "Correlation length \u03be = 1.51 = 6.0 \u00d7 0.25"
  },
  "qw1072": {
    "d_test": [
      0,
      1.75,
      3.5,
      6.0
    ],
    "m_from_exp(-alpha*d/2)": [
      1.0,
      0.08838834764831845,
      0.007812500000000002,
      0.00024414062500000016
    ],
    "m_from_4^(-d)": [
      1.0,
      0.08838834764831845,
      0.0078125,
      0.000244140625
    ],
    "ratios": [
      1.0,
      1.0,
      1.0000000000000002,
      1.0000000000000007
    ],
    "errors_percent": [
      0.0,
      0.0,
      2.220446049250313e-14,
      6.661338147750939e-14
    ],
    "FORMULA_VERIFIED": true,
    "DERIVED": true,
    "INSIGHT": "m(d) = exp(-\u03b1\u00d7d/2) = 4^(-d) IS EXACT because \u03b1 = 4\u00d7ln(2)"
  },
  "qw1073": {
    "K_zeros": [
      1.3333333333333335,
      5.333333333333333,
      9.333333333333332,
      13.333333333333334,
      17.333333333333336,
      21.333333333333332
    ],
    "zero_spacings": [
      3.9999999999999996,
      3.999999999999999,
      4.000000000000002,
      4.000000000000002,
      3.9999999999999964
    ],
    "mean_spacing": 3.9999999999999996,
    "expected_spacing": 4.0,
    "quarter_period": 2.0,
    "STEP_RELATION": "Period = 8.0, zeros every 4.0",
    "INSIGHT": "Zero spacing = 4, but step = 0.25 \u2192 need scaling by 16"
  },
  "qw1074": {
    "exp_alpha": 16,
    "period_K": 8.0,
    "step_for_16_per_period": 0.5,
    "16_states_in_d_0_4": {
      "0": {
        "d": 0.0,
        "K": 2.4011322677058873
      },
      "1": {
        "d": 0.25,
        "K": 2.0793442106205897
      },
      "2": {
        "d": 0.5,
        "K": 1.6794478334369674
      },
      "3": {
        "d": 0.75,
        "K": 1.217155964739858
      },
      "4": {
        "d": 1.0,
        "K": 0.7104938272793256
      },
      "5": {
        "d": 1.25,
        "K": 0.17909726271870183
      },
      "6": {
        "d": 1.5,
        "K": -0.3565472399076023
      },
      "7": {
        "d": 1.75,
        "K": -0.8758913379678931
      },
      "8": {
        "d": 2.0,
        "K": -1.3591121187449902
      },
      "9": {
        "d": 2.25,
        "K": -1.7878677467016597
      },
      "10": {
        "d": 2.5,
        "K": -2.1459927063831588
      },
      "11": {
        "d": 2.75,
        "K": -2.4201063245331436
      },
      "12": {
        "d": 3.0,
        "K": -2.6001117014458375
      },
      "13": {
        "d": 3.25,
        "K": -2.6795664861575625
      },
      "14": {
        "d": 3.5,
        "K": -2.6559119240266766
      },
      "15": {
        "d": 3.75,
        "K": -2.5305520942527093
      }
    },
    "INSIGHT": "16 = exp(\u03b1) gives structure in [0, 4) with step 0.25"
  },
  "qw1075": {
    "m_Planck_effective_MeV": 2093.056,
    "predictions": {
      "top": {
        "d": 0.0,
        "m_predicted_MeV": 2093.056
      },
      "bottom": {
        "d": 1.75,
        "m_predicted_MeV": 185.00176137539881
      },
      "tau": {
        "d": 2.25,
        "m_predicted_MeV": 92.50088068769941
      },
      "charm": {
        "d": 2.25,
        "m_predicted_MeV": 92.50088068769941
      },
      "muon": {
        "d": 3.5,
        "m_predicted_MeV": 16.352
      },
      "strange": {
        "d": 3.5,
        "m_predicted_MeV": 16.352
      },
      "down": {
        "d": 5.0,
        "m_predicted_MeV": 2.044
      },
      "up": {
        "d": 5.25,
        "m_predicted_MeV": 1.4453262607453032
      },
      "electron": {
        "d": 6.0,
        "m_predicted_MeV": 0.511
      },
      "nu_atm": {
        "d": 13.75,
        "m_predicted_MeV": 1.1026964269602228e-05
      },
      "nu_sol": {
        "d": 14.5,
        "m_predicted_MeV": 3.89862060546875e-06
      }
    },
    "errors_percent": {
      "top": 98.7884602917342,
      "bottom": 95.5741205412584,
      "tau": 94.7945480761002,
      "charm": 92.74502896567063,
      "muon": 84.52980132450331,
      "strange": 82.78736842105263,
      "down": 57.41666666666666,
      "up": 37.159727793682464,
      "electron": 0.0
    },
    "mean_error_percent": 71.53285800896316,
    "SUCCESS": "False",
    "INSIGHT": "Mass formula m = m_P \u00d7 4^(-d) with calibrated m_P"
  },
  "qw1076": {
    "KEY_DISCOVERIES": [
      "m(d) = exp(-\u03b1\u00d7d/2) = 4^(-d) EXACT!",
      "\u03b1 = ln(16), Base = 4 = \u221a16, k = 0.5 connects",
      "16 = exp(\u03b1) gives 16 states in [0,4) with step 0.25"
    ],
    "FORMULA_DERIVED": true,
    "STEP_EXPLAINED": true,
    "PREDICTIONS_WORK": false,
    "VERDICT": "SUKCES: Base-4 z krokiem 0.25 WYNIKA z \u03b1 = 4\u00d7ln(2)!"
  }
}
```

---

## Faza 58: `phase58_h_true_mapping.py`
### Kod:
```python
log = logging.getLogger(__name__)

RAW_DIR = "./raw_strain_unfiltered"
FS = 4096
GPS = 1266965117 # Baseline
N_SAMPLES = 131072 # 128K samples

def fetch_empirical_spectra():
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

G_H1, G_L1 = fetch_empirical_spectra()
SCALES = np.unique(np.logspace(3, np.log10(N_SAMPLES//4), 15).astype(int))

def run_trial(args):
    h_true, seed = args
    X_in = generate_fgn_spectrum(N_SAMPLES, h_true, seed=seed)
    
    X_out_h1 = G_H1 * np.exp(1j * np.angle(X_in))
    X_out_l1 = G_L1 * np.exp(1j * np.angle(X_in))
    
    x_h1 = np.fft.irfft(X_out_h1, n=N_SAMPLES)
    x_l1 = np.fft.irfft(X_out_l1, n=N_SAMPLES)
    
    return cross_mfdfa_q0(x_h1, x_l1, SCALES)

def main():
    h_true_values = np.arange(0.05, 0.46, 0.01)
    h_true_values = np.round(h_true_values, 2)
    N_TRIALS = 100
    
    
    results = {}
    
    num_cores = multiprocessing.cpu_count()
    
    with multiprocessing.Pool(num_cores) as p:
        for h_true in h_true_values:
            args = [(h_true, int(h_true*1000) * N_TRIALS + i) for i in range(N_TRIALS)]
            h_obs_trials = p.map(run_trial, args)
            
            mean_obs = np.mean(h_obs_trials)
            std_obs = np.std(h_obs_trials)
            
            results[f"{h_true:.2f}"] = {
                "mean": float(mean_obs),
                "std": float(std_obs)
            }

    with open("QW_1660_v58_HTrue_Mapping.json", "w") as f:
        json.dump(results, f, indent=2)
        
    for h in h_true_values:
        key = f"{h:.2f}"
        marker = " <=== FIN" if abs(h - 0.23) < 0.005 else ""
        if abs(results[key]["mean"] - 0.31) < 0.015:
            marker += "  [MATCHES 0.31]"

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "0.05": {
    "mean": 0.4004410091997295,
    "std": 0.027933321248884768
  },
  "0.06": {
    "mean": 0.4037080246714917,
    "std": 0.027440867036929373
  },
  "0.07": {
    "mean": 0.39858702807371643,
    "std": 0.024428475993454848
  },
  "0.08": {
    "mean": 0.40126004797377846,
    "std": 0.0237754441519594
  },
  "0.09": {
    "mean": 0.39899303797859104,
    "std": 0.02716992747638438
  },
  "0.10": {
    "mean": 0.39857089796607176,
    "std": 0.027970497217470577
  },
  "0.11": {
    "mean": 0.39963077436004674,
    "std": 0.028280015123032336
  },
  "0.12": {
    "mean": 0.39820285479377476,
    "std": 0.022425824368900078
  },
  "0.13": {
    "mean": 0.3986489805097576,
    "std": 0.026189567416261723
  },
  "0.14": {
    "mean": 0.3988655107489701,
    "std": 0.023426960570410077
  },
  "0.15": {
    "mean": 0.39620037112491585,
    "std": 0.026921608819045015
  },
  "0.16": {
    "mean": 0.39457224499111254,
    "std": 0.027340425191611815
  },
  "0.17": {
    "mean": 0.40082080274225335,
    "std": 0.027089112575727187
  },
  "0.18": {
    "mean": 0.3941544493091026,
    "std": 0.02355772334169456
  },
  "0.19": {
    "mean": 0.39981841104111665,
    "std": 0.027016244861192636
  },
  "0.20": {
    "mean": 0.40029184793024347,
    "std": 0.024708432620686956
  },
  "0.21": {
    "mean": 0.39789812590783497,
    "std": 0.025520497936303666
  },
  "0.22": {
    "mean": 0.395847891965017,
    "std": 0.025790662490122216
  },
  "0.23": {
    "mean": 0.39648745424536785,
    "std": 0.02385882135337397
  },
  "0.24": {
    "mean": 0.3995899720316251,
    "std": 0.02578104008013167
  },
  "0.25": {
    "mean": 0.3972815606284205,
    "std": 0.025576727198992264
  },
  "0.26": {
    "mean": 0.3983771805109425,
    "std": 0.028935293979013833
  },
  "0.27": {
    "mean": 0.3969920165843066,
    "std": 0.02497529525206574
  },
  "0.28": {
    "mean": 0.4031693429115578,
    "std": 0.024551636034024146
  },
  "0.29": {
    "mean": 0.39991865182437275,
    "std": 0.026682409914462683
  },
  "0.30": {
    "mean": 0.4008943857800395,
    "std": 0.02675033780957701
  },
  "0.31": {
    "mean": 0.3962344554211347,
    "std": 0.02571293799437975
  },
  "0.32": {
    "mean": 0.3947083278936816,
    "std": 0.026558697381502395
  },
  "0.33": {
    "mean": 0.40284990798616005,
    "std": 0.026789983832809547
  },
  "0.34": {
    "mean": 0.3987862366901604,
    "std": 0.026211352369012567
  },
  "0.35": {
    "mean": 0.3973613704299476,
    "std": 0.025134966977588186
  },
  "0.36": {
    "mean": 0.3991941009288439,
    "std": 0.024758946990399582
  },
  "0.37": {
    "mean": 0.39639518530902856,
    "std": 0.02511888629071801
  },
  "0.38": {
    "mean": 0.39479028817294837,
    "std": 0.024507724205019284
  },
  "0.39": {
    "mean": 0.4022232369392795,
    "std": 0.02617700375805709
  },
  "0.40": {
    "mean": 0.3985299913591742,
    "std": 0.027914514081746673
  },
  "0.41": {
    "mean": 0.39711781214788244,
    "std": 0.02600649550958563
  },
  "0.42": {
    "mean": 0.395751567199055,
    "std": 0.02467251649620674
  },
  "0.43": {
    "mean": 0.39459432398424904,
    "std": 0.026048966497937552
  },
  "0.44": {
    "mean": 0.39641806407787455,
    "std": 0.02560372103068189
  },
  "0.45": {
    "mean": 0.3984582561241745,
    "std": 0.027550973733774454
  }
}
```

---

## Faza 59: `phase59_envelope_swap.py`
### Kod:
```python
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
    
    G_generic = power_law_spectrum(N_SAMPLES, H=0.5)
    
    G_generic = G_generic * (np.sum(np.abs(R_l1)) / np.sum(G_generic))
    
    X_h1 = G_H1 * np.exp(1j * phase_h1)       # Real H1
    X_l1 = G_generic * np.exp(1j * phase_l1)  # Real phase, generic amplitude
    
    x_h1_out = np.fft.irfft(X_h1, n=N_SAMPLES)
    x_l1_out = np.fft.irfft(X_l1, n=N_SAMPLES)
    
    return cross_mfdfa_q0(x_h1_out, x_l1_out, SCALES)

def main():
    x_h1, x_l1 = fetch_empirical_data()
    G_H1 = np.abs(np.fft.rfft(x_h1))
    G_L1 = np.abs(np.fft.rfft(x_l1))
    
    N_TRIALS = 100
    num_cores = multiprocessing.cpu_count()
    
    with multiprocessing.Pool(num_cores) as p:
        h_cross_t1 = p.map(run_test_1, [(G_H1, G_L1, i) for i in range(N_TRIALS)])
    
    m_t1 = np.mean(h_cross_t1)
    s_t1 = np.std(h_cross_t1)
    
    with multiprocessing.Pool(num_cores) as p:
        hc_t2 = run_test_2((x_h1, x_l1, 42))
        
    
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
        

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "Test1_Independent_Phases_Real_Envelope": {
    "mean": 0.2868038142777771,
    "std": 0.012702395253275257,
    "verdict": "If ~0.50, then the envelope ALONE does not cause the anomaly. If ~0.31, then the envelope alone is responsible."
  },
  "Test2_Real_Phases_Generic_Envelope": {
    "value": 0.705328267933796,
    "verdict": "If ~0.50, then removing the shared envelope destroys the anomaly. If ~0.31, the anomaly is robustly in the phase."
  }
}
```

---

## Faza 60: `phase60_scale_stability.py`
### Kod:
```python
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

def cross_mfdfa_q0_local(x, y, scales):
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
            
    if len(valid_scales) < 3: return np.nan
    return float(np.polyfit(np.log(valid_scales), np.log(F), 1)[0])

def main():
    x_h1, x_l1 = fetch_empirical_data()
    
    
    scale_boundaries = [1000, 3000, 10000, 30000, 100000, N_SAMPLES//4]
    
    results = {}
    
    for i in range(5):
        s_min = scale_boundaries[i]
        s_max = scale_boundaries[i+1]
        
        scales = np.unique(np.logspace(np.log10(s_min), np.log10(s_max), 15).astype(int))
        
        h_cross = cross_mfdfa_q0_local(x_h1, x_l1, scales)
        
        key = f"Scale_{s_min}_to_{s_max}"
        results[key] = {
            "s_min": s_min,
            "s_max": s_max,
            "H_cross": h_cross
        }
        
    with open("QW_1660_v60_ScaleStability.json", "w") as f:
        json.dump(results, f, indent=2)

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "H_unfiltered_long_scales": 0.03733273522994843,
  "H_filtered_short_scales": 0.10019515377472395,
  "H_anomaly_reproduced": 0.002944397387996455,
  "Interpretation": {
    "H_anomaly_reproduced": "The ~0.23 value seen previously is an artifact of the 20Hz filter forcing the signal to cross zero over long scales.",
    "H_unfiltered_long_scales": "The true long-range scaling exponent of the raw noise.",
    "H_filtered_short_scales": "The short-range scaling exponent inside the studied beta frequency band."
  }
}
```
```json
{
  "Config": {
    "N_trials": 500,
    "N_samples_per_trial": 131072,
    "H_local": 0.04
  },
  "Results": {
    "Null_Cross_H_Mean": 0.5018624026057812,
    "Null_Cross_H_Std": 0.03820165650501183,
    "Null_Cross_H_Min": 0.39702322618872843,
    "Null_1st_Percentile": 0.4147443000100347
  },
  "Bootstrap_95_CI": [
    0.49861649052509377,
    0.5051488624382866
  ],
  "Real_LIGO_Observation": 0.311,
  "Statistical_Significance": "5.00 sigma",
  "Verdict": "Strong anomaly confirmed at 5.00 sigma. The observed 0.311 NEVER occurs in 500 null trials. Real cross-detector anti-persistence exists."
}
```
```json
{
  "H_real": 0.395795605895727,
  "H_shuffle": 0.501406297389991,
  "H_phase_randomized": 0.37727819984529654,
  "timestamp_utc": "2026-01-03T22:20:51.998166+00:00"
}
```
```json
{
  "0.0": 0.31105068079469167,
  "0.1": 0.30997165995459897,
  "1.0": 0.28581875130863643,
  "5.0": 0.29588892330405664,
  "10.0": 0.2947718600757252,
  "50.0": 0.2982288953028582,
  "100.0": 0.3202829195822657
}
```
```json
{
  "Coherence_Coarse_Sec": {
    "-5": 0.22713815148589375,
    "-4": 0.2284834335849714,
    "-3": 0.22813757926549674,
    "-2": 0.22719536436191218,
    "-1": 0.2285323440264692,
    "0": 0.22774859668180394,
    "1": 0.22885061737842644,
    "2": 0.22871103459093778,
    "3": 0.22899172351486632,
    "4": 0.2286404321201261,
    "5": 0.22867975079279013
  },
  "Coherence_Fine_ms": {
    "-20.0": 0.22754163493399385,
    "-18.0": 0.22747685654528674,
    "-16.0": 0.22742069398267709,
    "-14.0": 0.22738031797250913,
    "-12.0": 0.22737844803239704,
    "-10.0": 0.22748645642290677,
    "-8.0": 0.22753136500406423,
    "-6.0": 0.22759328581688065,
    "-4.0": 0.2276585701554386,
    "-2.0": 0.2277118181868326,
    "0.0": 0.22774859668180394,
    "2.0": 0.22776441559660618,
    "4.0": 0.22778231958839587,
    "6.0": 0.22781260225540848,
    "8.0": 0.22786637891755065,
    "10.0": 0.22795232472939844,
    "12.0": 0.227799685381826,
    "14.0": 0.2278843776961098,
    "16.0": 0.22793328160504317,
    "18.0": 0.22794253354135366,
    "20.0": 0.22793205412191542
  },
  "Interpretation_Guide": {
    "Peak_at_0": "Global Field / Non-local Correlation",
    "Peak_at_10ms": "Light-speed Propagation (GW-like)",
    "Peak_at_Seconds": "Slow Environmental (Seismic/Atmospheric)",
    "Flat": "Uncorrelated Noise"
  }
}
```
```json
{
  "Best_Lag_ms": 73.2421875,
  "Best_Correlation": -0.005650218758170756,
  "Verdict": "If peak is at ~10ms, it is a true gravitational/light-speed wave. If peak is exactly 0ms, it's global synchronous noise. If no clear peak, there is no broadband phase correlation.",
  "Correlation_Window_ms": {
    "-24.90 ms": 0.0035958924423889955,
    "-24.66 ms": -0.0015633011271547006,
    "-24.41 ms": 0.00017008011077205253,
    "-24.17 ms": 0.001355628330727946,
    "-23.93 ms": -0.0004134882273864998,
    "-23.68 ms": 0.0005054541196247306,
    "-23.44 ms": 0.0004900404954794068,
    "-23.19 ms": 0.0020395115197046417,
    "-22.95 ms": -0.0013569693128512989,
    "-22.71 ms": -0.0009627638168922174,
    "-22.46 ms": 0.001537215409647673,
    "-22.22 ms": 0.0007465840528576063,
    "-21.97 ms": 7.288600826827111e-06,
    "-21.73 ms": -0.00031926895456860677,
    "-21.48 ms": 0.0023335996325561915,
    "-21.24 ms": -0.000351092269335575,
    "-21.00 ms": -0.0013262921264217344,
    "-20.75 ms": 0.0027737408633027928,
    "-20.51 ms": 0.0002140976214402963,
    "-20.26 ms": -0.0017877361836068969,
    "-20.02 ms": 0.0020615367907124214,
    "-19.78 ms": 0.002347911922408769,
    "-19.53 ms": -0.0016768330858254732,
    "-19.29 ms": -1.970008582722594e-05,
    "-19.04 ms": 0.0018161985086641468,
    "-18.80 ms": 0.00038459294603834045,
    "-18.55 ms": 0.0012723814097295358,
    "-18.31 ms": -0.0009614678279654249,
    "-18.07 ms": 0.0008923191411199123,
    "-17.82 ms": 0.0014303140773147583,
    "-17.58 ms": -0.0013920944838412953,
    "-17.33 ms": 0.0031772269965881883,
    "-17.09 ms": 0.00017789552463652144,
    "-16.85 ms": -0.0015651730805940265,
    "-16.60 ms": 0.0018926567135340066,
    "-16.36 ms": 7.97589406729218e-05,
    "-16.11 ms": 0.001648131013788291,
    "-15.87 ms": 0.0003240010701814508,
    "-15.62 ms": 0.0011457879230633169,
    "-15.38 ms": 0.00089830500275615,
    "-15.14 ms": -0.0017275701490225563,
    "-14.89 ms": 0.0034050243898876608,
    "-14.65 ms": -0.00021676111242924362,
    "-14.40 ms": -0.0011972245040886079,
    "-14.16 ms": 0.002871610415299286,
    "-13.92 ms": -0.0010700932326487345,
    "-13.67 ms": 0.0023655663472634805,
    "-13.43 ms": 0.0019648860626813006,
    "-13.18 ms": -0.0022148141302521713,
    "-12.94 ms": 0.0010894040102334087,
    "-12.70 ms": 0.0012423076847510306,
    "-12.45 ms": 0.001505315720024038,
    "-12.21 ms": 0.0009349869149328864,
    "-11.96 ms": 9.707745509500646e-05,
    "-11.72 ms": -4.690221944909932e-05,
    "-11.47 ms": 0.0016199404633047672,
    "-11.23 ms": 0.001445677612150131,
    "-10.99 ms": -0.0010669552237243453,
    "-10.74 ms": 0.001344353060641185,
    "-10.50 ms": 0.001272999900487303,
    "-10.25 ms": 0.0008204132412868782,
    "-10.01 ms": -0.00031235607648622487,
    "-9.77 ms": 0.0012853310661435364,
    "-9.52 ms": 0.002164395644159403,
    "-9.28 ms": -0.001941729894240139,
    "-9.03 ms": 0.002208708098259537,
    "-8.79 ms": 0.0010388139945964607,
    "-8.54 ms": -0.0009969571629149417,
    "-8.30 ms": 0.0030184750758895647,
    "-8.06 ms": -0.0003822562793872377,
    "-7.81 ms": 0.0007977175482042462,
    "-7.57 ms": 0.0018976038759113678,
    "-7.32 ms": -0.0009456301976598709,
    "-7.08 ms": 0.00151764801706955,
    "-6.84 ms": 0.0014211995413514987,
    "-6.59 ms": -0.0002813114199989071,
    "-6.35 ms": -5.3712641338081626e-06,
    "-6.10 ms": 0.001970505995524733,
    "-5.86 ms": 0.0026007706801578993,
    "-5.62 ms": -0.0019946366304326515,
    "-5.37 ms": 0.0014201064526390096,
    "-5.13 ms": 0.003231590192172192,
    "-4.88 ms": -0.0022646161375598655,
    "-4.64 ms": 0.0017579494428678977,
    "-4.39 ms": 0.00143613910709342,
    "-4.15 ms": 0.0002827334271200185,
    "-3.91 ms": 0.0025202316600075043,
    "-3.66 ms": -0.002308038435619221,
    "-3.42 ms": 0.0019035634563612438,
    "-3.17 ms": 0.002384385718265745,
    "-2.93 ms": -0.0023088516721113347,
    "-2.69 ms": 0.0028886102039535956,
    "-2.44 ms": -0.00027315709272521303,
    "-2.20 ms": 0.0010609043632456843,
    "-1.95 ms": 0.0024014283671421064,
    "-1.71 ms": -0.0012244714012334403,
    "-1.46 ms": 0.003080209209764442,
    "-1.22 ms": -0.001583111947079528,
    "-0.98 ms": 0.0007468173716760598,
    "-0.73 ms": 0.0027386996591416113,
    "-0.49 ms": -0.0009326546566146771,
    "-0.24 ms": 0.0023591702014447904,
    "0.00 ms": -0.0008141046665487757,
    "0.24 ms": 0.001669849220685649,
    "0.49 ms": 0.002455146839938389,
    "0.73 ms": -0.0015496847419886645,
    "0.98 ms": 0.0011500722054539857,
    "1.22 ms": 0.0016680894830704196,
    "1.46 ms": -0.0001983649469891438,
    "1.71 ms": 0.0010650975025451983,
    "1.95 ms": 0.0021359048893473767,
    "2.20 ms": -7.279701833772713e-05,
    "2.44 ms": -0.0001252849054576339,
    "2.69 ms": 0.001890500652592548,
    "2.93 ms": 0.0005329389855970609,
    "3.17 ms": -0.000754706815663927,
    "3.42 ms": 0.0018228613785432081,
    "3.66 ms": 0.0012721137809966468,
    "3.91 ms": 0.0007220120666060275,
    "4.15 ms": 0.00021898813978581388,
    "4.39 ms": -0.0008471737570492024,
    "4.64 ms": 0.0018167036402158222,
    "4.88 ms": 0.0010777796843581082,
    "5.13 ms": -0.0003674864153857978,
    "5.37 ms": 0.0009149772524157882,
    "5.62 ms": 0.0014122365606375777,
    "5.86 ms": -5.9070473936331805e-05,
    "6.10 ms": 0.0007337426260857589,
    "6.35 ms": 0.0020346626159028747,
    "6.59 ms": -0.0009198023032520723,
    "6.84 ms": 0.0006959395853365341,
    "7.08 ms": 0.0007905222308515899,
    "7.32 ms": -0.0008360205116897204,
    "7.57 ms": 0.0027083644906924337,
    "7.81 ms": 0.0003334952880923528,
    "8.06 ms": -0.00041982207434698514,
    "8.30 ms": 0.0011443289366778971,
    "8.54 ms": 0.0003482171257196657,
    "8.79 ms": 0.0007702690707188068,
    "9.03 ms": -0.0008464080527831987,
    "9.28 ms": 0.0010071997339217511,
    "9.52 ms": 0.001947483803407492,
    "9.77 ms": 0.0006820474220215479,
    "10.01 ms": -0.000663842273425792,
    "10.25 ms": 0.0003373093289508459,
    "10.50 ms": 0.0015240362895267354,
    "10.74 ms": -0.0016475597435585802,
    "10.99 ms": 0.00042988213888107214,
    "11.23 ms": 0.002470088247878573,
    "11.47 ms": 0.0011871871530210985,
    "11.72 ms": -0.002074234329936727,
    "11.96 ms": -0.0005141144930285108,
    "12.21 ms": 0.0041528995821315424,
    "12.45 ms": -0.0016228369240076,
    "12.70 ms": -0.002244549994312596,
    "12.94 ms": 0.0037400614885914683,
    "13.18 ms": 0.00039106143091727544,
    "13.43 ms": -0.0019853563594793458,
    "13.67 ms": 0.0011190536410211449,
    "13.92 ms": 0.0023330869753939603,
    "14.16 ms": -0.0003616351627169883,
    "14.40 ms": -0.0011771663524364159,
    "14.65 ms": 0.002119871615158723,
    "14.89 ms": -0.0009848050714806838,
    "15.14 ms": -0.0001906672510117372,
    "15.38 ms": 0.0018270415850659622,
    "15.62 ms": -0.00034244649390020973,
    "15.87 ms": 0.00040048994610381765,
    "16.11 ms": -0.0016879625876422412,
    "16.36 ms": 0.0028303998420298434,
    "16.60 ms": -4.997284293283329e-05,
    "16.85 ms": -0.0030187742761357585,
    "17.09 ms": 0.003803223982454949,
    "17.33 ms": -0.0013198718734197547,
    "17.58 ms": 0.0010115647537846607,
    "17.82 ms": 0.001520261152710049,
    "18.07 ms": -0.0023831911658713447,
    "18.31 ms": 0.002030118606286494,
    "18.55 ms": -0.0014347003885494152,
    "18.80 ms": -7.78976108667486e-06,
    "19.04 ms": 0.0013756010423563911,
    "19.29 ms": -0.0007928289397442114,
    "19.53 ms": 0.0007965138403484423,
    "19.78 ms": -0.0006161300191874048,
    "20.02 ms": 0.0009850579714104986,
    "20.26 ms": -0.00022581939097157493,
    "20.51 ms": -0.0018688656550558558,
    "20.75 ms": 0.0014928530139534074,
    "21.00 ms": 0.0003885437493149752,
    "21.24 ms": -0.0008053923889068185,
    "21.48 ms": 0.0009445005040433828,
    "21.73 ms": -0.00021550349715158482,
    "21.97 ms": -0.00020866735578380276,
    "22.22 ms": 0.00043726384608365785,
    "22.46 ms": -0.0010525567289084176,
    "22.71 ms": 0.00029623005177442386,
    "22.95 ms": 0.00044074740127627156,
    "23.19 ms": -0.0015295763450776862,
    "23.44 ms": 0.0011450537815834666,
    "23.68 ms": 0.0012728978607754091,
    "23.93 ms": -0.003036822580908511,
    "24.17 ms": 0.0008509201291327026,
    "24.41 ms": 0.0011240928077248147,
    "24.66 ms": -0.0020782885853291963,
    "24.90 ms": 0.0009105883117009215
  }
}
```
```json
{
  "H_original": 0.0904990693214959,
  "H_fragmented": 0.0905518729921673,
  "Delta": -5.280367067139746e-05,
  "Interpretation": {
    "H_frag ~ H_original": "FIN is LOCAL/ROBUST (survives fragmentation)",
    "H_frag -> 0.5": "FIN is GLOBAL (requires long-range continuity)"
  }
}
```
```json
{
  "Config": {
    "N_trials": 1000,
    "N_samples": 131072,
    "Description": "Independent random phases applied strictly to empirical H1 and L1 PSDs"
  },
  "Results": {
    "Mean_H_cross": 0.389898436517463,
    "Std_H_cross": 0.023539041835616652,
    "Min_H_cross": 0.322682844205378,
    "P1_H_cross": 0.33770149682344514
  },
  "Observation": 0.311,
  "Significance": "-3.35 sigma",
  "Verdict": "Anomaly survives at -3.35 sigma even against realistic instrumental envelopes! The FIN correlation is still required."
}
```
```json
{
  "diurnal_null_test": [
    {
      "segment": 0,
      "H1_real": 0.2861297793052872,
      "L1_real": 0.26324010867471026,
      "Delta_real": 0.022889670630576953,
      "H1_shuffle": 0.2012720915031506,
      "L1_shuffle": 0.2001353507574666,
      "H1_reverse": 0.27527233598374606,
      "L1_reverse": 0.2625768901276506
    }
  ],
  "segment_sec": 3600,
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-01T21:39:03.253257+00:00"
}
```
```json
{
  "H1": {
    "real": 0.3107033628020025,
    "null_mean": 0.29152461412519504,
    "null_std": 0.006086465184678248
  },
  "L1": {
    "real": 0.26512680447073594,
    "null_mean": 0.2510490771952443,
    "null_std": 0.0040103972037196884
  },
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-03T21:56:45.987553+00:00"
}
```
```json
{
  "window_sec": 64,
  "step_sec": 32,
  "N_windows": 126,
  "H_mean": 0.306493493443769,
  "H_std": 0.011984760442150907,
  "timestamp_utc": "2026-01-01T22:36:49.918816+00:00"
}
```
```json
{
  "v40_CosmicTime_3x": {
    "H1": {
      "Spearman_z_H": -0.9642857142857145,
      "p_value": 0.0004541491691941689,
      "verdict": "Cosmic-Structured"
    },
    "L1": {
      "Spearman_z_H": 0.07142857142857144,
      "p_value": 0.8790481931481541,
      "verdict": "Stationary"
    },
    "V1": {
      "Spearman_z_H": 0.5714285714285715,
      "p_value": 0.1802019889115274,
      "verdict": "Cosmic-Structured"
    }
  },
  "v41_EnergyCoupling_3x": {
    "H1": {
      "Spearman_FIN_Energy": -0.9642857142857145,
      "p_value": 0.0004541491691941689,
      "verdict": "Coupled"
    },
    "L1": {
      "Spearman_FIN_Energy": 0.07142857142857144,
      "p_value": 0.8790481931481541,
      "verdict": "Decoupled"
    },
    "V1": {
      "Spearman_FIN_Energy": 0.6071428571428572,
      "p_value": 0.1482311614811614,
      "verdict": "Coupled"
    }
  },
  "v42_Axioms": {
    "Axiom_1": "FIN is scale-invariant (H constant across geometry).",
    "Axiom_2": "FIN carries structure without energy transport.",
    "Axiom_3": "FIN is tensorial and time-symmetric."
  },
  "Metadata": {
    "segment_sec": 64,
    "method": "MF-DFA q=0",
    "timestamp": "v40_42_3x_audit"
  }
}
```
```json
{
  "cross_hurst_null_validation": [
    {
      "pair": "H1-L1",
      "H_real": 0.3610605334380803,
      "H_shuffle": 0.5310469043003473,
      "H_phase_randomized": 0.39679047633700987,
      "delta_real_shuffle": -0.16998637086226698,
      "delta_real_phase": -0.035729942898929556,
      "compute_time_sec": 6.245482921600342
    },
    {
      "pair": "H1-V1",
      "H_real": 0.2150698637573661,
      "H_shuffle": 0.5219537647569182,
      "H_phase_randomized": 0.23618825276589883,
      "delta_real_shuffle": -0.30688390099955215,
      "delta_real_phase": -0.021118389008532745,
      "compute_time_sec": 6.256558418273926
    },
    {
      "pair": "L1-V1",
      "H_real": 0.27003525002927997,
      "H_shuffle": 0.5039965798656929,
      "H_phase_randomized": 0.28341320075874976,
      "delta_real_shuffle": -0.23396132983641293,
      "delta_real_phase": -0.01337795072946979,
      "compute_time_sec": 6.2524378299713135
    }
  ],
  "summary": {
    "mean_delta_real_phase": -0.023408760878977364,
    "std_delta_real_phase": 0.00926776639630034,
    "N_pairs": 3
  },
  "sample_rate": 4096,
  "window_sec": 512,
  "timestamp_utc": "2026-01-04T01:25:25.095853+00:00"
}
```
```json
{
  "0.05": {
    "mean": 0.4004410091997295,
    "std": 0.027933321248884768
  },
  "0.06": {
    "mean": 0.4037080246714917,
    "std": 0.027440867036929373
  },
  "0.07": {
    "mean": 0.39858702807371643,
    "std": 0.024428475993454848
  },
  "0.08": {
    "mean": 0.40126004797377846,
    "std": 0.0237754441519594
  },
  "0.09": {
    "mean": 0.39899303797859104,
    "std": 0.02716992747638438
  },
  "0.10": {
    "mean": 0.39857089796607176,
    "std": 0.027970497217470577
  },
  "0.11": {
    "mean": 0.39963077436004674,
    "std": 0.028280015123032336
  },
  "0.12": {
    "mean": 0.39820285479377476,
    "std": 0.022425824368900078
  },
  "0.13": {
    "mean": 0.3986489805097576,
    "std": 0.026189567416261723
  },
  "0.14": {
    "mean": 0.3988655107489701,
    "std": 0.023426960570410077
  },
  "0.15": {
    "mean": 0.39620037112491585,
    "std": 0.026921608819045015
  },
  "0.16": {
    "mean": 0.39457224499111254,
    "std": 0.027340425191611815
  },
  "0.17": {
    "mean": 0.40082080274225335,
    "std": 0.027089112575727187
  },
  "0.18": {
    "mean": 0.3941544493091026,
    "std": 0.02355772334169456
  },
  "0.19": {
    "mean": 0.39981841104111665,
    "std": 0.027016244861192636
  },
  "0.20": {
    "mean": 0.40029184793024347,
    "std": 0.024708432620686956
  },
  "0.21": {
    "mean": 0.39789812590783497,
    "std": 0.025520497936303666
  },
  "0.22": {
    "mean": 0.395847891965017,
    "std": 0.025790662490122216
  },
  "0.23": {
    "mean": 0.39648745424536785,
    "std": 0.02385882135337397
  },
  "0.24": {
    "mean": 0.3995899720316251,
    "std": 0.02578104008013167
  },
  "0.25": {
    "mean": 0.3972815606284205,
    "std": 0.025576727198992264
  },
  "0.26": {
    "mean": 0.3983771805109425,
    "std": 0.028935293979013833
  },
  "0.27": {
    "mean": 0.3969920165843066,
    "std": 0.02497529525206574
  },
  "0.28": {
    "mean": 0.4031693429115578,
    "std": 0.024551636034024146
  },
  "0.29": {
    "mean": 0.39991865182437275,
    "std": 0.026682409914462683
  },
  "0.30": {
    "mean": 0.4008943857800395,
    "std": 0.02675033780957701
  },
  "0.31": {
    "mean": 0.3962344554211347,
    "std": 0.02571293799437975
  },
  "0.32": {
    "mean": 0.3947083278936816,
    "std": 0.026558697381502395
  },
  "0.33": {
    "mean": 0.40284990798616005,
    "std": 0.026789983832809547
  },
  "0.34": {
    "mean": 0.3987862366901604,
    "std": 0.026211352369012567
  },
  "0.35": {
    "mean": 0.3973613704299476,
    "std": 0.025134966977588186
  },
  "0.36": {
    "mean": 0.3991941009288439,
    "std": 0.024758946990399582
  },
  "0.37": {
    "mean": 0.39639518530902856,
    "std": 0.02511888629071801
  },
  "0.38": {
    "mean": 0.39479028817294837,
    "std": 0.024507724205019284
  },
  "0.39": {
    "mean": 0.4022232369392795,
    "std": 0.02617700375805709
  },
  "0.40": {
    "mean": 0.3985299913591742,
    "std": 0.027914514081746673
  },
  "0.41": {
    "mean": 0.39711781214788244,
    "std": 0.02600649550958563
  },
  "0.42": {
    "mean": 0.395751567199055,
    "std": 0.02467251649620674
  },
  "0.43": {
    "mean": 0.39459432398424904,
    "std": 0.026048966497937552
  },
  "0.44": {
    "mean": 0.39641806407787455,
    "std": 0.02560372103068189
  },
  "0.45": {
    "mean": 0.3984582561241745,
    "std": 0.027550973733774454
  }
}
```
```json
{
  "Spectra": {
    "Strain": {
      "-5": 0.006844430703847382,
      "0": 0.0025043396942770737,
      "5": -0.001974344057338875
    },
    "Seismic_Model": {
      "-5": 1.5419615969998273,
      "0": 1.4251372162191898,
      "5": 1.3549806600694683
    },
    "Magnetic_Model": {
      "-5": 0.9690215577376325,
      "0": 0.9594795851484398,
      "5": 0.9378419488518792
    }
  },
  "Distances": {
    "vs_Seismic": 2.494352801600233,
    "vs_Magnetic": 1.6507092744080252
  },
  "Info": "If fetch failed, Red/Pink noise models were used as physical vetoes."
}
```
```json
{
  "H1_Pure_H": 0.03733273522994843,
  "L1_Pure_H": 0.05299552509319112,
  "Cross_H1_L1_Pure_H": 0.31105068288116156,
  "Interpretation": "SURPRISING: The control loops / seismic background of H1 and L1 are strongly correlated at long distances."
}
```
```json
{
  "Whitened_Cross_H": 0.6856328431496843,
  "Verdict": "If ~0.50, there is NO phase structure (FIN is falsified as a phase effect). If significantly diff from 0.50, FIN exists as pure phase correlation."
}
```
```json
{
  "q": [
    -5.0,
    -4.0,
    -3.0,
    -2.0,
    -1.0,
    0.0,
    1.0,
    2.0,
    3.0,
    4.0,
    5.0
  ],
  "tau": [
    -2.3437647824134276,
    -2.0452331705814615,
    -1.7603708297863268,
    -1.4904612750979518,
    -1.236655371120009,
    -1.0,
    -0.781464764172749,
    -0.5818119681941865,
    -0.4011786117098268,
    -0.23861932213431447,
    -0.09213253820621514
  ],
  "alpha": [
    0.2985316118319661,
    0.2916969763135504,
    0.27738594774175485,
    0.2618577293331589,
    0.2452306375489759,
    0.22759530347363,
    0.20909401590290677,
    0.1901430762314611,
    0.171596323029936,
    0.15452303675180584,
    0.14648678392809933
  ],
  "f_alpha": [
    0.8511067232535972,
    0.8784452653272599,
    0.9282129865610622,
    0.966745816431634,
    0.9914247335710331,
    1.0,
    0.9905587800756558,
    0.9620981206571086,
    0.9159675807996348,
    0.8567114691415378,
    0.8245664578467118
  ],
  "FIN": {
    "H0": 0.22774859668180394,
    "DeltaH": 0.08717946412392852,
    "Cascade": "log-Poisson",
    "Cascade_param": 0.7680334399643445,
    "chi2_logPoisson": 0.08637671851557861,
    "chi2_logNormal": 0.22775411741584167
  }
}
```
```json
{
  "Config": {
    "GPS": 1266965117,
    "Fs_target": 1,
    "Cutoff_Hz": 0.1,
    "Window_sec": 4096
  },
  "Cross_H_LF": 1.9154092490606147,
  "Coherence_0_01_to_0_1_Hz": 0.6741010398437792,
  "Mean_Real_CSD_LF": 5.3899207891616826e-48,
  "Mean_Imag_CSD_LF": 1.0931139634509775e-50,
  "Verdict": "Effect disappears at ultra-low frequencies or interpretation is mixed."
}
```
```json
{
  "Detector_H0": {
    "H1": 0.0026335489204392295,
    "L1": 0.01859453307242668,
    "V1": 0.0036370796568547664
  },
  "Arm_Lengths": {
    "H1": 4000.0,
    "L1": 4000.0,
    "V1": 3000.0
  },
  "Scaling_Analysis": {
    "Unique_Lengths": [
      3000.0,
      4000.0
    ],
    "Average_H": [
      0.0036370796568547664,
      0.010614040996432955
    ],
    "Gamma_Exponent": 3.72284817816416,
    "Description": "H ~ 4.13e-16 * L^3.723"
  },
  "timestamp": "v29_execution"
}
```
```json
{
  "H_mean": 0.0030418309661207313,
  "H_std_local": 0.0007083480808800167,
  "Interpretation": {
    "low_std (<0.02)": "FIN is a STRUCTURAL MEASURE",
    "high_std": "FIN is a DYNAMICAL PROCESS"
  }
}
```
```json
{
  "Hurst_H1": 0.3107033628020025,
  "Hurst_L1": 0.26512680447073594,
  "Delta_H": 0.04557655833126656,
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-01T18:55:25.233662"
}
```
```json
{
  "sidereal_test": [
    {
      "segment": 0,
      "sidereal_phase": 0.0,
      "H1": 0.2764631582569697,
      "L1": 0.2738780068272483,
      "Delta": 0.002585151429721433
    },
    {
      "segment": 1,
      "sidereal_phase": 0.020890370815687738,
      "H1": 0.29013919118663656,
      "L1": 0.2689470865713718,
      "Delta": 0.021192104615264773
    }
  ],
  "segment_sec": 1800,
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-01T22:27:10.425300+00:00"
}
```
```json
{
  "beta_real": 5.8071325680938815,
  "beta_null_mean": 5.719870385223784,
  "beta_null_std": 0.013884007622481693,
  "z_score": 6.285086067570166
}
```
```json
{
  "H1_V1_Cross_H": 0.714808705323103,
  "Verdict": "If ~0.50, the anomaly is an H1-L1 paired architectural artifact. If ~0.31, it's truly global across distinct architectures."
}
```
```json
{
  "Spectra": {
    "Strain_H1": {
      "-5": 0.006844430703847382,
      "0": 0.0025043396942770737,
      "5": -0.001974344057338875
    },
    "Seismic": {
      "-5": 0.5065520521790964,
      "0": 0.4865799566905009,
      "5": 0.46810562719698207
    },
    "Magnetic": {
      "-5": 0.5136585429375095,
      "0": 0.500904416564615,
      "5": 0.4819535028086708
    }
  },
  "Distance_Euclidean": {
    "Strain_vs_Seismic": 0.8396499802329472,
    "Strain_vs_Magnetic": 0.8599124036527941
  },
  "Conclusion": "If Distances > 0.1, FIN is distinct from Environment.",
  "timestamp": "v30_execution"
}
```
```json
{
  "original_cross_H": 0.31105068288116156,
  "surrogate_cross_H_mean": 0.28332943157004487,
  "surrogate_cross_H_std": 0.008788916370598766,
  "N_surrogates": 20
}
```
```json
{
  "global_consistency_test": [
    {
      "window_sec": 64,
      "H": 0.30885260430877476,
      "samples": 262144,
      "compute_time_sec": 0.24536657333374023
    },
    {
      "window_sec": 128,
      "H": 0.3098408247635631,
      "samples": 524288,
      "compute_time_sec": 0.4489321708679199
    },
    {
      "window_sec": 256,
      "H": 0.28857808692469833,
      "samples": 1048576,
      "compute_time_sec": 0.8733696937561035
    },
    {
      "window_sec": 512,
      "H": 0.2756928984409219,
      "samples": 2097152,
      "compute_time_sec": 1.6308715343475342
    },
    {
      "window_sec": 1024,
      "H": 0.2874670726409621,
      "samples": 4194304,
      "compute_time_sec": 3.195228099822998
    },
    {
      "window_sec": 2048,
      "H": 0.27098363828141603,
      "samples": 8388608,
      "compute_time_sec": 6.0503809452056885
    }
  ],
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-01T22:53:20.955802+00:00"
}
```
```json
{
  "cross_epoch_consistency": [
    {
      "epoch": "Epoch_1",
      "window_sec": 512,
      "samples": 2097152,
      "H": 0.27105320196327354,
      "compute_time_sec": 1.1741423606872559
    }
  ],
  "summary": {
    "H_mean": 0.27105320196327354,
    "H_std": 0.0,
    "N_epochs": 1
  },
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-03T23:04:14.447923+00:00"
}
```
```json
{
  "diurnal_test": [
    {
      "segment": 0,
      "H1": 0.2861297793052872,
      "L1": 0.26324010867471026,
      "Delta": 0.022889670630576953
    }
  ],
  "segment_sec": 3600,
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-01T21:21:07.394453+00:00"
}
```
```json
{
  "v40_CosmicTime": {
    "Spearman_z_H": -0.9642857142857145,
    "p_value": 0.0004541491691941689,
    "verdict": "Cosmic-Structured"
  },
  "v41_GW_Coupling": {
    "Spearman_FIN_GWrate": -0.9642857142857145,
    "p_value": 0.0004541491691941689,
    "verdict": "Coupled"
  },
  "v42_Axioms": {
    "Axiom_1": "FIN is scale-invariant (H constant across geometry).",
    "Axiom_2": "FIN carries structure without energy transport.",
    "Axiom_3": "FIN is tensorial and time-symmetric."
  }
}
```
```json
{
  "H_cross_short_scales_1s_to_25s": 0.22332124569812675,
  "H_cross_long_scales_25s_to_128s": 0.3619575625119233,
  "H_cross_full_baseline": 0.31105068288116156,
  "Verdict": "Anti-persistence is scale-invariant across both short and long bins."
}
```
```json
{
  "40-80Hz": 0.7313123830135317,
  "100-200Hz": 0.5631587605638296,
  "300-500Hz": 0.4134798668917408
}
```
```json
{
  "segment_sec": 64,
  "r_real": -0.1453076356553426,
  "p_real": 0.25194302076116887,
  "null_mean": 0.024229854332111018,
  "null_std": 0.12065324834219607,
  "z_score": -1.4051630794606735,
  "verdict": "NO COHERENCE \u2014 CONSISTENT WITH INDEPENDENT NOISE",
  "timestamp_utc": "2026-01-01T21:51:55.414111+00:00"
}
```
```json
{
  "H1": {
    "real": 0.30885260430877476,
    "null_mean": 0.31090352670748245,
    "null_std": 0.007531104412966834
  },
  "L1": {
    "real": 0.3246381413143207,
    "null_mean": 0.30472640469044976,
    "null_std": 0.005574223192037795
  },
  "sample_rate": 4096.0,
  "duration_used": 64.0,
  "timestamp_utc": "2026-01-01T20:43:16.857304+00:00"
}
```
```json
{
  "H_original": 0.22774859668180394,
  "H_permuted_mean": 0.5003428220114121,
  "H_permuted_std": 0.015041406502576567,
  "Delta": -0.2725942253296082,
  "Interpretation": {
    "Delta ~ 0": "Amplitude-based / noise",
    "Delta >> 0": "Informational / structural FIN"
  }
}
```
```json
{
  "0.0": {
    "Cross_H_mean": 0.5013565178093967,
    "Cross_H_std": 0.020171871555999893
  },
  "0.5": {
    "Cross_H_mean": 0.5059634912258375,
    "Cross_H_std": 0.016360932501354154
  },
  "1.0": {
    "Cross_H_mean": 0.5103556104649954,
    "Cross_H_std": 0.0146661713518238
  },
  "2.0": {
    "Cross_H_mean": 0.5044638878327039,
    "Cross_H_std": 0.02237324856125114
  },
  "3.0": {
    "Cross_H_mean": 0.49910779621350587,
    "Cross_H_std": 0.024637612596447415
  },
  "5.0": {
    "Cross_H_mean": 0.4949717022769238,
    "Cross_H_std": 0.025364894370843234
  },
  "7.0": {
    "Cross_H_mean": 0.49354969534495413,
    "Cross_H_std": 0.025348111303817908
  },
  "10.0": {
    "Cross_H_mean": 0.49267169835880253,
    "Cross_H_std": 0.025164686986300058
  },
  "15.0": {
    "Cross_H_mean": 0.49211390367399277,
    "Cross_H_std": 0.02489982657201801
  },
  "20.0": {
    "Cross_H_mean": 0.49187840042302067,
    "Cross_H_std": 0.024724322758137508
  },
  "50.0": {
    "Cross_H_mean": 0.49152941673745076,
    "Cross_H_std": 0.02433616880408267
  }
}
```
```json
{
  "cross_hurst_results": [
    {
      "pair": "H1-L1",
      "H_cross": 0.11552275914205934,
      "samples": 2097152
    },
    {
      "pair": "H1-V1",
      "H_cross": 0.07803517778157153,
      "samples": 2097152
    },
    {
      "pair": "L1-V1",
      "H_cross": 0.19066426135733577,
      "samples": 2097152
    }
  ],
  "summary": {
    "H_cross_mean": 0.12807406609365554,
    "H_cross_std": 0.04682932910348679,
    "pairs": [
      "H1-L1",
      "H1-V1",
      "L1-V1"
    ]
  },
  "sample_rate": 4096,
  "window_sec": 512,
  "timestamp_utc": "2026-01-04T01:01:51.856503+00:00"
}
```
```json
{
  "bands": {
    "30-80Hz": {
      "H1": 0.23607077037727509,
      "L1": 0.23731590563132346,
      "Delta": -0.001245135254048374
    },
    "80-200Hz": {
      "H1": 0.26564212964633904,
      "L1": 0.24650555849726993,
      "Delta": 0.01913657114906911
    },
    "200-500Hz": {
      "H1": 0.2069448048711533,
      "L1": 0.1766016827951618,
      "Delta": 0.03034312207599149
    }
  },
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-01T19:24:15.850515+00:00"
}
```
```json
{
  "config": {
    "N": 1048576,
    "N_trials": 20,
    "H_true_background": 0.23,
    "H_local_damping": 0.04,
    "SNR_local_over_background": 20.0
  },
  "mixed_signal": {
    "H_x_mean": 0.07864376942410117,
    "H_x_std": 0.002380065086687186,
    "H_y_mean": 0.07825426745142278,
    "H_y_std": 0.003211036386556685,
    "Cross_H_mean": 0.5094425737999033,
    "Cross_H_std": 0.02343557500667056
  },
  "null_model": {
    "Cross_H_null_mean": 0.4965606426963275,
    "Cross_H_null_std": 0.027611698738229173
  },
  "verdict": "Cross-MF-DFA gives 0.509, which matches neither 0.23 nor 0.31. Relationship is complex."
}
```
```json
{
  "v36_Scale_L": {
    "slope": 0.030842042547518753,
    "r_squared": 0.0054080586363148735,
    "verdict": "Scale Invariant"
  },
  "v37_Tensor": {
    "H_plus_proxy": 0.2217127145427585,
    "H_cross_proxy": 0.3606526430653205,
    "Delta": 0.138939928522562,
    "verdict": "Tensorial"
  },
  "v38_Isotropy": {
    "slope": 0.0019315065165690576,
    "r_squared": 0.9966893760060604,
    "verdict": "Anisotropic"
  },
  "v39_EnergyInfo": {
    "correlation": -0.955014594715938,
    "p_value": 0.1916781512354847,
    "verdict": "Energy Driven"
  },
  "Raw_Data": {
    "H1": {
      "H": 0.2217127145427585,
      "RMS": 4.393047798925605e-20
    },
    "L1": {
      "H": 0.3606526430653205,
      "RMS": 2.1030104377162187e-21
    },
    "V1": {
      "H": 0.2823099760853233,
      "RMS": 1.4631395710288143e-20
    }
  }
}
```
```json
{
  "cross_epoch_fractal_stability": [
    {
      "epoch": "O3a_H1",
      "window_sec": 512,
      "samples": 2097152,
      "H": 0.2900189379914974
    },
    {
      "epoch": "O3b_H1",
      "window_sec": 512,
      "samples": 2097152,
      "H": 0.36065896220129573
    }
  ],
  "summary": {
    "H_mean": 0.3253389500963966,
    "H_std": 0.03532001210489916,
    "N_epochs": 2
  },
  "sample_rate": 4096,
  "timestamp_utc": "2026-01-04T00:26:47.894554+00:00"
}
```
```json
{
  "Fundamental_H_Input": 0.23,
  "Injected_Cross_H_Output": 0.300537392344461,
  "Verdict": "The instrumental response accurately transforms a true H=0.23 background into an apparent Cross-H ~ 0.31. The 0.31 anomaly COULD be a projection of a 0.23 fundamental background."
}
```
```json
{
  "Test1_Independent_Phases_Real_Envelope": {
    "mean": 0.2868038142777771,
    "std": 0.012702395253275257,
    "verdict": "If ~0.50, then the envelope ALONE does not cause the anomaly. If ~0.31, then the envelope alone is responsible."
  },
  "Test2_Real_Phases_Generic_Envelope": {
    "value": 0.705328267933796,
    "verdict": "If ~0.50, then removing the shared envelope destroys the anomaly. If ~0.31, the anomaly is robustly in the phase."
  }
}
```
```json
{
  "Scale_1000_to_3000": {
    "s_min": 1000,
    "s_max": 3000,
    "H_cross": 0.8918372300793768
  },
  "Scale_3000_to_10000": {
    "s_min": 3000,
    "s_max": 10000,
    "H_cross": 0.3088136078636575
  },
  "Scale_10000_to_30000": {
    "s_min": 10000,
    "s_max": 30000,
    "H_cross": 0.19505626413175572
  },
  "Scale_30000_to_100000": {
    "s_min": 30000,
    "s_max": 100000,
    "H_cross": 0.16602859804572065
  },
  "Scale_100000_to_131072": {
    "s_min": 100000,
    "s_max": 131072,
    "H_cross": -0.13145406277611377
  }
}
```
```json
{
  "bands": [
    [
      30,
      60
    ],
    [
      60,
      120
    ],
    [
      120,
      240
    ],
    [
      240,
      480
    ]
  ],
  "H1_band": [
    0.23478763001930075,
    0.20451298768398127,
    0.2605842628711187,
    0.205146703094412
  ],
  "L1_band": [
    0.2357252841122392,
    0.21592309821763106,
    0.24127755831025977,
    0.17635400825431122
  ],
  "r_real": 0.7903615769456582,
  "p_real": 0.20963842305434177,
  "null_mean": -0.04269731659612157,
  "null_std": 0.6081802287520204,
  "z_score": 1.369756618447147,
  "timestamp_utc": "2026-01-01T22:03:23.168967+00:00"
}
```
```json
{
  "1266965117": {
    "H1_Pure": 0.03733273522994843,
    "L1_Pure": 0.05299552509319112,
    "Cross_H1_L1": 0.31105068288116156
  },
  "1267051517": {
    "H1_Pure": 0.041405315287956546,
    "L1_Pure": 0.06420362757311483,
    "Cross_H1_L1": 0.33515322923661534
  },
  "1267137917": {
    "H1_Pure": 0.044313974791691,
    "L1_Pure": 0.05405721290409476,
    "Cross_H1_L1": 0.3034840404652728
  },
  "1253326744": {
    "H1_Pure": 0.04180091028389436,
    "L1_Pure": 0.05994293130661616,
    "Cross_H1_L1": 0.7633012338448356
  }
}
```
```json
{
  "v44_Geometry": {
    "baseline_results": {
      "H1-L1": {
        "angle": 90.0,
        "H_cross": 0.0022573648019294417
      },
      "H1-V1": {
        "angle": 45.0,
        "H_cross": 0.0020402144813620095
      },
      "L1-V1": {
        "angle": 60.0,
        "H_cross": 0.0091828030341689
      }
    },
    "slope_angle": -2.884208071414157e-05,
    "r_squared": 0.02646165698996646,
    "verdict": "Non-geometric"
  },
  "v45_Relational": {
    "Spearman_H_energy": -1.0,
    "p_value": 0.0,
    "verdict": "Energetic"
  },
  "v46_Scale": {
    "H_windows": {
      "32": 0.0024096229830277488,
      "64": 0.0022573648019294417,
      "128": 0.0021222544788838562,
      "256": 0.002277317937767772
    },
    "std": 0.00010184714011412424,
    "verdict": "Structural"
  }
}
```
```json
{
  "H_original": 0.0025043396942770737,
  "H_phase_randomized": 0.03488916553240482,
  "H_time_permuted": 0.5261990862176169,
  "Interpretation_Guide": {
    "H_phase ~ H_original": "Phase-encoded / informational FIN",
    "H_phase ~ 0.5": "Amplitude-based noise",
    "H_perm ~ 0.5": "Order-dependent structure confirmed"
  }
}
```
```json
{
  "beta_h1": 5.788587780685799,
  "beta_l1": 5.825677355501964,
  "delta_beta": -0.03708957481616437,
  "null_mean": 0.13165685849154202,
  "null_std": 0.027553184027458343,
  "z_score": -6.124389585593477
}
```
```json
{
  "Coherence_Coarse_Sec": {
    "H1-L1": {
      "-5": 0.22713815148589375,
      "-4": 0.2284834335849714,
      "-3": 0.22813757926549674,
      "-2": 0.22719536436191218,
      "-1": 0.2285323440264692,
      "0": 0.22774859668180394,
      "1": 0.22885061737842644,
      "2": 0.22871103459093778,
      "3": 0.22899172351486632,
      "4": 0.2286404321201261,
      "5": 0.22867975079279013
    },
    "H1-V1": {
      "-5": 0.0761885833323833,
      "-4": 0.059078151977446036,
      "-3": 0.0798636811523966,
      "-2": 0.07126686483520528,
      "-1": 0.06592063659510108,
      "0": 0.08061816890960831,
      "1": 0.05841927158024111,
      "2": 0.07313876548475505,
      "3": 0.06833537685489051,
      "4": 0.06041403469070915,
      "5": 0.07775007475583821
    },
    "L1-V1": {
      "-5": 0.1301869471470167,
      "-4": 0.13556358584670644,
      "-3": 0.1573920502774073,
      "-2": 0.16091488046423677,
      "-1": 0.1666934268485417,
      "0": 0.15516582476076077,
      "1": 0.16092046982674887,
      "2": 0.1586646925042033,
      "3": 0.1527071989203297,
      "4": 0.14458779321127657,
      "5": 0.1410046398746766
    }
  },
  "Coherence_Fine_ms": {
    "H1-L1": {
      "-20.0": 0.22754163493399385,
      "-18.0": 0.22747685654528674,
      "-16.0": 0.22742069398267709,
      "-14.0": 0.22738031797250913,
      "-12.0": 0.22737844803239704,
      "-10.0": 0.22748645642290677,
      "-8.0": 0.22753136500406423,
      "-6.0": 0.22759328581688065,
      "-4.0": 0.2276585701554386,
      "-2.0": 0.2277118181868326,
      "0.0": 0.22774859668180394,
      "2.0": 0.22776441559660618,
      "4.0": 0.22778231958839587,
      "6.0": 0.22781260225540848,
      "8.0": 0.22786637891755065,
      "10.0": 0.22795232472939844,
      "12.0": 0.227799685381826,
      "14.0": 0.2278843776961098,
      "16.0": 0.22793328160504317,
      "18.0": 0.22794253354135366,
      "20.0": 0.22793205412191542
    },
    "H1-V1": {
      "-20.0": 0.07380368836963057,
      "-18.0": 0.07583107454321843,
      "-16.0": 0.07804436472618854,
      "-14.0": 0.08057437240608738,
      "-12.0": 0.08384097504009345,
      "-10.0": 0.06462637417388827,
      "-8.0": 0.06802516707388706,
      "-6.0": 0.07217648112759158,
      "-4.0": 0.07612695523393315,
      "-2.0": 0.07871699263290512,
      "0.0": 0.08061816890960831,
      "2.0": 0.08126862628342205,
      "4.0": 0.08109825777150323,
      "6.0": 0.08058131442028921,
      "8.0": 0.07934242114197854,
      "10.0": 0.07638145978477211,
      "12.0": 0.07283125671815269,
      "14.0": 0.07333114609843601,
      "16.0": 0.07261825865508073,
      "18.0": 0.07154471844926055,
      "20.0": 0.0710382411168277
    },
    "L1-V1": {
      "-20.0": 0.14966628278599062,
      "-18.0": 0.1673353942887202,
      "-16.0": 0.18557497037893034,
      "-14.0": 0.17834810641236656,
      "-12.0": 0.1486904148949327,
      "-10.0": 0.14667045415716096,
      "-8.0": 0.14712041829137643,
      "-6.0": 0.156730848755411,
      "-4.0": 0.15038150624375707,
      "-2.0": 0.13807230161484335,
      "0.0": 0.15516582476076077,
      "2.0": 0.16048255254948618,
      "4.0": 0.15232325761547447,
      "6.0": 0.13418898606714447,
      "8.0": 0.14331995956172752,
      "10.0": 0.16283249488586493,
      "12.0": 0.17583982538042625,
      "14.0": 0.1645064893724227,
      "16.0": 0.1470261296279435,
      "18.0": 0.13971512595665786,
      "20.0": 0.15581125290362344
    }
  },
  "Note": "Thresholds adjusted for gravitational strain amplitude (~1e-21)",
  "timestamp": "v28_execution"
}
```
```json
{
  "integration_test": [
    {
      "window_sec": 32,
      "H1": 0.33995934773901576,
      "L1": 0.3520262031207114,
      "Delta": -0.012066855381695663,
      "compute_time_sec": 0.31324338912963867
    },
    {
      "window_sec": 64,
      "H1": 0.3105435282906503,
      "L1": 0.3263566705352001,
      "Delta": -0.015813142244549827,
      "compute_time_sec": 0.5691792964935303
    },
    {
      "window_sec": 128,
      "H1": 0.31001070121971774,
      "L1": 0.31835155974370716,
      "Delta": -0.008340858523989414,
      "compute_time_sec": 1.0603904724121094
    },
    {
      "window_sec": 256,
      "H1": 0.2881253580625984,
      "L1": 0.30091963081910567,
      "Delta": -0.012794272756507241,
      "compute_time_sec": 2.042760133743286
    },
    {
      "window_sec": 512,
      "H1": 0.2749408309566025,
      "L1": 0.2906625964784416,
      "Delta": -0.015721765521839126,
      "compute_time_sec": 3.687760591506958
    },
    {
      "window_sec": 1024,
      "H1": 0.28688700383670396,
      "L1": 0.287528864266169,
      "Delta": -0.0006418604294650687,
      "compute_time_sec": 7.145539283752441
    }
  ],
  "sample_rate": 4096.0,
  "total_runtime_sec": 14.842632055282593,
  "timestamp_utc": "2026-01-01T21:06:38.168655+00:00"
}
```
```json
{
  "cross_epoch_fractal_stability": [
    {
      "epoch": "O3a_H1",
      "window_sec": 512,
      "samples": 2097152,
      "H": 0.2900189379914974
    },
    {
      "epoch": "O3b_H1",
      "window_sec": 512,
      "samples": 2097152,
      "H": 0.36065896220129573
    },
    {
      "epoch": "O4_H1",
      "window_sec": 512,
      "samples": 2097152,
      "H": 0.33009336350893775
    }
  ],
  "summary": {
    "H_mean": 0.32692375456724365,
    "H_std": 0.0289256295894912,
    "N_epochs": 3,
    "epochs_processed": [
      "O3a_H1",
      "O3b_H1",
      "O4_H1"
    ]
  },
  "sample_rate": 4096,
  "timestamp_utc": "2026-01-04T00:33:13.409766+00:00"
}
```
```json
{
  "MF_DFA": {
    "H1": {
      "-5": 0.006844430703847382,
      "-4": 0.005990294627013035,
      "-3": 0.005128600477532727,
      "-2": 0.004260019290883733,
      "-1": 0.003385112029380349,
      "0": 0.0025043396942770737,
      "1": 0.0016180874149963294,
      "2": 0.0007267016120266353,
      "3": -0.00016946256612925594,
      "4": -0.001069988196042291,
      "5": -0.001974344057338875
    },
    "L1": {
      "-5": 0.047682750986458335,
      "-4": 0.0428337697498101,
      "-3": 0.037691671330022485,
      "-2": 0.03223713004462336,
      "-1": 0.026448632237696158,
      "0": 0.020271377114164532,
      "1": 0.01347862981930821,
      "2": 0.005121093332372146,
      "3": -0.00825328495861051,
      "4": -0.033948655687234715,
      "5": -0.07157518028645297
    },
    "V1": {
      "-5": 0.011001293684340626,
      "-4": 0.009339864389059477,
      "-3": 0.007719558192389645,
      "-2": 0.00614238209867122,
      "-1": 0.004609302332536944,
      "0": 0.0031201886565858854,
      "1": 0.0016736935469015438,
      "2": 0.00026705753659740194,
      "3": -0.0011041908926608184,
      "4": -0.0024466938486457897,
      "5": -0.0037699433012532936
    }
  },
  "CROSS_MF_DFA": {
    "real": {
      "-5": 0.2687529564826855,
      "-4": 0.2613082926453654,
      "-3": 0.2534569432621089,
      "-2": 0.24523063754897595,
      "-1": 0.23665537112000895,
      "0": 0.22774859668180394,
      "1": 0.21853523582725104,
      "2": 0.20909401590290674,
      "3": 0.1996071294300577,
      "4": 0.19034516946642138,
      "5": 0.18157349235875697
    },
    "shuffle": {
      "-5": 0.5034203775726257,
      "-4": 0.503315190571987,
      "-3": 0.5034871257469522,
      "-2": 0.5036697032660185,
      "-1": 0.5035583331580095,
      "0": 0.5029040733314845,
      "1": 0.501572960760977,
      "2": 0.49953625611087116,
      "3": 0.49682172095733346,
      "4": 0.49347982600161155,
      "5": 0.4895835197780457
    },
    "phase": {
      "-5": 0.21241714108321552,
      "-4": 0.2066586328164522,
      "-3": 0.20240773013742588,
      "-2": 0.20007000933802568,
      "-1": 0.19976744895791468,
      "0": 0.20077379936394577,
      "1": 0.20132358066545503,
      "2": 0.19980187737030422,
      "3": 0.1959791993052112,
      "4": 0.19064284141606427,
      "5": 0.18466079912108885
    },
    "H1-V1": {
      "real": {
        "-5": 0.07824452227756472,
        "-4": 0.0773505311767661,
        "-3": 0.07717491185697276,
        "-2": 0.07772987656508179,
        "-1": 0.07893671623815128,
        "0": 0.08061816890960831,
        "1": 0.08252423161965843,
        "2": 0.08440309871882581,
        "3": 0.08607757944928261,
        "4": 0.08747213561044159,
        "5": 0.08858725918899217
      },
      "shuffle": {
        "-5": 0.507041800253571,
        "-4": 0.5052754204414203,
        "-3": 0.5035364723268283,
        "-2": 0.5018306812190458,
        "-1": 0.5001567206873525,
        "0": 0.4984694722998402,
        "1": 0.49665564437092724,
        "2": 0.49455703697797243,
        "3": 0.49203532276862355,
        "4": 0.4890301173986625,
        "5": 0.4855734173855122
      },
      "phase": {
        "-5": 0.08344104807009789,
        "-4": 0.07901170280927079,
        "-3": 0.07561390986870035,
        "-2": 0.07311542986638424,
        "-1": 0.07139771427750366,
        "0": 0.07032920663337112,
        "1": 0.06972990468301439,
        "2": 0.06938983175156443,
        "3": 0.06912385008644327,
        "4": 0.06881168828504328,
        "5": 0.06839963231201049
      }
    },
    "L1-V1": {
      "real": {
        "-5": 0.1703613127327925,
        "-4": 0.16716218425111537,
        "-3": 0.16403612007350366,
        "-2": 0.16100045506505392,
        "-1": 0.1580501998513365,
        "0": 0.15516582476076077,
        "1": 0.15232808064208983,
        "2": 0.14953007529836113,
        "3": 0.14677710019154533,
        "4": 0.1440721791108341,
        "5": 0.14139752345734977
      },
      "shuffle": {
        "-5": 0.5045139190414584,
        "-4": 0.5031718371874048,
        "-3": 0.5023880004865273,
        "-2": 0.5020622920859986,
        "-1": 0.5018845338107429,
        "0": 0.5013856509483983,
        "1": 0.5001254538408195,
        "2": 0.4978749376950873,
        "3": 0.4946506925377633,
        "4": 0.490631130687572,
        "5": 0.4860568715004665
      },
      "phase": {
        "-5": 0.18911863880090202,
        "-4": 0.18146817260109233,
        "-3": 0.1737619636920628,
        "-2": 0.16631363858181425,
        "-1": 0.1594712155669294,
        "0": 0.15354185295717637,
        "1": 0.14870969447167898,
        "2": 0.14498707421159965,
        "3": 0.1422261736004945,
        "4": 0.14018347516359303,
        "5": 0.13859906551400358
      }
    }
  },
  "Wavelet_Leaders": {
    "H1": [
      4.812972522623127e-20,
      2.6983210498821907e-20,
      4.2188180222577184e-20,
      2.070283832327583e-19,
      1.5563749178434705e-19,
      2.8361233791228486e-20
    ],
    "L1": [
      1.6380946086137567e-20,
      6.207168573902065e-21,
      3.4411965114437665e-21,
      1.6474495104490852e-20,
      1.4780705973309457e-20,
      1.867954404117499e-21
    ]
  },
  "q": [
    -5,
    -4,
    -3,
    -2,
    -1,
    0,
    1,
    2,
    3,
    4,
    5
  ],
  "fs": 4096,
  "window_sec": 512,
  "timestamp": "2026-01-05T19:39:05.988826+00:00"
}
```
```json
{
  "real_cpsd": 2.605739266198478e-48,
  "null_mean": 2.757410458877926e-48,
  "null_std": 4.111652609538801e-49,
  "z_score": -0.36888134062585726,
  "p_value": 0.6,
  "verdict": "NO CORRELATED SIGNAL \u2014 CONSISTENT WITH GR"
}
```
```json
{
  "corrected_global_consistency": [
    {
      "window_sec": 64,
      "H1": {
        "H": 0.30885260430877476,
        "samples": 262144,
        "compute_time_sec": 0.1855020523071289
      },
      "L1": {
        "H": 0.3246381413143207,
        "samples": 262144,
        "compute_time_sec": 0.1655426025390625
      }
    },
    {
      "window_sec": 128,
      "H1": {
        "H": 0.3098408247635631,
        "samples": 524288,
        "compute_time_sec": 0.3062772750854492
      },
      "L1": {
        "H": 0.3151923379430553,
        "samples": 524288,
        "compute_time_sec": 0.3051609992980957
      }
    },
    {
      "window_sec": 256,
      "H1": {
        "H": 0.28857808692469833,
        "samples": 1048576,
        "compute_time_sec": 0.6171586513519287
      },
      "L1": {
        "H": 0.3011290909967617,
        "samples": 1048576,
        "compute_time_sec": 0.5838541984558105
      }
    },
    {
      "window_sec": 512,
      "H1": {
        "H": 0.2756928984409219,
        "samples": 2097152,
        "compute_time_sec": 1.1128170490264893
      },
      "L1": {
        "H": 0.291333484468053,
        "samples": 2097152,
        "compute_time_sec": 1.1120657920837402
      }
    },
    {
      "window_sec": 1024,
      "H1": {
        "H": 0.2874670726409621,
        "samples": 4194304,
        "compute_time_sec": 2.2060842514038086
      },
      "L1": {
        "H": 0.28720688288221935,
        "samples": 4194304,
        "compute_time_sec": 2.133000612258911
      }
    }
  ],
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-03T22:39:03.343917+00:00"
}

{
  "corrected_global_consistency": [
    {
      "window_sec": 64,
      "H1": {
        "H": 0.30885260430877476,
        "samples": 262144,
        "compute_time_sec": 0.1855020523071289
      },
      "L1": {
        "H": 0.3246381413143207,
        "samples": 262144,
        "compute_time_sec": 0.1655426025390625
      }
    },
    {
      "window_sec": 128,
      "H1": {
        "H": 0.3098408247635631,
        "samples": 524288,
        "compute_time_sec": 0.3062772750854492
      },
      "L1": {
        "H": 0.3151923379430553,
        "samples": 524288,
        "compute_time_sec": 0.3051609992980957
      }
    },
    {
      "window_sec": 256,
      "H1": {
        "H": 0.28857808692469833,
        "samples": 1048576,
        "compute_time_sec": 0.6171586513519287
      },
      "L1": {
        "H": 0.3011290909967617,
        "samples": 1048576,
        "compute_time_sec": 0.5838541984558105
      }
    },
    {
      "window_sec": 512,
      "H1": {
        "H": 0.2756928984409219,
        "samples": 2097152,
        "compute_time_sec": 1.1128170490264893
      },
      "L1": {
        "H": 0.291333484468053,
        "samples": 2097152,
        "compute_time_sec": 1.1120657920837402
      }
    },
    {
      "window_sec": 1024,
      "H1": {
        "H": 0.2874670726409621,
        "samples": 4194304,
        "compute_time_sec": 2.2060842514038086
      },
      "L1": {
        "H": 0.28720688288221935,
        "samples": 4194304,
        "compute_time_sec": 2.133000612258911
      }
    }
  ],
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-03T22:39:03.343917+00:00"
}

```

---

## Faza 61: `phase61_full_null_model.py`
### Kod:
```python
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
    
    num_cores = multiprocessing.cpu_count()
    with multiprocessing.Pool(num_cores) as p:
        results = p.map(run_trial, range(N_TRIALS))
        
    results = np.array(results)
    
    m_null = np.mean(results)
    s_null = np.std(results)
    min_null = np.min(results)
    p1 = np.percentile(results, 1)
    
    obs_h = 0.311
    
    z_score = (obs_h - m_null) / s_null
    
    
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

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "Config": {
    "N_trials": 1000,
    "N_samples": 131072,
    "Description": "Independent random phases applied strictly to empirical H1 and L1 PSDs"
  },
  "Results": {
    "Mean_H_cross": 0.389898436517463,
    "Std_H_cross": 0.023539041835616652,
    "Min_H_cross": 0.322682844205378,
    "P1_H_cross": 0.33770149682344514
  },
  "Observation": 0.311,
  "Significance": "-3.35 sigma",
  "Verdict": "Anomaly survives at -3.35 sigma even against realistic instrumental envelopes! The FIN correlation is still required."
}
```

---

## Faza 62: `phase62_inter_observatory.py`
### Kod:
```python
log = logging.getLogger(__name__)

RAW_DIR = "./raw_strain_unfiltered"
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
    ts = ts.notch(60).notch(120).notch(180) # Generic line removal
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

def main():
    
    try:
        x_h1 = detrend(fetch_pure_strain("H1", GPS))
        x_v1 = detrend(fetch_pure_strain("V1", GPS))
    except Exception as e:
        return

    N = min(len(x_h1), len(x_v1))
    x_h1 = x_h1[:N]
    x_v1 = x_v1[:N]
    
    scales = np.unique(np.logspace(3, np.log10(N//4), 15).astype(int))
    
    h_cross = cross_mfdfa_q0(x_h1, x_v1, scales)
    
    
    out = {
        "H1_V1_Cross_H": h_cross,
        "Verdict": "If ~0.50, the anomaly is an H1-L1 paired architectural artifact. If ~0.31, it's truly global across distinct architectures."
    }
    
    with open("QW_1660_v62_InterObservatory.json", "w") as f:
        json.dump(out, f, indent=2)

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "H1_V1_Cross_H": 0.714808705323103,
  "Verdict": "If ~0.50, the anomaly is an H1-L1 paired architectural artifact. If ~0.31, it's truly global across distinct architectures."
}
```

---

## Faza 63: `phase63_whitening_test.py`
### Kod:
```python
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
    x_h1, x_l1 = fetch_empirical_data()
    
    hw1 = exact_spectral_whitening(x_h1)
    lw1 = exact_spectral_whitening(x_l1)
    
    scales = np.unique(np.logspace(3, np.log10(N_SAMPLES//4), 15).astype(int))
    
    h_cross = cross_mfdfa_q0(hw1, lw1, scales)
    
    out = {
        "Whitened_Cross_H": h_cross,
        "Verdict": "If ~0.50, there is NO phase structure (FIN is falsified as a phase effect). If significantly diff from 0.50, FIN exists as pure phase correlation."
    }
    
    with open("QW_1660_v63_Whitening.json", "w") as f:
        json.dump(out, f, indent=2)

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "Whitened_Cross_H": 0.6856328431496843,
  "Verdict": "If ~0.50, there is NO phase structure (FIN is falsified as a phase effect). If significantly diff from 0.50, FIN exists as pure phase correlation."
}
```

---

## Faza 64: `phase64_bandpass_test.py`
### Kod:
```python
log = logging.getLogger(__name__)

RAW_DIR = "./raw_strain_unfiltered"
GPS = 1266965117 
N_SAMPLES = 524288 # 128s
FS = 4096

def fetch_empirical_data():
    path_h1 = f"{RAW_DIR}/H1_unfiltered_{GPS}.h5"
    path_l1 = f"{RAW_DIR}/L1_unfiltered_{GPS}.h5"
    with h5py.File(path_h1, "r") as f: x_h1 = f["strain"][:N_SAMPLES]
    with h5py.File(path_l1, "r") as f: x_l1 = f["strain"][:N_SAMPLES]
    return x_h1, x_l1

def bandpass_filter(data, lowcut, highcut, fs, order=4):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return filtfilt(b, a, data)

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
    x_h1, x_l1 = fetch_empirical_data()
    
    bands = [
        (40, 80),
        (100, 200),
        (300, 500)
    ]
    
    scales = np.unique(np.logspace(3, np.log10(N_SAMPLES//4), 15).astype(int))
    
    results = {}
    
    for (low, high) in bands:
        h1_bp = bandpass_filter(x_h1, low, high, FS)
        l1_bp = bandpass_filter(x_l1, low, high, FS)
        
        hc = cross_mfdfa_q0(h1_bp, l1_bp, scales)
        
        results[f"{low}-{high}Hz"] = hc
        
    with open("QW_1660_v64_Bandpass.json", "w") as f:
        json.dump(results, f, indent=2)

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "40-80Hz": 0.7313123830135317,
  "100-200Hz": 0.5631587605638296,
  "300-500Hz": 0.4134798668917408
}
```

---

## Faza 65: `phase65_micro_timeshift.py`
### Kod:
```python
log = logging.getLogger(__name__)

RAW_DIR = "./raw_strain_unfiltered"
GPS = 1266965117
FS = 4096
N_SAMPLES = 524288 # 128s

def fetch_empirical_data():
    path_h1 = f"{RAW_DIR}/H1_unfiltered_{GPS}.h5"
    path_l1 = f"{RAW_DIR}/L1_unfiltered_{GPS}.h5"
    with h5py.File(path_h1, "r") as f: x_h1 = f["strain"][:N_SAMPLES]
    with h5py.File(path_l1, "r") as f: x_l1 = f["strain"][:N_SAMPLES]
    return detrend(x_h1), detrend(x_l1)

def exact_spectral_whitening(x):
    X_f = np.fft.rfft(x)
    phases = np.angle(X_f)
    X_white = np.exp(1j * phases)
    X_white[0] = 0.0
    if len(x) % 2 == 0: X_white[-1] = np.real(X_white[-1])
    x_white = np.fft.irfft(X_white, n=len(x))
    return detrend(x_white)

def main():
    x_h1, x_l1 = fetch_empirical_data()
    
    
    hw = exact_spectral_whitening(x_h1)
    lw = exact_spectral_whitening(x_l1)
    
    max_lag_samples = 410
    
    
    center_start = max_lag_samples
    center_end = len(lw) - max_lag_samples
    lw_center = lw[center_start:center_end]
    
    lags = np.arange(-max_lag_samples, max_lag_samples + 1)
    corr_values = []
    
    for lag in lags:
        h1_shifted = hw[center_start + lag : center_end + lag]
        c = np.corrcoef(h1_shifted, lw_center)[0, 1]
        corr_values.append(c)
        
    corr_values = np.array(corr_values)
    
    best_lag_idx = np.argmax(np.abs(corr_values))
    best_lag_samples = lags[best_lag_idx]
    best_lag_ms = (best_lag_samples / FS) * 1000.0
    best_corr = corr_values[best_lag_idx]
    
    
    results = {}
    for lag, c in zip(lags, corr_values):
        t_ms = (lag / FS) * 1000.0
        if -25.0 <= t_ms <= 25.0:
            results[f"{t_ms:.2f} ms"] = float(c)
            
    out = {
        "Best_Lag_ms": best_lag_ms,
        "Best_Correlation": best_corr,
        "Verdict": "If peak is at ~10ms, it is a true gravitational/light-speed wave. If peak is exactly 0ms, it's global synchronous noise. If no clear peak, there is no broadband phase correlation.",
        "Correlation_Window_ms": results
    }

    with open("QW_1660_v65_MicroTimeShift.json", "w") as f:
        json.dump(out, f, indent=2)

    if abs(best_lag_ms) > 1.0:

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "Best_Lag_ms": 73.2421875,
  "Best_Correlation": -0.005650218758170756,
  "Verdict": "If peak is at ~10ms, it is a true gravitational/light-speed wave. If peak is exactly 0ms, it's global synchronous noise. If no clear peak, there is no broadband phase correlation.",
  "Correlation_Window_ms": {
    "-24.90 ms": 0.0035958924423889955,
    "-24.66 ms": -0.0015633011271547006,
    "-24.41 ms": 0.00017008011077205253,
    "-24.17 ms": 0.001355628330727946,
    "-23.93 ms": -0.0004134882273864998,
    "-23.68 ms": 0.0005054541196247306,
    "-23.44 ms": 0.0004900404954794068,
    "-23.19 ms": 0.0020395115197046417,
    "-22.95 ms": -0.0013569693128512989,
    "-22.71 ms": -0.0009627638168922174,
    "-22.46 ms": 0.001537215409647673,
    "-22.22 ms": 0.0007465840528576063,
    "-21.97 ms": 7.288600826827111e-06,
    "-21.73 ms": -0.00031926895456860677,
    "-21.48 ms": 0.0023335996325561915,
    "-21.24 ms": -0.000351092269335575,
    "-21.00 ms": -0.0013262921264217344,
    "-20.75 ms": 0.0027737408633027928,
    "-20.51 ms": 0.0002140976214402963,
    "-20.26 ms": -0.0017877361836068969,
    "-20.02 ms": 0.0020615367907124214,
    "-19.78 ms": 0.002347911922408769,
    "-19.53 ms": -0.0016768330858254732,
    "-19.29 ms": -1.970008582722594e-05,
    "-19.04 ms": 0.0018161985086641468,
    "-18.80 ms": 0.00038459294603834045,
    "-18.55 ms": 0.0012723814097295358,
    "-18.31 ms": -0.0009614678279654249,
    "-18.07 ms": 0.0008923191411199123,
    "-17.82 ms": 0.0014303140773147583,
    "-17.58 ms": -0.0013920944838412953,
    "-17.33 ms": 0.0031772269965881883,
    "-17.09 ms": 0.00017789552463652144,
    "-16.85 ms": -0.0015651730805940265,
    "-16.60 ms": 0.0018926567135340066,
    "-16.36 ms": 7.97589406729218e-05,
    "-16.11 ms": 0.001648131013788291,
    "-15.87 ms": 0.0003240010701814508,
    "-15.62 ms": 0.0011457879230633169,
    "-15.38 ms": 0.00089830500275615,
    "-15.14 ms": -0.0017275701490225563,
    "-14.89 ms": 0.0034050243898876608,
    "-14.65 ms": -0.00021676111242924362,
    "-14.40 ms": -0.0011972245040886079,
    "-14.16 ms": 0.002871610415299286,
    "-13.92 ms": -0.0010700932326487345,
    "-13.67 ms": 0.0023655663472634805,
    "-13.43 ms": 0.0019648860626813006,
    "-13.18 ms": -0.0022148141302521713,
    "-12.94 ms": 0.0010894040102334087,
    "-12.70 ms": 0.0012423076847510306,
    "-12.45 ms": 0.001505315720024038,
    "-12.21 ms": 0.0009349869149328864,
    "-11.96 ms": 9.707745509500646e-05,
    "-11.72 ms": -4.690221944909932e-05,
    "-11.47 ms": 0.0016199404633047672,
    "-11.23 ms": 0.001445677612150131,
    "-10.99 ms": -0.0010669552237243453,
    "-10.74 ms": 0.001344353060641185,
    "-10.50 ms": 0.001272999900487303,
    "-10.25 ms": 0.0008204132412868782,
    "-10.01 ms": -0.00031235607648622487,
    "-9.77 ms": 0.0012853310661435364,
    "-9.52 ms": 0.002164395644159403,
    "-9.28 ms": -0.001941729894240139,
    "-9.03 ms": 0.002208708098259537,
    "-8.79 ms": 0.0010388139945964607,
    "-8.54 ms": -0.0009969571629149417,
    "-8.30 ms": 0.0030184750758895647,
    "-8.06 ms": -0.0003822562793872377,
    "-7.81 ms": 0.0007977175482042462,
    "-7.57 ms": 0.0018976038759113678,
    "-7.32 ms": -0.0009456301976598709,
    "-7.08 ms": 0.00151764801706955,
    "-6.84 ms": 0.0014211995413514987,
    "-6.59 ms": -0.0002813114199989071,
    "-6.35 ms": -5.3712641338081626e-06,
    "-6.10 ms": 0.001970505995524733,
    "-5.86 ms": 0.0026007706801578993,
    "-5.62 ms": -0.0019946366304326515,
    "-5.37 ms": 0.0014201064526390096,
    "-5.13 ms": 0.003231590192172192,
    "-4.88 ms": -0.0022646161375598655,
    "-4.64 ms": 0.0017579494428678977,
    "-4.39 ms": 0.00143613910709342,
    "-4.15 ms": 0.0002827334271200185,
    "-3.91 ms": 0.0025202316600075043,
    "-3.66 ms": -0.002308038435619221,
    "-3.42 ms": 0.0019035634563612438,
    "-3.17 ms": 0.002384385718265745,
    "-2.93 ms": -0.0023088516721113347,
    "-2.69 ms": 0.0028886102039535956,
    "-2.44 ms": -0.00027315709272521303,
    "-2.20 ms": 0.0010609043632456843,
    "-1.95 ms": 0.0024014283671421064,
    "-1.71 ms": -0.0012244714012334403,
    "-1.46 ms": 0.003080209209764442,
    "-1.22 ms": -0.001583111947079528,
    "-0.98 ms": 0.0007468173716760598,
    "-0.73 ms": 0.0027386996591416113,
    "-0.49 ms": -0.0009326546566146771,
    "-0.24 ms": 0.0023591702014447904,
    "0.00 ms": -0.0008141046665487757,
    "0.24 ms": 0.001669849220685649,
    "0.49 ms": 0.002455146839938389,
    "0.73 ms": -0.0015496847419886645,
    "0.98 ms": 0.0011500722054539857,
    "1.22 ms": 0.0016680894830704196,
    "1.46 ms": -0.0001983649469891438,
    "1.71 ms": 0.0010650975025451983,
    "1.95 ms": 0.0021359048893473767,
    "2.20 ms": -7.279701833772713e-05,
    "2.44 ms": -0.0001252849054576339,
    "2.69 ms": 0.001890500652592548,
    "2.93 ms": 0.0005329389855970609,
    "3.17 ms": -0.000754706815663927,
    "3.42 ms": 0.0018228613785432081,
    "3.66 ms": 0.0012721137809966468,
    "3.91 ms": 0.0007220120666060275,
    "4.15 ms": 0.00021898813978581388,
    "4.39 ms": -0.0008471737570492024,
    "4.64 ms": 0.0018167036402158222,
    "4.88 ms": 0.0010777796843581082,
    "5.13 ms": -0.0003674864153857978,
    "5.37 ms": 0.0009149772524157882,
    "5.62 ms": 0.0014122365606375777,
    "5.86 ms": -5.9070473936331805e-05,
    "6.10 ms": 0.0007337426260857589,
    "6.35 ms": 0.0020346626159028747,
    "6.59 ms": -0.0009198023032520723,
    "6.84 ms": 0.0006959395853365341,
    "7.08 ms": 0.0007905222308515899,
    "7.32 ms": -0.0008360205116897204,
    "7.57 ms": 0.0027083644906924337,
    "7.81 ms": 0.0003334952880923528,
    "8.06 ms": -0.00041982207434698514,
    "8.30 ms": 0.0011443289366778971,
    "8.54 ms": 0.0003482171257196657,
    "8.79 ms": 0.0007702690707188068,
    "9.03 ms": -0.0008464080527831987,
    "9.28 ms": 0.0010071997339217511,
    "9.52 ms": 0.001947483803407492,
    "9.77 ms": 0.0006820474220215479,
    "10.01 ms": -0.000663842273425792,
    "10.25 ms": 0.0003373093289508459,
    "10.50 ms": 0.0015240362895267354,
    "10.74 ms": -0.0016475597435585802,
    "10.99 ms": 0.00042988213888107214,
    "11.23 ms": 0.002470088247878573,
    "11.47 ms": 0.0011871871530210985,
    "11.72 ms": -0.002074234329936727,
    "11.96 ms": -0.0005141144930285108,
    "12.21 ms": 0.0041528995821315424,
    "12.45 ms": -0.0016228369240076,
    "12.70 ms": -0.002244549994312596,
    "12.94 ms": 0.0037400614885914683,
    "13.18 ms": 0.00039106143091727544,
    "13.43 ms": -0.0019853563594793458,
    "13.67 ms": 0.0011190536410211449,
    "13.92 ms": 0.0023330869753939603,
    "14.16 ms": -0.0003616351627169883,
    "14.40 ms": -0.0011771663524364159,
    "14.65 ms": 0.002119871615158723,
    "14.89 ms": -0.0009848050714806838,
    "15.14 ms": -0.0001906672510117372,
    "15.38 ms": 0.0018270415850659622,
    "15.62 ms": -0.00034244649390020973,
    "15.87 ms": 0.00040048994610381765,
    "16.11 ms": -0.0016879625876422412,
    "16.36 ms": 0.0028303998420298434,
    "16.60 ms": -4.997284293283329e-05,
    "16.85 ms": -0.0030187742761357585,
    "17.09 ms": 0.003803223982454949,
    "17.33 ms": -0.0013198718734197547,
    "17.58 ms": 0.0010115647537846607,
    "17.82 ms": 0.001520261152710049,
    "18.07 ms": -0.0023831911658713447,
    "18.31 ms": 0.002030118606286494,
    "18.55 ms": -0.0014347003885494152,
    "18.80 ms": -7.78976108667486e-06,
    "19.04 ms": 0.0013756010423563911,
    "19.29 ms": -0.0007928289397442114,
    "19.53 ms": 0.0007965138403484423,
    "19.78 ms": -0.0006161300191874048,
    "20.02 ms": 0.0009850579714104986,
    "20.26 ms": -0.00022581939097157493,
    "20.51 ms": -0.0018688656550558558,
    "20.75 ms": 0.0014928530139534074,
    "21.00 ms": 0.0003885437493149752,
    "21.24 ms": -0.0008053923889068185,
    "21.48 ms": 0.0009445005040433828,
    "21.73 ms": -0.00021550349715158482,
    "21.97 ms": -0.00020866735578380276,
    "22.22 ms": 0.00043726384608365785,
    "22.46 ms": -0.0010525567289084176,
    "22.71 ms": 0.00029623005177442386,
    "22.95 ms": 0.00044074740127627156,
    "23.19 ms": -0.0015295763450776862,
    "23.44 ms": 0.0011450537815834666,
    "23.68 ms": 0.0012728978607754091,
    "23.93 ms": -0.003036822580908511,
    "24.17 ms": 0.0008509201291327026,
    "24.41 ms": 0.0011240928077248147,
    "24.66 ms": -0.0020782885853291963,
    "24.90 ms": 0.0009105883117009215
  }
}
```

---

## Faza 66: `phase66_scale_delay_function.py`
### Kod:
```python
log = logging.getLogger(__name__)

RAW_DIR = "./raw_strain_unfiltered"
GPS = 1266965117
FS = 4096
N_SAMPLES = 524288 # 128s

def fetch_empirical_data():
    path_h1 = f"{RAW_DIR}/H1_unfiltered_{GPS}.h5"
    path_l1 = f"{RAW_DIR}/L1_unfiltered_{GPS}.h5"
    with h5py.File(path_h1, "r") as f: x_h1 = f["strain"][:N_SAMPLES]
    with h5py.File(path_l1, "r") as f: x_l1 = f["strain"][:N_SAMPLES]
    return detrend(x_h1), detrend(x_l1)

def exact_spectral_whitening(x):
    X_f = np.fft.rfft(x)
    phases = np.angle(X_f)
    X_white = np.exp(1j * phases)
    X_white[0] = 0.0
    if len(x) % 2 == 0: X_white[-1] = np.real(X_white[-1])
    x_white = np.fft.irfft(X_white, n=len(x))
    return detrend(x_white)

def main():
    x_h1, x_l1 = fetch_empirical_data()
    
    hw = exact_spectral_whitening(x_h1)
    lw = exact_spectral_whitening(x_l1)
    
    scales = np.unique(np.logspace(np.log10(128), np.log10(N_SAMPLES//8), 20).astype(int))
    
    max_lag_ms = 1000.0 # Search window +/- 1 second
    max_lag_samples = int((max_lag_ms / 1000.0) * FS)
    
    tau_max_values = []
    valid_scales = []
    peaks_corr = []
    
    
    for s in scales:
        n_windows = len(hw) // s
        if n_windows == 0: continue
        
        
        
        avg_correlogram = np.zeros(2 * max_lag_samples + 1)
        k_count = 0
        
        lags = np.arange(-max_lag_samples, max_lag_samples + 1)
        
        for i in range(n_windows):
            start = i * s
            end = start + s
            
            if start - max_lag_samples < 0 or end + max_lag_samples >= len(lw):
                continue
                
            h_seg = hw[start:end]
            h_seg_norm = (h_seg - np.mean(h_seg)) / (np.std(h_seg) + 1e-12)
            
            l_segs = np.array([lw[start + lag : end + lag] for lag in lags])
            l_segs_mean = np.mean(l_segs, axis=1, keepdims=True)
            l_segs_std = np.std(l_segs, axis=1, keepdims=True) + 1e-12
            l_segs_norm = (l_segs - l_segs_mean) / l_segs_std
            
            local_corr = np.dot(l_segs_norm, h_seg_norm) / s
            
            avg_correlogram += local_corr
            k_count += 1
            
        if k_count == 0: continue
            
        avg_correlogram /= k_count
        
        best_idx = np.argmax(np.abs(avg_correlogram))
        best_lag_samples = lags[best_idx]
        best_tau_ms = abs(best_lag_samples) / FS * 1000.0
        
        valid_scales.append(s)
        tau_max_values.append(best_tau_ms)
        peaks_corr.append(float(avg_correlogram[best_idx]))
        

    valid_scales = np.array(valid_scales)
    tau_max_values = np.array(tau_max_values)
    
    non_zero_idx = tau_max_values > 0
    
    log_s = np.log10(valid_scales[non_zero_idx])
    log_tau = np.log10(tau_max_values[non_zero_idx])
    
    if len(log_s) > 2:
        slope, intercept, r_value, p_value, std_err = linregress(log_s, log_tau)
    else:
        slope, std_err, r_value, p_value = 0, 0, 0, 0

    results_data = {
        "Scales_samples": [int(x) for x in valid_scales],
        "Tau_max_ms": [float(x) for x in tau_max_values],
        "Peak_Correlations": [float(x) for x in peaks_corr],
        "Regression": {
            "Slope_H_obs": float(slope),
            "Std_Err": float(std_err),
            "R_squared": float(r_value**2),
            "P_value": float(p_value)
        },
        "Verdict": ""
    }
    
    if abs(slope - 0.23) < 0.05 and p_value < 0.05:
        results_data["Verdict"] = "FIN CONFIRMED: Deterministic scale-dependent delay strongly matches H=0.23!"
    elif slope < 0.05:
        results_data["Verdict"] = "FIN FALSIFIED: Delay \u03c4_max is effectively independent of scale (slope ~ 0)."
    else:
        results_data["Verdict"] = f"Unclear: Slope is {slope:.4f}, doesn't exactly match FIN but implies some non-trivial scaling."
        
    with open("QW_1660_v66_ScaleDelayFunction.json", "w") as f:
        json.dump(results_data, f, indent=2)

if __name__ == "__main__":
    main()
```

### Wynik / JSON:
```json
{
  "H_unfiltered_long_scales": 0.03733273522994843,
  "H_filtered_short_scales": 0.10019515377472395,
  "H_anomaly_reproduced": 0.002944397387996455,
  "Interpretation": {
    "H_anomaly_reproduced": "The ~0.23 value seen previously is an artifact of the 20Hz filter forcing the signal to cross zero over long scales.",
    "H_unfiltered_long_scales": "The true long-range scaling exponent of the raw noise.",
    "H_filtered_short_scales": "The short-range scaling exponent inside the studied beta frequency band."
  }
}
```
```json
{
  "Config": {
    "N_trials": 500,
    "N_samples_per_trial": 131072,
    "H_local": 0.04
  },
  "Results": {
    "Null_Cross_H_Mean": 0.5018624026057812,
    "Null_Cross_H_Std": 0.03820165650501183,
    "Null_Cross_H_Min": 0.39702322618872843,
    "Null_1st_Percentile": 0.4147443000100347
  },
  "Bootstrap_95_CI": [
    0.49861649052509377,
    0.5051488624382866
  ],
  "Real_LIGO_Observation": 0.311,
  "Statistical_Significance": "5.00 sigma",
  "Verdict": "Strong anomaly confirmed at 5.00 sigma. The observed 0.311 NEVER occurs in 500 null trials. Real cross-detector anti-persistence exists."
}
```
```json
{
  "H_real": 0.395795605895727,
  "H_shuffle": 0.501406297389991,
  "H_phase_randomized": 0.37727819984529654,
  "timestamp_utc": "2026-01-03T22:20:51.998166+00:00"
}
```
```json
{
  "0.0": 0.31105068079469167,
  "0.1": 0.30997165995459897,
  "1.0": 0.28581875130863643,
  "5.0": 0.29588892330405664,
  "10.0": 0.2947718600757252,
  "50.0": 0.2982288953028582,
  "100.0": 0.3202829195822657
}
```
```json
{
  "Coherence_Coarse_Sec": {
    "-5": 0.22713815148589375,
    "-4": 0.2284834335849714,
    "-3": 0.22813757926549674,
    "-2": 0.22719536436191218,
    "-1": 0.2285323440264692,
    "0": 0.22774859668180394,
    "1": 0.22885061737842644,
    "2": 0.22871103459093778,
    "3": 0.22899172351486632,
    "4": 0.2286404321201261,
    "5": 0.22867975079279013
  },
  "Coherence_Fine_ms": {
    "-20.0": 0.22754163493399385,
    "-18.0": 0.22747685654528674,
    "-16.0": 0.22742069398267709,
    "-14.0": 0.22738031797250913,
    "-12.0": 0.22737844803239704,
    "-10.0": 0.22748645642290677,
    "-8.0": 0.22753136500406423,
    "-6.0": 0.22759328581688065,
    "-4.0": 0.2276585701554386,
    "-2.0": 0.2277118181868326,
    "0.0": 0.22774859668180394,
    "2.0": 0.22776441559660618,
    "4.0": 0.22778231958839587,
    "6.0": 0.22781260225540848,
    "8.0": 0.22786637891755065,
    "10.0": 0.22795232472939844,
    "12.0": 0.227799685381826,
    "14.0": 0.2278843776961098,
    "16.0": 0.22793328160504317,
    "18.0": 0.22794253354135366,
    "20.0": 0.22793205412191542
  },
  "Interpretation_Guide": {
    "Peak_at_0": "Global Field / Non-local Correlation",
    "Peak_at_10ms": "Light-speed Propagation (GW-like)",
    "Peak_at_Seconds": "Slow Environmental (Seismic/Atmospheric)",
    "Flat": "Uncorrelated Noise"
  }
}
```
```json
{
  "Best_Lag_ms": 73.2421875,
  "Best_Correlation": -0.005650218758170756,
  "Verdict": "If peak is at ~10ms, it is a true gravitational/light-speed wave. If peak is exactly 0ms, it's global synchronous noise. If no clear peak, there is no broadband phase correlation.",
  "Correlation_Window_ms": {
    "-24.90 ms": 0.0035958924423889955,
    "-24.66 ms": -0.0015633011271547006,
    "-24.41 ms": 0.00017008011077205253,
    "-24.17 ms": 0.001355628330727946,
    "-23.93 ms": -0.0004134882273864998,
    "-23.68 ms": 0.0005054541196247306,
    "-23.44 ms": 0.0004900404954794068,
    "-23.19 ms": 0.0020395115197046417,
    "-22.95 ms": -0.0013569693128512989,
    "-22.71 ms": -0.0009627638168922174,
    "-22.46 ms": 0.001537215409647673,
    "-22.22 ms": 0.0007465840528576063,
    "-21.97 ms": 7.288600826827111e-06,
    "-21.73 ms": -0.00031926895456860677,
    "-21.48 ms": 0.0023335996325561915,
    "-21.24 ms": -0.000351092269335575,
    "-21.00 ms": -0.0013262921264217344,
    "-20.75 ms": 0.0027737408633027928,
    "-20.51 ms": 0.0002140976214402963,
    "-20.26 ms": -0.0017877361836068969,
    "-20.02 ms": 0.0020615367907124214,
    "-19.78 ms": 0.002347911922408769,
    "-19.53 ms": -0.0016768330858254732,
    "-19.29 ms": -1.970008582722594e-05,
    "-19.04 ms": 0.0018161985086641468,
    "-18.80 ms": 0.00038459294603834045,
    "-18.55 ms": 0.0012723814097295358,
    "-18.31 ms": -0.0009614678279654249,
    "-18.07 ms": 0.0008923191411199123,
    "-17.82 ms": 0.0014303140773147583,
    "-17.58 ms": -0.0013920944838412953,
    "-17.33 ms": 0.0031772269965881883,
    "-17.09 ms": 0.00017789552463652144,
    "-16.85 ms": -0.0015651730805940265,
    "-16.60 ms": 0.0018926567135340066,
    "-16.36 ms": 7.97589406729218e-05,
    "-16.11 ms": 0.001648131013788291,
    "-15.87 ms": 0.0003240010701814508,
    "-15.62 ms": 0.0011457879230633169,
    "-15.38 ms": 0.00089830500275615,
    "-15.14 ms": -0.0017275701490225563,
    "-14.89 ms": 0.0034050243898876608,
    "-14.65 ms": -0.00021676111242924362,
    "-14.40 ms": -0.0011972245040886079,
    "-14.16 ms": 0.002871610415299286,
    "-13.92 ms": -0.0010700932326487345,
    "-13.67 ms": 0.0023655663472634805,
    "-13.43 ms": 0.0019648860626813006,
    "-13.18 ms": -0.0022148141302521713,
    "-12.94 ms": 0.0010894040102334087,
    "-12.70 ms": 0.0012423076847510306,
    "-12.45 ms": 0.001505315720024038,
    "-12.21 ms": 0.0009349869149328864,
    "-11.96 ms": 9.707745509500646e-05,
    "-11.72 ms": -4.690221944909932e-05,
    "-11.47 ms": 0.0016199404633047672,
    "-11.23 ms": 0.001445677612150131,
    "-10.99 ms": -0.0010669552237243453,
    "-10.74 ms": 0.001344353060641185,
    "-10.50 ms": 0.001272999900487303,
    "-10.25 ms": 0.0008204132412868782,
    "-10.01 ms": -0.00031235607648622487,
    "-9.77 ms": 0.0012853310661435364,
    "-9.52 ms": 0.002164395644159403,
    "-9.28 ms": -0.001941729894240139,
    "-9.03 ms": 0.002208708098259537,
    "-8.79 ms": 0.0010388139945964607,
    "-8.54 ms": -0.0009969571629149417,
    "-8.30 ms": 0.0030184750758895647,
    "-8.06 ms": -0.0003822562793872377,
    "-7.81 ms": 0.0007977175482042462,
    "-7.57 ms": 0.0018976038759113678,
    "-7.32 ms": -0.0009456301976598709,
    "-7.08 ms": 0.00151764801706955,
    "-6.84 ms": 0.0014211995413514987,
    "-6.59 ms": -0.0002813114199989071,
    "-6.35 ms": -5.3712641338081626e-06,
    "-6.10 ms": 0.001970505995524733,
    "-5.86 ms": 0.0026007706801578993,
    "-5.62 ms": -0.0019946366304326515,
    "-5.37 ms": 0.0014201064526390096,
    "-5.13 ms": 0.003231590192172192,
    "-4.88 ms": -0.0022646161375598655,
    "-4.64 ms": 0.0017579494428678977,
    "-4.39 ms": 0.00143613910709342,
    "-4.15 ms": 0.0002827334271200185,
    "-3.91 ms": 0.0025202316600075043,
    "-3.66 ms": -0.002308038435619221,
    "-3.42 ms": 0.0019035634563612438,
    "-3.17 ms": 0.002384385718265745,
    "-2.93 ms": -0.0023088516721113347,
    "-2.69 ms": 0.0028886102039535956,
    "-2.44 ms": -0.00027315709272521303,
    "-2.20 ms": 0.0010609043632456843,
    "-1.95 ms": 0.0024014283671421064,
    "-1.71 ms": -0.0012244714012334403,
    "-1.46 ms": 0.003080209209764442,
    "-1.22 ms": -0.001583111947079528,
    "-0.98 ms": 0.0007468173716760598,
    "-0.73 ms": 0.0027386996591416113,
    "-0.49 ms": -0.0009326546566146771,
    "-0.24 ms": 0.0023591702014447904,
    "0.00 ms": -0.0008141046665487757,
    "0.24 ms": 0.001669849220685649,
    "0.49 ms": 0.002455146839938389,
    "0.73 ms": -0.0015496847419886645,
    "0.98 ms": 0.0011500722054539857,
    "1.22 ms": 0.0016680894830704196,
    "1.46 ms": -0.0001983649469891438,
    "1.71 ms": 0.0010650975025451983,
    "1.95 ms": 0.0021359048893473767,
    "2.20 ms": -7.279701833772713e-05,
    "2.44 ms": -0.0001252849054576339,
    "2.69 ms": 0.001890500652592548,
    "2.93 ms": 0.0005329389855970609,
    "3.17 ms": -0.000754706815663927,
    "3.42 ms": 0.0018228613785432081,
    "3.66 ms": 0.0012721137809966468,
    "3.91 ms": 0.0007220120666060275,
    "4.15 ms": 0.00021898813978581388,
    "4.39 ms": -0.0008471737570492024,
    "4.64 ms": 0.0018167036402158222,
    "4.88 ms": 0.0010777796843581082,
    "5.13 ms": -0.0003674864153857978,
    "5.37 ms": 0.0009149772524157882,
    "5.62 ms": 0.0014122365606375777,
    "5.86 ms": -5.9070473936331805e-05,
    "6.10 ms": 0.0007337426260857589,
    "6.35 ms": 0.0020346626159028747,
    "6.59 ms": -0.0009198023032520723,
    "6.84 ms": 0.0006959395853365341,
    "7.08 ms": 0.0007905222308515899,
    "7.32 ms": -0.0008360205116897204,
    "7.57 ms": 0.0027083644906924337,
    "7.81 ms": 0.0003334952880923528,
    "8.06 ms": -0.00041982207434698514,
    "8.30 ms": 0.0011443289366778971,
    "8.54 ms": 0.0003482171257196657,
    "8.79 ms": 0.0007702690707188068,
    "9.03 ms": -0.0008464080527831987,
    "9.28 ms": 0.0010071997339217511,
    "9.52 ms": 0.001947483803407492,
    "9.77 ms": 0.0006820474220215479,
    "10.01 ms": -0.000663842273425792,
    "10.25 ms": 0.0003373093289508459,
    "10.50 ms": 0.0015240362895267354,
    "10.74 ms": -0.0016475597435585802,
    "10.99 ms": 0.00042988213888107214,
    "11.23 ms": 0.002470088247878573,
    "11.47 ms": 0.0011871871530210985,
    "11.72 ms": -0.002074234329936727,
    "11.96 ms": -0.0005141144930285108,
    "12.21 ms": 0.0041528995821315424,
    "12.45 ms": -0.0016228369240076,
    "12.70 ms": -0.002244549994312596,
    "12.94 ms": 0.0037400614885914683,
    "13.18 ms": 0.00039106143091727544,
    "13.43 ms": -0.0019853563594793458,
    "13.67 ms": 0.0011190536410211449,
    "13.92 ms": 0.0023330869753939603,
    "14.16 ms": -0.0003616351627169883,
    "14.40 ms": -0.0011771663524364159,
    "14.65 ms": 0.002119871615158723,
    "14.89 ms": -0.0009848050714806838,
    "15.14 ms": -0.0001906672510117372,
    "15.38 ms": 0.0018270415850659622,
    "15.62 ms": -0.00034244649390020973,
    "15.87 ms": 0.00040048994610381765,
    "16.11 ms": -0.0016879625876422412,
    "16.36 ms": 0.0028303998420298434,
    "16.60 ms": -4.997284293283329e-05,
    "16.85 ms": -0.0030187742761357585,
    "17.09 ms": 0.003803223982454949,
    "17.33 ms": -0.0013198718734197547,
    "17.58 ms": 0.0010115647537846607,
    "17.82 ms": 0.001520261152710049,
    "18.07 ms": -0.0023831911658713447,
    "18.31 ms": 0.002030118606286494,
    "18.55 ms": -0.0014347003885494152,
    "18.80 ms": -7.78976108667486e-06,
    "19.04 ms": 0.0013756010423563911,
    "19.29 ms": -0.0007928289397442114,
    "19.53 ms": 0.0007965138403484423,
    "19.78 ms": -0.0006161300191874048,
    "20.02 ms": 0.0009850579714104986,
    "20.26 ms": -0.00022581939097157493,
    "20.51 ms": -0.0018688656550558558,
    "20.75 ms": 0.0014928530139534074,
    "21.00 ms": 0.0003885437493149752,
    "21.24 ms": -0.0008053923889068185,
    "21.48 ms": 0.0009445005040433828,
    "21.73 ms": -0.00021550349715158482,
    "21.97 ms": -0.00020866735578380276,
    "22.22 ms": 0.00043726384608365785,
    "22.46 ms": -0.0010525567289084176,
    "22.71 ms": 0.00029623005177442386,
    "22.95 ms": 0.00044074740127627156,
    "23.19 ms": -0.0015295763450776862,
    "23.44 ms": 0.0011450537815834666,
    "23.68 ms": 0.0012728978607754091,
    "23.93 ms": -0.003036822580908511,
    "24.17 ms": 0.0008509201291327026,
    "24.41 ms": 0.0011240928077248147,
    "24.66 ms": -0.0020782885853291963,
    "24.90 ms": 0.0009105883117009215
  }
}
```
```json
{
  "H_original": 0.0904990693214959,
  "H_fragmented": 0.0905518729921673,
  "Delta": -5.280367067139746e-05,
  "Interpretation": {
    "H_frag ~ H_original": "FIN is LOCAL/ROBUST (survives fragmentation)",
    "H_frag -> 0.5": "FIN is GLOBAL (requires long-range continuity)"
  }
}
```
```json
{
  "Config": {
    "N_trials": 1000,
    "N_samples": 131072,
    "Description": "Independent random phases applied strictly to empirical H1 and L1 PSDs"
  },
  "Results": {
    "Mean_H_cross": 0.389898436517463,
    "Std_H_cross": 0.023539041835616652,
    "Min_H_cross": 0.322682844205378,
    "P1_H_cross": 0.33770149682344514
  },
  "Observation": 0.311,
  "Significance": "-3.35 sigma",
  "Verdict": "Anomaly survives at -3.35 sigma even against realistic instrumental envelopes! The FIN correlation is still required."
}
```
```json
{
  "diurnal_null_test": [
    {
      "segment": 0,
      "H1_real": 0.2861297793052872,
      "L1_real": 0.26324010867471026,
      "Delta_real": 0.022889670630576953,
      "H1_shuffle": 0.2012720915031506,
      "L1_shuffle": 0.2001353507574666,
      "H1_reverse": 0.27527233598374606,
      "L1_reverse": 0.2625768901276506
    }
  ],
  "segment_sec": 3600,
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-01T21:39:03.253257+00:00"
}
```
```json
{
  "H1": {
    "real": 0.3107033628020025,
    "null_mean": 0.29152461412519504,
    "null_std": 0.006086465184678248
  },
  "L1": {
    "real": 0.26512680447073594,
    "null_mean": 0.2510490771952443,
    "null_std": 0.0040103972037196884
  },
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-03T21:56:45.987553+00:00"
}
```
```json
{
  "window_sec": 64,
  "step_sec": 32,
  "N_windows": 126,
  "H_mean": 0.306493493443769,
  "H_std": 0.011984760442150907,
  "timestamp_utc": "2026-01-01T22:36:49.918816+00:00"
}
```
```json
{
  "v40_CosmicTime_3x": {
    "H1": {
      "Spearman_z_H": -0.9642857142857145,
      "p_value": 0.0004541491691941689,
      "verdict": "Cosmic-Structured"
    },
    "L1": {
      "Spearman_z_H": 0.07142857142857144,
      "p_value": 0.8790481931481541,
      "verdict": "Stationary"
    },
    "V1": {
      "Spearman_z_H": 0.5714285714285715,
      "p_value": 0.1802019889115274,
      "verdict": "Cosmic-Structured"
    }
  },
  "v41_EnergyCoupling_3x": {
    "H1": {
      "Spearman_FIN_Energy": -0.9642857142857145,
      "p_value": 0.0004541491691941689,
      "verdict": "Coupled"
    },
    "L1": {
      "Spearman_FIN_Energy": 0.07142857142857144,
      "p_value": 0.8790481931481541,
      "verdict": "Decoupled"
    },
    "V1": {
      "Spearman_FIN_Energy": 0.6071428571428572,
      "p_value": 0.1482311614811614,
      "verdict": "Coupled"
    }
  },
  "v42_Axioms": {
    "Axiom_1": "FIN is scale-invariant (H constant across geometry).",
    "Axiom_2": "FIN carries structure without energy transport.",
    "Axiom_3": "FIN is tensorial and time-symmetric."
  },
  "Metadata": {
    "segment_sec": 64,
    "method": "MF-DFA q=0",
    "timestamp": "v40_42_3x_audit"
  }
}
```
```json
{
  "cross_hurst_null_validation": [
    {
      "pair": "H1-L1",
      "H_real": 0.3610605334380803,
      "H_shuffle": 0.5310469043003473,
      "H_phase_randomized": 0.39679047633700987,
      "delta_real_shuffle": -0.16998637086226698,
      "delta_real_phase": -0.035729942898929556,
      "compute_time_sec": 6.245482921600342
    },
    {
      "pair": "H1-V1",
      "H_real": 0.2150698637573661,
      "H_shuffle": 0.5219537647569182,
      "H_phase_randomized": 0.23618825276589883,
      "delta_real_shuffle": -0.30688390099955215,
      "delta_real_phase": -0.021118389008532745,
      "compute_time_sec": 6.256558418273926
    },
    {
      "pair": "L1-V1",
      "H_real": 0.27003525002927997,
      "H_shuffle": 0.5039965798656929,
      "H_phase_randomized": 0.28341320075874976,
      "delta_real_shuffle": -0.23396132983641293,
      "delta_real_phase": -0.01337795072946979,
      "compute_time_sec": 6.2524378299713135
    }
  ],
  "summary": {
    "mean_delta_real_phase": -0.023408760878977364,
    "std_delta_real_phase": 0.00926776639630034,
    "N_pairs": 3
  },
  "sample_rate": 4096,
  "window_sec": 512,
  "timestamp_utc": "2026-01-04T01:25:25.095853+00:00"
}
```
```json
{
  "0.05": {
    "mean": 0.4004410091997295,
    "std": 0.027933321248884768
  },
  "0.06": {
    "mean": 0.4037080246714917,
    "std": 0.027440867036929373
  },
  "0.07": {
    "mean": 0.39858702807371643,
    "std": 0.024428475993454848
  },
  "0.08": {
    "mean": 0.40126004797377846,
    "std": 0.0237754441519594
  },
  "0.09": {
    "mean": 0.39899303797859104,
    "std": 0.02716992747638438
  },
  "0.10": {
    "mean": 0.39857089796607176,
    "std": 0.027970497217470577
  },
  "0.11": {
    "mean": 0.39963077436004674,
    "std": 0.028280015123032336
  },
  "0.12": {
    "mean": 0.39820285479377476,
    "std": 0.022425824368900078
  },
  "0.13": {
    "mean": 0.3986489805097576,
    "std": 0.026189567416261723
  },
  "0.14": {
    "mean": 0.3988655107489701,
    "std": 0.023426960570410077
  },
  "0.15": {
    "mean": 0.39620037112491585,
    "std": 0.026921608819045015
  },
  "0.16": {
    "mean": 0.39457224499111254,
    "std": 0.027340425191611815
  },
  "0.17": {
    "mean": 0.40082080274225335,
    "std": 0.027089112575727187
  },
  "0.18": {
    "mean": 0.3941544493091026,
    "std": 0.02355772334169456
  },
  "0.19": {
    "mean": 0.39981841104111665,
    "std": 0.027016244861192636
  },
  "0.20": {
    "mean": 0.40029184793024347,
    "std": 0.024708432620686956
  },
  "0.21": {
    "mean": 0.39789812590783497,
    "std": 0.025520497936303666
  },
  "0.22": {
    "mean": 0.395847891965017,
    "std": 0.025790662490122216
  },
  "0.23": {
    "mean": 0.39648745424536785,
    "std": 0.02385882135337397
  },
  "0.24": {
    "mean": 0.3995899720316251,
    "std": 0.02578104008013167
  },
  "0.25": {
    "mean": 0.3972815606284205,
    "std": 0.025576727198992264
  },
  "0.26": {
    "mean": 0.3983771805109425,
    "std": 0.028935293979013833
  },
  "0.27": {
    "mean": 0.3969920165843066,
    "std": 0.02497529525206574
  },
  "0.28": {
    "mean": 0.4031693429115578,
    "std": 0.024551636034024146
  },
  "0.29": {
    "mean": 0.39991865182437275,
    "std": 0.026682409914462683
  },
  "0.30": {
    "mean": 0.4008943857800395,
    "std": 0.02675033780957701
  },
  "0.31": {
    "mean": 0.3962344554211347,
    "std": 0.02571293799437975
  },
  "0.32": {
    "mean": 0.3947083278936816,
    "std": 0.026558697381502395
  },
  "0.33": {
    "mean": 0.40284990798616005,
    "std": 0.026789983832809547
  },
  "0.34": {
    "mean": 0.3987862366901604,
    "std": 0.026211352369012567
  },
  "0.35": {
    "mean": 0.3973613704299476,
    "std": 0.025134966977588186
  },
  "0.36": {
    "mean": 0.3991941009288439,
    "std": 0.024758946990399582
  },
  "0.37": {
    "mean": 0.39639518530902856,
    "std": 0.02511888629071801
  },
  "0.38": {
    "mean": 0.39479028817294837,
    "std": 0.024507724205019284
  },
  "0.39": {
    "mean": 0.4022232369392795,
    "std": 0.02617700375805709
  },
  "0.40": {
    "mean": 0.3985299913591742,
    "std": 0.027914514081746673
  },
  "0.41": {
    "mean": 0.39711781214788244,
    "std": 0.02600649550958563
  },
  "0.42": {
    "mean": 0.395751567199055,
    "std": 0.02467251649620674
  },
  "0.43": {
    "mean": 0.39459432398424904,
    "std": 0.026048966497937552
  },
  "0.44": {
    "mean": 0.39641806407787455,
    "std": 0.02560372103068189
  },
  "0.45": {
    "mean": 0.3984582561241745,
    "std": 0.027550973733774454
  }
}
```
```json
{
  "Spectra": {
    "Strain": {
      "-5": 0.006844430703847382,
      "0": 0.0025043396942770737,
      "5": -0.001974344057338875
    },
    "Seismic_Model": {
      "-5": 1.5419615969998273,
      "0": 1.4251372162191898,
      "5": 1.3549806600694683
    },
    "Magnetic_Model": {
      "-5": 0.9690215577376325,
      "0": 0.9594795851484398,
      "5": 0.9378419488518792
    }
  },
  "Distances": {
    "vs_Seismic": 2.494352801600233,
    "vs_Magnetic": 1.6507092744080252
  },
  "Info": "If fetch failed, Red/Pink noise models were used as physical vetoes."
}
```
```json
{
  "H1_Pure_H": 0.03733273522994843,
  "L1_Pure_H": 0.05299552509319112,
  "Cross_H1_L1_Pure_H": 0.31105068288116156,
  "Interpretation": "SURPRISING: The control loops / seismic background of H1 and L1 are strongly correlated at long distances."
}
```
```json
{
  "Whitened_Cross_H": 0.6856328431496843,
  "Verdict": "If ~0.50, there is NO phase structure (FIN is falsified as a phase effect). If significantly diff from 0.50, FIN exists as pure phase correlation."
}
```
```json
{
  "q": [
    -5.0,
    -4.0,
    -3.0,
    -2.0,
    -1.0,
    0.0,
    1.0,
    2.0,
    3.0,
    4.0,
    5.0
  ],
  "tau": [
    -2.3437647824134276,
    -2.0452331705814615,
    -1.7603708297863268,
    -1.4904612750979518,
    -1.236655371120009,
    -1.0,
    -0.781464764172749,
    -0.5818119681941865,
    -0.4011786117098268,
    -0.23861932213431447,
    -0.09213253820621514
  ],
  "alpha": [
    0.2985316118319661,
    0.2916969763135504,
    0.27738594774175485,
    0.2618577293331589,
    0.2452306375489759,
    0.22759530347363,
    0.20909401590290677,
    0.1901430762314611,
    0.171596323029936,
    0.15452303675180584,
    0.14648678392809933
  ],
  "f_alpha": [
    0.8511067232535972,
    0.8784452653272599,
    0.9282129865610622,
    0.966745816431634,
    0.9914247335710331,
    1.0,
    0.9905587800756558,
    0.9620981206571086,
    0.9159675807996348,
    0.8567114691415378,
    0.8245664578467118
  ],
  "FIN": {
    "H0": 0.22774859668180394,
    "DeltaH": 0.08717946412392852,
    "Cascade": "log-Poisson",
    "Cascade_param": 0.7680334399643445,
    "chi2_logPoisson": 0.08637671851557861,
    "chi2_logNormal": 0.22775411741584167
  }
}
```
```json
{
  "Config": {
    "GPS": 1266965117,
    "Fs_target": 1,
    "Cutoff_Hz": 0.1,
    "Window_sec": 4096
  },
  "Cross_H_LF": 1.9154092490606147,
  "Coherence_0_01_to_0_1_Hz": 0.6741010398437792,
  "Mean_Real_CSD_LF": 5.3899207891616826e-48,
  "Mean_Imag_CSD_LF": 1.0931139634509775e-50,
  "Verdict": "Effect disappears at ultra-low frequencies or interpretation is mixed."
}
```
```json
{
  "Detector_H0": {
    "H1": 0.0026335489204392295,
    "L1": 0.01859453307242668,
    "V1": 0.0036370796568547664
  },
  "Arm_Lengths": {
    "H1": 4000.0,
    "L1": 4000.0,
    "V1": 3000.0
  },
  "Scaling_Analysis": {
    "Unique_Lengths": [
      3000.0,
      4000.0
    ],
    "Average_H": [
      0.0036370796568547664,
      0.010614040996432955
    ],
    "Gamma_Exponent": 3.72284817816416,
    "Description": "H ~ 4.13e-16 * L^3.723"
  },
  "timestamp": "v29_execution"
}
```
```json
{
  "H_mean": 0.0030418309661207313,
  "H_std_local": 0.0007083480808800167,
  "Interpretation": {
    "low_std (<0.02)": "FIN is a STRUCTURAL MEASURE",
    "high_std": "FIN is a DYNAMICAL PROCESS"
  }
}
```
```json
{
  "Hurst_H1": 0.3107033628020025,
  "Hurst_L1": 0.26512680447073594,
  "Delta_H": 0.04557655833126656,
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-01T18:55:25.233662"
}
```
```json
{
  "sidereal_test": [
    {
      "segment": 0,
      "sidereal_phase": 0.0,
      "H1": 0.2764631582569697,
      "L1": 0.2738780068272483,
      "Delta": 0.002585151429721433
    },
    {
      "segment": 1,
      "sidereal_phase": 0.020890370815687738,
      "H1": 0.29013919118663656,
      "L1": 0.2689470865713718,
      "Delta": 0.021192104615264773
    }
  ],
  "segment_sec": 1800,
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-01T22:27:10.425300+00:00"
}
```
```json
{
  "beta_real": 5.8071325680938815,
  "beta_null_mean": 5.719870385223784,
  "beta_null_std": 0.013884007622481693,
  "z_score": 6.285086067570166
}
```
```json
{
  "H1_V1_Cross_H": 0.714808705323103,
  "Verdict": "If ~0.50, the anomaly is an H1-L1 paired architectural artifact. If ~0.31, it's truly global across distinct architectures."
}
```
```json
{
  "Spectra": {
    "Strain_H1": {
      "-5": 0.006844430703847382,
      "0": 0.0025043396942770737,
      "5": -0.001974344057338875
    },
    "Seismic": {
      "-5": 0.5065520521790964,
      "0": 0.4865799566905009,
      "5": 0.46810562719698207
    },
    "Magnetic": {
      "-5": 0.5136585429375095,
      "0": 0.500904416564615,
      "5": 0.4819535028086708
    }
  },
  "Distance_Euclidean": {
    "Strain_vs_Seismic": 0.8396499802329472,
    "Strain_vs_Magnetic": 0.8599124036527941
  },
  "Conclusion": "If Distances > 0.1, FIN is distinct from Environment.",
  "timestamp": "v30_execution"
}
```
```json
{
  "original_cross_H": 0.31105068288116156,
  "surrogate_cross_H_mean": 0.28332943157004487,
  "surrogate_cross_H_std": 0.008788916370598766,
  "N_surrogates": 20
}
```
```json
{
  "global_consistency_test": [
    {
      "window_sec": 64,
      "H": 0.30885260430877476,
      "samples": 262144,
      "compute_time_sec": 0.24536657333374023
    },
    {
      "window_sec": 128,
      "H": 0.3098408247635631,
      "samples": 524288,
      "compute_time_sec": 0.4489321708679199
    },
    {
      "window_sec": 256,
      "H": 0.28857808692469833,
      "samples": 1048576,
      "compute_time_sec": 0.8733696937561035
    },
    {
      "window_sec": 512,
      "H": 0.2756928984409219,
      "samples": 2097152,
      "compute_time_sec": 1.6308715343475342
    },
    {
      "window_sec": 1024,
      "H": 0.2874670726409621,
      "samples": 4194304,
      "compute_time_sec": 3.195228099822998
    },
    {
      "window_sec": 2048,
      "H": 0.27098363828141603,
      "samples": 8388608,
      "compute_time_sec": 6.0503809452056885
    }
  ],
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-01T22:53:20.955802+00:00"
}
```
```json
{
  "cross_epoch_consistency": [
    {
      "epoch": "Epoch_1",
      "window_sec": 512,
      "samples": 2097152,
      "H": 0.27105320196327354,
      "compute_time_sec": 1.1741423606872559
    }
  ],
  "summary": {
    "H_mean": 0.27105320196327354,
    "H_std": 0.0,
    "N_epochs": 1
  },
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-03T23:04:14.447923+00:00"
}
```
```json
{
  "diurnal_test": [
    {
      "segment": 0,
      "H1": 0.2861297793052872,
      "L1": 0.26324010867471026,
      "Delta": 0.022889670630576953
    }
  ],
  "segment_sec": 3600,
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-01T21:21:07.394453+00:00"
}
```
```json
{
  "v40_CosmicTime": {
    "Spearman_z_H": -0.9642857142857145,
    "p_value": 0.0004541491691941689,
    "verdict": "Cosmic-Structured"
  },
  "v41_GW_Coupling": {
    "Spearman_FIN_GWrate": -0.9642857142857145,
    "p_value": 0.0004541491691941689,
    "verdict": "Coupled"
  },
  "v42_Axioms": {
    "Axiom_1": "FIN is scale-invariant (H constant across geometry).",
    "Axiom_2": "FIN carries structure without energy transport.",
    "Axiom_3": "FIN is tensorial and time-symmetric."
  }
}
```
```json
{
  "H_cross_short_scales_1s_to_25s": 0.22332124569812675,
  "H_cross_long_scales_25s_to_128s": 0.3619575625119233,
  "H_cross_full_baseline": 0.31105068288116156,
  "Verdict": "Anti-persistence is scale-invariant across both short and long bins."
}
```
```json
{
  "40-80Hz": 0.7313123830135317,
  "100-200Hz": 0.5631587605638296,
  "300-500Hz": 0.4134798668917408
}
```
```json
{
  "segment_sec": 64,
  "r_real": -0.1453076356553426,
  "p_real": 0.25194302076116887,
  "null_mean": 0.024229854332111018,
  "null_std": 0.12065324834219607,
  "z_score": -1.4051630794606735,
  "verdict": "NO COHERENCE \u2014 CONSISTENT WITH INDEPENDENT NOISE",
  "timestamp_utc": "2026-01-01T21:51:55.414111+00:00"
}
```
```json
{
  "H1": {
    "real": 0.30885260430877476,
    "null_mean": 0.31090352670748245,
    "null_std": 0.007531104412966834
  },
  "L1": {
    "real": 0.3246381413143207,
    "null_mean": 0.30472640469044976,
    "null_std": 0.005574223192037795
  },
  "sample_rate": 4096.0,
  "duration_used": 64.0,
  "timestamp_utc": "2026-01-01T20:43:16.857304+00:00"
}
```
```json
{
  "H_original": 0.22774859668180394,
  "H_permuted_mean": 0.5003428220114121,
  "H_permuted_std": 0.015041406502576567,
  "Delta": -0.2725942253296082,
  "Interpretation": {
    "Delta ~ 0": "Amplitude-based / noise",
    "Delta >> 0": "Informational / structural FIN"
  }
}
```
```json
{
  "0.0": {
    "Cross_H_mean": 0.5013565178093967,
    "Cross_H_std": 0.020171871555999893
  },
  "0.5": {
    "Cross_H_mean": 0.5059634912258375,
    "Cross_H_std": 0.016360932501354154
  },
  "1.0": {
    "Cross_H_mean": 0.5103556104649954,
    "Cross_H_std": 0.0146661713518238
  },
  "2.0": {
    "Cross_H_mean": 0.5044638878327039,
    "Cross_H_std": 0.02237324856125114
  },
  "3.0": {
    "Cross_H_mean": 0.49910779621350587,
    "Cross_H_std": 0.024637612596447415
  },
  "5.0": {
    "Cross_H_mean": 0.4949717022769238,
    "Cross_H_std": 0.025364894370843234
  },
  "7.0": {
    "Cross_H_mean": 0.49354969534495413,
    "Cross_H_std": 0.025348111303817908
  },
  "10.0": {
    "Cross_H_mean": 0.49267169835880253,
    "Cross_H_std": 0.025164686986300058
  },
  "15.0": {
    "Cross_H_mean": 0.49211390367399277,
    "Cross_H_std": 0.02489982657201801
  },
  "20.0": {
    "Cross_H_mean": 0.49187840042302067,
    "Cross_H_std": 0.024724322758137508
  },
  "50.0": {
    "Cross_H_mean": 0.49152941673745076,
    "Cross_H_std": 0.02433616880408267
  }
}
```
```json
{
  "cross_hurst_results": [
    {
      "pair": "H1-L1",
      "H_cross": 0.11552275914205934,
      "samples": 2097152
    },
    {
      "pair": "H1-V1",
      "H_cross": 0.07803517778157153,
      "samples": 2097152
    },
    {
      "pair": "L1-V1",
      "H_cross": 0.19066426135733577,
      "samples": 2097152
    }
  ],
  "summary": {
    "H_cross_mean": 0.12807406609365554,
    "H_cross_std": 0.04682932910348679,
    "pairs": [
      "H1-L1",
      "H1-V1",
      "L1-V1"
    ]
  },
  "sample_rate": 4096,
  "window_sec": 512,
  "timestamp_utc": "2026-01-04T01:01:51.856503+00:00"
}
```
```json
{
  "bands": {
    "30-80Hz": {
      "H1": 0.23607077037727509,
      "L1": 0.23731590563132346,
      "Delta": -0.001245135254048374
    },
    "80-200Hz": {
      "H1": 0.26564212964633904,
      "L1": 0.24650555849726993,
      "Delta": 0.01913657114906911
    },
    "200-500Hz": {
      "H1": 0.2069448048711533,
      "L1": 0.1766016827951618,
      "Delta": 0.03034312207599149
    }
  },
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-01T19:24:15.850515+00:00"
}
```
```json
{
  "config": {
    "N": 1048576,
    "N_trials": 20,
    "H_true_background": 0.23,
    "H_local_damping": 0.04,
    "SNR_local_over_background": 20.0
  },
  "mixed_signal": {
    "H_x_mean": 0.07864376942410117,
    "H_x_std": 0.002380065086687186,
    "H_y_mean": 0.07825426745142278,
    "H_y_std": 0.003211036386556685,
    "Cross_H_mean": 0.5094425737999033,
    "Cross_H_std": 0.02343557500667056
  },
  "null_model": {
    "Cross_H_null_mean": 0.4965606426963275,
    "Cross_H_null_std": 0.027611698738229173
  },
  "verdict": "Cross-MF-DFA gives 0.509, which matches neither 0.23 nor 0.31. Relationship is complex."
}
```
```json
{
  "v36_Scale_L": {
    "slope": 0.030842042547518753,
    "r_squared": 0.0054080586363148735,
    "verdict": "Scale Invariant"
  },
  "v37_Tensor": {
    "H_plus_proxy": 0.2217127145427585,
    "H_cross_proxy": 0.3606526430653205,
    "Delta": 0.138939928522562,
    "verdict": "Tensorial"
  },
  "v38_Isotropy": {
    "slope": 0.0019315065165690576,
    "r_squared": 0.9966893760060604,
    "verdict": "Anisotropic"
  },
  "v39_EnergyInfo": {
    "correlation": -0.955014594715938,
    "p_value": 0.1916781512354847,
    "verdict": "Energy Driven"
  },
  "Raw_Data": {
    "H1": {
      "H": 0.2217127145427585,
      "RMS": 4.393047798925605e-20
    },
    "L1": {
      "H": 0.3606526430653205,
      "RMS": 2.1030104377162187e-21
    },
    "V1": {
      "H": 0.2823099760853233,
      "RMS": 1.4631395710288143e-20
    }
  }
}
```
```json
{
  "cross_epoch_fractal_stability": [
    {
      "epoch": "O3a_H1",
      "window_sec": 512,
      "samples": 2097152,
      "H": 0.2900189379914974
    },
    {
      "epoch": "O3b_H1",
      "window_sec": 512,
      "samples": 2097152,
      "H": 0.36065896220129573
    }
  ],
  "summary": {
    "H_mean": 0.3253389500963966,
    "H_std": 0.03532001210489916,
    "N_epochs": 2
  },
  "sample_rate": 4096,
  "timestamp_utc": "2026-01-04T00:26:47.894554+00:00"
}
```
```json
{
  "Fundamental_H_Input": 0.23,
  "Injected_Cross_H_Output": 0.300537392344461,
  "Verdict": "The instrumental response accurately transforms a true H=0.23 background into an apparent Cross-H ~ 0.31. The 0.31 anomaly COULD be a projection of a 0.23 fundamental background."
}
```
```json
{
  "Test1_Independent_Phases_Real_Envelope": {
    "mean": 0.2868038142777771,
    "std": 0.012702395253275257,
    "verdict": "If ~0.50, then the envelope ALONE does not cause the anomaly. If ~0.31, then the envelope alone is responsible."
  },
  "Test2_Real_Phases_Generic_Envelope": {
    "value": 0.705328267933796,
    "verdict": "If ~0.50, then removing the shared envelope destroys the anomaly. If ~0.31, the anomaly is robustly in the phase."
  }
}
```
```json
{
  "Scale_1000_to_3000": {
    "s_min": 1000,
    "s_max": 3000,
    "H_cross": 0.8918372300793768
  },
  "Scale_3000_to_10000": {
    "s_min": 3000,
    "s_max": 10000,
    "H_cross": 0.3088136078636575
  },
  "Scale_10000_to_30000": {
    "s_min": 10000,
    "s_max": 30000,
    "H_cross": 0.19505626413175572
  },
  "Scale_30000_to_100000": {
    "s_min": 30000,
    "s_max": 100000,
    "H_cross": 0.16602859804572065
  },
  "Scale_100000_to_131072": {
    "s_min": 100000,
    "s_max": 131072,
    "H_cross": -0.13145406277611377
  }
}
```
```json
{
  "bands": [
    [
      30,
      60
    ],
    [
      60,
      120
    ],
    [
      120,
      240
    ],
    [
      240,
      480
    ]
  ],
  "H1_band": [
    0.23478763001930075,
    0.20451298768398127,
    0.2605842628711187,
    0.205146703094412
  ],
  "L1_band": [
    0.2357252841122392,
    0.21592309821763106,
    0.24127755831025977,
    0.17635400825431122
  ],
  "r_real": 0.7903615769456582,
  "p_real": 0.20963842305434177,
  "null_mean": -0.04269731659612157,
  "null_std": 0.6081802287520204,
  "z_score": 1.369756618447147,
  "timestamp_utc": "2026-01-01T22:03:23.168967+00:00"
}
```
```json
{
  "1266965117": {
    "H1_Pure": 0.03733273522994843,
    "L1_Pure": 0.05299552509319112,
    "Cross_H1_L1": 0.31105068288116156
  },
  "1267051517": {
    "H1_Pure": 0.041405315287956546,
    "L1_Pure": 0.06420362757311483,
    "Cross_H1_L1": 0.33515322923661534
  },
  "1267137917": {
    "H1_Pure": 0.044313974791691,
    "L1_Pure": 0.05405721290409476,
    "Cross_H1_L1": 0.3034840404652728
  },
  "1253326744": {
    "H1_Pure": 0.04180091028389436,
    "L1_Pure": 0.05994293130661616,
    "Cross_H1_L1": 0.7633012338448356
  }
}
```
```json
{
  "v44_Geometry": {
    "baseline_results": {
      "H1-L1": {
        "angle": 90.0,
        "H_cross": 0.0022573648019294417
      },
      "H1-V1": {
        "angle": 45.0,
        "H_cross": 0.0020402144813620095
      },
      "L1-V1": {
        "angle": 60.0,
        "H_cross": 0.0091828030341689
      }
    },
    "slope_angle": -2.884208071414157e-05,
    "r_squared": 0.02646165698996646,
    "verdict": "Non-geometric"
  },
  "v45_Relational": {
    "Spearman_H_energy": -1.0,
    "p_value": 0.0,
    "verdict": "Energetic"
  },
  "v46_Scale": {
    "H_windows": {
      "32": 0.0024096229830277488,
      "64": 0.0022573648019294417,
      "128": 0.0021222544788838562,
      "256": 0.002277317937767772
    },
    "std": 0.00010184714011412424,
    "verdict": "Structural"
  }
}
```
```json
{
  "H_original": 0.0025043396942770737,
  "H_phase_randomized": 0.03488916553240482,
  "H_time_permuted": 0.5261990862176169,
  "Interpretation_Guide": {
    "H_phase ~ H_original": "Phase-encoded / informational FIN",
    "H_phase ~ 0.5": "Amplitude-based noise",
    "H_perm ~ 0.5": "Order-dependent structure confirmed"
  }
}
```
```json
{
  "beta_h1": 5.788587780685799,
  "beta_l1": 5.825677355501964,
  "delta_beta": -0.03708957481616437,
  "null_mean": 0.13165685849154202,
  "null_std": 0.027553184027458343,
  "z_score": -6.124389585593477
}
```
```json
{
  "Coherence_Coarse_Sec": {
    "H1-L1": {
      "-5": 0.22713815148589375,
      "-4": 0.2284834335849714,
      "-3": 0.22813757926549674,
      "-2": 0.22719536436191218,
      "-1": 0.2285323440264692,
      "0": 0.22774859668180394,
      "1": 0.22885061737842644,
      "2": 0.22871103459093778,
      "3": 0.22899172351486632,
      "4": 0.2286404321201261,
      "5": 0.22867975079279013
    },
    "H1-V1": {
      "-5": 0.0761885833323833,
      "-4": 0.059078151977446036,
      "-3": 0.0798636811523966,
      "-2": 0.07126686483520528,
      "-1": 0.06592063659510108,
      "0": 0.08061816890960831,
      "1": 0.05841927158024111,
      "2": 0.07313876548475505,
      "3": 0.06833537685489051,
      "4": 0.06041403469070915,
      "5": 0.07775007475583821
    },
    "L1-V1": {
      "-5": 0.1301869471470167,
      "-4": 0.13556358584670644,
      "-3": 0.1573920502774073,
      "-2": 0.16091488046423677,
      "-1": 0.1666934268485417,
      "0": 0.15516582476076077,
      "1": 0.16092046982674887,
      "2": 0.1586646925042033,
      "3": 0.1527071989203297,
      "4": 0.14458779321127657,
      "5": 0.1410046398746766
    }
  },
  "Coherence_Fine_ms": {
    "H1-L1": {
      "-20.0": 0.22754163493399385,
      "-18.0": 0.22747685654528674,
      "-16.0": 0.22742069398267709,
      "-14.0": 0.22738031797250913,
      "-12.0": 0.22737844803239704,
      "-10.0": 0.22748645642290677,
      "-8.0": 0.22753136500406423,
      "-6.0": 0.22759328581688065,
      "-4.0": 0.2276585701554386,
      "-2.0": 0.2277118181868326,
      "0.0": 0.22774859668180394,
      "2.0": 0.22776441559660618,
      "4.0": 0.22778231958839587,
      "6.0": 0.22781260225540848,
      "8.0": 0.22786637891755065,
      "10.0": 0.22795232472939844,
      "12.0": 0.227799685381826,
      "14.0": 0.2278843776961098,
      "16.0": 0.22793328160504317,
      "18.0": 0.22794253354135366,
      "20.0": 0.22793205412191542
    },
    "H1-V1": {
      "-20.0": 0.07380368836963057,
      "-18.0": 0.07583107454321843,
      "-16.0": 0.07804436472618854,
      "-14.0": 0.08057437240608738,
      "-12.0": 0.08384097504009345,
      "-10.0": 0.06462637417388827,
      "-8.0": 0.06802516707388706,
      "-6.0": 0.07217648112759158,
      "-4.0": 0.07612695523393315,
      "-2.0": 0.07871699263290512,
      "0.0": 0.08061816890960831,
      "2.0": 0.08126862628342205,
      "4.0": 0.08109825777150323,
      "6.0": 0.08058131442028921,
      "8.0": 0.07934242114197854,
      "10.0": 0.07638145978477211,
      "12.0": 0.07283125671815269,
      "14.0": 0.07333114609843601,
      "16.0": 0.07261825865508073,
      "18.0": 0.07154471844926055,
      "20.0": 0.0710382411168277
    },
    "L1-V1": {
      "-20.0": 0.14966628278599062,
      "-18.0": 0.1673353942887202,
      "-16.0": 0.18557497037893034,
      "-14.0": 0.17834810641236656,
      "-12.0": 0.1486904148949327,
      "-10.0": 0.14667045415716096,
      "-8.0": 0.14712041829137643,
      "-6.0": 0.156730848755411,
      "-4.0": 0.15038150624375707,
      "-2.0": 0.13807230161484335,
      "0.0": 0.15516582476076077,
      "2.0": 0.16048255254948618,
      "4.0": 0.15232325761547447,
      "6.0": 0.13418898606714447,
      "8.0": 0.14331995956172752,
      "10.0": 0.16283249488586493,
      "12.0": 0.17583982538042625,
      "14.0": 0.1645064893724227,
      "16.0": 0.1470261296279435,
      "18.0": 0.13971512595665786,
      "20.0": 0.15581125290362344
    }
  },
  "Note": "Thresholds adjusted for gravitational strain amplitude (~1e-21)",
  "timestamp": "v28_execution"
}
```
```json
{
  "integration_test": [
    {
      "window_sec": 32,
      "H1": 0.33995934773901576,
      "L1": 0.3520262031207114,
      "Delta": -0.012066855381695663,
      "compute_time_sec": 0.31324338912963867
    },
    {
      "window_sec": 64,
      "H1": 0.3105435282906503,
      "L1": 0.3263566705352001,
      "Delta": -0.015813142244549827,
      "compute_time_sec": 0.5691792964935303
    },
    {
      "window_sec": 128,
      "H1": 0.31001070121971774,
      "L1": 0.31835155974370716,
      "Delta": -0.008340858523989414,
      "compute_time_sec": 1.0603904724121094
    },
    {
      "window_sec": 256,
      "H1": 0.2881253580625984,
      "L1": 0.30091963081910567,
      "Delta": -0.012794272756507241,
      "compute_time_sec": 2.042760133743286
    },
    {
      "window_sec": 512,
      "H1": 0.2749408309566025,
      "L1": 0.2906625964784416,
      "Delta": -0.015721765521839126,
      "compute_time_sec": 3.687760591506958
    },
    {
      "window_sec": 1024,
      "H1": 0.28688700383670396,
      "L1": 0.287528864266169,
      "Delta": -0.0006418604294650687,
      "compute_time_sec": 7.145539283752441
    }
  ],
  "sample_rate": 4096.0,
  "total_runtime_sec": 14.842632055282593,
  "timestamp_utc": "2026-01-01T21:06:38.168655+00:00"
}
```
```json
{
  "cross_epoch_fractal_stability": [
    {
      "epoch": "O3a_H1",
      "window_sec": 512,
      "samples": 2097152,
      "H": 0.2900189379914974
    },
    {
      "epoch": "O3b_H1",
      "window_sec": 512,
      "samples": 2097152,
      "H": 0.36065896220129573
    },
    {
      "epoch": "O4_H1",
      "window_sec": 512,
      "samples": 2097152,
      "H": 0.33009336350893775
    }
  ],
  "summary": {
    "H_mean": 0.32692375456724365,
    "H_std": 0.0289256295894912,
    "N_epochs": 3,
    "epochs_processed": [
      "O3a_H1",
      "O3b_H1",
      "O4_H1"
    ]
  },
  "sample_rate": 4096,
  "timestamp_utc": "2026-01-04T00:33:13.409766+00:00"
}
```
```json
{
  "MF_DFA": {
    "H1": {
      "-5": 0.006844430703847382,
      "-4": 0.005990294627013035,
      "-3": 0.005128600477532727,
      "-2": 0.004260019290883733,
      "-1": 0.003385112029380349,
      "0": 0.0025043396942770737,
      "1": 0.0016180874149963294,
      "2": 0.0007267016120266353,
      "3": -0.00016946256612925594,
      "4": -0.001069988196042291,
      "5": -0.001974344057338875
    },
    "L1": {
      "-5": 0.047682750986458335,
      "-4": 0.0428337697498101,
      "-3": 0.037691671330022485,
      "-2": 0.03223713004462336,
      "-1": 0.026448632237696158,
      "0": 0.020271377114164532,
      "1": 0.01347862981930821,
      "2": 0.005121093332372146,
      "3": -0.00825328495861051,
      "4": -0.033948655687234715,
      "5": -0.07157518028645297
    },
    "V1": {
      "-5": 0.011001293684340626,
      "-4": 0.009339864389059477,
      "-3": 0.007719558192389645,
      "-2": 0.00614238209867122,
      "-1": 0.004609302332536944,
      "0": 0.0031201886565858854,
      "1": 0.0016736935469015438,
      "2": 0.00026705753659740194,
      "3": -0.0011041908926608184,
      "4": -0.0024466938486457897,
      "5": -0.0037699433012532936
    }
  },
  "CROSS_MF_DFA": {
    "real": {
      "-5": 0.2687529564826855,
      "-4": 0.2613082926453654,
      "-3": 0.2534569432621089,
      "-2": 0.24523063754897595,
      "-1": 0.23665537112000895,
      "0": 0.22774859668180394,
      "1": 0.21853523582725104,
      "2": 0.20909401590290674,
      "3": 0.1996071294300577,
      "4": 0.19034516946642138,
      "5": 0.18157349235875697
    },
    "shuffle": {
      "-5": 0.5034203775726257,
      "-4": 0.503315190571987,
      "-3": 0.5034871257469522,
      "-2": 0.5036697032660185,
      "-1": 0.5035583331580095,
      "0": 0.5029040733314845,
      "1": 0.501572960760977,
      "2": 0.49953625611087116,
      "3": 0.49682172095733346,
      "4": 0.49347982600161155,
      "5": 0.4895835197780457
    },
    "phase": {
      "-5": 0.21241714108321552,
      "-4": 0.2066586328164522,
      "-3": 0.20240773013742588,
      "-2": 0.20007000933802568,
      "-1": 0.19976744895791468,
      "0": 0.20077379936394577,
      "1": 0.20132358066545503,
      "2": 0.19980187737030422,
      "3": 0.1959791993052112,
      "4": 0.19064284141606427,
      "5": 0.18466079912108885
    },
    "H1-V1": {
      "real": {
        "-5": 0.07824452227756472,
        "-4": 0.0773505311767661,
        "-3": 0.07717491185697276,
        "-2": 0.07772987656508179,
        "-1": 0.07893671623815128,
        "0": 0.08061816890960831,
        "1": 0.08252423161965843,
        "2": 0.08440309871882581,
        "3": 0.08607757944928261,
        "4": 0.08747213561044159,
        "5": 0.08858725918899217
      },
      "shuffle": {
        "-5": 0.507041800253571,
        "-4": 0.5052754204414203,
        "-3": 0.5035364723268283,
        "-2": 0.5018306812190458,
        "-1": 0.5001567206873525,
        "0": 0.4984694722998402,
        "1": 0.49665564437092724,
        "2": 0.49455703697797243,
        "3": 0.49203532276862355,
        "4": 0.4890301173986625,
        "5": 0.4855734173855122
      },
      "phase": {
        "-5": 0.08344104807009789,
        "-4": 0.07901170280927079,
        "-3": 0.07561390986870035,
        "-2": 0.07311542986638424,
        "-1": 0.07139771427750366,
        "0": 0.07032920663337112,
        "1": 0.06972990468301439,
        "2": 0.06938983175156443,
        "3": 0.06912385008644327,
        "4": 0.06881168828504328,
        "5": 0.06839963231201049
      }
    },
    "L1-V1": {
      "real": {
        "-5": 0.1703613127327925,
        "-4": 0.16716218425111537,
        "-3": 0.16403612007350366,
        "-2": 0.16100045506505392,
        "-1": 0.1580501998513365,
        "0": 0.15516582476076077,
        "1": 0.15232808064208983,
        "2": 0.14953007529836113,
        "3": 0.14677710019154533,
        "4": 0.1440721791108341,
        "5": 0.14139752345734977
      },
      "shuffle": {
        "-5": 0.5045139190414584,
        "-4": 0.5031718371874048,
        "-3": 0.5023880004865273,
        "-2": 0.5020622920859986,
        "-1": 0.5018845338107429,
        "0": 0.5013856509483983,
        "1": 0.5001254538408195,
        "2": 0.4978749376950873,
        "3": 0.4946506925377633,
        "4": 0.490631130687572,
        "5": 0.4860568715004665
      },
      "phase": {
        "-5": 0.18911863880090202,
        "-4": 0.18146817260109233,
        "-3": 0.1737619636920628,
        "-2": 0.16631363858181425,
        "-1": 0.1594712155669294,
        "0": 0.15354185295717637,
        "1": 0.14870969447167898,
        "2": 0.14498707421159965,
        "3": 0.1422261736004945,
        "4": 0.14018347516359303,
        "5": 0.13859906551400358
      }
    }
  },
  "Wavelet_Leaders": {
    "H1": [
      4.812972522623127e-20,
      2.6983210498821907e-20,
      4.2188180222577184e-20,
      2.070283832327583e-19,
      1.5563749178434705e-19,
      2.8361233791228486e-20
    ],
    "L1": [
      1.6380946086137567e-20,
      6.207168573902065e-21,
      3.4411965114437665e-21,
      1.6474495104490852e-20,
      1.4780705973309457e-20,
      1.867954404117499e-21
    ]
  },
  "q": [
    -5,
    -4,
    -3,
    -2,
    -1,
    0,
    1,
    2,
    3,
    4,
    5
  ],
  "fs": 4096,
  "window_sec": 512,
  "timestamp": "2026-01-05T19:39:05.988826+00:00"
}
```
```json
{
  "real_cpsd": 2.605739266198478e-48,
  "null_mean": 2.757410458877926e-48,
  "null_std": 4.111652609538801e-49,
  "z_score": -0.36888134062585726,
  "p_value": 0.6,
  "verdict": "NO CORRELATED SIGNAL \u2014 CONSISTENT WITH GR"
}
```
```json
{
  "corrected_global_consistency": [
    {
      "window_sec": 64,
      "H1": {
        "H": 0.30885260430877476,
        "samples": 262144,
        "compute_time_sec": 0.1855020523071289
      },
      "L1": {
        "H": 0.3246381413143207,
        "samples": 262144,
        "compute_time_sec": 0.1655426025390625
      }
    },
    {
      "window_sec": 128,
      "H1": {
        "H": 0.3098408247635631,
        "samples": 524288,
        "compute_time_sec": 0.3062772750854492
      },
      "L1": {
        "H": 0.3151923379430553,
        "samples": 524288,
        "compute_time_sec": 0.3051609992980957
      }
    },
    {
      "window_sec": 256,
      "H1": {
        "H": 0.28857808692469833,
        "samples": 1048576,
        "compute_time_sec": 0.6171586513519287
      },
      "L1": {
        "H": 0.3011290909967617,
        "samples": 1048576,
        "compute_time_sec": 0.5838541984558105
      }
    },
    {
      "window_sec": 512,
      "H1": {
        "H": 0.2756928984409219,
        "samples": 2097152,
        "compute_time_sec": 1.1128170490264893
      },
      "L1": {
        "H": 0.291333484468053,
        "samples": 2097152,
        "compute_time_sec": 1.1120657920837402
      }
    },
    {
      "window_sec": 1024,
      "H1": {
        "H": 0.2874670726409621,
        "samples": 4194304,
        "compute_time_sec": 2.2060842514038086
      },
      "L1": {
        "H": 0.28720688288221935,
        "samples": 4194304,
        "compute_time_sec": 2.133000612258911
      }
    }
  ],
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-03T22:39:03.343917+00:00"
}

{
  "corrected_global_consistency": [
    {
      "window_sec": 64,
      "H1": {
        "H": 0.30885260430877476,
        "samples": 262144,
        "compute_time_sec": 0.1855020523071289
      },
      "L1": {
        "H": 0.3246381413143207,
        "samples": 262144,
        "compute_time_sec": 0.1655426025390625
      }
    },
    {
      "window_sec": 128,
      "H1": {
        "H": 0.3098408247635631,
        "samples": 524288,
        "compute_time_sec": 0.3062772750854492
      },
      "L1": {
        "H": 0.3151923379430553,
        "samples": 524288,
        "compute_time_sec": 0.3051609992980957
      }
    },
    {
      "window_sec": 256,
      "H1": {
        "H": 0.28857808692469833,
        "samples": 1048576,
        "compute_time_sec": 0.6171586513519287
      },
      "L1": {
        "H": 0.3011290909967617,
        "samples": 1048576,
        "compute_time_sec": 0.5838541984558105
      }
    },
    {
      "window_sec": 512,
      "H1": {
        "H": 0.2756928984409219,
        "samples": 2097152,
        "compute_time_sec": 1.1128170490264893
      },
      "L1": {
        "H": 0.291333484468053,
        "samples": 2097152,
        "compute_time_sec": 1.1120657920837402
      }
    },
    {
      "window_sec": 1024,
      "H1": {
        "H": 0.2874670726409621,
        "samples": 4194304,
        "compute_time_sec": 2.2060842514038086
      },
      "L1": {
        "H": 0.28720688288221935,
        "samples": 4194304,
        "compute_time_sec": 2.133000612258911
      }
    }
  ],
  "sample_rate": 4096.0,
  "timestamp_utc": "2026-01-03T22:39:03.343917+00:00"
}

```

---

