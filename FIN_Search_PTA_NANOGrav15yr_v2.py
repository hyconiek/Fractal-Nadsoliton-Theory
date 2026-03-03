"""
FIN Search in PTA: DFA + Cross-Hurst on NANOGrav 15yr Residuals
Phase 1 of looking for H_FIN ≈ 0.23 global temporal memory outside LIGO.
Optimized for Kaggle Execution.

Author: Antigravity & User | Date: Feb 2026
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
from tqdm import tqdm
import tarfile, requests, os
from pathlib import Path

sns.set_style("darkgrid")
plt.rcParams['figure.figsize']=(12,6)

# ===================== 2. SMART DOWNLOAD =====================
DATA = Path("./nano15")
DATA.mkdir(exist_ok=True)

url = "https://zenodo.org/records/16051178/files/NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz?download=1"
tar_path = DATA/"nano15.tar.gz"

if not tar_path.exists():
    print("Downloading NANOGrav 15-yr Data...")
    r = requests.get(url, stream=True)
    with open(tar_path,"wb") as f:
        for chunk in tqdm(r.iter_content(1024*1024)):
            f.write(chunk)

print("Download OK")

# ===================== 3. EXTRACT ONLY RESIDUALS =====================
extract_dir = DATA/"residuals"
extract_dir.mkdir(exist_ok=True)

if not list(extract_dir.glob("*.res")):
    print("Extracting residuals...")
    with tarfile.open(tar_path,"r:gz") as tar:
        members = [
            m for m in tar.getmembers()
            if "residuals/" in m.name
            and m.name.endswith(".res")
            and "avg.res" in m.name  # Using avg epoch residuals to reduce noise
        ]
        for m in tqdm(members):
            tar.extract(m, extract_dir)

print("Residuals extracted")

# ===================== 4. LOAD DATA =====================
residuals = {}
files = list(extract_dir.rglob("*.res"))

def format_df(f):
    df = pd.read_csv(
        f,
        sep=r"\s+",
        comment="#",
        names=["MJD","res_us","err_us"],
        usecols=[0,1,2]
    )
    df["res_s"]=df["res_us"]*1e-6
    return df.sort_values("MJD")

for f in tqdm(files, desc="Loading txt files"):
    psr = f.stem
    residuals[psr] = format_df(f)

print("Loaded pulsars:",len(residuals))

# ===================== INTERPOLATION FOR UNIFORM TIME SERIES =====================
# For strict DFA, we need uniform sampling. We resample to a common barycentric time grid for Cross-DFA.
def get_uniform_grid_and_interp(r1_df, r2_df=None, N_points=1000):
    t1 = r1_df['MJD'].values
    x1 = r1_df['res_s'].values
    
    if r2_df is not None:
        t2 = r2_df['MJD'].values
        x2 = r2_df['res_s'].values
        common_t = np.linspace(max(t1.min(), t2.min()), min(t1.max(), t2.max()), N_points)
        r1_interp = np.interp(common_t, t1, x1)
        r2_interp = np.interp(common_t, t2, x2)
        return r1_interp, r2_interp
    else:
        common_t = np.linspace(t1.min(), t1.max(), N_points)
        return np.interp(common_t, t1, x1)

# ===================== 5. DFA (POPRAWNA IMPLEMENTACJA) =====================
def dfa(signal, min_scale=20):
    x = signal - np.mean(signal)
    y = np.cumsum(x)
    
    # Needs to be at least 4*min_scale long
    if len(y) < 4*min_scale: return np.nan
    
    scales = np.unique(np.logspace(np.log10(min_scale), np.log10(len(y)//4), 20).astype(int))
    F = []

    for s in scales:
        n = len(y)//s
        if n < 2: continue

        rms = []
        for i in range(n):
            seg = y[i*s:(i+1)*s]
            t = np.arange(s)
            p = np.polyfit(t,seg,1)
            rms.append(np.sqrt(np.mean((seg - np.polyval(p,t))**2)))

        if len(rms) > 0 and np.mean(rms) > 0:
            F.append(np.exp(np.mean(np.log(np.array(rms) + 1e-300))))

    if len(F) < 3: return np.nan
    slope, _, r, _, _ = linregress(np.log(scales[:len(F)]), np.log(F))
    return slope

# ===================== 6. CROSS-DFA =====================
def cross_dfa(x, y, min_scale=20):
    n = min(len(x), len(y))
    x = x[:n] - np.mean(x[:n])
    y = y[:n] - np.mean(y[:n])
    
    z = np.cumsum(x * y) # THIS IS STILL CUMULANT OF PRODUCT.
    # User proposed the correct local covariance approach:
    # F^2_xy(s) = 1/s sum (X_s - ~X_s)(Y_s - ~Y_s)
    
    X = np.cumsum(x)
    Y = np.cumsum(y)
    
    if len(X) < 4*min_scale: return np.nan

    scales = np.unique(np.logspace(np.log10(min_scale), np.log10(n//4), 20).astype(int))
    F = []

    for s in scales:
        k = len(X)//s
        if k < 2: continue

        covs = []
        for i in range(k):
            xs = X[i*s:(i+1)*s]
            ys = Y[i*s:(i+1)*s]
            t = np.arange(s)
            px = np.polyfit(t, xs, 1)
            py = np.polyfit(t, ys, 1)
            
            dx = xs - np.polyval(px, t)
            dy = ys - np.polyval(py, t)
            
            # Local covariance
            val = np.mean(dx * dy)
            covs.append(val)

        # F_xy = sqrt(mean(|cov|)) for stable log-log plotting
        mean_cov = np.mean(np.abs(covs))
        if mean_cov > 0:
            F.append(np.sqrt(mean_cov))

    if len(F) < 3: return np.nan
    slope, _, _, _, _ = linregress(np.log(scales[:len(F)]), np.log(F))
    return slope

# ===================== 7. KLUCZOWE: SURROGATE TESTY =====================
def phase_randomize(ts):
    ft = np.fft.rfft(ts)
    amp = np.abs(ft)
    random_phase = np.exp(1j * np.random.uniform(0, 2*np.pi, len(ft)))
    new_ft = amp * random_phase
    return np.fft.irfft(new_ft, n=len(ts))

def time_shuffle(ts):
    return np.random.permutation(ts)

# ===================== 8. SINGLE H TEST =====================
psr_list = list(residuals.keys())[:25] # Selecting first 25 for quick execution

H_real = {}
H_phase = {}
H_shuffle = {}

print("Calculating Single H (Real, Phase-Rand, Shuffle)...")
for psr in tqdm(psr_list, desc="Single DFA"):
    # Resample to uniform grid first to make DFA valid
    ts_uniform = get_uniform_grid_and_interp(residuals[psr], N_points=1000)
    
    H_r = dfa(ts_uniform)
    H_p = dfa(phase_randomize(ts_uniform))
    H_s = dfa(time_shuffle(ts_uniform))
    
    if not np.isnan(H_r):
        H_real[psr] = H_r
        H_phase[psr] = H_p
        H_shuffle[psr] = H_s

# ===================== 9. FIN TEST (WIZUALIZACJA) =====================
if H_real:
    plt.figure()
    sns.kdeplot(list(H_real.values()), label="REAL")
    sns.kdeplot(list(H_phase.values()), label="Phase rand")
    sns.kdeplot(list(H_shuffle.values()), label="Shuffle")

    plt.axvline(0.23, ls="--", color="red", label="FIN prediction")

    plt.legend()
    plt.title("Does H≈0.23 survive PSD destruction? (Single Pulsar)")
    plt.savefig(DATA/"fin_pta_single_h.png")
    
    print("\n--- TEST SURROGATE (SINGLE PULSAR) ---")
    print(f"REAL:    {np.mean(list(H_real.values())):.4f}")
    print(f"PHASE:   {np.mean(list(H_phase.values())):.4f}")
    print(f"SHUFFLE: {np.mean(list(H_shuffle.values())):.4f}")

    print("\nNext step for full physical check: Cross-DFA H(theta) Hellings-Downs verification!")
