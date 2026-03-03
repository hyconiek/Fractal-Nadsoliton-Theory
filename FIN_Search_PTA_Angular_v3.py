"""
FIN Search in PTA Phase 69: FIN Hellings-Downs Analog Test
Testing if Fractal Memory (Cross-H) depends on Angular Separation.
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
from itertools import combinations
from astropy.coordinates import SkyCoord
import astropy.units as u
import warnings

warnings.filterwarnings('ignore')
sns.set_style("darkgrid")
plt.rcParams['figure.figsize']=(12,6)

DATA = Path("./nano15")
tar_path = DATA/"nano15.tar.gz"
extract_dir = DATA/"residuals"
par_dir = DATA/"parfiles"

# 1. Extract .res and .par files if needed
if not list(extract_dir.glob("*.res")) or not list(par_dir.glob("*.par")):
    print("Extracting residuals and par files...")
    extract_dir.mkdir(exist_ok=True)
    par_dir.mkdir(exist_ok=True)
    with tarfile.open(tar_path,"r:gz") as tar:
        for m in tqdm(tar.getmembers(), desc="Scanning tarball"):
            if "residuals/" in m.name and m.name.endswith(".res") and "avg.res" in m.name:
                tar.extract(m, extract_dir)
            if "par/" in m.name and m.name.endswith(".par"):
                tar.extract(m, par_dir)

# 2. Load Data
residuals = {}
for f in tqdm(list(extract_dir.rglob("*.res")), desc="Loading TXT"):
    psr = f.stem
    df = pd.read_csv(f, sep=r"\s+", comment="#", names=["MJD","res_us","err_us"], usecols=[0,1,2])
    df["res_s"]=df["res_us"]*1e-6
    residuals[psr] = df.sort_values("MJD")

# 3. Load Pulsar Positions
def load_pulsar_positions(parfile_dir):
    pos = {}
    for par in Path(parfile_dir).rglob("*.par"):
        name = par.stem
        # Extract base name without extra suffix if present
        base_name = name.split('_')[0]
        ra, dec, elong, elat = None, None, None, None
        with open(par) as f:
            for line in f:
                if line.startswith("RAJ"):
                    ra = line.split()[1]
                if line.startswith("DECJ"):
                    dec = line.split()[1]
                if line.startswith("ELONG"):
                    elong = line.split()[1]
                if line.startswith("ELAT"):
                    elat = line.split()[1]
        
        try:
            if ra and dec:
                pos[base_name] = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
            elif elong and elat:
                pos[base_name] = SkyCoord(lon=float(elong)*u.deg, lat=float(elat)*u.deg, frame='barycentrictrueecliptic').transform_to('icrs')
        except:
            pass
    return pos

positions = load_pulsar_positions(par_dir)
print(f"Loaded {len(positions)} pulsar positions.")

def angular_sep(psr1, psr2):
    p1_base = psr1.split('_')[0]
    p2_base = psr2.split('_')[0]
    if p1_base in positions and p2_base in positions:
        return positions[p1_base].separation(positions[p2_base]).deg
    return None

# 4. Interpolation & DFA Functions
def get_uniform_grid_and_interp(r1_df, r2_df=None, N_points=500):
    t1 = r1_df['MJD'].values
    x1 = r1_df['res_s'].values
    if r2_df is not None:
        t2 = r2_df['MJD'].values
        x2 = r2_df['res_s'].values
        common_t = np.linspace(max(t1.min(), t2.min()), min(t1.max(), t2.max()), N_points)
        return np.interp(common_t, t1, x1), np.interp(common_t, t2, x2)
    common_t = np.linspace(t1.min(), t1.max(), N_points)
    return np.interp(common_t, t1, x1)

def cross_dfa(x, y, min_scale=15):
    n = min(len(x), len(y))
    x, y = x[:n] - np.mean(x[:n]), y[:n] - np.mean(y[:n])
    X, Y = np.cumsum(x), np.cumsum(y)
    
    if len(X) < 4*min_scale: return np.nan

    scales = np.unique(np.logspace(np.log10(min_scale), np.log10(n//4), 15).astype(int))
    F = []
    for s in scales:
        k = len(X)//s
        if k < 2: continue
        covs = []
        for i in range(k):
             xs, ys = X[i*s:(i+1)*s], Y[i*s:(i+1)*s]
             t = np.arange(s)
             px, py = np.polyfit(t, xs, 1), np.polyfit(t, ys, 1)
             covs.append(np.mean((xs - np.polyval(px, t)) * (ys - np.polyval(py, t))))
        mean_cov = np.mean(np.abs(covs))
        if mean_cov > 0: F.append(np.sqrt(mean_cov))

    if len(F) < 3: return np.nan
    return linregress(np.log(scales[:len(F)]), np.log(F))[0]

def phase_randomize(ts):
    ft = np.fft.rfft(ts)
    new_ft = np.abs(ft) * np.exp(1j * np.random.uniform(0, 2*np.pi, len(ft)))
    return np.fft.irfft(new_ft, n=len(ts))

# 5. Execute Angular Cross-DFA
psr_list = list(residuals.keys())
pairs = list(combinations(psr_list, 2))

# Subsample for faster execution, e.g. 500 random pairs
np.random.seed(42)
if len(pairs) > 500:
    indices = np.random.choice(len(pairs), 500, replace=False)
    pairs = [pairs[i] for i in indices]

theta_vals = []
Hcross_real = []
Hcross_phase = []

for p1, p2 in tqdm(pairs, desc="Cross DFA pairs"):
    sep = angular_sep(p1, p2)
    if sep is None: continue

    x, y = get_uniform_grid_and_interp(residuals[p1], residuals[p2], N_points=500)
    
    Hxy = cross_dfa(x, y)
    Hxy_p = cross_dfa(phase_randomize(x), phase_randomize(y))

    if not np.isnan(Hxy) and not np.isnan(Hxy_p):
        theta_vals.append(sep)
        Hcross_real.append(Hxy)
        Hcross_phase.append(Hxy_p)

theta_vals = np.array(theta_vals)
Hcross_real = np.array(Hcross_real)
Hcross_phase = np.array(Hcross_phase)

# 6. Plotting
bins = np.linspace(0, 180, 12)
digitized = np.digitize(theta_vals, bins)

theta_mean, H_mean_real, H_mean_phase = [], [], []

for b in range(1, len(bins)):
    mask = digitized == b
    if np.sum(mask) > 3:
        theta_mean.append(np.mean(theta_vals[mask]))
        H_mean_real.append(np.mean(Hcross_real[mask]))
        H_mean_phase.append(np.mean(Hcross_phase[mask]))

plt.figure(figsize=(10,6))
plt.scatter(theta_vals, Hcross_real, alpha=0.15, color='blue')
plt.plot(theta_mean, H_mean_real, "o-", lw=3, color='blue', label="Real Residuals (FIN)")

plt.plot(theta_mean, H_mean_phase, "s--", lw=3, color='orange', label="Phase Rand (Red Noise Null)")

plt.axhline(0.5, color='black', linestyle=':', label="White Noise Baseline")

plt.xlabel("Angular separation (degrees)")
plt.ylabel("Cross Hurst exponent $H_{cross}$")
plt.title("FIN Hellings-Downs-like Test: $H_{cross}(\\theta)$")
plt.legend()
plt.tight_layout()
plt.savefig("phase69_angular_fin_test.png")
print("Saved phase69_angular_fin_test.png")

# Calculate flat line correlation (Pearson r) to check constraint
real_r = linregress(theta_mean, H_mean_real)[2]**2
phase_r = linregress(theta_mean, H_mean_phase)[2]**2

with open("phase69_angular_fin_results.json", "w") as f:
    import json
    json.dump({
        "R2_real_vs_angle": real_r,
        "R2_phase_vs_angle": phase_r,
        "real_mean": np.mean(H_mean_real),
        "phase_mean": np.mean(H_mean_phase)
    }, f, indent=4)
