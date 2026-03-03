"""
FIN Search in PTA Phase 70: Referee-Grade FIN Hellings-Downs Analog Test
Fixing Interpolation Bias, DCCA sign bias, Shared-Phase Null, and HD curve fitting.
Author: Antigravity & User | Date: Feb 2026
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
from scipy.optimize import curve_fit
from tqdm import tqdm
import tarfile, os
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

# 1. Load Data
residuals = {}
for f in tqdm(list(extract_dir.rglob("*.res")), desc="Loading TXT"):
    psr = f.stem
    df = pd.read_csv(f, sep=r"\s+", comment="#", names=["MJD","res_us","err_us"], usecols=[0,1,2])
    df["res_s"]=df["res_us"]*1e-6
    residuals[psr] = df.sort_values("MJD")

# 2. Load Pulsar Positions
def load_pulsar_positions(parfile_dir):
    pos = {}
    for par in Path(parfile_dir).rglob("*.par"):
        name = par.stem
        base_name = name.split('_')[0]
        ra, dec, elong, elat = None, None, None, None
        with open(par) as f:
            for line in f:
                if line.startswith("RAJ"): ra = line.split()[1]
                if line.startswith("DECJ"): dec = line.split()[1]
                if line.startswith("ELONG"): elong = line.split()[1]
                if line.startswith("ELAT"): elat = line.split()[1]
        try:
            if ra and dec:
                pos[base_name] = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
            elif elong and elat:
                pos[base_name] = SkyCoord(lon=float(elong)*u.deg, lat=float(elat)*u.deg, frame='barycentrictrueecliptic').transform_to('icrs')
        except: pass
    return pos

positions = load_pulsar_positions(par_dir)
print(f"Loaded {len(positions)} pulsar positions.")

def angular_sep(psr1, psr2):
    p1_base, p2_base = psr1.split('_')[0], psr2.split('_')[0]
    if p1_base in positions and p2_base in positions:
        return positions[p1_base].separation(positions[p2_base]).deg
    return None

# 3. Patch 1: Epoch Matching (No Interpolation)
def match_epochs(df1, df2, tol_days=30):
    t1 = df1["MJD"].values
    x1 = df1["res_s"].values
    t2 = df2["MJD"].values
    x2 = df2["res_s"].values

    matched_x = []
    matched_y = []
    j = 0
    for i in range(len(t1)):
        while j < len(t2)-1 and t2[j] < t1[i] - tol_days:
            j += 1
        if abs(t2[j] - t1[i]) <= tol_days:
            matched_x.append(x1[i])
            matched_y.append(x2[j])

    if len(matched_x) < 50:
        return None, None
    return np.array(matched_x), np.array(matched_y)

# 4. Patch 2: Correct DCCA Estimator
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
             
        mean_cov = np.mean(covs)
        if mean_cov != 0:
            F.append(np.sign(mean_cov) * np.sqrt(abs(mean_cov)))
        else:
            F.append(0.0)

    # Need positive values for log-log regression. 
    # If the general trend is negative, H is negative (anti-persistence).
    # Since H > 0 usually, F should be mostly positive. 
    # Let's take abs(F) for regression, but track the overall sign.
    F = np.array(F)
    valid = F != 0
    if np.sum(valid) < 3: return np.nan
    
    slope, _, _, _, _ = linregress(np.log(scales[valid]), np.log(np.abs(F[valid])))
    # Assign overall sign based on the mean covariance sign at large scales
    overall_sign = np.sign(np.mean(F))
    return slope * overall_sign

# 5. Patch 3: Shared-Phase Surrogate
def shared_phase_randomize(x, y):
    fx = np.fft.rfft(x)
    fy = np.fft.rfft(y)
    phase = np.exp(1j*np.random.uniform(0, 2*np.pi, len(fx)))
    xr = np.fft.irfft(np.abs(fx)*phase, n=len(x))
    yr = np.fft.irfft(np.abs(fy)*phase, n=len(y))
    return xr, yr

# 6. Execute Angular Cross-DFA
psr_list = list(residuals.keys())
pairs = list(combinations(psr_list, 2))

# Execute on all valid pairs instead of full random subselection
theta_vals = []
Hcross_real = []
Hcross_phase = []

print("Running Unbiased Cross-DFA on all possible pairs...")
for p1, p2 in tqdm(pairs, desc="Cross DFA"):
    sep = angular_sep(p1, p2)
    if sep is None: continue

    x, y = match_epochs(residuals[p1], residuals[p2], tol_days=30)
    if x is None: continue
    
    Hxy = cross_dfa(x, y)
    
    # Null Hypothesis: Shared Phase Red Noise
    xr, yr = shared_phase_randomize(x, y)
    Hxy_p = cross_dfa(xr, yr)

    if not np.isnan(Hxy) and not np.isnan(Hxy_p):
        theta_vals.append(sep)
        Hcross_real.append(Hxy)
        Hcross_phase.append(Hxy_p)

theta_vals = np.array(theta_vals)
Hcross_real = np.array(Hcross_real)
Hcross_phase = np.array(Hcross_phase)

# 7. Patch 5: Balanced Angle Sampling and Binning
def balanced_angle_sampling(theta, H, bins=10):
    edges = np.linspace(0,180,bins+1)
    idx = np.digitize(theta, edges)
    out_t, out_H = [], []
    # Find min count across bins that have at least some minimal threshold (e.g. 3)
    counts = [np.sum(idx==i) for i in range(1,bins+1)]
    valid_counts = [c for c in counts if c >= 3]
    if not valid_counts: return np.array([]), np.array([])
    nmin = min(valid_counts)

    for i in range(1,bins+1):
        mask = np.where(idx==i)[0]
        if len(mask) >= nmin:
            choose = np.random.choice(mask, nmin, replace=False)
            out_t.extend(theta[choose])
            out_H.extend(H[choose])
    return np.array(out_t), np.array(out_H)

# Sample to avoid bin bias
t_bal_real, H_bal_real = balanced_angle_sampling(theta_vals, Hcross_real, bins=12)
t_bal_phase, H_bal_phase = balanced_angle_sampling(theta_vals, Hcross_phase, bins=12)

# Binning for plotting and fitting
bins = np.linspace(0,180,12)
digitized = np.digitize(t_bal_real, bins)

theta_mean, H_mean_real, H_mean_phase = [], [], []

for b in range(1, len(bins)):
    mask = digitized == b
    if np.sum(mask) >= 3:
        theta_mean.append(np.mean(t_bal_real[mask]))
        H_mean_real.append(np.mean(H_bal_real[mask]))
        H_mean_phase.append(np.mean(H_bal_phase[mask]))

# 8. Patch 4: HD Functional Fit
def hd_curve(theta_deg, A, C):
    theta = np.radians(theta_deg)
    x = (1 - np.cos(theta)) / 2
    return A*(1.5*x*np.log(x+1e-10) - 0.25*x + 0.5) + C

try:
    popt_hd, _ = curve_fit(hd_curve, theta_mean, H_mean_real, p0=[0.1, np.mean(H_mean_real)])
    model_hd_real = hd_curve(np.array(theta_mean), *popt_hd)
    flat_real = np.full_like(theta_mean, np.mean(H_mean_real))
    rss_hd_real = np.sum((H_mean_real - model_hd_real)**2)
    rss_flat_real = np.sum((H_mean_real - flat_real)**2)
    delta_rss_real = rss_flat_real - rss_hd_real
except:
    delta_rss_real = 0

try:
    popt_phase, _ = curve_fit(hd_curve, theta_mean, H_mean_phase, p0=[0.1, np.mean(H_mean_phase)])
    model_hd_phase = hd_curve(np.array(theta_mean), *popt_phase)
    flat_phase = np.full_like(theta_mean, np.mean(H_mean_phase))
    rss_hd_phase = np.sum((H_mean_phase - model_hd_phase)**2)
    rss_flat_phase = np.sum((H_mean_phase - flat_phase)**2)
    delta_rss_phase = rss_flat_phase - rss_hd_phase
except:
    delta_rss_phase = 0

print(f"ΔRSS REAL (Flat - HD): {delta_rss_real:.6f}")
print(f"ΔRSS PHASE (Flat - HD): {delta_rss_phase:.6f}")

# 9. Plotting
plt.figure()
plt.scatter(t_bal_real, H_bal_real, alpha=0.15, label="Balanced Pairs (Real)", color='blue')
plt.plot(theta_mean, H_mean_real, "o-", lw=3, label="FIN Cross-H (Real)", color='blue')

if 'model_hd_real' in locals():
    smooth_theta = np.linspace(0, 180, 100)
    plt.plot(smooth_theta, hd_curve(smooth_theta, *popt_hd), "b--", label="HD Curve Fit (Real)")

plt.plot(theta_mean, H_mean_phase, "s-", lw=2, label="FIN Cross-H (Shared Phase)", color='orange')

plt.axhline(np.mean(H_mean_real), color='black', linestyle=':', label="Flat Mean")

plt.xlabel("Angular separation (deg)")
plt.ylabel("Cross Hurst exponent $H_{cross}$")
plt.title(f"FIN Phase 70 HD-Test (ΔRSS Real: {delta_rss_real:.4f})")
plt.legend()
plt.tight_layout()
plt.savefig(DATA/"phase70_hd_test.png")

with open(DATA/"phase70_results.json", "w") as f:
    import json
    json.dump({
        "delta_rss_real": delta_rss_real,
        "delta_rss_phase": delta_rss_phase,
        "mean_real": np.mean(H_mean_real),
        "mean_phase": np.mean(H_mean_phase),
        "hd_params": popt_hd.tolist() if 'popt_hd' in locals() else []
    }, f, indent=4)
