"""
FIN Search in PTA Phase 75: Log-Poisson Multifractal Scaling (MF-DFA q=0)
Uses the true Formal Definition of FIN: H(MF-DFA_q=0[C(x,y)])
Author: Antigravity & User | Date: Feb 2026
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
from scipy.special import logsumexp
from tqdm import tqdm
from pathlib import Path
from itertools import combinations
from astropy.coordinates import SkyCoord
import astropy.units as u
import warnings
import json

warnings.filterwarnings('ignore')
np.random.seed(42)

DATA = Path("./nano15")
extract_dir = DATA/"residuals"
par_dir = DATA/"parfiles"

print("Loading Data...")
residuals = {}
for f in tqdm(list(extract_dir.rglob("*.res")), desc="Loading TXT"):
    psr = f.stem
    df = pd.read_csv(f, sep=r"\s+", comment="#", names=["MJD","res_us","err_us"], usecols=[0,1,2])
    df["res_s"]=df["res_us"]*1e-6
    residuals[psr] = df.sort_values("MJD")

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

def angular_sep(psr1, psr2):
    p1_base, p2_base = psr1.split('_')[0], psr2.split('_')[0]
    if p1_base in positions and p2_base in positions:
        return positions[p1_base].separation(positions[p2_base]).deg
    return None

def match_epochs(df1, df2, tol_days=30):
    t1, x1 = df1["MJD"].values, df1["res_s"].values
    t2, x2 = df2["MJD"].values, df2["res_s"].values
    matched_x, matched_y = [], []
    j = 0
    for i in range(len(t1)):
        while j < len(t2)-1 and t2[j] < t1[i] - tol_days: j += 1
        if abs(t2[j] - t1[i]) <= tol_days:
            matched_x.append(x1[i])
            matched_y.append(x2[j])
    if len(matched_x) < 50: return None, None
    return np.array(matched_x), np.array(matched_y)

def mf_dfa_q0_cross(x, y, min_scale=15):
    """
    Formal FIN Functional: H( MF-DFA_{q=0} [ C(x,y) ] )
    Where C(x,y) is the normalized relational product.
    """
    n = min(len(x), len(y))
    x, y = x[:n], y[:n]
    
    # 1. Non-energetic relational operator: Standardized Product
    sx = np.std(x)
    sy = np.std(y)
    if sx == 0 or sy == 0: return np.nan
    
    z = ((x - np.mean(x)) / sx) * ((y - np.mean(y)) / sy)
    
    # 2. Cumulative profile
    Z = np.cumsum(z - np.mean(z))
    if len(Z) < 4*min_scale: return np.nan

    scales = np.unique(np.logspace(np.log10(min_scale), np.log10(len(Z)//4), 15).astype(int))
    F_q0 = []
    
    for s in scales:
        k = len(Z)//s
        if k < 2: continue
        variances = []
        for i in range(k):
             zs = Z[i*s:(i+1)*s]
             t = np.arange(s)
             pz = np.polyfit(t, zs, 1)
             var = np.mean((zs - np.polyval(pz, t))**2)
             if var > 1e-15:  # Avoid log(0)
                 variances.append(var)
        
        if len(variances) < 2:
            F_q0.append(np.nan)
        else:
            # 3. Log-Poisson Multifractal limit (q=0): exp( mean( ln(var) ) / 2 )
            mean_ln_var = np.mean(np.log(variances))
            F_q0.append(np.exp(mean_ln_var / 2.0))

    F_q0 = np.array(F_q0)
    valid = ~np.isnan(F_q0) & (F_q0 > 0)
    if np.sum(valid) < 3: return np.nan
    
    slope, _, _, _, _ = linregress(np.log(scales[valid]), np.log(F_q0[valid]))
    return slope

psr_list = list(residuals.keys())
pairs = list(combinations(psr_list, 2))

print("\n--- Phase 75: Extracting Relational Functional (q=0) ---")
theta_vals = []
Hcross_real = []
pairs_used = []

for p1, p2 in tqdm(pairs, desc="MF-DFA(q=0) Real Pairs"):
    sep = angular_sep(p1, p2)
    if sep is None: continue
    x, y = match_epochs(residuals[p1], residuals[p2], tol_days=30)
    if x is None: continue
    
    Hxy_q0 = mf_dfa_q0_cross(x, y)
    if not np.isnan(Hxy_q0):
        theta_vals.append(sep)
        Hcross_real.append(Hxy_q0)
        pairs_used.append((p1, p2))

theta_vals = np.array(theta_vals)
Hcross_real = np.array(Hcross_real)

sigma_q0 = np.std(Hcross_real)
mean_q0 = np.mean(Hcross_real)
print(f"\n[q=0 Estimator Stats] Mean: {mean_q0:.4f} | Noise Floor (Sigma): {sigma_q0:.4f}")

# =====================================================
# Bayesian Analysis (Phase 71 Re-run)
# =====================================================
def hellings_downs(theta_deg):
    theta = np.radians(theta_deg)
    x = (1 - np.cos(theta)) / 2
    return 1.5*x*np.log(x + 1e-12) - 0.25*x + 0.5

def loglike_flat(theta, H, sigma, C):
    model = np.full_like(theta, C)
    return -0.5*np.sum(((H-model)/sigma)**2)

def loglike_fin(theta, H, sigma, A, C):
    model = A*hellings_downs(theta) + C
    return -0.5*np.sum(((H-model)/sigma)**2)

def logprior_fin(A, C):
    if -2 < A < 2 and -1 < C < 2: return 0.0
    return -np.inf

def logprior_flat(C):
    if -1 < C < 2: return 0.0
    return -np.inf

def evidence_flat(theta, H, sigma, N=20000):
    samples = np.random.uniform(-1, 2, N)
    logL = []
    for C in samples:
        lp = logprior_flat(C)
        if np.isfinite(lp):
            logL.append(loglike_flat(theta, H, sigma, C))
    return logsumexp(logL) - np.log(len(logL))

def evidence_fin(theta, H, sigma, N=40000):
    A = np.random.uniform(-2, 2, N)
    C = np.random.uniform(-1, 2, N)
    logL = []
    for a, c in zip(A, C):
        lp = logprior_fin(a, c)
        if np.isfinite(lp):
            logL.append(loglike_fin(theta, H, sigma, a, c))
    return logsumexp(logL) - np.log(len(logL))

print("\n--- Re-running Phase 71 (Bayesian Evidence) with q=0 ---")
logZ_flat = evidence_flat(theta_vals, Hcross_real, sigma_q0)
logZ_fin = evidence_fin(theta_vals, Hcross_real, sigma_q0)
logB_real = logZ_fin - logZ_flat

print(f"Log Evidence FLAT (q=0): {logZ_flat:.4f}")
print(f"Log Evidence FIN  (q=0): {logZ_fin:.4f}")
print(f"log Bayes Factor (FIN vs Flat) = {logB_real:.4f}")

# =====================================================
# FIN Injection Recovery Test (Phase 74 Re-run)
# =====================================================
print("\n--- Re-running Phase 74 (Injection Recovery) with q=0 ---")
A_inj = 0.35
def inject_fin(theta, H, sigma_noise):
    fin_signal = A_inj * hellings_downs(theta)
    noise = np.random.normal(0, sigma_noise, len(theta))
    return H + fin_signal + noise

H_injected = inject_fin(theta_vals, Hcross_real, sigma_q0)
sigma_inj = np.std(H_injected)

logZ_flat_inj = evidence_flat(theta_vals, H_injected, sigma_inj)
logZ_fin_inj = evidence_fin(theta_vals, H_injected, sigma_inj)

logB_inj = logZ_fin_inj - logZ_flat_inj

print(f"Injection Log Evidence FLAT: {logZ_flat_inj:.4f}")
print(f"Injection Log Evidence FIN:  {logZ_fin_inj:.4f}")
print(f"Injection log Bayes (FIN vs Flat) = {logB_inj:.4f}")

# Plotting the raw q=0 profile
plt.figure()
plt.scatter(theta_vals, Hcross_real, alpha=0.3, label="MF-DFA(q=0) Pairs")
# Fit Hellings-Downs just to visualize HD shape on this data
from scipy.optimize import curve_fit
def hd_fit(theta, A, C): return A*hellings_downs(theta) + C
try:
    popt, _ = curve_fit(hd_fit, theta_vals, Hcross_real)
    t_plot = np.linspace(0, 180, 100)
    plt.plot(t_plot, hd_fit(t_plot, *popt), 'r-', linewidth=2, label=f"HD Fit (A={popt[0]:.3f})")
except: pass
plt.axhline(mean_q0, color='k', linestyle='--', label=f"Mean={mean_q0:.3f}")
plt.xlabel("Angular Separation (degrees)")
plt.ylabel(r"Cross-Hurst $H_{cross}(q=0)$")
plt.title(f"FIN Phase 75: MF-DFA(q=0) Relational Functional\n$\sigma={sigma_q0:.4f}$, logB={logB_real:.2f}")
plt.legend()
plt.savefig(DATA/"phase75_mdfa_q0.png", dpi=200, bbox_inches='tight')
plt.close()

with open(DATA/"phase75_results.json", "w") as f:
    json.dump({
        "sigma_q0": sigma_q0,
        "logB_real": logB_real,
        "logB_inj": logB_inj
    }, f, indent=4)
