"""
FIN Search in PTA Phases 72-74: Final Falsification Framework
Includes PPC, Sky Scramble, and Injection Recovery
Author: Antigravity & User | Date: Feb 2026
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
from scipy.special import logsumexp
from tqdm import tqdm
import tarfile, os
from pathlib import Path
from itertools import combinations
from astropy.coordinates import SkyCoord
import astropy.units as u
import warnings

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
        if mean_cov != 0: F.append(np.sign(mean_cov) * np.sqrt(abs(mean_cov)))
        else: F.append(0.0)

    F = np.array(F)
    valid = F != 0
    if np.sum(valid) < 3: return np.nan
    slope, _, _, _, _ = linregress(np.log(scales[valid]), np.log(np.abs(F[valid])))
    return slope * np.sign(np.mean(F))

psr_list = list(residuals.keys())
pairs = list(combinations(psr_list, 2))

theta_vals = []
Hcross_real = []
pairs_used = []

for p1, p2 in tqdm(pairs, desc="Cross DFA Real Pairs"):
    sep = angular_sep(p1, p2)
    if sep is None: continue
    x, y = match_epochs(residuals[p1], residuals[p2], tol_days=30)
    if x is None: continue
    
    Hxy = cross_dfa(x, y)
    if not np.isnan(Hxy):
        theta_vals.append(sep)
        Hcross_real.append(Hxy)
        pairs_used.append((p1, p2))

theta_vals = np.array(theta_vals)
Hcross_real = np.array(Hcross_real)

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

sigma = np.std(Hcross_real)

# =====================================================
# Phase 72 — Posterior Predictive Check
# =====================================================
print("\n[Phase 72] Running Posterior Predictive Check (PPC)...")

def rss_stat(theta, H, A, C):
    model = A*hellings_downs(theta) + C
    return np.sum((H - model)**2)

Nrep = 2000
rss_rep = []

A_samples = np.random.uniform(-2,2,40000)
C_samples = np.random.uniform(-1,2,40000)

logL_all = np.array([
    loglike_fin(theta_vals, Hcross_real, sigma, a, c)
    for a,c in zip(A_samples, C_samples)
])

weights = np.exp(logL_all - logsumexp(logL_all))

idx = np.random.choice(
    np.arange(len(weights)),
    size=Nrep,
    p=weights
)

for i in tqdm(idx, desc="PPC simulations"):
    A = A_samples[i]
    C = C_samples[i]
    H_sim = A*hellings_downs(theta_vals) + C + np.random.normal(0, sigma, len(theta_vals))
    rss_rep.append(rss_stat(theta_vals, H_sim, A, C))

rss_rep = np.array(rss_rep)
rss_obs = rss_stat(
    theta_vals,
    Hcross_real,
    np.mean(A_samples[idx]),
    np.mean(C_samples[idx])
)

p_ppc = np.mean(rss_rep >= rss_obs)
print(f"PPC p-value = {p_ppc:.4f}")

# =====================================================
# Phase 73 — Sky Scramble Null Test
# =====================================================
print("\n[Phase 73] Running Sky Scramble Null Test...")

def compute_rss(theta, H):
    A = 1.0 # Or use fixed mean amplitude instead
    C = np.mean(H)
    model = A*hellings_downs(theta) + C
    return np.sum((H-model)**2)

rss_real_scramble = compute_rss(theta_vals, Hcross_real)

Nscramble = 500
rss_scrambled = []

psr_names = list(positions.keys())
coords = list(positions.values())

for _ in tqdm(range(Nscramble), desc="Sky scrambles"):
    shuffled = np.random.permutation(coords)
    fake_pos = dict(zip(psr_names, shuffled))
    
    theta_fake = []
    H_fake = []
    for (p1,p2), H in zip(pairs_used, Hcross_real):
        b1 = p1.split('_')[0]
        b2 = p2.split('_')[0]
        if b1 in fake_pos and b2 in fake_pos:
            sep = fake_pos[b1].separation(fake_pos[b2]).deg
            theta_fake.append(sep)
            H_fake.append(H)
    
    if len(theta_fake) > 0:
        rss_scrambled.append(compute_rss(np.array(theta_fake), np.array(H_fake)))

rss_scrambled = np.array(rss_scrambled)
p_null = np.mean(rss_scrambled <= rss_real_scramble)

print(f"REAL RSS = {rss_real_scramble:.4f}")
print(f"Null probability = {p_null:.4f}")

# =====================================================
# Phase 74 — FIN Injection Recovery Test
# =====================================================
print("\n[Phase 74] Running FIN Injection Test...")

A_inj = 0.35   # injected FIN strength

def inject_fin(theta, H, sigma_noise):
    fin_signal = A_inj * hellings_downs(theta)
    noise = np.random.normal(0, sigma_noise, len(theta))
    return H + fin_signal + noise

H_injected = inject_fin(theta_vals, Hcross_real, sigma)

# recompute evidences for injection
sigma_inj = np.std(H_injected)

logZ_flat_inj = evidence_flat(theta_vals, H_injected, sigma_inj)
logZ_fin_inj = evidence_fin(theta_vals, H_injected, sigma_inj)

logB_inj = logZ_fin_inj - logZ_flat_inj

print("\nInjection Results:")
print(f"logZ_flat = {logZ_flat_inj:.4f}")
print(f"logZ_fin  = {logZ_fin_inj:.4f}")
print(f"log Bayes (FIN vs Flat) = {logB_inj:.4f}")

with open(DATA/"phases72_74_results.txt", "w") as f:
    f.write(f"Phase 72 PPC p-value = {p_ppc:.4f}\n")
    f.write(f"Phase 73 REAL RSS = {rss_real_scramble:.4f}\n")
    f.write(f"Phase 73 Null p-value = {p_null:.4f}\n")
    f.write(f"Phase 74 Injection logB = {logB_inj:.4f}\n")
