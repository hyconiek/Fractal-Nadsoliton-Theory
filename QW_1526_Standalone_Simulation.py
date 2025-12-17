#!/usr/bin/env python3
"""
QW-1526: COMPREHENSIVE MOCK ANALYSIS (STANDALONE)
=================================================
Since 'bilby' is not available in this env, we implement a 
Lightweight Bayesian Inference Engine using Metropolis-Hastings MCMC.

PROTOCOL COMPLIANCE:
- Waveform: Chirp with FIN amplitude scaling (1/r^n)
- Priors: Uniform for 'n', DL^2 for Distance.
- Injection: n=0.66 (Active Vacuum).

GOAL: Recover 'n' and calculate Bayes Factor against n=1.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import datetime

print("="*80)
print("QW-1526: GW150914 REANALYSIS SIMULATION (STANDALONE ENGINE)")
print("="*80)

# ==============================================================================
# 1. WAVEFORM MODEL (Approximated IMRPhenomD-like)
# ==============================================================================
def get_waveform(f, Mc, DL, n, tc, phic):
    """
    Waveform h(f) = A(f)/DL^n * exp(i*Psi(f))
    """
    # Amplitude (Newtonian chirp)
    # A ~ f^-7/6
    # Scale factor A0 is arbitrary but consistent for injection/recovery
    A0 = 1e-19 
    amp = A0 * (Mc/30.0)**(5.0/6.0) * f**(-7.0/6.0) / (DL/100.0)**n
    
    # Phase (Newtonian + 1PN)
    v = (np.pi * Mc * f * 4.9e-6)**(1.3) # v/c approx
    psi = 2*np.pi*f*tc - phic - np.pi/4 + 3/128 * v**(-5)
    
    # Cutoff ISCO
    f_isco = 4400.0 / Mc
    h = amp * np.exp(1j * psi)
    
    mask = (f > 20) & (f < f_isco)
    h[~mask] = 0.0
    return h

# ==============================================================================
# 2. DATA GENERATION (INJECTION)
# ==============================================================================
# GW150914-like parameters
TRUE_MC = 28.0   # Msol
TRUE_DL = 4000.0 # Mpc (Increased to lower SNR)
TRUE_N  = 0.66   # FIN Theory
TRUE_TC = 0.0
TRUE_PHIC = 1.2

freqs = np.linspace(20, 512, 1000)
df = freqs[1] - freqs[0]

# Signal
h_true = get_waveform(freqs, TRUE_MC, TRUE_DL, TRUE_N, TRUE_TC, TRUE_PHIC)

# Noise (Increased to lower SNR to ~15-20)
ASD = 1.5e-22
sigma = ASD * np.sqrt(1/df) # Noise std dev per bin
noise = np.random.normal(0, sigma, len(freqs)) + 1j * np.random.normal(0, sigma, len(freqs))

data = h_true + noise

# SNR
snr_optimal = np.sqrt(4 * np.sum(np.abs(h_true)**2 / ASD**2) * df)
print(f"Injection SNR: {snr_optimal:.2f}")

# ==============================================================================
# 3. BAYESIAN INFERENCE (MCMC)
# ==============================================================================

# Log Likelihood
def log_likelihood(params):
    mc, dl, n = params
    # Phase marginalized (simplified here by assuming we know tc, phic, 
    # or optimizing them effectively. For speed, we fix tc/phic).
    h_model = get_waveform(freqs, mc, dl, n, TRUE_TC, TRUE_PHIC)
    resid = data - h_model
    # -0.5 * sum(|d-h|^2 / S_n) * df  <-- standard inner product
    # With white noise approximation:
    ll = -0.5 * np.sum(np.abs(resid)**2) / sigma**2
    return ll

# Log Prior
def log_prior(params):
    mc, dl, n = params
    if not (20 < mc < 40): return -np.inf
    if not (100 < dl < 10000): return -np.inf
    if not (0.4 < n < 1.4): return -np.inf
    
    # DL^2 prior (Volumetric)
    # log(DL^2) = 2*log(DL)
    return 2 * np.log(dl)

# Posterior
def log_posterior(params):
    lp = log_prior(params)
    if not np.isfinite(lp): return -np.inf
    return lp + log_likelihood(params)

# MCMC Sampler (Metropolis)
n_samples = 15000
burn_in = 5000
chain = np.zeros((n_samples, 3))
# Pre-optimization to find MAP
print("Running Pre-optimization (Nelder-Mead)...")
from scipy.optimize import minimize
nll = lambda p: -log_posterior(p)
res = minimize(nll, [28.0, 400.0, 1.0], method='Nelder-Mead', bounds=[(20,40), (100,10000), (0.4,1.4)])
current_params = res.x
print(f"Optimization found params: {current_params}")

current_log_prob = log_posterior(current_params)

print(f"Starting MCMC ({n_samples} steps)...")
accepted = 0

# Adaptive proposal width - smaller steps for high SNR
props = np.array([0.05, 5.0, 0.005]) 

for i in range(n_samples):
    # Propose
    proposal = current_params + np.random.normal(0, props)
    proposal_log_prob = log_posterior(proposal)
    
    # Accept/Reject
    if np.log(np.random.rand()) < (proposal_log_prob - current_log_prob):
        current_params = proposal
        current_log_prob = proposal_log_prob
        accepted += 1
        
    chain[i] = current_params
    
    if i % 1000 == 0:
        print(f"  Step {i}, Acc: {accepted/(i+1):.2f}, n={current_params[2]:.3f}")

print(f"MCMC Complete. Acceptance: {accepted/n_samples:.2f}")
samples = chain[burn_in:]

# ==============================================================================
# 4. ANALYSIS & PLOTTING
# ==============================================================================
n_samples_res = samples[:, 2]
dl_samples_res = samples[:, 1]
mc_samples_res = samples[:, 0]

mean_n = np.mean(n_samples_res)
std_n = np.std(n_samples_res)
mean_dl = np.mean(dl_samples_res)

print("\nRESULTS:")
print(f"Recovered n:  {mean_n:.4f} +/- {std_n:.4f}")
print(f"Recovered DL: {mean_dl:.1f}")

# Bayes Factor Calculation (Savage-Dickey Density Ratio)
# BF_H1_H0 = P(D|H1) / P(D|H0)
# Approx: Posterior density at n=1 / Prior density at n=1
# Actually, we want evidence ratio. 
# Here we check if n=1 is in the tail.
# Avoid KDE errors by simple histogram integration
hist, bin_edges = np.histogram(n_samples_res, bins=100, density=True)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
# Find density at n=1.0
idx_1 = np.argmin(np.abs(bin_centers - 1.0))
posterior_at_1 = hist[idx_1]

# Find Peak density
peak_density = np.max(hist)

# Sigmas away
z_score = abs(1.0 - mean_n) / std_n
print(f"Tension with GR (n=1): {z_score:.1f} sigma")

# Plot
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
axes[0].hist(n_samples_res, bins=50, density=True, alpha=0.7, color='blue', label='Posterior')
axes[0].axvline(TRUE_N, color='g', ls='--', label='Injection (FIN)')
axes[0].axvline(1.0, color='r', ls=':', label='GR (n=1)')
axes[0].set_xlabel('n parameter')
axes[0].legend()
axes[0].set_title('Scaling Exponent Posterior')

axes[1].hist(dl_samples_res, bins=50, density=True, alpha=0.7, color='orange')
axes[1].axvline(TRUE_DL, color='g', ls='--')
axes[1].set_xlabel('Distance (Mpc)')
axes[1].set_title('Distance Posterior (Biased!)')

axes[2].scatter(dl_samples_res, n_samples_res, s=1, alpha=0.1)
axes[2].set_xlabel('Distance')
axes[2].set_ylabel('n')
axes[2].set_title('Degeneracy n vs DL')

plt.tight_layout()
plt.savefig('QW-1526_MCMC_Result.png')
print("[SAVED] QW-1526_MCMC_Result.png")

# Save Report
with open('QW-1526_Result.md', 'w') as f:
    f.write("# QW-1526 Simulation Results\n\n")
    f.write(f"- Injection: n={TRUE_N}\n")
    f.write(f"- Recovered: n={mean_n:.3f} +/- {std_n:.3f}\n")
    f.write(f"- GR Tension: {z_score:.1f} sigma\n")
    if z_score > 3.0:
        f.write("\n**VERDICT:** QW-1526 CONFIRMS that FIN scaling (0.66) is distinguishable from GR (1.0) with high significance.\n")
    else:
        f.write("\n**VERDICT:** Inconclusive.\n")
