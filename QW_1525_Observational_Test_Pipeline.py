#!/usr/bin/env python3
"""
QW-1525: OBSERVATIONAL TEST PIPELINE (PROOF OF CONCEPT)
=========================================================
Demonstration of recovering the amplitude scaling scaling exponent 'n'
from Gravitational Wave data.

We test the model: h ~ 1 / r^n

If FIN Theory is correct, we should recover n ~ 0.66.
If GR is correct, we should recover n ~ 1.00.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import datetime

# ==============================================================================
# 1. WAVEFORM MODEL (frequency domain)
# ==============================================================================
# Simplified TaylorF2-like chirp (Newtonian amplitude, PN phase)

def get_chirp_waveform(f, Mc, DL, n_scaling, tc, phic):
    """
    Generate frequency domain waveform h(f).
    
    Arguments:
    - Mc: Chirp mass [Msol]
    - DL: Luminosity distance [Mpc]
    - n_scaling: Amplitude scaling exponent (1.0 = GR, 0.66 = FIN)
    """
    # Constants
    G = 4.30e-9  # Mpc * (km/s)^2 / Msol (approx, strictly need units)
    # Using natural units for simplicity in PoC, scaled to sensible Hz/Strain
    
    # Amplitude A(f) ~ f^(-7/6) * Mc^(5/6) / DL^n
    # Note: Standard GR has 1/DL. We use 1/DL^n.
    
    # Arbitrary prefactor for PoC scale
    A0 = 1e-21 
    
    # Frequency cutoff (ISCO)
    f_isco = 4400.0 / Mc # Approx for solar mass binaries
    
    amp = A0 * (Mc/30.0)**(5.0/6.0) * f**(-7.0/6.0) / (DL/100.0)**n_scaling
    
    # Phase Psi(f) (Newtonian)
    # Psi(f) = 2pi*f*tc - phic - pi/4 + 3/128 * (pi*Mc*f)^(-5/3)
    # Note: Mc in phase is well constrained by f_dot.
    
    psi = 2 * np.pi * f * tc - phic - np.pi/4 + (f/100.0)**(-5.0/3.0) 
    # (Simplified phase term for stability in PoC optimization)
    
    h = amp * np.exp(1j * psi)
    
    # Cutoff
    mask = (f > 20.0) & (f < f_isco)
    h[~mask] = 0.0
    
    return h

# ==============================================================================
# 2. INJECTION (SIMULATING NATURE)
# ==============================================================================

print("="*80)
print("QW-1525: GW OBSERVATIONAL TEST PROOF-OF-CONCEPT")
print("="*80)

# True Parameters (Nature = FIN Theory assumption for this test)
TRUE_N = 0.66      # Active Vacuum
TRUE_DL = 500.0    # Mpc
TRUE_MC = 30.0     # Solar masses
TRUE_TC = 0.0
TRUE_PHIC = 0.0

print(f"Hypothesis to Test: n_scaling")
print(f"Injection (Nature): n = {TRUE_N} (FIN Theory)")
print(f"Distance:           D = {TRUE_DL} Mpc")
print(f"Chirp Mass:         M = {TRUE_MC} Msol")

# Frequency grid (LIGO-like)
freqs = np.linspace(20, 1024, 1000)

# Generate Signal
h_true = get_chirp_waveform(freqs, TRUE_MC, TRUE_DL, TRUE_N, TRUE_TC, TRUE_PHIC)

# Generate Noise
np.random.seed(42)
asd_noise = 1e-23 # Flat ASD for simplicity
noise_real = np.random.normal(0, asd_noise, len(freqs))
noise_imag = np.random.normal(0, asd_noise, len(freqs))
noise = noise_real + 1j * noise_imag

# Data = Signal + Noise
data = h_true + noise

# Calculate SNR
snr = np.sqrt(4 * np.sum(np.abs(h_true)**2 / asd_noise**2) * (freqs[1]-freqs[0]))
print(f"Injection SNR: {snr:.2f} (Loud detection)")

# ==============================================================================
# 3. INFERENCE (RECOVERY)
# ==============================================================================

# Likelihood Function
def log_likelihood(params):
    # Unpack
    mc, dl, n, tc, phic = params
    
    # Constraints/Priors (Simple bounds)
    if not (20 < mc < 40): return -np.inf
    if not (100 < dl < 2000): return -np.inf
    if not (0.2 < n < 2.0): return -np.inf
    
    # Template
    h_temp = get_chirp_waveform(freqs, mc, dl, n, tc, phic)
    
    # Residual
    diff = data - h_temp
    
    # LogL (Gaussian noise)
    # -0.5 * sum(|d-h|^2 / sigma^2)
    ll = -0.5 * np.sum(np.abs(diff)**2 / asd_noise**2)
    return ll

# Negative Log Likelihood for minimizer
def nll(params):
    return -log_likelihood(params)

print("\nRunning Parameter Estimation...")

# Search starting point (perturbed)
x0 = [28.0, 400.0, 1.0, 0.01, 0.1] # Start closer to GR (n=1)

# Minimize
# We freeze tc, phic for simplicity to focus on Mc, DL, n correlation
# In real life, MCMC handles all. Here we do optimization.
# Let's optimize [Mc, DL, n] only
def nll_reduced(p):
    return nll([p[0], p[1], p[2], TRUE_TC, TRUE_PHIC])

res = minimize(nll_reduced, [28.0, 400.0, 1.0], 
               bounds=[(20,40), (100,2000), (0.1, 2.0)], 
               method='Nelder-Mead')

print("\nRecovered Maximum Likelihood Estimation:")
print(f"Mc: {res.x[0]:.4f} (True: {TRUE_MC})")
print(f"DL: {res.x[1]:.4f} (True: {TRUE_DL})")
print(f"n:  {res.x[2]:.4f} (True: {TRUE_N})")

# ==============================================================================
# 4. MCMC SCAN (Simulating Posterior)
# ==============================================================================
# To see if n=1 is excluded, we scan the likelihood over 'n' 
# while marginalizing (optimizing) over DL and Mc at each step.

print("\nScanning Likelihood over 'n' parameter...")
n_grid = np.linspace(0.4, 1.2, 50)
log_evidences = []

for n_val in n_grid:
    # Profile likelihood: maximize over Mc, DL for fixed n
    def nll_profile(p): # p = [Mc, DL]
        return nll([p[0], p[1], n_val, TRUE_TC, TRUE_PHIC])
    
    r = minimize(nll_profile, [TRUE_MC, TRUE_DL], bounds=[(20,40), (100,2000)], method='L-BFGS-B')
    log_evidences.append(-r.fun)

log_evidences = np.array(log_evidences)
# Normalize to probability
probs = np.exp(log_evidences - np.max(log_evidences))
probs /= np.sum(probs) * (n_grid[1]-n_grid[0])

# ==============================================================================
# 5. VISUALIZATION
# ==============================================================================
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(n_grid, probs, 'b-', lw=2, label='Posterior P(n|data)')
ax.axvline(TRUE_N, color='g', linestyle='--', label=f'Injection (FIN): n={TRUE_N}')
ax.axvline(1.0, color='r', linestyle=':', label='GR Expectation: n=1.0')

ax.set_xlabel('Amplitude Scaling Exponent $n$')
ax.set_ylabel('Posterior Probability')
ax.set_title(f'QW-1525 PoC: Recovering FIN Scaling\nInjection: n={TRUE_N}, SNR={snr:.1f}')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('QW-1525_Result_Posterior.png')
print("\n[SAVED] QW-1525_Result_Posterior.png")

# ==============================================================================
# 6. VERDICT
# ==============================================================================

# Confidence interval
max_idx = np.argmax(probs)
n_best = n_grid[max_idx]
print(f"\nBest Fit n = {n_best:.3f}")

# Bayes Factor (approximate ratio of peaks/likelihoods)
# Compare Likelihood at n=0.66 vs n=1.0
L_FIN = np.max(log_evidences[np.abs(n_grid - 0.66) < 0.05])
L_GR = np.max(log_evidences[np.abs(n_grid - 1.0) < 0.05])
bayes_factor = np.exp(L_FIN - L_GR)

print(f"Log Likelihood FIN (n~0.66): {L_FIN:.1f}")
print(f"Log Likelihood GR (n=1.0):   {L_GR:.1f}")
print(f"Bayes Factor (FIN/GR):       {bayes_factor:.2e}")

with open("QW-1525_PoC_Result.md", "w") as f:
    f.write(f"# QW-1525 Proof of Concept\n")
    f.write(f"Injection n: {TRUE_N}\n")
    f.write(f"Recovered n: {n_best:.3f}\n")
    f.write(f"Bayes Factor FIN/GR: {bayes_factor:.2e}\n")
    if bayes_factor > 100:
        f.write("\n**CONCLUSION:** The pipeline SUCCESSFULLY recovered the non-standard scaling.\n")
        f.write("A deviation n=0.66 is potentially detectable with SNR~20 events.\n")
    else:
        f.write("\n**CONCLUSION:** The pipeline failed to distinguish scaling.\n")

print("\nDone.")
