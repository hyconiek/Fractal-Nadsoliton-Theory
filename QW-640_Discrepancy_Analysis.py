#!/usr/bin/env python3
"""
QW-640: EXACT MISSING FACTOR ANALYSIS
Purpose: Calculate precise discrepancy for Muon and Tau to identify geometric source.
Date: 2025-12-06
"""

import numpy as np
from scipy.linalg import eigh
from scipy.stats import entropy

# Constants
ALPHA = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA = 0.01
M_PLANCK = 1.2209e19
C_STAB = 12.027675

# Exp Data
M_E = 0.511
M_MU = 105.66
M_TAU = 1776.86

def compute_mass_raw(N_oct, N_layer):
    # Same physics engine as Electron
    kappa = ALPHA / (OMEGA * PHI)
    amp_oct = kappa ** (N_oct/12)
    
    H = np.zeros((12, 12))
    def K(d): return ALPHA * np.cos(OMEGA * d + PHI) / (1 + BETA * d)
    for i in range(12):
        for j in range(12): H[i, j] = -K(abs(i - j))
    evals, evecs = eigh(H)
    psi = np.zeros(12); psi[N_oct] = 1.0
    A_res = abs(np.dot(psi, evecs[:, 0]))
    
    amp_lay = BETA ** N_layer
    
    psi_dist = np.exp(-0.5 * ((np.arange(12) - N_oct) / 0.8)**2)
    psi_dist /= np.sum(psi_dist)
    I_proc = (entropy(psi_dist, base=2) * 0.1 / C_STAB)
    
    return M_PLANCK * 1000 * amp_oct * A_res * amp_lay * I_proc

print("="*60)
print("QW-640: DISCREPANCY ANALYSIS")
print("="*60)

# 1. Electron Baseline
m_e_pred = compute_mass_raw(1, 10.0)
print(f"Electron (O1, L10): {m_e_pred:.4f} MeV [Exp: {M_E}] Diff: {abs(m_e_pred-M_E)/M_E*100:.2f}%")

# 2. Muon Analysis
m_mu_pred = compute_mass_raw(4, 9.0)
missing_mu = M_MU / m_mu_pred
print(f"\nMuon (O4, L9):     {m_mu_pred:.4f} MeV [Exp: {M_MU}]")
print(f"  -> Missing Factor: {missing_mu:.6f}")

# 3. Tau Analysis
m_tau_pred = compute_mass_raw(7, 8.5)
missing_tau = M_TAU / m_tau_pred
print(f"\nTau (O7, L8.5):    {m_tau_pred:.4f} MeV [Exp: {M_TAU}]")
print(f"  -> Missing Factor: {missing_tau:.6f}")

# 4. Hypothesis Testing
print("\n" + "-"*60)
print("GEOMETRIC CANDIDATES for Missing Factors:")
print("-" * 60)

candidates = {
    "4/pi": 4/np.pi,
    "pi/4": np.pi/4,
    "sqrt(2)": np.sqrt(2),
    "1/sin(phi) (2.0)": 1/np.sin(PHI),
    "1/cos(phi) (1.15)": 1/np.cos(PHI),
    "cos(phi) (0.866)": np.cos(PHI),
    "ln(2)": np.log(2),
    "1/ln(2)": 1/np.log(2)
}

print(f"{'Candidate':<20} {'Value':<10} {'Diff Mu (%)':<12} {'Diff Tau (%)':<12}")
for name, val in candidates.items():
    d_mu = abs(missing_mu - val)/missing_mu*100
    d_tau = abs(missing_tau - val)/missing_tau*100
    print(f"{name:<20} {val:.4f}     {d_mu:.1f}%          {d_tau:.1f}%")

print("\nSpecific Check: Muon needs factor 0.84? (Inverse?)")
inv_missing_mu = 1/missing_mu
print(f"Inverse Missing Mu: {inv_missing_mu:.4f}")
# Maybe relativistic gamma?
# Or maybe the layer step is not EXACTLY 1.0 but depends on phi?
