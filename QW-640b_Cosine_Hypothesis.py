#!/usr/bin/env python3
"""
QW-640b: COSINE PROJECTION HYPOTHESIS
Purpose: Test if M_mu and M_tau are corrected by geometric projection factors cos(phi).
Hypothesis: Higher generation particles are 'tilted' in 12D space relative to 3D slice.
            Mass we measure is a projection: m_eff = m_raw * cos(tilt).
Date: 2025-12-06
"""

import numpy as np
from scipy.linalg import eigh
from scipy.stats import entropy

# Constants
ALPHA = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6 # 30 degrees
BETA = 0.01
M_PLANCK = 1.2209e19
C_STAB = 12.027675

# Exp Data
M_E = 0.511
M_MU = 105.66
M_TAU = 1776.86

def compute_mass_projected(N_oct, N_layer, projection_power=0):
    # Base Raw Mass
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
    
    m_raw = M_PLANCK * 1000 * amp_oct * A_res * amp_lay * I_proc
    
    # PROJECTION GEN: m_eff = m_raw * (cos(PHI))^P
    # Why cos(PHI)? phi=30deg is the "twist" of the lattice.
    factor = (np.cos(PHI)) ** projection_power
    
    return m_raw * factor

print("="*60)
print("QW-640b: COSINE PROJECTION TEST")
print(f"Angle PHI = 30 deg, cos(PHI) = {np.cos(PHI):.4f}")
print("="*60)

# Electron: Gen 1, Tilt 0?
m_e = compute_mass_projected(1, 10.0, projection_power=0)
print(f"Electron (Raw): {m_e:.4f} MeV [Exp: {M_E}] Diff: {abs(m_e-M_E)/M_E*100:.2f}%")

# Muon: Gen 2, Tilt 1?
m_mu_1 = compute_mass_projected(4, 9.0, projection_power=1)
print(f"\nMuon (Raw * cos^1): {m_mu_1:.4f} MeV [Exp: {M_MU}] Diff: {abs(m_mu_1-M_MU)/M_MU*100:.2f}%")

m_mu_2 = compute_mass_projected(4, 9.0, projection_power=2)
print(f"Muon (Raw * cos^2): {m_mu_2:.4f} MeV [Exp: {M_MU}] Diff: {abs(m_mu_2-M_MU)/M_MU*100:.2f}%")

# Tau: Gen 3, Tilt 2?
m_tau_1 = compute_mass_projected(7, 8.5, projection_power=1)
print(f"\nTau (Raw * cos^1):  {m_tau_1:.4f} MeV [Exp: {M_TAU}] Diff: {abs(m_tau_1-M_TAU)/M_TAU*100:.2f}%")

m_tau_2 = compute_mass_projected(7, 8.5, projection_power=2)
print(f"Tau (Raw * cos^2):  {m_tau_2:.4f} MeV [Exp: {M_TAU}] Diff: {abs(m_tau_2-M_TAU)/M_TAU*100:.2f}%") # 0.866^2 = 0.75

print("\n" + "="*60)
print("OPTIMIZED TILT SEARCH")
print("="*60)
# Find exact power P such that factor = cos(phi)^P
ratio_mu = M_MU / compute_mass_projected(4, 9.0, 0)
P_mu = np.log(ratio_mu) / np.log(np.cos(PHI))
print(f"Exact Power for Muon: P = {P_mu:.4f}")

ratio_tau = M_TAU / compute_mass_projected(7, 8.5, 0)
P_tau = np.log(ratio_tau) / np.log(np.cos(PHI))
print(f"Exact Power for Tau:  P = {P_tau:.4f}")

print("\nChecking Integer/Half-Integer Hypothesis:")
print("Is P close to 1, 1.5, 2?")
