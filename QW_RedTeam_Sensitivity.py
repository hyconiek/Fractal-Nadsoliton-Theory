#!/usr/bin/env python3
"""
QW-640_RedTeam_Sensitivity.py
Purpose: Analyze the sensitivity of the mass derivation model in QW-640c 
         to small variations in input parameters (N_layer, C_STAB, PHI).
         This is a Red Team tool to detect overfitting.
"""

import numpy as np
from scipy.linalg import eigh
from scipy.stats import entropy
import matplotlib.pyplot as plt

# --- Original Logic from QW-640c ---

ALPHA = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6 # 30 degrees
BETA = 0.01
M_PLANCK = 1.2209e19
# Original C_STAB from QW-640
C_STAB_ORIGINAL = 12.027675

def compute_mass_redteam(N_oct, N_layer, c_stab_val=C_STAB_ORIGINAL, phi_val=PHI, is_excited=False):
    # Reconstruct mass calculation with variable parameters
    kappa = ALPHA / (OMEGA * phi_val)
    amp_oct = kappa ** (N_oct/12.0)
    
    H = np.zeros((12, 12))
    def K(d): return ALPHA * np.cos(OMEGA * d + phi_val) / (1 + BETA * d)
    for i in range(12):
        for j in range(12): H[i, j] = -K(abs(i - j))
    
    # Simple check for stability of diagonalization
    try:
        evals, evecs = eigh(H)
    except:
        return 0.0
        
    psi = np.zeros(12); psi[N_oct] = 1.0
    A_res = abs(np.dot(psi, evecs[:, 0]))
    
    amp_lay = BETA ** N_layer
    
    psi_dist = np.exp(-0.5 * ((np.arange(12) - N_oct) / 0.8)**2)
    psi_dist /= np.sum(psi_dist)
    I_proc = (entropy(psi_dist, base=2) * 0.1 / c_stab_val)
    
    m_raw = M_PLANCK * 1000 * amp_oct * A_res * amp_lay * I_proc
    
    if is_excited:
        factor = np.cos(phi_val)
        m_final = m_raw * factor
    else:
        m_final = m_raw
        
    return m_final

# --- Analysis ---

print("="*60)
print("RED TEAM SENSITIVITY ANALYSIS: QW-640c Model")
print("="*60)

# Target values
M_TAU_EXP = 1776.86
M_MU_EXP = 105.66

print(f"Baseline C_STAB: {C_STAB_ORIGINAL}")
print(f"Baseline PHI:    {PHI} ({np.degrees(PHI)} deg)")
print("-" * 60)

# 1. Parameter Sensitivity for TAU (N_oct=7, N_layer=8.5)
print("\n[TEST 1] TAU MASS SENSITIVITY (Target: 1776.86 MeV)")
print("Nominal parameters: N_layer=8.5, C_STAB=12.027675")

base_tau = compute_mass_redteam(7, 8.5, is_excited=True)
print(f"Base Tau Mass: {base_tau:.4f} MeV")

# Vary N_layer
print("\nSensitivity to N_layer (Layer Selection):")
for delta in [-0.5, -0.1, -0.01, 0.01, 0.1, 0.5]:
    val = 8.5 + delta
    m = compute_mass_redteam(7, val, is_excited=True)
    diff_pct = (m - base_tau) / base_tau * 100
    print(f"  N_layer = {val:<6} -> Mass = {m:.4f} MeV (Change: {diff_pct:+.2f}%)")

# Vary C_STAB
print("\nSensitivity to C_STAB (Normalization Constant):")
for delta_pct in [-5, -1, 1, 5]:
    val = C_STAB_ORIGINAL * (1 + delta_pct/100)
    m = compute_mass_redteam(7, 8.5, c_stab_val=val, is_excited=True)
    diff_pct = (m - base_tau) / base_tau * 100
    print(f"  C_STAB = {val:.4f} ({delta_pct:+}%) -> Mass = {m:.4f} MeV (Change: {diff_pct:+.2f}%)")

# 2. Logic Check: Layer Quantization
print("\n[TEST 2] LOGIC CHECK - Is N_layer=8.5 justified?")
print("The theory claims layers are quantized integers usually. Why 8.5?")
print("Checking N=8 and N=9 for Tau:")
m8 = compute_mass_redteam(7, 8.0, is_excited=True)
m9 = compute_mass_redteam(7, 9.0, is_excited=True)
print(f"  N=8.0 -> {m8:.4f} MeV (vs Exp {M_TAU_EXP}: {(m8-M_TAU_EXP)/M_TAU_EXP*100:+.2f}%)")
print(f"  N=9.0 -> {m9:.4f} MeV (vs Exp {M_TAU_EXP}: {(m9-M_TAU_EXP)/M_TAU_EXP*100:+.2f}%)")
print("  => Heavy reliance on half-integer layer 8.5 for fit.")

# 3. Parameter Sensitivity for MUON (N_oct=4, N_layer=9.0)
print("\n[TEST 3] MUON MASS SENSITIVITY (Target: 105.66 MeV)")
print("Nominal parameters: N_layer=9.0")
base_mu = compute_mass_redteam(4, 9.0, is_excited=True)
print(f"Base Muon Mass: {base_mu:.4f} MeV")

print("\nSensitivity to N_layer:")
for delta in [-0.1, 0.1]:
    val = 9.0 + delta
    m = compute_mass_redteam(4, val, is_excited=True)
    diff_pct = (m - base_mu) / base_mu * 100
    print(f"  N_layer = {val:<6} -> Mass = {m:.4f} MeV (Change: {diff_pct:+.2f}%)")

# 4. Tautology Check on C_STAB
print("\n[TEST 4] TAUTOLOGY CHECK")
# Does C_STAB just cancel out an error?
# Formula: m ~ 1/C_STAB
# If we chose C_STAB to fix Electron mass, does it hold for others naturally?
# Let's see what happens if C_STAB was 1.0 (natural)
m_e_natural = compute_mass_redteam(1, 10.0, c_stab_val=1.0, is_excited=False)
ratio = m_e_natural / 0.511
print(f"Mass of Electron if C_STAB=1.0: {m_e_natural:.4f} MeV")
print(f"Factor needed to fix Electron: {ratio:.4f}")
print(f"Used C_STAB: {C_STAB_ORIGINAL}")
print(f"Difference: {abs(ratio - C_STAB_ORIGINAL):.6f}")
print("=> C_STAB is exactly calibrating the electron mass. It works for others only if the ratios are perfect.")

print("="*60)
