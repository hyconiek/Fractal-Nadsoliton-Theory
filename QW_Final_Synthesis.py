#!/usr/bin/env python3
"""
QW_Final_Synthesis.py
Purpose: Calculate precise mass predictions for the Topological-Fractal Model.
         Formula: M(n) = M_scale * |W_n| * (d_n)^Alpha
         Where d_n are stable equilibria of the Kernel.
"""

import numpy as np
import scipy.optimize as optimize

# --- Constants ---
ALPHA = 4 * np.log(2)  # 2.772588...
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA = 0.01

print("="*60)
print("QW FINAL SYNTHESIS: TOPOLOGICAL-FRACTAL MODEL")
print("="*60)

# --- 1. Kernel Equilibria (Recalculating for precision) ---
def K(d):
    return ALPHA * np.cos(OMEGA * d + PHI) / (1 + BETA * d)

roots = []
# Find first 3 stable roots (Attractors: Slope < 0)
x_scan = np.linspace(0, 25, 500)
for i in range(len(x_scan)-1):
    if K(x_scan[i]) > 0 and K(x_scan[i+1]) < 0: # Zero crossing downwards
        root = optimize.brentq(K, x_scan[i], x_scan[i+1])
        roots.append(root)

if len(roots) < 3:
    print("Error: Could not find 3 stable roots.")
    exit()

d1 = roots[0]
d2 = roots[1]
d3 = roots[2]

print(f"Stable Radii (d_n):")
print(f"  Gen 1: {d1:.6f}")
print(f"  Gen 2: {d2:.6f}")
print(f"  Gen 3: {d3:.6f}")

# --- 2. Model Parameters ---
# Winding Numbers (Hypothesis)
W1 = 1
W2 = 1
W3 = 3

# Exponent (Hypothesis)
Df = ALPHA # 2.7726

print(f"\nModel Parameters:")
print(f"  Fractal Dimension: Df = Alpha = {Df:.6f}")
print(f"  Winding Numbers: W=[{W1}, {W2}, {W3}]")

# --- 3. Predictions ---
# Calibrate scale M0 using Electron Mass
M_E_EXP = 0.511
M_MU_EXP = 105.66
M_TAU_EXP = 1776.86

# Formula: M = M0 * W * d^Df
# M0 = M_E / (W1 * d1^Df)
gamma_1 = W1 * (d1 ** Df)
M0 = M_E_EXP / gamma_1

print(f"\nCalibration (normalized to Electron):")
print(f"  Gamma_1 (W={W1}, d={d1:.4f})  = {gamma_1:.6f}")
print(f"  Scale M0 = {M_E_EXP} / {gamma_1:.6f} = {M0:.6f} MeV")

# Predict Muon
gamma_2 = W2 * (d2 ** Df)
M_mu_pred = M0 * gamma_2
err_mu = abs(M_mu_pred - M_MU_EXP) / M_MU_EXP * 100

print(f"\nGen 2 (Muon) Prediction:")
print(f"  W={W2}, d={d2:.4f} -> Gamma={gamma_2:.6f}")
print(f"  M_pred = {M0:.6f} * {gamma_2:.6f} = {M_mu_pred:.4f} MeV")
print(f"  Exp    = {M_MU_EXP:.4f} MeV")
print(f"  Error  = {err_mu:.2f}%")

# Predict Tau
gamma_3 = W3 * (d3 ** Df)
M_tau_pred = M0 * gamma_3
err_tau = abs(M_tau_pred - M_TAU_EXP) / M_TAU_EXP * 100

print(f"\nGen 3 (Tau) Prediction:")
print(f"  W={W3}, d={d3:.4f} -> Gamma={gamma_3:.6f}")
print(f"  M_pred = {M0:.6f} * {gamma_3:.6f} = {M_tau_pred:.4f} MeV")
print(f"  Exp    = {M_TAU_EXP:.4f} MeV")
print(f"  Error  = {err_tau:.2f}%")

# --- 4. Sensitivity to Alpha ---
print(f"\nSensitivity Analysis (Alpha +/- 1%):")
for factor in [0.99, 1.01]:
    d_f_sens = Df * factor
    # Recalculate M0
    g1 = W1 * (d1 ** d_f_sens)
    m0_sens = M_E_EXP / g1
    
    # Predict others
    m_mu = m0_sens * W2 * (d2 ** d_f_sens)
    m_tau = m0_sens * W3 * (d3 ** d_f_sens)
    
    e_mu = abs(m_mu - M_MU_EXP)/M_MU_EXP*100
    e_tau = abs(m_tau - M_TAU_EXP)/M_TAU_EXP*100
    
    print(f"  Alpha * {factor:.2f}: Mu Err={e_mu:.1f}%, Tau Err={e_tau:.1f}%")

# --- 5. Conclusion ---
print("\n" + "="*60)
avg_err = (err_mu + err_tau) / 2
if avg_err < 5.0:
    print(f"✅ FINAL VERDICT: SUCCESS (Avg Error {avg_err:.2f}%)")
    print("   The Topological-Fractal Model explains the hierarchy!")
    print(f"   Formula: M_n = M_0 * W_n * (d_n)^Alpha")
else:
    print(f"❌ FINAL VERDICT: FAILURE (Avg Error {avg_err:.2f}%)")

print("="*60)
