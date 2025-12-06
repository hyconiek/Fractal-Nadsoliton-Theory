#!/usr/bin/env python3
"""
QW-647_Vortex_Volume.py
Purpose: Verify the Vortex Volume Hypothesis: Mass ~ Radius^Df
         using stable equilibrium points from QW-646.
"""

import numpy as np
import scipy.optimize as optimize

print("="*60)
print("QW-647: VORTEX VOLUME SCALING ANALYSIS")
print("="*60)

# --- Data from QW-646 (Stable Equilibria) ---
# Kernel: 2.77 * cos(0.79*d + 0.52) / (1 + 0.01*d)
D1 = 1.3333
D2 = 9.3333
D3 = 17.3333

print(f"Stable Radii (from QW-646):")
print(f"  r1 (Gen 1) = {D1:.4f}")
print(f"  r2 (Gen 2) = {D2:.4f}")
print(f"  r3 (Gen 3) = {D3:.4f}")

# --- Experimental Mass Ratios ---
RATIO_MU_E = 105.66 / 0.511      # 206.77
RATIO_TAU_MU = 1776.86 / 105.66  # 16.82

print(f"\nTarget Mass Ratios:")
print(f"  Mu/E   = {RATIO_MU_E:.2f}")
print(f"  Tau/Mu = {RATIO_TAU_MU:.2f}")

# --- 1. Calculate Fractal Dimension Df ---
# Hypothesis: M ~ r^Df
# M2 / M1 = (r2 / r1)^Df
# 206.77  = (9.3333 / 1.3333)^Df
# 206.77  = (7.0)^Df
# Df = log(206.77) / log(7.0)

ratio_r2_r1 = D2 / D1
Df_calc = np.log(RATIO_MU_E) / np.log(ratio_r2_r1)

print(f"\n1. Calculating Fractal Dimension (from Gen 1 -> Gen 2):")
print(f"   Radial Ratio r2/r1 = {ratio_r2_r1:.4f}")
print(f"   Mass Ratio M2/M1   = {RATIO_MU_E:.4f}")
print(f"   implied Df = log({RATIO_MU_E:.2f}) / log({ratio_r2_r1:.2f})")
print(f"   Df = {Df_calc:.6f}")

# --- 2. Predict Gen 3 (Tau) ---
# Prediction: M3 / M2 = (r3 / r2)^Df
ratio_r3_r2 = D3 / D2
pred_ratio_tau_mu = ratio_r3_r2 ** Df_calc

print(f"\n2. Predicting Gen 3 (Tau/Mu Ratio):")
print(f"   Radial Ratio r3/r2 = {ratio_r3_r2:.4f}")
print(f"   Using Df = {Df_calc:.4f}")
print(f"   Predicted Ratio = {ratio_r3_r2:.4f}^{Df_calc:.4f} = {pred_ratio_tau_mu:.4f}")
print(f"   Actual Ratio    = {RATIO_TAU_MU:.4f}")

error_tau = abs(pred_ratio_tau_mu - RATIO_TAU_MU) / RATIO_TAU_MU * 100
print(f"   Error: {error_tau:.2f}%")

# --- 3. Is Df = 3 (Volumetric)? ---
# Or is it Df = Alpha?
print(f"\n3. Interpretation of Df = {Df_calc:.4f}")
print(f"   Is it 3.0? (Euclidean Volume): Error = {abs(Df_calc - 3.0)/3.0*100:.1f}%")
print(f"   Is it Alpha? (Information Capacity 2.772): Error = {abs(Df_calc - 2.772)/2.772*100:.1f}%")
print(f"   Is it e? (2.718): Error = {abs(Df_calc - np.e)/np.e*100:.1f}%")

# --- 4. Global Fit Optimization ---
# Find best Df that minimizes total error
def total_error(df):
    p1 = (D2/D1)**df
    p2 = (D3/D2)**df
    e1 = abs(p1 - RATIO_MU_E)/RATIO_MU_E
    e2 = abs(p2 - RATIO_TAU_MU)/RATIO_TAU_MU
    return e1 + e2

res = optimize.minimize_scalar(total_error, bounds=(2, 4), method='bounded')
best_Df = res.x
min_err = res.fun

print(f"\n4. Optimized Global Df:")
print(f"   Best Df = {best_Df:.6f}")
print(f"   Total Relative Error metric: {min_err:.4f}")

# Check predictions with Best Df
pred_mu_rat = (D2/D1)**best_Df
pred_tau_rat = (D3/D2)**best_Df
print(f"   Pred Mu/E Ratio: {pred_mu_rat:.2f} (Target {RATIO_MU_E:.2f}, Err {abs(pred_mu_rat-RATIO_MU_E)/RATIO_MU_E*100:.1f}%)")
print(f"   Pred Tau/Mu Ratio: {pred_tau_rat:.2f} (Target {RATIO_TAU_MU:.2f}, Err {abs(pred_tau_rat-RATIO_TAU_MU)/RATIO_TAU_MU*100:.1f}%)")

# --- 5. Physical Meaning ---
print("\n" + "="*60)
print("CONCLUSION")
print("="*60)

if min_err < 0.5: # < 50% combined error suggests viability
    print("✅ HYPOTHESIS PLAUSIBLE")
    print(f"   Mass scales as M ~ r^{best_Df:.2f}")
    if abs(best_Df - 2.772) < 0.1:
        print("   AMAZING MATCH: Df ≈ Alpha_Geo (4ln2 = 2.772)!")
        print("   Theory: Mass scales with Information Dimension (Alpha)!")
    elif abs(best_Df - 3.0) < 0.2:
        print("   MATCH: Df ≈ 3 (Volumetric Scaling)")
    else:
        print(f"   Df is a fractional dimension ~ {best_Df:.2f}")
else:
    print("❌ HYPOTHESIS FAILED")
    print("   Scaling law does not hold across generations.")

print("="*60)
