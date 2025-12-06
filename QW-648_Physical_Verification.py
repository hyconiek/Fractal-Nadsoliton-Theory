#!/usr/bin/env python3
"""
QW-648_Physical_Verification.py
Purpose: Critical Physical Verification of the Topological-Fractal Model.
         1. Check consistency with QM (Mass-Size relation).
         2. Calculate Tunneling Probabilities (Decay Lifetimes).
"""

import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

# --- Constants ---
ALPHA = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA = 0.01

M_E = 0.511
M_MU = 105.66
M_TAU = 1776.86

# --- Kernel and Potential ---
def K(d):
    return ALPHA * np.cos(OMEGA * d + PHI) / (1 + BETA * d)

# Potential V(d) = - Integral K(x) dx
# We integrate numerically for the Potential Landscape
xs = np.linspace(0, 25, 1000)
ks = K(xs)
# Cumulative trapezoidal integration
vs = -integrate.cumulative_trapezoid(ks, xs, initial=0)

# Shift potential so gen 1 is at 0? Or just leave relative.
# Let's verify depths.

# --- 1. Identify Wells and Barriers ---
d1_approx = 1.333
d2_approx = 9.333
d3_approx = 17.333

def get_potential(d_val):
    idx = np.abs(xs - d_val).argmin()
    return vs[idx]

V1 = get_potential(d1_approx)
V2 = get_potential(d2_approx)
V3 = get_potential(d3_approx)

print("="*60)
print("QW-648/649: PHYSICAL VERIFICATION")
print("="*60)

print("Potential Wells (Depth):")
print(f"  Gen 1 (d={d1_approx:.3f}): V = {V1:.4f}")
print(f"  Gen 2 (d={d2_approx:.3f}): V = {V2:.4f}")
print(f"  Gen 3 (d={d3_approx:.3f}): V = {V3:.4f}")

# Find Barriers (peaks between wells)
# Between d1 and d2
mask12 = (xs > d1_approx) & (xs < d2_approx)
barrier12_h = np.max(vs[mask12])
barrier12_x = xs[mask12][np.argmax(vs[mask12])]

# Between d2 and d3
mask23 = (xs > d2_approx) & (xs < d3_approx)
barrier23_h = np.max(vs[mask23])
barrier23_x = xs[mask23][np.argmax(vs[mask23])]

print(f"\nPotential Barriers:")
print(f"  Barrier 1-2 (at x={barrier12_x:.3f}): Height = {barrier12_h:.4f}")
print(f"  Barrier 2-3 (at x={barrier23_x:.3f}): Height = {barrier23_h:.4f}")

print(f"\nBarrier Heights (Relative to 'decaying' state):")
# Decay 2 -> 1 (Muon decay)
# Barrier height relative to V2
H_mu = barrier12_h - V2
print(f"  Muon Escape Barrier (2->1): {H_mu:.4f}")

# Decay 3 -> 2 (Tau decay)
# Barrier height relative to V3
H_tau = barrier23_h - V3
print(f"  Tau Escape Barrier (3->2): {H_tau:.4f}")

# --- 2. Tunneling Probability (WKB Approximation) ---
# T ~ exp( -2 * Integral sqrt(2m(V(x)-E)) dx / hbar )
# Here we treat 'd' as reaction coordinate. 'm' is 'effective mass' of the reaction coordinate?
# Let's calculate the "WKB Factor" integral: Integral sqrt(V(x) - E_start)
# Assume E_start is slightly above the bottom (Ground state energy ~ hbar*omega/2)
# For rough scaling, assume E_start approx V_bottom.
# So Integral sqrt(V(x) - V_bottom).

def wkb_integral(x_start, x_end, V_base):
    # Integrate sqrt(V(x) - V_base) only where V(x) > V_base
    indices = (xs >= x_start) & (xs <= x_end)
    x_seg = xs[indices]
    v_seg = vs[indices]
    
    # Integrand
    integrand = np.sqrt(np.maximum(0, v_seg - V_base))
    
    integral = np.trapz(integrand, x_seg)
    return integral

# Calculate Action for Muon decay (2->1)
# Tunneling from d2 towards d1 through barrier12
# Find turning points where V(x) = V2 roughly?
# Actually barrier is between 1 and 2.
# Tunneling from 2 to 1:
Action_Mu = wkb_integral(d1_approx, d2_approx, V2)

# Calculate Action for Tau decay (3->2)
Action_Tau = wkb_integral(d2_approx, d3_approx, V3)

print(f"\nTunneling Actions (WKB Integral):")
print(f"  Action_Tau (3->2): {Action_Tau:.4f}")
print(f"  Action_Mu  (2->1): {Action_Mu:.4f}")

# Lifetime ratio ~ exp(2 * (Action_Mu - Action_Tau)) ?
# Actually Gammas ~ exp(-2S).
# Lifetime Tau / Lifetime Mu ~ Gamma_Mu / Gamma_Tau ~ exp(-2 S_Mu) / exp(-2 S_Tau)
# Ratio ~ exp(2 * (Action_Tau - Action_Mu))
# Wait. Longer lifetime means LARGER barrier action.
# Muon lives much longer (10^-6) than Tau (10^-13).
# So Action_Mu should be MUCH LARGER than Action_Tau.

diff_action = Action_Mu - Action_Tau
print(f"\nLifetime Prediction Check:")
print(f"  Difference (Mu - Tau): {diff_action:.4f}")
print(f"  Sign check: {diff_action} > 0? {'✅ YES' if diff_action > 0 else '❌ NO'}")

if diff_action > 0:
    # Estimate magnitude
    # Lifetime ratio ~ 10^7
    # log(10^7) ~ 16.
    # We need 2 * C * diff_action ~ 16.
    # What is the effective coupling constant C?
    print(f"  Muon lives longer than Tau.")
    print(f"  Log Lifetime Ratio (Exp): ln(2.2e-6 / 2.9e-13) = ln(7.6e6) approx 15.8")
    print(f"  Action Difference: {diff_action:.4f}")
    print(f"  Required prefactor 'k': 15.8 / (2 * {diff_action:.4f}) = {15.8 / (2*diff_action):.2f}")
    print(f"  This prefactor relates to effective mass of the field coordinate.")
else:
    print(f"  ❌ MODEL FAILURE: Predicts Tau lives longer than Muon?")
    print(f"     (Action_Mu={Action_Mu} vs Action_Tau={Action_Tau})")
    print(f"     Barrier between 1-2 should be thicker/higher than 2-3.")
    # Check Barrier Heights again
    print(f"     H_mu ({H_mu:.4f}) vs H_tau ({H_tau:.4f})")
    
# --- 3. QM Mass-Size Critique ---
print("\n" + "="*60)
print("QM MASS-SIZE RELATION CRITIQUE")
print("="*60)
print("Model claims M ~ d^2.77")
print("Standard QM: M ~ 1/L (Compton)")
print("Conflict?")

# Interpretation: d is NOT box size L.
# d is "Information Path Length" (Vortex Filament Length).
# A vortex filament has Tension T. Energy E = T * Length.
# So M ~ L_filament.
# But what is the physical SIZE (Radius R) of the particle?
# If the filament is wound up in a ball of radius R...
# Does R correlate with d?

print("Hypothesis: d = Length of Filament.")
print("If scaling is Fractally crumpled:")
print("  Length ~ R^Df (where R is physical radius)")
print("  So M ~ Length ~ R^Df")
print("  But QM says M ~ 1/R")
print("  Contradiction: R^Df vs 1/R")
print("  UNLESS: The 'Stable Radius d' in Kernel space corresponds to MOMENTUM scale k!")

print("\nAlternative Interpretation check:")
print("  Suppose 'd' in Kernel is MOMENTUM MAGNITUDE (k).")
print("  Then Equilibrium at k1, k2, k3.")
print("  Mass ~ k.")
print("  Does mass ratio match k ratios?")
print(f"  k ratios: {d2_approx/d1_approx:.2f}, {d3_approx/d2_approx:.2f}")
print(f"  Mass ratios: {M_MU/M_E:.2f}, {M_TAU/M_MU:.2f}")
print("  Linear match (M ~ k)?")
print(f"  Gen 2: 7.0 vs 206 (No, off by factor 30)")
print(f"  Gen 3: 1.8 vs 16 (No, off by factor 10)")
print("  => 'd' is NOT momentum k linearly.")

print("\nRetaining Vortex Filament Interpretation:")
print("  M ~ d^2.77 (Filament Mass)")
print("  Physical Size R ~ 1/M (Compton)")
print("  The particle is a TIGHTLY WOUND ball.")
print("  Heavier = Tighter = Smaller R = Longer Filament packed inside.")
print("  Consistent if: Packing density increases enormously.")

print("="*60)
