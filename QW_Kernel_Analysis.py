#!/usr/bin/env python3
"""
QW_Kernel_Analysis.py
Purpose: Analyze the structure of the Coupling Kernel K(d) to find natural 
         resonance points (Octaves) that might correspond to particles.
         Also test Lyapunov stability of these points.
"""

import numpy as np
import matplotlib.pyplot as plt

# --- Constants ---
ALPHA = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA = 0.01

def K(d):
    return ALPHA * np.cos(OMEGA * d + PHI) / (1 + BETA * d)

print("="*60)
print("KERNEL RESONANCE ANALYSIS")
print("="*60)

# 1. Analyze K(d) for d in [0, 24]
d_vals = np.linspace(0, 24, 241)
k_vals = K(d_vals)

print(f"Kernel Function: K(d) = {ALPHA:.2f} * cos({OMEGA:.2f}*d + {PHI:.2f}) / (1 + {BETA}*d)")

# Find Peaks and Valleys
peaks = []
valleys = []

for i in range(1, len(d_vals)-1):
    if k_vals[i-1] < k_vals[i] > k_vals[i+1]:
        peaks.append((d_vals[i], k_vals[i]))
    if k_vals[i-1] > k_vals[i] < k_vals[i+1]:
        valleys.append((d_vals[i], k_vals[i]))

print("\nNatural Resonances (Peaks of Interaction):")
print("d (Octave dist) | Strength (K)")
print("----------------|-------------")
for d, val in peaks:
    print(f"{d:6.2f}          | {val:6.4f}")

print("\nNatural Anti-Resonances (Valleys):")
print("d (Octave dist) | Strength (K)")
print("----------------|-------------")
for d, val in valleys:
    print(f"{d:6.2f}          | {val:6.4f}")

# 2. Map to Generations?
# Electron is Octave 1? Or Octave 0?
# Previous model used N=1, 4, 7.
# Let's see if there are peaks near 1, 4, 7.
# Peak 1: ~7.3? No.

# Let's check Integer Interactions
print("\nInteger Interactions K(n):")
integers = np.arange(0, 13)
k_integers = K(integers)
for n in integers:
    print(f"n={n:<2}: {k_integers[n]:.4f}")

# 3. Mass Hypothesis: M ~ 1 / |K(d)| ?
# If interaction is strong (resonance), the particle is "light" (easy to excite)?
# Or "heavy" (strongly bound)?
# Usually Binding Energy ~ K.
# Free Particle Mass?
# In Solid State: Effective Mass m* ~ 1 / curvature of band.
# Curvature ~ Interaction Strength K.
# So m* ~ 1/K.

print("\nTesting Hypothesis: Mass ~ 1 / |K(n)|")
ref_mass = 0.511 # Electron
# Assuming Electron is at n=1?
k1 = abs(K(1))
m_scale = ref_mass * k1
print(f"Scale derived from Electron (n=1): {m_scale:.4f}")

for n in [1, 2, 3, 4, 5, 6, 7, 8]:
    kn = abs(K(n))
    if kn > 1e-4:
        m_pred = m_scale / kn
        print(f"n={n}: K={kn:.4f} -> M={m_pred:.4f} MeV")
    else:
        print(f"n={n}: K={kn:.4f} -> M=INF")

# 4. Testing Hypothesis: Mass ~ exp(d) * K(d)?
# From QW-641: m ~ beta^d?
# Let's try to find if ANY integer n matches Muon (105) and Tau (1776).
print("\nReverse Search for Muon (105) and Tau (1776):")
target_mu = 105.66
target_tau = 1776.86

for n in range(1, 20):
    kn = abs(K(n))
    # Test m ~ 1/K
    m1 = m_scale / kn
    
    # Test m ~ 1/K^2
    m2 = m_scale / (kn**2) * k1 # re-normalize
    
    # Test m ~ sum(K) up to n? (Cumulative)
    
    print(f"n={n}: m(1/K)={m1:.1f}, m(1/K^2)={m2:.1f}")
    
    if abs(m1 - target_mu)/target_mu < 0.2: print(f"  >>> MATCH MUON (1/K) at n={n}")
    if abs(m2 - target_mu)/target_mu < 0.2: print(f"  >>> MATCH MUON (1/K^2) at n={n}")
    
    if abs(m1 - target_tau)/target_tau < 0.2: print(f"  >>> MATCH TAU (1/K) at n={n}")
    if abs(m2 - target_tau)/target_tau < 0.2: print(f"  >>> MATCH TAU (1/K^2) at n={n}")

# 5. Lyapunov Stability (Mockup)
# Real Lyapunov requires time evolution.
# Heuristic: Stability ~ K(d).
# If K(d) is positive -> Attractor? Negative -> Repulsor?
# Local minima of Potential V(d) ~ -Integral K(x) dx
print("\nPotential landscape V(d):")
print("(Assuming Force F ~ K(d))")
# Simple integration
V = [0]
for i in range(len(d_vals)-1):
    V_next = V[-1] - k_vals[i] * (d_vals[1]-d_vals[0]) # Minus because F = -dV/dx and assume K is force-like
    V.append(V_next)

# Find local minima of V (Stable Orbits)
minima = []
for i in range(1, len(d_vals)-1):
    if V[i-1] > V[i] < V[i+1]:
        minima.append((d_vals[i], V[i]))

print("Stable Orbits (Minima of V):")
for d, val in minima:
    print(f"d = {d:.2f}")

print("="*60)
