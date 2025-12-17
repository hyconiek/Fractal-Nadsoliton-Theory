#!/usr/bin/env python3
"""
QW-1524: 3D NEURO-GRAVITY WAVE SIMULATION
==========================================
Testing the "Gravity as Weight Update Propagation" hypothesis in full 3D.

HYPOTHESIS:
- FIN Vacuum is a Neural Network.
- Mass is a topological defect (active region).
- Gravity is the propagation of weight updates (learning).
- In 1D (QW-1523), we saw amplification (n < 0).
- In 3D, geometric spreading (1/r) should compete with neural amplification.

Does 3D spreading win? Do we get A ~ 1/r?

METHODOLOGY:
1. 3D Grid of neurons (Simulating vacuum).
2. Central disturbance (Merger).
3. Neural Wave Equation (damped wave eq for potential u).
4. Measure amplitude vs distance.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import laplace
import datetime

print("="*80)
print("QW-1524: 3D NEURO-GRAVITY WAVE SIMULATION")
print("="*80)

# Parameters
N = 64  # 3D Grid size (64^3 = 262k neurons)
L = 40.0
dx = L / N
dt = 0.05
t_max = 100.0
n_steps = int(t_max / dt)

# FIN Parameters
ALPHA_GEO = 4 * np.log(2)
BETA_TORS = 0.01  # Damping

# Neural Wave Speed (Topology limit)
c_neural = 1.0

print(f"Grid: {N}x{N}x{N}")
print(f"Neural Speed c={c_neural}")
print(f"Damping β={BETA_TORS}")

# Fields
u = np.zeros((N, N, N))
u_dot = np.zeros((N, N, N))

# Source
center = N // 2
f_gw = 0.15

# Detectors (Indices relative to center)
radii_idx = [3, 6, 9, 12, 15, 18, 21]
det_histories = {r: [] for r in radii_idx}

print("\nRunning 3D Simulation...")

for step in range(n_steps):
    t = step * dt
    
    # 1. Source (Center)
    source = np.zeros((N, N, N))
    source[center, center, center] = np.sin(2 * np.pi * f_gw * t)
    
    # 2. Laplacian (3D Connectivity)
    # Using scipy.ndimage for speed
    lap = laplace(u) / dx**2
    
    # 3. Neural Wave Equation
    # u_tt = c^2 Δu - β u_t + S
    u_ddot = c_neural**2 * lap - BETA_TORS * u_dot + source
    
    u_dot += u_ddot * dt
    u += u_dot * dt
    
    # 4. Record at detectors
    # Average over spherical shells to reduce grid noise
    if step % 2 == 0: # Optimize recording
        for r_idx in radii_idx:
            # Simple point detector along axis for speed
             # (Averaging shells is expensive in loop)
            val = u[center + r_idx, center, center]
            det_histories[r_idx].append(val)
            
    if step % (n_steps // 10) == 0:
        print(f"  Step {step}/{n_steps}, t={t:.1f}")

print("Simulation Complete.")

# Analysis
print("\n[AMPLITUDE ANALYSIS]")
amplitudes = []
radii_phys = []

for r_idx in radii_idx:
    r_phys = r_idx * dx
    sig = np.array(det_histories[r_idx])
    # Cut transient
    sig = sig[int(len(sig)/2):]
    amp = np.std(sig)
    
    amplitudes.append(amp)
    radii_phys.append(r_phys)
    print(f"r={r_phys:.2f}: Amp={amp:.6f}")

amplitudes = np.array(amplitudes)
radii_phys = np.array(radii_phys)

# Fit
valid = amplitudes > 1e-10
if np.sum(valid) >= 3:
    log_r = np.log(radii_phys[valid])
    log_A = np.log(amplitudes[valid])
    coeffs = np.polyfit(log_r, log_A, 1)
    n_fit = -coeffs[0]
    print(f"\nScaling: A ~ 1/r^{n_fit:.2f}")
else:
    n_fit = 0
    print("Insufficient data")

# Visualization
fig, ax = plt.subplots(figsize=(8, 6))
ax.loglog(radii_phys[valid], amplitudes[valid], 'bo-', label='Measured (3D)')
if len(radii_phys) > 0:
    r_ref = radii_phys[valid]
    A_ref = amplitudes[valid][0]
    ax.loglog(r_ref, A_ref * (r_ref[0]/r_ref)**1.0, 'g--', label='1/r (Ideal GW)')
    ax.loglog(r_ref, A_ref * (r_ref[0]/r_ref)**0.5, 'r:', label='1/sqrt(r)')

ax.set_title(f'QW-1524: 3D Neuro-Gravity Scaling (n={n_fit:.2f})')
ax.set_xlabel('Distance')
ax.set_ylabel('Amplitude')
ax.legend()
ax.grid(True)

plt.savefig('QW-1524_3D_NeuroGravity.png')
print("[SAVED] QW-1524_3D_NeuroGravity.png")

# Verdict
print("\nVERDICT for QW-1524:")
if abs(n_fit - 1.0) < 0.2:
    print("✅ SUCCESS: 3D Geometry restores 1/r scaling!")
    print("   The Theory is saved by dimensionality.")
elif n_fit < 0.5:
    print("❌ FAILURE: Still anomalous scaling.")
    print("   The vacuum is too active/amplifying.")
else:
    print(f"🟡 PARTIAL: n={n_fit:.2f}")

# Report
with open('QW-1524_3D_NeuroGravity_Result.md', 'w') as f:
    f.write(f"# QW-1524 Result\n\nScale n = {n_fit:.2f}\n")
    if abs(n_fit - 1.0) < 0.2:
        f.write("Verdict: SUCCESS (1/r recovered in 3D)")
    else:
        f.write("Verdict: FAILURE")
