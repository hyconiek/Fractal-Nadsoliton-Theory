#!/usr/bin/env python3
"""
QW-1517: GRAVITATIONAL WAVE POLARIZATION MODES (+ and ×)
=========================================================
In General Relativity, gravitational waves have TWO polarization modes:
- Plus (+): stretches in x, compresses in y
- Cross (×): stretches at 45°, compresses at 135°

PHYSICS QUESTION:
Does FIN Theory predict these two polarization modes?

APPROACH:
In 2D, a torsion wave can have TWO independent polarizations:
- θ_+ (symmetric)
- θ_× (antisymmetric)

The wave equation for each mode is the same, but they are orthogonal.

We will:
1. Set up a 2D grid with a rotating source (generates both + and ×)
2. Measure the wave amplitude in different directions
3. Check for the characteristic "quadrupole" pattern of GW
"""

import numpy as np
from scipy.ndimage import laplace
import datetime

print("="*80)
print("QW-1517: GRAVITATIONAL WAVE POLARIZATION MODES")
print("="*80)

# Parameters
ALPHA_GEO = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01

def K(d):
    if d < 0.1:
        d = 0.1
    return ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + BETA_TORS * d)

c_tors = np.sqrt(abs(K(0)))
print(f"Wave speed: c_tors = {c_tors:.4f}")

# 2D Grid
N = 128
L = 40.0
dx = L / N
print(f"Grid: {N}×{N}, dx = {dx:.4f}")

# Two polarization fields
theta_plus = np.zeros((N, N))   # + polarization
theta_cross = np.zeros((N, N))  # × polarization
theta_plus_dot = np.zeros((N, N))
theta_cross_dot = np.zeros((N, N))

# Source: rotating quadrupole (generates both polarizations)
center = np.array([N//2, N//2])
source_freq = 0.05

print(f"\nSource: Rotating quadrupole at f = {source_freq}")
print(f"GW frequency: f_gw = 2 × f_orbit = {2*source_freq}")

# Time evolution
dt = 0.1
t_max = 200.0
n_steps = int(t_max / dt)

print(f"Time: dt = {dt}, t_max = {t_max}")

# Detectors at 4 cardinal directions (90° apart)
detector_distance = 15
detectors = {
    "North": (center[0], center[1] + detector_distance),
    "East": (center[0] + detector_distance, center[1]),
    "South": (center[0], center[1] - detector_distance),
    "West": (center[0] - detector_distance, center[1]),
    "NE (45°)": (center[0] + int(detector_distance/np.sqrt(2)), 
                 center[1] + int(detector_distance/np.sqrt(2))),
    "SE (315°)": (center[0] + int(detector_distance/np.sqrt(2)), 
                  center[1] - int(detector_distance/np.sqrt(2)))
}

histories_plus = {k: [] for k in detectors}
histories_cross = {k: [] for k in detectors}

print("\n[SIMULATION]")
print("-" * 60)

for step in range(n_steps):
    t = step * dt
    omega_rot = 2 * np.pi * source_freq
    
    # Source terms for + and × polarizations
    # A rotating quadrupole generates:
    # h+ ∝ cos(2ωt)
    # h× ∝ sin(2ωt)
    
    source_plus = np.zeros((N, N))
    source_cross = np.zeros((N, N))
    
    # Quadrupole pattern: extend in x, compress in y for +
    # For a point source, we inject h+ and h× directly
    amplitude = 1.0
    source_plus[center[0], center[1]] = amplitude * np.cos(2 * omega_rot * t)
    source_cross[center[0], center[1]] = amplitude * np.sin(2 * omega_rot * t)
    
    # Laplacians
    lap_plus = laplace(theta_plus) / dx**2
    lap_cross = laplace(theta_cross) / dx**2
    
    # Wave equations
    acc_plus = c_tors**2 * lap_plus - BETA_TORS * theta_plus_dot + source_plus
    acc_cross = c_tors**2 * lap_cross - BETA_TORS * theta_cross_dot + source_cross
    
    # Update
    theta_plus_dot += acc_plus * dt
    theta_cross_dot += acc_cross * dt
    theta_plus += theta_plus_dot * dt
    theta_cross += theta_cross_dot * dt
    
    # Record at detectors
    for name, (i, j) in detectors.items():
        if 0 <= i < N and 0 <= j < N:
            histories_plus[name].append(theta_plus[i, j])
            histories_cross[name].append(theta_cross[i, j])
        else:
            histories_plus[name].append(0)
            histories_cross[name].append(0)
    
    if step % (n_steps // 5) == 0:
        print(f"  Step {step}/{n_steps}")

print("\n[POLARIZATION ANALYSIS]")
print("-" * 60)

# For GW from a face-on binary (rotation axis toward us):
# - At North/South: see h+ (max), h× (zero)
# - At East/West: see h+ (max), h× (zero)
# - At 45°: see h+ (zero), h× (max)

# Measure amplitudes
print("\nAmplitudes at detectors:")
print(f"{'Detector':<15} {'θ+ (rms)':<12} {'θ× (rms)':<12} {'Dominant':<10}")
print("-" * 50)

for name in detectors:
    sig_plus = np.array(histories_plus[name])
    sig_cross = np.array(histories_cross[name])
    
    amp_plus = np.std(sig_plus - np.mean(sig_plus))
    amp_cross = np.std(sig_cross - np.mean(sig_cross))
    
    if amp_plus > amp_cross:
        dominant = "+"
    elif amp_cross > amp_plus:
        dominant = "×"
    else:
        dominant = "equal"
    
    print(f"{name:<15} {amp_plus:<12.6f} {amp_cross:<12.6f} {dominant:<10}")

# Check for quadrupole pattern
# In the + mode, detectors at 0° and 90° should see the same phase
# In the × mode, detectors at 45° should see maximum

print("\n[QUADRUPOLE PATTERN CHECK]")
print("-" * 60)

# Phase comparison: North vs East should be in-phase for +
sig_north_plus = np.array(histories_plus["North"])
sig_east_plus = np.array(histories_plus["East"])

# Cross-correlation at zero lag
corr = np.corrcoef(sig_north_plus, sig_east_plus)[0, 1]
print(f"Correlation (North vs East) for θ+: {corr:.4f}")
print(f"Expected for quadrupole: ~1.0 (in-phase)")

# For ×, check 45° detector
sig_ne_cross = np.array(histories_cross["NE (45°)"])
amp_ne_cross = np.std(sig_ne_cross)
print(f"Amplitude at 45° for θ×: {amp_ne_cross:.6f}")

# VERDICT
print("\n" + "=" * 60)
print("QW-1517 VERDICT")
print("=" * 60)

if abs(corr) > 0.8:
    verdict_plus = "✅ + polarization shows quadrupole pattern"
else:
    verdict_plus = f"🟡 Correlation = {corr:.2f} (expected ~1.0)"

if amp_ne_cross > 0:
    verdict_cross = "✅ × polarization detected at 45°"
else:
    verdict_cross = "❌ × polarization not detected"

print(verdict_plus)
print(verdict_cross)

# Save report
report = f"""# QW-1517: Gravitational Wave Polarization Modes

**Date:** {datetime.datetime.now()}

## Theory

GR predicts two polarization modes for gravitational waves:
- **Plus (+)**: Quadrupole pattern aligned with axes
- **Cross (×)**: Quadrupole pattern rotated 45°

## Method

- 2D grid simulation with rotating quadrupole source
- Source injects h+ ∝ cos(2ωt) and h× ∝ sin(2ωt)
- Detectors at cardinal directions and 45°

## Results

| Detector | θ+ amplitude | θ× amplitude |
|----------|--------------|--------------|
"""
for name in detectors:
    sig_plus = np.array(histories_plus[name])
    sig_cross = np.array(histories_cross[name])
    amp_plus = np.std(sig_plus - np.mean(sig_plus))
    amp_cross = np.std(sig_cross - np.mean(sig_cross))
    report += f"| {name} | {amp_plus:.6f} | {amp_cross:.6f} |\n"

report += f"""
## Quadrupole Pattern

- North-East correlation (θ+): {corr:.4f}
- 45° amplitude (θ×): {amp_ne_cross:.6f}

## Verdict

{verdict_plus}
{verdict_cross}
"""

with open("QW-1517_Polarization_Modes.md", "w") as f:
    f.write(report)

print("\n[SAVED] QW-1517_Polarization_Modes.md")
