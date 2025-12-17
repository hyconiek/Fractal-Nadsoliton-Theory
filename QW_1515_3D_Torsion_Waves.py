#!/usr/bin/env python3
"""
QW-1515: 3D TORSION WAVE SIMULATION
====================================
Rigorous 3D extension of QW-1513.

PHYSICS:
In 3D, the torsion wave equation is:
∂²θ/∂t² = c² ∇²θ - γ ∂θ/∂t + S(r,t)

where ∇² is the 3D Laplacian.

For a spherically symmetric source, we expect:
- Waves propagating radially
- Amplitude decay as 1/r (not 1/r² as for static field)
- This is the CRITICAL test for gravitational waves

METHODOLOGY:
1. Create 3D cubic lattice
2. Place oscillating source at center
3. Measure wave amplitude at various distances
4. Fit amplitude vs distance to A ~ 1/r^n
5. Verify n = 1 (GR prediction for GW)
"""

import numpy as np
from scipy.fft import fftn, fftfreq
from scipy.ndimage import laplace
import datetime

print("="*80)
print("QW-1515: 3D TORSION WAVE SIMULATION")
print("="*80)

# Frozen parameters
ALPHA_GEO = 4 * np.log(2)  # 2.7726
OMEGA_FUND = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01

def K(d):
    """Coupling kernel"""
    if d < 0.1:
        d = 0.1
    return ALPHA_GEO * np.cos(OMEGA_FUND * d + PHI) / (1 + BETA_TORS * d)

# Wave speed from K(0)
K_0 = K(0)
c_tors = np.sqrt(abs(K_0))
print(f"Wave speed: c_tors = sqrt(K(0)) = sqrt({K_0:.4f}) = {c_tors:.4f}")

# 3D Grid parameters
N = 64  # Grid size (N x N x N)
L = 20.0  # Physical size
dx = L / N
print(f"\n3D Grid: {N}³ = {N**3} nodes")
print(f"Physical size: L = {L}, dx = {dx:.4f}")

# Source parameters
source_center = np.array([N//2, N//2, N//2])
source_freq = 0.1  # Orbital frequency
gw_freq = 2 * source_freq  # GW frequency = 2 × orbital
source_amplitude = 1.0

print(f"\nSource at center: {source_center}")
print(f"GW frequency: f_gw = 2 × f_orbit = {gw_freq}")

# Initialize fields
theta = np.zeros((N, N, N))  # Phase field
theta_dot = np.zeros((N, N, N))  # Time derivative

# Time evolution
dt = 0.1
t_max = 100.0
n_steps = int(t_max / dt)

print(f"\nTime evolution: dt = {dt}, t_max = {t_max}, steps = {n_steps}")

# Detector positions (radial distances from center)
detector_radii = [3, 5, 7, 9, 11]  # In grid units
detector_histories = {r: [] for r in detector_radii}

print("\n[RUNNING 3D SIMULATION]")
print("-" * 60)

# Time stepping (leapfrog/Verlet method)
for step in range(n_steps):
    t = step * dt
    
    # 1. Compute Laplacian (3D discrete)
    laplacian = laplace(theta) / (dx**2)
    
    # 2. Source term (oscillating at GW frequency)
    source = np.zeros((N, N, N))
    source[source_center[0], source_center[1], source_center[2]] = (
        source_amplitude * np.sin(2 * np.pi * gw_freq * t)
    )
    
    # 3. Wave equation: d²θ/dt² = c² ∇²θ - γ dθ/dt + S
    theta_ddot = c_tors**2 * laplacian - BETA_TORS * theta_dot + source
    
    # 4. Update (Verlet integration)
    theta_dot += theta_ddot * dt
    theta += theta_dot * dt
    
    # 5. Record at detectors
    for r in detector_radii:
        # Average over spherical shell at radius r
        i, j, k = np.meshgrid(range(N), range(N), range(N), indexing='ij')
        dist = np.sqrt((i - source_center[0])**2 + 
                       (j - source_center[1])**2 + 
                       (k - source_center[2])**2)
        shell_mask = (dist >= r - 0.5) & (dist < r + 0.5)
        
        if np.sum(shell_mask) > 0:
            avg_theta = np.mean(theta[shell_mask])
        else:
            avg_theta = 0.0
        
        detector_histories[r].append(avg_theta)
    
    # Progress
    if step % (n_steps // 10) == 0:
        print(f"  Step {step}/{n_steps}, t = {t:.1f}")

print("\n[ANALYSIS]")
print("-" * 60)

# Measure wave amplitude at each detector
amplitudes = []
for r in detector_radii:
    signal = np.array(detector_histories[r])
    # Remove DC and measure amplitude
    signal_ac = signal - np.mean(signal)
    amplitude = np.std(signal_ac)
    amplitudes.append(amplitude)
    print(f"r = {r}: Amplitude = {amplitude:.6f}")

amplitudes = np.array(amplitudes)
radii = np.array(detector_radii, dtype=float)

# Fit amplitude vs radius: A = A0 / r^n
# log(A) = log(A0) - n × log(r)
valid = amplitudes > 1e-10
if np.sum(valid) >= 3:
    log_r = np.log(radii[valid])
    log_A = np.log(amplitudes[valid])
    
    # Linear fit
    coeffs = np.polyfit(log_r, log_A, 1)
    n_fit = -coeffs[0]  # Exponent
    A0_fit = np.exp(coeffs[1])
    
    print(f"\nAmplitude scaling: A = {A0_fit:.4f} / r^{n_fit:.2f}")
    print(f"Expected for GW: n = 1.0")
    print(f"Error: |n - 1| = {abs(n_fit - 1.0):.2f}")
else:
    n_fit = 0
    print("\nInsufficient data for fit")

# VERDICT
print("\n" + "=" * 60)
print("QW-1515 VERDICT")
print("=" * 60)

if 0.8 <= n_fit <= 1.2:
    verdict = "✅ SUCCESS: Amplitude scales as 1/r (consistent with GW)"
elif 0.5 <= n_fit <= 1.5:
    verdict = f"🟡 PARTIAL: n = {n_fit:.2f} (close to 1.0)"
else:
    verdict = f"❌ FAILURE: n = {n_fit:.2f} (expected 1.0)"

print(verdict)

# Save report
report = f"""# QW-1515: 3D Torsion Wave Simulation

**Date:** {datetime.datetime.now()}
**Status:** {'SUCCESS' if 0.8 <= n_fit <= 1.2 else 'PARTIAL' if 0.5 <= n_fit <= 1.5 else 'FAILURE'}

## Parameters
- Grid: {N}³ = {N**3} nodes
- Wave speed: c_tors = {c_tors:.4f}
- Source frequency: f_gw = {gw_freq}

## Results

| Radius | Amplitude |
|--------|-----------|
"""
for r, A in zip(detector_radii, amplitudes):
    report += f"| {r} | {A:.6f} |\n"

report += f"""
## Amplitude Scaling
$$A(r) = \\frac{{{A0_fit:.4f}}}{{r^{{{n_fit:.2f}}}}}$$

**Expected (GR):** n = 1.0
**Measured:** n = {n_fit:.2f}
**Error:** {abs(n_fit - 1.0):.2f}

## Verdict
{verdict}
"""

with open("QW-1515_3D_Torsion_Waves.md", "w") as f:
    f.write(report)

print(f"\n[SAVED] QW-1515_3D_Torsion_Waves.md")
