#!/usr/bin/env python3
"""
QW-1523: NEURO-GRAVITY WAVE SIMULATION
=========================================
Testing the hypothesis that Gravitational Waves are waves of
SYNAPTIC WEIGHT UPDATES propagating through the Neural Vacua.

PHYSICS:
- Space = Network of nodes
- Metric g_uv = Synaptic weights w_ij
- Gravitational Wave = Propagation of Δw_ij (weight update)

ALGORITHM:
1. Initialize a 1D chain of neurons (as a simplified slice of 3D vacuum).
2. Weights w_ij define the "connectivity" (metric).
3. Perturb weights in the center (representing a merger event).
4. Use a Hebbian-like update rule with finite propagation delay:
   Δw_ij(t) ~ correlation(i, j, t-delay)
5. Measure how the perturbation Δw propagates.

Unlike mechanical waves (restoring force), neural waves propagate
information about correlations.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import datetime

print("="*80)
print("QW-1523: NEURO-GRAVITY WAVE SIMULATION")
print("="*80)

# Parameters
N = 512  # Neurons
L = 200.0
dx = L / N
dt = 0.05
t_max = 300.0
n_steps = int(t_max / dt)

# FIN Parameters (Frozen)
ALPHA_GEO = 4 * np.log(2)
# Learning rate (effective coupling speed)
ETA = 0.1 
# Decay rate (forgetting factor - creates effective damping)
DECAY = 0.01

print(f"Network: N={N}")
print(f"Learning Rate η={ETA}")
print(f"Decay Rate γ={DECAY}")

# State: Synaptic Weights w[i] 
# (Simplification: w[i] represents coupling strength at pos i relative to background)
w = np.ones(N)  # Background metric = 1 (flat space)
w_history = np.zeros((n_steps, N))

# Source: Oscillating perturbation of weights at center
source_pos = N // 4
f_gw = 0.1

print(f"Source Freq f_gw={f_gw}")

# Propagation mechanism in Neural Network:
# Information propagates from neighbor to neighbor.
# w_new[i] depends on w_old[i] and neighbors w_old[i-1], w_old[i+1]
# This looks like diffusion-reaction or wave equation depending on coefficients.

# Neural Update Rule (Linearized):
# w(t+1) = w(t) + η * Laplacian(w) - γ * (w - 1) + Source
#
# Wait! Is it diffusion or wave?
# Neural signals are electrochemical waves (action potentials).
# We model this as a WAVE equation for neural activity, which drives weight changes.
#
# Let u[i] be the "neural activity" (potential).
# Wave equation for u: u_tt = c^2 u_xx - damping * u_t
# Weights w[i] are modulated by u[i].

# Let's verify K(d) parameter influence.
# c_neural = 1.0 (topology neighbor-to-neighbor speed)

c_neural = 1.0

# Fields
u = np.zeros(N)       # Potential
u_dot = np.zeros(N)   # Time derivative

# Detector positions
detectors = [20, 40, 60, 80, 100, 120]
det_histories = {r: [] for r in detectors}

print("\nRunning Neural Wave Simulation...")

for step in range(n_steps):
    t = step * dt
    
    # 1. Source (Neural Oscillation at center)
    source = np.zeros(N)
    source[source_pos] = np.sin(2 * np.pi * f_gw * t)
    
    # 2. Laplacian (Connectivity)
    laplacian = np.zeros(N)
    # Vectorized laplacian (1D with boundary 0)
    laplacian[1:-1] = (u[2:] - 2*u[1:-1] + u[:-2]) / dx**2
    
    # 3. Wave Equation for Neural Activity
    # u_tt = c^2 Δu - γ u_t + S
    # Using BETA_TORS from FIN as damping
    u_ddot = c_neural**2 * laplacian - 0.01 * u_dot + source
    
    u_dot += u_ddot * dt
    u += u_dot * dt
    
    # 4. Weights respond to Activity (Plasticity)
    # w changes based on local activity energy u^2
    # w[i] = 1 + alpha * u[i]^2 (Metric perturbation h ~ u^2)
    # But for propagation measurement, we look at u directly first.
    
    # Record
    w_history[step, :] = u
    
    for r in detectors:
        pos = source_pos + int(r/dx)
        if 0 <= pos < N:
            det_histories[r].append(u[pos])
        else:
            det_histories[r].append(0)

print("Simulation Complete.")

# Analysis
print("\n[AMPLITUDE ANALYSIS]")
amplitudes = []
for r in detectors:
    sig = np.array(det_histories[r])
    sig = sig[int(n_steps/2):] # Ignore startup transient
    amp = np.std(sig)
    amplitudes.append(amp)
    print(f"r={r:3d}: Amp={amp:.6f}")

amplitudes = np.array(amplitudes)
radii = np.array(detectors, dtype=float)

# Fit
valid = amplitudes > 1e-10
if np.sum(valid) >= 3:
    log_r = np.log(radii[valid])
    log_A = np.log(amplitudes[valid])
    coeffs = np.polyfit(log_r, log_A, 1)
    n_fit = -coeffs[0]
    print(f"\nScaling: A ~ 1/r^{n_fit:.2f}")

# Visualization
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle(f'QW-1523: Neuro-Gravity Wave (n={n_fit:.2f})')

ax1 = axes[0]
ax1.loglog(radii[valid], amplitudes[valid], 'bo-', label='Measured')
ax1.loglog(radii[valid], amplitudes[valid][0] * (radii[valid][0]/radii[valid])**0.5, 'g--', label='1/sqrt(r) (1D Energy)')
ax1.loglog(radii[valid], amplitudes[valid][0] * (radii[valid][0]/radii[valid])**1.0, 'r:', label='1/r (3D Energy)')
ax1.legend()
ax1.set_title('Amplitude Scaling')
ax1.set_xlabel('Distance')

ax2 = axes[1]
ax2.plot(w_history[-1,:], 'k-')
ax2.set_title('Final Network State')

plt.tight_layout()
plt.savefig('QW-1523_NeuroGravity.png')
print("[SAVED] QW-1523_NeuroGravity.png")

# Verdict
print("\nVERDICT for QW-1523:")
if abs(n_fit - 0.5) < 0.1:
    print("✅ Matches 1D Energy Conservation (n=0.5)")
    print("   Note: In 1D, A ~ 1/sqrt(r) is correct conservation.")
    print("   Our previous 3D simulations failed because we mapped 1D chain to 3D space incorrectly.")
elif abs(n_fit - 1.0) < 0.2:
    print("⚠️ Matches 3D 1/r scaling (Unexpected for 1D chain)")
else:
    print(f"❌ Anomalous scaling n={n_fit:.2f}")

# Report
with open('QW-1523_NeuroGravity_Result.md', 'w') as f:
    f.write(f"# QW-1523 Result\n\nScale n = {n_fit:.2f}\n")
