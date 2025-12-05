#!/usr/bin/env python3
# QW-605b: PHASE DECOHERENCE DIAGNOSIS
# Purpose: Understand why phase decreases with multiple braids
# Hypothesis: Attractor dynamics (gamma) causes phase decay
# Date: 2025-12-05

import numpy as np
from scipy.ndimage import convolve

print("="*80)
print("QW-605b: PHASE DECOHERENCE DIAGNOSIS")
print("="*80)
print("Test: Does gamma (attractor) cause phase decay?")
print("="*80)

GRID_SIZE = 32
DT = 0.01
STEPS_PER_BRAID = 400
BETA_TORS = 0.01

# Test with different gamma values
GAMMA_VALUES = [0.0, 0.1, 0.3]  # 0.3 was used in QW-605

print(f"Testing gamma values: {GAMMA_VALUES}")
print("-" * 40)

def hopfion_field(grid_size, center, R=3.0, winding=1):
    x = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    y = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    z = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    X, Y, Z = X - center[0], Y - center[1], Z - center[2]
    rho = np.sqrt(X**2 + Y**2)
    rho[rho == 0] = 1e-10
    eta = np.arctan2(Z, rho - R)
    xi = np.arctan2(Y, X)
    phase = winding * (xi + eta)
    dist_to_ring = np.sqrt((rho - R)**2 + Z**2)
    amplitude = np.tanh(dist_to_ring / 1.5)
    return amplitude * np.exp(1j * phase)

laplacian_kernel = np.zeros((3,3,3))
laplacian_kernel[1,1,1] = -6
for i in [(0,1,1), (2,1,1), (1,0,1), (1,2,1), (1,1,0), (1,1,2)]:
    laplacian_kernel[i] = 1

def laplacian(field):
    return convolve(field, laplacian_kernel, mode='wrap')

results = []

for gamma in GAMMA_VALUES:
    print(f"\nGamma: {gamma}")
    
    # Initialize
    psi_A = hopfion_field(GRID_SIZE, (0, -6, 0), R=3.0, winding=+1)
    psi_B = hopfion_field(GRID_SIZE, (0, +6, 0), R=3.0, winding=+1)
    psi = 0.7 * (psi_A/3 + psi_B/3)
    
    initial_phase = np.sum(np.angle(psi))
    
    # Single braid
    x = np.arange(GRID_SIZE) - GRID_SIZE/2
    y = np.arange(GRID_SIZE) - GRID_SIZE/2
    z = np.arange(GRID_SIZE) - GRID_SIZE/2
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    for t in range(STEPS_PER_BRAID):
        rho = np.abs(psi)**2
        attractor = gamma * (1.0 - rho) * psi
        kin = 1j * laplacian(psi)
        back = -1j * BETA_TORS * rho * psi
        
        theta_force = 2 * np.pi * t / STEPS_PER_BRAID
        force_strength = 0.03
        
        force_A_x = force_strength * np.cos(theta_force) * np.exp(-(X**2 + (Y+6)**2 + Z**2)/16)
        force_A_y = force_strength * np.sin(theta_force) * np.exp(-(X**2 + (Y+6)**2 + Z**2)/16)
        force_B_x = -force_strength * np.cos(theta_force) * np.exp(-(X**2 + (Y-6)**2 + Z**2)/16)
        force_B_y = -force_strength * np.sin(theta_force) * np.exp(-(X**2 + (Y-6)**2 + Z**2)/16)
        
        force_total = 1j * (force_A_x + force_B_x) * X * psi + 1j * (force_A_y + force_B_y) * Y * psi
        
        dpsi_dt = attractor + kin + back + force_total
        psi += DT * dpsi_dt
    
    final_phase = np.sum(np.angle(psi))
    delta_phase = ((final_phase - initial_phase) + np.pi) % (2*np.pi) - np.pi
    theta = abs(delta_phase) / 2
    
    results.append({'gamma': gamma, 'theta': theta})
    print(f"  θ = {theta:.3f}")

print("\n" + "="*80)
print("DIAGNOSIS")
print("="*80)

print("\n| Gamma | θ (single braid) |")
print("|-------|------------------|")
for r in results:
    print(f"| {r['gamma']:.1f}   | {r['theta']:.3f}            |")

# Check if theta changes with gamma
theta_values = [r['theta'] for r in results]
theta_std = np.std(theta_values)

print(f"\nTheta standard deviation: {theta_std:.3f}")

if theta_std < 0.1:
    print("\n✅ Gamma NIE wpływa na θ")
    print("   Phase decay w QW-605 nie jest przez attractor")
else:
    print("\n⚠️ Gamma WPŁYWA na θ")
    print("   Attractor może powodować phase decay")

with open("raport_qw605b_diagnosis.md", "w") as f:
    f.write("# QW-605b: Phase Decoherence Diagnosis\n\n")
    f.write("**Test:** Czy gamma (attractor) powoduje phase decay?\n\n")
    f.write("| Gamma | θ |\n|-------|---|\n")
    for r in results:
        f.write(f"| {r['gamma']:.1f} | {r['theta']:.3f} |\n")
    f.write(f"\n**Std Dev:** {theta_std:.3f}\n\n")
    if theta_std < 0.1:
        f.write("**Wniosek:** Gamma nie jest przyczyną phase decay.\n")
    else:
        f.write("**Wniosek:** Gamma wpływa na fazę - może być przyczyną decay.\n")

print("="*80)
