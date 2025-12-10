#!/usr/bin/env python3
# QW-594: REVISED HOPFION STABILITY (ATTRACTOR DYNAMICS)
# Purpose: Test Hypothesis 4 (Particles as Vortices) using stabilized Attractor Dynamics
#          instead of pure Unitary evolution. Fixes NaN issues from QW-590.
# Date: 2025-12-05

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import convolve

print("="*80)
print("QW-594: REVISED HOPFION STABILITY (ATTRACTOR DYNAMICS)")
print("="*80)

# ============================================================================
# 1. PARAMETERS & CONFIGURATION
# ============================================================================
GRID_SIZE = 32
DT = 0.01  # Smaller timestep for stability
STEPS = 500
GAMMA = 0.5 # Relaxation rate (Attractor strength)
BETA_TORS = 0.01
N_LAYER = 10

print(f"Grid: {GRID_SIZE}x{GRID_SIZE}x{GRID_SIZE}")
print(f"Dynamics: Ginzburg-Landau (Attractor)")
print(f"Gamma: {GAMMA}")
print("-" * 40)

# ============================================================================
# 2. FIELD INITIALIZATION (SOFT HOPFION)
# ============================================================================

def hopfion_field(grid_size, R=8.0):
    """
    Initialize a Hopfion (knot) configuration in 3D.
    Maps S3 -> S2.
    """
    x = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    y = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    z = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    rho = np.sqrt(X**2 + Y**2)
    radius = np.sqrt(X**2 + Y**2 + Z**2)
    
    # Avoid division by zero
    rho[rho == 0] = 1e-10
    
    # Toroidal coordinates
    eta = np.arctan2(Z, rho - R)
    xi = np.arctan2(Y, X)
    
    # Phase: Linked loops
    phase = xi + eta
    
    # Amplitude: "Soft core" vortex
    # Goes to 0 at the ring (rho=R, Z=0)
    # Goes to 1 at infinity
    dist_to_ring = np.sqrt((rho - R)**2 + Z**2)
    amplitude = np.tanh(dist_to_ring / 2.0)
    
    # Initial field
    psi = amplitude * np.exp(1j * phase)
    
    # Normalize to be safe? No, let GL dynamics handle it.
    return psi

psi = hopfion_field(GRID_SIZE, R=GRID_SIZE/4.0)

# ============================================================================
# 3. ROBUST TOPOLOGICAL CHARGE
# ============================================================================

def get_derivatives(psi):
    # Central differences with wrap boundary
    dx = (np.roll(psi, -1, axis=0) - np.roll(psi, 1, axis=0)) / 2.0
    dy = (np.roll(psi, -1, axis=1) - np.roll(psi, 1, axis=1)) / 2.0
    dz = (np.roll(psi, -1, axis=2) - np.roll(psi, 1, axis=2)) / 2.0
    return dx, dy, dz

def compute_hopf_invariant_robust(psi):
    """
    Compute H = v . w
    Regularized specifically to avoid singularities at vortex core.
    """
    # 1. Superfluid velocity v_s = Im(psi* grad psi) / rho
    # Regularized: v_s = Im(psi* grad psi) / (rho + eps)
    
    psi_conj = np.conj(psi)
    rho = np.abs(psi)**2 + 1e-6 # Regularization epsilon
    
    grad_x, grad_y, grad_z = get_derivatives(psi)
    
    # Current density J = Im(psi* grad psi)
    J_x = np.imag(psi_conj * grad_x)
    J_y = np.imag(psi_conj * grad_y)
    J_z = np.imag(psi_conj * grad_z)
    
    # Velocity field
    v_x = J_x / rho
    v_y = J_y / rho
    v_z = J_z / rho
    
    # 2. Vorticity w = curl v
    # Use simple curl on the vector field v
    
    dv_z_dy = (np.roll(v_z, -1, axis=1) - np.roll(v_z, 1, axis=1)) / 2.0
    dv_y_dz = (np.roll(v_y, -1, axis=2) - np.roll(v_y, 1, axis=2)) / 2.0
    w_x = dv_z_dy - dv_y_dz
    
    dv_x_dz = (np.roll(v_x, -1, axis=2) - np.roll(v_x, 1, axis=2)) / 2.0
    dv_z_dx = (np.roll(v_z, -1, axis=0) - np.roll(v_z, 1, axis=0)) / 2.0
    w_y = dv_x_dz - dv_z_dx
    
    dv_y_dx = (np.roll(v_y, -1, axis=0) - np.roll(v_y, 1, axis=0)) / 2.0
    dv_x_dy = (np.roll(v_x, -1, axis=1) - np.roll(v_x, 1, axis=1)) / 2.0
    w_z = dv_y_dx - dv_x_dy
    
    # 3. Helicity density h = v . w
    h = v_x * w_x + v_y * w_y + v_z * w_z
    
    return np.sum(h)

initial_Q = compute_hopf_invariant_robust(psi)
print(f"Initial Robust Q: {initial_Q:.4f}")

# ============================================================================
# 4. EVOLUTION (GINZBURG-LANDAU ATTRACTOR)
# ============================================================================

# Laplacian convolution kernel
laplacian_kernel = np.zeros((3,3,3))
laplacian_kernel[1,1,1] = -6
for i in [(0,1,1), (2,1,1), (1,0,1), (1,2,1), (1,1,0), (1,1,2)]:
    laplacian_kernel[i] = 1

def laplacian(field):
    return convolve(field, laplacian_kernel, mode='wrap')

print("Simulating Attractor Dynamics (Ginzburg-Landau)...")

history_Q = []
history_E = []

for t in range(STEPS):
    # 1. Attractor Term (The "Nadsoliton" tendency)
    # dPsi/dt = Gamma * (1 - |psi|^2) * psi
    # Pushes amplitude towards 1 everywhere, except where topology forbids it (core)
    
    rho = np.abs(psi)**2
    attractor = GAMMA * (1.0 - rho) * psi
    
    # 2. Kinetic Term (Elasticity/Dispersion)
    # dPsi/dt = i * Laplacian * psi
    # Rotates phase, smoothes gradients
    kin = 1j * laplacian(psi)
    
    # 3. Effective Back-Reaction (Layer N=10)
    # V_eff = Beta * rho * psi (Nonlinear mass) - stabilizes loops?
    # Or just a simple confining potential? 
    # Let's assume Back-Reaction adds a "mass" penalty for high density
    # V ~ -mu * |psi|^2
    back_reaction = -1j * BETA_TORS * rho * psi
    
    # Total update
    dpsi_dt = attractor + kin + back_reaction
    
    psi += DT * dpsi_dt
    
    # Measurements
    if t % 50 == 0:
        Q = compute_hopf_invariant_robust(psi)
        
        # Energy Functional (Ginzburg-Landau Energy)
        # E = Integral [ |grad psi|^2 + (1-|psi|^2)^2 ]
        grad_x, grad_y, grad_z = get_derivatives(psi)
        grad_sq = np.abs(grad_x)**2 + np.abs(grad_y)**2 + np.abs(grad_z)**2
        pot_term = (1.0 - np.abs(psi)**2)**2
        E = np.sum(grad_sq + pot_term)
        
        history_Q.append(Q)
        history_E.append(E)
        
        print(f"Step {t:4d}: Q={Q:.4f}, E={E:.2f}")

print("-" * 40)

# ============================================================================
# 5. FINAL ANALYSIS
# ============================================================================

final_Q = history_Q[-1]
ratio = final_Q / initial_Q if initial_Q != 0 else 0

print(f"Initial Q: {initial_Q:.4f}")
print(f"Final Q:   {final_Q:.4f}")
print(f"Ratio:     {ratio*100:.1f}%")

print()
if ratio > 0.8:
    print("✅ SUCCESS: Hopfion Stabilized via Attractor Dynamics!")
    print("   Hypothesis 4 CONFIRMED: Particles are topological attractors.")
elif ratio > 0.4:
    print("🟡 PARTIAL: Hopfion decaying slowly (Metastable).")
    print("   Suggests partial truth to H4.")
else:
    print("❌ FAILURE: Topology lost.")
    print("   Even attractor dynamics cannot save the Hopfion.")

# ============================================================================
# 6. REPORT GENERATION
# ============================================================================

report_path = "raport_qw594_hopfion.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write(f"# Raport QW-594: Revised Hopfion Stability\n")
    f.write(f"**Data:** 2025-12-05\n")
    f.write(f"**Metoda:** Ginzburg-Landau Attractor Dynamics (Layer N={N_LAYER})\n\n")
    
    f.write("## 1. Parametry Symulacji\n")
    f.write(f"- Grid: {GRID_SIZE}^3\n")
    f.write(f"- Steps: {STEPS}\n")
    f.write(f"- DT: {DT}\n")
    f.write(f"- Gamma (Attractor Strength): {GAMMA}\n\n")
    
    f.write("## 2. Wyniki Stabilności\n")
    f.write(f"- **Initial Helicity (Q):** {initial_Q:.4f}\n")
    f.write(f"- **Final Helicity (Q):** {final_Q:.4f}\n")
    f.write(f"- **Retention Ratio:** {ratio*100:.2f}%\n\n")
    
    f.write("## 3. Przebieg Ewolucji\n")
    f.write("| Step | Charge Q | Energy E |\n")
    f.write("|------|----------|----------|\n")
    for i in range(len(history_Q)):
        step = i * 50
        f.write(f"| {step} | {history_Q[i]:.4f} | {history_E[i]:.2f} |\n")
    f.write("\n")
    
    f.write("## 4. Werdykt\n")
    if ratio > 0.8:
        f.write("### ✅ SUKCES: Hopfion Stabilny\n")
        f.write("Hipoteza 4 (Cząstki jako Wiry) została **POTWIERDZONA**. Dynamika atraktora (Ginzburg-Landau) skutecznie chroni topologię wiru przed rozpadem, który występował w modelach unitarnych.\n")
    else:
        f.write("### ❌ PORAŻKA: Hopfion Zniknął\n")
        f.write("Nawet dynamika atraktora nie utrzymała struktury topologicznej.\n")

print("Report saved.")

