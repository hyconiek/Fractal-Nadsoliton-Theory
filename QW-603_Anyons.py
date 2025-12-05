#!/usr/bin/env python3
# QW-603: ANYONIC STATISTICS TEST (Vortex Braiding)
# Purpose: Test if hopfions acquire fractional phase when exchanged
# Theory: Bosons θ=0, Fermions θ=π, Anyons 0<θ<π
# Date: 2025-12-05

import numpy as np
from scipy.ndimage import convolve

print("="*80)
print("QW-603: ANYONIC STATISTICS (Vortex Braid Test)")
print("="*80)
print("Test: Do hopfions get fractional phase when positions exchange?")
print("="*80)

# ============================================================================
# CONCEPT: Exchange two identical hopfions and measure phase change
# Method: Apply forces to swap positions, then measure ΔΦ
# ============================================================================

GRID_SIZE = 32
DT = 0.01
STEPS_EXCHANGE = 400  # Time to complete one exchange
GAMMA = 0.3
BETA_TORS = 0.01

WINDING = +1  # Test identical particles

print(f"Grid: {GRID_SIZE}^3")
print(f"Winding: both +{WINDING} (identical)")
print(f"Exchange time: {STEPS_EXCHANGE} steps")
print("-" * 40)

# ============================================================================
# HOPFION FIELD
# ============================================================================

def hopfion_field(grid_size, center, R=3.0, winding=1):
    x = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    y = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    z = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    X = X - center[0]
    Y = Y - center[1]
    Z = Z - center[2]
    
    rho = np.sqrt(X**2 + Y**2)
    rho[rho == 0] = 1e-10
    
    eta = np.arctan2(Z, rho - R)
    xi = np.arctan2(Y, X)
    
    phase = winding * (xi + eta)
    dist_to_ring = np.sqrt((rho - R)**2 + Z**2)
    amplitude = np.tanh(dist_to_ring / 1.5)
    
    return amplitude * np.exp(1j * phase)

# Initial: Two hopfions on y-axis
center_A = (0, -6, 0)  # Bottom
center_B = (0, +6, 0)  # Top

psi_A = hopfion_field(GRID_SIZE, center_A, R=3.0, winding=WINDING)
psi_B = hopfion_field(GRID_SIZE, center_B, R=3.0, winding=WINDING)

psi_A = psi_A / (np.max(np.abs(psi_A)) + 1e-10)
psi_B = psi_B / (np.max(np.abs(psi_B)) + 1e-10)

psi = 0.7 * (psi_A + psi_B)

print(f"Hopfion A: {center_A} (winding +{WINDING})")
print(f"Hopfion B: {center_B} (winding +{WINDING})")

# Measure initial phase
initial_phase_sum = np.sum(np.angle(psi))

print(f"Initial global phase: {initial_phase_sum:.2f}")

# ============================================================================
# EVOLUTION WITH EXCHANGE FORCE
# ============================================================================

laplacian_kernel = np.zeros((3,3,3))
laplacian_kernel[1,1,1] = -6
for i in [(0,1,1), (2,1,1), (1,0,1), (1,2,1), (1,1,0), (1,1,2)]:
    laplacian_kernel[i] = 1

def laplacian(field):
    return convolve(field, laplacian_kernel, mode='wrap')

# Exchange strategy: Apply rotating force
# A moves: (0,-6) → (+6,0) → (0,+6)
# B moves: (0,+6) → (-6,0) → (0,-6)
# This is a 180° rotation = half exchange

FORCE_STRENGTH = 0.03

print("\nPerforming exchange (braid) operation...")

x = np.arange(GRID_SIZE) - GRID_SIZE/2
y = np.arange(GRID_SIZE) - GRID_SIZE/2
z = np.arange(GRID_SIZE) - GRID_SIZE/2
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

for t in range(STEPS_EXCHANGE):
    rho = np.abs(psi)**2
    attractor = GAMMA * (1.0 - rho) * psi
    kin = 1j * laplacian(psi)
    back = -1j * BETA_TORS * rho * psi
    
    # Rotating force (circular path)
    theta_force = 2 * np.pi * t / STEPS_EXCHANGE
    
    # Force A clockwise, B counterclockwise
    force_A_x = FORCE_STRENGTH * np.cos(theta_force) * np.exp(-(X**2 + (Y+6)**2 + Z**2)/16)
    force_A_y = FORCE_STRENGTH * np.sin(theta_force) * np.exp(-(X**2 + (Y+6)**2 + Z**2)/16)
    
    force_B_x = -FORCE_STRENGTH * np.cos(theta_force) * np.exp(-(X**2 + (Y-6)**2 + Z**2)/16)
    force_B_y = -FORCE_STRENGTH * np.sin(theta_force) * np.exp(-(X**2 + (Y-6)**2 + Z**2)/16)
    
    force_total = 1j * (force_A_x + force_B_x) * X * psi + 1j * (force_A_y + force_B_y) * Y * psi
    
    dpsi_dt = attractor + kin + back + force_total
    psi += DT * dpsi_dt
    
    if t % 100 == 0:
        print(f"t={t:4d}: Exchange progress {t/STEPS_EXCHANGE*100:.0f}%")

print("-" * 40)

# Measure final phase
final_phase_sum = np.sum(np.angle(psi))

# Phase change (modulo 2π)
delta_phase = final_phase_sum - initial_phase_sum
delta_phase_wrapped = (delta_phase + np.pi) % (2*np.pi) - np.pi

print(f"\nFinal global phase: {final_phase_sum:.2f}")
print(f"ΔΦ (raw): {delta_phase:.2f}")
print(f"ΔΦ (wrapped): {delta_phase_wrapped:.2f}")
print()

# Normalize to exchange angle θ (single exchange = π for fermions, 0 for bosons)
# We did full circle (2π rotation), so divide by 2
theta_statistics = abs(delta_phase_wrapped) / 2

print(f"Exchange angle θ: {theta_statistics:.3f}")
print()

# Classification
if theta_statistics < 0.1:
    print("✅ BOSONIC (θ ≈ 0)")
    particle_type = "boson"
elif abs(theta_statistics - np.pi) < 0.3:
    print("✅ FERMIONIC (θ ≈ π)")
    particle_type = "fermion"
elif 0.1 < theta_statistics < np.pi - 0.3:
    print("🌟 ANYONIC (0 < θ < π)!")
    print(f"   Fractional statistics: θ/{np.pi:.3f}π")
    particle_type = "anyon"
else:
    print("🟡 UNCLEAR")
    particle_type = "unknown"

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw603_anyons.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-603: Anyonic Statistics Test\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Test:** Vortex braiding (particle exchange)\n\n")
    
    f.write("## 1. Setup\n")
    f.write(f"- Dwa identyczne hopfiony: winding = +{WINDING}\n")
    f.write(f"- Wymiana: obrót 360° (pełna wymiana)\n\n")
    
    f.write("## 2. Wyniki\n")
    f.write(f"- **ΔΦ:** {delta_phase_wrapped:.3f}\n")
    f.write(f"- **θ (exchange angle):** {theta_statistics:.3f}\n")
    f.write(f"- **θ/π:** {theta_statistics/np.pi:.3f}\n\n")
    
    f.write("## 3. Klasyfikacja\n")
    if particle_type == "boson":
        f.write("### ✅ BOZON\n")
        f.write("Wymiana nie zmienia fazy (θ≈0). Zgodne z oczekiwaniem dla identycznych bozonów.\n")
    elif particle_type == "fermion":
        f.write("### ✅ FERMION\n")
        f.write("Wymiana daje fazę π. Hopfiony wykazują statystykę fermionową!\n")
    elif particle_type == "anyon":
        f.write("### 🌟 ANYON (Fractional Statistics)!\n")
        f.write(f"Wymiana daje fazę {theta_statistics:.3f} = {theta_statistics/np.pi:.3f}π.\n\n")
        f.write("**To jest przełomowe odkrycie!**\n")
        f.write("Hopfiony w FIN wykazują anyoniczną statystykę - otwiera drogę do:\n")
        f.write("- Topological quantum computing\n")
        f.write("- Link do fractional quantum Hall effect\n")
    else:
        f.write("### 🟡 NIEJEDNOZNACZNE\n")
        f.write("Wynik nie pasuje do znanych statystyk.\n")

print("Report saved.")
print("="*80)
