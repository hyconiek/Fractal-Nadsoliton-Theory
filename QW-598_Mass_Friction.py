#!/usr/bin/env python3
# QW-598: MASS AS VACUUM FRICTION
# Purpose: Test if hopfion "mass" (inertia) comes from friction with vacuum field
# Hypothesis H5: Masa = Opór Eteru
# Date: 2025-12-05

import numpy as np
from scipy.ndimage import convolve

print("="*80)
print("QW-598: MASS AS VACUUM FRICTION (Test H5)")
print("="*80)

# ============================================================================
# CONCEPT: Apply external "force" to hopfion and measure acceleration
# If m_eff = F/a, and m_eff depends on vacuum coupling beta_tors,
# then we confirm mass = vacuum friction
# ============================================================================

GRID_SIZE = 32
DT = 0.01
STEPS = 400
GAMMA = 0.3

# Test different vacuum coupling strengths
BETA_VALUES = [0.005, 0.01, 0.02, 0.03]

print(f"Grid: {GRID_SIZE}^3")
print(f"Testing {len(BETA_VALUES)} vacuum coupling strengths")
print("-" * 40)

# ============================================================================
# HOPFION INITIALIZATION
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

# ============================================================================
# CORE TRACKING
# ============================================================================

def find_core_position(psi):
    density = np.abs(psi)**2
    threshold = 0.15 * np.mean(density)
    core_mask = density < threshold
    
    x = np.arange(GRID_SIZE) - GRID_SIZE/2
    y = np.arange(GRID_SIZE) - GRID_SIZE/2
    z = np.arange(GRID_SIZE) - GRID_SIZE/2
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    if np.sum(core_mask) > 10:
        x_c = np.sum(X[core_mask]) / np.sum(core_mask)
        y_c = np.sum(Y[core_mask]) / np.sum(core_mask)
        z_c = np.sum(Z[core_mask]) / np.sum(core_mask)
        return np.array([x_c, y_c, z_c])
    else:
        return np.array([0, 0, 0])

# ============================================================================
# EVOLUTION WITH EXTERNAL FORCE
# ============================================================================

laplacian_kernel = np.zeros((3,3,3))
laplacian_kernel[1,1,1] = -6
for i in [(0,1,1), (2,1,1), (1,0,1), (1,2,1), (1,1,0), (1,1,2)]:
    laplacian_kernel[i] = 1

def laplacian(field):
    return convolve(field, laplacian_kernel, mode='wrap')

results = []

print("\nTesting effective mass at different vacuum couplings...")

for beta_tors in BETA_VALUES:
    print(f"\nBeta_tors: {beta_tors:.3f}")
    
    # Initialize hopfion at center
    psi = hopfion_field(GRID_SIZE, center=(0, 0, 0), R=3.0, winding=+1)
    
    # Apply constant "force" in +x direction
    # This is done by adding a phase gradient: F ~ grad(phase)
    FORCE_STRENGTH = 0.05
    
    positions = []
    velocities = []
    
    for t in range(STEPS):
        rho = np.abs(psi)**2
        attractor = GAMMA * (1.0 - rho) * psi
        kin = 1j * laplacian(psi)
        back = -1j * beta_tors * rho * psi
        
        # External force: phase gradient in x direction
        x = np.arange(GRID_SIZE) - GRID_SIZE/2
        y = np.arange(GRID_SIZE) - GRID_SIZE/2
        z = np.arange(GRID_SIZE) - GRID_SIZE/2
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        
        force = 1j * FORCE_STRENGTH * X * psi  # Force proportional to x
        
        dpsi_dt = attractor + kin + back + force
        psi += DT * dpsi_dt
        
        # Track position
        if t % 20 == 0:
            pos = find_core_position(psi)
            positions.append(pos[0])  # x-coordinate only
    
    # Compute acceleration from position data
    if len(positions) > 3:
        # Fit x(t) = x0 + v0*t + 0.5*a*t^2
        times = np.arange(len(positions)) * 20 * DT
        
        # Simple finite difference for acceleration
        velocities_calc = np.diff(positions) / (20 * DT)
        if len(velocities_calc) > 1:
            acceleration = np.mean(np.diff(velocities_calc) / (20 * DT))
        else:
            acceleration = 0
        
        # Effective mass: m_eff = F / a
        if abs(acceleration) > 1e-6:
            m_eff = FORCE_STRENGTH / abs(acceleration)
        else:
            m_eff = float('inf')
        
        results.append({
            'beta': beta_tors,
            'acceleration': acceleration,
            'm_eff': m_eff,
            'final_pos': positions[-1]
        })
        
        print(f"  Acceleration: {acceleration:.4f}")
        print(f"  Effective mass: {m_eff:.2f}")
        print(f"  Final position: {positions[-1]:.2f}")

print("-" * 40)

# ============================================================================
# ANALYSIS
# ============================================================================

print("\nMass-Friction Relationship:")
print("| Beta (friction) | m_eff | Ratio m/β |")
print("|-----------------|-------|-----------|")
for r in results:
    ratio = r['m_eff'] / r['beta'] if r['beta'] > 0 else 0
    print(f"| {r['beta']:15.3f} | {r['m_eff']:5.1f} | {ratio:9.1f} |")

# Check if m_eff ∝ beta (linear relationship)
betas = [r['beta'] for r in results]
masses = [r['m_eff'] for r in results if r['m_eff'] < 1000]  # Filter infinities

if len(masses) >= 2:
    # Linear fit
    correlation = np.corrcoef(betas[:len(masses)], masses)[0, 1]
    
    print(f"\nCorrelation (m_eff vs beta): {correlation:.3f}")
    
    if correlation > 0.8:
        print("\n✅ H5 POTWIERDZONE:")
        print("   Masa hopfiona ROŚNIE z tarciem próżni (beta)!")
        print("   m_eff ∝ beta_tors")
        confirmed = True
    else:
        print("\n🟡 SŁABA KORELACJA")
        confirmed = False
else:
    confirmed = False

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw598_mass_friction.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-598: Mass as Vacuum Friction\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Hipoteza:** H5 - Masa to Opór Eteru\n\n")
    
    f.write("## 1. Koncepcja\n")
    f.write("Zastosowano zewnętrzną 'siłę' do hopfiona i zmierzono przyspieszenie.\n")
    f.write("Jeśli masa efektywna m_eff = F/a zależy od sprzężenia z próżnią (beta_tors),\n")
    f.write("to potwierdza że masa = opór (tarcie) próżni.\n\n")
    
    f.write("## 2. Wyniki\n")
    f.write("| Beta (tarcie) | Przyspieszenie | m_eff |\n")
    f.write("|---------------|----------------|-------|\n")
    for r in results:
        f.write(f"| {r['beta']:.3f} | {r['acceleration']:.4f} | {r['m_eff']:.1f} |\n")
    f.write("\n")
    
    if confirmed:
        f.write("## 3. Wnioski\n")
        f.write("### ✅ HIPOTEZA H5 POTWIERDZONA\n")
        f.write(f"Korelacja m_eff vs beta: {correlation:.3f} > 0.8\n\n")
        f.write("**Masa hopfiona jest proporcjonalna do sprzężenia z próżnią!**\n\n")
        f.write("To oznacza że:\n")
        f.write("- Masa nie jest wewnętrzną właściwością cząstki\n")
        f.write("- Masa wynika z OPORU (tarcia) przy poruszaniu się przez pole próżni\n")
        f.write("- Mechanizm Higgsa w FIN = topologiczne tarcie!\n")
    else:
        f.write("## 3. Wnioski\n")
        f.write("### 🟡 WYNIKI NIEPRZEKONUJĄCE\n")
        f.write("Brak wyraźnej korelacji. Możliwe że:\n")
        f.write("- Potrzebna większa sieć\n")
        f.write("- Siła jest za słaba lub za silna\n")
        f.write("- Efekty graniczne dominują\n")

print("Report saved.")
print("="*80)
