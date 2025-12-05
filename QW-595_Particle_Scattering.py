#!/usr/bin/env python3
# QW-595: PARTICLE-PARTICLE INTERACTIONS (HOPFION SCATTERING)
# Purpose: Test if forces emerge from topological dynamics between stable hopfions
# Date: 2025-12-05

import numpy as np
from scipy.ndimage import convolve

print("="*80)
print("QW-595: PARTICLE-PARTICLE INTERACTIONS (HOPFION SCATTERING)")
print("="*80)

# ============================================================================
# 1. PARAMETERS
# ============================================================================
GRID_SIZE = 48  # Larger grid to fit two hopfions
DT = 0.01
STEPS = 800
GAMMA = 0.5
BETA_TORS = 0.01

# Hopfion parameters
WINDING_A = +1  # First hopfion
WINDING_B = +1  # Second hopfion (same sign = repulsion, opposite = attraction?)
SEPARATION = 12.0  # Initial separation

print(f"Grid: {GRID_SIZE}^3")
print(f"Winding: A={WINDING_A:+d}, B={WINDING_B:+d}")
print(f"Initial separation: {SEPARATION}")
print("-" * 40)

# ============================================================================
# 2. INITIALIZE TWO HOPFIONS
# ============================================================================

def hopfion_field(grid_size, center, R=4.0, winding=1):
    """
    Create a hopfion centered at 'center' with given winding number.
    """
    x = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    y = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    z = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    # Shift to center
    X = X - center[0]
    Y = Y - center[1]
    Z = Z - center[2]
    
    rho = np.sqrt(X**2 + Y**2)
    rho[rho == 0] = 1e-10
    
    # Toroidal coordinates
    eta = np.arctan2(Z, rho - R)
    xi = np.arctan2(Y, X)
    
    # Phase with winding
    phase = winding * (xi + eta)
    
    # Amplitude
    dist_to_ring = np.sqrt((rho - R)**2 + Z**2)
    amplitude = np.tanh(dist_to_ring / 2.0)
    
    return amplitude * np.exp(1j * phase)

# Create two hopfions
center_A = (-SEPARATION/2, 0, 0)
center_B = (+SEPARATION/2, 0, 0)

psi_A = hopfion_field(GRID_SIZE, center_A, R=4.0, winding=WINDING_A)
psi_B = hopfion_field(GRID_SIZE, center_B, R=4.0, winding=WINDING_B)

# Superposition (linear combination - will evolve nonlinearly)
psi = psi_A + psi_B

print(f"Hopfion A center: {center_A}")
print(f"Hopfion B center: {center_B}")

# ============================================================================
# 3. VORTEX CORE TRACKING
# ============================================================================

def find_vortex_cores(psi):
    """
    Find vortex cores by locating regions where |psi| is minimal.
    Returns list of core positions.
    """
    density = np.abs(psi)**2
    
    # Smooth to avoid noise
    from scipy.ndimage import gaussian_filter
    density_smooth = gaussian_filter(density, sigma=2.0)
    
    # Find local minima (cores are where amplitude → 0)
    # Simple approach: divide grid into two halves and find min in each
    
    mid = GRID_SIZE // 2
    
    # Left half (Hopfion A)
    left_half = density_smooth[:mid, :, :]
    idx_A = np.unravel_index(np.argmin(left_half), left_half.shape)
    pos_A = np.array(idx_A) - GRID_SIZE/2
    
    # Right half (Hopfion B)
    right_half = density_smooth[mid:, :, :]
    idx_B = np.unravel_index(np.argmin(right_half), right_half.shape)
    pos_B = np.array(idx_B) + np.array([mid, 0, 0]) - GRID_SIZE/2
    
    return pos_A, pos_B

initial_pos_A, initial_pos_B = find_vortex_cores(psi)
print(f"Detected cores: A={initial_pos_A}, B={initial_pos_B}")

# ============================================================================
# 4. EVOLUTION WITH INTERACTION
# ============================================================================

laplacian_kernel = np.zeros((3,3,3))
laplacian_kernel[1,1,1] = -6
for i in [(0,1,1), (2,1,1), (1,0,1), (1,2,1), (1,1,0), (1,1,2)]:
    laplacian_kernel[i] = 1

def laplacian(field):
    return convolve(field, laplacian_kernel, mode='wrap')

print("\nSimulating scattering...")

history_distance = []
history_energy = []
history_charge = []

for t in range(STEPS):
    # 1. Attractor term
    rho = np.abs(psi)**2
    attractor = GAMMA * (1.0 - rho) * psi
    
    # 2. Kinetic
    kin = 1j * laplacian(psi)
    
    # 3. Self-interaction (back-reaction)
    back = -1j * BETA_TORS * rho * psi
    
    # Total
    dpsi_dt = attractor + kin + back
    psi += DT * dpsi_dt
    
    # Measurements
    if t % 50 == 0:
        pos_A, pos_B = find_vortex_cores(psi)
        distance = np.linalg.norm(pos_A - pos_B)
        
        # Total energy
        grad_x = (np.roll(psi, -1, axis=0) - np.roll(psi, 1, axis=0)) / 2.0
        grad_y = (np.roll(psi, -1, axis=1) - np.roll(psi, 1, axis=1)) / 2.0
        grad_z = (np.roll(psi, -1, axis=2) - np.roll(psi, 1, axis=2)) / 2.0
        grad_sq = np.abs(grad_x)**2 + np.abs(grad_y)**2 + np.abs(grad_z)**2
        pot = (1.0 - rho)**2
        E = np.sum(grad_sq + pot)
        
        # Topological charge (simplified - just total phase winding)
        phase = np.angle(psi)
        Q_approx = np.sum(np.abs(np.gradient(phase)[0])) / (2*np.pi)
        
        history_distance.append(distance)
        history_energy.append(E)
        history_charge.append(Q_approx)
        
        print(f"t={t:4d}: d={distance:.2f}, E={E:.1f}, Q~{Q_approx:.1f}")
        
        if distance < 3.0:
            print("⚠️ Hopfions very close - possible collision!")

print("-" * 40)

# ============================================================================
# 5. ANALYSIS
# ============================================================================

initial_distance = history_distance[0]
final_distance = history_distance[-1]
delta_distance = final_distance - initial_distance

initial_energy = history_energy[0]
final_energy = history_energy[-1]
delta_energy = final_energy - initial_energy

print(f"\nInitial distance: {initial_distance:.2f}")
print(f"Final distance:   {final_distance:.2f}")
print(f"Change:           {delta_distance:+.2f}")
print()

if delta_distance < -1.0:
    print("✅ ATTRACTION: Hopfions moved closer!")
    print(f"   Energy change: {delta_energy:+.1f} (should decrease if attractive)")
elif delta_distance > 1.0:
    print("✅ REPULSION: Hopfions moved apart!")
    print(f"   Energy change: {delta_energy:+.1f} (should increase if repulsive)")
else:
    print("🟡 NEUTRAL: No significant interaction.")

print()
print("Topological charge conservation:")
print(f"  Initial Q: {history_charge[0]:.2f}")
print(f"  Final Q:   {history_charge[-1]:.2f}")
if abs(history_charge[-1] - history_charge[0]) < 1.0:
    print("  ✅ Charge conserved")
else:
    print("  ❌ Charge not conserved (topology broken)")

# ============================================================================
# 6. REPORT GENERATION
# ============================================================================

report_path = "raport_qw595_scattering.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-595: Particle-Particle Interactions\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Metoda:** Hopfion Scattering (Ginzburg-Landau)\n\n")
    
    f.write("## 1. Setup\n")
    f.write(f"- Grid: {GRID_SIZE}^3\n")
    f.write(f"- Winding Numbers: A={WINDING_A:+d}, B={WINDING_B:+d}\n")
    f.write(f"- Initial Separation: {SEPARATION}\n\n")
    
    f.write("## 2. Wyniki Dynamiki\n")
    f.write(f"- **Initial Distance:** {initial_distance:.2f}\n")
    f.write(f"- **Final Distance:** {final_distance:.2f}\n")
    f.write(f"- **Change:** {delta_distance:+.2f}\n\n")
    
    f.write("## 3. Energia\n")
    f.write(f"- Initial: {initial_energy:.1f}\n")
    f.write(f"- Final: {final_energy:.1f}\n")
    f.write(f"- Change: {delta_energy:+.1f}\n\n")
    
    f.write("## 4. Trajektoria\n")
    f.write("| Step | Distance | Energy |\n")
    f.write("|------|----------|--------|\n")
    for i, (d, e) in enumerate(zip(history_distance, history_energy)):
        step = i * 50
        f.write(f"| {step} | {d:.2f} | {e:.1f} |\n")
    f.write("\n")
    
    f.write("## 5. Wnioski\n")
    if delta_distance < -1.0:
        f.write("### ✅ PRZYCIĄGANIE\n")
        f.write("Hopfiony zbliżyły się, co sugeruje emergencję siły przyciągającej z topologii.\n")
    elif delta_distance > 1.0:
        f.write("### ✅ ODPYCHANIE\n")
        f.write("Hopfiony oddaliły się, co sugeruje emergencję siły odpychającej z topologii.\n")
    else:
        f.write("### 🟡 BRAK INTERAKCJI\n")
        f.write("Nie zaobserwowano znaczącej siły między hopfionami.\n")

print("Report saved.")
