#!/usr/bin/env python3
# QW-595 REVISED: PARTICLE-PARTICLE INTERACTIONS
# Purpose: Test force emergence between hopfions with proper core tracking
# Date: 2025-12-05

import numpy as np
from scipy.ndimage import convolve

print("="*80)
print("QW-595 REVISED: PARTICLE-PARTICLE INTERACTIONS")
print("="*80)

# ============================================================================
# 1. PARAMETERS
# ============================================================================
GRID_SIZE = 40
DT = 0.01
STEPS = 600
GAMMA = 0.3  # Lower for stability
BETA_TORS = 0.01

# Test configurations
TEST_CONFIG = "same_winding"  # Options: "same_winding", "opposite_winding"

if TEST_CONFIG == "same_winding":
    WINDING_A = +1
    WINDING_B = +1
    SEPARATION = 16.0
    print("Test: SAME WINDING (+1, +1) - Expect REPULSION")
elif TEST_CONFIG == "opposite_winding":
    WINDING_A = +1
    WINDING_B = -1
    SEPARATION = 16.0
    print("Test: OPPOSITE WINDING (+1, -1) - Expect ATTRACTION")

print(f"Grid: {GRID_SIZE}^3")
print(f"Initial separation: {SEPARATION}")
print("-" * 40)

# ============================================================================
# 2. INITIALIZE TWO HOPFIONS
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

# Create hopfions on x-axis
center_A = (-SEPARATION/2, 0, 0)
center_B = (+SEPARATION/2, 0, 0)

psi_A = hopfion_field(GRID_SIZE, center_A, R=3.0, winding=WINDING_A)
psi_B = hopfion_field(GRID_SIZE, center_B, R=3.0, winding=WINDING_B)

# Normalize each separately to avoid initial clash
psi_A = psi_A / (np.max(np.abs(psi_A)) + 1e-10)
psi_B = psi_B / (np.max(np.abs(psi_B)) + 1e-10)

# Superposition with weighting to avoid overlap
psi = 0.7 * psi_A + 0.7 * psi_B

print(f"Hopfion A: center={center_A}, winding={WINDING_A:+d}")
print(f"Hopfion B: center={center_B}, winding={WINDING_B:+d}")

# ============================================================================
# 3. IMPROVED VORTEX CORE TRACKING
# ============================================================================

def find_vortex_cores_improved(psi):
    """
    Find cores using center-of-mass of low-density regions.
    """
    density = np.abs(psi)**2
    
    # Threshold: cores are where density < threshold
    threshold = 0.1 * np.mean(density)
    core_mask = density < threshold
    
    # Grid coordinates
    x = np.arange(GRID_SIZE) - GRID_SIZE/2
    y = np.arange(GRID_SIZE) - GRID_SIZE/2
    z = np.arange(GRID_SIZE) - GRID_SIZE/2
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    # Separate by hemisphere
    left_mask = (X < 0) & core_mask
    right_mask = (X >= 0) & core_mask
    
    if np.sum(left_mask) > 0:
        x_A = np.sum(X[left_mask]) / np.sum(left_mask)
        y_A = np.sum(Y[left_mask]) / np.sum(left_mask)
        z_A = np.sum(Z[left_mask]) / np.sum(left_mask)
        pos_A = np.array([x_A, y_A, z_A])
    else:
        pos_A = np.array([-SEPARATION/2, 0, 0])
    
    if np.sum(right_mask) > 0:
        x_B = np.sum(X[right_mask]) / np.sum(right_mask)
        y_B = np.sum(Y[right_mask]) / np.sum(right_mask)
        z_B = np.sum(Z[right_mask]) / np.sum(right_mask)
        pos_B = np.array([x_B, y_B, z_B])
    else:
        pos_B = np.array([+SEPARATION/2, 0, 0])
    
    return pos_A, pos_B

initial_pos_A, initial_pos_B = find_vortex_cores_improved(psi)
initial_distance = np.linalg.norm(initial_pos_A - initial_pos_B)
print(f"Initial core positions: A={initial_pos_A}, B={initial_pos_B}")
print(f"Initial separation: {initial_distance:.2f}")

# ============================================================================
# 4. EVOLUTION
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
history_pos_A = []
history_pos_B = []

for t in range(STEPS):
    rho = np.abs(psi)**2
    attractor = GAMMA * (1.0 - rho) * psi
    kin = 1j * laplacian(psi)
    back = -1j * BETA_TORS * rho * psi
    
    dpsi_dt = attractor + kin + back
    psi += DT * dpsi_dt
    
    # Measure every 50 steps
    if t % 50 == 0:
        pos_A, pos_B = find_vortex_cores_improved(psi)
        distance = np.linalg.norm(pos_A - pos_B)
        
        # Energy
        grad_x = (np.roll(psi, -1, axis=0) - np.roll(psi, 1, axis=0)) / 2.0
        grad_y = (np.roll(psi, -1, axis=1) - np.roll(psi, 1, axis=1)) / 2.0
        grad_z = (np.roll(psi, -1, axis=2) - np.roll(psi, 1, axis=2)) / 2.0
        grad_sq = np.abs(grad_x)**2 + np.abs(grad_y)**2 + np.abs(grad_z)**2
        pot = (1.0 - rho)**2
        E = np.sum(grad_sq + pot)
        
        history_distance.append(distance)
        history_energy.append(E)
        history_pos_A.append(pos_A.copy())
        history_pos_B.append(pos_B.copy())
        
        print(f"t={t:4d}: d={distance:.2f}, E={E:.1f}, A=({pos_A[0]:.1f}, {pos_A[1]:.1f}, {pos_A[2]:.1f})")

print("-" * 40)

# ============================================================================
# 5. ANALYSIS
# ============================================================================

initial_d = history_distance[0]
final_d = history_distance[-1]
delta_d = final_d - initial_d

initial_E = history_energy[0]
final_E = history_energy[-1]
delta_E = final_E - initial_E

print(f"\nInitial distance: {initial_d:.2f}")
print(f"Final distance:   {final_d:.2f}")
print(f"Change:           {delta_d:+.2f}")
print(f"Energy change:    {delta_E:+.1f}")
print()

# Determine interaction type
if delta_d < -2.0:
    result = "ATTRACTION"
    symbol = "✅"
elif delta_d > 2.0:
    result = "REPULSION"
    symbol = "✅"
else:
    result = "NEUTRAL"
    symbol = "🟡"

print(f"{symbol} {result}")

# Check prediction
if TEST_CONFIG == "same_winding" and result == "REPULSION":
    print("   ✅ Prediction CORRECT: Same winding → Repulsion")
elif TEST_CONFIG == "opposite_winding" and result == "ATTRACTION":
    print("   ✅ Prediction CORRECT: Opposite winding → Attraction")
else:
    print("   ❌ Prediction FAILED or inconclusive")

# ============================================================================
# 6. REPORT
# ============================================================================

report_path = "raport_qw595_revised.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-595: Particle-Particle Interactions (REVISED)\n")
    f.write("**Data:** 2025-12-05\n")
    f.write(f"**Test:** {TEST_CONFIG}\n\n")
    
    f.write("## 1. Setup\n")
    f.write(f"- Grid: {GRID_SIZE}^3\n")
    f.write(f"- Windings: A={WINDING_A:+d}, B={WINDING_B:+d}\n")
    f.write(f"- Initial Separation: {SEPARATION}\n\n")
    
    f.write("## 2. Wyniki\n")
    f.write(f"- **Initial Distance:** {initial_d:.2f}\n")
    f.write(f"- **Final Distance:** {final_d:.2f}\n")
    f.write(f"- **Change:** {delta_d:+.2f}\n")
    f.write(f"- **Result:** {symbol} **{result}**\n\n")
    
    f.write("## 3. Trajektoria Wirów\n")
    f.write("| Step | Distance | Energy | A_x | B_x |\n")
    f.write("|------|----------|--------|-----|-----|\n")
    for i in range(len(history_distance)):
        step = i * 50
        f.write(f"| {step} | {history_distance[i]:.2f} | {history_energy[i]:.1f} | {history_pos_A[i][0]:.1f} | {history_pos_B[i][0]:.1f} |\n")
    f.write("\n")
    
    f.write("## 4. Wnioski\n")
    if result == "REPULSION":
        f.write("### ✅ ODPYCHANIE POTWIERDZONE\n")
        f.write("Hopfiony o tym samym windings się odpychają. Siła emerge z topologii bez dodatkowych założeń!\n")
    elif result == "ATTRACTION":
        f.write("### ✅ PRZYCIĄGANIE POTWIERDZONE\n")
        f.write("Hopfiony o przeciwnych windings się przyciągają (antypartner). Emergencja sił potwierdzona!\n")
    else:
        f.write("### 🟡 BRAK WYRAŹNEJ SIŁY\n")
        f.write("Interakcja zbyt słaba lub wymaga dłuższej symulacji.\n")

print("Report saved.")
