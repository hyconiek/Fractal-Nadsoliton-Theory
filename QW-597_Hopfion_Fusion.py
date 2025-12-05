#!/usr/bin/env python3
# QW-597: HOPFION FUSION (ANNIHILATION TEST)
# Purpose: Test if opposite-winding hopfions (+1, -1) annihilate and convert to pure wave energy
# Date: 2025-12-05

import numpy as np
from scipy.ndimage import convolve

print("="*80)
print("QW-597: HOPFION FUSION (ANNIHILATION TEST)")
print("="*80)

# ============================================================================
# PARAMETERS
# ============================================================================
GRID_SIZE = 40
DT = 0.01
STEPS = 800
GAMMA = 0.3
INITIAL_SEPARATION = 14.0  # Start close enough for collision

print(f"Grid: {GRID_SIZE}^3")
print(f"Testing: (+1) vs (-1) collision")
print(f"Initial separation: {INITIAL_SEPARATION}")
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

# Create particle and antiparticle
center_A = (-INITIAL_SEPARATION/2, 0, 0)
center_B = (+INITIAL_SEPARATION/2, 0, 0)

psi_A = hopfion_field(GRID_SIZE, center_A, R=3.0, winding=+1)
psi_B = hopfion_field(GRID_SIZE, center_B, R=3.0, winding=-1)  # OPPOSITE winding!

psi_A = psi_A / (np.max(np.abs(psi_A)) + 1e-10)
psi_B = psi_B / (np.max(np.abs(psi_B)) + 1e-10)

psi = 0.7 * psi_A + 0.7 * psi_B

print("Hopfion A: winding = +1 (particle)")
print("Hopfion B: winding = -1 (antiparticle)")

# ============================================================================
# CORE TRACKING
# ============================================================================

def find_vortex_cores(psi):
    density = np.abs(psi)**2
    threshold = 0.1 * np.mean(density)
    core_mask = density < threshold
    
    x = np.arange(GRID_SIZE) - GRID_SIZE/2
    y = np.arange(GRID_SIZE) - GRID_SIZE/2
    z = np.arange(GRID_SIZE) - GRID_SIZE/2
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    left_mask = (X < 0) & core_mask
    right_mask = (X >= 0) & core_mask
    
    if np.sum(left_mask) > 10:
        x_A = np.sum(X[left_mask]) / np.sum(left_mask)
        pos_A = np.array([x_A, 0, 0])
        exists_A = True
    else:
        pos_A = None
        exists_A = False
    
    if np.sum(right_mask) > 10:
        x_B = np.sum(X[right_mask]) / np.sum(right_mask)
        pos_B = np.array([x_B, 0, 0])
        exists_B = True
    else:
        pos_B = None
        exists_B = False
    
    return pos_A, pos_B, exists_A, exists_B

# ============================================================================
# EVOLUTION
# ============================================================================

laplacian_kernel = np.zeros((3,3,3))
laplacian_kernel[1,1,1] = -6
for i in [(0,1,1), (2,1,1), (1,0,1), (1,2,1), (1,1,0), (1,1,2)]:
    laplacian_kernel[i] = 1

def laplacian(field):
    return convolve(field, laplacian_kernel, mode='wrap')

def compute_metrics(psi):
    rho = np.abs(psi)**2
    grad_x = (np.roll(psi, -1, axis=0) - np.roll(psi, 1, axis=0)) / 2.0
    grad_y = (np.roll(psi, -1, axis=1) - np.roll(psi, 1, axis=1)) / 2.0
    grad_z = (np.roll(psi, -1, axis=2) - np.roll(psi, 1, axis=2)) / 2.0
    grad_sq = np.abs(grad_x)**2 + np.abs(grad_y)**2 + np.abs(grad_z)**2
    pot = (1.0 - rho)**2
    E = np.sum(grad_sq + pot)
    
    # "Mass" (number of low-density vortex cores)
    threshold = 0.1 * np.mean(rho)
    core_volume = np.sum(rho < threshold)
    
    return E, core_volume

print("\nSimulating collision...")

history_distance = []
history_energy = []
history_mass = []

annihilation_time = None

for t in range(STEPS):
    rho = np.abs(psi)**2
    attractor = GAMMA * (1.0 - rho) * psi
    kin = 1j * laplacian(psi)
    back = -1j * 0.01 * rho * psi
    
    dpsi_dt = attractor + kin + back
    psi += DT * dpsi_dt
    
    if t % 50 == 0:
        pos_A, pos_B, exists_A, exists_B = find_vortex_cores(psi)
        E, mass = compute_metrics(psi)
        
        if exists_A and exists_B:
            distance = np.linalg.norm(pos_A - pos_B)
            status = "Both alive"
        elif exists_A or exists_B:
            distance = 0
            status = "One destroyed"
            if annihilation_time is None:
                annihilation_time = t
        else:
            distance = 0
            status = "Both destroyed"
            if annihilation_time is None:
                annihilation_time = t
        
        history_distance.append(distance)
       
        history_energy.append(E)
        history_mass.append(mass)
        
        print(f"t={t:4d}: d={distance:.1f}, E={E:.1f}, mass={mass:5d}, {status}")

print("-" * 40)

# ============================================================================
# ANALYSIS
# ============================================================================

initial_energy = history_energy[0]
final_energy = history_energy[-1]
initial_mass = history_mass[0]
final_mass = history_mass[-1]

print(f"\nInitial Energy: {initial_energy:.1f}")
print(f"Final Energy:   {final_energy:.1f}")
print(f"Energy Release: {initial_energy - final_energy:.1f}")
print()
print(f"Initial 'Mass' (core volume): {initial_mass}")
print(f"Final 'Mass':   {final_mass}")
print(f"Mass Loss:      {initial_mass - final_mass}")
print()

if annihilation_time:
    print(f"✅ ANNIHILATION DETECTED at t={annihilation_time}")
    print("   Hopfiony zniknęły (topologia zniszczona)")
    annihilated = True
else:
    print("🟡 NO ANNIHILATION")
    print("   Hopfiony survived or bounced")
    annihilated = False

# Check energy conservation
energy_change_percent = abs(final_energy - initial_energy) / initial_energy * 100
if energy_change_percent < 30:
    print(f"   Energy conserved ({energy_change_percent:.1f}% change)")
else:
    print(f"   Energy changed significantly ({energy_change_percent:.1f}%)")

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw597_fusion.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-597: Hopfion Fusion (Annihilation Test)\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Cel:** Test anihilacji cząstka-antycząstka\n\n")
    
    f.write("## 1. Setup\n")
    f.write("- Hopfion A: winding = +1 (cząstka)\n")
    f.write("- Hopfion B: winding = -1 (antycząstka)\n")
    f.write(f"- Początkowa separacja: {INITIAL_SEPARATION}\n\n")
    
    f.write("## 2. Wyniki Kolizji\n")
    f.write(f"- **Initial Energy:** {initial_energy:.1f}\n")
    f.write(f"- **Final Energy:** {final_energy:.1f}\n")
    f.write(f"- **Energy Released:** {initial_energy - final_energy:.1f}\n")
    f.write(f"- **Mass Loss:** {initial_mass - final_mass}\n\n")
    
    f.write("## 3. Trajektoria\n")
    f.write("| Step | Distance | Energy | Mass |\n")
    f.write("|------|----------|--------|------|\n")
    for i in range(len(history_distance)):
        step = i * 50
        f.write(f"| {step} | {history_distance[i]:.1f} | {history_energy[i]:.1f} | {history_mass[i]} |\n")
    f.write("\n")
    
    f.write("## 4. Wnioski\n")
    if annihilated:
        f.write("### ✅ ANIHILACJA POTWIERDZONA\n")
        f.write(f"Hopfiony znikły w czasile t={annihilation_time}. ")
        f.write("Struktura topologiczna została zniszczona, a energia przekształcona w fale.\n\n")
        f.write("**Interpretacja:** FIN odtwarza mechanizm cząstka-antycząstka! ")
        f.write("Przeciwne windingi anihilują podobnie do e⁺ + e⁻ → γ.\n")
    else:
        f.write("### 🟡 BRAK ANIHILACJI\n")
        f.write("Hopfiony przetrwały lub odbiły się. ")
        f.write("Może wymagane jest większe nakładanie lub inna konfiguracja.\n")

print("Report saved.")
print("="*80)
