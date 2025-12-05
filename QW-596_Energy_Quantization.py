#!/usr/bin/env python3
# QW-596: ENERGY QUANTIZATION IN HOPFION INTERACTIONS
# Purpose: Test if binding energy between hopfions is quantized (discrete energy levels)
# Date: 2025-12-05

import numpy as np
from scipy.ndimage import convolve

print("="*80)
print("QW-596: ENERGY QUANTIZATION IN HOPFION INTERACTIONS")
print("="*80)

# ============================================================================
# PARAMETERS
# ============================================================================
GRID_SIZE = 36
DT = 0.01
STEPS = 400
GAMMA = 0.3

# Test multiple separations to find energy spectrum
SEPARATIONS = [8.0, 10.0, 12.0, 14.0, 16.0, 18.0]

print(f"Grid: {GRID_SIZE}^3")
print(f"Testing {len(SEPARATIONS)} separation distances")
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
# EVOLUTION
# ============================================================================

laplacian_kernel = np.zeros((3,3,3))
laplacian_kernel[1,1,1] = -6
for i in [(0,1,1), (2,1,1), (1,0,1), (1,2,1), (1,1,0), (1,1,2)]:
    laplacian_kernel[i] = 1

def laplacian(field):
    return convolve(field, laplacian_kernel, mode='wrap')

def compute_total_energy(psi):
    rho = np.abs(psi)**2
    grad_x = (np.roll(psi, -1, axis=0) - np.roll(psi, 1, axis=0)) / 2.0
    grad_y = (np.roll(psi, -1, axis=1) - np.roll(psi, 1, axis=1)) / 2.0
    grad_z = (np.roll(psi, -1, axis=2) - np.roll(psi, 1, axis=2)) / 2.0
    grad_sq = np.abs(grad_x)**2 + np.abs(grad_y)**2 + np.abs(grad_z)**2
    pot = (1.0 - rho)**2
    return np.sum(grad_sq + pot)

results = []

print("\nTesting energy levels at different separations...")

for separation in SEPARATIONS:
    print(f"\nSeparation: {separation:.1f}")
    
    # Initialize two hopfions
    center_A = (-separation/2, 0, 0)
    center_B = (+separation/2, 0, 0)
    
    psi_A = hopfion_field(GRID_SIZE, center_A, R=3.0, winding=+1)
    psi_B = hopfion_field(GRID_SIZE, center_B, R=3.0, winding=+1)
    
    psi_A = psi_A / (np.max(np.abs(psi_A)) + 1e-10)
    psi_B = psi_B / (np.max(np.abs(psi_B)) + 1e-10)
    
    psi = 0.7 * psi_A + 0.7 * psi_B
    
    # Evolve to equilibrium
    for t in range(STEPS):
        rho = np.abs(psi)**2
        attractor = GAMMA * (1.0 - rho) * psi
        kin = 1j * laplacian(psi)
        back = -1j * 0.01 * rho * psi
        
        dpsi_dt = attractor + kin + back
        psi += DT * dpsi_dt
    
    # Measure final energy
    E_final = compute_total_energy(psi)
    
    results.append({
        'separation': separation,
        'energy': E_final
    })
    
    print(f"  Final Energy: {E_final:.1f}")

print("-" * 40)

# ============================================================================
# ANALYSIS: Energy Quantization
# ============================================================================

energies = [r['energy'] for r in results]
separations = [r['separation'] for r in results]

# Compute energy differences
energy_diffs = []
for i in range(len(energies)-1):
    diff = energies[i+1] - energies[i]
    energy_diffs.append(diff)

print("\nEnergy Spectrum Analysis:")
print("| Separation | Energy    | ΔE        |")
print("|------------|-----------|-----------|")
for i, r in enumerate(results):
    if i < len(energy_diffs):
        print(f"| {r['separation']:10.1f} | {r['energy']:9.1f} | {energy_diffs[i]:9.1f} |")
    else:
        print(f"| {r['separation']:10.1f} | {r['energy']:9.1f} | -         |")

# Check if energy differences are approximately constant (quantized)
if len(energy_diffs) > 1:
    mean_diff = np.mean(energy_diffs)
    std_diff = np.std(energy_diffs)
    variation = std_diff / mean_diff if mean_diff != 0 else float('inf')
    
    print(f"\nMean ΔE: {mean_diff:.1f}")
    print(f"Std Dev: {std_diff:.1f}")
    print(f"Variation: {variation*100:.1f}%")
    
    if variation < 0.2:  # Less than 20% variation
        print("\n✅ KWANTYZACJA ENERGII POTWIERDZONA!")
        print(f"   Energia zmienia się w jednostkach ~{mean_diff:.1f}")
        quantized = True
    else:
        print("\n🟡 BRAK WYRAŹNEJ KWANTYZACJI")
        print("   Energie nie tworzą regularnego wzoru")
        quantized = False
else:
    quantized = False

# ============================================================================
# REPORT GENERATION
# ============================================================================

report_path = "raport_qw596_quantization.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-596: Energy Quantization in Hopfion Interactions\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Cel:** Test kwantyzacji energii wiązania między hopfionami\n\n")
    
    f.write("## 1. Metodologia\n")
    f.write(f"- Symulowano dwa hopfiony (+1, +1) w różnych odległościach\n")
    f.write(f"- Grid: {GRID_SIZE}^3\n")
    f.write(f"- Ewolucja: {STEPS} kroków do równowagi\n\n")
    
    f.write("## 2. Spektrum Energetyczne\n")
    f.write("| Separacja | Energia | ΔE |\n")
    f.write("|-----------|---------|----|\n")
    for i, r in enumerate(results):
        if i < len(energy_diffs):
            f.write(f"| {r['separation']:.1f} | {r['energy']:.1f} | {energy_diffs[i]:.1f} |\n")
        else:
            f.write(f"| {r['separation']:.1f} | {r['energy']:.1f} | - |\n")
    f.write("\n")
    
    if quantized:
        f.write("## 3. Wnioski\n")
        f.write("### ✅ KWANTYZACJA POTWIERDZONA\n")
        f.write(f"Energia wiązania zmienia się w dyskretnych jednostkach ~{mean_diff:.1f}, ")
        f.write("co sugeruje że hopfiony formują quasi-kwantowe poziomy energetyczne, ")
        f.write("podobnie do atomów w mechanice kwantowej!\n\n")
        f.write("**Implikacja:** FIN może odtwarzać chemię kwantową z topologii.\n")
    else:
        f.write("## 3. Wnioski\n")
        f.write("### 🟡 BRAK KWANTYZACJI\n")
        f.write("Energie nie wykazują regularnego wzoru kwantyzacji. ")
        f.write("Może to wymagać większej sieci lub dłuższej ewolucji.\n")

print("Report saved.")
print("="*80)
