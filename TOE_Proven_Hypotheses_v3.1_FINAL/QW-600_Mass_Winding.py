#!/usr/bin/env python3
# QW-600: MASS FROM TOPOLOGICAL WINDING (Correct H5 Test)
# Purpose: Test if hopfion mass correlates with WINDING NUMBER, not beta_tors
# Based on FIN_Theory_Paper.tex: m = m0 × |w| × A
# Date: 2025-12-05

import numpy as np
from scipy.ndimage import convolve

print("="*80)
print("QW-600: MASS FROM TOPOLOGICAL WINDING (Revised H5)")
print("="*80)
print("Based on FIN Theory: m = m0 × |winding| × Amplification")
print("="*80)

# ============================================================================
# CONCEPT: Create hopfions with different winding numbers and measure inertia
# If theory correct: m_eff ∝ |winding|
# ============================================================================

GRID_SIZE = 32
DT = 0.01
STEPS = 400
GAMMA = 0.3
BETA_TORS = 0.01  # Fixed (not varying!)

# Test different winding numbers
WINDINGS = [1, 2, 3]

print(f"Grid: {GRID_SIZE}^3")
print(f"Testing windings: {WINDINGS}")
print(f"Beta_tors: {BETA_TORS} (FIXED)")
print("-" * 40)

# ============================================================================
# HOPFION WITH VARIABLE WINDING
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
    
    # WINDING enters here!
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
# EVOLUTION WITH FORCE
# ============================================================================

laplacian_kernel = np.zeros((3,3,3))
laplacian_kernel[1,1,1] = -6
for i in [(0,1,1), (2,1,1), (1,0,1), (1,2,1), (1,1,0), (1,1,2)]:
    laplacian_kernel[i] = 1

def laplacian(field):
    return convolve(field, laplacian_kernel, mode='wrap')

results = []
FORCE_STRENGTH = 0.05

print("\nTesting effective mass vs winding number...")

for winding in WINDINGS:
    print(f"\nWinding: {winding}")
    
    psi = hopfion_field(GRID_SIZE, center=(0, 0, 0), R=3.0, winding=winding)
    
    positions = []
    
    for t in range(STEPS):
        rho = np.abs(psi)**2
        attractor = GAMMA * (1.0 - rho) * psi
        kin = 1j * laplacian(psi)
        back = -1j * BETA_TORS * rho * psi
        
        # External force in x direction
        x = np.arange(GRID_SIZE) - GRID_SIZE/2
        y = np.arange(GRID_SIZE) - GRID_SIZE/2
        z = np.arange(GRID_SIZE) - GRID_SIZE/2
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        
        force = 1j * FORCE_STRENGTH * X * psi
        
        dpsi_dt = attractor + kin + back + force
        psi += DT * dpsi_dt
        
        if t % 20 == 0:
            pos = find_core_position(psi)
            positions.append(pos[0])
    
    # Compute acceleration
    if len(positions) > 3:
        times = np.arange(len(positions)) * 20 * DT
        velocities = np.diff(positions) / (20 * DT)
        if len(velocities) > 1:
            acceleration = np.mean(np.diff(velocities) / (20 * DT))
        else:
            acceleration = 0
        
        if abs(acceleration) > 1e-6:
            m_eff = FORCE_STRENGTH / abs(acceleration)
        else:
            m_eff = float('inf')
        
        results.append({
            'winding': winding,
            'acceleration': acceleration,
            'm_eff': m_eff,
            'final_pos': positions[-1]
        })
        
        print(f"  Acceleration: {acceleration:.4f}")
        print(f"  Effective mass: {m_eff:.3f}")
        print(f"  Final position: {positions[-1]:.2f}")

print("-" * 40)

# ============================================================================
# ANALYSIS
# ============================================================================

print("\nMass-Winding Relationship:")
print("| Winding |w| | m_eff | Ratio m/|w| |")
print("|-------------|-------|-------------|")
for r in results:
    ratio = r['m_eff'] / r['winding'] if r['winding'] > 0 else 0
    print(f"| {r['winding']:11d} | {r['m_eff']:5.3f} | {ratio:11.3f} |")

# Check if m_eff ∝ winding (linear relationship)
windings = [r['winding'] for r in results]
masses = [r['m_eff'] for r in results if r['m_eff'] < 1000]

if len(masses) >= 2:
    correlation = np.corrcoef(windings[:len(masses)], masses)[0, 1]
    
    print(f"\nCorrelation (m_eff vs |winding|): {correlation:.3f}")
    
    if correlation > 0.8:
        print("\n✅ H5 POTWIERDZONE (POPRAWNIE):")
        print("   Masa hopfiona ROŚNIE z liczbą wirową |w|!")
        print("   m ∝ |winding| (zgodnie z FIN Theory Paper)")
        confirmed = True
    else:
        print("\n🟡 SŁABA KORELACJA")
        confirmed = False
else:
    confirmed = False
    correlation = 0

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw600_mass_winding.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-600: Mass from Topological Winding\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Hipoteza:** H5 (REVISED) - Masa = Topologia × Rezonans\n\n")
    
    f.write("## 1. Korekta Założeń\n")
    f.write("**QW-598 testowało błędną hipotezę!**\n\n")
    f.write("Według FIN_Theory_Paper.tex, masa NIE pochodzi z beta_tors, ale z:\n")
    f.write("```\n")
    f.write("m = m0 × |winding| × Amplification\n")
    f.write("```\n\n")
    f.write("Beta_tors pojawia się tylko jako korekta tłumienia, nie źródło masy.\n\n")
    
    f.write("## 2. Nowa Metodologia\n")
    f.write("Testowano hopfiony o różnych **liczbach wirowych** (winding = 1, 2, 3)\n")
    f.write("przy **stałym** beta_tors = 0.01.\n\n")
    
    f.write("## 3. Wyniki\n")
    f.write("| Winding | m_eff | m/w |\n")
    f.write("|---------|-------|-----|\n")
    for r in results:
        ratio = r['m_eff'] / r['winding']
        f.write(f"| {r['winding']} | {r['m_eff']:.3f} | {ratio:.3f} |\n")
    f.write("\n")
    
    if confirmed:
        f.write("## 4. Wnioski\n")
        f.write("### ✅ HIPOTEZA H5 POTWIERDZONA (Poprawiona Wersja)\n")
        f.write(f"Korelacja: {correlation:.3f} > 0.8\n\n")
        f.write("**Masa hopfiona jest proporcjonalna do liczby wirowej (topologii)!**\n\n")
        f.write("To oznacza że:\n")
        f.write("- Masa wynika z STRUKTURY TOPOLOGICZNEJ cząstki (winding number)\n")
        f.write("- Beta_tors NIE jest źródłem masy (to był błąd w QW-598)\n")
        f.write("- Zgodność z FIN Theory Paper: m ∝ |w|\n\n")
        f.write("**Rewizja H5:**\n") 
        f.write("- Stara (błędna): Masa = Opór (beta_tors)\n")
        f.write("- Nowa (poprawna): Masa = Topologia (|winding|) × Amplifikacja (rezonans)\n")
    else:
        f.write("## 4. Wnioski\n")
        f.write("### 🟡 WYNIKI NIEPRZEKONUJĄCE\n")
        f.write("Korelacja za słaba. Możliwe przyczyny:\n")
        f.write("- Wyższe windingi wymagają większej siatki\n")
        f.write("- Amplifikacja (rezonans) nie jest uwzględniona\n")
        f.write("- Struktura hopfiona +2, +3 może być niestabilna\n")

print("Report saved.")
print("="*80)
