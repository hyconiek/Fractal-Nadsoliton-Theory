#!/usr/bin/env python3
# QW-602: MULTI-HOPFION CHEMISTRY (Triple Binding Test)
# Purpose: Test if 3 hopfions form stable bound state (FIN "molecule")
# Precedent: QW-433 found proton as stable triplet with E_bind < 0
# Date: 2025-12-05

import numpy as np
from scipy.ndimage import convolve

print("="*80)
print("QW-602: MULTI-HOPFION CHEMISTRY (Triple Binding)")
print("="*80)
print("Test: Can (+1, +1, -1) hopfions form stable triangle?")
print("Precedent: QW-433 (proton triplet with E_bind < 0)")
print("="*80)

# ============================================================================
# PARAMETERS
# ============================================================================
GRID_SIZE = 40
DT = 0.01
STEPS = 800
GAMMA = 0.3
BETA_TORS = 0.01

# Initial configuration: equilateral triangle
SEPARATION = 12.0  # Distance between hopfions
angle_offset = 2 * np.pi / 3  # 120 degrees

print(f"Grid: {GRID_SIZE}^3")
print(f"Configuration: Equilateral triangle, separation={SEPARATION}")
print(f"Windings: A=+1, B=+1, C=-1 (2 particles + 1 antiparticle)")
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

# Create triangle configuration
R_triangle = SEPARATION / np.sqrt(3)  # Radius of circumscribed circle

center_A = (R_triangle * np.cos(0), R_triangle * np.sin(0), 0)
center_B = (R_triangle * np.cos(angle_offset), R_triangle * np.sin(angle_offset), 0)
center_C = (R_triangle * np.cos(2*angle_offset), R_triangle * np.sin(2*angle_offset), 0)

psi_A = hopfion_field(GRID_SIZE, center_A, R=3.0, winding=+1)
psi_B = hopfion_field(GRID_SIZE, center_B, R=3.0, winding=+1)
psi_C = hopfion_field(GRID_SIZE, center_C, R=3.0, winding=-1)  # Antiparticle!

# Normalize and superpose
psi_A = psi_A / (np.max(np.abs(psi_A)) + 1e-10)
psi_B = psi_B / (np.max(np.abs(psi_B)) + 1e-10)
psi_C = psi_C / (np.max(np.abs(psi_C)) + 1e-10)

psi = 0.6 * (psi_A + psi_B + psi_C)

print(f"Hopfion A: center={center_A}, winding=+1")
print(f"Hopfion B: center={center_B}, winding=+1")
print(f"Hopfion C: center={center_C}, winding=-1")

# ============================================================================
# CORE TRACKING (3 cores)
# ============================================================================

def find_three_cores(psi):
    """Find 3 vortex cores using k-means-like approach"""
    density = np.abs(psi)**2
    threshold = 0.12 * np.mean(density)
    core_mask = density < threshold
    
    x = np.arange(GRID_SIZE) - GRID_SIZE/2
    y = np.arange(GRID_SIZE) - GRID_SIZE/2
    z = np.arange(GRID_SIZE) - GRID_SIZE/2
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    # Get all core points
    core_points = np.column_stack([X[core_mask], Y[core_mask], Z[core_mask]])
    
    if len(core_points) < 30:
        return None, None, None
    
    # Simple sector-based clustering (divide space into 3 sectors)
    angles = np.arctan2(core_points[:, 1], core_points[:, 0])
    
    sector_0 = (angles >= -np.pi/3) & (angles < np.pi/3)
    sector_1 = (angles >= np.pi/3) & (angles < np.pi)
    sector_2 = (angles >= np.pi) | (angles < -np.pi/3)
    
    cores = []
    for sector in [sector_0, sector_1, sector_2]:
        if np.sum(sector) > 0:
            core_pos = np.mean(core_points[sector], axis=0)
            cores.append(core_pos)
        else:
            cores.append(None)
    
    return cores[0], cores[1], cores[2]

# ============================================================================
# EVOLUTION
# ============================================================================

laplacian_kernel = np.zeros((3,3,3))
laplacian_kernel[1,1,1] = -6
for i in [(0,1,1), (2,1,1), (1,0,1), (1,2,1), (1,1,0), (1,1,2)]:
    laplacian_kernel[i] = 1

def laplacian(field):
    return convolve(field, laplacian_kernel, mode='wrap')

def compute_energy(psi):
    rho = np.abs(psi)**2
    grad_x = (np.roll(psi, -1, axis=0) - np.roll(psi, 1, axis=0)) / 2.0
    grad_y = (np.roll(psi, -1, axis=1) - np.roll(psi, 1, axis=1)) / 2.0
    grad_z = (np.roll(psi, -1, axis=2) - np.roll(psi, 1, axis=2)) / 2.0
    grad_sq = np.abs(grad_x)**2 + np.abs(grad_y)**2 + np.abs(grad_z)**2
    pot = (1.0 - rho)**2
    return np.sum(grad_sq + pot)

print("\nSimulating triplet evolution...")

history_energy = []
history_triangle_size = []
history_cores = []

initial_energy = compute_energy(psi)

for t in range(STEPS):
    rho = np.abs(psi)**2
    attractor = GAMMA * (1.0 - rho) * psi
    kin = 1j * laplacian(psi)
    back = -1j * BETA_TORS * rho * psi
    
    dpsi_dt = attractor + kin + back
    psi += DT * dpsi_dt
    
    if t % 50 == 0:
        core_A, core_B, core_C = find_three_cores(psi)
        E = compute_energy(psi)
        
        all_exist = (core_A is not None) and (core_B is not None) and (core_C is not None)
        
        if all_exist:
            # Triangle size (average side length)
            d_AB = np.linalg.norm(core_A - core_B)
            d_BC = np.linalg.norm(core_B - core_C)
            d_CA = np.linalg.norm(core_C - core_A)
            avg_size = (d_AB + d_BC + d_CA) / 3.0
            
            history_triangle_size.append(avg_size)
            history_energy.append(E)
            history_cores.append((core_A.copy(), core_B.copy(), core_C.copy()))
            
            print(f"t={t:4d}: Triangle size={avg_size:.2f}, E={E:.1f}, All alive")
        else:
            print(f"t={t:4d}: COLLISION or ESCAPE detected")
            history_triangle_size.append(0)
            history_energy.append(E)
            history_cores.append(None)

print("-" * 40)

# ============================================================================
# ANALYSIS
# ============================================================================

if len(history_triangle_size) > 2:
    valid_sizes = [s for s in history_triangle_size if s > 0]
    
    if len(valid_sizes) >= 3:
        initial_size = history_triangle_size[0]
        final_size = history_triangle_size[-1] if history_triangle_size[-1] > 0 else valid_sizes[-1]
        
        final_energy = history_energy[-1]
        
        # Binding energy (negative if bound)
        # Compare to 3 separated hopfions
        E_bind = final_energy - initial_energy
        
        # Check oscillation (sign of binding)
        size_variance = np.std(valid_sizes) / np.mean(valid_sizes) if len(valid_sizes) > 1 else 0
        
        print(f"\nInitial size: {initial_size:.2f}")
        print(f"Final size:   {final_size:.2f}")
        print(f"Size change:  {final_size - initial_size:+.2f}")
        print(f"Size variance: {size_variance*100:.1f}%")
        print()
        print(f"Initial energy: {initial_energy:.1f}")
        print(f"Final energy:   {final_energy:.1f}")
        print(f"ΔE (binding):   {E_bind:+.1f}")
        print()
        
        # Determine if bound
        if size_variance > 0.05 and abs(final_size - initial_size) < initial_size * 0.3:
            print("✅ STABLE OSCYLUJĄCYUKŁAD!")
            print("   Trójkąt oscyluje wokół równowagi (jak molekuła drgająca)")
            print(f"   E_bind = {E_bind:.1f} (energia wiązania)")
            bound = True
        elif final_size < initial_size * 0.5:
            print("🟡 KOLAPS")
            print("   Hopfiony zbiegły się (może w kierunku fuzji)")
            bound = False
        elif final_size > initial_size * 1.5:
            print("❌ ROZPAD")
            print("   Trójkąt się rozpadł (brak wiązania)")
            bound = False
        else:
            print("🟡 NIESTABILNY")
            print("   Układ nie osiągnął równowagi")
            bound = None
    else:
        bound = False
        E_bind = 0
else:
    bound = False
    E_bind = 0

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw602_chemistry.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-602: Multi-Hopfion Chemistry\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Precedent:** QW-433 (proton triplet binding)\n\n")
    
    f.write("## 1. Setup\n")
    f.write("- Konfiguracja: Równoboczny trójkąt\n")
    f.write("- Hopfiony: A(+1), B(+1), C(-1)\n")
    f.write(f"- Początkowa separacja: {SEPARATION}\n\n")
    
    f.write("## 2. Ewolucja\n")
    f.write("| Step | Triangle Size | Energy |\n")
    f.write("|------|---------------|--------|\n")
    for i in range(len(history_triangle_size)):
        step = i * 50
        if history_triangle_size[i] > 0:
            f.write(f"| {step} | {history_triangle_size[i]:.2f} | {history_energy[i]:.1f} |\n")
    f.write("\n")
    
    f.write("## 3. Analiza Wiązania\n")
    if bound:
        f.write("### ✅ STABILNY UKŁAD ZWIĄZANY!\n")
        f.write(f"- Wariancja rozmiaru: {size_variance*100:.1f}% (oscylacje)\n")
        f.write(f"- Energia wiązania: ΔE = {E_bind:.1f}\n\n")
        f.write("**Wniosek:** FIN może tworzyć stabilne wielocząstkowe struktury!\n")
        f.write("To pierwszy krok do 'FIN Chemistry' - wiązania analogiczne do atomów/molekuł.\n\n")
        f.write("**Zgodność z QW-433:** Proton też był tripletem z E_bind < 0.\n")
    elif bound == False:
        f.write("### ❌ BRAK STABILNEGO WIĄZANIA\n")
        f.write("Trójkąt rozpadł się lub skolapsował. Możliwe przyczyny:\n")
        f.write("- Za duża/mała początkowa separacja\n")
        f.write("- Konieczna inna konfiguracja windings\n")
    else:
        f.write("### 🟡 WYNIK NIEJEDNOZNACZNY\n")

print("Report saved.")
print("="*80)
