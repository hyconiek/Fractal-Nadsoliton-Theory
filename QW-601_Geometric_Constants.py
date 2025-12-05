#!/usr/bin/env python3
# QW-601: GEOMETRIC ORIGIN OF COUPLING CONSTANTS
# Purpose: Test if fundamental constants emerge from network topology/geometry
# Hypothesis H7: Stałe Fizyczne to Liczby Geometryczne
# Date: 2025-12-05

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.linalg import eigh

print("="*80)
print("QW-601: GEOMETRIC ORIGIN OF COUPLING CONSTANTS")
print("="*80)
print("Testing H7: Czy stałe wynikają z geometrii sieci?")
print("="*80)

# ============================================================================
# CONCEPT: Build network from pure geometry, measure emergent constants
# If H7 true: alpha_geo, omega, phi emerge from topology
# If H7 false: they are free parameters
# ============================================================================

print("\n1. TESTING: Network Eigenvalue Spectrum")
print("-" * 40)

# Create different network topologies
TOPOLOGIES = {
    'random': 'Random Erdős-Rényi',
    'regular': 'Regular lattice',
    'fractal': 'Hierarchical fractal',
    'fcc': 'Face-centered cubic'
}

N_NODES = 100

results = {}

for topo_name, topo_desc in TOPOLOGIES.items():
    print(f"\nTopology: {topo_desc}")
    
    if topo_name == 'random':
        # Erdős-Rényi random graph
        np.random.seed(601)
        positions = np.random.randn(N_NODES, 3)
    
    elif topo_name == 'regular':
        # Regular cubic lattice
        n = int(N_NODES**(1/3))
        x = np.linspace(0, 1, n)
        X, Y, Z = np.meshgrid(x, x, x)
        positions = np.column_stack([X.ravel()[:N_NODES], 
                                      Y.ravel()[:N_NODES], 
                                      Z.ravel()[:N_NODES]])
    
    elif topo_name == 'fractal':
        # Hierarchical  fractal (recursive subdivision)
        np.random.seed(602)
        positions = []
        for level in range(3):
            n_per_level = N_NODES // 3
            scale = 2**(-level)
            pos_level = np.random.randn(n_per_level, 3) * scale
            positions.append(pos_level)
        positions = np.vstack(positions)[:N_NODES]
    
    elif topo_name == 'fcc':
        # FCC lattice
        n = int((N_NODES/4)**(1/3)) + 1
        positions = []
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    positions.append([i, j, k])
                    if len(positions) < N_NODES:
                        positions.append([i+0.5, j+0.5, k])
                    if len(positions) < N_NODES:
                        positions.append([i+0.5, j, k+0.5])
                    if len(positions) < N_NODES:
                        positions.append([i, j+0.5, k+0.5])
                    if len(positions) >= N_NODES:
                        break
                if len(positions) >= N_NODES:
                    break
            if len(positions) >= N_NODES:
                break
        positions = np.array(positions)[:N_NODES]
    
    # Distance matrix
    dist_matrix = squareform(pdist(positions))
    
    # Connection matrix (K_geo only - pure geometric decay)
    # K(d) = alpha * exp(-d/lambda)
    # Choose lambda = 1 (characteristic length)
    K_matrix = np.exp(-dist_matrix)
    np.fill_diagonal(K_matrix, 0)
    
    # Eigenvalues
    eigenvalues = eigh(K_matrix, eigvals_only=True)
    eigenvalues = np.sort(eigenvalues)[::-1]  # Descending
    
    # Key ratios (looking for geometric constants)
    lambda_1 = eigenvalues[0]
    lambda_2 = eigenvalues[1]
    lambda_3 = eigenvalues[2]
    
    ratio_21 = lambda_2 / lambda_1
    ratio_32 = lambda_3 / lambda_2
    
    # Spectral gap
    gap = lambda_1 - lambda_2
    
    results[topo_name] = {
        'lambda_1': lambda_1,
        'lambda_2': lambda_2,
        'lambda_3': lambda_3,
        'ratio_21': ratio_21,
        'ratio_32': ratio_32,
        'gap': gap
    }
    
    print(f"  λ₁ = {lambda_1:.4f}")
    print(f"  λ₂ = {lambda_2:.4f}")
    print(f"  λ₃ = {lambda_3:.4f}")
    print(f"  λ₂/λ₁ = {ratio_21:.4f}")
    print(f"  Gap = {gap:.4f}")

print("\n" + "="*80)
print("2. SEARCHING FOR UNIVERSAL RATIOS")
print("="*80)

# Check if any ratio is topology-independent (universal)
ratios_21 = [r['ratio_21'] for r in results.values()]
ratios_32 = [r['ratio_32'] for r in results.values()]

std_21 = np.std(ratios_21)
std_32 = np.std(ratios_32)
mean_21 = np.mean(ratios_21)
mean_32 = np.mean(ratios_32)

print(f"\nλ₂/λ₁ across topologies:")
print(f"  Mean: {mean_21:.4f}")
print(f"  Std:  {std_21:.4f}")
print(f"  Variation: {std_21/mean_21*100:.1f}%")

print(f"\nλ₃/λ₂ across topologies:")
print(f"  Mean: {mean_32:.4f}")
print(f"  Std:  {std_32:.4f}")
print(f"  Variation: {std_32/mean_32*100:.1f}%")

# If variation < 10%, it's universal (topology-independent)
if std_21/mean_21 < 0.1:
    print(f"\n✅ UNIVERSAL RATIO FOUND: λ₂/λ₁ ≈ {mean_21:.4f}")
    universal_found = True
else:
    print(f"\n❌ λ₂/λ₁ is TOPOLOGY-DEPENDENT (variation {std_21/mean_21*100:.1f}%)")
    universal_found = False

print("\n" + "="*80)
print("3. COMPARISON WITH FIN CONSTANTS")
print("="*80)

# Known FIN constants
ALPHA_GEO = 4 * np.log(2)  # 2.772588...
OMEGA = np.pi / 4           # 0.785398...
PHI = np.pi / 6             # 0.523599...
KAPPA = 6.74               # From QW-481

print(f"\nKnown FIN constants:")
print(f"  α_geo = {ALPHA_GEO:.6f}")
print(f"  ω = {OMEGA:.6f}")
print(f"  φ = {PHI:.6f}")
print(f"  κ = {KAPPA:.2f}")

print(f"\nDo network ratios match FIN constants?")

# Check if any eigenvalue ratio matches known constants
for const_name, const_value in [('α_geo', ALPHA_GEO), ('ω', OMEGA), 
                                  ('φ', PHI), ('κ', KAPPA)]:
    closest_ratio = min(ratios_21 + ratios_32, key=lambda x: abs(x - const_value))
    error = abs(closest_ratio - const_value) / const_value
    
    if error < 0.1:  # Within 10%
        print(f"  {const_name}: MATCH with ratio {closest_ratio:.6f} (error {error*100:.1f}%)")
    else:
        print(f"  {const_name}: NO MATCH (closest {closest_ratio:.6f}, error {error*100:.1f}%)")

print("\n" + "="*80)
print("4. CONCLUSION")
print("="*80)

if universal_found:
    print("\n✅ CZĘŚCIOWE POTWIERDZENIE H7:")
    print(f"   Found topology-independent ratio: {mean_21:.4f}")
    print("   This suggests SOME constants emerge from pure geometry")
    confirmed = "partial"
else:
    print("\n❌ H7 NIEPOTWIERDZONY:")
    print("   Eigenvalue ratios are topology-dependent")
    print("   Constants α_geo, ω, φ nie wynikają z czystej geometrii")
    confirmed = "no"

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw601_geometric_constants.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-601: Geometric Origin of Coupling Constants\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Hipoteza:** H7 - Stałe Fizyczne = Liczby Geometryczne\n\n")
    
    f.write("## 1. Metodologia\n")
    f.write("Zbadano 4 różne topologie sieci:\n")
    for name, desc in TOPOLOGIES.items():
        f.write(f"- {desc}\n")
    f.write("\nObliczono spektrum własne macierzy połączeń K(d) = exp(-d).\n\n")
    
    f.write("## 2. Wyniki Spektralne\n")
    f.write("| Topologia | λ₁ | λ₂ | λ₂/λ₁ | Gap |\n")
    f.write("|-----------|-----|-----|--------|-----|\n")
    for name, r in results.items():
        f.write(f"| {TOPOLOGIES[name]} | {r['lambda_1']:.3f} | {r['lambda_2']:.3f} | ")
        f.write(f"{r['ratio_21']:.3f} | {r['gap']:.3f} |\n")
    f.write("\n")
    
    f.write("## 3. Analiza Uniwersalności\n")
    f.write(f"- **λ₂/λ₁ wariancja:** {std_21/mean_21*100:.1f}%\n")
    if std_21/mean_21 < 0.1:
        f.write(f"- **Wniosek:** Stosunek jest uniwersalny (niezależny od topologii)\n\n")
    else:
        f.write(f"- **Wniosek:** Stosunek zależy od topologii (nie jest uniwersalny)\n\n")
    
    f.write("## 4. Werdykt\n")
    if confirmed == "partial":
        f.write("### 🟡 CZĘŚCIOWE POTWIERDZENIE\n")
        f.write("Niektóre stosunki spektralne są uniwersalne, co sugeruje że część stałych\n")
        f.write("może wynikać z geometrii. Jednak α_geo, ω, φ nie pasują do znalezionych wartości.\n")
    else:
        f.write("### ❌ H7 NIEPOTWIERDZONY\n")
        f.write("Stałe FIN (α_geo=4ln2, ω=π/4, φ=π/6) NIE wyłaniają się naturalnie\n")
        f.write("ze spektrum prostych topologii geometrycznych.\n\n")
        f.write("**Implikacja:** Są to parametry fenomenologiczne, nie fundamentalne.\n")

print("Report saved.")
print("="*80)
