#!/usr/bin/env python3
# QW-620: PROTON STRUCTURE (3-MODE BOUND STATE)
# Purpose: Test if 3 modes (quarks) bind in 12-octave network
# Model: Tensor product of 3 single-particle spaces (12x12x12 = 1728 states)
# Interaction: Pairwise attractive potentials V_12 + V_23 + V_31
# Date: 2025-12-05

import numpy as np
from scipy.linalg import eigh
from itertools import product

print("="*80)
print("QW-620: PROTON STRUCTURE (3-MODE BOUND STATE)")
print("="*80)
print("Test: Czy 3 mody (kwarki) tworzą stabilny proton w sieci oktaw?")
print("Hypothesis: Strong binding for triplet configuration")
print("="*80)

# ============================================================================
# PARAMETERS
# ============================================================================
N_OCTAVES = 12
ALPHA_GEO = 4 * np.log(2)
BETA_TORS = 0.01
COUPLING_STRENGTH = 5.0 # Same as QW-619

print(f"Network: {N_OCTAVES} octaves")
print(f"State space: 12^3 = {N_OCTAVES**3} basis states")
print("-" * 40)

# ============================================================================
# HAMILTONIAN CONSTRUCTION
# ============================================================================

def K(d):
    return ALPHA_GEO / (1 + BETA_TORS * d)

# Single particle energies (Vacuum)
H1 = np.zeros(N_OCTAVES)
for i in range(N_OCTAVES):
    # Self-energy or interaction with background?
    # Let's use simple diagonal energy for "mass" of localized mode
    # Assuming quarks have intrinsic mass from octave confinement
    H1[i] = 1.0 # Base mass

# Full basis indices
basis = list(product(range(N_OCTAVES), repeat=3)) # (i, j, k) tuples
N_STATES = len(basis)
idx_map = {state: i for i, state in enumerate(basis)}

# Build Sparse-like Hamiltonian (by iterating)
# H = H1 + H2 + H3 + V12 + V23 + V31

diag_H = np.zeros(N_STATES)
off_diag = {} # (row, col) -> value

print("Building Hamiltonian...")

for idx, (i, j, k) in enumerate(basis):
    # Kinetic/Mass Energy
    E_kin = H1[i] + H1[j] + H1[k]
    
    # Interaction Potentials (Pairwise)
    # V(d) = -g * K(d)
    V12 = -COUPLING_STRENGTH * K(abs(i - j))
    V23 = -COUPLING_STRENGTH * K(abs(j - k))
    V31 = -COUPLING_STRENGTH * K(abs(k - i))
    
    # Hopping / Tunneling?
    # FIN Theory: transition between octaves via K(d)
    # <i'|H|i> ~ -K(|i'-i|)
    
    # Diagonal term: Potential + Mass
    diag_H[idx] = E_kin + V12 + V23 + V31
    
    # Off-diagonal hopping (Tunneling) - essential for binding dynamics?
    # Let's add hopping for each particle individually
    # Particle 1 hops i -> i'
    # For now, let's stick to POTENTIAL binding first (Static Model)
    # If we add hopping, matrix becomes dense.
    
    # Actually, proper H_vac in QW-619 was off-diagonal coupling.
    # Let's include simplified hopping: only nearest octaves?
    # Or full K(d) hopping.
    
    # To keep it computable for 1728 states, let's limit hopping to nearest neighbors
    # or just analyze the Potential Landscape (Static Binding)
    pass

# Interpretation: Finding the lowest energy CONFIGURATION (i,j,k)
# This corresponds to the ground state if kinetic energy is small.

min_energy = np.min(diag_H)
min_idx = np.argmin(diag_H)
best_config = basis[min_idx]

print(f"\nLowest Potential Energy Configuration: {best_config}")
print(f"Energy: {min_energy:.4f}")

# Compare to isolated
E_isolated = 3 * 1.0 # 3 masses
print(f"Isolated Energy: {E_isolated:.4f}")
print(f"Binding Energy (Static): {min_energy - E_isolated:.4f}")

# ============================================================================
# FULL QUANTUM SOLUTION (Small Subspace)
# ============================================================================
# Solve full H with hopping for a subset of octaves relevant to best_config
# E.g. restrict to octaves 0-5

print("\n" + "="*80)
print("FULL QUANTUM SOLUTION (Subspace 6^3 = 216 states)")
print("="*80)

N_SUB = 6
H_sub = np.zeros((N_SUB**3, N_SUB**3))
basis_sub = list(product(range(N_SUB), repeat=3))

for idx_row, (i1, j1, k1) in enumerate(basis_sub):
    for idx_col, (i2, j2, k2) in enumerate(basis_sub):
        
        val = 0.0
        
        # Diagonal elements (if same state)
        if idx_row == idx_col:
            V12 = -COUPLING_STRENGTH * K(abs(i1 - j1))
            V23 = -COUPLING_STRENGTH * K(abs(j1 - k1))
            V31 = -COUPLING_STRENGTH * K(abs(k1 - i1))
            val += (H1[i1] + H1[j1] + H1[k1]) + (V12 + V23 + V31)
            
        # Off-diagonal (Hopping)
        # Particle 1 hops: (i1 != i2) but (j1==j2, k1==k2)
        if (j1==j2 and k1==k2) and (i1 != i2):
            val += -K(abs(i1 - i2))
            
        # Particle 2 hops
        if (i1==i2 and k1==k2) and (j1 != j2):
            val += -K(abs(j1 - j2))
            
        # Particle 3 hops
        if (i1==i2 and j1==j2) and (k1 != k2):
            val += -K(abs(k1 - k2))
            
        H_sub[idx_row, idx_col] = val

# Diagonalize
evals_sub, evecs_sub = eigh(H_sub)

E_proton_ground = evals_sub[0]
E_isolated_sub = 3 * (-K(0)) # Roughly ground state of single particle in hopping H (approx)
# Or calculate single particle ground state in N_SUB
H1_single = np.zeros((N_SUB, N_SUB))
for i in range(N_SUB):
    for j in range(N_SUB):
        if i==j: H1_single[i,j] = H1[i]
        else: H1_single[i,j] = -K(abs(i-j))
e1, _ = eigh(H1_single)
E_single = e1[0]

print(f"Single Particle Ground: {E_single:.4f}")
print(f"3 Isolated Particles: {3*E_single:.4f}")
print(f"Proton Ground State: {E_proton_ground:.4f}")
Binding_Energy_Proton = E_proton_ground - 3*E_single

print(f"\nProton Binding Energy: {Binding_Energy_Proton:.4f}")

if Binding_Energy_Proton < -1.0:
    print("\n✅ STABLE PROTON STRUCTURE")
    print("   3-mode bound state exists!")
    verdict = "stable"
else:
    print("\n❌ UNSTABLE")
    verdict = "unstable"

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw620_proton_structure.md"
with open(report_path, "w") as f:
    f.write("# Raport QW-620: Proton Structure (3-Mode)\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Test:** Stabilność układu 3 modów w sieci oktaw\n\n")
    
    f.write("## Wyniki\n")
    f.write(f"- Isolated Energy: {3*E_single:.4f}\n")
    f.write(f"- Proton Energy: {E_proton_ground:.4f}\n")
    f.write(f"- **Binding Energy:** {Binding_Energy_Proton:.4f}\n\n")
    
    if verdict == "stable":
        f.write("### ✅ STABILNY PROTON\n")
        f.write("Potwierdzono istnienie stanu związanego 3 modów (uud).\n")
    else:
        f.write("### ❌ BRAK WIĄZANIA\n")

print("Report saved.")
