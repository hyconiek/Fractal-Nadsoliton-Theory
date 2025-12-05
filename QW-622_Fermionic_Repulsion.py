#!/usr/bin/env python3
# QW-622: FERMIONIC REPULSION (SPIN GAP)
# Purpose: Test if Pauli Exclusion (anti-symmetry) allows stable binding
#          while preventing particles from occupying the same state.
# Model: 2 Fermions in 12-Octave Network.
#        Wavefunction must be Anti-Symmetric: psi(i,j) = -psi(j,i)
# Date: 2025-12-05

import numpy as np
from scipy.linalg import eigh

print("="*80)
print("QW-622: FERMIONIC REPULSION (SPIN GAP)")
print("="*80)
print("Test: Czy zakaz Pauliego (anty-symetria) pozwala na wiązanie?")
print("Hypothesis: Binding exists but excludes i=j states (Fermi hole)")
print("="*80)

# ============================================================================
# PARAMETERS
# ============================================================================
N_OCTAVES = 12
ALPHA_GEO = 4 * np.log(2)
BETA_TORS = 0.01
g_coupling = 5.0

print(f"Network: {N_OCTAVES} octaves")
print("Particles: 2 Identical Fermions (e.g. 2 Electrons)")
print("-" * 40)

# ============================================================================
# BASIS CONSTUCTION (FERMIONIC)
# ============================================================================

# For fermions, basis states are |i,j> with i < j.
# Total dim = N(N-1)/2. For N=12, dim=66.
# If we use full tensor product N*N=144, we must projector onto anti-symmetric subspace.

# Let's use the full N*N approach and check the symmetry of the ground state.
# But verify if Hamiltonian preserves symmetry.

def K(d):
    return ALPHA_GEO / (1 + BETA_TORS * d)

# Single particle Hamiltonian (e.g. 2 electrons in potential well)
H_single = np.zeros((N_OCTAVES, N_OCTAVES))
for i in range(N_OCTAVES):
    H_single[i, i] = 1.0 + 0.5 * (i - 1)**2  # Potential well at O1
    # Add hopping
    for j in range(N_OCTAVES):
        if i!=j: H_single[i,j] = -1.0 * np.exp(-abs(i-j))

# 2-Particle Hamiltonian (Distinguishable first)
dim = N_OCTAVES * N_OCTAVES
H_total = np.zeros((dim, dim))

for i in range(N_OCTAVES):
    for j in range(N_OCTAVES):
        idx = i * N_OCTAVES + j
        
        # Diagonal E
        E_site = H_single[i,i] + H_single[j,j]
        
        # Interaction V(|i-j|) - Repulsive or Attractive?
        # Electrons repel via Coulomb. But here we test "Binding" mechanism.
        # Let's assume we maintain the 'Octave Attraction' V = -g*K (e.g. pairing mechanism)
        # But distinguish it from Coulomb.
        # If we stick to QW-619/621 model, interaction was attractive.
        # Let's test attractive interaction with Fermi statistics (Cooper Pair analog).
        interaction = -g_coupling * K(abs(i - j))
        if i == j: 
            # On-site interaction. For fermions, usually U (Hubbard).
            # But strictly, two spinless fermions cannot be at same site i.
            # So interaction term for i=j is irrelevant if physics is correct.
            pass
            
        H_total[idx, idx] = E_site + interaction
        
        # Hopping P1
        for next_i in range(N_OCTAVES):
            if i != next_i:
                idx_target = next_i * N_OCTAVES + j
                H_total[idx, idx_target] += H_single[i, next_i]
                
        # Hopping P2
        for next_j in range(N_OCTAVES):
            if j != next_j:
                idx_target = i * N_OCTAVES + next_j
                H_total[idx, idx_target] += H_single[j, next_j]

# ============================================================================
# SYMMETRY PROJECTION
# ============================================================================

# We want Anti-Symmetric Eigenstates.
# H commutes with Exchange Operator P_ij?
# If H is symmetric in particles 1 and 2 (H(1,2) = H(2,1)), then eigenstates are either Sym or Anti-Sym.
# Our H is symmetric (identical particles).

print("Diagonalizing...")
evals, evecs = eigh(H_total)

print("Filtering for Anti-Symmetric States (Fermions)...")
fermionic_states = []

for k in range(len(evals)):
    psi = evecs[:, k]
    psi_matrix = psi.reshape((N_OCTAVES, N_OCTAVES))
    
    # Check symmetry: psi(j,i) = -psi(i,j)
    psi_transpose = psi_matrix.T
    
    # Measure Anti-Symmetry
    # If Anti: psi + psi.T = 0
    # If Sym: psi - psi.T = 0
    
    diff_anti = np.linalg.norm(psi_matrix + psi_transpose) # Should be 0 for Anti
    diff_sym = np.linalg.norm(psi_matrix - psi_transpose)  # Should be 0 for Sym
    
    if diff_anti < 1e-5:
        fermionic_states.append((evals[k], psi_matrix, "Fermion"))
    elif diff_sym < 1e-5:
        # fermionic_states.append((evals[k], psi_matrix, "Boson")) # Ignore bosons
        pass

if len(fermionic_states) > 0:
    E_ground_fermi = fermionic_states[0][0]
    psi_ground_fermi = fermionic_states[0][1]
    
    print(f"\nFermionic Ground State Energy: {E_ground_fermi:.4f}")
    
    # Check binding
    # Compare to 2 isolated fermions in lowest DIFFERENT states
    # Single particle energies
    e_single, _ = eigh(H_single)
    E1 = e_single[0]
    E2 = e_single[1] # Next lowest state
    E_isolated_fermi = E1 + E2
    
    print(f"Isolated Fermions (E1+E2): {E_isolated_fermi:.4f}")
    
    Binding_Energy = E_ground_fermi - E_isolated_fermi
    print(f"Binding Energy: {Binding_Energy:.4f}")
    
    # Analyze Structure (Fermi Hole)
    prob_matrix = np.abs(psi_ground_fermi)**2
    diag_prob = np.sum(np.diag(prob_matrix))
    
    print(f"\nProbability of being in same octave (i=j): {diag_prob:.6f}")
    
    if diag_prob < 1e-6:
        print("✅ PAULI EXCLUSION CONFIRMED")
        print("   Particles completely avoid same-state (Fermi hole).")
    else:
        print("❌ PAOLI VIOLATION")
    
    if Binding_Energy < -0.1:
        print("✅ STABLE FERMIONIC BINDING")
        print("   Cooper-like pairing in octave network confirmed.")
        verdict = "success"
    else:
        print("❌ NO BINDING for Fermions")
        verdict = "fail"
else:
    print("❌ No Fermionic States found?? Check Hamiltonian symmetry.")
    verdict = "error"

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw622_fermionic_repulsion.md"
with open(report_path, "w") as f:
    f.write("# Raport QW-622: Fermionic Repulsion (Spin Gap)\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Test:** Stabilność wiązania fermionów (anty-symetria)\n\n")
    
    if verdict == "success":
        f.write(f"## Wyniki\n")
        f.write(f"- Fermionic Ground Energy: {E_ground_fermi:.4f}\n")
        f.write(f"- Binding Energy: {Binding_Energy:.4f}\n")
        f.write(f"- Probability(i=j): {diag_prob:.6f} (Fermi Hole check)\n\n")
        f.write("### ✅ SUKCES\n")
        f.write("Model oktawowy obsługuje statystykę Fermiego.\n")
        f.write("Pauli Exclusion działa (dziura Fermiego na diagonali).\n")
        f.write("Wiązanie (Cooper-like) jest stabilne mimo odpychania Pauliego.\n")
    else:
        f.write("### ❌ PORAŻKA\n")

print("Report saved.")
