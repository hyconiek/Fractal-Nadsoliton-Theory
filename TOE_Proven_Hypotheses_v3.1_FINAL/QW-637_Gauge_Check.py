#!/usr/bin/env python3
# QW-637: GAUGE INVARIANCE CHECK (THE FINAL CANNON)
# Purpose: Test if FIN Theory Hamiltonian is locally Gauge Invariant.
#          Scepticus: "If NO, you have no Photons/Charge."
# System: 1D Chain with Complex Hopping (Hopfion Model from QW-636).
# Date: 2025-12-05

import numpy as np
import matplotlib.pyplot as plt

print("="*80)
print("QW-637: GAUGE INVARIANCE CHECK")
print("="*80)

# System
N = 10
state = np.random.rand(N) + 1j * np.random.rand(N)
state /= np.linalg.norm(state)

# Hamiltonian (Hopfion Phase Model from QW-636)
# H = Sum t * c_i+ c_j + h.c.
# t_ij = |t| * exp(i * phi_ij)
# We assume fixed background t_ij (Static Lattice).

def energy(psi, phases_ij):
    E = 0.0
    for i in range(N):
        # Neighbor j = (i+1)%N
        j = (i+1)%N
        
        # Link Phase (fixed background or dynamical?)
        # Standard Matter Hamiltonian assumes fixed U_ij.
        # But we want to check invariance of the MATTER sector under Matter Gauge transform.
        # If we rotate psi_i -> psi_i * e^i alpha_i
        # The term psi_i^* U_ij psi_j transforms as e^-i alpha_i * U * e^i alpha_j
        # It is ONLY invariant if U_ij transforms too: U -> e^i alpha_i U e^-i alpha_j
        
        # In FIN Theory, do we have a U field?
        # The theory says "Space is Information". Couplings are geometrical.
        # If geometry is fixed, U is fixed.
        
        val = np.vdot(psi[i], psi[j]) # Hopping
        # Add phase from link
        phase = phases_ij[i]
        term = val * np.exp(1j * phase)
        E += term
        
    return - (E + np.conj(E)).real

# 1. Setup Static Background
link_phases = np.zeros(N) # Vacuum

E_original = energy(state, link_phases)
print(f"Original Energy: {E_original:.6f}")

# 2. Apply Local Gauge Transform (Random Phases for each node)
alphas = np.random.rand(N) * 2 * np.pi
state_transformed = np.zeros(N, dtype=complex)
for i in range(N):
    state_transformed[i] = state[i] * np.exp(1j * alphas[i])

# 3. Check Energy of Transformed State 
# (Without changing link variables - naive check)
E_gauge = energy(state_transformed, link_phases)
print(f"Transformed Energy: {E_gauge:.6f}")

diff = abs(E_original - E_gauge)
print(f"\nEnergy Difference: {diff:.6f}")

if diff > 1e-6:
    print("❌ GAUGE SYMMETRY BROKEN (Naive)")
    print("   The Hamiltonian with fixed links is NOT invariant.")
    print("   Scepticus is right: Matter coupling is sensitive to local phase.")
    print("   This implies Charge is NOT conserved locally unless...")
else:
    print("✅ GAUGE SYMMETRY HELD")

# 4. The Defense: Emergent Gauge Field?
# If we allow link_phases to update aka "The geometry adapts".
# Transform U_ij -> U_ij * exp(i * (alpha_i - alpha_j))
print("\nApplying Compensating Link Transform (Emergent Gauge Field)...")

link_phases_prime = np.zeros(N)
for i in range(N):
    j = (i+1)%N
    # U' = U * e^(i(alpha_i - alpha_j)) ?
    # Term is psi_i^* * U * psi_j
    # psi_i^* -> psi_i^* e^-i alpha_i
    # psi_j -> psi_j e^i alpha_j
    # Total phase shift: alpha_j - alpha_i
    # We need U to shift by alpha_i - alpha_j to cancel.
    link_phases_prime[i] = link_phases[i] + (alphas[i] - alphas[j])

E_invariant = energy(state_transformed, link_phases_prime)
print(f"Corrected Energy: {E_invariant:.6f}")
diff_inv = abs(E_original - E_invariant)

if diff_inv < 1e-6:
    print("✅ GAUGE SYMMETRY RESTORED by Dynamic Geometry.")
    print("   Conclusion: The 'Link Variables' (Geometry) MUST match the Matter Phase.")
    print("   This means Geometry IS the Gauge Field.")
    print("   Spin Waves in Geomery = Gauge Bosons.")

# ============================================================================
# REPORT
# ============================================================================
with open("raport_qw637_gauge.md", "w") as f:
    f.write("# Raport QW-637: Gauge Invariance\n")
    f.write(f"Energy Diff (Naive): {diff:.6f}\n")
    if diff > 1e-6:
        f.write("### ❌ SCEPTICUS MA RACJĘ (Częściowo)\n")
        f.write("Przy sztywnej geometrii nie ma symetrii cechowania.\n")
    
    f.write(f"Energy Diff (Dynamic Geometry): {diff_inv:.6f}\n")
    if diff_inv < 1e-6:
        f.write("### ✅ OBRONA: EMERGENCE\n")
        f.write("Jeśli geometria sieci (fazy wiązań) reaguje na fazę materii, symetria jest zachowana.\n")
        f.write("Oznacza to, że Geometria Sieci = Pole Cechowania ($A_\\mu$).\n")

print("Report saved.")
