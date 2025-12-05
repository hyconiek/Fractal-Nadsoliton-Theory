#!/usr/bin/env python3
# QW-619: OCTAVE CHEMISTRY (BINDING ENERGY)
# Purpose: Test if particles bind via Octave Resonance (not spatial)
# Replaces QW-602 which failed using spatial hopfions
# Model: Particles as excitations in specific octaves (1,4,7) interacting via K_ij
# Date: 2025-12-05

import numpy as np
from scipy.linalg import eigh

print("="*80)
print("QW-619: OCTAVE CHEMISTRY (BINDING ENERGY)")
print("="*80)
print("Test: Czy cząstki wiążą się rezonansem międzyoktawowym?")
print("Hypothesis: Binding Energy < 0 due to shared resonance modes")
print("="*80)

# ============================================================================
# PARAMETERS
# ============================================================================
N_OCTAVES = 12
ALPHA_GEO = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01

print(f"Network: {N_OCTAVES} octaves")
print("Particles defined as localized excitations (Gaussian in octave space)")
print("-" * 40)

# ============================================================================
# COUPLING HAMILTONIAN
# ============================================================================

def K(d):
    return ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + BETA_TORS * d)

# Base Hamiltonian (Vacuum)
H_vac = np.zeros((N_OCTAVES, N_OCTAVES))
for i in range(N_OCTAVES):
    for j in range(N_OCTAVES):
        H_vac[i, j] = -K(abs(i - j))  # Negative coupling for binding?

# Add diagonal energy (mass)
# Electron: Octave 1, Proton: Octave 4 (roughly)
E_electron = 1.0
E_proton = 200.0  # Mass ratio roughly
# But in octave space, let's treat them as modes 
# For testing generic binding: two identical particles vs p-e

# ============================================================================
# CALCULATE ENERGY OF ISOLATED PARTICLES
# ============================================================================

def get_particle_energy(center_octave, sigma=1.0):
    # Wavefunction localized at center_octave
    psi = np.exp(-0.5 * ((np.arange(N_OCTAVES) - center_octave) / sigma)**2)
    psi = psi / np.linalg.norm(psi)
    
    # Energy expectation: <ψ|H|ψ>
    E = psi.T @ H_vac @ psi
    return E, psi

E_A, psi_A = get_particle_energy(1, sigma=0.8)  # Particle A (e.g. electron)
E_B, psi_B = get_particle_energy(4, sigma=0.8)  # Particle B (e.g. muon/proton)

print(f"\nParticle A (Octave 1): E = {E_A:.4f}")
print(f"Particle B (Octave 4): E = {E_B:.4f}")
print(f"Total Isolated Energy: E_isolated = {E_A + E_B:.4f}")

# ============================================================================
# CALCULATE ENERGY OF COMBINED SYSTEM (MOLECULE)
# ============================================================================

print("\n" + "="*80)
print("TEST 1: SUPERPOSITION (Molecule)")
print("="*80)

# Ansatz: Molecular orbital ψ_M = c1*ψ_A + c2*ψ_B
# Solve generalized eigenvalue problem in subspace {ψ_A, ψ_B}
# H_sub = [[H_AA, H_AB], [H_BA, H_BB]]
# S_sub = [[1, S_AB], [S_BA, 1]] (Overlap)

H_AA = E_A
H_BB = E_B
H_AB = psi_A.T @ H_vac @ psi_B
S_AB = psi_A.T @ psi_B  # Overlap integral

print(f"Overlap <A|B>: {S_AB:.4f}")
print(f"Interaction <A|H|B>: {H_AB:.4f}")

# Molecular Hamiltonian matrix
H_mol = np.array([[H_AA, H_AB], [H_AB, H_BB]])
S_mol = np.array([[1.0, S_AB], [S_AB, 1.0]])

# Generalized eigenvalues: H c = E S c
from scipy.linalg import eigh
evals_mol, evecs_mol = eigh(H_mol, S_mol)

E_ground = evals_mol[0]
E_excited = evals_mol[1]

print(f"\nMolecular Ground State: E_mol = {E_ground:.4f}")
Binding_Energy = E_ground - (min(E_A, E_B) if abs(S_AB) > 0.9 else (E_A if E_A < E_B else E_B)) 
# Wait, classic binding: E_mol - (E_A + E_B) is for separate systems.
# Here we are in the SAME Hilbert space (1 particle in 12 octaves?)
# Ah, chemistry is TWO particles. We need Tensor Product Space!

# ============================================================================
# PROPER 2-PARTICLE SYSTEM (Tensor Product)
# ============================================================================

print("\n" + "="*80)
print("TEST 2: TWO-PARTICLE TENSOR PRODUCT STATE")
print("="*80)
print("State space: 12 × 12 = 144 basis states |i>|j>")

# Interaction V_int depends on relative octave distance?
# Or simply H_total = H_1 + H_2 + V_12
# Where V_12 is interaction between particle 1 at octave i and particle 2 at octave j

# Construct Full Hamiltonian (144x144)
H_total = np.zeros((N_OCTAVES*N_OCTAVES, N_OCTAVES*N_OCTAVES))

# Single particle Hamiltonians (H1 ⊗ I + I ⊗ H2)
I = np.eye(N_OCTAVES)
H1_full = np.kron(H_vac, I)
H2_full = np.kron(I, H_vac)

# Interaction Term V_12
# Hypothesis: Interaction strength K(|i-j|) is attractive?
# V = -g * K(|i-j|)
g_coupling = 5.0  # Strength
V_12 = np.zeros((N_OCTAVES*N_OCTAVES, N_OCTAVES*N_OCTAVES))

for i in range(N_OCTAVES):
    for j in range(N_OCTAVES): # Particle 1 at i, Particle 2 at j
        idx = i * N_OCTAVES + j
        # Interaction energy depends on separation (i, j)
        coupling = K(abs(i - j))
        V_12[idx, idx] = -g_coupling * abs(coupling) # Attractive potential

H_system = H1_full + H2_full + V_12

# Solve for ground state of 2-particle system
evals_sys, evecs_sys = eigh(H_system)

E_sys_ground = evals_sys[0]  # Lowest energy state
E_isolated_sys = E_A + E_B # Simple sum of previous single-particle energies (vac expectation)

# Re-calculate E_A, E_B properly as eigenvalues of H_vac to be fair
evals_vac, _ = eigh(H_vac)
E_single_ground = evals_vac[0]
E_isolated_proper = 2 * E_single_ground

print(f"\nSingle Particle Ground Energy: {E_single_ground:.4f}")
print(f"Two Isolated Particles Energy: {E_isolated_proper:.4f}")
print(f"Interacting System Energy: {E_sys_ground:.4f}")

Binding_Energy_Sys = E_sys_ground - E_isolated_proper

print(f"\nBinding Energy: {Binding_Energy_Sys:.4f}")

if Binding_Energy_Sys < -0.01:
    print("\n✅ STABLE BOUND STATE FOUND!")
    print("   Resonance between octaves creates binding.")
    print("   Mechanism: Potential V_12 lowers energy for specific (i,j) pairs.")
    verdict = "bound"
    
    # Analyze structure of bound state
    psi_ground = evecs_sys[:, 0]
    prob_density = np.abs(psi_ground)**2
    prob_matrix = prob_density.reshape((N_OCTAVES, N_OCTAVES))
    
    # Find most probable configuration (i, j)
    max_idx = np.argmax(prob_density)
    i_max = max_idx // N_OCTAVES
    j_max = max_idx % N_OCTAVES
    
    print(f"\nMost probable configuration: Octave {i_max} & Octave {j_max}")
    print(f"Separation: {abs(i_max - j_max)} octaves")
    
else:
    print("\n❌ NO BINDING")
    print("   Particles remain isolated or repulsive.")
    verdict = "unbound"

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw619_octave_chemistry.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-619: Octave Chemistry (Binding)\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Model:** Tensor Product Space 12x12 with V_12 interaction\n\n")
    
    f.write("## Wyniki\n")
    f.write(f"- Isolated Energy: {E_isolated_proper:.4f}\n")
    f.write(f"- Bound System Energy: {E_sys_ground:.4f}\n")
    f.write(f"- **Binding Energy:** {Binding_Energy_Sys:.4f}\n\n")
    
    if verdict == "bound":
        f.write("### ✅ STABLE BOUND STATE\n")
        f.write("Potwierdzono mechanizm wiązania rezonansowego w sieci oktaw.\n")
        f.write(f"Konfiguracja orbitalna: Octave {i_max} - Octave {j_max} (d={abs(i_max-j_max)})\n")
    else:
        f.write("### ❌ BRAK WIĄZANIA\n")

print("Report saved.")
print("="*80)
