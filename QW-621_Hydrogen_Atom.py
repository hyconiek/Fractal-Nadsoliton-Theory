#!/usr/bin/env python3
# QW-621: HYDROGEN ATOM (PROTON + ELECTRON)
# Purpose: Test binding of Electron (Octave 1) and Proton (Octave 7)
# Model: 2-Particle Tensor Product Space in 12-Octave Network
# Interaction: Attractive V(d) derived from coupling
# Question: Do we get a stable Hydrogen atom? What is the binding energy?
# Date: 2025-12-05

import numpy as np
from scipy.linalg import eigh

print("="*80)
print("QW-621: HYDROGEN ATOM (p + e)")
print("="*80)
print("Test: Wiązanie Elektronu (O1) z Protonem (O7)")
print("Hypothesis: Strong stable binding, localized relative state")
print("="*80)

# ============================================================================
# PARAMETERS
# ============================================================================
N_OCTAVES = 12
ALPHA_GEO = 4 * np.log(2)
BETA_TORS = 0.01

# Coupling strength
g_coupling = 5.0

print(f"Network: {N_OCTAVES} octaves")
print("-" * 40)

# ============================================================================
# HAMILTONIAN
# ============================================================================

def K(d):
    return ALPHA_GEO / (1 + BETA_TORS * d)

# Single particle Hamiltonians (Diagonal Mass/Energy terms)
# Electron in Octave 1 (low mass effective?)
# Proton in Octave 7 (high mass effective?)
# In FIN, "Octave" is a mode index. Mass emerges from topology/winding.
# Here we assign intrinsic energy to the mode to simulate mass difference.

H_electron = np.zeros((N_OCTAVES, N_OCTAVES))
H_proton = np.zeros((N_OCTAVES, N_OCTAVES))

for i in range(N_OCTAVES):
    # Electron prefers Octave 1 (Potential well)
    H_electron[i, i] = 1.0 + 0.5 * (i - 1)**2  # Harmonic trap around 1
    
    # Proton prefers Octave 7 (Potential well)
    H_proton[i, i] = 200.0 + 10.0 * (i - 7)**2 # Stiffer trap (heavier) around 7

# Tensor Product Basis: |i>|j> (Electron at i, Proton at j)
dim = N_OCTAVES * N_OCTAVES
H_total = np.zeros((dim, dim))

print("Building 2-Particle Hamiltonian...")

for i_e in range(N_OCTAVES):
    for i_p in range(N_OCTAVES):
        idx = i_e * N_OCTAVES + i_p
        
        # Kinetic / Potential Energy of individual particles
        E_site = H_electron[i_e, i_e] + H_proton[i_p, i_p]
        
        # Interaction V_ep depends on distance |i_e - i_p|
        # Attractive force mediated by network coupling
        interaction = -g_coupling * K(abs(i_e - i_p))
        
        H_total[idx, idx] = E_site + interaction
        
        # Hopping (Tunneling)
        # Allow electron to hop (light)
        # Allow proton to hop (heavy -> small hopping)
        
        hopping_e = 1.0 # Amplitude
        hopping_p = 0.05 # Much heavier
        
        # Add off-diagonal terms
        # Electron hop i_e -> j_e
        for j_e in range(N_OCTAVES):
            if i_e != j_e:
                idx_target = j_e * N_OCTAVES + i_p
                dist_hop = abs(i_e - j_e)
                # Hopping decreases with distance
                amp = -hopping_e * np.exp(-dist_hop)
                H_total[idx, idx_target] += amp
                
        # Proton hop i_p -> j_p
        for j_p in range(N_OCTAVES):
            if i_p != j_p:
                idx_target = i_e * N_OCTAVES + j_p
                dist_hop = abs(i_p - j_p)
                amp = -hopping_p * np.exp(-dist_hop)
                H_total[idx, idx_target] += amp

# ============================================================================
# SOLVE
# ============================================================================

print("Diagonalizing...")
evals, evecs = eigh(H_total)

E_ground = evals[0]
psi_ground = evecs[:, 0]

# ============================================================================
# ANALYSIS
# ============================================================================

# Isolated Energies (Approximate lowest energy of H_electron + H_proton without interaction)
E_iso_e = H_electron[1, 1] - hopping_e # Approx ground
E_iso_p = H_proton[7, 7] - hopping_p # Approx ground
E_isolated = E_iso_e + E_iso_p

# Re-calculate accurate isolated
H_e_iso = np.zeros((N_OCTAVES, N_OCTAVES))
for i in range(N_OCTAVES):
    H_e_iso[i,i] = H_electron[i,i]
    for j in range(N_OCTAVES):
        if i!=j: H_e_iso[i,j] = -hopping_e * np.exp(-abs(i-j))
e_val, _ = eigh(H_e_iso)
real_E_iso_e = e_val[0]

H_p_iso = np.zeros((N_OCTAVES, N_OCTAVES))
for i in range(N_OCTAVES):
    H_p_iso[i,i] = H_proton[i,i]
    for j in range(N_OCTAVES):
        if i!=j: H_p_iso[i,j] = -hopping_p * np.exp(-abs(i-j))
p_val, _ = eigh(H_p_iso)
real_E_iso_p = p_val[0]

Real_Isolated_Total = real_E_iso_e + real_E_iso_p

print(f"\nElectron Isolated E: {real_E_iso_e:.4f}")
print(f"Proton Isolated E: {real_E_iso_p:.4f}")
print(f"Total Isolated E: {Real_Isolated_Total:.4f}")
print(f"Hydrogen Ground E: {E_ground:.4f}")

Binding_Energy = E_ground - Real_Isolated_Total
print(f"\nBinding Energy: {Binding_Energy:.4f}")

# Structure analysis
prob_density = np.abs(psi_ground)**2
max_idx = np.argmax(prob_density)
i_e_max = max_idx // N_OCTAVES
i_p_max = max_idx % N_OCTAVES

print(f"\nMost probable configuration:")
print(f"  Electron Octave: {i_e_max}")
print(f"  Proton Octave: {i_p_max}")
print(f"  Separation: {abs(i_e_max - i_p_max)} octaves")

# Check Physics
if Binding_Energy < -1.0:
    print("\n✅ STABLE HYDROGEN ATOM")
    print("   Electron bound to Proton potential well.")
    verdict = "stable"
    
    # Check if electron is pulled towards proton
    if i_e_max != 1:
        print(f"   Note: Electron shifted from O1 to O{i_e_max} due to attraction!")
    else:
        print("   Electron remains in O1 (Balance of attraction vs trap)")
else:
    print("\n❌ UNSTABLE")
    verdict = "unstable"

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw621_hydrogen_atom.md"
with open(report_path, "w") as f:
    f.write("# Raport QW-621: Hydrogen Atom (p+e)\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Test:** Czy p(O7) i e(O1) tworzą stabilny atom?\n\n")
    
    f.write("## Wyniki\n")
    f.write(f"- Isolated Energy: {Real_Isolated_Total:.4f}\n")
    f.write(f"- Hydrogen Energy: {E_ground:.4f}\n")
    f.write(f"- **Binding Energy:** {Binding_Energy:.4f}\n\n")
    
    f.write("## Struktura\n")
    f.write(f"- Electron Octave: {i_e_max}\n")
    f.write(f"- Proton Octave: {i_p_max}\n")
    f.write(f"- Separation: {abs(i_e_max - i_p_max)}\n\n")
    
    if verdict == "stable":
        f.write("### ✅ STABILNY WODÓR\n")
        f.write("Potwierdzono stan związany protonu i elektronu w sieci oktaw.\n")
    else:
        f.write("### ❌ BRAK WIĄZANIA\n")

print("Report saved.")
