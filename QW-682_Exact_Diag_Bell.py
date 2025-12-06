#!/usr/bin/env python3
"""
QW-682: EXACT DIAGONALIZATION BELL TEST
=========================================
Purpose: Test Bell inequality with FULL quantum mechanics (not mean-field)

Key difference from QW-680:
  - QW-680: Mean-field approximation traci korelacje kwantowe
  - QW-682: Exact Diagonalization zachowuje pełne splątanie

Method:
  1. Small cluster (N=8) - exact diagonalizable
  2. Heisenberg antiferromagnet on triangle
  3. Full Hilbert space H = ⊗_i C² (2^N dimensional)
  4. Find ground state exactly
  5. Measure CHSH from quantum correlations

Output: RAPORT_QW682_EXACT_DIAG.md
"""

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import kron, eye
from scipy.sparse.linalg import eigsh
import datetime

print("="*80)
print("QW-682: EXACT DIAGONALIZATION BELL TEST")
print("="*80)
print("Testing Bell inequality with FULL quantum mechanics (Exact Diagonalization)")
print()

# --- Parameters ---
N_SPINS = 8  # 2^8 = 256 dimensional Hilbert space (manageable)
J_COUPLING = -1.0  # Antiferromagnetic (frustration source)

# Pauli Matrices (sparse)
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
sigma_0 = np.eye(2, dtype=complex)

REPORT_FILE = "RAPORT_QW682_EXACT_DIAG.md"

# --- BUILD FULL HAMILTONIAN ---
print("1. Building exact Hamiltonian...")
print(f"   N = {N_SPINS} spins → Hilbert space dim = {2**N_SPINS}")

def build_spin_operator(op, site, n_sites):
    """Build full operator for single site: I ⊗ ... ⊗ op ⊗ ... ⊗ I"""
    result = 1.0
    for i in range(n_sites):
        if i == site:
            result = np.kron(result, op)
        else:
            result = np.kron(result, sigma_0)
    return result

def build_heisenberg_hamiltonian(n_sites, J, topology='ring'):
    """
    Build Heisenberg Hamiltonian: H = J Σ_{<i,j>} S_i · S_j
    
    topology:
      - 'ring': nearest neighbor only
      - 'triangle': triangular lattice (frustration)
    """
    dim = 2**n_sites
    H = np.zeros((dim, dim), dtype=complex)
    
    # Build S^x, S^y, S^z for each site
    Sx = [build_spin_operator(0.5 * sigma_x, i, n_sites) for i in range(n_sites)]
    Sy = [build_spin_operator(0.5 * sigma_y, i, n_sites) for i in range(n_sites)]
    Sz = [build_spin_operator(0.5 * sigma_z, i, n_sites) for i in range(n_sites)]
    
    # Build bonds based on topology
    if topology == 'ring':
        bonds = [(i, (i+1) % n_sites) for i in range(n_sites)]
    elif topology == 'triangle':
        # Triangular lattice: add frustration with diagonal bonds
        bonds = [(i, (i+1) % n_sites) for i in range(n_sites)]
        # Add some long-range bonds for frustration
        for i in range(0, n_sites, 2):
            bonds.append((i, (i+2) % n_sites))
    
    # Add Heisenberg interaction for each bond
    for i, j in bonds:
        # S_i · S_j = Sx_i Sx_j + Sy_i Sy_j + Sz_i Sz_j
        H += J * (Sx[i] @ Sx[j] + Sy[i] @ Sy[j] + Sz[i] @ Sz[j])
    
    return H, Sx, Sy, Sz, bonds

# Build Hamiltonian
H, Sx, Sy, Sz, bonds = build_heisenberg_hamiltonian(N_SPINS, J_COUPLING, 'triangle')

print(f"   Topology: triangular (frustrated)")
print(f"   Number of bonds: {len(bonds)}")
print(f"   J = {J_COUPLING} (antiferromagnetic)")

# --- FIND GROUND STATE ---
print("\n2. Finding ground state (exact diagonalization)...")

eigenvalues, eigenvectors = eigh(H)
ground_state = eigenvectors[:, 0]
E_ground = eigenvalues[0]
E_gap = eigenvalues[1] - eigenvalues[0]

print(f"   Ground state energy: E_0 = {E_ground:.6f}")
print(f"   Energy gap: ΔE = {E_gap:.6f}")
print(f"   Ground state norm: {np.linalg.norm(ground_state):.6f}")

# --- MEASURE CHSH ---
print("\n3. Measuring CHSH parameter...")

def measure_correlation(psi, obs_A, obs_B):
    """
    Measure <psi| A ⊗ B |psi> exactly
    """
    AB = obs_A @ obs_B
    return np.real(psi.conj().T @ AB @ psi)

def spin_operator_at_angle(theta, phi, site, n_sites, Sx, Sy, Sz):
    """
    Spin measurement in direction (theta, phi)
    S_n = sin(theta)cos(phi) Sx + sin(theta)sin(phi) Sy + cos(theta) Sz
    """
    nx = np.sin(theta) * np.cos(phi)
    ny = np.sin(theta) * np.sin(phi)
    nz = np.cos(theta)
    return nx * Sx[site] + ny * Sy[site] + nz * Sz[site]

# Pick two distant spins (maximally separated in ring)
spin_A = 0
spin_B = N_SPINS // 2

print(f"   Measuring spins: A = {spin_A}, B = {spin_B}")

# CHSH optimal angles (in XZ plane, phi=0)
theta_a = 0
theta_a_prime = np.pi / 2
theta_b = np.pi / 4
theta_b_prime = 3 * np.pi / 4

# Build measurement operators
A = 2 * spin_operator_at_angle(theta_a, 0, spin_A, N_SPINS, Sx, Sy, Sz)
A_prime = 2 * spin_operator_at_angle(theta_a_prime, 0, spin_A, N_SPINS, Sx, Sy, Sz)
B = 2 * spin_operator_at_angle(theta_b, 0, spin_B, N_SPINS, Sx, Sy, Sz)
B_prime = 2 * spin_operator_at_angle(theta_b_prime, 0, spin_B, N_SPINS, Sx, Sy, Sz)

# Measure correlations
E_ab = measure_correlation(ground_state, A, B)
E_ab_prime = measure_correlation(ground_state, A, B_prime)
E_a_prime_b = measure_correlation(ground_state, A_prime, B)
E_a_prime_b_prime = measure_correlation(ground_state, A_prime, B_prime)

# CHSH parameter
S_chsh = np.abs(E_ab - E_ab_prime + E_a_prime_b + E_a_prime_b_prime)

print(f"\n   Correlations:")
print(f"   E(a,b) = {E_ab:.6f}")
print(f"   E(a,b') = {E_ab_prime:.6f}")
print(f"   E(a',b) = {E_a_prime_b:.6f}")
print(f"   E(a',b') = {E_a_prime_b_prime:.6f}")
print(f"\n   CHSH Parameter S = {S_chsh:.4f}")
print(f"   Classical bound: S ≤ 2.0")
print(f"   Tsirelson bound: S ≤ 2√2 ≈ 2.828")

if S_chsh > 2.0:
    bell_status = f"✅ QUANTUM! Bell violation S = {S_chsh:.4f} > 2.0"
    success = True
elif S_chsh > 1.5:
    bell_status = f"🟡 PARTIAL: S = {S_chsh:.4f} (approaching quantum regime)"
    success = False
else:
    bell_status = f"❌ CLASSICAL: S = {S_chsh:.4f} ≤ 2.0"
    success = False

print(f"\n   {bell_status}")

# --- TEST MULTIPLE SPIN PAIRS ---
print("\n4. Testing multiple spin pairs...")

results = []
for spin_A in range(N_SPINS // 2):
    for spin_B in range(N_SPINS // 2, N_SPINS):
        if spin_A == spin_B:
            continue
        
        A = 2 * spin_operator_at_angle(theta_a, 0, spin_A, N_SPINS, Sx, Sy, Sz)
        A_prime = 2 * spin_operator_at_angle(theta_a_prime, 0, spin_A, N_SPINS, Sx, Sy, Sz)
        B = 2 * spin_operator_at_angle(theta_b, 0, spin_B, N_SPINS, Sx, Sy, Sz)
        B_prime = 2 * spin_operator_at_angle(theta_b_prime, 0, spin_B, N_SPINS, Sx, Sy, Sz)
        
        E1 = measure_correlation(ground_state, A, B)
        E2 = measure_correlation(ground_state, A, B_prime)
        E3 = measure_correlation(ground_state, A_prime, B)
        E4 = measure_correlation(ground_state, A_prime, B_prime)
        
        S = np.abs(E1 - E2 + E3 + E4)
        results.append({'A': spin_A, 'B': spin_B, 'S': S})

# Find best pair
best = max(results, key=lambda x: x['S'])
print(f"   Best pair: A={best['A']}, B={best['B']} → S = {best['S']:.4f}")

# Statistics
S_values = [r['S'] for r in results]
print(f"   Max S: {max(S_values):.4f}")
print(f"   Mean S: {np.mean(S_values):.4f}")
print(f"   Min S: {min(S_values):.4f}")
num_violation = sum(1 for s in S_values if s > 2.0)
print(f"   Pairs with S > 2: {num_violation}/{len(results)}")

# --- ENTANGLEMENT ENTROPY ---
print("\n5. Computing entanglement entropy...")

# Partition system into A (first half) and B (second half)
n_A = N_SPINS // 2
dim_A = 2**n_A
dim_B = 2**(N_SPINS - n_A)

# Reshape ground state as matrix (Schmidt decomposition)
psi_matrix = ground_state.reshape(dim_A, dim_B)

# Reduced density matrix
rho_A = psi_matrix @ psi_matrix.conj().T

# Von Neumann entropy
eigenvalues_rho = np.linalg.eigvalsh(rho_A)
eigenvalues_rho = eigenvalues_rho[eigenvalues_rho > 1e-12]
S_entanglement = -np.sum(eigenvalues_rho * np.log(eigenvalues_rho))

print(f"   Partition: A = spins 0-{n_A-1}, B = spins {n_A}-{N_SPINS-1}")
print(f"   S_vN (entanglement entropy) = {S_entanglement:.4f}")

if S_entanglement > 0.5:
    print(f"   ✅ System is ENTANGLED!")
else:
    print(f"   ❌ Weak or no entanglement")

# --- WRITE REPORT ---
print("\n6. Writing report...")

with open(REPORT_FILE, "w") as f:
    f.write("# RAPORT: QW-682 EXACT DIAGONALIZATION BELL TEST\n")
    f.write(f"**Data:** {datetime.datetime.now()}\n")
    f.write("**Cel:** Test Bell inequality z PEŁNĄ mechaniką kwantową (nie mean-field)\n\n")
    
    f.write("## 1. METODOLOGIA\n\n")
    f.write(f"- **N spins:** {N_SPINS} (Hilbert space dim = {2**N_SPINS})\n")
    f.write(f"- **Topology:** Triangular (frustrated)\n")
    f.write(f"- **Coupling:** J = {J_COUPLING} (antiferromagnetic)\n")
    f.write(f"- **Method:** Full exact diagonalization (scipy.linalg.eigh)\n\n")
    
    f.write("## 2. GROUND STATE\n\n")
    f.write(f"- **Energy:** E_0 = {E_ground:.6f}\n")
    f.write(f"- **Gap:** ΔE = {E_gap:.6f}\n")
    f.write(f"- **Entanglement entropy:** S_vN = {S_entanglement:.4f}\n\n")
    
    f.write("## 3. CHSH TEST\n\n")
    f.write("| Para | S (CHSH) | Violation? |\n")
    f.write("|------|----------|------------|\n")
    f.write(f"| (0, {N_SPINS//2}) | {S_chsh:.4f} | {'✅' if S_chsh > 2 else '❌'} |\n")
    f.write(f"| Best ({best['A']}, {best['B']}) | {best['S']:.4f} | {'✅' if best['S'] > 2 else '❌'} |\n")
    f.write(f"| Mean | {np.mean(S_values):.4f} | - |\n\n")
    
    f.write("**Classical bound:** S ≤ 2.0\n")
    f.write("**Tsirelson bound:** S ≤ 2√2 ≈ 2.828\n\n")
    
    f.write("## 4. KORELACJE (para 0-4)\n\n")
    f.write(f"- E(a,b) = {E_ab:.6f}\n")
    f.write(f"- E(a,b') = {E_ab_prime:.6f}\n")
    f.write(f"- E(a',b) = {E_a_prime_b:.6f}\n")
    f.write(f"- E(a',b') = {E_a_prime_b_prime:.6f}\n\n")
    
    f.write("## 5. PODSUMOWANIE\n\n")
    if max(S_values) > 2.0:
        f.write("### ✅ SUKCES: BELL VIOLATION DETECTED!\n\n")
        f.write(f"Max S = {max(S_values):.4f} > 2.0\n\n")
        f.write("**Exact Diagonalization** zachowuje pełne korelacje kwantowe.\n")
        f.write("W przeciwieństwie do mean-field (QW-680), ED daje prawdziwe splątanie.\n\n")
        f.write("**Wniosek:** Sieć FIN z frustrated Heisenberg coupling **MOŻE** łamać Bell inequality.\n")
    else:
        f.write("### ❌ PORAŻKA: NO BELL VIOLATION\n\n")
        f.write(f"Max S = {max(S_values):.4f} ≤ 2.0\n\n")
        f.write("**Przyczyna:** Ground state nie jest maksymalnie splątany.\n")
        f.write("Frustracja tworzy spin liquid, ale nie singlet-like correlations.\n\n")
        f.write("**Rekomendacja:**\n")
        f.write("1. Zwiększyć N (więcej sprzężeń)\n")
        f.write("2. Próbować innych topologii (kagome, triangular 2D)\n")
        f.write("3. Dodać anisotropię (XXZ model)\n")
    
    f.write("\n## 6. PORÓWNANIE Z POPRZEDNIMI\n\n")
    f.write("| Badanie | Metoda | Max S | Status |\n")
    f.write("|---------|--------|-------|--------|\n")
    f.write("| QW-680 | Mean-field | 0.17 | ❌ |\n")
    f.write(f"| **QW-682** | **Exact Diag** | **{max(S_values):.2f}** | {'✅' if max(S_values) > 2 else '❌'} |\n")

print(f"   Report saved to: {REPORT_FILE}")

print("\n" + "="*80)
print("QW-682 COMPLETE")
print("="*80)
print(f"\nFINAL RESULT: Max S = {max(S_values):.4f}")
print(bell_status)
print(f"Entanglement entropy: S_vN = {S_entanglement:.4f}")
