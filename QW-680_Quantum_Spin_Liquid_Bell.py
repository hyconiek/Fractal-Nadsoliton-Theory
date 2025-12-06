#!/usr/bin/env python3
"""
QW-680: QUANTUM SPIN LIQUID BELL TEST
======================================
Purpose: Test Bell inequality in a FRUSTRATED spin network to achieve S > 2.0
         
Key Insight from PLAN_BADAWCZY:
  - QW-237, QW-545: S < 2 because FIN uses SCALAR fields (classical LHV)
  - QW-592, QW-678: S > 2 but use STANDARD QM (Pauli singlet), not FIN
  - Need: Bell violation from FIN NETWORK ITSELF
  
Approach:
  1. Use triangular lattice (geometric frustration)
  2. Heisenberg antiferromagnetic coupling (J < 0)
  3. Tune temperature near critical point (T_c)
  4. Measure CHSH from network correlations
  
Theory:
  In Quantum Spin Liquid (QSL), spins don't freeze - they fluctuate quantum mechanically.
  This allows long-range entanglement even at T > 0.
  
Output: RAPORT_QW680_BELL_QSL.md
"""

import numpy as np
from scipy.linalg import expm
from scipy.spatial import Delaunay
import datetime

print("="*80)
print("QW-680: QUANTUM SPIN LIQUID BELL TEST")
print("="*80)
print("Testing Bell inequality with frustrated triangular spin network")
print()

# --- Parameters ---
N_NODES = 50  # Smaller for numerical accuracy
J_COUPLING = -1.0  # ANTIFERROMAGNETIC (negative = frustration)
DT = 0.02
STEPS = 500
BETA = 0.1  # Inverse temperature (higher = colder)

# Pauli Matrices
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

REPORT_FILE = "RAPORT_QW680_BELL_QSL.md"

# --- CREATE TRIANGULAR LATTICE ---
print("1. Creating frustrated triangular lattice...")

# Random positions in 2D, then Delaunay triangulation
np.random.seed(680)
positions = np.random.randn(N_NODES, 2) * 3.0

# Delaunay gives triangular mesh
tri = Delaunay(positions)

# Build adjacency matrix from triangulation edges
adj_matrix = np.zeros((N_NODES, N_NODES))
for simplex in tri.simplices:
    for i in range(3):
        for j in range(i+1, 3):
            a, b = simplex[i], simplex[j]
            adj_matrix[a, b] = 1.0
            adj_matrix[b, a] = 1.0

avg_degree = np.sum(adj_matrix > 0) / N_NODES
print(f"   Nodes: {N_NODES}, Avg degree: {avg_degree:.1f}")
print(f"   Coupling J = {J_COUPLING} (antiferromagnetic = frustrated)")

# --- INITIALIZE SPINORS ---
print("\n2. Initializing random spinors...")
spinors = np.random.randn(N_NODES, 2) + 1j * np.random.randn(N_NODES, 2)
spinors /= np.linalg.norm(spinors, axis=1, keepdims=True)

def get_spin_vector(psi):
    """Calculate <psi|sigma|psi>"""
    sx = np.real(psi.conj().T @ sigma_x @ psi)
    sy = np.real(psi.conj().T @ sigma_y @ psi)
    sz = np.real(psi.conj().T @ sigma_z @ psi)
    return np.array([sx, sy, sz])

# --- EVOLVE WITH FRUSTRATED HEISENBERG ---
print("\n3. Evolving with frustrated Heisenberg dynamics...")
print(f"   Steps: {STEPS}, dt = {DT}")

spin_vectors = np.zeros((N_NODES, 3))
for i in range(N_NODES):
    spin_vectors[i] = get_spin_vector(spinors[i])

history_total_spin = []
history_correlations = []

for step in range(STEPS):
    new_spinors = np.zeros_like(spinors)
    
    for i in range(N_NODES):
        neighbors = np.where(adj_matrix[i] > 0)[0]
        
        if len(neighbors) == 0:
            new_spinors[i] = spinors[i]
            continue
        
        # Mean field from neighbors (frustrated = antiferromagnetic)
        B_eff = np.zeros(3)
        for j in neighbors:
            B_eff += adj_matrix[i, j] * spin_vectors[j]
        
        B_eff *= J_COUPLING  # Negative J = frustration
        
        # Local Hamiltonian H_i = - B_eff . sigma
        H_local = -(B_eff[0]*sigma_x + B_eff[1]*sigma_y + B_eff[2]*sigma_z)
        
        # Add thermal noise (simulate finite temperature)
        noise = np.random.randn(1) * np.sqrt(BETA) * 0.1
        H_noise = noise * sigma_z
        
        # Time evolution
        U = expm(-1j * (H_local + H_noise) * DT)
        new_spinors[i] = U @ spinors[i]
    
    spinors = new_spinors
    spinors /= np.linalg.norm(spinors, axis=1, keepdims=True)
    
    # Update spin vectors
    total_spin = np.zeros(3)
    for i in range(N_NODES):
        spin_vectors[i] = get_spin_vector(spinors[i])
        total_spin += spin_vectors[i]
    
    history_total_spin.append(np.linalg.norm(total_spin) / N_NODES)
    
    # Sample correlations periodically
    if step % 50 == 0:
        # Calculate nearest-neighbor correlations
        corr_sum = 0
        count = 0
        for i in range(N_NODES):
            for j in np.where(adj_matrix[i] > 0)[0]:
                if j > i:  # Avoid double counting
                    corr_sum += np.dot(spin_vectors[i], spin_vectors[j])
                    count += 1
        if count > 0:
            history_correlations.append(corr_sum / count)

final_magnetization = history_total_spin[-1]
print(f"   Final magnetization: {final_magnetization:.4f}")

# --- CHSH MEASUREMENT ---
print("\n4. Measuring CHSH parameter from FIN network...")

def measure_fin_correlation(spinors, adj, axis_A, axis_B):
    """
    Measure <S_A · n_A> <S_B · n_B> correlation between connected pairs.
    This is the FIN network version - NOT standard QM tensor product!
    """
    correlations = []
    
    for i in range(len(spinors)):
        for j in np.where(adj[i] > 0)[0]:
            if j > i:  # Avoid double counting
                # Measure spin_i along axis_A
                op_A = axis_A[0]*sigma_x + axis_A[1]*sigma_y + axis_A[2]*sigma_z
                s_A = np.real(spinors[i].conj().T @ op_A @ spinors[i])
                
                # Measure spin_j along axis_B
                op_B = axis_B[0]*sigma_x + axis_B[1]*sigma_y + axis_B[2]*sigma_z
                s_B = np.real(spinors[j].conj().T @ op_B @ spinors[j])
                
                correlations.append(s_A * s_B)
    
    return np.mean(correlations) if correlations else 0

# Standard CHSH angles
theta_a = 0
theta_a_prime = np.pi / 4
theta_b = np.pi / 8
theta_b_prime = -np.pi / 8

axis_a = np.array([np.sin(theta_a), 0, np.cos(theta_a)])
axis_a_prime = np.array([np.sin(theta_a_prime), 0, np.cos(theta_a_prime)])
axis_b = np.array([np.sin(theta_b), 0, np.cos(theta_b)])
axis_b_prime = np.array([np.sin(theta_b_prime), 0, np.cos(theta_b_prime)])

# Measure correlations
E_ab = measure_fin_correlation(spinors, adj_matrix, axis_a, axis_b)
E_ab_prime = measure_fin_correlation(spinors, adj_matrix, axis_a, axis_b_prime)
E_a_prime_b = measure_fin_correlation(spinors, adj_matrix, axis_a_prime, axis_b)
E_a_prime_b_prime = measure_fin_correlation(spinors, adj_matrix, axis_a_prime, axis_b_prime)

# CHSH parameter
S_fin = np.abs(E_ab - E_ab_prime + E_a_prime_b + E_a_prime_b_prime)

print(f"   E(a,b) = {E_ab:.4f}")
print(f"   E(a,b') = {E_ab_prime:.4f}")
print(f"   E(a',b) = {E_a_prime_b:.4f}")
print(f"   E(a',b') = {E_a_prime_b_prime:.4f}")
print(f"\n   CHSH Parameter S = {S_fin:.4f}")

if S_fin > 2.0:
    bell_status = "✅ QUANTUM VIOLATION! S > 2.0"
else:
    bell_status = f"❌ CLASSICAL: S = {S_fin:.4f} ≤ 2.0"

print(f"\n   {bell_status}")

# --- COMPARE WITH FERROMAGNETIC ---
print("\n5. Control test: Ferromagnetic coupling (J > 0)...")

# Reset spinors
spinors_ferro = np.random.randn(N_NODES, 2) + 1j * np.random.randn(N_NODES, 2)
spinors_ferro /= np.linalg.norm(spinors_ferro, axis=1, keepdims=True)

J_FERRO = +1.0  # Ferromagnetic

for step in range(STEPS):
    new_spinors = np.zeros_like(spinors_ferro)
    spin_vecs = np.array([get_spin_vector(s) for s in spinors_ferro])
    
    for i in range(N_NODES):
        neighbors = np.where(adj_matrix[i] > 0)[0]
        if len(neighbors) == 0:
            new_spinors[i] = spinors_ferro[i]
            continue
        
        B_eff = np.sum([spin_vecs[j] for j in neighbors], axis=0) * J_FERRO
        H_local = -(B_eff[0]*sigma_x + B_eff[1]*sigma_y + B_eff[2]*sigma_z)
        U = expm(-1j * H_local * DT)
        new_spinors[i] = U @ spinors_ferro[i]
    
    spinors_ferro = new_spinors
    spinors_ferro /= np.linalg.norm(spinors_ferro, axis=1, keepdims=True)

E_ab_ferro = measure_fin_correlation(spinors_ferro, adj_matrix, axis_a, axis_b)
E_ab_prime_ferro = measure_fin_correlation(spinors_ferro, adj_matrix, axis_a, axis_b_prime)
E_a_prime_b_ferro = measure_fin_correlation(spinors_ferro, adj_matrix, axis_a_prime, axis_b)
E_a_prime_b_prime_ferro = measure_fin_correlation(spinors_ferro, adj_matrix, axis_a_prime, axis_b_prime)

S_ferro = np.abs(E_ab_ferro - E_ab_prime_ferro + E_a_prime_b_ferro + E_a_prime_b_prime_ferro)
print(f"   S (ferromagnetic) = {S_ferro:.4f}")

# --- WRITE REPORT ---
print("\n6. Writing report...")

with open(REPORT_FILE, "w") as f:
    f.write("# RAPORT: QW-680 QUANTUM SPIN LIQUID BELL TEST\n")
    f.write(f"**Data:** {datetime.datetime.now()}\n")
    f.write("**Cel:** Uzyskać S > 2.0 z SIECI FIN (nie standard QM)\n\n")
    
    f.write("## 1. METODOLOGIA\n\n")
    f.write("- **Sieć:** Delaunay triangulation (N=50 nodes)\n")
    f.write(f"- **Sprzężenie:** J = {J_COUPLING} (antiferromagnetyczne = frustracja)\n")
    f.write(f"- **Dynamika:** Heisenberg mean-field + thermal noise\n")
    f.write(f"- **Ewolucja:** {STEPS} steps, dt = {DT}\n\n")
    
    f.write("## 2. WYNIKI CHSH\n\n")
    f.write("| Test | Coupling | S | Status |\n")
    f.write("|------|----------|---|--------|\n")
    f.write(f"| Frustrated (QSL) | J = -1 | {S_fin:.4f} | {'✅' if S_fin > 2 else '❌'} |\n")
    f.write(f"| Ferromagnetic | J = +1 | {S_ferro:.4f} | {'✅' if S_ferro > 2 else '❌'} |\n")
    f.write(f"| Classical bound | - | 2.0 | - |\n")
    f.write(f"| Tsirelson (QM max) | - | 2.83 | - |\n\n")
    
    f.write("## 3. KORELACJE\n\n")
    f.write("### Frustrated Network:\n")
    f.write(f"- E(a,b) = {E_ab:.4f}\n")
    f.write(f"- E(a,b') = {E_ab_prime:.4f}\n")
    f.write(f"- E(a',b) = {E_a_prime_b:.4f}\n")
    f.write(f"- E(a',b') = {E_a_prime_b_prime:.4f}\n\n")
    
    f.write("### Ferromagnetic Network:\n")
    f.write(f"- E(a,b) = {E_ab_ferro:.4f}\n")
    f.write(f"- E(a,b') = {E_ab_prime_ferro:.4f}\n")
    f.write(f"- E(a',b) = {E_a_prime_b_ferro:.4f}\n")
    f.write(f"- E(a',b') = {E_a_prime_b_prime_ferro:.4f}\n\n")
    
    f.write("## 4. ANALIZA\n\n")
    if S_fin > 2.0:
        f.write("### ✅ SUKCES: FRUSTRACJA DAJE ŁAMANIE BELLA!\n\n")
        f.write("Quantum Spin Liquid phase generuje prawdziwe splątanie kwantowe.\n")
        f.write("To dowodzi że sieć FIN może być kwantowa przy odpowiedniej topologii.\n")
    else:
        f.write("### ❌ PORAŻKA: FRUSTRACJA NIE WYSTARCZA\n\n")
        f.write(f"S = {S_fin:.4f} < 2.0 (granica klasyczna)\n\n")
        f.write("**Przyczyny:**\n")
        f.write("1. Mean-field approximation traci korelacje kwantowe\n")
        f.write("2. Thermal noise niszczy splątanie\n")
        f.write("3. Potrzeba pełnej kwantowej dynamiki (nie mean-field)\n\n")
        f.write("**Rekomendacja:**\n")
        f.write("- Użyć Exact Diagonalization małego klastra (N~16)\n")
        f.write("- Lub DMRG dla większych układów\n")
        f.write("- Lub wprowadzić Resonating Valence Bond (RVB) state\n")
    
    f.write("\n## 5. WNIOSEK\n\n")
    if S_fin > 2.0:
        f.write("**STATUS:** ✅ FIN może być KWANTOWA przy frustracji geometrycznej.\n")
        f.write("H12 (Partial Quantumness) → FULL QUANTUMNESS\n")
    else:
        f.write(f"**STATUS:** ❌ Mean-field FIN pozostaje klasyczna (S = {S_fin:.4f})\n")
        f.write("H12 pozostaje 'Partial' - potrzeba pełnej kwantyzacji.\n")

print(f"   Report saved to: {REPORT_FILE}")

print("\n" + "="*80)
print("QW-680 COMPLETE")
print("="*80)
print(f"\nFINAL RESULT: S = {S_fin:.4f}")
print(bell_status)
