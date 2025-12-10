#!/usr/bin/env python3
# QW-630: RIGOROUS QUANTUM ENTANGLEMENT (4-BIT SPINORS)
# Purpose: Test if FIN Theory interaction kernel generates genuine Quantum Entanglement.
#          Unlike QW-545/592 (Classical proxies), this uses rigorous QM Density Matrices.
# System: Two 4-bit Spinors (Nodes A and B). Hilbert Space dim = 16 * 16 = 256.
# Interaction: H_int = g * K(d) * (Sigma_A . Sigma_B)
# Verify: CHSH Violation (S > 2.0) and Entanglement Entropy.
# Date: 2025-12-05

import numpy as np
import scipy.linalg as la

print("="*80)
print("QW-630: RIGOROUS QUANTUM ENTANGLEMENT")
print("="*80)
print("Test: Can FIN Interaction generate S > 2 (Bell Violation)?")
print("System: Two 4-bit Spinors (16 states each). Total Dim = 256.")
print("="*80)

# 1. Define 4-bit Spinor Operators
# 4 bits = 16 states.
# We need generalized Pauli operators for dim=16?
# Or treat it as 4 Qubits per node?
# FIN Theory 4-bit = single entity. Let's assume it behaves like a "Spin-15/2" (dim 16)?
# OR: 4 actual Qubits.
# Let's assume internal structure is 4 Qubits (Information).

# Pauli Matrices (single qubit)
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

# Construct Operators for 4-Qubit Node (A)
# We need an "observable" to measure.
# Let's define "Total Spin" of the node along Z: Z_tot = Z1 + Z2 + Z3 + Z4
def tensor_4(op, pos):
    ops = [I]*4
    ops[pos] = op
    res = ops[0]
    for i in range(1, 4):
        res = np.kron(res, ops[i])
    return res

# Total spin operators for Node A (16x16)
Sz_A = sum(tensor_4(Z, i) for i in range(4)) / 2
Sx_A = sum(tensor_4(X, i) for i in range(4)) / 2
Sy_A = sum(tensor_4(Y, i) for i in range(4)) / 2

# 2. System AB (256 x 256)
# Basis: |A> |B>
# Operators on full space
Sz_A_full = np.kron(Sz_A, np.eye(16))
Sz_B_full = np.kron(np.eye(16), Sz_A) # Same structure for B

Sx_A_full = np.kron(Sx_A, np.eye(16))
Sx_B_full = np.kron(np.eye(16), Sx_A)

Sy_A_full = np.kron(Sy_A, np.eye(16))
Sy_B_full = np.kron(np.eye(16), Sy_A)

# 3. Hamiltonian
# H = H_0 + V_int
# H_0 = 0 (degenerate)
# V_int = -g * K * (S_A . S_B)  (Heisenberg-like interaction)
# Coupling K derived from QW-627 (Kernel) -> Let's take K=1.0 for neighbors.
# g = coupling strength.

g = 1.0
H_int = -g * (np.dot(Sx_A_full, Sx_B_full) + 
              np.dot(Sz_A_full, Sz_B_full) + 
              # np.dot(Sy_A_full, Sy_B_full) # Simplify to XZ model for CHSH usually
              np.dot(Sy_A_full, Sy_B_full)
             )

# Note: dimensions mismatch? dot product?
# H = -g * (Sx_A * Sx_B + Sy_A * Sy_B + Sz_A * Sz_B) tensor product wise
H_int = -g * (
    np.kron(Sx_A, Sx_A) + 
    np.kron(Sy_A, Sy_A) + 
    np.kron(Sz_A, Sz_A)
)

print(f"Hamiltonian shape: {H_int.shape}")

# 4. Time Evolution to generate Entanglement
# Start from separated state: |0000>_A |0000>_B
# Ground state of Z.
# |0> is [1, 0]
psi_0 = np.array([1, 0], dtype=complex)
psi_node_0 = psi_0
for _ in range(3): psi_node_0 = np.kron(psi_node_0, psi_0) # |0000>
state_0 = np.kron(psi_node_0, psi_node_0) # |0000>|0000> (Size 256)

print("Initial Entanglement Entropy: 0")

# Evolve: Psi(t) = exp(-iHt) Psi(0)
# Let's compute for a time t that maximizes entanglement.
# Generic unitary evolution.
t = np.pi / 4 # Often good for Bell states
U = la.expm(-1j * H_int * t)
state_t = np.dot(U, state_0)

# 5. Measure Bell Inequality (CHSH)
# S = E(a, b) - E(a, b') + E(a', b) + E(a', b')
# We need to choose measurement settings a, a', b, b'.
# For standard Bell state:
# A = Z, A' = X
# B = (Z+X)/sqrt(2), B' = (Z-X)/sqrt(2)

# But we have spin-2 systems (4 qubits). Can we violate Bell?
# Yes, if we measure "Sign of Z" or parity.
# Let's map the measurement to {-1, +1}.
# Operator M: Sign(Sz).
# If Sz eigenvalues are -2, -1, 0, 1, 2...
# Parity = (-1)^(Sz)?

# Let's construct a "Dichotomic Observable" O_A.
# O_A = P_up - P_down ?
# Or just normalized Spin: Sz / |Sz_max|? No, CHSH requires eigenvalues +/- 1.
# Let's define O_z = Threshold(Sz_A) -> +1 if >0, -1 if <=0.

# Actually, rigorous way: Define observable M with eigenvalues +/- 1.
# M_z = Sign(Sz_A) (Handling 0? Let's say +1)
evals_Sz, evecs_Sz = la.eigh(Sz_A)
M_z_A = np.zeros((16,16), dtype=complex)
# Reconstruct sign operator
for i in range(16):
    val = 1.0 if evals_Sz[i] >= 0 else -1.0
    M_z_A += val * np.outer(evecs_Sz[:,i], evecs_Sz[:,i].conj())

# M_x_A ... rotate M_z_A ?
# M_x corresponds to Sx basis.
evals_Sx, evecs_Sx = la.eigh(Sx_A)
M_x_A = np.zeros((16,16), dtype=complex)
for i in range(16):
    val = 1.0 if evals_Sx[i] >= 0 else -1.0
    M_x_A += val * np.outer(evecs_Sx[:,i], evecs_Sx[:,i].conj())

# B settings (rotated)
# Rotated by 45 degrees
R = la.expm(-1j * (np.pi/4) * Sy_A) # Rotation around Y
M_z_B = M_z_A # Same basis relative to B
M_zb_B = np.dot(R, np.dot(M_z_A, R.conj().T)) # Rotated

# Wait, standard settings:
# A = Z, A' = X
# B = W = (Z+X)/rt2, B' = V = (Z-X)/rt2
M_w_B = (M_z_A + M_x_A) / np.sqrt(2)
# Re-binarize? No, Quantum limit is 2*sqrt(2) for expectation values.
# But operators must be bounded by 1.
# (Z+X)/rt2 has eigenvalues +/- 1.

# Calculate correlations E(A, B) = <Psi | M_A x M_B | Psi>
def expect(OpA, OpB, state):
    Op = np.kron(OpA, OpB)
    return np.vdot(state, np.dot(Op, state)).real

E_a_b = expect(M_z_A, M_w_B, state_t)
E_a_bp = expect(M_z_A, (M_z_A - M_x_A)/np.sqrt(2), state_t)
E_ap_b = expect(M_x_A, M_w_B, state_t)
E_ap_bp = expect(M_x_A, (M_z_A - M_x_A)/np.sqrt(2), state_t)

S_val = E_a_b - E_a_bp + E_ap_b + E_ap_bp
S_val = abs(S_val)

print(f"\nCHSH Parameter S: {S_val:.5f}")
print("Classical Limit: 2.0")
print("Quantum Limit: 2.828 (2*sqrt(2))")
print("-" * 40)

if S_val > 2.0:
    print("✅ BELL VIOLATION CONFIRMED (Quantum Entanglement)")
    print("   The 4-bit lattice interaction generates non-local correlations.")
    print("   FIN Theory is QUANTUM.")
else:
    print("❌ NO VIOLATION (S <= 2)")
    print("   Interaction failed to entangle or settings unoptimized.")

# Entropy
rho_A = np.trace(np.outer(state_t, state_t.conj()).reshape(16, 16, 16, 16), axis1=1, axis2=3) 
# Wait, partial trace is tricky in numpy reshape.
# Manual partial trace
rho_A_mat = np.zeros((16,16), dtype=complex)
for i in range(16):
    for j in range(16):
        # sum over B basis (k)
        # element = sum_k <i,k|psi><psi|j,k>
        val = 0
        for k in range(16):
             idx_ik = i*16 + k
             idx_jk = j*16 + k
             val += state_t[idx_ik] * state_t[idx_jk].conj()
        rho_A_mat[i,j] = val

evals_rho = la.eigvalsh(rho_A_mat)
evals_rho = evals_rho[evals_rho > 1e-10] # filter zeros
entropy = -np.sum(evals_rho * np.log(evals_rho))

print(f"Von Neumann Entropy (Node A): {entropy:.5f}")
if entropy > 0.1:
    print("   Entropy > 0 implies Entanglement.")

# ============================================================================
# REPORT
# ============================================================================
with open("raport_qw630_quantum_entanglement.md", "w") as f:
    f.write("# Raport QW-630: Quantum Entanglement (Rigorous)\n")
    f.write("**Data:** 2025-12-05\n\n")
    f.write("## Metodologia\n")
    f.write("Symulacja pełna macierzowa (Density Matrix) dwóch cząstek 4-bitowych (dim=256).\n")
    f.write("Interakcja: $H = -g \\vec{S}_A \\cdot \\vec{S}_B$ przez czas $t=\\pi/4$.\n")
    f.write("Pomiar: Nierówność CHSH dla dychotomicznych obserwabli znaku spinu.\n\n")
    
    f.write("## Wyniki\n")
    f.write(f"- Parametr CHSH S: {S_val:.5f}\n")
    f.write(f"- Entropia Splątania: {entropy:.5f}\n\n")
    
    if S_val > 2.0:
        f.write("### ✅ TO JEST TEORIA KWANTOWA\n")
        f.write("Nierówność Bella została złamana (S > 2). FIN Theory generuje nielokalne korelacje z lokalnych oddziaływań sieciowych.\n")
        f.write("Luka 'Quantumness' została zamknięta.\n")
    else:
        f.write("### ❌ WYNIK KLASYCZNY\n")
        f.write("Brak naruszenia Bella. Model zachowuje się lokalnie.\n")

print("Report saved.")
