#!/usr/bin/env python3
# QW-631: ENTANGLEMENT DEBUG SCAN
# Purpose: Find conditions where Entanglement Entropy > 0.
#          Scan Initial States and Time steps.
#          Debug Operator Definitions.
# Date: 2025-12-05

import numpy as np
import scipy.linalg as la

print("="*80)
print("QW-631: ENTANGLEMENT DEBUG")
print("="*80)

# Re-define Operators (carefully)
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

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

# Full Space Operators (256x256)
Sx_A_full = np.kron(Sx_A, np.eye(16))
Sx_B_full = np.kron(np.eye(16), Sx_A)
Sy_A_full = np.kron(Sy_A, np.eye(16))
Sy_B_full = np.kron(np.eye(16), Sy_A)
Sz_A_full = np.kron(Sz_A, np.eye(16))
Sz_B_full = np.kron(np.eye(16), Sz_A)

# Interaction: Heisenberg
g = 1.0
H_int = -g * (
    np.dot(Sx_A_full, Sx_B_full) + 
    np.dot(Sy_A_full, Sy_B_full) + 
    np.dot(Sz_A_full, Sz_B_full)
)

print(f"H Norm: {np.linalg.norm(H_int)}")

# Scan Initial States
# State 1: |0000>|0000> (All Up)
# This is eigenstate of SzSz, but NOT SxSx?
# SxSx flips spins: |0000> -> |1000> etc.
# So it SHOULD evolve.

# State 2: Superposition |+>|+>...
plus = np.array([1, 1], dtype=complex)/np.sqrt(2)
psi_plus = plus
for _ in range(3): psi_plus = np.kron(psi_plus, plus)
state_plus = np.kron(psi_plus, psi_plus)

def get_entropy(state):
    # Reduced rho A
    # Manual partial trace
    rho_A = np.zeros((16,16), dtype=complex)
    for i in range(16):
        for j in range(16):
            val = 0
            for k in range(16): # Trace out B
                 idx_ik = i*16 + k
                 idx_jk = j*16 + k
                 val += state[idx_ik] * state[idx_jk].conj()
            rho_A[i,j] = val
    
    evals = la.eigvalsh(rho_A)
    evals = evals[evals > 1e-10]
    return -np.sum(evals * np.log(evals))

print("\nRunning Time Scan...")
times = np.linspace(0, np.pi, 20)
max_S = 0.0

psi_0 = np.array([1, 0], dtype=complex)
psi_node = psi_0
for _ in range(3): psi_node = np.kron(psi_node, psi_0)
state_vac = np.kron(psi_node, psi_node)

print(f"Checking Vacuum State Evolution...")
for t in times:
    U = la.expm(-1j * H_int * t)
    st = np.dot(U, state_vac) # Changed to simple dot via expon
    # Check if state changed
    overlap = abs(np.vdot(state_vac, st))
    ent = get_entropy(st)
    if ent > max_S: max_S = ent
    # print(f"t={t:.2f} Overlap={overlap:.4f} Entropy={ent:.4f}")

print(f"Max Entropy (Vacuum Start): {max_S:.5f}")

print(f"Checking Superposition State Evolution...")
max_S_plus = 0.0
for t in times:
    U = la.expm(-1j * H_int * t)
    st = np.dot(U, state_plus)
    ent = get_entropy(st)
    if ent > max_S_plus: max_S_plus = ent

print(f"Max Entropy (Plus Start): {max_S_plus:.5f}")

if max_S_plus > 0.1:
    print("✅ ENTANGLEMENT GENERATED (in superposition)")
else:
    print("❌ STILL NO ENTANGLEMENT")
    # Debug: Check commutator [H, State]
    comm = np.dot(H_int, state_plus)
    print(f"H|psi> norm: {np.linalg.norm(comm)}")

# ============================================================================
# REPORT
# ============================================================================
with open("raport_qw631_debug.md", "w") as f:
    f.write("# Raport QW-631: Entanglement Debug\n")
    f.write(f"Max Entropy (Vac): {max_S}\n")
    f.write(f"Max Entropy (Plus): {max_S_plus}\n")

print("Report saved.")
