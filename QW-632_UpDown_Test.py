#!/usr/bin/env python3
# QW-632: UP-DOWN STATE ENTANGLEMENT
# Purpose: Test Entanglement generation starting from |Up>_A |Down>_B.
#          This state is NOT an eigenstate of Heisenberg H (which flip-flops).
#          This MUST generate entanglement.
# Date: 2025-12-05

import numpy as np
import scipy.linalg as la

print("="*80)
print("QW-632: UP-DOWN ENTANGLEMENT TEST")
print("="*80)

# Re-define Operators
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

Sz_A = sum(tensor_4(Z, i) for i in range(4)) / 2
Sx_A = sum(tensor_4(X, i) for i in range(4)) / 2
Sy_A = sum(tensor_4(Y, i) for i in range(4)) / 2

Sx_A_full = np.kron(Sx_A, np.eye(16))
Sx_B_full = np.kron(np.eye(16), Sx_A)
Sy_A_full = np.kron(Sy_A, np.eye(16))
Sy_B_full = np.kron(np.eye(16), Sy_A)
Sz_A_full = np.kron(Sz_A, np.eye(16))
Sz_B_full = np.kron(np.eye(16), Sz_A)

# Interaction H = - (SxSx + SySy + SzSz)
# Note: Sy is complex.
H_int = -(
    np.dot(Sx_A_full, Sx_B_full) + 
    np.dot(Sy_A_full, Sy_B_full) + 
    np.dot(Sz_A_full, Sz_B_full)
).real # H must be Hermitian/Real? No, Hermitian. complex?
# Sy*Sy involves i*i = -1. Should be real?
# Let's keep it complex but ensure Hermiticity.
H_int = -(
    np.dot(Sx_A_full, Sx_B_full) + 
    np.dot(Sy_A_full, Sy_B_full) + 
    np.dot(Sz_A_full, Sz_B_full)
)

print(f"H Hermitian? {np.allclose(H_int, H_int.conj().T)}")

# Initial State: |Up>_A |Down>_B
# Node A: |0000> (Spin +2)
psi_0 = np.array([1, 0], dtype=complex)
psi_A = psi_0
for _ in range(3): psi_A = np.kron(psi_A, psi_0)

# Node B: |1111> (Spin -2)
psi_1 = np.array([0, 1], dtype=complex)
psi_B = psi_1
for _ in range(3): psi_B = np.kron(psi_B, psi_1)

state_0 = np.kron(psi_A, psi_B)
print("Initial state |Up>|Down> prepared.")

# Evolve
times = np.linspace(0, 2, 10)
print("\nScanning time evolution...")

def get_entropy(state):
    rho_A = np.zeros((16,16), dtype=complex)
    for i in range(16):
        for j in range(16):
            val = 0
            for k in range(16): 
                 idx_ik = i*16 + k
                 idx_jk = j*16 + k
                 val += state[idx_ik] * state[idx_jk].conj()
            rho_A[i,j] = val
    evals = la.eigvalsh(rho_A)
    evals = evals[evals > 1e-10]
    return -np.sum(evals * np.log(evals))

for t in times:
    U = la.expm(-1j * H_int * t)
    st = np.dot(U, state_0)
    ent = get_entropy(st)
    print(f"t={t:.2f} S={ent:.4f}")

    if ent > 0.5:
        # Calculate CHSH for entangled state
        # ... logic as before ...
        pass

# ============================================================================
# REPORT
# ============================================================================
with open("raport_qw632_interaction.md", "w") as f:
    f.write("# Raport QW-632: Up-Down Interaction\n")
    f.write("Test state: |Up>|Down> (Sz eigenvalues +2, -2).\n")
    f.write("Evolution shows entropy generation?\n")

print("Report saved.")
