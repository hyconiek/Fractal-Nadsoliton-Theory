
import numpy as np
from scipy.linalg import eigh
import datetime

print("="*80)
print("QW-692: RIGOROUS SCIENTIFIC VERIFICATION")
print("="*80)

# Pauli Matrices
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
sigma_0 = np.eye(2, dtype=complex)

def build_op(op, site, n_sites):
    result = 1
    for i in range(n_sites):
        if i == site:
            result = np.kron(result, op)
        else:
            result = np.kron(result, sigma_0)
    return result

def solve_system(n_per_chain, j_vertical, j_pair):
    n_total = 2 * n_per_chain
    dim_total = 2**n_total
    
    # Operators
    Sx = [build_op(sigma_x, i, n_total) for i in range(n_total)]
    Sz = [build_op(sigma_z, i, n_total) for i in range(n_total)]
    
    # Hamiltonian
    # Chain A: 0 to n-1
    # Chain B: n to 2n-1
    H = np.zeros((dim_total, dim_total), dtype=complex)
    
    # Vertical Coupling (Intra-particle)
    for i in range(n_per_chain - 1):
        # Chain A
        H += -j_vertical * Sz[i] @ Sz[i+1] # Ising vertical
        H += -1.0 * Sx[i] # Transverse field
        
        # Chain B
        j = i + n_per_chain
        H += -j_vertical * Sz[j] @ Sz[j+1]
        H += -1.0 * Sx[j]
        
    # Last sites transverse
    H += -1.0 * Sx[n_per_chain-1]
    H += -1.0 * Sx[n_total-1]
        
    # Entanglement (Layer 0 to Layer 0)
    # Target singlet-like state: Minimize energy with anti-alignment or specific correlation
    # To get S ~ 2.82, we need a Bell state like (|00> + |11>)/sqrt(2) or similar.
    # Standard Heisenberg interaction or Ising?
    # Let's use the one from the original file: +J * (XX + ZZ)
    H += j_pair * (Sx[0] @ Sx[n_per_chain] + Sz[0] @ Sz[n_per_chain])
    
    # Solve
    evals, evecs = eigh(H)
    psi = evecs[:, 0]
    
    return psi, Sx, Sz

def measure_S(psi, Sx, Sz, n_per_chain, mode="lab"):
    if mode == "lab":
        # Measure Layer 0
        idxA = 0
        idxB = n_per_chain
        opsA = [(Sx[idxA], Sz[idxA], 1.0)]
        opsB = [(Sx[idxB], Sz[idxB], 1.0)]
    else:
        # Natural: Average over all layers
        opsA = []
        opsB = []
        weight = 1.0 / n_per_chain
        for l in range(n_per_chain):
            opsA.append((Sx[l], Sz[l], weight))
            opsB.append((Sx[n_per_chain+l], Sz[n_per_chain+l], weight))
            
    # Standard Bell Angles
    angles = [0, np.pi/2, np.pi/4, 3*np.pi/4]
    theta_a, theta_ap, theta_b, theta_bp = angles
    
    def get_corr(thA, thB):
        # Construct Operator A
        OpA_mat = 0
        for sx, sz, w in opsA:
            OpA_mat += w * (np.sin(thA)*sx + np.cos(thA)*sz)
            
        # Construct Operator B
        OpB_mat = 0
        for sx, sz, w in opsB:
            OpB_mat += w * (np.sin(thB)*sx + np.cos(thB)*sz)
            
        Expectation = psi.conj().T @ (OpA_mat @ OpB_mat) @ psi
        return np.real(Expectation)

    E1 = get_corr(theta_a, theta_b)
    E2 = get_corr(theta_a, theta_bp)
    E3 = get_corr(theta_ap, theta_b)
    E4 = get_corr(theta_ap, theta_bp)
    
    return np.abs(E1 - E2 + E3 + E4)

# --- RUN SWEEP ---
results = []
# Vary N
n_values = [2, 3, 4, 5, 6]
j_pair_values = [1.0, 2.0, 5.0, 10.0]

print(f"{'N':<5} | {'J_pair':<8} | {'S_lab':<8} | {'S_nat':<8} | {'Gain':<8}")
print("-" * 50)

best_s_lab = 0
best_s_nat = 100

for n in n_values:
    # Scale vertical coupling? 
    # If layers are "fractal", maybe coupling is constant.
    j_vert = 1.0 
    
    for j_p in j_pair_values:
        try:
            psi, Sx, Sz = solve_system(n, j_vert, j_p)
            s_lab = measure_S(psi, Sx, Sz, n, "lab")
            s_nat = measure_S(psi, Sx, Sz, n, "natural")
            
            print(f"{n:<5} | {j_p:<8.1f} | {s_lab:<8.4f} | {s_nat:<8.4f} | {s_lab-s_nat:<8.4f}")
            results.append((n, j_p, s_lab, s_nat))
            
            if s_lab > best_s_lab: best_s_lab = s_lab
            if s_nat < best_s_nat: best_s_nat = s_nat
            
        except Exception as e:
            print(f"Error for N={n}: {e}")

print("-" * 50)
print(f"Best Lab Quantumness: {best_s_lab:.4f}")
print(f"Lowest Natural Class: {best_s_nat:.4f}")
