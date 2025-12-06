#!/usr/bin/env python3
"""
QW-686: DYNAMIC OBSERVER BELL TEST
===================================
Purpose: Test if perceived quantumness fluctuates over time for a DYNAMIC observer.

Hypothesis:
  "The Observer is a PROCESS, not a static entity."
  Entanglement/Quantumness (S) should fluctuate as the observer process evolves.

Paradigm:
  Nadsoliton = Fundamental Process
  Observer = Emergent Dynamic Sub-process (Octaves 0-3)
  Observed = Rest of System (Octaves 4-7)

Method:
  1. Initialize system in Ground State (or slight excitation)
  2. Evolve full state |ψ(t)⟩ = e^{-iHt} |ψ(0)⟩
  3. At each step t:
     - Trace out observer: ρ_obs(t) = Tr_observer(|ψ(t)⟩⟨ψ(t)|)
     - Measure S(t) from internal perspective using ρ_obs(t)
     
Output: RAPORT_QW686_DYNAMIC_OBSERVER.md
"""

import numpy as np
from scipy.linalg import eigh, expm
import datetime

print("="*80)
print("QW-686: DYNAMIC OBSERVER BELL TEST")
print("="*80)
print("Testing: How does perceived quantumness S(t) evolve for a DYNAMIC observer?")
print()

# --- Parameters ---
N_TOTAL = 8
N_OBSERVER = 4
N_OBSERVED = 4
ALPHA_GEO = 4 * np.log(2)
BETA_TORS = 0.01
OMEGA = np.pi / 4
PHI = np.pi / 6
STEPS = 50
DT = 0.1

# Pauli Matrices
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
sigma_0 = np.eye(2, dtype=complex)

REPORT_FILE = "RAPORT_QW686_DYNAMIC_OBSERVER.md"

# --- K(d) KERNEL ---
def K_coupling(d):
    return ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + BETA_TORS * d)

# --- BUILD OPERATORS ---
def build_spin_operator(op, site, n_sites):
    result = 1
    for i in range(n_sites):
        if i == site:
            result = np.kron(result, op)
        else:
            result = np.kron(result, sigma_0)
    return result

print("1. Building Hamiltonian...")
dim_total = 2**N_TOTAL
dim_obs = 2**N_OBSERVER
dim_meas = 2**N_OBSERVED

Sx = [build_spin_operator(0.5*sigma_x, i, N_TOTAL) for i in range(N_TOTAL)]
Sy = [build_spin_operator(0.5*sigma_y, i, N_TOTAL) for i in range(N_TOTAL)]
Sz = [build_spin_operator(0.5*sigma_z, i, N_TOTAL) for i in range(N_TOTAL)]

H = np.zeros((dim_total, dim_total), dtype=complex)
for i in range(N_TOTAL):
    for j in range(i+1, N_TOTAL):
        d = abs(i - j)
        K = K_coupling(d)
        H += K * (Sx[i] @ Sx[j] + Sy[i] @ Sy[j] + Sz[i] @ Sz[j])

# --- INITIAL STATE ---
print("\n2. Preparing initial state...")
# Start with ground state + slight perturbation (to induce dynamics)
eigenvalues, eigenvectors = eigh(H)
psi_ground = eigenvectors[:, 0]
psi_excited = eigenvectors[:, 1]

# Mix states: |ψ(0)⟩ = 0.8|GS⟩ + 0.6|ES⟩
psi_t = 0.8 * psi_ground + 0.6 * psi_excited
psi_t /= np.linalg.norm(psi_t)

E_0 = np.real(psi_t.conj().T @ H @ psi_t)
print(f"   E_avg = {E_0:.4f}")

# --- INTERNAL MEASUREMENT FUNCTIONS ---
Sx_sub = [build_spin_operator(0.5*sigma_x, i, N_OBSERVED) for i in range(N_OBSERVED)] # Actually careful here
# Need proper subsystem operators for trace calculation
# Actually, for trace, we just reshape and contract

def get_reduced_rho(psi, dim_A, dim_B):
    """Trace out system A (first dim_A states), keep system B (dim_B states)"""
    psi_reshaped = psi.reshape(dim_A, dim_B)
    # rho_B = sum_i <i_A|psi><psi|i_A>
    # psi_reshaped[i, j] = coeff for |i>_A |j>_B
    # We want to sum over i (trace out A)
    # rho_B[j, k] = sum_i psi[i, j] * conj(psi[i, k])
    # This corresponds to: psi.T @ psi.conj()  <-- Wait.
    # (dim_B, dim_A) @ (dim_A, dim_B) -> (dim_B, dim_B)
    rho_B = psi_reshaped.T @ psi_reshaped.conj()
    return rho_B

# Subsystem operators (for observed system only)
def build_sub_op(op, site, n_sites):
    result = 1
    for i in range(n_sites):
        if i == site:
            result = np.kron(result, op)
        else:
            result = np.kron(result, sigma_0)
    return result

Sx_meas = [build_sub_op(0.5*sigma_x, i, N_OBSERVED) for i in range(N_OBSERVED)]
Sz_meas = [build_sub_op(0.5*sigma_z, i, N_OBSERVED) for i in range(N_OBSERVED)]

def spin_angle_meas(theta, site):
    return np.sin(theta)*Sx_meas[site] + np.cos(theta)*Sz_meas[site]

# CHSH angles
theta_a = 0
theta_a_prime = np.pi / 2
theta_b = np.pi / 4
theta_b_prime = 3 * np.pi / 4

# Measure between spin 0 and 2 of OBSERVED system (which are spins 4 and 6 of total)
spin_A = 0
spin_B = 2

A = 2 * spin_angle_meas(theta_a, spin_A)
Ap = 2 * spin_angle_meas(theta_a_prime, spin_A)
B = 2 * spin_angle_meas(theta_b, spin_B)
Bp = 2 * spin_angle_meas(theta_b_prime, spin_B)

def measure_S_internal(rho):
    def trace_val(Op):
        return np.real(np.trace(rho @ Op))
    
    E1 = trace_val(A @ B)
    E2 = trace_val(A @ Bp)
    E3 = trace_val(Ap @ B)
    E4 = trace_val(Ap @ Bp)
    return np.abs(E1 - E2 + E3 + E4)

def get_entanglement_entropy(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-12]
    return -np.sum(evals * np.log(evals))

# --- EVOLUTION LOOP ---
print("\n3. Running time evolution...")
print(f"   Steps: {STEPS}, dt = {DT}")

# Time evolution operator U = e^{-iHdt}
U = expm(-1j * H * DT)

history_S = []
history_SvN = []
history_time = []

for t in range(STEPS):
    # 1. Get reduced density matrix
    # Partition: A=Observer (trace out), B=Observed (keep)
    rho_observed = get_reduced_rho(psi_t, dim_obs, dim_meas)
    
    # 2. Measure S and S_vN
    S_val = measure_S_internal(rho_observed)
    S_vN = get_entanglement_entropy(rho_observed)
    
    history_S.append(S_val)
    history_SvN.append(S_vN)
    history_time.append(t * DT)
    
    if t % 10 == 0:
        print(f"   t={t*DT:.1f}: S={S_val:.4f}, S_vN={S_vN:.4f}")
    
    # 3. Evolve
    psi_t = U @ psi_t
    psi_t /= np.linalg.norm(psi_t)

# --- ANALYSIS ---
print("\n4. Analysis...")
S_min, S_max = min(history_S), max(history_S)
S_mean = np.mean(history_S)
S_std = np.std(history_S)

print(f"   S(t) Range: [{S_min:.4f}, {S_max:.4f}]")
print(f"   Mean S: {S_mean:.4f} ± {S_std:.4f}")

violation_fraction = sum(1 for s in history_S if s > 2.0) / len(history_S)
print(f"   Fraction of time S > 2.0: {violation_fraction*100:.1f}%")

# --- WRITE REPORT ---
print("\n5. Writing report...")

with open(REPORT_FILE, "w") as f:
    f.write("# RAPORT: QW-686 DYNAMIC OBSERVER BELL TEST\n")
    f.write(f"**Data:** {datetime.datetime.now()}\n")
    f.write("**Hipoteza:** Kwantowość S(t) fluktuuje w czasie dla emergentnego obserwatora (procesu)\n\n")
    
    f.write("## 1. METODOLOGIA\n\n")
    f.write("- **System:** N=8 (Observer=4, Observed=4)\n")
    f.write("- **Dynamics:** Perturbed Ground State evolution\n")
    f.write("- **Perspective:** Internal (Reduced Density Matrix)\n\n")
    
    f.write("## 2. WYNIKI DYNAMICZNE\n\n")
    f.write(f"- **S_max:** {S_max:.4f}\n")
    f.write(f"- **S_min:** {S_min:.4f}\n")
    f.write(f"- **Mean:** {S_mean:.4f}\n")
    f.write(f"- **Violation Time:** {violation_fraction*100:.1f}%\n\n")
    
    f.write("### Przebieg (próbka):\n")
    f.write("| t | S(t) | S_vN(t) |\n")
    f.write("|---|------|---------|\n")
    for i in range(0, len(history_S), 5):
        f.write(f"| {history_time[i]:.1f} | {history_S[i]:.4f} | {history_SvN[i]:.4f} |\n")
    
    f.write("\n## 3. INTERPRETACJA\n\n")
    if S_std > 0.05:
        f.write("### ✅ POTWIERDZENIE: KWANTOWOŚĆ JEST DYNAMICZNA\n\n")
        f.write("Wartość S(t) fluktuuje znacząco w czasie.\n")
        f.write("To oznacza, że 'kwantowość' nie jest stałą cechą, ale ZMIENNĄ fazą procesu.\n\n")
        if violation_fraction > 0:
            f.write(f"**Obserwator widzi łamanie Bella przez {violation_fraction*100:.1f}% czasu!**\n")
    else:
        f.write("### ❌ BRAK DYNAMIKI\n\n")
        f.write("S(t) jest stałe. Perturbacja nie wpłynęła na korelacje.\n")

print(f"   Report saved to: {REPORT_FILE}")

print("\n" + "="*80)
print("QW-686 COMPLETE")
print("="*80)
print(f"\nWYNIKI KOŃCOWE:")
print(f"  S_max = {S_max:.4f}")
print(f"  Violation fraction: {violation_fraction*100:.1f}%")
if S_std > 0.05:
    print("  ✅ DYNAMIC FLUCTUATIONS CONFIRMED")
else:
    print("  ❌ STATIC CORRELATIONS")
