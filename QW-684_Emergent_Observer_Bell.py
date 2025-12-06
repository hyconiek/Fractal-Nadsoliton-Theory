#!/usr/bin/env python3
"""
QW-684: EMERGENT OBSERVER BELL TEST
====================================
Purpose: Test if Bell violation appears from INTERNAL observer perspective

User Hypothesis:
  "Obserwator emergentny może mieć znaczenie"
  
  If observer is EMERGENT from nadsoliton:
  1. We cannot observe from "outside" 
  2. Bell test is self-referential (nadsoliton measuring itself)
  3. Quantum correlations may be perspectival

Method:
  1. Partition system: Observer (octaves 0-3) + Observed (octaves 4-7)
  2. Observer has access only to REDUCED density matrix
  3. Compare: External Bell (full state) vs Internal Bell (observer's view)
  
Key Insight:
  External: |ψ⟩ → S_external
  Internal: Tr_obs(|ψ⟩⟨ψ|) = ρ_measured → S_internal
  
  Hypothesis: S_internal > S_external (observer sees more "quantumness")

Output: RAPORT_QW684_EMERGENT_OBSERVER.md
"""

import numpy as np
from scipy.linalg import eigh, expm
import datetime

print("="*80)
print("QW-684: EMERGENT OBSERVER BELL TEST")
print("="*80)
print("Testing if Bell violation appears from INTERNAL observer perspective")
print()

# --- Parameters ---
N_TOTAL = 8
N_OBSERVER = 4  # Observer is first 4 octaves
N_OBSERVED = 4  # Observed is last 4 octaves
ALPHA_GEO = 4 * np.log(2)
BETA_TORS = 0.01
OMEGA = np.pi / 4
PHI = np.pi / 6

# Pauli Matrices
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
sigma_0 = np.eye(2, dtype=complex)

REPORT_FILE = "RAPORT_QW684_EMERGENT_OBSERVER.md"

# --- K(d) KERNEL ---
def K_coupling(d):
    return ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + BETA_TORS * d)

# --- BUILD OPERATORS ---
dim_total = 2**N_TOTAL
dim_obs = 2**N_OBSERVER
dim_meas = 2**N_OBSERVED

def build_spin_operator(op, site, n_sites):
    result = 1
    for i in range(n_sites):
        if i == site:
            result = np.kron(result, op)
        else:
            result = np.kron(result, sigma_0)
    return result

print("1. Building system...")
print(f"   Total: {N_TOTAL} spins (dim = {dim_total})")
print(f"   Observer: spins 0-{N_OBSERVER-1} (dim = {dim_obs})")
print(f"   Observed: spins {N_OBSERVER}-{N_TOTAL-1} (dim = {dim_meas})")

# Build full operators
Sx = [build_spin_operator(0.5*sigma_x, i, N_TOTAL) for i in range(N_TOTAL)]
Sy = [build_spin_operator(0.5*sigma_y, i, N_TOTAL) for i in range(N_TOTAL)]
Sz = [build_spin_operator(0.5*sigma_z, i, N_TOTAL) for i in range(N_TOTAL)]

# Build Hamiltonian
H = np.zeros((dim_total, dim_total), dtype=complex)
for i in range(N_TOTAL):
    for j in range(i+1, N_TOTAL):
        d = abs(i - j)
        K = K_coupling(d)
        H += K * (Sx[i] @ Sx[j] + Sy[i] @ Sy[j] + Sz[i] @ Sz[j])

# --- FIND GROUND STATE ---
print("\n2. Finding ground state...")
eigenvalues, eigenvectors = eigh(H)
psi_full = eigenvectors[:, 0]
E_ground = eigenvalues[0]
print(f"   E_0 = {E_ground:.4f}")

# --- EXTERNAL BELL TEST (FULL STATE) ---
print("\n3. EXTERNAL Bell test (from 'outside')...")

def measure_correlation_full(psi, obs_A, obs_B):
    AB = obs_A @ obs_B
    return np.real(psi.conj().T @ AB @ psi)

theta_a = 0
theta_a_prime = np.pi / 2
theta_b = np.pi / 4
theta_b_prime = 3 * np.pi / 4

def spin_at_angle(theta, site):
    return np.sin(theta)*Sx[site] + np.cos(theta)*Sz[site]

# External: measure spins 0 (in observer) and 4 (in observed)
A_ext = 2 * spin_at_angle(theta_a, 0)
Ap_ext = 2 * spin_at_angle(theta_a_prime, 0)
B_ext = 2 * spin_at_angle(theta_b, N_OBSERVER)
Bp_ext = 2 * spin_at_angle(theta_b_prime, N_OBSERVER)

E1_ext = measure_correlation_full(psi_full, A_ext, B_ext)
E2_ext = measure_correlation_full(psi_full, A_ext, Bp_ext)
E3_ext = measure_correlation_full(psi_full, Ap_ext, B_ext)
E4_ext = measure_correlation_full(psi_full, Ap_ext, Bp_ext)

S_external = np.abs(E1_ext - E2_ext + E3_ext + E4_ext)
print(f"   S_external = {S_external:.4f}")

# --- INTERNAL BELL TEST (FROM OBSERVER'S PERSPECTIVE) ---
print("\n4. INTERNAL Bell test (from observer's perspective)...")

# Observer's reduced density matrix (trace out observer, keep observed)
# ρ_observed = Tr_observer(|ψ⟩⟨ψ|)
psi_matrix = psi_full.reshape(dim_obs, dim_meas)
rho_observed = psi_matrix.conj().T @ psi_matrix  # dim_meas × dim_meas

# Build operators on observed subsystem only
def build_spin_operator_subsystem(op, site, n_sites):
    result = 1
    for i in range(n_sites):
        if i == site:
            result = np.kron(result, op)
        else:
            result = np.kron(result, sigma_0)
    return result

Sx_obs = [build_spin_operator_subsystem(0.5*sigma_x, i, N_OBSERVED) for i in range(N_OBSERVED)]
Sy_obs = [build_spin_operator_subsystem(0.5*sigma_y, i, N_OBSERVED) for i in range(N_OBSERVED)]
Sz_obs = [build_spin_operator_subsystem(0.5*sigma_z, i, N_OBSERVED) for i in range(N_OBSERVED)]

def spin_at_angle_obs(theta, site):
    return np.sin(theta)*Sx_obs[site] + np.cos(theta)*Sz_obs[site]

# Internal: observer measures TWO spins in observed subsystem
# (spins 0 and 2 in observed = spins 4 and 6 in total)
A_int = 2 * spin_at_angle_obs(theta_a, 0)
Ap_int = 2 * spin_at_angle_obs(theta_a_prime, 0)
B_int = 2 * spin_at_angle_obs(theta_b, 2)
Bp_int = 2 * spin_at_angle_obs(theta_b_prime, 2)

def measure_correlation_internal(rho, obs_A, obs_B):
    """Observer's view: Tr(ρ A⊗B)"""
    AB = obs_A @ obs_B
    return np.real(np.trace(rho @ AB))

E1_int = measure_correlation_internal(rho_observed, A_int, B_int)
E2_int = measure_correlation_internal(rho_observed, A_int, Bp_int)
E3_int = measure_correlation_internal(rho_observed, Ap_int, B_int)
E4_int = measure_correlation_internal(rho_observed, Ap_int, Bp_int)

S_internal = np.abs(E1_int - E2_int + E3_int + E4_int)
print(f"   S_internal = {S_internal:.4f}")

# --- COMPARISON ---
print("\n5. COMPARISON...")
print(f"   S_external (outside view): {S_external:.4f}")
print(f"   S_internal (observer view): {S_internal:.4f}")
print(f"   Difference: {S_internal - S_external:.4f}")

if S_internal > S_external:
    comparison = "✅ S_internal > S_external: Observer sees MORE quantumness!"
elif S_internal > 2.0:
    comparison = "✅ S_internal > 2: Observer sees QUANTUM correlations!"
else:
    comparison = "❌ No enhancement from internal perspective"
    
print(f"\n   {comparison}")

# --- ENTANGLEMENT ANALYSIS ---
print("\n6. Entanglement analysis...")

# Observer's entanglement with observed (from perspective of observer)
eigenvalues_rho = np.linalg.eigvalsh(rho_observed)
eigenvalues_rho = eigenvalues_rho[eigenvalues_rho > 1e-12]
S_entanglement = -np.sum(eigenvalues_rho * np.log(eigenvalues_rho))
print(f"   S_vN (observer-observed entanglement): {S_entanglement:.4f}")

# Purity of observed state (from observer's perspective)
purity = np.real(np.trace(rho_observed @ rho_observed))
print(f"   Purity of observed state: {purity:.4f}")
print(f"   (Purity = 1 means pure, < 1 means mixed = observer entangled)")

# --- WRITE REPORT ---
print("\n7. Writing report...")

with open(REPORT_FILE, "w") as f:
    f.write("# RAPORT: QW-684 EMERGENT OBSERVER BELL TEST\n")
    f.write(f"**Data:** {datetime.datetime.now()}\n")
    f.write("**Hipoteza:** Emergentny obserwator widzi inną kwantowość niż widok zewnętrzny\n\n")
    
    f.write("## 1. KONCEPCJA\n\n")
    f.write("- **Observer:** Oktawy 0-3 (emergentna część nadsolitona)\n")
    f.write("- **Observed:** Oktawy 4-7 (reszta systemu)\n")
    f.write("- **External view:** Pełna funkcja falowa |ψ⟩\n")
    f.write("- **Internal view:** Reduced density matrix Tr_obs(|ψ⟩⟨ψ|)\n\n")
    
    f.write("## 2. WYNIKI BELL\n\n")
    f.write("| Perspektywa | S (CHSH) | Violation? |\n")
    f.write("|-------------|----------|------------|\n")
    f.write(f"| External | {S_external:.4f} | {'✅' if S_external > 2 else '❌'} |\n")
    f.write(f"| Internal (observer) | {S_internal:.4f} | {'✅' if S_internal > 2 else '❌'} |\n\n")
    
    f.write("## 3. SPLĄTANIE\n\n")
    f.write(f"- **S_vN (obs↔meas):** {S_entanglement:.4f}\n")
    f.write(f"- **Purity:** {purity:.4f}\n\n")
    
    f.write("## 4. INTERPRETACJA\n\n")
    if S_internal > S_external:
        f.write("### ✅ EMERGENTNA PERSPEKTYWA WZMACNIA KORELACJE\n\n")
        f.write(f"Obserwator (część systemu) widzi S = {S_internal:.4f} > {S_external:.4f} (zewnętrzne).\n")
        f.write("To sugeruje, że kwantowość jest **perspektywalna** - zależy od tego, czy mierzymy z wewnątrz.\n")
    elif S_internal > 2.0:
        f.write("### ✅ TYLKO OBSERWATOR WIDZI ŁAMANIE BELLA!\n\n")
        f.write("Z zewnątrz: S < 2 (klasyczne)\n")
        f.write("Z wewnątrz: S > 2 (kwantowe)\n\n")
        f.write("**WNIOSEK:** Kwantowość wyłania się z perspektywy emergentnego obserwatora!\n")
    else:
        f.write("### ❌ BRAK WZMOCNIENIA Z PERSPEKTYWY WEWNĘTRZNEJ\n\n")
        f.write(f"S_external = {S_external:.4f}, S_internal = {S_internal:.4f}\n")
        f.write("Obie perspektywy dają S < 2.\n")
    
    f.write("\n## 5. PORÓWNANIE Z POPRZEDNIMI\n\n")
    f.write("| Badanie | S (CHSH) | Perspektywa |\n")
    f.write("|---------|----------|-------------|\n")
    f.write("| QW-680 | 0.17 | zewnętrzna |\n")
    f.write("| QW-682 | 0.66 | zewnętrzna |\n")
    f.write("| QW-683 | 0.43 | zewnętrzna |\n")
    f.write(f"| **QW-684** | **{S_internal:.2f}** | **wewnętrzna** |\n")

print(f"   Report saved to: {REPORT_FILE}")

print("\n" + "="*80)
print("QW-684 COMPLETE")
print("="*80)
print(f"\nKLUCZOWE WYNIKI:")
print(f"  S_external: {S_external:.4f}")
print(f"  S_internal: {S_internal:.4f}")
print(f"  Różnica: {S_internal - S_external:+.4f}")
print(f"\n  {comparison}")
