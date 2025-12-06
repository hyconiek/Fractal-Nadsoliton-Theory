#!/usr/bin/env python3
"""
QW-687: OBSERVER HIERARCHY TEST
================================
Purpose: Test how observer SIZE affects perceived quantumness

Hypothesis: Different "scales" of emergent observers see different physics.
  - Small observer (1-2 octaves): Sees strong quantum effects
  - Large observer (6-7 octaves): Sees classical physics (measuring self)

Paradigm:
  "Nadsoliton is the ONLY fundamental. Everything else is emergent,
   INCLUDING the observer."

Method:
  1. Fix total system: 10 octaves
  2. Vary observer size: 1, 2, 3, 4, 5 octaves
  3. Observed = remaining octaves
  4. Measure S (CHSH) from observer's internal perspective
  
Expected: S decreases as observer size increases (classical limit)

Output: RAPORT_QW687_OBSERVER_HIERARCHY.md
"""

import numpy as np
from scipy.linalg import eigh
import datetime

print("="*80)
print("QW-687: OBSERVER HIERARCHY TEST")
print("="*80)
print("Testing: How does observer SIZE affect perceived quantumness?")
print()

# --- Parameters ---
N_TOTAL = 10  # Total system size
ALPHA_GEO = 4 * np.log(2)
BETA_TORS = 0.01
OMEGA = np.pi / 4
PHI = np.pi / 6

# Pauli Matrices
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
sigma_0 = np.eye(2, dtype=complex)

REPORT_FILE = "RAPORT_QW687_OBSERVER_HIERARCHY.md"

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

def compute_entanglement_entropy(psi, n_sites, partition):
    """Compute S_vN for bipartition at 'partition' """
    dim_A = 2**partition
    dim_B = 2**(n_sites - partition)
    psi_matrix = psi.reshape(dim_A, dim_B)
    rho_A = psi_matrix @ psi_matrix.conj().T
    eigenvalues = np.linalg.eigvalsh(rho_A)
    eigenvalues = eigenvalues[eigenvalues > 1e-12]
    return -np.sum(eigenvalues * np.log(eigenvalues)) if len(eigenvalues) > 0 else 0

# --- BUILD HAMILTONIAN ---
print("1. Building full Hamiltonian...")
dim_total = 2**N_TOTAL
print(f"   N = {N_TOTAL} octaves, dim = {dim_total}")

Sx = [build_spin_operator(0.5*sigma_x, i, N_TOTAL) for i in range(N_TOTAL)]
Sy = [build_spin_operator(0.5*sigma_y, i, N_TOTAL) for i in range(N_TOTAL)]
Sz = [build_spin_operator(0.5*sigma_z, i, N_TOTAL) for i in range(N_TOTAL)]

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

# --- OBSERVER HIERARCHY TEST ---
print("\n3. Testing observer hierarchy...")

# CHSH angles
theta_a = 0
theta_a_prime = np.pi / 2
theta_b = np.pi / 4
theta_b_prime = 3 * np.pi / 4

results = []
observer_sizes = [1, 2, 3, 4, 5]

for n_observer in observer_sizes:
    n_observed = N_TOTAL - n_observer
    
    if n_observed < 2:
        print(f"   Observer {n_observer}: Skipping (need at least 2 observed)")
        continue
    
    dim_obs = 2**n_observer
    dim_meas = 2**n_observed
    
    print(f"\n   Observer size: {n_observer} octaves, Observed: {n_observed} octaves")
    
    # Get reduced density matrix (trace out observer, keep observed)
    psi_matrix = psi_full.reshape(dim_obs, dim_meas)
    rho_observed = psi_matrix.conj().T @ psi_matrix
    
    # Build operators on observed subsystem
    def build_spin_op_sub(op, site, n_sites):
        result = 1
        for i in range(n_sites):
            if i == site:
                result = np.kron(result, op)
            else:
                result = np.kron(result, sigma_0)
        return result
    
    Sx_obs = [build_spin_op_sub(0.5*sigma_x, i, n_observed) for i in range(n_observed)]
    Sy_obs = [build_spin_op_sub(0.5*sigma_y, i, n_observed) for i in range(n_observed)]
    Sz_obs = [build_spin_op_sub(0.5*sigma_z, i, n_observed) for i in range(n_observed)]
    
    def spin_at_angle_obs(theta, site):
        return np.sin(theta)*Sx_obs[site] + np.cos(theta)*Sz_obs[site]
    
    def measure_correlation_internal(rho, obs_A, obs_B):
        AB = obs_A @ obs_B
        return np.real(np.trace(rho @ AB))
    
    # Measure CHSH between first and last spin in observed subsystem
    spin_A = 0
    spin_B = n_observed - 1
    
    A = 2 * spin_at_angle_obs(theta_a, spin_A)
    Ap = 2 * spin_at_angle_obs(theta_a_prime, spin_A)
    B = 2 * spin_at_angle_obs(theta_b, spin_B)
    Bp = 2 * spin_at_angle_obs(theta_b_prime, spin_B)
    
    E1 = measure_correlation_internal(rho_observed, A, B)
    E2 = measure_correlation_internal(rho_observed, A, Bp)
    E3 = measure_correlation_internal(rho_observed, Ap, B)
    E4 = measure_correlation_internal(rho_observed, Ap, Bp)
    
    S_internal = np.abs(E1 - E2 + E3 + E4)
    
    # Compute entanglement within observed subsystem
    if n_observed >= 2:
        # Entanglement between first and second half of observed
        S_vN = compute_entanglement_entropy(psi_full, N_TOTAL, n_observer + n_observed//2)
    else:
        S_vN = 0
    
    # Purity of observed state
    purity = np.real(np.trace(rho_observed @ rho_observed))
    
    results.append({
        'n_observer': n_observer,
        'n_observed': n_observed,
        'S_chsh': S_internal,
        'S_vN': S_vN,
        'purity': purity
    })
    
    print(f"   S (CHSH) = {S_internal:.4f}, S_vN = {S_vN:.4f}, Purity = {purity:.4f}")

# --- ANALYSIS ---
print("\n4. Analysis...")

# Check if S decreases with observer size
observer_sizes_done = [r['n_observer'] for r in results]
S_values = [r['S_chsh'] for r in results]

if len(S_values) >= 2:
    corr = np.corrcoef(observer_sizes_done, S_values)[0, 1]
    print(f"   Correlation (observer_size vs S): r = {corr:.4f}")
    
    if corr < -0.5:
        trend = "✅ S DECREASES with observer size (classical limit emerges)"
    elif corr > 0.5:
        trend = "⚠️ S INCREASES with observer size (unexpected!)"
    else:
        trend = "❌ No clear trend"
else:
    corr = 0
    trend = "Not enough data"

print(f"   {trend}")

# --- WRITE REPORT ---
print("\n5. Writing report...")

with open(REPORT_FILE, "w") as f:
    f.write("# RAPORT: QW-687 OBSERVER HIERARCHY TEST\n")
    f.write(f"**Data:** {datetime.datetime.now()}\n")
    f.write("**Hipoteza:** Rozmiar obserwatora wpływa na postrzeganą kwantowość\n\n")
    
    f.write("## 1. PARADYGMAT\n\n")
    f.write("> Nadsoliton jest JEDYNYM fundamentem.\n")
    f.write("> WSZYSTKO inne jest emergentne, w tym obserwator.\n\n")
    
    f.write("## 2. WYNIKI\n\n")
    f.write("| Observer | Observed | S (CHSH) | S_vN | Purity |\n")
    f.write("|----------|----------|----------|------|--------|\n")
    for r in results:
        f.write(f"| {r['n_observer']} | {r['n_observed']} | {r['S_chsh']:.4f} | {r['S_vN']:.2f} | {r['purity']:.4f} |\n")
    
    f.write(f"\n**Korelacja (observer_size ↔ S):** r = {corr:.4f}\n\n")
    
    f.write("## 3. INTERPRETACJA\n\n")
    if corr < -0.5:
        f.write("### ✅ POTWIERDZENIE: EMERGENCJA KLASYCZNOŚCI\n\n")
        f.write("Gdy obserwator jest większą częścią systemu:\n")
        f.write("- S maleje → system wygląda bardziej klasycznie\n")
        f.write("- Limit: obserwator = cały system → S = 0 (mierzy sam siebie)\n\n")
        f.write("**WNIOSEK:** Klasyczność jest perspektywą dużego obserwatora.\n")
    else:
        f.write("### ❌ BRAK OCZEKIWANEGO TRENDU\n\n")
        f.write(f"Korelacja r = {corr:.4f} nie potwierdza hipotezy.\n")
        f.write("Możliwe przyczyny:\n")
        f.write("- System zbyt mały (N=10)\n")
        f.write("- K(d) zbyt słabe\n")
    
    f.write("\n## 4. SCHEMAT\n\n")
    f.write("```\n")
    f.write("Mały obserwator (1 oktawa):\n")
    f.write("  [O] [====OBSERVED====]\n")
    f.write("  S = wysokie (widzi dużo kwantowości)\n\n")
    f.write("Średni obserwator (5 oktaw):\n")
    f.write("  [===OBSERVER===] [OBSERVED]\n")
    f.write("  S = średnie\n\n")
    f.write("Duży obserwator (9 oktaw):\n")
    f.write("  [=======OBSERVER=======] [O]\n")
    f.write("  S = niskie (prawie mierzy siebie)\n")
    f.write("```\n")

print(f"   Report saved to: {REPORT_FILE}")

print("\n" + "="*80)
print("QW-687 COMPLETE")
print("="*80)
print(f"\nWYNIKI:")
for r in results:
    print(f"  Observer {r['n_observer']}: S = {r['S_chsh']:.4f}")
print(f"\nKorelacja: r = {corr:.4f}")
print(f"{trend}")
