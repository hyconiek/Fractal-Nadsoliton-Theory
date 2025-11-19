#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Krzysztof Żuchowski

"""
================================================================================
BADANIE 114 v2: ZAAWANSOWANE MAPOWANIE GENERATORÓW NA OBSERWABLE
================================================================================

PROBLEM W V1:
Naiwne mapowanie eigenvalues → masy dało błędy ~97%.
To sygnalizuje, że generatory to nie bezpośrednie wartości własne,
lecz STRUKTURY ALGEBRAICZNE (komutatory, reprezentacje SU(N)).

ROZWIĄZANIE V2:
1. Analiza macierzy komutatora [G_i, G_j] (algebra Liego)
2. Ekstraktuj STRUKTURALNE INVARIANTY (Casimir operators)
3. Stosunki mas z HERMITOWSKIEGO TRACJU generatorów
4. Coupling ratios z NORMY operatorów
5. CKM angles z STRUKTURY REPREZENTACJI

KLUCZOWA IDEA:
W algebrze Liego, masy cząstek odpowiadają EIGENVALUOM CASIMIRA,
nie eigenvaluom macierzy S bezpośrednio!

================================================================================
"""

import numpy as np
import json
from datetime import datetime
from scipy.linalg import eigh, norm
import warnings

warnings.filterwarnings('ignore')

# Parametry kernelu
ALPHA_GEO = 2.7715
BETA_TORS = 0.01
OMEGA = 0.7854
PHI = 0.0

# Kalibracja
V_HIGGS = 246.0  # GeV
TOP_MASS = 172.76  # GeV

# Obserwable eksperymentalne
EXPERIMENTAL = {
    'electron_mass_MeV': 0.511,
    'muon_mass_MeV': 105.66,
    'tau_mass_MeV': 1776.86,
    'W_mass_GeV': 80.377,
    'Z_mass_GeV': 91.188,
    'Higgs_mass_GeV': 125.10,
    'sin2_theta_W': 0.2223,
    'alpha_s_Mz': 0.118,
    'theta_12_rad': 0.2261,
    'theta_23_rad': 0.0403,
    'theta_13_rad': 0.00365,
}

DEFAULT_NS = [12, 16, 20, 24, 28, 32]

# ============================================================================
# FUNKCJE JĄDRA
# ============================================================================

def kernel_K(d, s=1.0):
    """Uniwersalne jądro sprzężeń."""
    d_scaled = d * s
    return ALPHA_GEO * s * np.cos(OMEGA * d_scaled + PHI) / (1.0 + BETA_TORS * d_scaled)

def build_resonance_matrix(N, s=1.0):
    """Zbuduj macierz rezonansową S_ij."""
    S = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            d = abs(i - j)
            S[i, j] = kernel_K(d, s)
    return S

def get_eigenmodes(S):
    """Diagonalizuj S."""
    evals, evecs = eigh(S)
    idx = np.argsort(-evals)
    return evals[idx], evecs[:, idx]

# ============================================================================
# ANALIZA ALGEBRY LIEGO
# ============================================================================

def extract_generators_as_matrices(N, s_values):
    """
    Wygeneruj 11 generatorów algebry Liego poprzez PCA
    macierzy S na różnych skalach.
    
    Zwraca:
      generators: lista N×N macierzy (generatory)
      singular_values: wartości singularne (energia)
    """
    matrices_stack = []
    
    for s in s_values:
        S = build_resonance_matrix(N, s=s)
        matrices_stack.append(S.flatten())
    
    matrices_stack = np.array(matrices_stack)
    
    # SVD na stacku macierzy
    U, sing_vals, Vt = np.linalg.svd(matrices_stack, full_matrices=False)
    
    # Każdy wektor singularny V to "generator" (kierunek w przestrzeni macierzy)
    generators = []
    for i in range(min(11, len(Vt))):
        gen_matrix = Vt[i].reshape((N, N))
        gen_matrix = (gen_matrix + gen_matrix.T) / 2  # Hermitian
        generators.append(gen_matrix)
    
    return generators, sing_vals

def compute_commutator(G_i, G_j):
    """Komutator [G_i, G_j] = G_i @ G_j - G_j @ G_i"""
    return G_i @ G_j - G_j @ G_i

def compute_casimir_invariant(N, generators):
    """
    Casimir invariant C = Σ G_i²
    Odzwierciedla łączną "masę" algebry.
    """
    C = np.zeros((N, N))
    for G in generators:
        C += G @ G
    return C

def compute_generator_norms(generators):
    """Norma Frobeniusa każdego generatora (energia)."""
    norms = [norm(G, 'fro') for G in generators]
    return np.array(norms)

def extract_mass_hierarchy_from_casimir(N, generators):
    """
    Eigenvalues Casimira indukują hierarchię mas.
    
    IDEA:
    m_i ∝ λ_i(Casimir) — większy eigenvalue → większa masa
    """
    C = compute_casimir_invariant(N, generators)
    evals_casimir, _ = eigh(C)
    evals_casimir = np.sort(-evals_casimir)  # Descending
    
    return evals_casimir[:min(12, len(evals_casimir))]

# ============================================================================
# ZADANIE 1 v2: Algebra Characteristics
# ============================================================================

def task_1_algebra_characteristics(N):
    """ZADANIE 1 v2: Charakterystyka algebry generatorów"""
    s_values = np.linspace(0.5, 2.5, 11)
    generators, sing_vals = extract_generators_as_matrices(N, s_values)
    
    # Norms generatorów
    gen_norms = compute_generator_norms(generators)
    total_norm = np.sum(gen_norms**2)
    frac_top3 = np.sum(gen_norms[:3]**2) / total_norm if total_norm > 0 else 0
    
    # Casimir eigenvalues (masa charakterystyczna algebry)
    mass_eigenvalues = extract_mass_hierarchy_from_casimir(N, generators)
    
    return {
        'task': 'Algebra Characteristics',
        'N': N,
        'n_generators': len(generators),
        'generator_norms_top5': [float(n) for n in gen_norms[:5]],
        'casimir_eigenvalues_top5': [float(e) for e in mass_eigenvalues[:5]],
        'energy_concentration': float(frac_top3),
        'casimir_trace': float(np.trace(compute_casimir_invariant(N, generators)))
    }

# ============================================================================
# ZADANIE 2 v2: Lepton Mass Ratios (via Casimir)
# ============================================================================

def task_2_lepton_ratios_v2(N):
    """ZADANIE 2 v2: Stosunki mas leptonów z Casimira"""
    s_values = np.linspace(0.5, 2.5, 11)
    generators, _ = extract_generators_as_matrices(N, s_values)
    mass_evals = extract_mass_hierarchy_from_casimir(N, generators)
    
    # Top-3 eigenvalues Casimira → masy leptonów (e, μ, τ)
    # ZAŁOŻENIE: większy eigenvalue → większa masa
    # Normalizujemy do m_e = 1 (bezwymiarowe)
    
    if len(mass_evals) < 3:
        return {
            'error': f'Not enough Casimir eigenvalues (only {len(mass_evals)} found)',
            'N': N
        }
    
    m_e_proxy = mass_evals[2]  # Najmniejszy (electron)
    m_mu_proxy = mass_evals[1]  # Średni (muon)
    m_tau_proxy = mass_evals[0]  # Największy (tau)
    
    # Stosunki
    theory_mu_e = m_mu_proxy / m_e_proxy if m_e_proxy != 0 else 0
    theory_tau_e = m_tau_proxy / m_e_proxy if m_e_proxy != 0 else 0
    
    # Eksperyment
    exp_mu_e = EXPERIMENTAL['muon_mass_MeV'] / EXPERIMENTAL['electron_mass_MeV']
    exp_tau_e = EXPERIMENTAL['tau_mass_MeV'] / EXPERIMENTAL['electron_mass_MeV']
    
    error_mu_e = abs(theory_mu_e - exp_mu_e) / exp_mu_e * 100 if exp_mu_e != 0 else 100
    error_tau_e = abs(theory_tau_e - exp_tau_e) / exp_tau_e * 100 if exp_tau_e != 0 else 100
    
    return {
        'task': 'Lepton Mass Ratios (Casimir)',
        'N': N,
        'theory': {
            'mu_e_ratio': float(theory_mu_e),
            'tau_e_ratio': float(theory_tau_e),
        },
        'experiment': {
            'mu_e_ratio': float(exp_mu_e),
            'tau_e_ratio': float(exp_tau_e),
        },
        'errors_percent': {
            'mu_e_error': float(error_mu_e),
            'tau_e_error': float(error_tau_e),
            'mean_error': float((error_mu_e + error_tau_e) / 2),
        }
    }

# ============================================================================
# ZADANIE 3 v2: Boson Mass Ratios
# ============================================================================

def task_3_boson_ratios_v2(N):
    """ZADANIE 3 v2: M_Z/M_W z algebry"""
    s_values = np.linspace(0.5, 2.5, 11)
    generators, sing_vals = extract_generators_as_matrices(N, s_values)
    
    # Singular values mogą reprezentować energie bozonów
    if len(sing_vals) < 2:
        return {'error': 'Not enough singular values', 'N': N}
    
    # Ratio z singular values
    theory_z_w = sing_vals[1] / sing_vals[0] if sing_vals[0] != 0 else 0
    
    # Eksperyment
    exp_z_w = EXPERIMENTAL['Z_mass_GeV'] / EXPERIMENTAL['W_mass_GeV']
    error_z_w = abs(theory_z_w - exp_z_w) / exp_z_w * 100 if exp_z_w != 0 else 100
    
    return {
        'task': 'Boson Mass Ratios',
        'N': N,
        'theory_z_w': float(theory_z_w),
        'experiment_z_w': float(exp_z_w),
        'error_percent': float(error_z_w),
    }

# ============================================================================
# ZADANIE 4 v2: Coupling Ratios z norm generatorów
# ============================================================================

def task_4_coupling_ratios_v2(N):
    """ZADANIE 4 v2: Sprzężenia gauge z norm generatorów"""
    s_values = np.linspace(0.5, 2.5, 11)
    generators, _ = extract_generators_as_matrices(N, s_values)
    
    gen_norms = compute_generator_norms(generators)
    
    if len(gen_norms) < 3:
        return {'error': 'Not enough generators', 'N': N}
    
    # Top-3 normy → sprzężenia (g₃ > g₂ > g₁)
    g3_proxy = gen_norms[0]  # Największy
    g2_proxy = gen_norms[1]
    g1_proxy = gen_norms[2]
    
    # Normalize do rozsądnego zakresu
    g3_theory = g3_proxy
    g2_theory = g2_proxy
    g1_theory = g1_proxy
    
    # Eksperyment (rough)
    exp_g3 = 1.22
    exp_g2 = 0.653
    exp_g1 = 0.364
    
    ratio_theory_g3_g2 = g3_theory / g2_theory if g2_theory != 0 else 0
    ratio_theory_g2_g1 = g2_theory / g1_theory if g1_theory != 0 else 0
    
    ratio_exp_g3_g2 = exp_g3 / exp_g2
    ratio_exp_g2_g1 = exp_g2 / exp_g1
    
    error_g3_g2 = abs(ratio_theory_g3_g2 - ratio_exp_g3_g2) / ratio_exp_g3_g2 * 100
    error_g2_g1 = abs(ratio_theory_g2_g1 - ratio_exp_g2_g1) / ratio_exp_g2_g1 * 100
    
    return {
        'task': 'Coupling Ratios',
        'N': N,
        'theory': {
            'g3_g2_ratio': float(ratio_theory_g3_g2),
            'g2_g1_ratio': float(ratio_theory_g2_g1),
        },
        'experiment': {
            'g3_g2_ratio': float(ratio_exp_g3_g2),
            'g2_g1_ratio': float(ratio_exp_g2_g1),
        },
        'errors_percent': {
            'g3_g2_error': float(error_g3_g2),
            'g2_g1_error': float(error_g2_g1),
            'mean_error': float((error_g3_g2 + error_g2_g1) / 2),
        }
    }

# ============================================================================
# ZADANIE 5 v2: CKM Angles z komutatora
# ============================================================================

def task_5_ckm_angles_v2(N):
    """ZADANIE 5 v2: CKM angles z struktury komutatora"""
    s_values = np.linspace(0.5, 2.5, 11)
    generators, _ = extract_generators_as_matrices(N, s_values)
    
    if len(generators) < 3:
        return {'error': 'Not enough generators for CKM', 'N': N}
    
    # Komutatory [G_i, G_j]
    comm_01 = compute_commutator(generators[0], generators[1])
    comm_12 = compute_commutator(generators[1], generators[2])
    comm_02 = compute_commutator(generators[0], generators[2])
    
    # Normy komutatora (miara sprzężenia między generatorami)
    norm_comm_01 = norm(comm_01, 'fro')
    norm_comm_12 = norm(comm_12, 'fro')
    norm_comm_02 = norm(comm_02, 'fro')
    
    # Kąty CKM z normalizowanych norm
    # ZAŁOŻENIE: norm(commutator) ∝ sin(θ_CKM)
    sin_theta_12 = norm_comm_01 / (norm_comm_01 + norm_comm_12 + norm_comm_02 + 1e-10)
    sin_theta_23 = norm_comm_12 / (norm_comm_01 + norm_comm_12 + norm_comm_02 + 1e-10)
    sin_theta_13 = norm_comm_02 / (norm_comm_01 + norm_comm_12 + norm_comm_02 + 1e-10)
    
    # Konwertuj na radians (rough mapping)
    theta_12 = np.arcsin(np.clip(sin_theta_12, 0, 1))
    theta_23 = np.arcsin(np.clip(sin_theta_23, 0, 1))
    theta_13 = np.arcsin(np.clip(sin_theta_13, 0, 1))
    
    exp_theta_12 = EXPERIMENTAL['theta_12_rad']
    exp_theta_23 = EXPERIMENTAL['theta_23_rad']
    exp_theta_13 = EXPERIMENTAL['theta_13_rad']
    
    return {
        'task': 'CKM Angles (Commutator)',
        'N': N,
        'theory_rad': {
            'theta_12': float(theta_12),
            'theta_23': float(theta_23),
            'theta_13': float(theta_13),
        },
        'experiment_rad': {
            'theta_12': float(exp_theta_12),
            'theta_23': float(exp_theta_23),
            'theta_13': float(exp_theta_13),
        },
        'note': 'Commutator structure provides complementary CKM prediction'
    }

# ============================================================================
# ZADANIE 6 v2: Synthesis
# ============================================================================

def task_6_synthesis_v2(all_results):
    """ZADANIE 6 v2: Oceń zgodność z eksperymentem"""
    errors_collected = []
    
    if 'task_2' in all_results:
        for res in all_results['task_2']:
            if 'errors_percent' in res:
                mean_err = res['errors_percent'].get('mean_error', 0)
                errors_collected.append(mean_err)
    
    if 'task_3' in all_results:
        for res in all_results['task_3']:
            if 'error_percent' in res:
                errors_collected.append(res['error_percent'])
    
    if 'task_4' in all_results:
        for res in all_results['task_4']:
            if 'errors_percent' in res:
                mean_err = res['errors_percent'].get('mean_error', 0)
                errors_collected.append(mean_err)
    
    synthesis = {
        'timestamp': datetime.now().isoformat(),
        'methodology': 'Casimir invariants + Commutator analysis',
    }
    
    if errors_collected:
        mean_error_total = np.mean(errors_collected)
        synthesis['mean_error_percent'] = float(mean_error_total)
        synthesis['assessment'] = (
            'EXCELLENT!' if mean_error_total < 10 else
            'VERY GOOD' if mean_error_total < 20 else
            'GOOD' if mean_error_total < 40 else
            'MODERATE - needs refinement' if mean_error_total < 60 else
            'Poor fit - different theoretical framework needed'
        )
    
    return synthesis

# ============================================================================
# GŁÓWNA
# ============================================================================

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='BADANIE 114 v2: Advanced Generator Mapping')
    parser.add_argument('--Ns', type=int, nargs='+', default=DEFAULT_NS)
    
    args = parser.parse_args()
    
    print("=" * 80)
    print("BADANIE 114 v2: ZAAWANSOWANE MAPOWANIE GENERATORÓW NA OBSERWABLE")
    print("=" * 80)
    print(f"Metoda: Casimir Invariants + Commutator Analysis")
    print(f"Ensemble: N ∈ {args.Ns}")
    print("=" * 80)
    print()
    
    all_results = {
        'metadata': {
            'created': datetime.now().isoformat(),
            'script': '114_GENERATOR_OBSERVABLE_MAPPING_v2.py',
            'Ns': args.Ns,
            'methodology': 'Casimir invariants + Commutator structure',
        }
    }
    
    # TASK 1
    print("[TASK 1] Algebra Characteristics...")
    task_1_results = []
    for N in args.Ns:
        result = task_1_algebra_characteristics(N)
        task_1_results.append(result)
        print(f"  N={N}: {result['n_generators']} generators, energy concentration {result['energy_concentration']:.1%}")
    all_results['task_1_algebra_characteristics'] = task_1_results
    
    # TASK 2
    print("\n[TASK 2] Lepton Mass Ratios (Casimir)...")
    task_2_results = []
    for N in args.Ns:
        result = task_2_lepton_ratios_v2(N)
        task_2_results.append(result)
        if 'errors_percent' in result:
            print(f"  N={N}: m_μ/m_e error {result['errors_percent']['mu_e_error']:.1f}%, m_τ/m_e error {result['errors_percent']['tau_e_error']:.1f}%")
        else:
            print(f"  N={N}: {result.get('error', 'unknown error')}")
    all_results['task_2_lepton_ratios'] = task_2_results
    
    # TASK 3
    print("\n[TASK 3] Boson Mass Ratios...")
    task_3_results = []
    for N in args.Ns:
        result = task_3_boson_ratios_v2(N)
        task_3_results.append(result)
        if 'error_percent' in result:
            print(f"  N={N}: M_Z/M_W error {result['error_percent']:.1f}%")
        else:
            print(f"  N={N}: {result.get('error', 'unknown error')}")
    all_results['task_3_boson_ratios'] = task_3_results
    
    # TASK 4
    print("\n[TASK 4] Coupling Ratios...")
    task_4_results = []
    for N in args.Ns:
        result = task_4_coupling_ratios_v2(N)
        task_4_results.append(result)
        if 'errors_percent' in result:
            print(f"  N={N}: g₃/g₂ error {result['errors_percent']['g3_g2_error']:.1f}%, g₂/g₁ error {result['errors_percent']['g2_g1_error']:.1f}%")
        else:
            print(f"  N={N}: {result.get('error', 'unknown error')}")
    all_results['task_4_coupling_ratios'] = task_4_results
    
    # TASK 5
    print("\n[TASK 5] CKM Angles...")
    task_5_results = []
    for N in args.Ns:
        result = task_5_ckm_angles_v2(N)
        task_5_results.append(result)
        print(f"  N={N}: CKM angles computed")
    all_results['task_5_ckm_angles'] = task_5_results
    
    # TASK 6
    print("\n[TASK 6] Synthesis...")
    synthesis = task_6_synthesis_v2(all_results)
    all_results['task_6_synthesis'] = synthesis
    print(f"  Assessment: {synthesis.get('assessment', 'unknown')}")
    if 'mean_error_percent' in synthesis:
        print(f"  Mean error: {synthesis['mean_error_percent']:.2f}%")
    
    # SAVE
    output_file = 'report_114_v2_advanced_mapping.json'
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(all_results, f, indent=2)
    
    print(f"\n✅ Raport: {output_file}")
    print("=" * 80)

if __name__ == "__main__":
    main()
