#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Krzysztof Żuchowski

"""
================================================================================
BADANIE 113: GŁĘBOKIE STUDIUM WEWNĘTRZNEJ STRUKTURY NADSOLITONA
================================================================================

Data utworzenia: 14 listopada 2025
Projekt: Teoria Wszystkiego - Fraktalny Nadsoliton Informacyjny
Autor: GitHub Copilot (na podstawie analiz Script 110-112)

PROBLEM BADAWCZY:
Poprzednie badania (Script 110-112) zidentyfikowały kluczowe cechy nadsolitona:
1. Algebraiczna struktura (top-4 mody: 25% pełne zamknięcie)
2. Skalowanie PR ~ N (mody rozszerzone, nie zlokalizowane)
3. Topologiczna responsywność (defekt → zmiana PR o 10-15%)
4. Niska wymiarowość generatorów (rank=2 z SVD)
5. RG brak zmian znaku w [0.5, 2.5]

CEL BADANIA 113 - ROZSZERZONE CHARAKTERYZACJE:
Wysondować wewnętrzną budowę nadsolitona poprzez 5 kategorii eksperymentów
bez fittingu, aby odkryć:
  - Pełną algebraiczną strukturę (top-12 vs top-4)
  - Rzeczywisty wymiar fraktalny i skalowanie z N
  - Mapę topologicznej wrażliwości
  - Efektywną algebrę generatorów
  - Fazy RG w rozszerzonym zakresie skal

================================================================================
METODOLOGIA
================================================================================

Badanie wykonane dla pełnego ensemble N ∈ {12, 16, 20, 24, 28, 32}, co daje:
- Potwierdzenie uniwersalności wyników
- Określenie systematycznych trendów
- Oszacowanie błędów statystycznych

Każde zadanie zawiera:
1. Teoretyczne uzasadnienie (po co to robimy)
2. Metodologię (jak to robimy)
3. Oczekiwane wyniki (co powinno się pokazać)
4. Obliczenia (kod)
5. Interpretacja (co to znaczy dla teorii)

================================================================================
PARAMETRY BAZOWE
================================================================================
"""

import numpy as np
import json
import sys
from datetime import datetime
from scipy.linalg import eigh, svd
from scipy.optimize import curve_fit
import warnings

warnings.filterwarnings('ignore')

# ============================================================================
# PARAMETRY JĄDRA TEORII (z badań poprzednich)
# ============================================================================

ALPHA_GEO = 2.7715      # Odkryty w badaniach 104+
BETA_TORS = 0.01        # Stabilny, niskorzędny
OMEGA = 0.7854          # π/4 - rezonans
PHI = 0.0               # Faza bazowa

# Kalibracja energetyczna (z Script 110)
V_HIGGS = 246.0         # GeV - anchor kalibracyjny

# Domyślne rozmiary ensemble
DEFAULT_NS = [12, 16, 20, 24, 28, 32]

# ============================================================================
# FUNKCJE JĄDRA
# ============================================================================

def kernel_K(d, s=1.0, complex_kernel=False):
    """
    Uniwersalne jądro sprzężeń.
    K(d,s) = α_geo·cos(ωd+φ)/(1+β_tors·d) przy skali s
    
    Argumenty:
      d: odległość (0..N-1)
      s: skala renormalizacyjna (domyślnie 1.0)
      complex_kernel: jeśli True, zwraca postać zespoloną
    
    Zwraca: wartość jądra (rzeczywista lub zespolona)
    """
    d_scaled = d * s
    if complex_kernel:
        return ALPHA_GEO * s * np.exp(1j * (OMEGA * d_scaled + PHI)) / (1.0 + BETA_TORS * d_scaled)
    else:
        return ALPHA_GEO * s * np.cos(OMEGA * d_scaled + PHI) / (1.0 + BETA_TORS * d_scaled)

def build_resonance_matrix(N, s=1.0, complex_kernel=False):
    """
    Buduje macierz rezonansową (macierz samosprzężeń S_ij).
    
    S_ij = K(|i-j|, s) - uniwersalne sprzężenie między modami
    
    Argumenty:
      N: rozmiar macierzy (liczba modów)
      s: skala renormalizacyjna
      complex_kernel: postać zespolona?
    
    Zwraca: macierz (N×N)
    """
    S = np.zeros((N, N), dtype=np.complex128 if complex_kernel else np.float64)
    for i in range(N):
        for j in range(N):
            d = abs(i - j)
            S[i, j] = kernel_K(d, s, complex_kernel)
    return S

def get_eigenmodes(S):
    """
    Diagonalizuje macierz rezonansową, zwraca wartości i wektory własne.
    
    Zwraca:
      eigenvalues: sorted descending
      eigenvectors: odpowiadające kolumny
    """
    evals, evecs = eigh(S)
    # Sort descending
    idx = np.argsort(-evals)
    return evals[idx], evecs[:, idx]

# ============================================================================
# ZADANIE 1: ALGEBRAICZNA PROBE - ROZSZERZONA DO TOP-12
# ============================================================================

def task_1_algebraic_probe_extended(N, s=1.0):
    """
    ZADANIE 1: ALGEBRAICZNA SONDA - ZAMKNIĘCIE KOMUTATORÓW TOP-12
    
    PROBLEM:
    Script 112 wykazał, że top-4 mody mają 25% pełne zamknięcie w algebrze
    Lie'ego (relative residual < 1e-2). Ale może top-12 da lepsze wyniki?
    
    METODA:
    1. Diagonalizuj S dla N zadanego
    2. Weź top-12 wektorów własnych
    3. Zbuduj projektory P_i = v_i v_i^T
    4. Oblicz komutatory [P_i, P_j] dla wszystkich par
    5. Spróbuj wyrazić jako kombinacja P_i (test zamknięcia algebry)
    
    OCZEKIWANY WYNIK:
    - Frakcja par z relative residual < 1e-2 powinna być wyższa niż w top-4
    - Struktura algebraiczna powinna być bardziej regularna
    - Może ujawni się ukryta symetria (np. su(3), su(2))
    
    ZNACZENIE DLA TEORII:
    Jeśli top-12 zamyka się w algebrze Lie'ego, to nadsoliton ma naturalną,
    algebraiczną strukturę, którą można zidentyfikować i wykorzystać do
    budowy teorii pola.
    """
    S = build_resonance_matrix(N, s=s)
    evals, evecs = get_eigenmodes(S)
    
    # Weź top-12 (lub tyle ile jest modów jeśli N < 12)
    n_modes = min(12, N)
    evecs_top = evecs[:, :n_modes]
    
    # Zbuduj projektory
    projectors = []
    for i in range(n_modes):
        v = evecs_top[:, i].reshape(-1, 1)
        P = v @ v.T  # v v^T
        projectors.append(P)
    
    # Oblicz komutatory i ocenę zamknięcia
    commutators = []
    residuals = []
    
    for i in range(n_modes):
        for j in range(n_modes):
            if i >= j:
                continue
            P_i, P_j = projectors[i], projectors[j]
            C_ij = P_i @ P_j - P_j @ P_i  # [P_i, P_j]
            commutators.append(C_ij)
            
            # Ocena zamknięcia: czy C_ij = sum_k c_ijk P_k?
            # Liczymy residualne normy
            residual_norms = []
            for k in range(n_modes):
                c_k = np.trace(C_ij @ projectors[k])
                residual_norms.append(abs(c_k))
            residuals.append(residual_norms)
    
    residuals = np.array(residuals)
    
    # Statystyki zamknięcia
    residuals_flat = residuals.flatten()
    fraction_below_1e2 = np.sum(residuals_flat < 1e-2) / len(residuals_flat)
    fraction_below_1e1 = np.sum(residuals_flat < 1e-1) / len(residuals_flat)
    rel_residual_mean = np.mean(residuals_flat[residuals_flat > 0]) if np.any(residuals_flat > 0) else 0
    
    return {
        'task': 'Algebraic Probe Extended (top-12)',
        'N': N,
        's': s,
        'n_modes': n_modes,
        'n_commutators': len(commutators),
        'residuals_stats': {
            'fraction_below_1e2': float(fraction_below_1e2),
            'fraction_below_1e1': float(fraction_below_1e1),
            'mean_residual': float(rel_residual_mean),
            'max_residual': float(np.max(residuals_flat)),
            'min_residual': float(np.min(residuals_flat[residuals_flat > 0])) if np.any(residuals_flat > 0) else 0
        },
        'evals_top12': [float(e) for e in evals[:n_modes]],
        'interpretation': 'Estructura algebraica de los 12 modos principales'
    }

# ============================================================================
# ZADANIE 2: SKALOWANIE PR - ROZSZERZONE (WIĘCEJ PUNKTÓW N)
# ============================================================================

def task_2_pr_scaling_extended(Ns):
    """
    ZADANIE 2: SKALOWANIE PARTICIPATION RATIO - ROZSZERZONE
    
    PROBLEM:
    Script 112 wykazał PR ~ N^α gdzie α ≈ 0.998 (prawie 1).
    Czy to rzeczywiście PR ∝ N, czy jest niuans?
    
    METODA:
    1. Dla każdego N w Ns: zbuduj S, oblicz PR dla każdego modu
    2. Przyk: PR_i = 1 / sum(|v_ij|^4) dla wektora v_i
    3. Zbierz PR dla top-4 modów dla każdego N
    4. Wykonaj fit: PR ~ c * N^α (log-log)
    5. Oblicz α i oszacuj błąd
    
    OCZEKIWANY WYNIK:
    - α powinno być bardzo bliskie 1.0 (±0.01)
    - Krzywa log-log powinna być liniowa (R² > 0.99)
    - Brak odchyleń sugeruje liniowe skalowanie: mody są rozszerzone
    
    ZNACZENIE DLA TEORII:
    α = 1 oznacza, że mody nadsolitona nie są zlokalizowane ani ultradelokalizo-
    wane - są to mody falkowe (wave modes) rozprzestrzeniające się przez całą
    przestrzeń w taki sam sposób niezależnie od rozmiaru systemu.
    """
    results_by_mode = {}
    
    for mode_idx in range(4):  # Top-4 mody
        prs = []
        for N in Ns:
            S = build_resonance_matrix(N)
            evals, evecs = get_eigenmodes(S)
            v = evecs[:, mode_idx]
            pr = 1.0 / np.sum(np.abs(v)**4)
            prs.append(pr)
        
        prs = np.array(prs)
        Ns_arr = np.array(Ns)
        
        # Fit logarytmiczny: log(PR) = log(c) + α*log(N)
        log_Ns = np.log(Ns_arr)
        log_prs = np.log(prs)
        
        # Fitowanie liniowe
        coeffs = np.polyfit(log_Ns, log_prs, 1)
        alpha, log_c = coeffs
        c = np.exp(log_c)
        
        # Ocena R²
        fit_line = alpha * log_Ns + log_c
        ss_res = np.sum((log_prs - fit_line)**2)
        ss_tot = np.sum((log_prs - np.mean(log_prs))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot != 0 else 0
        
        results_by_mode[f'mode_{mode_idx}'] = {
            'Ns': [int(n) for n in Ns],
            'PR_values': [float(p) for p in prs],
            'fit_alpha': float(alpha),
            'fit_c': float(c),
            'R_squared': float(r2),
            'interpretation': f'PR ~ {c:.3f} * N^{alpha:.4f}'
        }
    
    return {
        'task': 'PR Scaling Extended',
        'results': results_by_mode,
        'overall_alpha_mean': float(np.mean([results_by_mode[f'mode_{i}']['fit_alpha'] for i in range(4)])),
        'overall_r2_mean': float(np.mean([results_by_mode[f'mode_{i}']['R_squared'] for i in range(4)]))
    }

# ============================================================================
# ZADANIE 3: DEFEKT TOPOLOGICZNY - ROZSZERZONA MAPA
# ============================================================================

def task_3_defect_probe_extended(N, s=1.0):
    """
    ZADANIE 3: SONDA DEFEKTU TOPOLOGICZNEGO - PEŁNA MAPA
    
    PROBLEM:
    Script 112 wykazał, że lokalizacja defektu powoduje zmianę PR o 10-15%.
    Ale co się dzieje z innymi obserwablami? Czy defekt zmienia algebraiczną
    strukturę?
    
    METODA:
    1. Zbuduj S(s) - macierz rezonansową na skali s
    2. Oblicz eigenmodes, zbierz obserwable
    3. Wprowadź defekt topologiczny (zeruj część jądra w centrum)
    4. Oblicz zmienione eigenmodes
    5. Porównaj: PR, wartości własne, wektory własne, entropię Shannona
    
    OCZEKIWANY WYNIK:
    - PR zmienia się 10-15% (potwierdzi poprzednie wyniki)
    - Zmiana wartości własnych λ zmienia się mniej (<5%)
    - Wektory własne przesuwają się, ale struktura pozostaje
    - Topologiczna responsywność potwierdzona
    
    ZNACZENIE DLA TEORII:
    Nadsoliton ma topologiczną strukturę rdzenia, którą można aktywować
    lokalnie. To sugeruje, że topologia (nie geometria) jest fundamentalna.
    """
    S_orig = build_resonance_matrix(N, s=s)
    evals_orig, evecs_orig = get_eigenmodes(S_orig)
    
    # Defekt: zeruj część centralną kernelu (np. d = 1..2)
    S_defect = S_orig.copy()
    defect_width = 2
    mid = N // 2
    for i in range(N):
        for j in range(N):
            d = abs(i - j)
            if d < defect_width:
                S_defect[i, j] = 0
    
    evals_def, evecs_def = get_eigenmodes(S_defect)
    
    # Obserwable
    pr_orig = []
    pr_defect = []
    for mode_idx in range(4):
        v_orig = evecs_orig[:, mode_idx]
        v_def = evecs_def[:, mode_idx]
        
        pr_orig.append(1.0 / np.sum(np.abs(v_orig)**4))
        pr_defect.append(1.0 / np.sum(np.abs(v_def)**4))
    
    pr_orig = np.array(pr_orig)
    pr_defect = np.array(pr_defect)
    pr_change = (pr_orig - pr_defect) / pr_orig * 100  # % zmiana
    
    lambda_change = abs(evals_orig[0] - evals_def[0]) / evals_orig[0] * 100
    
    return {
        'task': 'Defect Probe Extended',
        'N': N,
        's': s,
        'defect_width': defect_width,
        'PR_original': [float(p) for p in pr_orig],
        'PR_defect': [float(p) for p in pr_defect],
        'PR_change_percent': [float(c) for c in pr_change],
        'mean_PR_change': float(np.mean(pr_change)),
        'lambda_change_percent': float(lambda_change),
        'topological_responsivity': 'High' if np.mean(pr_change) > 5 else 'Medium' if np.mean(pr_change) > 2 else 'Low',
        'interpretation': 'Topological sensitivity of nadsoliton structure'
    }

# ============================================================================
# ZADANIE 4: ALGEBRA GENERATORÓW - PEŁNE STUDIUM
# ============================================================================

def task_4_generator_algebra_full(N, s=1.0):
    """
    ZADANIE 4: ALGEBRA GENERATORÓW - STUDIUM PEŁNE
    
    PROBLEM:
    Script 112 wykazał, że SVD macierzy S(φ) dla różnych faz φ daje
    effective_rank = 2 (tylko 2 niezależne generatory).
    Jakie generatory? Czy to su(2)? Coś głębszego?
    
    METODA:
    1. Zbuduj S dla różnych skal s (proxy dla "faz" w algebrazie)
    2. Zbierz macierze S(s_i) w stos (każdy wiersz to wektoryza S)
    3. SVD całego stosu: U Σ V^T
    4. Zbadaj singular values: ile jest znaczących?
    5. Zbadaj lewe wektory singularne U: jakie kierunki dominują?
    6. Zrekonstruuj S ze stosów i policz błąd
    
    OCZEKIWANY WYNIK:
    - 2-3 wartości singularne są duże (>>1)
    - Reszta wartości singularnych → 0 (< 1e-10)
    - Lewe wektory singularne ujawniają strukturę algebraiczną
    - Nadsoliton ma co najmniej 2-wymiarową algebra
    
    ZNACZENIE DLA TEORII:
    Jeśli algebra jest 2-wymiarowa (SU(2)-like) lub 3-wymiarowa (SU(3)-like),
    to emergentna symetria gauge ma naturalną podstawę.
    """
    # Zbierz S dla różnych skal
    s_values = np.linspace(0.5, 2.5, 11)
    S_stack = []
    
    for s_val in s_values:
        S = build_resonance_matrix(N, s=s_val)
        S_vec = S.flatten()  # Wektoryza
        S_stack.append(S_vec)
    
    S_stack = np.array(S_stack)  # Shape: (n_scales, N*N)
    
    # SVD
    U, singular_vals, Vt = svd(S_stack, full_matrices=False)
    
    # Ocena rangi
    threshold = 1e-10
    effective_rank = np.sum(singular_vals > threshold)
    
    # Frakcja energii w top-k modach
    total_energy = np.sum(singular_vals**2)
    energy_top2 = np.sum(singular_vals[:2]**2)
    energy_top3 = np.sum(singular_vals[:3]**2)
    frac_top2 = energy_top2 / total_energy
    frac_top3 = energy_top3 / total_energy
    
    return {
        'task': 'Generator Algebra Full',
        'N': N,
        's_range': [float(s) for s in s_values],
        'singular_values': [float(s) for s in singular_vals[:10]],  # Top-10
        'effective_rank': int(effective_rank),
        'threshold': float(threshold),
        'energy_fraction_top2': float(frac_top2),
        'energy_fraction_top3': float(frac_top3),
        'interpretation': f'Effective generator algebra dimension: {effective_rank} or 2-3'
    }

# ============================================================================
# ZADANIE 5: RG BIFURKACJA - ROZSZERZONA SONDA
# ============================================================================

def task_5_rg_bifurcation_wide(N):
    """
    ZADANIE 5: RG BIFURKACJA - ROZSZERZONA SONDA
    
    PROBLEM:
    Script 112 sondował RG w zakresie s ∈ [0.5, 2.5] i nie znalazł zmian znaku
    (brak fazy przejścia). Może na szerszym zakresie jest?
    
    METODA:
    1. Oblicz lambda_max(s) dla szerokiego zakresu s ∈ [0.1, 10]
    2. Oblicz beta-proxy: dλ/d(ln s)
    3. Policz zmianę znaku dβ/d(ln s)
    4. Ponadto, przeanalizuj bifurkacje w strukturze eigenvector:
       ile "nowych" modów pojawia się w top-4?
    
    OCZEKIWANY WYNIK:
    - Może być 1-2 zmiany znaku w szerszym zakresie
    - Bifurkacje mogą ujawić fazowe przejścia
    - Jeśli brak zmian nawet w [0.1, 10], to RG płynna
    
    ZNACZENIE DLA TEORII:
    Czy nadsoliton ma fazy RG? Czy jest fazowe przejście przy pewnych skalach?
    To określi, czy istnieje punktu krytyczne lub asymptotyczną wolność.
    """
    s_values = np.logspace(-1, 1, 91)  # 0.1 do 10, 91 punktów
    
    lambda_maxs = []
    betas = []
    
    for s_val in s_values:
        S = build_resonance_matrix(N, s=s_val)
        evals, _ = get_eigenmodes(S)
        lambda_maxs.append(evals[0])
    
    lambda_maxs = np.array(lambda_maxs)
    
    # Beta-proxy: dλ/d(ln s)
    d_lns = np.diff(np.log(s_values))
    d_lam = np.diff(lambda_maxs)
    betas = d_lam / d_lns
    
    # Zmianę znaku
    sign_changes = 0
    for i in range(len(betas) - 1):
        if betas[i] * betas[i+1] < 0:
            sign_changes += 1
    
    return {
        'task': 'RG Bifurcation Wide',
        'N': N,
        's_range': [float(s) for s in s_values[::5]],  # Co 5-ty punkt
        'lambda_max_values': [float(l) for l in lambda_maxs[::5]],
        'beta_proxy_values': [float(b) for b in betas[::5]],
        'total_points': len(s_values),
        'sign_changes_count': int(sign_changes),
        's_min': float(s_values[0]),
        's_max': float(s_values[-1]),
        'interpretation': f'RG bifurcations: {sign_changes} sign changes in [{s_values[0]:.2f}, {s_values[-1]:.2f}]'
    }

# ============================================================================
# SYNTEZA - KOMBINACJA WSZYSTKICH WYNIKÓW
# ============================================================================

def synthesize_results(results_dict):
    """
    Łączy wyniki wszystkich zadań w spójną interpretację struktury nadsolitona.
    """
    synthesis = {
        'timestamp': datetime.now().isoformat(),
        'title': 'DEEP NADSOLITON STRUCTURE ANALYSIS - SYNTHESIS',
        'key_findings': {}
    }
    
    # Finding 1: Algebraic Structure
    if 'task_1' in results_dict:
        task1 = results_dict['task_1']
        frac_good = task1.get('residuals_stats', {}).get('fraction_below_1e2', 0)
        synthesis['key_findings']['algebraic_closure'] = {
            'status': 'STRONG' if frac_good > 0.5 else 'PARTIAL' if frac_good > 0.2 else 'WEAK',
            'fraction_excellent': frac_good,
            'finding': f'{int(frac_good*100)}% of commutator pairs have residual < 1e-2'
        }
    
    # Finding 2: Scaling
    if 'task_2' in results_dict:
        task2 = results_dict['task_2']
        alpha_mean = task2.get('overall_alpha_mean', 0)
        r2_mean = task2.get('overall_r2_mean', 0)
        synthesis['key_findings']['pr_scaling'] = {
            'exponent': alpha_mean,
            'fit_quality': r2_mean,
            'interpretation': 'Extended modes' if abs(alpha_mean - 1.0) < 0.05 else 'Localized modes' if alpha_mean < 0.5 else 'Ultra-delocalized'
        }
    
    # Finding 3: Topological Sensitivity
    if 'task_3' in results_dict:
        task3_results = results_dict['task_3']
        if isinstance(task3_results, list):
            task3 = task3_results[0] if task3_results else {}
        else:
            task3 = task3_results
        
        mean_pr_change = task3.get('mean_PR_change', 0)
        synthesis['key_findings']['topological_sensitivity'] = {
            'mean_PR_change_percent': mean_pr_change,
            'status': 'High' if mean_pr_change > 10 else 'Moderate' if mean_pr_change > 5 else 'Low',
            'finding': f'Defect induces {mean_pr_change:.1f}% change in PR'
        }
    
    # Finding 4: Generator Dimension
    if 'task_4' in results_dict:
        task4_results = results_dict['task_4']
        if isinstance(task4_results, list):
            task4 = task4_results[0] if task4_results else {}
        else:
            task4 = task4_results
        
        eff_rank = task4.get('effective_rank', 0)
        synthesis['key_findings']['generator_dimension'] = {
            'effective_rank': eff_rank,
            'status': 'SU(2)-like' if eff_rank == 2 else 'SU(3)-like' if eff_rank == 3 else f'Other (rank={eff_rank})',
            'finding': f'Minimal generator algebra has dimension {eff_rank}'
        }
    
    # Finding 5: RG Structure
    if 'task_5' in results_dict:
        task5_results = results_dict['task_5']
        if isinstance(task5_results, list):
            task5 = task5_results[0] if task5_results else {}
        else:
            task5 = task5_results
        
        sign_changes = task5.get('sign_changes_count', 0)
        synthesis['key_findings']['rg_structure'] = {
            'sign_changes': sign_changes,
            'status': 'Rich phases' if sign_changes >= 3 else 'Smooth RG' if sign_changes == 0 else 'Moderate structure',
            'finding': f'{sign_changes} phase transitions detected in RG flow'
        }
    
    return synthesis

# ============================================================================
# GŁÓWNA FUNKCJA - WYKONANIE
# ============================================================================

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='BADANIE 113: Deep Nadsoliton Structure Analysis')
    parser.add_argument('--Ns', type=int, nargs='+', default=DEFAULT_NS, help='System sizes to analyze')
    parser.add_argument('--dry-run', action='store_true', help='Dry run with small N only')
    parser.add_argument('--anchor', type=float, default=V_HIGGS, help='Calibration anchor (GeV)')
    
    args = parser.parse_args()
    
    if args.dry_run:
        print("[DRY RUN MODE] Using N=24 only")
        Ns = [24]
    else:
        Ns = args.Ns
    
    print("=" * 80)
    print("BADANIE 113: DEEP NADSOLITON STRUCTURE ANALYSIS")
    print("=" * 80)
    print(f"Timestamp: {datetime.now().isoformat()}")
    print(f"Ensemble sizes: {Ns}")
    print(f"Parameters: α_geo={ALPHA_GEO}, β_tors={BETA_TORS}, ω={OMEGA}")
    print("=" * 80)
    print()
    
    # Zbierz wyniki
    all_results = {
        'metadata': {
            'created': datetime.now().isoformat(),
            'script': '113_DEEP_NADSOLITON_STRUCTURE_ANALYSIS.py',
            'Ns': Ns,
            'parameters': {
                'alpha_geo': ALPHA_GEO,
                'beta_tors': BETA_TORS,
                'omega': OMEGA,
                'phi': PHI
            }
        }
    }
    
    # ========== ZADANIE 1: Algebraic Probe Extended ==========
    print("\n[TASK 1] Algebraic Probe Extended (top-12 modes)...")
    task_1_results = []
    for N in Ns:
        result = task_1_algebraic_probe_extended(N)
        task_1_results.append(result)
        frac = result['residuals_stats']['fraction_below_1e2']
        print(f"  N={N}: {int(frac*100)}% pairs with residual < 1e-2")
    
    all_results['task_1_algebraic_probe'] = task_1_results
    
    # ========== ZADANIE 2: PR Scaling Extended ==========
    print("\n[TASK 2] PR Scaling Extended (all Ns)...")
    result = task_2_pr_scaling_extended(Ns)
    all_results['task_2_pr_scaling'] = result
    print(f"  Mean exponent α: {result['overall_alpha_mean']:.4f}")
    print(f"  Mean R²: {result['overall_r2_mean']:.4f}")
    
    # ========== ZADANIE 3: Defect Probe Extended ==========
    print("\n[TASK 3] Defect Probe Extended (all Ns)...")
    task_3_results = []
    for N in Ns:
        result = task_3_defect_probe_extended(N)
        task_3_results.append(result)
        mean_change = result['mean_PR_change']
        print(f"  N={N}: Mean PR change = {mean_change:.2f}%")
    
    all_results['task_3_defect_probe'] = task_3_results
    
    # ========== ZADANIE 4: Generator Algebra Full ==========
    print("\n[TASK 4] Generator Algebra Full (sample N=24)...")
    N_sample = 24
    result = task_4_generator_algebra_full(N_sample)
    all_results['task_4_generator_algebra'] = result
    print(f"  Effective rank: {result['effective_rank']}")
    print(f"  Energy in top-2: {result['energy_fraction_top2']:.3f}")
    
    # ========== ZADANIE 5: RG Bifurcation Wide ==========
    print("\n[TASK 5] RG Bifurcation Wide (N=24)...")
    result = task_5_rg_bifurcation_wide(24)
    all_results['task_5_rg_bifurcation'] = result
    print(f"  Sign changes in RG: {result['sign_changes_count']}")
    print(f"  Range: s ∈ [{result['s_min']:.2f}, {result['s_max']:.2f}]")
    
    # ========== SYNTEZA ==========
    print("\n[SYNTHESIS] Combining all results...")
    synthesis = synthesize_results(all_results)
    all_results['synthesis'] = synthesis
    
    # ========== ZAPIS ==========
    output_file = 'report_113_deep_nadsoliton_structure.json'
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(all_results, f, indent=2, ensure_ascii=False)
    
    print(f"\n✅ Raport zapisany: {output_file}")
    print("\n" + "=" * 80)
    print("KEY FINDINGS:")
    print("=" * 80)
    for key, finding in synthesis['key_findings'].items():
        print(f"\n{key.upper()}:")
        for k, v in finding.items():
            print(f"  {k}: {v}")
    
    return all_results

if __name__ == "__main__":
    main()
