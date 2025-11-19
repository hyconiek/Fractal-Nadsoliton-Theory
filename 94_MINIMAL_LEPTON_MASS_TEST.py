# Author: Krzysztof Żuchowski


import numpy as np
import json
from datetime import datetime

# ============================================================================
# 4 Fundamentalne Parametry (zdefiniowane w teorii)
# ============================================================================
ALPHA_GEO = 1.0
BETA_TORS = 0.1
OMEGA = 0.7854  # rad
PHI = 0.5236    # rad

# 8 Efektywnych Oktaw (z Badania 66/116)
OCTAVES_EFFECTIVE = np.array([1, 3, 4, 6, 7, 9, 10, 12])

# SM Lepton Masses (GeV) for comparison
SM_MASSES = {
    'electron': 0.511e-3,
    'muon': 105.7e-3,
    'tau': 1776.9e-3,
}

# ============================================================================
# FUNKCJE
# ============================================================================

def kernel_K(d: int) -> float:
    """Universal coupling kernel K(d) - podstawowy element interakcji."""
    numerator = ALPHA_GEO * np.cos(OMEGA * d + PHI)
    denominator = 1.0 + BETA_TORS * d
    return numerator / denominator

def construct_static_hamiltonian_matrix(n_modes: int) -> np.ndarray:
    """
    Konstrukcja efektywnego Hamiltonianu macierzowego H_Lepton.
    
    H_ij = E_i δ_ij + K_ij (1 - δ_ij)
    - Diagonalne E_i: Energia własna oktawy (odpowiada masie m_o²).
    - Pozadiagonalne K_ij: Sprzężenie międzyoktawowe.
    """
    H = np.zeros((n_modes, n_modes))
    
    # 1. Elementy diagonalne (Masa Własna Oktawy)
    # Używamy K(d) jako heurystycznej miary energii/masy dla danej oktawy.
    # Przez brak VEV/Higgs, E_i jest proporcjonalne do √|K(d)|.
    for i, d in enumerate(OCTAVES_EFFECTIVE):
        # Używamy K(d) jako bazowej miary energii/masy, co jest ekwiwalentem
        # masy renormalizowanej w teorii pola.
        H[i, i] = np.abs(kernel_K(d))
    
    # 2. Elementy pozadiagonalne (Sprzężenia Międzyoktawowe)
    for i in range(n_modes):
        for j in range(i + 1, n_modes):
            dist = abs(OCTAVES_EFFECTIVE[i] - OCTAVES_EFFECTIVE[j])
            K_ij = kernel_K(dist)
            # Sprzężenie jest symetryczne i wnosi wkład do masy (elementy nidiagonalne)
            H[i, j] = H[j, i] = K_ij * 0.5 # Skalowanie ze względu na symetrię

    return H

def diagonalize_and_analyze(H: np.ndarray, octave_labels: list):
    """Diagonalizuje macierz H i analizuje wartości własne (masy)."""
    
    eigenvalues, _ = np.linalg.eigh(H)
    
    # Wybieramy 3 najmniejsze dodatnie wartości własne jako kandydatów na masy leptonów (e, μ, τ).
    # W teorii powinny to być 3 najmniejsze pozycje, które odpowiadają fermionom.
    
    # Odrzucamy ujemne (niestabilność) lub bardzo małe (zero-modes)
    positive_eigenvalues = sorted([e for e in eigenvalues if e > 1e-8])
    
    if len(positive_eigenvalues) < 3:
        print("❌ Błąd: Nie znaleziono co najmniej 3 stabilnych, dodatnich modów.")
        return None
        
    # Przyjmujemy 3 najlżejsze mody
    m_e_pred = positive_eigenvalues[0]
    m_mu_pred = positive_eigenvalues[1]
    m_tau_pred = positive_eigenvalues[2]
    
    print("\n--- WYNIKI ANALIZY SPEKTRALNEJ ---")
    print(f"Wartości własne (proporcjonalne do masy):")
    for i, val in enumerate(positive_eigenvalues[:6]):
        print(f"  Mod {i+1}: {val:.6f}")
        
    print("\n--- PORÓWNANIE HIERARCHII LEPTODÓW (NA PODSTAWIE K(d)) ---")
    
    print(f"SM Ratios: m_μ/m_e = {SM_MASSES['muon']/SM_MASSES['electron']:.2f} | m_τ/m_μ = {SM_MASSES['tau']/SM_MASSES['muon']:.2f}")
    
    # Predykcja bez dodatkowych wzmocnień A_i (czysty wynik K(d))
    if m_e_pred > 1e-9:
        pred_ratio_mu_e = m_mu_pred / m_e_pred
        pred_ratio_tau_mu = m_tau_pred / m_mu_pred
        
        print(f"PRED. Ratios (K(d) ONLY): m_μ/m_e = {pred_ratio_mu_e:.3f} | m_τ/m_μ = {pred_ratio_tau_mu:.3f}")
        
        print("\n--- WNIOSEK Z CZYSTEGO TESTU ---")
        print("Hierarchia jest nadal zdominowana przez błąd strukturalny.")
        if pred_ratio_mu_e < 10.0:
            print("❌ Różnica rzędu 100x (SM: 207) wskazuje, że K(d) SAMODZIELNIE NIE GENERUJE HIERARCHII.")
        else:
            print("✅ Hierarchia jest bliższa SM niż w poprzednich testach (sprawdzić dlaczego).")
            
    return m_e_pred, m_mu_pred, m_tau_pred

# ============================================================================
# MAIN
# ============================================================================

def main():
    print("======================================================================")
    print(" TEST MINIMALNY: HIERARCHIA MAS TYLKO Z JĄDRA K(d)")
    print("======================================================================")
    
    n = len(OCTAVES_EFFECTIVE)
    H_matrix = construct_static_hamiltonian_matrix(n)
    
    print(f"✓ Skonstruowano macierz {n}x{n} H_Lepton.")
    
    # Przykładowe wartości na diagonali (aby zorientować się w skali)
    print(f"Przykładowa diagonala H (Energia/Masa bazowa): {H_matrix[0,0]:.4f} (d=1) -> {H_matrix[-1,-1]:.4f} (d=12)")
    
    diagonalize_and_analyze(H_matrix, OCTAVES_EFFECTIVE)
    
    print("\n======================================================================")
    print(" TEST ZAKOŃCZONY. Wynik analizy czystego K(d) jest miarą potencjału teorii.")
    print("======================================================================")

if __name__ == "__main__":
    main()
