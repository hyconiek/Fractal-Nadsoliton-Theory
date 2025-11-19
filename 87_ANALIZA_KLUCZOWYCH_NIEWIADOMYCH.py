# -*- coding: utf-8 -*-
# Author: Krzysztof Żuchowski

"""
BADANIE 87: ANALIZA KLUCZOWYCH NIEWIADOMYCH I TESTOWANIE NOWYCH PARADYGMATÓW

Data utworzenia: 14 listopada 2025
Autor: GitHub Copilot (na podstawie analizy SWOT i dyskusji)

CEL:
Ten skrypt stanowi framework badawczy do systematycznego badania kluczowych
problemów zidentyfikowanych w teorii fraktalnego nadsolitona. Zamiast próbować
"naprawiać" istniejące mechanizmy, skupia się na testowaniu nowych, fundamentalnych
hipotez (zaproponowanych w QW-V87 do QW-V91), które mogą rozwiązać problemy u ich
źródeł.

METODOLOGIA:
1.  Ustanowienie solidnej bazy: Implementacja 4-parametrowego rdzenia teorii.
2.  Zdefiniowanie funkcji dla każdego z 4 kluczowych problemów (niewiadomych).
3.  Implementacja "przełączników" do testowania starych vs. nowych hipotez.
4.  Skupienie na analizie strukturalnej i jakościowej, a nie na precyzyjnym dopasowaniu.
"""

import numpy as np
from scipy.linalg import eig
import os
import markdown
import pandas as pd
import sys
from io import StringIO
from datetime import datetime

# =============================================================================
# CZĘŚĆ 1: RDZEŃ TEORETYCZNY (ODKRYCIA Z QW-V46-V50)
# =============================================================================
# Tutaj definiujemy fundamenty, które okazały się solidne i predykcyjne.

# 4 fundamentalne parametry minimalne
ALPHA_GEO = 1.0
BETA_TORS = 0.1
OMEGA = 0.7854  # rad
PHI = 0.5236    # rad

# 8 efektywnych oktaw (odkrycie z QW-V40)
EFFECTIVE_OCTAVES = np.array([1, 3, 4, 6, 7, 9, 10, 12])

def get_coupling_kernel(d, complex_kernel=False):
    """
    Uniwersalne jądro sprzężeń K(d).
    Argument `complex_kernel` pozwala na testowanie hipotezy z QW-V88.
    """
    if complex_kernel:
        return ALPHA_GEO * np.exp(1j * (OMEGA * d + PHI)) / (1 + BETA_TORS * d)
    else:
        return ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + BETA_TORS * d)

def get_self_coupling_matrix(complex_kernel=False):
    """
    Generuje macierz samosprzężeń S_ij dla 8 efektywnych oktaw.
    """
    n = len(EFFECTIVE_OCTAVES)
    s_matrix = np.zeros((n, n), dtype=np.complex128 if complex_kernel else np.float64)
    for i in range(n):
        for j in range(n):
            d = abs(EFFECTIVE_OCTAVES[i] - EFFECTIVE_OCTAVES[j])
            if d > 0:
                s_matrix[i, j] = get_coupling_kernel(d, complex_kernel)
    return s_matrix

# Obliczamy bazową (rzeczywistą) macierz S_ij
S_ij_real = get_self_coupling_matrix(complex_kernel=False)

print("--- RDZEŃ TEORETYCZNY ZAINICJOWANY ---")
print(f"4 parametry fundamentalne: alpha_geo={ALPHA_GEO}, beta_tors={BETA_TORS}, omega={OMEGA}, phi={PHI}")
print(f"Macierz samosprzężeń S_ij (8x8) została obliczona.")
print("-" * 40, "\n")


# =============================================================================
# CZĘŚĆ 2: ANALIZA NIEWIADOMYCH
# =============================================================================

# --- NIEWIADOMA 1: HIERARCHIA MAS LEPTONÓW (Problem symetrii) ---
# Hipoteza (QW-V87): Hierarchia wynika z anharmonicznych przesunięć w geometrii
# oktaw, które łamią idealną symetrię.

def analyze_mass_hierarchy(use_anharmonic_shift=False, epsilon=1/137.0):
    """
    Analizuje siłę rezonansową dla każdej oktawy.
    - `use_anharmonic_shift=False`: Odtwarza problem symetrii (QW-V72, QW-V78).
    - `use_anharmonic_shift=True`: Testuje hipotezę złamania symetrii (QW-V87).
    """
    print("--- Analiza 1: Hierarchia Mas Leptonów ---")
    
    # Metoda bazowa: Suma sił sprzężeń (pokazuje problem symetrii)
    symmetric_strength = np.sum(np.abs(S_ij_real), axis=1)
    
    if not use_anharmonic_shift:
        print("Metoda bazowa (symetryczna): Siła sprzężenia dla każdej z 8 oktaw.")
        print(np.round(symmetric_strength, 4))
        print("Wniosek: Siły są bardzo podobne, co uniemożliwia generację hierarchii O(100).")
        return symmetric_strength
    else:
        print(f"Test hipotezy (anharmonicznej) z epsilon = {epsilon:.4f}")
        n = len(EFFECTIVE_OCTAVES)
        anharmonic_strength = np.zeros(n)
        for i in range(n):
            total_strength = 0
            for j in range(n):
                if i == j: continue
                
                # Anharmoniczne przesunięcie odległości efektywnej
                d_eff = abs(EFFECTIVE_OCTAVES[i] - EFFECTIVE_OCTAVES[j]) * (1 + epsilon * S_ij_real[i, j])
                
                # Siła sprzężenia zależy od tej nowej, złamanej odległości
                total_strength += np.abs(get_coupling_kernel(d_eff))
            anharmonic_strength[i] = total_strength
            
        print("Metoda anharmoniczna: Siła sprzężenia dla każdej z 8 oktaw.")
        print(np.round(anharmonic_strength, 4))
        
        ratio = np.max(anharmonic_strength) / np.min(anharmonic_strength)
        print(f"Stosunek max/min siły: {ratio:.4f}")
        print("Wniosek: Anharmoniczność łamie symetrię. Należy zbadać, czy ten efekt jest wystarczająco silny.")
        return anharmonic_strength

# --- NIEWIADOMA 2: DYNAMIKA FAZOWA I KĄTY CKM (Problem zespolony) ---
# Hipoteza (QW-V88): Kąty mieszania i rodziny cząstek wynikają z dynamicznej
# interferencji w zespolonej sieci przepływu informacji.

def analyze_complex_information_flow():
    """
    Bada strukturę przepływu informacji w zespolonej macierzy S_ij.
    """
    print("\n--- Analiza 2: Dynamika Fazowa i Kąty CKM ---")
    s_matrix_complex = get_self_coupling_matrix(complex_kernel=True)
    
    # Szukamy ścieżek interferencji (np. z oktawy i do j przez k)
    # Im(S_ik * S_kj) != 0 oznacza, że faza zależy od ścieżki
    interference_matrix = np.zeros_like(s_matrix_complex, dtype=float)
    n = s_matrix_complex.shape[0]
    for i in range(n):
        for j in range(n):
            # Sumujemy wkłady od wszystkich możliwych pośredników k
            interference_term = np.sum([np.imag(s_matrix_complex[i, k] * s_matrix_complex[k, j]) 
                                        for k in range(n) if k != i and k != j])
            interference_matrix[i, j] = interference_term
            
    print("Macierz interferencji Im(S_ik * S_kj):")
    print(np.round(interference_matrix, 4))
    print("Wniosek: Niezerowe wartości wskazują na złożoną dynamikę fazową.")
    print("Należy zbadać, czy wzorce w tej macierzy grupują oktawy w rodziny.")
    return interference_matrix

def cluster_octaves_by_interference(interference_matrix):
    """
    Używa klasteryzacji (K-Means), aby sprawdzić, czy oktawy grupują się w rodziny.
    """
    try:
        from sklearn.cluster import KMeans
        print("\n--- Analiza 2a: Klasteryzacja Oktaw (Hipoteza Rodzin) ---")
        
        # Używamy wartości bezwględnych interferencji jako miary "podobieństwa"
        kmeans = KMeans(n_clusters=3, random_state=42, n_init=10)
        clusters = kmeans.fit_predict(np.abs(interference_matrix))
        
        print("Wyniki klasteryzacji (przypisanie oktawy do rodziny):")
        for i, octave in enumerate(EFFECTIVE_OCTAVES):
            print(f"  Oktawa {octave:2d}: Rodzina {clusters[i]}")
            
        print("\nWniosek: Wynik pokazuje, czy struktura interferencji w naturalny sposób")
        print("dzieli oktawy na 3 grupy, co mogłoby odpowiadać trzem generacjom fermionów.")

    except ImportError:
        print("\n--- Analiza 2a: Klasteryzacja Oktaw (Pominięto) ---")
        print("Proszę zainstalować scikit-learn, aby uruchomić tę analizę: pip install scikit-learn")

# --- NIEWIADOMA 3: EMERGENTNA GRAWITACJA (Problem mechanizmu) ---
# Hipoteza (QW-V90): Grawitacja wyłania się nie z gęstości pola, ale z jego
# dynamiki, jako metryka akustyczna w "płynie informacyjnym".

def analyze_acoustic_metric_gravity():
    """
    Szkic koncepcyjny dla nowego paradygmatu grawitacji.
    """
    print("\n--- Analiza 3: Emergentna Grawitacja (Metryka Akustyczna) ---")
    print("Problem: Mapowanie g_μν(ρ) zawiodło (korelacja G~T ≈ 0 lub ujemna).")
    print("Hipoteza: Grawitacja jako metryka akustyczna w płynie informacyjnym.")
    
    # To jest funkcja koncepcyjna, ponieważ u_mu i c_s nie są jeszcze zdefiniowane
    def conceptual_acoustic_metric(u_mu, c_s):
        eta_munu = np.diag([-1, 1, 1, 1])
        # g_μν_eff ∝ η_μν + (1 - c_s²) u_μ u_ν
        # Ta formuła pokazuje, jak fluktuacje w płynie (zmiany c_s) mogą generować krzywiznę.
        pass
        
    print("Koncepcja: g_μν_eff ∝ η_μν + (1 - c_s²) u_μ u_ν")
    print("Wniosek: To podejście oddziela metrykę od T_μν, unikając tautologii.")
    print("Kolejny krok: Zdefiniować czteroprędkość (u_μ) i prędkość dźwięku (c_s) dla pola nadsolitona.")

# --- NIEWIADOMA 4: SYMETRIE GAUGE (Problem emergencji) ---
# Hipoteza (QW-V91): Symetrie cechowania nie wynikają z arbitralnego mapowania
# odległości 'd', ale z wewnętrznej algebry sieci rezonansowej.

def analyze_emergent_symmetries():
    """
    Bada wektory własne macierzy rezonansowej w poszukiwaniu struktur algebraicznych.
    """
    print("\n--- Analiza 4: Emergencja Symetrii Gauge ---")
    print("Problem: Mapowanie d=1,2,3 na SU(3),SU(2),U(1) jest ansatzem i zawodzi ilościowo.")
    print("Hipoteza: Symetrie są emergentną właściwością sieci rezonansowej.")
    
    # Krok 1: Zbuduj macierz rezonansową M_ij (uproszczona miara)
    # M_ij = siła rezonansu między oktawami i oraz j.
    # Użyjemy |S_ij| jako prostej miary.
    resonance_matrix = np.abs(S_ij_real)
    
    # Krok 2: Znajdź wektory własne
    eigenvalues, eigenvectors = eig(resonance_matrix)
    
    print("Wartości własne macierzy rezonansowej:")
    print(np.round(np.real(eigenvalues), 4))
    
    # Wektory własne tworzą bazę, w której operuje dynamika.
    # Należy zbadać ich komutatory, aby znaleźć generatory algebry.
    print("Wniosek: Wektory własne definiują 'naturalne' mody oscylacji systemu.")
    print("Kolejny krok: Zbudować generatory z wektorów własnych i zbadać ich algebrę komutacyjną.")
    return eigenvalues, eigenvectors


# =============================================================================
# CZĘŚĆ 3: FUNKCJE DO ZAPISYWANIA WYNIKÓW
# =============================================================================

def save_results_to_md(results_dict, filename='results_87_analysis.md'):
    """
    Zapisuje wyniki analiz do pliku Markdown z formatowaniem.
    """
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(f"# Analiza Kluczowych Niewiadomych - Raport\n")
        f.write(f"**Data wygenerowania:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"**Plik:** `87_ANALIZA_KLUCZOWYCH_NIEWIADOMYCH.py`\n\n")
        f.write("---\n\n")
        
        for section, content in results_dict.items():
            f.write(f"## {section}\n\n")
            f.write(f"{content}\n\n")
            f.write("---\n\n")

def capture_function_output(func):
    """
    Przechwytuje wyjście funkcji (print statements) i zwraca je jako string.
    """
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        func()
        output = sys.stdout.getvalue()
    finally:
        sys.stdout = old_stdout
    return output

# =============================================================================
# CZĘŚĆ 3: URUCHOMIENIE ANALIZ
# =============================================================================

if __name__ == "__main__":
    results = {}
    
    print("Uruchamianie analiz i zapisywanie wyników...\n")
    
    # Analiza 1: Hierarchia Mas
    print("=== Analiza 1: Hierarchia Mas Leptonów ===")
    output1_base = capture_function_output(lambda: analyze_mass_hierarchy(use_anharmonic_shift=False))
    output1_new = capture_function_output(lambda: analyze_mass_hierarchy(use_anharmonic_shift=True))
    results['Analiza 1: Hierarchia Mas Leptonów'] = f"{output1_base}\n{output1_new}"
    print(output1_base)
    print("\n--- Testowanie nowej hipotezy ---")
    print(output1_new)
    
    # Analiza 2: Dynamika CKM
    print("\n=== Analiza 2: Dynamika Fazowa i Kąty CKM ===")
    output2 = capture_function_output(lambda: analyze_complex_information_flow())
    results['Analiza 2: Dynamika Fazowa'] = output2
    print(output2)
    
    # Analiza 2a: Klasteryzacja
    print("=== Analiza 2a: Klasteryzacja Oktaw ===")
    interference_matrix = analyze_complex_information_flow()
    output2a = capture_function_output(lambda: cluster_octaves_by_interference(interference_matrix))
    results['Analiza 2a: Klasteryzacja Oktaw'] = output2a
    print(output2a)
    
    # Analiza 3: Grawitacja
    print("\n=== Analiza 3: Emergentna Grawitacja ===")
    output3 = capture_function_output(lambda: analyze_acoustic_metric_gravity())
    results['Analiza 3: Emergentna Grawitacja'] = output3
    print(output3)
    
    # Analiza 4: Symetrie Gauge
    print("\n=== Analiza 4: Emergencja Symetrii Gauge ===")
    output4 = capture_function_output(lambda: analyze_emergent_symmetries())
    results['Analiza 4: Emergencja Symetrii Gauge'] = output4
    print(output4)

    print("\n--- ZAKOŃCZONO ANALIZĘ KLUCZOWYCH NIEWIADOMYCH ---")
    print("Zapisywanie wyników do pliku Markdown...\n")
    
    # Zapisz wyniki do pliku
    save_results_to_md(results, filename='results_87_analysis.md')
    print(f"✅ Wyniki zapisane do pliku: results_87_analysis.md")
