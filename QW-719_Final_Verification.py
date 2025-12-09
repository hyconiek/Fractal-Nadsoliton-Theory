#!/usr/bin/env python3
"""
QW-719: FINAL MASS THEORY VERIFICATION
======================================
Cel: Ostateczna weryfikacja spójności Modelu Masy.
     Sprawdzamy, czy stany topologiczne (W=1, 2, 3) mają poprawne energie
     w PEŁNYM Hamiltonianie ZTP (a nie tylko w modelu geometrycznym).

Metoda:
1. Konstruujemy pełny Hamiltonian ZTP (12 oktaw, K(d), człony nieliniowe).
2. Konstruujemy stany kwantowe odpowiadające topologiom W=1, 2, 3.
   (Uwaga: Hamiltonian jest 1D w przestrzeni oktaw, więc musimy zrzutować
    efektywną energię wiru 3D na potencjał w przestrzeni oktaw).
3. Obliczamy energię własną tych stanów.
4. Stosujemy formułę 4-bitową: M = M_e * 2^(4 * E_rel).

Oczekiwanie:
  E(W=2)/E(W=1) ≈ 1.86 (żeby dać masę mionu ~200x)
  E(W=3)/E(W=1) ≈ 2.79 (żeby dać masę tau/protonu ~2000x)

Date: 2025-12-08
"""

import numpy as np
from scipy.linalg import eigh
import datetime

print("="*80)
print("QW-719: OSTATECZNA WERYFIKACJA MODELU MASY")
print("="*80)

# ===========================================================================
# 1. PEŁNY HAMILTONIAN ZTP (z QW-701)
# ===========================================================================

N_OCTAVES = 12
ALPHA = 4 * np.log(2)  # 2.7726 (Kluczowa stała!)
BETA = 0.01

def K_kernel(d):
    if d == 0: return ALPHA
    return ALPHA * np.cos(np.pi/4 * d + np.pi/6) / (1 + BETA * d)

# Hamiltonian macierzowy
H_ZTP = np.zeros((N_OCTAVES, N_OCTAVES))
for i in range(N_OCTAVES):
    for j in range(N_OCTAVES):
        H_ZTP[i, j] = -K_kernel(abs(i - j))

# Diagonalna poprawka (Self-energy rośnie z oktawą - "potencjał studni")
# To symuluje fakt, że wyższe oktawy są bardziej energetyczne
for i in range(N_OCTAVES):
    H_ZTP[i, i] += 0.1 * (i**2)  # Harmonic confinement (standard w ZTP)

print("Hamiltonian ZTP zbudowany.")

# ===========================================================================
# 2. EFEKTYWNY WKŁAD OD TOPOLOGII (WINDING NUMBER)
# ===========================================================================

print("""
JAK TOPOLOGIA W=1, 2, 3 WPŁYWA NA HAMILTONIAN 1D?

Mamy hipotezę z QW-718:
Energia kinetyczna wiru E_kin skaluje się jak geometria (W).
W przestrzeni oktaw przejawia się to jako DODATKOWY POTENCJAŁ ODŚRODKOWY.

V_topo(W, n) = C * W² / R(n)²

Gdzie R(n) to "promień" w przestrzeni oktaw (R ~ n + 1).
""")

def get_hamiltonian_with_topology(W, coupling_strength=0.05):
    H = H_ZTP.copy()
    # Dodajemy potencjał topologiczny (odśrodkowy)
    for n in range(N_OCTAVES):
        # V_centrifugal ~ W^2 / r^2
        # r w przestrzeni oktaw to (n+1)
        V_topo = coupling_strength * (W**2) / ((n + 1)**2)
        H[n, n] += V_topo
    return H

# ===========================================================================
# 3. OBLICZANIE ENERGII STANU PODSTAWOWEGO DLA W=1, 2, 3
# ===========================================================================

print("\nObliczanie stanów podstawowych dla różnych topologii...")
print(f"{'W':<5} {'Energia E0':>12} {'Delta E':>12} {'Ratio E/E(1)':>15}")
print("-"*70)

results = {}
base_energy = 0

# Parametr coupling_strength trzeba skalibrować, żeby W=1 pasowało do fizyki
# (To jedyny wolny parametr - siła sprzężenia topologii z oktawami)
# Ale sprawdźmy czy istnieje "naturalne" C, np. C=ALPHA?

C_topo = ALPHA # Testujemy hipotezę "wszystko jest ALPHA"

for W in [1, 2, 3]:
    H_W = get_hamiltonian_with_topology(W, coupling_strength=C_topo)
    evals, evecs = eigh(H_W)
    E0 = evals[0] # Energia stanu podstawowego
    
    if W == 1:
        base_energy = E0
    
    # Przesuwamy energię tak, by E(W=1) była naszą bazą (lub zerem, zależy od modelu)
    # W modelu 4-bitowym liczy się E_relatywne.
    # Użyjmy E_rel = E0(W) / E0(1) ? Nie, energia jest ujemna w wiązaniu.
    # Użyjmy E_binding = |E0|
    
    E_bind = abs(E0)
    base_bind = abs(base_energy)
    
    results[W] = E_bind
    
    ratio = E_bind / base_bind if base_bind > 0 else 0
    
    print(f"{W:<5} {E0:>12.4f} {abs(E0-base_energy):>12.4f} {ratio:>15.4f}")

# ===========================================================================
# 4. PREDYKCJA MASY (FORMUŁA 4-BITOWA)
# ===========================================================================

print("\n" + "="*80)
print("PREDYKCJA MASY Z FORMUŁY 4-BITOWEJ (H21)")
print("M = M_e * 2^(4 * E_skalowne)")
print("="*80)

# Jak zmapować E_binding na "E_skalowalne"?
# Z QW-718 wiemy, że E_kin(W)/E_kin(1) dawało alpha.
# Zobaczmy stosunki tutaj.

E1 = results[1]
E2 = results[2]
E3 = results[3]

# Obliczamy "efektywny n_skok" analogicznie do QW-718
# n_skok = (E(W) - E(1)) * k? Nie.
# Spróbujmy po prostu stosunku Binding Energy.

scaling_factor = 4.0 # Bo 4 bity

print(f"Stosunek E(2)/E(1) = {E2/E1:.4f}")
print(f"Stosunek E(3)/E(1) = {E3/E1:.4f}")

# Hipoteza: E_scaling = (E(W)/E(1)) * ALPHA ?
# Nie, sprawdźmy co pasuje.

# Mion: M = 207 * M_e -> log2(207) = 7.7 bitów
# Tau: M = 3477 * M_e -> log2(3477) = 11.8 bitów
# Proton: M = 1836 * M_e -> log2(1836) = 10.8 bitów

# Mamy 4 bity bazowe.
# n_bits_excess(Mion) = 7.7 - 0 = 7.7
# n_bits_excess(Tau) = 11.8 - 0 = 11.8

# Zobaczmy czy (E(W) - E(1)) pasuje do tych bitów.

dE2 = E1 - E2 # (Energia wiązania mniejsza czy większa? E1 jest najgłębsza?)
# W Hamiltonianie dodaliśmy V_topo > 0.
# Więc E0 idzie W GÓRĘ (mniej ujemne, słabsze wiązanie) dla wyższych W.
# E_ground = -Deep + Positive_Topo
# E(1) = -10 (np)
# E(2) = -8 (wyżej)
# Diff = 2

diff_21 = abs(results[2] - results[1])
diff_31 = abs(results[3] - results[1])

print(f"Różnica energii ΔE(2-1) = {diff_21:.4f}")
print(f"Różnica energii ΔE(3-1) = {diff_31:.4f}")

# Sprawdźmy stałą konwersji: k * ΔE = bity
# Dla mionu: k * diff_21 = 7.7
# k = 7.7 / diff_21

k_calib = 7.7 / diff_21
print(f"\nKalibracja na mionie: k = {k_calib:.4f} bitów/jednostkę energii")

# Predykcja dla Tau/Protonu (W=3)
bits_3 = k_calib * diff_31
mass_3_pred = 2**bits_3

print(f"\nPREDYKCJA DLA W=3 (Używając kalibracji z W=2):")
print(f"Liczba bitów: {bits_3:.2f}")
print(f"Przewidywany stosunek mas: {mass_3_pred:.2f}")
print(f"Cel (Tau): 3477")
print(f"Cel (Proton): 1836")

error_tau = abs(mass_3_pred - 3477)/3477 * 100
error_proton = abs(mass_3_pred - 1836)/1836 * 100

print(f"Błąd (Tau): {error_tau:.1f}%")
print(f"Błąd (Proton): {error_proton:.1f}%")

# Save report
with open("raport_qw719_final_verification.md", "w") as f:
    f.write("# RAPORT QW-719: Ostateczna Weryfikacja Modelu Masy\\n")
    f.write(f"**Date:** {datetime.datetime.now()}\\n\\n")
    
    f.write("## Wyniki Energii ZTP z Topologią\\n")
    f.write("| W | Energia E0 | ΔE od W=1 |\\n")
    f.write("|---|---|---|\\n")
    for W in [1, 2, 3]:
        diff = abs(results[W] - results[1])
        f.write(f"| {W} | {results[W]:.4f} | {diff:.4f} |\\n")
        
    f.write(f"\\n## Kalibracja i Predykcja\\n")
    f.write(f"Stała sprzężenia energii z bitami (z mionu): k = {k_calib:.4f} bit/E\\n")
    f.write(f"Predykcja bitów dla W=3: {bits_3:.2f} bitów\\n")
    f.write(f"Predykcja masy dla W=3: {mass_3_pred:.2f} M_e\\n")
    f.write(f"Błąd vs Tau: {error_tau:.1f}%\\n")
    f.write(f"Błąd vs Proton: {error_proton:.1f}%\\n")

print("\\nReport saved to raport_qw719_final_verification.md")
print("="*80)
