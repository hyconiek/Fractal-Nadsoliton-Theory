#!/usr/bin/env python3
"""
QW-720: MAPOWANIE OKTAW → PRZESTRZEŃ DLA WIDMA ATOMOWEGO
========================================================
Purpose: Wyprowadzić związek między odległością oktawową d a odległością przestrzenną r
         dla poprawnego modelowania widma atomowego.

Problem:
  - Oktawy: d = 1-12 (dyskretne, bez jednostek)
  - Przestrzeń: r ~ 0.5-10 Å (ciągłe, fizyczne jednostki)
  - QW-699/700 używały d bezpośrednio jako r → błąd 50-110%

Hipotezy do przetestowania:
  1. r ∝ d^α (potęgowe)
  2. r ∝ exp(d/α) (eksponencjalne)
  3. r ∝ 1/|K(d)| (z teorii: "r_ij ∝ 1/|K(d_ij)|")
  4. r ∝ d × a_B (Bohr radius jako skala)

Metoda:
  1. Wyprowadzić mapowanie z wymiarów fizycznych
  2. Kalibrować do znanych skali atomowych
  3. Zweryfikować na widmie wodoru
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.linalg import eigh
import datetime

print("="*80)
print("QW-720: MAPOWANIE OKTAW → PRZESTRZEŃ")
print("="*80)

# --- Constants ---
ALPHA = 4 * np.log(2)  # 2.7726
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA = 0.1

# Physical constants
A_BOHR = 0.529177  # Å (Bohr radius)
R_BOHR = A_BOHR * 1e-10  # m
ALPHA_EM = 1/137.036  # Fine structure constant
M_E = 9.109e-31  # kg (electron mass)
E_RYDBERG = 13.6  # eV

# Stable orbits (from theory)
D1 = 1.3333  # Electron
D2 = 9.3333  # Muon
D3 = 17.3333  # Tau

def K(d):
    """K(d) kernel"""
    if d < 0.1:
        d = 0.1
    return ALPHA * np.cos(OMEGA * d + PHI) / (1 + BETA * d)

# --- HYPOTHESIS 1: r ∝ d^α ---
print("\n[1] HYPOTEZA 1: r ∝ d^α (potęgowe)")
print("-" * 60)

# Calibration: Electron at d=1.33 should be at r ≈ a_B
# r(d) = r_0 × d^α
# For d=1.33, r ≈ a_B = 0.529 Å
r_0_power = A_BOHR / (D1 ** ALPHA)

def r_power(d):
    return r_0_power * (d ** ALPHA)

print(f"  Kalibracja: r(d=1.33) = a_B = {A_BOHR:.4f} Å")
print(f"  r_0 = {r_0_power:.6f} Å")
print(f"  r(d=1.33) = {r_power(D1):.4f} Å")
print(f"  r(d=9.33) = {r_power(D2):.4f} Å")
print(f"  r(d=17.33) = {r_power(D3):.4f} Å")

# --- HYPOTHESIS 2: r ∝ exp(d/α) ---
print("\n[2] HYPOTEZA 2: r ∝ exp(d/α) (eksponencjalne)")
print("-" * 60)

# r(d) = r_0 × exp(d/α)
# For d=1.33, r ≈ a_B
r_0_exp = A_BOHR / np.exp(D1 / ALPHA)

def r_exp(d):
    return r_0_exp * np.exp(d / ALPHA)

print(f"  r_0 = {r_0_exp:.6f} Å")
print(f"  r(d=1.33) = {r_exp(D1):.4f} Å")
print(f"  r(d=9.33) = {r_exp(D2):.4f} Å")
print(f"  r(d=17.33) = {r_exp(D3):.4f} Å")

# --- HYPOTHESIS 3: r ∝ 1/|K(d)| (z teorii) ---
print("\n[3] HYPOTEZA 3: r ∝ 1/|K(d)| (z teorii)")
print("-" * 60)

# From FIN_Theory_Paper.tex: "r_ij ∝ 1/|K(d_ij)|"
# Calibration: r(d=1.33) = a_B
K_d1 = abs(K(D1))
r_0_k = A_BOHR * K_d1

def r_k(d):
    K_d = abs(K(d))
    if K_d < 1e-10:
        K_d = 1e-10  # Avoid division by zero
    return r_0_k / K_d

print(f"  K(d=1.33) = {K_d1:.6f}")
print(f"  r_0 = {r_0_k:.6f} Å")
print(f"  r(d=1.33) = {r_k(D1):.4f} Å")
print(f"  r(d=9.33) = {r_k(D2):.4f} Å")
print(f"  r(d=17.33) = {r_k(D3):.4f} Å")

# --- HYPOTHESIS 4: r ∝ d × a_B (liniowe) ---
print("\n[4] HYPOTEZA 4: r ∝ d × a_B (liniowe)")
print("-" * 60)

# r(d) = d × a_B / d_0
# For d=1.33, r ≈ a_B
d_0_linear = D1

def r_linear(d):
    return d * A_BOHR / d_0_linear

print(f"  r(d=1.33) = {r_linear(D1):.4f} Å")
print(f"  r(d=9.33) = {r_linear(D2):.4f} Å")
print(f"  r(d=17.33) = {r_linear(D3):.4f} Å")

# --- PHYSICAL CONSTRAINTS ---
print("\n[5] OGRANICZENIA FIZYCZNE")
print("-" * 60)

# Hydrogen orbital radii (Bohr model)
# r_n = n² × a_B
r_hydrogen = {
    1: 1 * A_BOHR,      # n=1: a_B
    2: 4 * A_BOHR,      # n=2: 4a_B
    3: 9 * A_BOHR,      # n=3: 9a_B
    4: 16 * A_BOHR,     # n=4: 16a_B
}

print("  Promienie orbitali wodoru (Bohr):")
for n, r_n in r_hydrogen.items():
    print(f"    n={n}: r = {r_n:.4f} Å")

# Map octaves to quantum numbers
# Hypothesis: d_n = d_1 + 8(n-1) for n=1,2,3,4
d_hydrogen = {
    1: D1,              # n=1: d=1.33
    2: D1 + 8,          # n=2: d=9.33
    3: D1 + 16,         # n=3: d=17.33
    4: D1 + 24,         # n=4: d=25.33
}

print("\n  Mapowanie oktaw → liczby kwantowe:")
for n, d_n in d_hydrogen.items():
    print(f"    n={n}: d = {d_n:.2f}")

# Test each hypothesis
print("\n[6] TEST MAPOWAŃ")
print("-" * 60)

hypotheses = {
    "Potęgowe (d^α)": r_power,
    "Eksponencjalne (exp(d/α))": r_exp,
    "Odwrotność K(d)": r_k,
    "Liniowe (d×a_B)": r_linear,
}

results = {}

for name, func in hypotheses.items():
    errors = []
    print(f"\n  {name}:")
    for n in [1, 2, 3, 4]:
        d_n = d_hydrogen[n]
        r_pred = func(d_n)
        r_exp = r_hydrogen[n]
        err = abs(r_pred - r_exp) / r_exp * 100
        errors.append(err)
        print(f"    n={n}: r_pred={r_pred:.4f} Å, r_exp={r_exp:.4f} Å, błąd={err:.1f}%")
    
    mean_err = np.mean(errors)
    results[name] = mean_err
    print(f"    Średni błąd: {mean_err:.1f}%")

# Find best hypothesis
best_hyp = min(results, key=results.get)
best_err = results[best_hyp]

print(f"\n[7] NAJLEPSZA HYPOTEZA: {best_hyp}")
print(f"    Średni błąd: {best_err:.1f}%")

# --- DERIVATION FROM FIRST PRINCIPLES ---
print("\n[8] DERIWACJA Z PIERWSZYCH ZASAD")
print("-" * 60)

print("""
Z teorii FIN:
  1. Odległość przestrzenna wynika z korelacji informacyjnej
  2. r_ij ∝ 1/|K(d_ij)| (z FIN_Theory_Paper.tex)
  3. K(d) = α·cos(ωd+φ)/(1+βd)
  
Dla małych d: K(d) ≈ α·cos(φ) ≈ α·0.866 ≈ 2.4
Dla dużych d: K(d) ≈ α·cos(ωd+φ)/(βd) ~ 1/d

Więc: r(d) ∝ 1/|K(d)| ≈ d/α dla dużych d
Ale dla małych d: r(d) ≈ const

Lepsze przybliżenie: r(d) = r_0 / |K(d)|
""")

# Refined hypothesis: r(d) = r_0 / |K(d)| with calibration
K_d1_refined = abs(K(D1))
r_0_refined = A_BOHR * K_d1_refined

def r_refined(d):
    K_d = abs(K(d))
    if K_d < 1e-10:
        K_d = 1e-10
    return r_0_refined / K_d

print(f"  Refined: r(d) = {r_0_refined:.6f} / |K(d)|")
print(f"  r(d=1.33) = {r_refined(D1):.4f} Å")
print(f"  r(d=9.33) = {r_refined(D2):.4f} Å")
print(f"  r(d=17.33) = {r_refined(D3):.4f} Å")

# Test refined
errors_refined = []
for n in [1, 2, 3, 4]:
    d_n = d_hydrogen[n]
    r_pred = r_refined(d_n)
    r_exp = r_hydrogen[n]
    err = abs(r_pred - r_exp) / r_exp * 100
    errors_refined.append(err)

mean_err_refined = np.mean(errors_refined)
print(f"  Średni błąd (refined): {mean_err_refined:.1f}%")

# --- FINAL MAPPING FUNCTION ---
print("\n[9] FINALNA FUNKCJA MAPOWANIA")
print("-" * 60)

# Use the best hypothesis
if mean_err_refined < best_err:
    mapping_func = r_refined
    mapping_name = "r(d) = r_0 / |K(d)| (refined)"
    final_err = mean_err_refined
else:
    mapping_func = hypotheses[best_hyp]
    mapping_name = best_hyp
    final_err = best_err

print(f"  Wybrana funkcja: {mapping_name}")
print(f"  Średni błąd promieni: {final_err:.1f}%")

# Save mapping function for use in QW-721
def octave_to_space(d):
    """Map octave distance d to spatial distance r in Å"""
    return mapping_func(d)

# --- REPORT ---
report_file = "raport_qw720_octave_space_mapping.md"
print(f"\n[10] ZAPISYWANIE RAPORTU: {report_file}")

with open(report_file, "w") as f:
    f.write("# RAPORT QW-720: MAPOWANIE OKTAW → PRZESTRZEŃ\n")
    f.write(f"**Data:** {datetime.datetime.now()}\n\n")
    
    f.write("## 1. Problem\n")
    f.write("Oktawy (d=1-12) muszą być zmapowane na przestrzeń (r ~ Å) dla widma atomowego.\n\n")
    
    f.write("## 2. Testowane Hipotezy\n")
    f.write("| Hipoteza | Formuła | Średni błąd |\n")
    f.write("|----------|---------|-------------|\n")
    for name, err in results.items():
        f.write(f"| {name} | - | {err:.1f}% |\n")
    f.write(f"| **Refined (1/|K(d)|)** | **r = r₀/|K(d)|** | **{mean_err_refined:.1f}%** |\n\n")
    
    f.write("## 3. Finalna Funkcja Mapowania\n")
    f.write(f"**{mapping_name}**\n\n")
    f.write("```python\n")
    f.write("def octave_to_space(d):\n")
    f.write(f"    # {mapping_name}\n")
    f.write("    return mapping_func(d)\n")
    f.write("```\n\n")
    
    f.write("## 4. Weryfikacja na Promieniach Orbitali\n")
    f.write("| n | d | r_pred (Å) | r_exp (Å) | Błąd |\n")
    f.write("|---|---|------------|-----------|------|\n")
    for n in [1, 2, 3, 4]:
        d_n = d_hydrogen[n]
        r_pred = mapping_func(d_n)
        r_exp = r_hydrogen[n]
        err = abs(r_pred - r_exp) / r_exp * 100
        f.write(f"| {n} | {d_n:.2f} | {r_pred:.4f} | {r_exp:.4f} | {err:.1f}% |\n")
    
    f.write(f"\n**Średni błąd:** {final_err:.1f}%\n\n")
    
    f.write("## 5. Wnioski\n")
    if final_err < 20:
        f.write("✅ Mapowanie działa dobrze. Można użyć do widma atomowego.\n")
    elif final_err < 50:
        f.write("🟡 Mapowanie częściowo działa. Wymaga dalszych badań.\n")
    else:
        f.write("❌ Mapowanie nie działa. Wymaga fundamentalnego przeformułowania.\n")

print("Done.")
print(f"Final mapping: {mapping_name}")
print(f"Error: {final_err:.1f}%")
