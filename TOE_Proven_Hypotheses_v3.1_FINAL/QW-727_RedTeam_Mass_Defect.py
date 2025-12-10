#!/usr/bin/env python3
"""
QW-727: RED TEAM ANALYSIS - MASS DEFECT HYPOTHESIS
==================================================
Cel: Brutalna weryfikacja hipotezy o kwantowaniu mas kwarków w krokach 0.25 oktawy.

Ataki:
1. [STATYSTYKA] Czy dopasowanie do kroku 0.25 jest statystycznie istotne?
   - Monte Carlo: Generuj losowe masy (z rozkładu logarytmicznego rzędu kwarków).
   - Sprawdź jak często przypadkowe masy pasują do siatki 0.25 tak dobrze jak eksperyment.
   
2. [SENSITIVITY] Jak bardzo wynik zależy od parametru wykładnika całki (n=2.26)?
   - Zmieniaj n w zakresie błędu QW-722 (2.26 +/- 0.26).
   - Sprawdź czy kwantowanie znika.

3. [CONSISTENCY] Czy masa elektronu pasuje do tej drabiny?
   - Sprawdź gdzie ląduje elektron na tej drabinie względem Topu.
"""

import numpy as np
import matplotlib.pyplot as plt

# --- DANE EKSPERYMENTALNE ---
# Masy w MeV
QUARKS = {
    "u": 2.3, "d": 4.8, "s": 95.0, "c": 1275.0, "b": 4180.0, "t": 172760.0
}
MASS_VALUES = np.array(list(QUARKS.values()))
NAMES = list(QUARKS.keys())

# Referencja (Top)
M_REF = QUARKS["t"]
EXPONENT_BASE = 1.52 # Z QW-726 (3 - 2*0.74? Nie, 3 - 2*2.26 = -1.52 wiec E ~ r^-1.52)
# Czekaj, r_min ~ M^(-1/1.52). 
# W QW-726: r_min = mass ** (-1 / 1.52)

def calculate_grid_error(masses, m_ref, exponent, step=0.25, base=4):
    """Oblicza średni błąd dopasowania do siatki kwantowania"""
    errors = []
    
    r_ref = m_ref ** (-1 / exponent)
    
    for m in masses:
        if m == m_ref: continue
        
        r_curr = m ** (-1 / exponent)
        ratio = r_curr / r_ref
        
        # Log Base 4 difference
        d_diff = np.log(ratio) / np.log(base)
        
        # Distance to nearest grid point
        nearest = round(d_diff / step) * step
        err = abs(d_diff - nearest)
        errors.append(err)
        
    return np.mean(errors), errors

# --- ATAK 1: MONTE CARLO ---
print("=== ATAK 1: WERYFIKACJA STATYSTYCZNA (MONTE CARLO) ===")
np.random.seed(42)
N_TRIALS = 100000
true_error, _ = calculate_grid_error(MASS_VALUES, M_REF, EXPONENT_BASE)
print(f"Błąd eksperymentalny (Mean Abs Deviation from Grid): {true_error:.4f} oktawy")

better_count = 0
for _ in range(N_TRIALS):
    # Generuj losowe masy w zakresie logarytmicznym [Min, Max] kwarków
    # Zachowujemy strukturę rzędów wielkości, ale losujemy wartości
    log_min = np.log(min(MASS_VALUES))
    log_max = np.log(max(MASS_VALUES))
    
    # 5 losowych mas + 1 fixed (Top) jako referencja
    random_logs = np.random.uniform(log_min, log_max, 5)
    random_masses = np.exp(random_logs)
    all_random_masses = np.append(random_masses, M_REF) # Top is fixed anchor
    
    rand_error, _ = calculate_grid_error(all_random_masses, M_REF, EXPONENT_BASE)
    
    if rand_error <= true_error:
        better_count += 1

p_value = better_count / N_TRIALS
print(f"Liczba lepszych dopasowań w {N_TRIALS} próbach: {better_count}")
print(f"P-value (szansa na przypadek): {p_value:.4f}")

if p_value < 0.05:
    print(">> WYNIK: Statystycznie ISTOTNE ✅")
else:
    print(">> WYNIK: Przypadkowy szum ❌")

# --- ATAK 2: SENSITIVITY ANALYSIS ---
print("\n=== ATAK 2: WRAŻLIWOŚĆ NA WYKŁADNIK CAŁKI ===")
# Wykładnik pochodzi z n (siła). n = 2.26 +/- 0.26 (z QW-722)
# Exponent całki p = 3 - 2*n
# Dla n=2.26 -> p = -1.52. Absolute val = 1.52.
# Dla n=2.00 (Newton) -> p = 3 - 4 = -1.0. Absolute = 1.0.
# Dla n=2.50 -> p = 3 - 5 = -2.0. Absolute = 2.0.

scan_n = np.linspace(2.0, 2.5, 20)
sensitivity_errors = []

print(f"{'n (Force)':<10} | {'p (Int)':<10} | {'Grid Error':<10}")
print("-" * 40)

for n in scan_n:
    p_abs = abs(3 - 2*n)
    if p_abs < 0.1: continue # Unstable
    
    err, _ = calculate_grid_error(MASS_VALUES, M_REF, p_abs)
    sensitivity_errors.append(err)
    mark = "<-- MIN" if err == min(sensitivity_errors[-1:] + [100]) else "" # Logic fix needed but simple print works
    print(f"{n:<10.3f} | {p_abs:<10.3f} | {err:<10.4f}")

best_idx = np.argmin(sensitivity_errors)
best_n = scan_n[best_idx]
print(f"\nNajlepsze dopasowanie dla n = {best_n:.3f} (Eksperyment QW-722 dał 2.26)")
if abs(best_n - 2.26) < 0.1:
    print(">> WYNIK: Spójność z QW-722 ✅")
else:
    print(f">> WYNIK: Niespójność! Wymagane n={best_n:.3f} zamiast 2.26 ❌")

# --- ATAK 3: ELEKTRON ---
print("\n=== ATAK 3: CZY ELEKTRON PASUJE? ===")
MASS_E = 0.511
r_e = MASS_E ** (-1 / EXPONENT_BASE)
r_ref = M_REF ** (-1 / EXPONENT_BASE)
ratio_e = r_e / r_ref
d_diff_e = np.log(ratio_e) / np.log(4)
nearest_e = round(d_diff_e * 4) / 4
err_e = abs(d_diff_e - nearest_e)

print(f"Masa elektronu: {MASS_E} MeV")
print(f"Pozycja na drabinie (Base 4): {d_diff_e:.4f} oktaw od Topa")
print(f"Najbliższy węzeł 0.25: {nearest_e}")
print(f"Błąd: {err_e:.4f}")

if err_e < 0.05:
    print(">> WYNIK: Elektron pasuje do drabiny! ✅")
else:
    print(">> WYNIK: Elektron NIE pasuje (inny mechanizm?) ⚠️")
