#!/usr/bin/env python3
"""
QW-725: PREDYKCJA MAS KWARKÓW Z GEOMETRII OKTAWOWEJ
====================================================
Purpose: Przewidzieć masy kwarków używając struktury oktawowej i geometrii.

Hipoteza:
  - Kwarki mają strukturę podobną do leptonów, ale z dodatkowym czynnikiem koloru
  - Masa kwarka = M_lepton × f_color × f_generation
  - f_color: czynnik związany z SU(3) color (3 kolory)
  - f_generation: czynnik generacji (jak dla leptonów)

Experimental masses (MeV):
  u: 2.3, d: 4.8, s: 95, c: 1270, b: 4180, t: 172760

Output: raport_qw725_quark_mass_prediction.md
"""

import numpy as np
import datetime

print("="*80)
print("QW-725: PREDYKCJA MAS KWARKÓW")
print("="*80)

# Constants
ALPHA = 4 * np.log(2)  # 2.7726
D1 = 1.3333
D2 = 9.3333
D3 = 17.3333
M0 = 0.23015  # MeV (from lepton model)
W1, W2, W3 = 1, 1, 3

# Experimental quark masses (MeV)
M_U_EXP = 2.3
M_D_EXP = 4.8
M_S_EXP = 95
M_C_EXP = 1270
M_B_EXP = 4180
M_T_EXP = 172760

# Lepton masses (for reference)
M_E = 0.511
M_MU = 105.66
M_TAU = 1776.86

print("\n[1] MASY EKSPERYMENTALNE")
print("-" * 60)
print("| Kwark | Masa (MeV) | Generacja |\n")
print("|-------|------------|-----------|")
print(f"| u | {M_U_EXP:.1f} | 1 |")
print(f"| d | {M_D_EXP:.1f} | 1 |")
print(f"| s | {M_S_EXP:.1f} | 1 |")
print(f"| c | {M_C_EXP:.1f} | 2 |")
print(f"| b | {M_B_EXP:.1f} | 3 |")
print(f"| t | {M_T_EXP:.1f} | 3 |")

# Hypothesis 1: Quark = Lepton × Color Factor
print("\n[2] HYPOTEZA 1: Kwark = Lepton × Czynnik Koloru")
print("-" * 60)

# Color factor: SU(3) has 3 colors, so maybe factor of 3?
# But u/d are lighter than electron, so maybe 1/3?

# Generation mapping
GENERATION_MAP = {
    'u': 1, 'd': 1, 's': 1,  # Gen 1
    'c': 2,  # Gen 2
    'b': 3, 't': 3,  # Gen 3
}

# Try different color factors
color_factors = [1/3, 1/2, 1, 2, 3]

best_factor = None
best_error = 1000

for f_color in color_factors:
    errors = []
    predictions = {}
    
    for quark, gen in GENERATION_MAP.items():
        if gen == 1:
            M_base = M_E
            d = D1
            W = W1
        elif gen == 2:
            M_base = M_MU
            d = D2
            W = W2
        else:  # gen == 3
            M_base = M_TAU
            d = D3
            W = W3
        
        M_pred = M_base * f_color
        
        if quark == 'u':
            M_exp = M_U_EXP
        elif quark == 'd':
            M_exp = M_D_EXP
        elif quark == 's':
            M_exp = M_S_EXP
        elif quark == 'c':
            M_exp = M_C_EXP
        elif quark == 'b':
            M_exp = M_B_EXP
        else:  # t
            M_exp = M_T_EXP
        
        err = abs(M_pred - M_exp) / M_exp * 100
        errors.append(err)
        predictions[quark] = (M_pred, M_exp, err)
    
    mean_err = np.mean(errors)
    if mean_err < best_error:
        best_error = mean_err
        best_factor = f_color
    
    print(f"  f_color = {f_color:.2f}: średni błąd = {mean_err:.1f}%")

print(f"\n  Najlepszy czynnik: f_color = {best_factor:.2f} (błąd: {best_error:.1f}%)")

# Hypothesis 2: Quark mass from octave structure directly
print("\n[3] HYPOTEZA 2: Masa z Struktury Oktawowej")
print("-" * 60)

# Map quarks to octaves
# u, d: octaves 0-2 (light, gen 1)
# s: octaves 3-5 (strange, gen 1 but heavier)
# c: octaves 6-8 (charm, gen 2)
# b: octaves 9-11 (bottom, gen 3)
# t: octaves 12+? (top, gen 3)

QUARK_OCTAVES = {
    'u': 0.5,  # Average of 0-1
    'd': 1.5,  # Average of 1-2
    's': 4.0,  # Average of 3-5
    'c': 7.0,  # Average of 6-8
    'b': 10.0,  # Average of 9-11
    't': 13.0,  # Beyond 12
}

# Use lepton formula: M = M0 × W × d^α
predictions_octave = {}

for quark, d_oct in QUARK_OCTAVES.items():
    gen = GENERATION_MAP[quark]
    W = W1 if gen == 1 else (W2 if gen == 2 else W3)
    
    M_pred = M0 * W * (d_oct ** ALPHA)
    
    if quark == 'u':
        M_exp = M_U_EXP
    elif quark == 'd':
        M_exp = M_D_EXP
    elif quark == 's':
        M_exp = M_S_EXP
    elif quark == 'c':
        M_exp = M_C_EXP
    elif quark == 'b':
        M_exp = M_B_EXP
    else:  # t
        M_exp = M_T_EXP
    
    err = abs(M_pred - M_exp) / M_exp * 100
    predictions_octave[quark] = (M_pred, M_exp, err, d_oct)

print("  | Kwark | d_oct | M_pred (MeV) | M_exp (MeV) | Błąd (%) |")
print("  |-------|-------|--------------|-------------|----------|")
for quark, (M_pred, M_exp, err, d_oct) in predictions_octave.items():
    print(f"  | {quark} | {d_oct:.1f} | {M_pred:.1f} | {M_exp:.1f} | {err:.1f}% |")

mean_err_octave = np.mean([p[2] for p in predictions_octave.values()])
print(f"\n  Średni błąd: {mean_err_octave:.1f}%")

# Hypothesis 3: Top quark from W/Z/Higgs relation
print("\n[4] HYPOTEZA 3: Top Quark z Relacji W/Z/Higgs")
print("-" * 60)

# From QW-651: M_W = M_τ × 4^α
# Maybe: M_t = M_W × something?

M_W = 80.379 * 1000  # GeV → MeV
M_Z = 91.1876 * 1000
M_H = 125.10 * 1000

# Try: M_t = M_W × 2×4^α (from QW-667)
factor_t = 2 * (4 ** ALPHA)
M_t_pred = M_TAU * factor_t
err_t = abs(M_t_pred - M_T_EXP) / M_T_EXP * 100

print(f"  M_t = M_τ × 2×4^α = {M_TAU:.1f} × {factor_t:.1f} = {M_t_pred:.1f} MeV")
print(f"  M_t (exp) = {M_T_EXP:.1f} MeV")
print(f"  Błąd: {err_t:.1f}%")

# Verdict
print("\n[5] WERDYKT")
print("=" * 60)

if mean_err_octave < 50:
    verdict = "🟡 CZĘŚCIOWY SUKCES"
    conclusion = f"Struktura oktawowa daje predykcje z błędem {mean_err_octave:.1f}%."
elif mean_err_octave < 100:
    verdict = "🟠 PRZYBLIŻONA PREDYKCJA"
    conclusion = f"Model daje rząd wielkości, ale błąd {mean_err_octave:.1f}% jest duży."
else:
    verdict = "❌ PORAŻKA"
    conclusion = f"Model nie działa dla kwarków (błąd {mean_err_octave:.1f}%)."

print(f"\n{verdict}")
print(f"{conclusion}")

# Report
report_file = "raport_qw725_quark_mass_prediction.md"
with open(report_file, "w") as f:
    f.write("# RAPORT QW-725: PREDYKCJA MAS KWARKÓW\n")
    f.write(f"**Data:** {datetime.datetime.now()}\n\n")
    
    f.write("## 1. Masowe Eksperymentalne\n")
    f.write("| Kwark | Masa (MeV) | Generacja |\n")
    f.write("|-------|------------|-----------|\n")
    f.write(f"| u | {M_U_EXP:.1f} | 1 |\n")
    f.write(f"| d | {M_D_EXP:.1f} | 1 |\n")
    f.write(f"| s | {M_S_EXP:.1f} | 1 |\n")
    f.write(f"| c | {M_C_EXP:.1f} | 2 |\n")
    f.write(f"| b | {M_B_EXP:.1f} | 3 |\n")
    f.write(f"| t | {M_T_EXP:.1f} | 3 |\n\n")
    
    f.write("## 2. Predykcje z Struktury Oktawowej\n")
    f.write("| Kwark | d_oct | M_pred (MeV) | M_exp (MeV) | Błąd (%) |\n")
    f.write("|-------|-------|--------------|-------------|----------|\n")
    for quark, (M_pred, M_exp, err, d_oct) in predictions_octave.items():
        f.write(f"| {quark} | {d_oct:.1f} | {M_pred:.1f} | {M_exp:.1f} | {err:.1f}% |\n")
    f.write(f"\n**Średni błąd:** {mean_err_octave:.1f}%\n\n")
    
    f.write("## 3. Top Quark\n")
    f.write(f"M_t = M_τ × 2×4^α = {M_t_pred:.1f} MeV\n")
    f.write(f"M_t (exp) = {M_T_EXP:.1f} MeV\n")
    f.write(f"**Błąd:** {err_t:.1f}%\n\n")
    
    f.write("## 4. Wnioski\n")
    f.write(f"### {verdict}\n")
    f.write(f"{conclusion}\n")

print(f"\nRaport zapisany: {report_file}")
