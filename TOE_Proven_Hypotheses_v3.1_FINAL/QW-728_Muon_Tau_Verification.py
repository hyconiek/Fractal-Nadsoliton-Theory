#!/usr/bin/env python3
"""
QW-728: WERYFIKACJA MIONU I TAU (Moment Prawdy)
===============================================
Cel: Sprawdzić czy leptony pośrednie (Mion, Tau) pasują do drabiny kwarkowej (Base 4, step 0.25).
     Jeśli tak -> Unifikacja potwierdzona.
     Jeśli nie -> Teoria ma problem.

Dane:
  - Top Quark (Referencja): 172760 MeV
  - Tau Lepton: 1776.86 MeV
  - Muon Lepton: 105.66 MeV

Model:
  - Masa M to energia pola defektu: M ~ r_core^-1.52
  - Promień rdzenia r_core jest kwantowany w Bazie 4: r ~ 4^d
  - Oczekujemy d = k * 0.25

Obliczenia:
  1. Ratio r/r_top = (M_top / M)^(1/1.52)
  2. Octave diff d = log4(ratio)
  3. Sprawdzenie bliskości do gridu 0.25
"""

import numpy as np

# Masy w MeV
M_TOP = 172760.0
M_TAU = 1776.86
M_MUON = 105.66
M_ELECTRON = 0.510999

# Wykladnik z grawitacji (QW-722/726)
EXPONENT = 1.52 # pochodzi z n=2.26 (F~r^-2.26 => E~r^-1.52)

particles = [
    ("Top (Ref)", M_TOP),
    ("Tau", M_TAU),
    ("Muon", M_MUON),
    ("Electron", M_ELECTRON)
]

print(f"QW-728: WERYFIKACJA MIONU I TAU")
print(f"Exponent (Defect Energy Integral): {EXPONENT}")
print(f"Base: 4 (Information)")
print("-" * 65)
print(f"{'Particle':<12} | {'Mass (MeV)':<10} | {'d (calc)':<10} | {'Nearest 0.25':<12} | {'Error':<10} | {'Status'}")
print("-" * 65)

results = []

for name, mass in particles:
    if mass == M_TOP:
        d_val = 0.0
        nearest = 0.0
        err = 0.0
        status = "REF"
    else:
        # r/r_top = (M_top / M)^(1/1.52)
        # d = log4 (r/r_top)
        #   = log4 ( (M_top/M)^(1/1.52) )
        #   = (1/1.52) * log4(M_top/M)
        
        ratio_mass = M_TOP / mass
        d_val = (1.0 / EXPONENT) * (np.log(ratio_mass) / np.log(4))
        
        nearest = round(d_val * 4) / 4
        err = abs(d_val - nearest)
        
        status = "✅" if err < 0.05 else "⚠️" if err < 0.1 else "❌"
    
    print(f"{name:<12} | {mass:<10.2f} | {d_val:<10.4f} | {nearest:<12.2f} | {err:<10.4f} | {status}")
    results.append((name, err))

# Analiza Tau
tau_res = next(r for r in results if r[0] == "Tau")
muon_res = next(r for r in results if r[0] == "Muon")

print("-" * 65)
print("\nANALIZA WYNIKÓW:")

if tau_res[1] < 0.05:
    print(f"1. Tau: SUKCES (Błąd {tau_res[1]:.4f} oktawy)")
else:
    print(f"1. Tau: PORAŻKA (Błąd {tau_res[1]:.4f} oktawy)")

if muon_res[1] < 0.05:
    print(f"2. Mion: SUKCES (Błąd {muon_res[1]:.4f} oktawy)")
else:
    print(f"2. Mion: PORAŻKA (Błąd {muon_res[1]:.4f} oktawy)")

# Sprawdzenie kroku między Muon a Electron
# d_muon ~ 4.75? zależy od wyniku
# d_electron ~ 6.00
