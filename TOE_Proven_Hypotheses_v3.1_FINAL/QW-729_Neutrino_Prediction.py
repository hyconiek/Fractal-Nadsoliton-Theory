#!/usr/bin/env python3
"""
QW-729: NEUTRINO MASS PREDICTION
================================
Cel: Sprawdzić, czy drabina geometryczna (Base 4, d ~ 0.25) przewidujeasy neutrin.

Model:
  - M(d) = M_top * 4^(-1.52 * d)
  - M_top = 172760 MeV
  - Szukamy d, dla których M(d) < 1 eV (1e-6 MeV).

Oczekiwania:
  - Jeśli neutrina są na tej samej drabinie, powinniśmy znaleźć stabilne węzły (n * 0.25)
    w zakresie d odpowiadającym masom ~0.05 eV.
  - Z wstępnych szacunków d ~ 13-14.

"""

import numpy as np

M_TOP_MEV = 172760.0
EXPONENT = 1.52
BASE = 4

def mass_from_d(d):
    # M = M_ref * Base^(-exponent * d)
    # r/r_ref = Base^d
    # M/M_ref = (r/r_ref)^-exponent = (Base^d)^-exponent = Base^(-exponent * d)
    return M_TOP_MEV * (BASE ** (-EXPONENT * d))

print("QW-729: NEUTRINO MASS PREDICTION")
print("================================")
print(f"Top Mass: {M_TOP_MEV} MeV")
print(f"Formula: M(d) = M_top * 4^(-1.52 * d)")
print("-" * 60)
print(f"{'d (Octave)':<10} | {'Mass (MeV)':<15} | {'Mass (eV)':<15} | {'Note'}")
print("-" * 60)

# Scan range where mass passes through eV scale
# 0.1 eV = 1e-7 MeV.
# M_top ~ 1e5 MeV. Ratio ~ 1e12.
# 4^(1.5 * d) ~ 1e12 => 4^20 ~ 1e12. 1.5d ~ 20 => d ~ 13.

scan_range = np.arange(12.0, 16.0, 0.25)

candidates = []

for d in scan_range:
    m_mev = mass_from_d(d)
    m_ev = m_mev * 1e6
    
    note = ""
    if 0.001 < m_ev < 1.0:
        note = "✅ NEUTRINO RANGE"
        candidates.append((d, m_ev))
    elif m_ev < 0.001:
        note = "Too light"
    else:
        note = "Too heavy"
        
    print(f"{d:<10.2f} | {m_mev:<15.4e} | {m_ev:<15.4f} | {note}")

print("-" * 60)
print("\nWNIOSKI:")

if len(candidates) > 0:
    print(f"Znaleziono {len(candidates)} kandydatów w zakresie neutrinowym (0.001 - 1 eV):")
    for d, m_ev in candidates:
        print(f"  d={d:.2f}: {m_ev:.4f} eV")
        
    # Sprawdzenie Delta m^2 (Atmospheric / Solar)
    # Atm splitting: ~0.05 eV (mass difference or mass scale)
    # Solar splitting: ~0.009 eV
    
    # Czy różnica między szczeblami pasuje do eksperymentów?
    if len(candidates) >= 2:
        m1 = candidates[-1][1] # Lżejsze
        m2 = candidates[-2][1] # Cięższe
        print(f"\nPrzykładowa różnica między d={candidates[-2][0]} i d={candidates[-1][0]}:")
        print(f"  dm = |{m1:.4f} - {m2:.4f}| = {abs(m1-m2):.4f} eV")
        print(f"  dm^2 = |{m1**2:.6f} - {m2**2:.6f}| = {abs(m1**2 - m2**2):.6f} eV^2")
        print("  (Atmospheric target: ~2.5e-3 eV^2)")

else:
    print("Nie znaleziono kandydatów w zakresie neutrinowym (0.001 - 1 eV).")
    print("Neutrina mogą wymagać innego mechanizmu (Seesaw?).")
