#!/usr/bin/env python3
"""
QW-642: TRANSFORMACJA W GÓRĘ - Elektron z N=10 do N=20
Purpose: Sprawdzić czy transformacja w górę wymaga wzmocnienia
Date: 2025-12-06
"""

import numpy as np

BETA_TORS = 0.01
M_PLANCK_GeV = 1.2209e19
M_ELECTRON_EXP_MeV = 0.511

print("="*80)
print("QW-642: TRANSFORMACJA W GÓRĘ (N=10 → N=20)")
print("="*80)

# Elektron ISTNIEJE w N=10
m_e_intrinsic_GeV = M_PLANCK_GeV * (BETA_TORS ** 10)

print(f"\n📊 MASY:")
print(f"   Elektron ISTNIEJE w N=10:")
print(f"   m_e(intrinsic) = m_Planck × β^10 = {m_e_intrinsic_GeV:.4e} GeV")
print(f"   = {m_e_intrinsic_GeV * 1000:.4e} MeV")

# Gdy widzimy go w N=20, transformacja może być:
# 1. Tłumienie: m_observed = m_intrinsic × β^10
# 2. Wzmocnienie: m_observed = m_intrinsic × β^(-10)
# 3. Bez zmian: m_observed = m_intrinsic

print(f"\n🔀 OPCJE TRANSFORMACJI:")
print(f"   1. Tłumienie: m_obs = m_int × β^10 = {m_e_intrinsic_GeV * (BETA_TORS**10) * 1000:.4e} MeV")
print(f"   2. Wzmocnienie: m_obs = m_int × β^(-10) = {m_e_intrinsic_GeV * (BETA_TORS**(-10)) * 1000:.4e} MeV")
print(f"   3. Bez zmian: m_obs = m_int = {m_e_intrinsic_GeV * 1000:.4e} MeV")
print(f"\n   Eksperyment: m_e = {M_ELECTRON_EXP_MeV:.4f} MeV")

# Sprawdźmy, która transformacja daje poprawny wynik
factor_needed = M_ELECTRON_EXP_MeV / (m_e_intrinsic_GeV * 1000)
print(f"\n🔍 CZYNNIK POTRZEBNY:")
print(f"   m_exp / m_intrinsic = {factor_needed:.4e}")
print(f"   β^(-10) = {BETA_TORS**(-10):.4e}")
print(f"   Czy factor ≈ β^(-10)? {abs(factor_needed - BETA_TORS**(-10)) / BETA_TORS**(-10) * 100:.2f}% różnicy")

if abs(factor_needed - BETA_TORS**(-10)) / BETA_TORS**(-10) < 0.1:
    print(f"\n✅ TAK! Transformacja w górę WZMACNIA sygnał!")
    print(f"   m_observed = m_intrinsic × β^(-10)")
    print(f"   To znaczy: sygnał z N=10 jest WZMACNIANY gdy idzie w górę do N=20")
else:
    print(f"\n❌ NIE! Transformacja nie jest prostym β^(-10)")

print("\n" + "="*80)

