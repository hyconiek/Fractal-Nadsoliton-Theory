#!/usr/bin/env python3
"""
QW-643: INTERPRETACJA WARSTW - Gdzie mierzymy masę?
Purpose: Sprawdzić czy masa jest właściwością warstwy gdzie ISTNIEJE czy gdzie WIDZIMY
Date: 2025-12-06
"""

import numpy as np

BETA_TORS = 0.01
M_PLANCK_GeV = 1.2209e19
M_ELECTRON_EXP_MeV = 0.511

print("="*80)
print("QW-643: INTERPRETACJA WARSTW - Gdzie mierzymy masę?")
print("="*80)

print(f"\n📐 Sytuacja:")
print(f"   Elektron ISTNIEJE w N=10 (skala atomowa)")
print(f"   My WIDZIMY go w N=20 (nasza skala makroskopowa)")
print(f"   Eksperymentalna masa: m_e = {M_ELECTRON_EXP_MeV:.4f} MeV")

# Opcja 1: Masa jest właściwością warstwy gdzie ISTNIEJE (N=10)
m_scale_10 = M_PLANCK_GeV * (BETA_TORS ** 10)
print(f"\n🔍 OPCJA 1: Masa jest właściwością warstwy gdzie ISTNIEJE (N=10)")
print(f"   Skala masy: m_Planck × β^10 = {m_scale_10:.4e} GeV = {m_scale_10*1000:.4e} MeV")
print(f"   To jest skala masy w warstwie N=10")

# Opcja 2: Masa jest właściwością warstwy gdzie WIDZIMY (N=20)
m_scale_20 = M_PLANCK_GeV * (BETA_TORS ** 20)
print(f"\n🔍 OPCJA 2: Masa jest właściwością warstwy gdzie WIDZIMY (N=20)")
print(f"   Skala masy: m_Planck × β^20 = {m_scale_20:.4e} GeV = {m_scale_20*1000:.4e} MeV")
print(f"   To jest skala masy w warstwie N=20")

# Sprawdźmy, która skala jest bliższa eksperymentowi
ratio_10 = M_ELECTRON_EXP_MeV / (m_scale_10 * 1000)
ratio_20 = M_ELECTRON_EXP_MeV / (m_scale_20 * 1000)

print(f"\n📊 PORÓWNANIE:")
print(f"   m_exp / m_scale(N=10) = {ratio_10:.4e}")
print(f"   m_exp / m_scale(N=20) = {ratio_20:.4e}")

if ratio_10 < 1 and ratio_20 > 1:
    print(f"\n✅ WNIOSEK:")
    print(f"   m_exp jest POMIĘDZY m_scale(N=10) i m_scale(N=20)")
    print(f"   To sugeruje, że masa jest właściwością warstwy gdzie ISTNIEJE (N=10)")
    print(f"   Ale jest 'przeskalowana' przez transformację do N=20")
elif abs(ratio_10 - 1) < abs(ratio_20 - 1):
    print(f"\n✅ WNIOSEK:")
    print(f"   m_exp jest bliżej m_scale(N=10)")
    print(f"   To sugeruje, że masa jest właściwością warstwy gdzie ISTNIEJE (N=10)")
else:
    print(f"\n✅ WNIOSEK:")
    print(f"   m_exp jest bliżej m_scale(N=20)")
    print(f"   To sugeruje, że masa jest właściwością warstwy gdzie WIDZIMY (N=20)")

# Sprawdźmy transformację
transform_factor = M_ELECTRON_EXP_MeV / (m_scale_10 * 1000)
print(f"\n🔀 TRANSFORMACJA:")
print(f"   Jeśli używamy skali N=10 (gdzie ISTNIEJE):")
print(f"   m_observed = m_intrinsic × {transform_factor:.4e}")
print(f"   = m_Planck × β^10 × {transform_factor:.4e}")

print("\n" + "="*80)

