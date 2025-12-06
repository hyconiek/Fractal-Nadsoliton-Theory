#!/usr/bin/env python3
"""
QW-645: TRANSFORMACJA ODWROTNA - Masa w N=10 vs N=20
Purpose: Sprawdzić czy transformacja wzmacnia czy tłumi masę
Date: 2025-12-06
"""

import numpy as np

BETA_TORS = 0.01
M_PLANCK_GeV = 1.2209e19
M_ELECTRON_EXP_MeV = 0.511  # Masa zmierzona na N=20

print("="*80)
print("QW-645: TRANSFORMACJA ODWROTNA - Masa w N=10 vs N=20")
print("="*80)

print(f"\n📐 Sytuacja:")
print(f"   Elektron ISTNIEJE w N=10")
print(f"   My MIERZYMY masę w N=20: m_e(observed) = {M_ELECTRON_EXP_MeV:.4f} MeV")
print(f"   To jest odczyt na NASZEJ warstwie!")

# Opcja 1: Transformacja tłumi (idzie w górę)
# m_e(observed at N=20) = m_e(intrinsic at N=10) × β^10
m_e_intrinsic_damped = M_ELECTRON_EXP_MeV / (BETA_TORS ** 10)
print(f"\n🔍 OPCJA 1: Transformacja TŁUMI (idzie w górę):")
print(f"   m_e(intrinsic at N=10) = m_e(observed) / β^10")
print(f"   = {M_ELECTRON_EXP_MeV:.4f} / {BETA_TORS**10:.4e} = {m_e_intrinsic_damped:.4e} MeV")

# Opcja 2: Transformacja wzmacnia (idzie w górę)
# m_e(observed at N=20) = m_e(intrinsic at N=10) × β^(-10)
m_e_intrinsic_amplified = M_ELECTRON_EXP_MeV * (BETA_TORS ** 10)
print(f"\n🔍 OPCJA 2: Transformacja WZMACNIA (idzie w górę):")
print(f"   m_e(intrinsic at N=10) = m_e(observed) × β^10")
print(f"   = {M_ELECTRON_EXP_MeV:.4f} × {BETA_TORS**10:.4e} = {m_e_intrinsic_amplified:.4e} MeV")

# Sprawdźmy, która opcja daje sensowną masę w N=10
# Masa w N=10 powinna być w skali m_Planck × β^10 = 122.09 MeV
m_scale_10 = M_PLANCK_GeV * (BETA_TORS ** 10) * 1000  # MeV

print(f"\n📊 PORÓWNANIE ZE SKALĄ N=10:")
print(f"   Skala masy w N=10: m_Planck × β^10 = {m_scale_10:.4f} MeV")
print(f"   Opcja 1 (tłumienie): {m_e_intrinsic_damped:.4e} MeV")
print(f"   Opcja 2 (wzmocnienie): {m_e_intrinsic_amplified:.4e} MeV")

# Sprawdźmy, jaki czynnik transformacji daje sensowną masę w N=10
# m_e(intrinsic) = m_e(observed) × factor
# factor = m_scale_10 / m_e(observed)
factor_needed = m_scale_10 / M_ELECTRON_EXP_MeV
print(f"\n🔍 CZYNNIK TRANSFORMACJI:")
print(f"   Aby masa w N=10 była w skali m_Planck × β^10:")
print(f"   factor = m_scale(N=10) / m_e(observed) = {factor_needed:.4e}")
print(f"   β^(-10) = {BETA_TORS**(-10):.4e}")
print(f"   Czy factor ≈ β^(-10)? {abs(factor_needed - BETA_TORS**(-10)) / BETA_TORS**(-10) * 100:.2f}% różnicy")

if abs(factor_needed - BETA_TORS**(-10)) / BETA_TORS**(-10) < 0.1:
    print(f"\n✅ TAK! Transformacja WZMACNIA o β^(-10)!")
    print(f"   m_e(intrinsic at N=10) = m_e(observed at N=20) × β^(-10)")
    print(f"   = {M_ELECTRON_EXP_MeV:.4f} × {BETA_TORS**(-10):.4e} = {M_ELECTRON_EXP_MeV * BETA_TORS**(-10):.4f} MeV")
    transform_type = "amplification"
else:
    print(f"\n❓ Transformacja nie jest prostym β^(-10)")
    print(f"   Może jest bardziej skomplikowana?")
    transform_type = "complex"

print(f"\n📐 FORMULA:")
if transform_type == "amplification":
    print(f"   m_e(observed at N=20) = m_Planck × β^20 × (czynniki)")
    print(f"   m_e(intrinsic at N=10) = m_e(observed) × β^(-10)")
    print(f"   = m_Planck × β^10 × (czynniki)")
    print(f"   Transformacja WZMACNIA sygnał gdy idzie w górę!")
else:
    print(f"   Transformacja jest bardziej skomplikowana niż proste β^10")

print("\n" + "="*80)

