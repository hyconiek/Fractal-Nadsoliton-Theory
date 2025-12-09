#!/usr/bin/env python3
"""
QW-644: TRANSFORMACJA MASY między warstwami
Purpose: Sprawdzić jak masa transformuje się z N=10 (gdzie istnieje) do N=20 (gdzie mierzymy)
Date: 2025-12-06
"""

import numpy as np

BETA_TORS = 0.01
M_PLANCK_GeV = 1.2209e19
M_ELECTRON_EXP_MeV = 0.511  # To jest masa zmierzona na NASZEJ warstwie N=20!

print("="*80)
print("QW-644: TRANSFORMACJA MASY między warstwami")
print("="*80)

print(f"\n📐 Sytuacja:")
print(f"   Elektron ISTNIEJE w N=10 (skala atomowa)")
print(f"   My MIERZYMY masę w N=20 (nasza skala makroskopowa)")
print(f"   Zmierzona masa: m_e(observed at N=20) = {M_ELECTRON_EXP_MeV:.4f} MeV")
print(f"   To jest odczyt na NASZEJ warstwie!")

# Skale masy w różnych warstwach
m_scale_10 = M_PLANCK_GeV * (BETA_TORS ** 10)  # N=10 (gdzie istnieje)
m_scale_20 = M_PLANCK_GeV * (BETA_TORS ** 20)  # N=20 (gdzie mierzymy)

print(f"\n📊 SKALE MASY:")
print(f"   N=10 (gdzie ISTNIEJE): m_scale = {m_scale_10:.4e} GeV = {m_scale_10*1000:.4e} MeV")
print(f"   N=20 (gdzie MIERZYMY): m_scale = {m_scale_20:.4e} GeV = {m_scale_20*1000:.4e} MeV")

# Jeśli masa zmierzona jest w skali N=20, to:
# m_e(observed) = m_Planck × β^20 × (czynniki geometryczne)

# Sprawdźmy, jakie czynniki geometryczne dają 0.511 MeV
# 0.511 MeV = m_Planck × β^20 × factors
# factors = 0.511 MeV / (m_Planck × β^20)

factors_needed = M_ELECTRON_EXP_MeV / (m_scale_20 * 1000)
print(f"\n🔍 CZYNNIKI GEOMETRYCZNE POTRZEBNE:")
print(f"   Jeśli używamy skali N=20 (gdzie mierzymy):")
print(f"   factors = m_exp / (m_Planck × β^20) = {factors_needed:.4e}")

# Sprawdźmy też, jaka byłaby masa w N=10
# m_e(intrinsic at N=10) = m_e(observed at N=20) × (β^10 / β^20) = m_e(observed) / β^10
m_e_intrinsic = M_ELECTRON_EXP_MeV / (BETA_TORS ** 10)
print(f"\n🔍 MASY W RÓŻNYCH WARSTWACH:")
print(f"   m_e(observed at N=20) = {M_ELECTRON_EXP_MeV:.4f} MeV")
print(f"   m_e(intrinsic at N=10) = m_obs / β^10 = {m_e_intrinsic:.4e} MeV")
print(f"   Stosunek: m_intrinsic / m_observed = {m_e_intrinsic / M_ELECTRON_EXP_MeV:.4e} = β^(-10)")

# Sprawdźmy, czy formuła z β^20 daje sensowny wynik
# m = m_Planck × β^20 × |W| × κ^(1/12) × A_res × I_proc
# Musimy znaleźć I_proc, które daje 0.511 MeV

# Z QW-639 wiemy:
W = 1
kappa = 4 * np.log(2) / (np.pi/4 * np.pi/6)
octave_amp = kappa ** (1/12)
# A_res i I_proc musimy obliczyć lub oszacować

print(f"\n📐 FORMULA Z SKALĄ N=20:")
print(f"   m_e = m_Planck × β^20 × |W| × κ^(1/12) × A_res × I_proc")
print(f"   Jeśli m_e = {M_ELECTRON_EXP_MeV:.4f} MeV:")
print(f"   I_proc × A_res = m_e / (m_Planck × β^20 × |W| × κ^(1/12))")
print(f"   = {factors_needed:.4e}")

# Porównajmy z obecną wartością I_proc
I_proc_current = 0.016033  # Z QW-639
A_res_current = 0.267821   # Z QW-639
product_current = I_proc_current * A_res_current

print(f"\n📊 PORÓWNANIE:")
print(f"   Obecne: I_proc × A_res = {product_current:.6f}")
print(f"   Potrzebne: I_proc × A_res = {factors_needed:.4e}")
print(f"   Różnica: {abs(product_current - factors_needed) / factors_needed * 100:.2f}%")

if abs(product_current - factors_needed) / factors_needed < 0.5:
    print(f"\n✅ TAK! Formuła z β^20 jest POPRAWNA!")
    print(f"   Masa jest mierzona na NASZEJ warstwie N=20")
    print(f"   Więc używamy skali m_Planck × β^20")
else:
    print(f"\n❌ NIE! Formuła z β^20 wymaga innych wartości I_proc/A_res")

print("\n" + "="*80)


