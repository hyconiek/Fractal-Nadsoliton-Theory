#!/usr/bin/env python3
"""
QW-733/734: HIGGS & SPLITTING REBUTTAL
======================================
Cel: Odpowiedzieć na zarzuty Recenzenta (Prof. Kepticus).

1. QW-733: Gdzie jest Higgs?
   - M_Top = 172.76 GeV
   - M_Higgs = 125.25 GeV
   - Hipoteza: Higgs jest przesunięty o 0.25 oktawy (Rezonans boczny Topa).
   - Expected Ratio for 0.25 octave: 4^(1.52 * 0.25) = 4^0.38 = 1.69? Nie.
   - Pamiętajmy: M ~ 4^(-1.52 * d).
   - d = log4 (M_top / M_higgs) / 1.52

2. QW-734: Analiza Rozszczepienia (Splitting)
   - Tau (d_calc = 2.17) vs Node 2.25
   - Charm (d_calc = 2.33) vs Node 2.25
   - Czy średnia (Center of Mass) trafia w 2.25?
   
"""

import numpy as np

M_TOP = 172760.0
M_HIGGS = 125250.0 # 125.25 GeV
M_TAU = 1776.86
M_CHARM = 1275.0

EXPONENT = 1.52
BASE = 4

print("QW-733/734: PEER REVIEW REBUTTAL DATA")
print("-" * 60)

# --- 1. HIGGS ---
ratio_h = M_TOP / M_HIGGS
d_higgs = (1.0 / EXPONENT) * (np.log(ratio_h) / np.log(BASE))
print(f"1. HIGGS BOSON ANALYSIS")
print(f"   M_Top:   {M_TOP/1000:.2f} GeV")
print(f"   M_Higgs: {M_HIGGS/1000:.2f} GeV")
print(f"   Ratio:   {ratio_h:.4f}")
print(f"   Calculated d: {d_higgs:.4f}")

# Check distance to nearest node
grid = np.arange(-1.0, 1.0, 0.25)
nearest_h = grid[np.argmin(np.abs(grid - d_higgs))]
err_h = abs(d_higgs - nearest_h)
print(f"   Nearest Node: {nearest_h:.2f}")
print(f"   Error: {err_h:.4f}")

if err_h < 0.05:
    print("   STATUS: ✅ HIGGS FITS THE LADDER!")
else:
    print("   STATUS: ❌ HIGGS DOES NOT FIT DIRECTLY.")


# --- 2. SPLITTING (TAU / CHARM) ---
print(f"\n2. SPLITTING ANALYSIS (Level 2.25)")
# Recalculate d
d_tau = (1.0 / EXPONENT) * (np.log(M_TOP / M_TAU) / np.log(BASE))
d_charm = (1.0 / EXPONENT) * (np.log(M_TOP / M_CHARM) / np.log(BASE))

print(f"   Target Node: 2.25")
print(f"   d_Tau:   {d_tau:.4f} (Diff: {d_tau - 2.25:.4f})")
print(f"   d_Charm: {d_charm:.4f} (Diff: {d_charm - 2.25:.4f})")

# Symmetry check
avg_d = (d_tau + d_charm) / 2
print(f"   Average d: {avg_d:.4f}")
print(f"   Dist from 2.25: {avg_d - 2.25:.4f}")

if abs(avg_d - 2.25) < 0.02:
    print("   STATUS: ✅ PERFECT SYMMETRIC SPLITTING")
    print("   (Tau and Charm are split pairs around the 2.25 node)")
else:
    print("   STATUS: ❌ ASYMMETRIC")
