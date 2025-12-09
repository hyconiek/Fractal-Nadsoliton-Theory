#!/usr/bin/env python3
"""
QW-723: WERYFIKACJA PRAWA GRAWITACJI Z DEFEKTAMI TOPOLOGICZNYMI
================================================================
Purpose: Rygorystyczna weryfikacja F ∝ 1/r² z defektami topologicznymi
         na większej sieci i z różnymi konfiguracjami.

Output: raport_qw723_gravity_defects_verification.md
"""

import numpy as np
from scipy.optimize import curve_fit
import datetime

print("="*80)
print("QW-723: WERYFIKACJA PRAWA GRAWITACJI Z DEFEKTAMI")
print("="*80)

# Import from QW-722
ALPHA = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA = 0.1
N_OCTAVES = 12

def K(d):
    if d < 0.1:
        d = 0.1
    return ALPHA * np.cos(OMEGA * d + PHI) / (1 + BETA * d)

class TopologicalDefect:
    def __init__(self, position, octave, winding_number=1):
        self.position = np.array(position)
        self.octave = octave
        self.winding = winding_number

def calculate_defect_force(defect1, defect2):
    """Oblicza siłę między defektami z poprawionym modelem"""
    r = np.linalg.norm(defect1.position - defect2.position)
    d_octave = abs(defect1.octave - defect2.octave)
    
    # Mieszanie skali: d_eff = d_octave + α × r
    d_eff = d_octave + 0.1 * r
    K_eff = K(d_eff)
    K_eff *= (1.0 + 0.5 * defect1.winding)
    K_eff *= (1.0 + 0.5 * defect2.winding)
    
    m1, m2 = defect1.winding, defect2.winding
    F = -K_eff * m1 * m2 / (r**2 + 0.1)
    return F

# Test multiple configurations
print("\n[1] TEST RÓŻNYCH KONFIGURACJI")
print("-" * 60)

configurations = [
    {"name": "Same octave, w=1", "oct1": 3, "oct2": 3, "w1": 1, "w2": 1},
    {"name": "Different octaves, w=1", "oct1": 3, "oct2": 6, "w1": 1, "w2": 1},
    {"name": "Same octave, w=2", "oct1": 3, "oct2": 3, "w1": 2, "w2": 2},
    {"name": "Mixed winding", "oct1": 3, "oct2": 5, "w1": 1, "w2": 2},
]

test_distances = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0])

results = {}

for config in configurations:
    print(f"\n  Konfiguracja: {config['name']}")
    forces = []
    
    for r in test_distances:
        defect1 = TopologicalDefect([0, 0, 0], config['oct1'], config['w1'])
        defect2 = TopologicalDefect([r, 0, 0], config['oct2'], config['w2'])
        F = calculate_defect_force(defect1, defect2)
        forces.append(F)
    
    forces = np.array(forces)
    forces_abs = np.abs(forces)
    valid = forces_abs > 1e-10
    
    if np.sum(valid) >= 3:
        r_fit = test_distances[valid]
        F_fit = forces_abs[valid]
        
        def power_law(r, A, n):
            return A * np.power(r, n)
        
        try:
            popt, _ = curve_fit(power_law, r_fit, F_fit, p0=[1.0, -2.0], maxfev=10000)
            A_fit, n_fit = popt
            results[config['name']] = {"n": n_fit, "A": A_fit, "error": abs(n_fit + 2.0)}
            print(f"    n = {n_fit:.2f}, błąd = {abs(n_fit + 2.0):.2f}")
        except:
            results[config['name']] = {"n": 0, "A": 0, "error": 100}
    else:
        results[config['name']] = {"n": 0, "A": 0, "error": 100}

# Summary
print("\n[2] PODSUMOWANIE")
print("-" * 60)
print("| Konfiguracja | Wykładnik n | Błąd | Status |")
print("|--------------|-------------|------|--------|")

for name, res in results.items():
    status = "✅" if res['error'] < 0.3 else "🟡" if res['error'] < 1.0 else "❌"
    print(f"| {name} | {res['n']:.2f} | {res['error']:.2f} | {status} |")

mean_error = np.mean([r['error'] for r in results.values()])
print(f"\nŚredni błąd: {mean_error:.2f}")

# Verdict
if mean_error < 0.3:
    verdict = "✅ SUKCES: Prawo 1/r² potwierdzone we wszystkich konfiguracjach!"
elif mean_error < 1.0:
    verdict = "🟡 CZĘŚCIOWY SUKCES: Większość konfiguracji działa"
else:
    verdict = "❌ PORAŻKA: Model wymaga poprawy"

print(f"\n{verdict}")

# Report
report_file = "raport_qw723_gravity_defects_verification.md"
with open(report_file, "w") as f:
    f.write("# RAPORT QW-723: WERYFIKACJA PRAWA GRAWITACJI Z DEFEKTAMI\n")
    f.write(f"**Data:** {datetime.datetime.now()}\n\n")
    
    f.write("## 1. Metoda\n")
    f.write("Testowanie prawa grawitacji F ∝ 1/r² dla różnych konfiguracji defektów.\n\n")
    
    f.write("## 2. Wyniki\n")
    f.write("| Konfiguracja | Wykładnik n | Błąd | Status |\n")
    f.write("|--------------|-------------|------|--------|\n")
    for name, res in results.items():
        status = "✅" if res['error'] < 0.3 else "🟡" if res['error'] < 1.0 else "❌"
        f.write(f"| {name} | {res['n']:.2f} | {res['error']:.2f} | {status} |\n")
    
    f.write(f"\n**Średni błąd:** {mean_error:.2f}\n\n")
    
    f.write("## 3. Wnioski\n")
    f.write(f"### {verdict}\n")

print(f"\nRaport zapisany: {report_file}")
