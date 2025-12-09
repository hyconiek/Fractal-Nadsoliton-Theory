#!/usr/bin/env python3
"""
QW-721: WIDMO WODORU Z MAPOWANIEM OKTAW → PRZESTRZEŃ
====================================================
Purpose: Użyć mapowania z QW-720 do poprawnego modelowania widma wodoru.

Kluczowe odkrycie QW-720:
  - Liniowe mapowanie: r(d) = d × a_B / d_0 daje błąd 34.6%
  - Ale: promienie orbitali to r_n = n² × a_B
  - Więc: d_n = d_1 + 8(n-1) → r_n powinno być n² × a_B

Nowe podejście:
  - Użyć bezpośrednio n (liczba kwantowa) zamiast d
  - r_n = n² × a_B (Bohr model)
  - Mapować n → d przez: d_n = d_1 + 8(n-1)
  - Użyć r_n w Hamiltonianie zamiast d_n
"""

import numpy as np
from scipy.linalg import eigh
import datetime

print("="*80)
print("QW-721: WIDMO WODORU Z MAPOWANIEM OKTAW → PRZESTRZEŃ")
print("="*80)

# --- Constants ---
ALPHA = 4 * np.log(2)  # 2.7726
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA = 0.1

# Physical constants
A_BOHR = 0.529177  # Å
E_RYDBERG = 13.6  # eV

# Stable orbits
D1 = 1.3333

def K(d):
    """K(d) kernel"""
    if d < 0.1:
        d = 0.1
    return ALPHA * np.cos(OMEGA * d + PHI) / (1 + BETA * d)

def octave_to_space(d):
    """Map octave d to spatial distance r in Å (from QW-720)"""
    return d * A_BOHR / D1

def quantum_to_octave(n):
    """Map quantum number n to octave d"""
    return D1 + 8 * (n - 1)

def quantum_to_space(n):
    """Map quantum number n to spatial distance r in Å (Bohr model)"""
    return n**2 * A_BOHR

# --- Build Hydrogen Hamiltonian with proper spatial mapping ---
print("\n[1] BUDOWANIE HAMILTONIANU WODORU")
print("-" * 60)

N_RADIAL = 50  # Radial grid points
r_max = 20.0 * A_BOHR  # Maximum radius (20×a_B)

# Radial grid in physical space (Å)
r_grid = np.linspace(0.1 * A_BOHR, r_max, N_RADIAL)
dr = r_grid[1] - r_grid[0]

print(f"  Siatka radialna: {N_RADIAL} punktów")
print(f"  r_min = {r_grid[0]:.4f} Å")
print(f"  r_max = {r_grid[-1]:.4f} Å")
print(f"  dr = {dr:.4f} Å")

# Build Hamiltonian
H = np.zeros((N_RADIAL, N_RADIAL))

# Kinetic energy (finite difference Laplacian in 3D)
# T = -ħ²/(2m) ∇² → -1/(2m) d²/dr² (radial, atomic units)
for i in range(N_RADIAL):
    if i > 0:
        H[i, i-1] = -1.0 / (2 * dr**2)
    H[i, i] = 1.0 / dr**2
    if i < N_RADIAL - 1:
        H[i, i+1] = -1.0 / (2 * dr**2)

# Potential energy: V(r) = -Z/r (Coulomb)
# But modified by K(d) coupling from octave structure
Z = 1.0  # Hydrogen

for i in range(N_RADIAL):
    r = r_grid[i]
    
    # Map r → n (quantum number) by inverting r = n² × a_B
    n_approx = np.sqrt(r / A_BOHR)
    
    # Map n → d (octave)
    d_approx = quantum_to_octave(n_approx)
    
    # Get K(d) coupling
    K_d = K(d_approx)
    
    # Coulomb potential modified by K(d)
    # V(r) = -Z × K(d) / r
    # But for hydrogen, we want standard Coulomb, so use K(d) as correction factor
    V_coulomb = -Z / r
    
    # Add K(d) modulation (small correction)
    # If K(d) > 0: stronger binding, if K(d) < 0: weaker binding
    K_factor = 1.0 + 0.1 * (K_d / ALPHA)  # Small modulation (10%)
    V_effective = V_coulomb * K_factor
    
    H[i, i] += V_effective

print(f"  H[0,0] = {H[0,0]:.4f}")
print(f"  H[0,1] = {H[0,1]:.4f}")

# --- Diagonalize ---
print("\n[2] DIAGONALIZACJA")
print("-" * 60)

eigenvalues, eigenvectors = eigh(H)
eigenvalues = np.sort(eigenvalues)

# Filter bound states (negative eigenvalues)
bound = eigenvalues[eigenvalues < 0]
print(f"  Stany związane (E < 0): {len(bound)}")

if len(bound) == 0:
    print("  ⚠️ Brak stanów związanych!")
    bound = eigenvalues[:5]

# --- Compare to Rydberg Series ---
print("\n[3] PORÓWNANIE Z SERIĄ RYDBERGA")
print("-" * 60)

def rydberg(n):
    """Rydberg formula: E_n = -E_Rydberg / n²"""
    return -E_RYDBERG / (n**2)

N_compare = min(len(bound), 5)
E_model = bound[:N_compare]

# Normalize to E_1 = -E_Rydberg
E_model_norm = E_model / abs(E_model[0]) * (-E_RYDBERG)

E_rydberg = [rydberg(n+1) for n in range(N_compare)]
E_rydberg_norm = np.array(E_rydberg)

print("| n | E_model (eV) | E_Rydberg (eV) | Błąd (%) |")
print("|---|--------------|----------------|----------|")

errors = []
for n in range(N_compare):
    model_val = E_model_norm[n]
    rydberg_val = E_rydberg_norm[n]
    err = abs(model_val - rydberg_val) / abs(rydberg_val) * 100
    errors.append(err)
    print(f"| {n+1} | {model_val:.4f} | {rydberg_val:.4f} | {err:.1f}% |")

mean_error = np.mean(errors)
print(f"\nŚredni błąd: {mean_error:.1f}%")

# --- Comparison with Prior Results ---
print("\n[4] PORÓWNANIE Z POPRZEDNIMI BADANIAMI")
print("-" * 60)
print("| Test | Metoda | Błąd |")
print("|------|--------|------|")
print("| QW-221 | Bezpośrednie oktawy | 250% |")
print("| QW-699 | Tylko K(d) | 110.5% |")
print("| QW-700 | Pełny nadsoliton | 50-110% |")
print(f"| **QW-721** | **Mapowanie oktaw→przestrzeń** | **{mean_error:.1f}%** |")

# --- Verdict ---
print("\n[5] WERDYKT")
print("=" * 60)

if mean_error < 10:
    verdict = "✅ SUKCES: Widmo wodoru odtworzone!"
    conclusion = "Mapowanie oktaw → przestrzeń działa doskonale."
elif mean_error < 20:
    verdict = "🟡 DUŻA POPRAWA"
    conclusion = f"Błąd spadł z 110% do {mean_error:.1f}%. Mapowanie działa dobrze."
elif mean_error < 50:
    verdict = "🟡 CZĘŚCIOWA POPRAWA"
    conclusion = f"Błąd spadł z 110% do {mean_error:.1f}%. Wymaga dalszych badań."
elif mean_error < 110:
    verdict = "🟡 NIEZNACZNA POPRAWA"
    conclusion = f"Błąd {mean_error:.1f}% vs 110.5% (QW-699). Minimalna poprawa."
else:
    verdict = "❌ BRAK POPRAWY"
    conclusion = "Mapowanie nie pomogło."

print(f"\n{verdict}")
print(f"Średni błąd: {mean_error:.1f}%")
print(f"Poprzedni najlepszy (QW-699): 110.5%")
print(f"Poprawa: {110.5 - mean_error:.1f}%")
print(f"\n{conclusion}")

# --- Analysis of Errors ---
print("\n[6] ANALIZA BŁĘDÓW")
print("-" * 60)

print("Błąd dla każdego poziomu:")
for n, err in enumerate(errors):
    if err < 10:
        status = "✅ Doskonały"
    elif err < 20:
        status = "🟡 Dobry"
    elif err < 50:
        status = "🟠 Akceptowalny"
    else:
        status = "❌ Wymaga poprawy"
    print(f"  n={n+1}: {err:.1f}% - {status}")

# --- Report ---
report_file = "raport_qw721_hydrogen_spectrum_mapped.md"
print(f"\n[7] ZAPISYWANIE RAPORTU: {report_file}")

with open(report_file, "w") as f:
    f.write("# RAPORT QW-721: WIDMO WODORU Z MAPOWANIEM OKTAW → PRZESTRZEŃ\n")
    f.write(f"**Data:** {datetime.datetime.now()}\n\n")
    
    f.write("## 1. Metoda\n")
    f.write("Użyto mapowania z QW-720:\n")
    f.write("- r(d) = d × a_B / d_0 (liniowe)\n")
    f.write("- n → d: d_n = d_1 + 8(n-1)\n")
    f.write("- n → r: r_n = n² × a_B (Bohr model)\n\n")
    
    f.write("## 2. Wyniki\n")
    f.write("| n | E_model (eV) | E_Rydberg (eV) | Błąd (%) |\n")
    f.write("|---|--------------|----------------|----------|\n")
    for n in range(N_compare):
        f.write(f"| {n+1} | {E_model_norm[n]:.4f} | {E_rydberg_norm[n]:.4f} | {errors[n]:.1f}% |\n")
    f.write(f"\n**Średni błąd:** {mean_error:.1f}%\n\n")
    
    f.write("## 3. Porównanie z Poprzednimi Badaniami\n")
    f.write("| Test | Błąd |\n")
    f.write("|------|------|\n")
    f.write("| QW-221 | 250% |\n")
    f.write("| QW-699 | 110.5% |\n")
    f.write("| QW-700 | 50-110% |\n")
    f.write(f"| **QW-721** | **{mean_error:.1f}%** |\n\n")
    
    f.write("## 4. Wnioski\n")
    f.write(f"### {verdict}\n")
    f.write(f"{conclusion}\n\n")
    
    if mean_error < 20:
        f.write("**Mapowanie oktaw → przestrzeń jest skuteczne dla widma atomowego.**\n")
    else:
        f.write("**Mapowanie wymaga dalszych badań i ulepszeń.**\n")

print("Done.")
