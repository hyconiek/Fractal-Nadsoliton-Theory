#!/usr/bin/env python3
"""
QW-676_Hydrogen_Spectrum.py
Purpose: RIGOROUS TEST of Hydrogen spectrum in FIN Theory.
         Does the model produce E_n = -13.6/n² eV?

Previous Results:
- QW-221/515: 250% error in Balmer series

Method:
1. Model proton-electron system in octave network
2. Solve for bound state energies
3. Compare ratio E_2/E_1 with theoretical 0.25

Success Criterion: E_2/E_1 = 0.25 ± 0.05
"""

import numpy as np
from scipy.linalg import eigh
import datetime

# --- Constants ---
N_OCTAVES = 12
ALPHA_GEO = 4 * np.log(2)
BETA = 0.1
OMEGA = np.pi / 4
PHI = np.pi / 6

# Hydrogen binding energy
E_RYDBERG = 13.6  # eV

REPORT_FILE = "RAPORT_HYDROGEN_SPECTRUM.md"

print(f"QW-676: HYDROGEN SPECTRUM TEST - Output: {REPORT_FILE}")

def K(d):
    """Effective coupling kernel"""
    return ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + BETA * d)

def build_hydrogen_hamiltonian(n_radial=20, proton_oct=7, electron_oct=1, coupling_strength=5.0):
    """
    Build Hamiltonian for hydrogen-like system.
    Use radial grid with Coulomb-like potential modified by K(d).
    """
    n = n_radial
    
    # Radial grid
    r = np.linspace(0.5, 10.0, n)
    dr = r[1] - r[0]
    
    # Kinetic energy (finite difference Laplacian)
    T = np.zeros((n, n))
    for i in range(n):
        T[i, i] = 2 / dr**2
        if i > 0:
            T[i, i-1] = -1 / dr**2
        if i < n-1:
            T[i, i+1] = -1 / dr**2
    
    # Potential: V(r) = -Z/r modified by octave coupling K(|oct_e - oct_p|)
    d_oct = abs(electron_oct - proton_oct)
    K_factor = K(d_oct)
    
    V = np.diag(-coupling_strength * K_factor / r)
    
    H = T + V
    
    return H, r

# ===================================================================
# MAIN EXECUTION
# ===================================================================

with open(REPORT_FILE, "w") as f:
    f.write(f"# RAPORT: WIDMO WODORU (QW-676)\n")
    f.write(f"**Data:** {datetime.datetime.now()}\n")
    f.write("**Cel:** Czy FIN produkuje $E_n \\propto 1/n^2$?\n\n")

    # Build Hamiltonian
    print("Building hydrogen Hamiltonian...")
    H, r = build_hydrogen_hamiltonian(n_radial=50, coupling_strength=10.0)
    
    # Diagonalize
    eigenvalues, eigenvectors = eigh(H)
    
    # Find bound states (negative energy)
    bound_mask = eigenvalues < 0
    bound_energies = eigenvalues[bound_mask]
    
    f.write(f"## 1. SETUP\n")
    f.write(f"- Siatka radialna: 50 punktów\n")
    f.write(f"- Proton oktawa: 7\n")
    f.write(f"- Elektron oktawa: 1\n")
    f.write(f"- K(d=6) = {K(6):.4f}\n\n")
    
    f.write(f"## 2. STANY ZWIĄZANE\n")
    f.write(f"- Liczba stanów związanych: {len(bound_energies)}\n\n")
    
    if len(bound_energies) >= 3:
        E1 = bound_energies[0]
        E2 = bound_energies[1]
        E3 = bound_energies[2]
        
        f.write("| Stan n | Energia E_n | E_n/E_1 | Teoria (1/n²) |\n")
        f.write("|--------|-------------|---------|---------------|\n")
        f.write(f"| 1 | {E1:.4f} | 1.000 | 1.000 |\n")
        f.write(f"| 2 | {E2:.4f} | {E2/E1:.4f} | 0.250 |\n")
        f.write(f"| 3 | {E3:.4f} | {E3/E1:.4f} | 0.111 |\n\n")
        
        print(f"E1 = {E1:.4f}, E2 = {E2:.4f}, E3 = {E3:.4f}")
        print(f"Ratio E2/E1 = {E2/E1:.4f} (should be 0.25)")
        
        # Fit to 1/n² law
        n_levels = np.array([1, 2, 3])
        E_levels = np.array([E1, E2, E3])
        E_normalized = E_levels / E1
        E_theory = 1 / n_levels**2
        
        error = np.mean(np.abs(E_normalized - E_theory))
        
        f.write(f"## 3. ZGODNOŚĆ Z $1/n^2$\n")
        f.write(f"- Średni błąd: {error*100:.1f}%\n")
        f.write(f"- E_2/E_1: {E2/E1:.4f} (teoria: 0.25)\n")
        f.write(f"- E_3/E_1: {E3/E1:.4f} (teoria: 0.111)\n\n")
        
        if error < 0.1:
            result = "✅ **SUCCESS:** Widmo zgodne z $1/n^2$!"
        elif error < 0.3:
            result = f"🟡 **PARTIAL:** Przybliżona struktura powłokowa (błąd {error*100:.1f}%)"
        else:
            result = f"❌ **FAIL:** Brak struktury $1/n^2$ (błąd {error*100:.1f}%)"
            
    elif len(bound_energies) >= 1:
        E1 = bound_energies[0]
        f.write(f"- Stan podstawowy E_1 = {E1:.4f}\n")
        f.write(f"- Za mało stanów wzbudzonych do analizy spektralnej.\n\n")
        result = "🟡 **PARTIAL:** Tylko stan podstawowy znaleziony."
    else:
        f.write("**BRAK stanów związanych!**\n\n")
        result = "❌ **FAIL:** Brak stanów związanych."
    
    f.write(result + "\n\n")
    print(result)
    
    # ===================================================================
    # SUMMARY
    # ===================================================================
    f.write(f"## 4. PODSUMOWANIE\n\n")
    
    if len(bound_energies) >= 3:
        f.write("| Test | Wynik | Status |\n")
        f.write("|------|-------|--------|\n")
        f.write(f"| Stany związane | {len(bound_energies)} | ✅ |\n")
        f.write(f"| E_2/E_1 | {E2/E1:.4f} | {'✅' if abs(E2/E1 - 0.25) < 0.05 else '❌'} |\n")
        f.write(f"| Zgodność $1/n^2$ | {(1-error)*100:.1f}% | {'✅' if error < 0.1 else '❌'} |\n\n")

print(f"\nReport written to {REPORT_FILE}")
