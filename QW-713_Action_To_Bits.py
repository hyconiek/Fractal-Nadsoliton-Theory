#!/usr/bin/env python3
"""
QW-713: ACTION TO BITS (MASA Z CZYNNOŚCI LAGRANGE'A)
=====================================================
Cel: Wyprowadzić masę wprost z Lagrangianu L_ZTP, traktując
     Czynność (Action S) jako miarę intensywności obliczeniowej.

Hipoteza:
  Masa ∝ Action S (w jednostkach 4-bitowych)

Metoda:
1. Zdefiniować pełny Lagrangian L_ZTP dla pola skalarnego Ψ.
2. Znaleźć rozwiązania solitonowe (paczki falowe) dla e, μ, τ.
3. Obliczyć Czynność S = ∫ L dt dla jednego okresu.
4. Przeliczyć S na bity (dzieląc przez stałą sprzężenia?).
5. Sprawdzić czy S_τ / S_e ≈ 3477 itd.

Date: 2025-12-08
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
import datetime

print("="*80)
print("QW-713: ACTION TO BITS - MASA Z LAGRANGIANU")
print("="*80)

# ===========================================================================
# 1. PEŁNY LAGRANGIAN L_ZTP (MODEL UPROSZCZONY)
# ===========================================================================

print("""
LAGRANGIAN L_ZTP:
L = ∫ [ (∂Ψ/∂t)² - (∂Ψ/∂x)² - V(Ψ) - H_int ] dx

Gdzie:
  H_int = energia interakcji między oktawami (przez K(d))
  V(Ψ) = potencjał nieliniowy (nadsoliton)
""")

# Parametry
ALPHA = 4 * np.log(2)  # 2.7726
BETA = 0.01
N_OCTAVES = 12

def K_kernel(d):
    """Kernel oddziaływania między oktawami."""
    if d == 0: return ALPHA
    return ALPHA * np.cos(np.pi/4 * d + np.pi/6) / (1 + BETA * d)

# Modelujemy cząstkę jako Gaussa w przestrzeni oktaw
def psi_packet(x, center, width, amplitude):
    """Profil pola w przestrzeni oktaw (x = numer oktawy)."""
    return amplitude * np.exp(-0.5 * ((x - center) / width)**2)

# Energia potencjalna (z K(d) i self-interaction)
def potential_energy(center, width, amplitude):
    # 1. Self-energy (lokalna nieliniowość V(Ψ) ~ Ψ⁴)
    # Całka z Ψ⁴
    norm_factor = np.sqrt(np.pi) * width
    E_self = 0.5 * amplitude**4 * norm_factor  # przykładowy człon V
    
    # 2. Interaction energy (K(d) coupling)
    # E_int = ∫∫ Ψ(x) K(|x-y|) Ψ(y) dx dy
    # Uproszczenie: dla wąskich paczek to głównie K(0) + poprawki
    
    # Dyskretyzacja na oktawy
    E_int = 0
    x_range = np.arange(N_OCTAVES)
    psi_vals = psi_packet(x_range, center, width, amplitude)
    
    for i in range(N_OCTAVES):
        for j in range(N_OCTAVES):
            E_int += psi_vals[i] * K_kernel(abs(i-j)) * psi_vals[j]
            
    return E_self - 0.5 * E_int  # Minus bo wiązanie obniża energię

# Energia kinetyczna (gradient w przestrzeni oktaw)
def kinetic_energy(center, width, amplitude):
    # E_kin = ∫ (∂Ψ/∂x)² dx
    # Dla Gaussa: pochodna to wielomian × Gauss
    # Całka z (dΨ/dx)² ≈ amplitude² / width
    return 0.5 * amplitude**2 / width

# Całkowita energia (Hamiltonian)
def total_energy(center, width, amplitude):
    return kinetic_energy(center, width, amplitude) + potential_energy(center, width, amplitude)

# ===========================================================================
# 2. ZNAJDOWANIE "MASY" JAKO ENERGII ROZWIĄZANIA STABILNEGO
# ===========================================================================

print("\nSzukanie stabilnych paczek (minimów energii) dla e, μ, τ...")

# Cząstki są zdefiniowane przez CENTRUM (oktawę).
# Szukamy optymalnej SZEROKOŚCI i AMPLITUDY dla każdego centrum.

particles = {
    'electron': {'center': 1.33},  # Orbit z QW-703
    'muon': {'center': 9.33},
    'tau': {'center': 17.33},      # Tau (poza zakresem 12? Traktujemy cyklicznie lub rozszerzamy)
    'proton': {'center': 13.46}
}

# Rozszerzamy zakres oktaw dla obliczeń (żeby tau się zmieściło)
N_CALC = 20

def solve_particle(center):
    """Znajdź stabilną konfigurację (min E) dla danego centrum."""
    
    def objective(params):
        width, amp = params
        if width < 0.1 or amp < 0.1: return 1e9
        
        # Oblicz energię
        # Uwaga: używamy rozszerzonego K_kernel
        E_int = 0
        x_vals = np.arange(N_CALC)
        psi = psi_packet(x_vals, center, width, amp)
        
        # Self-energy
        E_self = 0.5 * np.sum(psi**4)
        
        # Interaction (K(d))
        for i in range(N_CALC):
            for j in range(N_CALC):
                # Efektywny K(d) zanika
                if abs(i-j) < 10:
                    E_int += psi[i] * K_kernel(abs(i-j)) * psi[j]
        
        # Kinetic (gradient)
        dx_psi = np.gradient(psi)
        E_kin = 0.5 * np.sum(dx_psi**2)
        
        # Frequency term (Ψ ~ e^(-iωt))
        # E_total = E_kin + E_self - E_int
        # Lagranżian L = T - V
        # Dla stacjonarnego solitonu: E ≈ H = T + V
        
        return E_kin + E_self - 0.5 * E_int

    # Prosta optymalizacja (grid search dla stabilności)
    best_E = 1e9
    best_params = (0,0)
    
    for w in np.linspace(0.5, 3.0, 10):
        for a in np.linspace(0.5, 3.0, 10):
            E = objective((w, a))
            if E < best_E:
                best_E = E
                best_params = (w, a)
                
    return best_E, best_params

print(f"{'Cząstka':<10} {'Centrum':>8} {'Szerokość':>10} {'Amplituda':>10} {'ENERGIA':>10}")
print("-"*60)

energies = {}
for name, p in particles.items():
    E, (w, a) = solve_particle(p['center'])
    energies[name] = E
    print(f"{name:<10} {p['center']:>8.2f} {w:>10.2f} {a:>10.2f} {E:>10.4f}")

# ===========================================================================
# 3. INTERPRETACJA WYNIKÓW: CZY HIERARCHIA WYNIKA Z ENERGII?
# ===========================================================================

print("\n" + "="*80)
print("CZY ENERGIA POLA DAJE HIERARCHIĘ MAS?")
print("="*80)

E_e = energies['electron']

print(f"{'Cząstka':<10} {'Energia E':>10} {'Stosunek E/E_e':>15} {'Cel (Exp)':>10}")
print("-"*60)

for name, E in energies.items():
    ratio = E / E_e
    
    if name == 'electron': target = 1.0
    elif name == 'muon': target = 206.77
    elif name == 'tau': target = 3477.22
    elif name == 'proton': target = 1836.14
    
    print(f"{name:<10} {E:>10.4f} {ratio:>15.4f} {target:>10.1f}")

# ===========================================================================
# 4. PRZELICZENIE NA BITY (HIPOTEZA QW-711)
# ===========================================================================

print("\n" + "="*80)
print("PRZELICZENIE LAGRANGIANU NA BITY")
print("="*80)

print("""
HIPOTEZA:
Lagranżian (Czynność) to nie energia w dżulach, ale INTENSYWNOŚĆ w bitach.
Jeśli α = 4 bity, to czynność S powinna skalować się jak 2^(S/4)?

NIE. Jeśli L to intensywność (bity/czas), to:
Masa = całka z intensywności? Nie, masa to sama intensywność.

Sprawdźmy skalowanie EKSPONENCJALNE:
Czy Masa ∝ exp(Action)?

M = M_0 * exp(E / E_scale)
""")

# Sprawdźmy czy log(M) koreluje z E
print("Test: czy log(Masa) jest liniowe względem Energii Pola?")

import math
log_masses = {
    'electron': 0,
    'muon': np.log(206.77),
    'tau': np.log(3477.22),
    'proton': np.log(1836.14)
}

print(f"{'Cząstka':<10} {'Energia E':>10} {'log(M/M_e)':>12} {'Korelacja?':>12}")
print("-"*60)
for name, E in energies.items():
    lm = log_masses[name]
    print(f"{name:<10} {E:>10.4f} {lm:>12.4f}")

# Prosty fit: log(M) = a * E + b
# Sprawdzamy electron i muon
E_e = energies['electron']
E_m = energies['muon']
lm_e = 0
lm_m = np.log(206.77)

# a = (lm_m - lm_e) / (E_m - E_e)
if E_m != E_e:
    slope = (lm_m - lm_e) / (E_m - E_e)
    print(f"\nNachylenie (Slope) z e-μ: {slope:.4f}")
    
    # Przewidywanie dla tau i protona
    for name in ['tau', 'proton']:
        pred_log_m = slope * (energies[name] - E_e)
        pred_m_ratio = np.exp(pred_log_m)
        
        if name == 'tau': target = 3477.22
        else: target = 1836.14
        
        print(f"Predykcja {name}: ratio = {pred_m_ratio:.2f} (Cel: {target:.1f})")
else:
    print("\nBrak różnicy energii e-μ (model nie rozróżnia orbit?)")
    
# Save raw energies for analysis
with open("raport_qw713_action_to_bits.md", "w") as f:
    f.write("# RAPORT QW-713: Masa z Lagrangianu (Energia Pola)\n")
    f.write(f"**Date:** {datetime.datetime.now()}\n\n")
    
    f.write("## Wyniki Energii Pola (rozwiązania solitonowe)\n")
    f.write("| Cząstka | Centrum | Energia E |\n")
    f.write("|---------|---------|-----------|\n")
    for name, E in energies.items():
        f.write(f"| {name} | {particles[name]['center']} | {E:.4f} |\n")
        
    f.write("\n## Hipoteza Eksponencjalna: M ~ exp(E)\n")
    if E_m != E_e:
        f.write(f"Nachylenie a = {slope:.4f}\n")
        
print("\nReport saved to raport_qw713_action_to_bits.md")
print("="*80)
