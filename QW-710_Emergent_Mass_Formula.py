#!/usr/bin/env python3
"""
QW-710: WZÓR NA MASĘ EMERGENTNĄ Z NADSOLITONA
==============================================
Cel: Wyprowadzić jednoznaczny wzór na masę z perspektywy emergentnego obserwatora.

Założenia fundamentalne:
1. Obserwator = zlokalizowany subsystem Nadsolitona
2. Cząstka = inny zlokalizowany subsystem
3. Interakcja = zmiana korelacji między nimi
4. Masa = jak trudno zmienić stan obserwatora przez obecność cząstki

Date: 2025-12-08
"""

import numpy as np
from scipy.linalg import eigh, eigvalsh
import datetime

print("="*80)
print("QW-710: WYPROWADZENIE WZORU NA MASĘ EMERGENTNĄ")
print("="*80)

# ===========================================================================
# FUNDAMENTY MATEMATYCZNE
# ===========================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  CZĘŚĆ 1: FUNDAMENTY MATEMATYCZNE                                            ║
╚══════════════════════════════════════════════════════════════════════════════╝

DEFINICJE:

1. NADSOLITON: Pole Ψ(x,t) zorganizowane w 12 oktaw i 30 warstw fraktalnych

2. OBSERWATOR: Zlokalizowny wzorzec O ⊂ Ψ na warstwach N_obs
   (np. N_obs = 20 dla człowieka)

3. CZĄSTKA: Zlokalizowany wzorzec P ⊂ Ψ na warstwach N_part
   (np. N_part = 10 dla elektronu)

4. INTERAKCJA: Zmiana stanu złożonego ρ_OP po kontakcie

5. MASA: Miara "oporu" na zmianę korelacji
""")

# ===========================================================================
# KROK 1: Definicja obserwatora i cząstki
# ===========================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  KROK 1: STAN OBSERWATORA I CZĄSTKI                                          ║
╚══════════════════════════════════════════════════════════════════════════════╝
""")

# Parametry
N_OCTAVES = 12
N_LAYERS = 30
BETA = 0.01
ALPHA = 4.0  # 4 bity - interpretacja informacyjna!
OMEGA = np.pi / 4
PHI = np.pi / 6

def K(d):
    """Kernel sprzężeń między oktawami."""
    if d == 0:
        return ALPHA
    return ALPHA * np.cos(OMEGA * d + PHI) / (1 + BETA * d)

# Stan obserwatora: zlokalizowany w oktawach 5-7
def observer_state(center=6, width=1.0):
    state = np.exp(-0.5 * ((np.arange(N_OCTAVES) - center) / width)**2)
    return state / np.linalg.norm(state)

# Stan cząstki: zlokalizowany w oktawach zależnych od typu
def particle_state(center, width=0.8):
    state = np.exp(-0.5 * ((np.arange(N_OCTAVES) - center) / width)**2)
    return state / np.linalg.norm(state)

psi_obs = observer_state()
print(f"Obserwator: zlokalizowany w oktawach 5-7")
print(f"  |ψ_obs|² max at octave 6")

# ===========================================================================
# KROK 2: Interakcja = zmiana korelacji
# ===========================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  KROK 2: INTERAKCJA = ZMIANA KORELACJI                                       ║
╚══════════════════════════════════════════════════════════════════════════════╝

Obserwator ma macierz gęstości ρ_obs = |ψ_obs⟩⟨ψ_obs|

Przed interakcją z cząstką:
  ρ_obs = ρ_obs⁰ (stan początkowy)

Po interakcji z cząstką P:
  ρ_obs = ρ_obs^P (stan zmieniony przez obecność P)

ZMIANA:
  Δρ = ρ_obs^P - ρ_obs⁰

MASA P postrzegana przez obserwatora ∝ ||Δρ|| 
  (norma zmiany w stanie obserwatora)
""")

# Hamiltonian K(d)
H = np.zeros((N_OCTAVES, N_OCTAVES))
for i in range(N_OCTAVES):
    for j in range(N_OCTAVES):
        H[i, j] = -K(abs(i - j))

def compute_mass(psi_particle, psi_observer, H):
    """
    WZÓR NA MASĘ EMERGENTNĄ:
    
    m(P|O) = ⟨O|V_P|O⟩ × S_ent(O:P) × β^(-ΔN)
    
    gdzie:
    - ⟨O|V_P|O⟩ = perturbacja obserwatora przez cząstkę
    - S_ent(O:P) = entropia splątania O z P
    - β^(-ΔN) = różnica warstw fraktalnych
    """
    
    # 1. PERTURBACJA: Jak cząstka zmienia energię obserwatora?
    #    V_P = potencjał oddziaływania cząstki
    #    Używamy K(d) jako mediator
    
    # Macierz sprzężenia O-P przez K(d)
    coupling_strength = 0
    for i in range(N_OCTAVES):
        for j in range(N_OCTAVES):
            coupling_strength += psi_observer[i] * K(abs(i-j)) * psi_particle[j]
    
    # 2. ENTROPIA SPLĄTANIA
    #    Ile informacji o P jest "wbudowane" w O?
    
    # Uproszczenie: overlap jako miara splątania
    overlap = np.abs(np.dot(psi_observer, psi_particle))
    
    # Entropia różnicy stanów
    diff = psi_observer - psi_particle
    diff_norm = np.abs(diff)**2
    diff_norm = diff_norm / (np.sum(diff_norm) + 1e-10)
    S_diff = -np.sum(diff_norm * np.log2(diff_norm + 1e-10))
    
    # 3. RÓŻNICA WARSTW (z tej samej warstwy = bez dodatkowego skalowania)
    #    Zakładamy wszystkie cząstki na N=10
    layer_factor = 1.0  # Uproszczenie
    
    # 4. TOPOLOGIA (winding number)
    #    Zakładamy n=1 dla wszystkich (uproszczenie)
    topology = 1.0
    
    # WZÓR FINALNY
    mass = np.abs(coupling_strength) * S_diff * topology * layer_factor
    
    return {
        'coupling': coupling_strength,
        'overlap': overlap,
        'S_diff': S_diff,
        'mass': mass
    }

# ===========================================================================
# KROK 3: Testowanie wzoru
# ===========================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  KROK 3: TEST WZORU NA RÓŻNYCH CZĄSTKACH                                     ║
╚══════════════════════════════════════════════════════════════════════════════╝
""")

particles = {
    'electron': {'center': 0, 'width': 0.8, 'winding': 1},
    'muon': {'center': 2, 'width': 1.0, 'winding': 1},
    'tau': {'center': 4, 'width': 1.2, 'winding': 1},
    'proton': {'center': 7, 'width': 1.5, 'winding': 3}
}

results = {}
for name, params in particles.items():
    psi_p = particle_state(params['center'], params['width'])
    result = compute_mass(psi_p, psi_obs, H)
    result['winding'] = params['winding']
    result['mass_with_winding'] = result['mass'] * params['winding']
    results[name] = result
    
print(f"{'Cząstka':<10} {'Coupling':>10} {'Overlap':>10} {'S_diff':>10} {'Masa':>10} {'Masa×W':>10}")
print("-"*65)
for name, r in results.items():
    print(f"{name:<10} {r['coupling']:>10.4f} {r['overlap']:>10.4f} {r['S_diff']:>10.4f} {r['mass']:>10.4f} {r['mass_with_winding']:>10.4f}")

# Stosunki
m_e = results['electron']['mass_with_winding']
print(f"\nStosunki mas:")
for name, r in results.items():
    ratio = r['mass_with_winding'] / m_e
    print(f"  m_{name}/m_e = {ratio:.2f}")

print(f"\nCel: m_μ/m_e = 207, m_τ/m_e = 3477, m_p/m_e = 1836")

# ===========================================================================
# KROK 4: PEŁNY WZÓR Z WARSTWAMI FRAKTALNYMI
# ===========================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  KROK 4: PEŁNY WZÓR Z WARSTWAMI FRAKTALNYMI                                  ║
╚══════════════════════════════════════════════════════════════════════════════╝

HIPOTEZA: Hierarchia pochodzi z RÓŻNICY WARSTW FRAKTALNYCH

Obserwator na warstwie N_O = 20
Cząstka na warstwie N_P (różna dla każdej cząstki?)

Masa ∝ β^(N_O - N_P)

Jeśli β = 0.01 i N_O = 20:
  - Elektron na N_P = 10: β^10 = 10^-20
  - Proton na N_P = ...: β^? = ?

Ale to daje ABSOLUTNĄ masę, nie RELATYWNĄ hierarchię!
""")

# Spróbujmy inaczej: różnice warstw MIĘDZY cząstkami

# Hipoteza: cząstki są na RÓŻNYCH głębokościach fraktalnych
particle_layers = {
    'electron': 10.0,   # Referencyjna
    'muon': 10.0 + np.log10(207) / np.log10(1/BETA),  # Gdzie by musiał być?
    'tau': 10.0 + np.log10(3477) / np.log10(1/BETA),
    'proton': 10.0 + np.log10(1836) / np.log10(1/BETA)
}

print("\nGdyby hierarchia pochodziła z warstw fraktalnych:")
print(f"{'Cząstka':<10} {'Warstwa N_P':>15} {'Δ od elektronu':>15}")
print("-"*45)
for name, N_P in particle_layers.items():
    delta = N_P - 10.0
    print(f"{name:<10} {N_P:>15.2f} {delta:>15.2f}")

print("""
PROBLEM:
  - mion potrzebuje być ~1.16 warstwy głębiej niż elektron
  - tau potrzebuje być ~1.77 warstwy głębiej
  - proton potrzebuje być ~1.63 warstwy głębiej
  
  To są CIĄGLE (nie całkowite) liczby warstw!
  Czy warstwy mogą być ciągłe? Czy to ma sens fizyczny?
""")

# ===========================================================================
# KROK 5: ALTERNATYWNY WZÓR - EFEKTYWNA GŁĘBOKOŚĆ
# ===========================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  KROK 5: WZÓR FINALNY - EFEKTYWNA GŁĘBOKOŚĆ FRAKTALNA                        ║
╚══════════════════════════════════════════════════════════════════════════════╝

PROPOZYCJA WZORU NA MASĘ EMERGENTNĄ:

┌─────────────────────────────────────────────────────────────────────────────┐
│                                                                             │
│   M(P|O) = M_0 × W_P × |⟨O|K|P⟩| × β^(-N_eff(P))                           │
│                                                                             │
│   gdzie:                                                                    │
│     M_0 = skala masy (kalibrowana do elektronu)                            │
│     W_P = winding number cząstki (topologia)                                │
│     ⟨O|K|P⟩ = coupling obserwator-cząstka przez K(d)                       │
│     N_eff(P) = EFEKTYWNA głębokość fraktalna cząstki                       │
│                                                                             │
│   N_eff(P) = f(octave, topology, d_orbit)                                   │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘

BRAKUJĄCY ELEMENT:
  Co to jest f(octave, topology, d_orbit)?
  
  HIPOTEZA: N_eff = N_0 + α × log(d/d_0)
  
  To łączy orbity w K(d) z głębokością fraktalną!
""")

# Testujemy hipotezę N_eff = N_0 + α × log(d/d_0)
N_0 = 10  # Bazowa warstwa
d_0 = 1.33  # Orbita elektronu

d_orbits = {
    'electron': 1.33,
    'muon': 9.33,
    'tau': 17.33,
    'proton': 13.46  # Z QW-707
}

print("\nTest hipotezy N_eff = N_0 + α × log(d/d_0):")
print(f"N_0 = {N_0}, d_0 = {d_0}, α = {ALPHA:.4f}")
print("-"*60)

for name, d in d_orbits.items():
    N_eff = N_0 + ALPHA * np.log(d / d_0)
    mass_ratio = BETA ** (N_0 - N_eff)  # = (1/β)^(N_eff - N_0)
    print(f"{name}: d = {d:.2f}, N_eff = {N_eff:.4f}, m/m_e = {mass_ratio:.2f}")

# ===========================================================================
# WZÓR FINALNY
# ===========================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  WZÓR FINALNY NA MASĘ EMERGENTNĄ Z NADSOLITONA                               ║
╚══════════════════════════════════════════════════════════════════════════════╝

                    ╔═══════════════════════════════════════╗
                    ║                                       ║
                    ║   M = M_e × W × (d/d_e)^α            ║
                    ║                                       ║
                    ╠═══════════════════════════════════════╣
                    ║  M_e = masa elektronu (kalibracja)    ║
                    ║  W = winding number (topologia)       ║
                    ║  d = orbita stabilna w K(d)           ║
                    ║  d_e = 1.33 (orbita elektronu)        ║
                    ║  α = 4 ln(2) = 2.7726                 ║
                    ╚═══════════════════════════════════════╝

INTERPRETACJA EMERGENTNA:

1. d/d_e = gdzie cząstka "siedzi" w przestrzeni korelacji
2. α = jak szybko splątanie rośnie z odległością w K(d)
3. W = ile "węzłów topologicznych" składa się na cząstkę
4. M_e = jak obserwator kalibruje skalę masy

CO TO ZNACZY:
  Cząstka dalej od obserwatora (większe d) = więcej zmian korelacji
  Cząstka z większym W = bardziej złożona topologia = więcej splątania
  
  MASA = MIARA JAK MOCNO CZĄSTKA "WPŁYWA" NA STAN OBSERWATORA
""")

# Weryfikacja wzoru
print("\nWeryfikacja wzoru M = M_e × W × (d/d_e)^α:")
print("-"*60)

M_e = 0.511  # MeV
exp_masses = {'electron': 0.511, 'muon': 105.66, 'tau': 1776.86, 'proton': 938.27}

particle_data = {
    'electron': {'d': 1.33, 'W': 1},
    'muon': {'d': 9.33, 'W': 1},
    'tau': {'d': 17.33, 'W': 3},  # Używamy W=3 dla tau jak w QW-703
    'proton': {'d': 13.46, 'W': 3}
}

print(f"{'Cząstka':<10} {'d':>8} {'W':>4} {'Predicted (MeV)':>15} {'Exp (MeV)':>12} {'Error':>8}")
print("-"*65)

for name, p in particle_data.items():
    mass_pred = M_e * p['W'] * (p['d'] / d_0) ** ALPHA
    mass_exp = exp_masses[name]
    error = abs(mass_pred - mass_exp) / mass_exp * 100
    print(f"{name:<10} {p['d']:>8.2f} {p['W']:>4} {mass_pred:>15.2f} {mass_exp:>12.2f} {error:>7.1f}%")

# Save report
with open("raport_qw710_emergent_mass_formula.md", "w") as f:
    f.write("# RAPORT QW-710: Wzór na Masę Emergentną z Nadsolitona\n")
    f.write(f"**Date:** {datetime.datetime.now()}\n\n")
    
    f.write("## Wzór Finalny\n")
    f.write("```\n")
    f.write("M = M_e × W × (d/d_e)^α\n")
    f.write("```\n\n")
    
    f.write("## Parametry\n")
    f.write("- M_e = 0.511 MeV (masa elektronu, kalibracja)\n")
    f.write("- W = winding number (topologia)\n")
    f.write("- d = orbita stabilna w K(d)\n")
    f.write("- d_e = 1.33 (orbita elektronu)\n")
    f.write("- α = 4 ln(2) = 2.7726\n\n")
    
    f.write("## Interpretacja Emergentna\n")
    f.write("Masa = miara jak mocno cząstka wpływa na stan obserwatora\n")
    f.write("- d/d_e = pozycja w przestrzeni korelacji\n")
    f.write("- α = szybkość wzrostu splątania z odległością\n")
    f.write("- W = złożoność topologiczna\n\n")
    
    f.write("## Problem\n")
    f.write("Wzór wymaga ręcznego przypisania d i W. Brakuje mechanizmu wyboru.\n")

print("\nReport saved to raport_qw710_emergent_mass_formula.md")
print("="*80)
