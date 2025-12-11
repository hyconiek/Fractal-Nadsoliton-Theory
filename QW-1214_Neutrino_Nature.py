import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# ==============================================================================
# QW-1214: NATURA NEUTRINA I FALE UDERZENIOWE (SHOCK WAVES)
# ==============================================================================
#
# CEL: Zbadać naturę neutrin w teorii FIN.
#      Czy są to bezmasowe fale (Q=0) czy masywne solitony?
#
# METODA: 
# 1. Analiza równania dyspersji dla Master Equation (QW-499) w granicy małych amplitud.
#    i dPsi/dt = - Laplacian Psi + V(Psi)
# 2. Poszukiwanie rozwiązań solitonowych typu "Breather" lub "Pulse" bez topologii węzła.
# 3. Sprawdzenie, czy takie fale mają masę spoczynkową (gap w widmie).
#
# ==============================================================================

REPORT_FILE = "RAPORT_QW1214_NEUTRINO_NATURE.md"
md = []

def log(msg):
    print(msg)
    md.append(msg)

log("=" * 80)
log("QW-1214: NATURA NEUTRINA - FALA CZY SOLITON?")
log("=" * 80)

def dispersion_relation(k):
    # Master Equation (uproszczone 1D):
    # i dt Psi = - d2 Psi + beta * Psi (beta - tłumienie lub masa efektywna)
    # W FIN, człon -beta_tors * Psi działa jak masa w równaniu Diraca?
    
    # Równanie (z dokumentacji):
    # dt psi = i(H0 + ...) - beta*psi ...
    # Ale dla swobodnej fali w próżni (z dala od węzłów):
    # E = sqrt(p^2 + m^2)
    # Skąd bierze się m dla neutrina?
    
    # Hipoteza: Neutrino to fala w ośrodku o lepkości beta.
    # omega(k) = c*k + i*beta (tłumiona)
    # LUB
    # omega(k) = sqrt(c^2 k^2 + omega_cut^2) (masa)
    
    # Sprawdźmy 'omega_cut' wynikającą z dyskretyzacji (siatki fraktalnej).
    # K(d) ma skalę 'omega = pi/4'. 
    # To narzuca gap energetyczny? E_min = h * omega_min.
    
    omega_fund = np.pi / 4.0
    return np.sqrt(k**2 + omega_fund**2)

log("\n[1] ANALIZA DYSPERSJI")
k_vals = np.linspace(0, 5, 20)
omega_vals = dispersion_relation(k_vals)

log(f"{'k':<10} {'omega':<10} {'v_group':<10}")
log("-" * 40)

for i in range(len(k_vals)):
    k = k_vals[i]
    w = omega_vals[i]
    # Prędkość grupowa dw/dk = k / sqrt(k^2 + m^2) -> k/w
    vg = k / w if w > 0 else 0
    
    log(f"{k:<10.2f} {w:<10.2f} {vg:<10.4f}")

log("\n[2] WNIOSKI O MASIE NEUTRINA")
log("-" * 80)
m_eff = np.pi / 4.0
log(f"Efektywna masa (Gap): m_eff = pi/4 = {m_eff:.4f} (jednostek naturalnych)")

# Przeliczenie na jednostki masy elektronu?
# Skala jest wyznaczona przez alpha_geo = 4ln2.
# Jeśli m_eff pochodzi z K(d), to neutrina POWINNY mieć masę.
# Ale neutrina są superlekkie (bezwzględnie).
# Dlaczego?

log("""
PARADOKS MASY:
Jeśli Nadsoliton ma strukturę dyskretną (omega = pi/4), to każda fala
powinna mieć masę (gap).
Dlaczego neutrina są prawie bezmasowe (m < 0.1 eV)?

HIPOTEZA ROZWIĄZANIA (CHIRAL WAVE):
Neutrina nie są wibracjami sieci (fononami), ale wibracjami SKRĘCENIA sieci (torsion waves).
Fale torsyjne w ciele stałym mogą być bezszczelinowe (gapless), jeśli symetria cechowania na to pozwala.

Wymagana symetria: Chiralna.
W FIN mamy beta_tors = 0.01. To małe łamanie symetrii chiralnej.
Masa neutrina ~ beta_tors * M_scale?
Jeśli M_scale ~ M_Planck, to za dużo.
Jeśli M_scale ~ M_electron...
m_nu ~ 0.01 * 0.5 MeV = 5 keV. (Wciąż za dużo, eksperyment m < 0.8 eV).

Wniosek: Neutrino musi być **Goldstone Boson** (lub Fermion)
związanym z łamaniem symetrii, która jest PRAWIE zachowana.
""")

# Zapis raportu
with open(REPORT_FILE, 'w') as f:
    f.write("# QW-1214: Natura Neutrina\n\n")
    f.write(f"**Data:** {datetime.now()}\n\n")
    f.write("```\n")
    f.write("\n".join(md))
    f.write("\n```\n")

print(f"Raport zapisano w {REPORT_FILE}")
