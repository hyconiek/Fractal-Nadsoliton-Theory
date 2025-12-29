#!/usr/bin/env python3
"""
QW-1613: MASS GENERATIONS FROM T(21,3) DEFORMATION
===================================================
Status: WYSOKIE RYZYKO (Optional)

Cel: Odtworzenie stosunków mas leptonów z energii węzłów torusowych

WAŻNE: Używamy JEDNEJ stałej kalibracyjnej (m_e = 0.511 MeV)
Sprawdzamy stosunki mas, nie wartości absolutne
"""

import numpy as np
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

REPORT_FILE = "RAPORT_QW1613_MASS_GENERATIONS.md"

# =============================================================================
# STAŁE
# =============================================================================
ALPHA_GEO = 4 * np.log(2)
GAMMA_MASS = 1.52  # Wykładnik masowy z poprzednich badań

# Masy eksperymentalne [MeV]
M_ELECTRON = 0.511
M_MUON = 105.66
M_TAU = 1776.86

md = []
def log(msg=""):
    print(msg)
    md.append(msg)

log("=" * 80)
log("QW-1613: MASS GENERATIONS FROM T(21,3) DEFORMATION")
log("=" * 80)
log(f"Data: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log("")

# =============================================================================
# MODEL MAS Z WĘZŁÓW TORUSOWYCH
# =============================================================================
def knot_energy(p, q):
    """
    Energia węzła torusowego T(p,q)
    
    E ∝ p² + q² (liczba przecięć)
    
    Stabilność wymaga gcd(p,q) = 1
    """
    return p**2 + q**2

def knot_topological_charge(p, q):
    """
    Ładunek topologiczny Q = p + q
    """
    return p + q

def mass_from_octave(d, gamma=GAMMA_MASS, m_ref=M_ELECTRON, d_ref=6.0):
    """
    Masa z pozycji oktawowej d
    
    m(d) = m_ref * exp(-γ * (d - d_ref))
    
    gdzie d_ref = 6 dla elektronu
    """
    return m_ref * np.exp(-gamma * (d - d_ref))

# =============================================================================
# WĘZŁY DLA LEPTONÓW
# =============================================================================
log("[1] WĘZŁY TORUSOWE DLA LEPTONÓW")
log("-" * 60)

# Elektron: T(21, 3) z QW-1201
p_e, q_e = 21, 3
d_e = 6.0  # Pozycja oktawowa

# Mion: T(13, 8) lub T(8, 5)?
# Z QW-1201: d_muon = 3.5, Q = 14
p_mu, q_mu = 8, 5  # Q = 13, close to 14
d_mu = 3.5

# Tau: T(5, 3) lub T(3, 2)?
# d_tau ≈ 2.25
p_tau, q_tau = 5, 3  # Q = 8
d_tau = 2.25

log(f"Elektron: T({p_e},{q_e}), Q = {p_e+q_e}, d = {d_e}")
log(f"Mion:     T({p_mu},{q_mu}), Q = {p_mu+q_mu}, d = {d_mu}")
log(f"Tau:      T({p_tau},{q_tau}), Q = {p_tau+q_tau}, d = {d_tau}")

# =============================================================================
# OBLICZENIE MAS
# =============================================================================
log("")
log("[2] PRZEWIDYWANE MASY (KALIBRACJA: m_e)")
log("-" * 60)

# Kalibracja: m_e = 0.511 MeV przy d = 6.0
m_e_pred = M_ELECTRON  # Punkt kalibracyjny

# Przewidywane masy
m_mu_pred = mass_from_octave(d_mu)
m_tau_pred = mass_from_octave(d_tau)

log(f"m_e (kalibracja):  {m_e_pred:.4f} MeV")
log(f"m_μ (predykcja):   {m_mu_pred:.4f} MeV")
log(f"m_τ (predykcja):   {m_tau_pred:.4f} MeV")

# =============================================================================
# STOSUNKI MAS
# =============================================================================
log("")
log("[3] STOSUNKI MAS")
log("-" * 60)

# Przewidywane stosunki
ratio_mu_e_pred = m_mu_pred / m_e_pred
ratio_tau_e_pred = m_tau_pred / m_e_pred
ratio_tau_mu_pred = m_tau_pred / m_mu_pred

# Eksperymentalne stosunki
ratio_mu_e_exp = M_MUON / M_ELECTRON
ratio_tau_e_exp = M_TAU / M_ELECTRON
ratio_tau_mu_exp = M_TAU / M_MUON

log(f"{'Stosunek':<15} {'Predykcja':>12} {'Eksperyment':>12} {'Błąd %':>10}")
log("-" * 50)

err_mu_e = abs(ratio_mu_e_pred - ratio_mu_e_exp) / ratio_mu_e_exp * 100
err_tau_e = abs(ratio_tau_e_pred - ratio_tau_e_exp) / ratio_tau_e_exp * 100
err_tau_mu = abs(ratio_tau_mu_pred - ratio_tau_mu_exp) / ratio_tau_mu_exp * 100

log(f"{'m_μ/m_e':<15} {ratio_mu_e_pred:>12.2f} {ratio_mu_e_exp:>12.2f} {err_mu_e:>10.2f}%")
log(f"{'m_τ/m_e':<15} {ratio_tau_e_pred:>12.2f} {ratio_tau_e_exp:>12.2f} {err_tau_e:>10.2f}%")
log(f"{'m_τ/m_μ':<15} {ratio_tau_mu_pred:>12.2f} {ratio_tau_mu_exp:>12.2f} {err_tau_mu:>10.2f}%")

# =============================================================================
# WERDYKT
# =============================================================================
log("")
log("[4] WERDYKT KOŃCOWY")
log("=" * 60)

avg_error = (err_mu_e + err_tau_e + err_tau_mu) / 3

if avg_error < 10:
    log(f"✅ Średni błąd stosunków: {avg_error:.1f}%")
    log("   Dobra zgodność fenomenologiczna!")
    overall_status = "VERIFIED"
elif avg_error < 30:
    log(f"⚠️ Średni błąd stosunków: {avg_error:.1f}%")
    log("   Jakościowa zgodność, wymaga udoskonalenia modelu")
    overall_status = "PARTIAL"
else:
    log(f"❌ Średni błąd stosunków: {avg_error:.1f}%")
    log("   Model wymaga fundamentalnej rewizji")
    overall_status = "FAILED"

log("")
log("UWAGA: To jest 'phenomenological consistency check', nie dowód TOE.")
log(f"OVERALL STATUS: {overall_status}")

# =============================================================================
# GENEROWANIE RAPORTU
# =============================================================================
with open(REPORT_FILE, 'w', encoding='utf-8') as f:
    f.write("# QW-1613: Mass Generations from T(21,3) Deformation\n\n")
    f.write(f"**Data:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    f.write(f"**Status:** {overall_status}\n\n")
    
    f.write("## Model\n")
    f.write("m(d) = m_e × exp(-γ × (d - 6.0))\n\n")
    f.write(f"gdzie γ = {GAMMA_MASS}\n\n")
    
    f.write("## Wyniki\n\n")
    f.write("| Stosunek | Predykcja | Eksperyment | Błąd |\n")
    f.write("|----------|-----------|-------------|------|\n")
    f.write(f"| m_μ/m_e | {ratio_mu_e_pred:.2f} | {ratio_mu_e_exp:.2f} | {err_mu_e:.1f}% |\n")
    f.write(f"| m_τ/m_e | {ratio_tau_e_pred:.2f} | {ratio_tau_e_exp:.2f} | {err_tau_e:.1f}% |\n")
    f.write(f"| m_τ/m_μ | {ratio_tau_mu_pred:.2f} | {ratio_tau_mu_exp:.2f} | {err_tau_mu:.1f}% |\n\n")
    
    f.write("## Raw Log\n```\n" + '\n'.join(md) + "\n```\n")

print(f"\n✅ Report saved to: {REPORT_FILE}")
