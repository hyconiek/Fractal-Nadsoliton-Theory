#!/usr/bin/env python3
"""
QW-1614: G-FACTOR FROM KNOT GEOMETRY
====================================
Test BARDZO MOCNY - emergencja Diraca z geometrii

Cel: g ≈ 2 + O(10⁻³) z czystej geometrii splotu

Metodologia:
1. Oblicz moment magnetyczny μ z pętli prądowej na splocie
2. Oblicz moment pędu L klasycznie
3. g = 2μ/L
4. Dodaj korekty krzywizną
"""

import numpy as np
from scipy.integrate import quad, dblquad
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

REPORT_FILE = "RAPORT_QW1614_G_FACTOR.md"

# =============================================================================
# STAŁE
# =============================================================================
ALPHA_GEO = 4 * np.log(2)
BETA_TORS = 0.01

md = []
def log(msg=""):
    print(msg)
    md.append(msg)

log("=" * 80)
log("QW-1614: G-FACTOR FROM KNOT GEOMETRY")
log("=" * 80)
log(f"Data: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log("")

# =============================================================================
# GEOMETRIA WĘZŁA TORUSOWEGO T(p,q)
# =============================================================================
def torus_knot_parametrization(t, p, q, R=2.0, r=1.0):
    """
    Parametryzacja węzła torusowego T(p,q)
    
    x(t) = (R + r*cos(q*t)) * cos(p*t)
    y(t) = (R + r*cos(q*t)) * sin(p*t)
    z(t) = r * sin(q*t)
    
    t ∈ [0, 2π]
    """
    x = (R + r * np.cos(q * t)) * np.cos(p * t)
    y = (R + r * np.cos(q * t)) * np.sin(p * t)
    z = r * np.sin(q * t)
    
    return np.array([x, y, z])

def torus_knot_tangent(t, p, q, R=2.0, r=1.0):
    """
    Wektor styczny do węzła T(p,q)
    """
    dx_dt = -p * (R + r*np.cos(q*t)) * np.sin(p*t) - q*r*np.sin(q*t)*np.cos(p*t)
    dy_dt = p * (R + r*np.cos(q*t)) * np.cos(p*t) - q*r*np.sin(q*t)*np.sin(p*t)
    dz_dt = q * r * np.cos(q*t)
    
    return np.array([dx_dt, dy_dt, dz_dt])

def compute_writhe(p, q, n_samples=500):
    """
    Oblicza writhe węzła (miara nie-płaskości)
    
    Wr = (1/4π) ∫∫ (r1 - r2) · (dr1 × dr2) / |r1 - r2|³
    """
    t_vals = np.linspace(0, 2*np.pi, n_samples)
    dt = t_vals[1] - t_vals[0]
    
    writhe = 0.0
    
    for i in range(n_samples):
        for j in range(n_samples):
            if i == j:
                continue
            
            t1, t2 = t_vals[i], t_vals[j]
            r1 = torus_knot_parametrization(t1, p, q)
            r2 = torus_knot_parametrization(t2, p, q)
            
            dr1 = torus_knot_tangent(t1, p, q) * dt
            dr2 = torus_knot_tangent(t2, p, q) * dt
            
            diff = r1 - r2
            dist = np.linalg.norm(diff) + 1e-10
            
            cross = np.cross(dr1, dr2)
            writhe += np.dot(diff, cross) / dist**3
    
    return writhe / (4 * np.pi)

# =============================================================================
# MOMENT MAGNETYCZNY Z PĘTLI PRĄDOWEJ
# =============================================================================
def compute_magnetic_moment(p, q, I=1.0, n_samples=200):
    """
    Moment magnetyczny pętli prądowej na węźle
    
    μ = (I/2c) ∮ r × dl
    
    Dla jednostek naturalnych c = 1
    """
    t_vals = np.linspace(0, 2*np.pi, n_samples)
    dt = t_vals[1] - t_vals[0]
    
    mu = np.zeros(3)
    
    for t in t_vals:
        r = torus_knot_parametrization(t, p, q)
        dr = torus_knot_tangent(t, p, q) * dt
        
        mu += 0.5 * I * np.cross(r, dr)
    
    return mu

def compute_angular_momentum(p, q, m=1.0, omega=1.0, n_samples=200):
    """
    Moment pędu dla orbity na węźle
    
    L = m ∮ r × v dl
    
    Zakładamy v ∝ tangent
    """
    t_vals = np.linspace(0, 2*np.pi, n_samples)
    dt = t_vals[1] - t_vals[0]
    
    L = np.zeros(3)
    
    for t in t_vals:
        r = torus_knot_parametrization(t, p, q)
        v = torus_knot_tangent(t, p, q) * omega
        
        L += m * np.cross(r, v) * dt
    
    return L

# =============================================================================
# OBLICZENIE G-FACTOR
# =============================================================================
log("[1] ANALIZA WĘZŁA ELEKTRONU T(21,3)")
log("-" * 60)

# Elektron jako T(21, 3)
p_e, q_e = 21, 3

log(f"Węzeł elektronu: T({p_e}, {q_e})")
log(f"Q = p + q = {p_e + q_e}")

# Oblicz moment magnetyczny
mu = compute_magnetic_moment(p_e, q_e)
mu_mag = np.linalg.norm(mu)
log(f"Moment magnetyczny |μ| = {mu_mag:.6f}")

# Oblicz moment pędu
L = compute_angular_momentum(p_e, q_e)
L_mag = np.linalg.norm(L)
log(f"Moment pędu |L| = {L_mag:.6f}")

# G-factor: g = 2|μ|/|L| (dla jednostek gdzie e = ℏ = 1)
g_classical = 2 * mu_mag / L_mag if L_mag > 0 else 0
log(f"g-factor (klasyczny) = {g_classical:.6f}")

# =============================================================================
# KOREKTY KRZYWIZNĄ
# =============================================================================
log("")
log("[2] KOREKTY GEOMETRYCZNE")
log("-" * 60)

# Writhe węzła (przybliżenie)
# Dla T(p,q): Wr ≈ (p-1)(q-1) / 2 (dla gcd(p,q)=1)
writhe_approx = (p_e - 1) * (q_e - 1) / 2
log(f"Writhe (przybliżenie): Wr ≈ {writhe_approx:.1f}")

# Korekta z torsji β
torsion_correction = BETA_TORS * (1 + 1/p_e)
log(f"Korekta torsyjna: δg_tors = {torsion_correction:.6f}")

# Korekta z informacji α
info_correction = ALPHA_GEO / (2 * np.pi * L_mag) if L_mag > 0 else 0
log(f"Korekta informacyjna: δg_info = {info_correction:.6f}")

# =============================================================================
# PORÓWNANIE Z WARTOŚCIĄ DIRACOWĄ
# =============================================================================
log("")
log("[3] PORÓWNANIE Z TEORIĄ DIRACA")
log("-" * 60)

# Teoretyczna wartość Diracowa
g_dirac = 2.0

# Eksperymentalna wartość z QED (g-2)/2 = anomalous magnetic moment
g_exp = 2.00231930436256  # Najbardziej precyzyjnie zmierzona stała w fizyce

# Nasza wartość z poprawkami
g_total = g_classical + torsion_correction + info_correction

log(f"g (Dirac, teoria):    {g_dirac:.10f}")
log(f"g (QED, eksperyment): {g_exp:.10f}")
log(f"g (FIN, geometria):   {g_total:.10f}")

# Anomalny moment
a_fin = (g_total - 2) / 2
a_qed = (g_exp - 2) / 2

log(f"")
log(f"Anomalny moment (g-2)/2:")
log(f"  QED:  a = {a_qed:.10f}")
log(f"  FIN:  a = {a_fin:.10f}")

# =============================================================================
# WERDYKT
# =============================================================================
log("")
log("[4] WERDYKT KOŃCOWY")
log("=" * 60)

# Kryterium: g ≈ 2 + O(10⁻³)
g_error = abs(g_total - g_dirac)

if g_error < 0.01:
    log(f"✅ g = {g_total:.6f} ≈ 2 + O(10⁻³)")
    log("")
    log("SUKCES: g-factor emerguje z geometrii węzła!")
    log("        Bez wprowadzania spinu jako fundamentalnego.")
    overall_status = "VERIFIED"
elif g_error < 0.1:
    log(f"⚠️ g = {g_total:.6f}, błąd = {g_error:.4f}")
    log("   Zbliżone do 2, ale wymaga dokładniejszych poprawek.")
    overall_status = "PARTIAL"
else:
    log(f"❌ g = {g_total:.6f}, błąd = {g_error:.4f}")
    log("   Znaczne odchylenie od wartości Diracowej.")
    overall_status = "FAILED"

log("")
log(f"OVERALL STATUS: {overall_status}")

# =============================================================================
# GENEROWANIE RAPORTU
# =============================================================================
with open(REPORT_FILE, 'w', encoding='utf-8') as f:
    f.write("# QW-1614: g-factor from Knot Geometry\n\n")
    f.write(f"**Data:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    f.write(f"**Status:** {overall_status}\n\n")
    
    f.write("## Hipoteza\n")
    f.write("g-factor elektronu (= 2 w teorii Diraca) emerguje z geometrii\n")
    f.write("węzła torusowego T(21,3) bez wprowadzania spinu jako fundamentalnego.\n\n")
    
    f.write("## Metodologia\n")
    f.write("1. μ = (1/2) ∮ r × J dl (moment magnetyczny pętli)\n")
    f.write("2. L = m ∮ r × v dl (moment pędu)\n")
    f.write("3. g = 2|μ|/|L|\n")
    f.write("4. Korekty: torsja (β), informacja (α)\n\n")
    
    f.write("## Wyniki\n\n")
    f.write("| Wielkość | Wartość |\n")
    f.write("|----------|--------|\n")
    f.write(f"| |μ| | {mu_mag:.6f} |\n")
    f.write(f"| |L| | {L_mag:.6f} |\n")
    f.write(f"| g (klasyczny) | {g_classical:.6f} |\n")
    f.write(f"| g (z poprawkami) | {g_total:.6f} |\n")
    f.write(f"| g (Dirac) | {g_dirac:.6f} |\n")
    f.write(f"| Błąd | {g_error:.6f} |\n\n")
    
    f.write("## Werdykt\n")
    if overall_status == "VERIFIED":
        f.write("> **SUKCES:** g ≈ 2 emerguje z geometrii węzła.\n")
        f.write("> Potwierdza to emergentną naturę spinu w FIN.\n")
    
    f.write("\n## Raw Log\n```\n" + '\n'.join(md) + "\n```\n")

print(f"\n✅ Report saved to: {REPORT_FILE}")
