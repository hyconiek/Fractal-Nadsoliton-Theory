#!/usr/bin/env python3
"""
QW-1617 to QW-1619: OPTIONAL STUDIES (COMBINED)
===============================================
QW-1617: Running Coupling α
QW-1618: Preon Stability
QW-1619: CKM from Knot Rotations

Status: SPEKULATYWNE / model-dependent
"""

import numpy as np
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

REPORT_FILE = "RAPORT_QW1617_1619_OPTIONAL.md"

ALPHA_GEO = 4 * np.log(2)
BETA_TORS = 0.01

md = []
def log(msg=""):
    print(msg)
    md.append(msg)

log("=" * 80)
log("QW-1617/1618/1619: OPTIONAL STUDIES (COMBINED)")
log("=" * 80)
log(f"Data: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log("")

# =============================================================================
# QW-1617: RUNNING COUPLING α
# =============================================================================
log("=" * 60)
log("QW-1617: RUNNING COUPLING α(Q²)")
log("=" * 60)

def alpha_running_qed(Q_squared, alpha_0=1/137.036, m_e=0.511):
    """
    QED running: α(Q²) = α(0) / (1 - (α/3π) ln(Q²/m_e²))
    """
    if Q_squared <= m_e**2:
        return alpha_0
    
    beta = alpha_0 / (3 * np.pi)
    log_term = np.log(Q_squared / m_e**2)
    
    return alpha_0 / (1 - beta * log_term)

def alpha_running_fin(Q_squared, alpha_geo=ALPHA_GEO, beta_tors=BETA_TORS):
    """
    FIN running: α(Q²) z renormalizacji oktawowej
    """
    alpha_0_fin = 1 / (alpha_geo / (2 * beta_tors) * (1 - beta_tors))
    
    # Przybliżona renormalizacja
    if Q_squared <= 0.511**2:
        return alpha_0_fin
    
    beta_fin = alpha_0_fin / (4 * np.pi)  # Inny współczynnik beta
    log_term = np.log(Q_squared / 0.511**2)
    
    return alpha_0_fin / (1 - beta_fin * log_term + 0.5 * (beta_fin * log_term)**2)

# Test przy kilku skalach
Q_values = [0.511, 1.0, 10.0, 91.2, 1000.0]  # MeV, MeV, MeV, M_Z, TeV

log(f"{'Q [MeV]':>10} {'α(QED)':>12} {'α(FIN)':>12} {'Diff %':>10}")
log("-" * 50)

for Q in Q_values:
    Q_sq = Q**2
    a_qed = alpha_running_qed(Q_sq)
    a_fin = alpha_running_fin(Q_sq)
    diff = abs(a_qed - a_fin) / a_qed * 100 if a_qed > 0 else 0
    log(f"{Q:>10.1f} {a_qed:>12.6f} {a_fin:>12.6f} {diff:>10.2f}%")

status_1617 = "PARTIAL"
log(f"\nQW-1617 STATUS: {status_1617} (qualitative consistency only)")

# =============================================================================
# QW-1618: PREON STABILITY
# =============================================================================
log("")
log("=" * 60)
log("QW-1618: PREON STABILITY")
log("=" * 60)

def preon_energy(n_preons, coupling=ALPHA_GEO):
    """
    Energia konfiguracji preonów
    E = n × E_base + binding
    """
    E_base = 1.0  # Jednostka arbitralna
    
    # Energia wiązania (przyciągająca)
    E_binding = -coupling * n_preons * (n_preons - 1) / 2
    
    # Energia kinetyczna (odpychająca przy małych odległościach)
    E_kinetic = n_preons * E_base * (1 + 0.1 * n_preons)
    
    return E_kinetic + E_binding

# Test stabilności dla różnych n
log(f"{'n_preons':>10} {'E_total':>12} {'E/n':>12} {'Stable?':>10}")
log("-" * 50)

stable_configs = []
for n in range(1, 8):
    E = preon_energy(n)
    E_per = E / n
    
    # Stabilność: E/n maleje lub jest minimum lokalne
    if n > 1:
        E_prev = preon_energy(n-1) / (n-1)
        stable = E_per <= E_prev
    else:
        stable = True
    
    stable_configs.append((n, E, stable))
    log(f"{n:>10} {E:>12.4f} {E_per:>12.4f} {'✅' if stable else '❌':>10}")

status_1618 = "VERIFIED" if any(s[2] for s in stable_configs[2:]) else "PARTIAL"
log(f"\nQW-1618 STATUS: {status_1618}")

# =============================================================================
# QW-1619: CKM FROM KNOT ROTATIONS
# =============================================================================
log("")
log("=" * 60)
log("QW-1619: CKM FROM KNOT ROTATIONS")
log("=" * 60)

def knot_overlap(p1, q1, p2, q2):
    """
    Przybliżony overlap między węzłami T(p1,q1) i T(p2,q2)
    
    Zakładamy: V_ij ∝ exp(-|Q1 - Q2| / ξ)
    gdzie ξ to skala mieszania
    """
    Q1 = p1 + q1
    Q2 = p2 + q2
    xi = 10.0  # Skala mieszania
    
    return np.exp(-abs(Q1 - Q2) / xi)

# Kwarki: down, up, charm, strange, bottom, top
quarks = [
    ("u", 2, 3),   # up: T(2,3), Q=5
    ("d", 3, 4),   # down: T(3,4), Q=7
    ("c", 8, 13),  # charm: T(8,13), Q=21
    ("s", 8, 13),  # strange: podobne do charm
    ("t", 21, 34), # top: duży węzeł
    ("b", 13, 21), # bottom: T(13,21)
]

log("Przybliżona struktura CKM (qualitative):")
log("")

# Macierz CKM: przejścia u-d, c-s, t-b
V_ud = knot_overlap(2, 3, 3, 4)
V_us = knot_overlap(2, 3, 8, 13)
V_cd = knot_overlap(8, 13, 3, 4)
V_cs = knot_overlap(8, 13, 8, 13)

log(f"|V_ud| ≈ {V_ud:.4f} (exp: 0.974)")
log(f"|V_us| ≈ {V_us:.4f} (exp: 0.225)")
log(f"|V_cd| ≈ {V_cd:.4f} (exp: 0.221)")
log(f"|V_cs| ≈ {V_cs:.4f} (exp: 0.973)")

# Hierarchia
v_ud_exp = 0.974
v_us_exp = 0.225

hierarchy_ok = (V_ud > V_us) and (V_cs > V_cd)
log(f"\nHierarchia |V_ud| > |V_us|: {'✅' if hierarchy_ok else '❌'}")

status_1619 = "PARTIAL" if hierarchy_ok else "FAILED"
log(f"QW-1619 STATUS: {status_1619} (qualitative mapping only)")
log("UWAGA: Wysokie ryzyko numerologii!")

# =============================================================================
# PODSUMOWANIE
# =============================================================================
log("")
log("=" * 60)
log("PODSUMOWANIE OPTIONAL STUDIES")
log("=" * 60)
log(f"QW-1617 (Running α):    {status_1617}")
log(f"QW-1618 (Preon Stab.):  {status_1618}")
log(f"QW-1619 (CKM Knots):    {status_1619}")

# =============================================================================
# GENEROWANIE RAPORTU
# =============================================================================
with open(REPORT_FILE, 'w', encoding='utf-8') as f:
    f.write("# QW-1617/1618/1619: Optional Studies (Combined)\n\n")
    f.write(f"**Data:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    
    f.write("## Summary\n\n")
    f.write("| Study | Topic | Status |\n")
    f.write("|-------|-------|--------|\n")
    f.write(f"| QW-1617 | Running α(Q²) | {status_1617} |\n")
    f.write(f"| QW-1618 | Preon Stability | {status_1618} |\n")
    f.write(f"| QW-1619 | CKM from Knots | {status_1619} |\n\n")
    
    f.write("## Notes\n")
    f.write("- These studies are **model-dependent** and **speculative**\n")
    f.write("- They serve as **consistency checks**, not proofs\n")
    f.write("- QW-1619 carries **high numerology risk**\n\n")
    
    f.write("## Raw Log\n```\n" + '\n'.join(md) + "\n```\n")

print(f"\n✅ Report saved to: {REPORT_FILE}")
