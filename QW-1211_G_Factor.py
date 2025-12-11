import numpy as np
from datetime import datetime

# ==============================================================================
# QW-1211: MOMENT MAGNETYCZNY SPLOTU (G-FACTOR)
# ==============================================================================
#
# CEL: Obliczyć gyromagnetyczny czynnik 'g' dla splotu elektronowego T(21,3).
#
# FIZYKA:
# Klasyczna pętla z prądem ma g = 1.
# Elektron Diraca ma g = 2.
# Anomalie QED dają g = 2.002...
#
# Hipoteza: Złożona geometria T(21,3) (splot 3 pętli) może zmieniać
# stosunek momentu magnetycznego (mu) do momentu pędu (L).
#
# METODA:
# 1. Zdefiniuj krzywą T(21,3).
# 2. Załóż przepływ "masy i ładunku" wzdłuż krzywej (v = stała wzdłuż ds).
# 3. Oblicz L = sum(r x v * dm)  (gdzie dm ~ dl)
# 4. Oblicz mu = 1/2 * sum(r x v * dq) (gdzie dq ~ dl)
#
# Jeśli rozkład masy i ładunku jest tożsamy (dm ~ dq), to g=1 klasycznie.
# ALE: W teorii FIN masa wynika z krzywizny (E ~ k^2), a ładunek z topologii.
# Gęstość liniowa masy lambda_m może nie być stała!
# Gęstość masy: lambda_m ~ (krzywizna)^2 ? (Energia zginania)
#
# Sprawdzimy g dla modelu: dm ~ kappa^2 * dl, dq ~ dl.
# Jeśli masa jest skoncentrowana w obszarach dużego zagięcia, L zmaleje?
#
# ==============================================================================

REPORT_FILE = "RAPORT_QW1211_G_FACTOR.md"
md = []

def log(msg):
    print(msg)
    md.append(msg)

log("=" * 80)
log("QW-1211: OBLICZENIE CZYNNIKA G DLA SPLOTU T(21,3)")
log("=" * 80)

def knot_curve(p, q, N=2000):
    t = np.linspace(0, 2*np.pi, N, endpoint=False)
    R, r = 2.0, 1.0
    x = (R + r * np.cos(q * t)) * np.cos(p * t)
    y = (R + r * np.cos(q * t)) * np.sin(p * t)
    z = r * np.sin(q * t)
    return np.column_stack((x, y, z)), t

def curvature(r_vec):
    """Oblicza krzywiznę w każdym punkcie krzywej numerycznej."""
    # Pierwsza pochodna (styczna)
    dr = np.gradient(r_vec, axis=0)
    ds = np.linalg.norm(dr, axis=1)
    T = dr / ds[:, np.newaxis]
    
    # Druga pochodna (normalna)
    dT = np.gradient(T, axis=0)
    # Kappa = |dT/ds|
    kappa = np.linalg.norm(dT, axis=1) / ds
    return kappa, dr, ds

# Parametry
p, q = 21, 3
r_vec, t = knot_curve(p, q, N=5000)
kappa, dr, ds = curvature(r_vec)

# Model 1: Klasyczny (gęstość masy = stała, gęstość ładunku = stała)
# dm = dq = dl
log("\n[1] MODEL KLASYCZNY (dm = dq = stała)")
L_vec_1 = np.sum(np.cross(r_vec, dr), axis=0)
mu_vec_1 = 0.5 * np.sum(np.cross(r_vec, dr), axis=0)

# Projekcja na oś Z (dominującą dla torusa p >> q)
Lz_1 = L_vec_1[2]
muz_1 = mu_vec_1[2]
g_1 = (muz_1 / Lz_1) * 2 # Definicja g: mu = g * (q/2m) * L. Tutaj q=m=1 => mu = g/2 L => g = 2*mu/L
log(f"g_factor: {g_1:.4f} (Oczekiwane 1.0)")


# Model 2: Topologiczny (Masa z energii krzywizny)
# dm ~ kappa^2 * dl (energia elastyczna węzła)
# dq ~ dl (ładunek jest topologiczny, rozłożony równomiernie wzdłuż włókna)
log("\n[2] MODEL TOPOLOGICZNY (Masa ~ krzywizna^2)")

# Normalizacja masy całkowitej do 1
dm = kappa**2 * ds
dm = dm / np.sum(dm) * np.sum(ds) # Skalowanie by total mass = total length

dq = ds # Ładunek jednorodny

# Moment pędu L = sum (r x v * dm) = sum (r x (dr/dt) * dm)
# Zakładamy v ~ dr (ruch wzdłuż krzywej)
L_vec_2 = np.cross(r_vec, dr)
L_vec_2 = np.sum(L_vec_2 * dm[:, np.newaxis] / ds[:, np.newaxis], axis=0) 
# Uwaga: dr zawiera już długość kroku. v * dm ~ (dr/dt) * dm.
# Uproszczenie: Operujemy na całkach geometrycznych.
# L ~ integral (r x dl) * rho_m
# mu ~ 1/2 integral (r x dl) * rho_q

mu_vec_2 = 0.5 * np.sum(np.cross(r_vec, dr) * dq[:, np.newaxis] / ds[:, np.newaxis], axis=0)

Lz_2 = L_vec_2[2]
muz_2 = mu_vec_2[2]

ratio = muz_2 / Lz_2
g_2 = ratio * 2 # Zakładając całkowitą masę = całkowity ładunek = 1

log(f"Lz (ważone masą):   {Lz_2:.4f}")
log(f"muz (ważone ład.):  {muz_2:.4f}")
log(f"g_factor:           {g_2:.6f}")

log("\n[3] WNIOSKI")
log("-" * 80)
diff = g_2 - 2.0
log(f"Odchylenie od g=2: {diff:.6f}")

if abs(g_2 - 2.0) < 0.2:
    log("SUKCES: Mechanizm koncentracji masy w zakrzywieniach podbija g w kierunku 2!")
    log("Geometryczne wyjaśnienie anomalnego momentu magnetycznego.")
else:
    log("WYNIK NEUTRALNY: Sama krzywiznwa nie wystarcza by dać g=2.")

# Zapis raportu
with open(REPORT_FILE, 'w') as f:
    f.write("# QW-1211: Moment Magnetyczny Splotu\n\n")
    f.write(f"**Data:** {datetime.now()}\n\n")
    f.write("```\n")
    f.write("\n".join(md))
    f.write("\n```\n")

print(f"Raport zapisano w {REPORT_FILE}")
