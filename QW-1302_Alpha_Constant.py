import numpy as np
from scipy.integrate import simpson
from datetime import datetime

# ==============================================================================
# QW-1302: OBLICZENIE STAŁEJ STRUKTURY SUBTELNEJ (ALPHA)
# ==============================================================================
#
# CEL: Spróbować wyprowadzić wartość alpha ~ 1/137 z czystej geometrii Skyrmiona.
#
# DANE z QW-1301:
# Stała oddziaływania C = E_int * D = 85.6 (jednostek geometrycznych).
#
# HIPOTEZA:
# Może alpha jest powiązana z odwrotnością stałej C lub stosunkiem C do E_rest?
#
# METODA:
# 1. Obliczyć E_rest (Energię Własną) dla tego samego profilu f(r).
#    E_rest = integral (grad f)^2 d^3x = 4pi * integral (f')^2 r^2 dr.
# 2. Sprawdzić różne kombinacje bezwymiarowe, np. alpha = 1/C, alpha = E_rest/C...
#
# ==============================================================================

REPORT_FILE = "RAPORT_QW1302_ALPHA_CALC.md"
md = []

def log(msg):
    print(msg)
    md.append(msg)

log("=" * 80)
log("QW-1302: POSZUKIWANIE GEOMETRYCZNEJ ALFY")
log("=" * 80)

def profile_deriv_sq(r, R_size=1.0):
    # f(r) = 2 * arctan(R/r)
    # f'(r) = 2 * (1 / (1 + (R/r)^2)) * (-R/r^2)
    # f'(r) = -2R / (r^2 + R^2)
    # f'(r)^2 = 4R^2 / (r^2 + R^2)^2
    return (2.0 * R_size / (r**2 + R_size**2))**2

def compute_self_energy(R_size=1.0):
    # E = 4pi * integral_0^inf (f')^2 * r^2 dr
    r = np.linspace(0, 1000, 100000)
    integrand = profile_deriv_sq(r, R_size) * r**2
    return 4 * np.pi * simpson(integrand, r)

log("\n[1] OBLICZENIA ENERGII")
E_rest = compute_self_energy()
C_coulomb = 85.655 # Średnia z QW-1301 (dla dużych D)

log(f"Energia Spoczynkowa E_rest: {E_rest:.4f}")
log(f"Stała Coulomba C (E*D):     {C_coulomb:.4f}")

log("\n[2] TESTOWANIE HIPOTEZ NA ALPHĘ")
log("-" * 60)

# Hipoteza 1: alpha = 1 / C
alpha_1 = 1.0 / C_coulomb
log(f"H1: alpha = 1/C             = {alpha_1:.6f} (Exp: 0.007297)")
log(f"    Różnica: {abs(alpha_1 - 1/137.036)/(1/137.036)*100:.2f}%")

# Hipoteza 2: alpha = 4pi / C (Normalizacja sferyczna pola)
# C odpowiada całce 4pi. Jeśli e^2 jest znormalizowane...
alpha_2 = 4 * np.pi / C_coulomb
log(f"H2: alpha = 4pi/C           = {alpha_2:.6f}")

# Hipoteza 3: Stosunek E_rest / C ?
# Masa / Ładunek?
alpha_3 = E_rest / C_coulomb
log(f"H3: alpha = E_rest/C        = {alpha_3:.6f}")

# Hipoteza 4: alpha = 1 / (E_rest) ?
alpha_4 = 1.0 / E_rest
log(f"H4: alpha = 1/E_rest        = {alpha_4:.6f}")

# Hipoteza 5 (Zasada holograficzna?): alpha = 1 / (C + E_rest)?
alpha_5 = 1.0 / (C_coulomb + E_rest)
log(f"H5: alpha = 1/(C+E)         = {alpha_5:.6f}")

log("\n[3] ANALIZA WYNIKÓW")
log("-" * 80)

# Szukamy wartości bliskiej 0.00729
candidates = [alpha_1, alpha_2, alpha_3, alpha_4, alpha_5]
best = min(candidates, key=lambda x: abs(x - 1/137.036))

log(f"Najlepszy kandydat: {best:.6f}")
if abs(best - 1/137.036) < 0.001:
    log("SUKCES! Znaleziono korelację geometryczną.")
else:
    log("PORAŻKA: Żadna prosta kombinacja nie daje 1/137.")
    log("Wniosek: Alpha wymaga ekranowania kwantowego (polaryzacji próżni),")
    log("sama klasyczna geometria daje zbyt silne sprzężenie.")
    
# Zapis raportu
with open(REPORT_FILE, 'w') as f:
    f.write("# QW-1302: Obliczenie Stałej Alpha\n\n")
    f.write(f"**Data:** {datetime.now()}\n\n")
    f.write("```\n")
    f.write("\n".join(md))
    f.write("\n```\n")

print(f"Raport zapisano w {REPORT_FILE}")
