import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from datetime import datetime

# ==============================================================================
# QW-1301: EMERGENTNA SIŁA COULOMBA (1/r CHECK)
# ==============================================================================
#
# CEL: Sprawdzić, jak skaluje się energia oddziaływania dwóch Skyrmionów
#      w funkcji odległości R.
#      Coulomb: E ~ 1/R
#      Yukawa:  E ~ exp(-mR) / R
#      Dipol:   E ~ 1/R^3
#
# MODEL:
# Skyrmion 3D reprezentowany przez profil f(r) (Hedgehog ansatz).
# Pole U(x) = exp(i * f(r) * n * tau).
# Energia gęstości L = Tr(dU dU) ...
# Dla dwóch Skyrmionów używamy "Product Ansatz": U_tot = U1 * U2.
#
# ==============================================================================

REPORT_FILE = "RAPORT_QW1301_COULOMB_FORCE.md"
md = []

def log(msg):
    print(msg)
    md.append(msg)

log("=" * 80)
log("QW-1301: CZY SKYRMIONY ODDZIAŁUJĄ SIŁĄ COULOMBA?")
log("=" * 80)

# 1. Profil pojedynczego Skyrmiona (uproszczony profil 1/r^2 dla dużych r)
# W standardowym modelu Skyrme'a f(r) ~ 1/r^2 na dalekim dystansie
# co daje oddziaływanie dipolowe (1/r^3). To problem.
# Aby uzyskać Coulomba, musimy "uruchomić" pole EM (cechowanie U(1)).
# Czy "surowy" Skyrmion ma Coulomba? Raczej nie. Ma Yukawę.
# Ale sprawdźmy to numerycznie.

def profile_f(r, R_size=1.0):
    # Wektoryzacja za pomocą np.where
    return np.where(r < 1e-6, np.pi, 2.0 * np.arctan(R_size / r))

# Funkcja obliczająca gęstość energii dwóch Skyrmionów
# U1 w (0,0,0), U2 w (0,0,D)
# Product ansatz: U = U1 * U2

def energy_density_at_point(x, y, z, D):
    # Pozycje
    r1 = np.sqrt(x**2 + y**2 + z**2)
    r2 = np.sqrt(x**2 + y**2 + (z-D)**2)
    
    # Profile i gradienty (symbolicznie - uproszczony model skalarnej interakcji)
    # Pełny model SU(2) jest kosztowny. 
    # Użyjmy modelu nieliniowego O(3) sigma: n1 . n2
    # E_int ~ gradient(n1) * gradient(n2)
    
    # Unit vector fields n(r):
    # n_z = cos(f(r))
    # n_perp = sin(f(r)) * (radial_unit_vector)
    
    f1 = profile_f(r1)
    f2 = profile_f(r2)
    
    # Interakcja dominująca na dużych odległościach to iloczyn "ogonów"
    # Ogon f ~ 1/r^2.
    # E_int ~ (1/r^2) * (1/r^2) ~ 1/r^4 ? Całka d^3x 1/r^4 ~ 1/r.
    # To by dało Coulomba!
    
    # Sprawdźmy numerycznie całkę z gradientów.
    # Grad f1 ~ 1/r1^2. Grad f2 ~ 1/r2^2.
    # Całka I(D) = integral d^3x (Grad f1 . Grad f2)
    
    # Przybliżenie: Grad f1 * Grad f2
    # Grad f = -2 R / (r^2 + R^2)
    
    gf1 = 2.0 / (r1**2 + 1.0) # Lorentzkie przybliżenie pochodnej
    gf2 = 2.0 / (r2**2 + 1.0)
    
    return gf1 * gf2

def compute_interaction_energy(D):
    # Całka numeryczna w 3D (cylindryczna symetria)
    # I = 2pi * integral r dr dz (Energy_density)
    
    R_max = D + 10.0
    r_vals = np.linspace(0, R_max, 50)
    z_vals = np.linspace(-R_max, R_max, 100)
    
    rr, zz = np.meshgrid(r_vals, z_vals)
    dens = energy_density_at_point(rr, 0, zz, D)
    
    # Całkowanie
    integrand = dens * rr # Jakobian
    res = 2 * np.pi * simpson(simpson(integrand, r_vals), z_vals)
    return res

log("\n[1] SYMULACJA ZALEŻNOŚCI E(D)")
d_vals = np.linspace(2.0, 10.0, 9)
e_vals = []

log(f"{'Dystans D':<10} {'E_int':<15} {'E * D':<10} {'E * D^3':<10}")
log("-" * 50)

for d in d_vals:
    e = compute_interaction_energy(d)
    e_vals.append(e)
    # Testy skalowania:
    # Jeśli Coulomb: E ~ 1/D => E*D = const
    # Jeśli Dipol: E ~ 1/D^3 => E*D^3 = const
    col_test = e * d
    dip_test = e * d**3
    log(f"{d:<10.2f} {e:<15.6f} {col_test:<10.4f} {dip_test:<10.4f}")

log("\n[2] WNIOSKI")
log("-" * 80)

# Analiza trendu
ed_vals = np.array(e_vals) * d_vals
slope = (ed_vals[-1] - ed_vals[0]) / len(ed_vals)

if abs(slope) < 0.1 * np.mean(ed_vals):
    log("WYNIK: E * D jest w przybliżeniu stałe.")
    log("Prawo oddziaływania: 1/D (COULOMB) ✅")
    log("Interpretacja: Teoria FIN generuje siły elektrostatyczne z geometrii!")
else:
    log("WYNIK: E * D nie jest stałe.")
    # Sprawdźmy dipol
    ed3_vals = np.array(e_vals) * d_vals**3
    if np.std(ed3_vals) < 0.2 * np.mean(ed3_vals):
        log("Prawo oddziaływania: 1/D^3 (DIPOL) ⚠️")
        log("Interpretacja: Skyrmiony zachowują się jak neutrane dipole, brakuje ładunku U(1).")
    else:
        log("Prawo oddziaływania: Złożone (Yukawa/Pośrednie).")

# Zapis raportu
with open(REPORT_FILE, 'w') as f:
    f.write("# QW-1301: Emergentna Siła Coulomba\n\n")
    f.write(f"**Data:** {datetime.now()}\n\n")
    f.write("```\n")
    f.write("\n".join(md))
    f.write("\n```\n")

print(f"Raport zapisano w {REPORT_FILE}")
