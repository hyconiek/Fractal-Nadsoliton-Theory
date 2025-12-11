import numpy as np
from scipy.integrate import tplquad
from datetime import datetime

# ==============================================================================
# QW-1208: STABILNOŚĆ SPLOTÓW (LINK STABILITY)
# ==============================================================================
#
# CEL: Sprawdzić, czy splot 3 pętli T(7,1) (model elektronu) jest stanem
#      związanym.
#
# FIZYKA: W teorii FIN linie topologiczne zachowują się jak prądy (vortex filaments).
#      Dymamika wirów: równoległe wiry o tym samym znaku cyrkulacji PRZYCIĄGAJĄ SIĘ
#      (analogia do prawa Ampere'a: F/L = 2k I1 I2 / r).
#
# METODA: 
# 1. Zdefiniować geometrię 3 pętli T(7,1) splecionych w torusie.
# 2. Obliczyć energię oddziaływania (Interaction Energy) między pętlami.
#    E_int = - (mu/4pi) * I1*I2 * integral (dl1 . dl2) / r
# 3. Jeśli E_int < 0 (zysk energetyczny), splot jest stabilny.
#
# ==============================================================================

REPORT_FILE = "RAPORT_QW1208_LINK_STABILITY.md"
md = []

def log(msg):
    print(msg)
    md.append(msg)

log("=" * 80)
log("QW-1208: STABILNOŚĆ SPLOTU ELEKTRONOWEGO (T(21,3))")
log("=" * 80)

def loop_points(p, q, phase_offset, N=1000):
    """Generuje punkty pojedynczej pętli T(p,q) z przesunięciem fazowym."""
    t = np.linspace(0, 2*np.pi, N, endpoint=False)
    # Geometria torusa
    R, r = 2.0, 0.5 # r mniejsze, żeby pętle były blisko siebie
    
    # Przesunięcie fazowe przesuwa pętlę "po rurze torusa"
    # Kąt toroidalny to q*t + phase
    phi = q*t + phase_offset
    theta = p*t
    
    x = (R + r * np.cos(phi)) * np.cos(theta)
    y = (R + r * np.cos(phi)) * np.sin(theta)
    z = r * np.sin(phi)
    
    return np.column_stack((x, y, z))

def compute_interaction_energy(loop1, loop2):
    """
    Oblicza energię oddziaływania magnetycznego/wirowego dwóch pętli.
    E_int ~ - integral (dl1 . dl2) / |r1 - r2|
    (Znak minus, bo równoległe prądy się przyciągają -> energia maleje)
    """
    N = len(loop1)
    
    # Wektory styczne (dl)
    dl1 = np.diff(loop1, axis=0, append=loop1[:1])
    dl2 = np.diff(loop2, axis=0, append=loop2[:1])
    
    E_int = 0.0
    
    # Podwójna pętla (uproszczona O(N^2))
    # Można przyspieszyć biorąc co k-ty punkt
    step = 5 
    for i in range(0, N, step):
        r1 = loop1[i]
        d1 = dl1[i]
        for j in range(0, N, step):
            r2 = loop2[j]
            d2 = dl2[j]
            
            diff = r1 - r2
            dist = np.linalg.norm(diff)
            if dist < 0.01: dist = 0.01 # Cutoff
            
            dot_prod = np.dot(d1, d2)
            E_int -= dot_prod / dist
            
    return E_int

# KONFIGURACJA ELEKTRONU: 3 pętle T(7,1)
# Przesunięte o 2pi/3 względem siebie w kącie toroidalnym
log("\n[1] KONFIGURACJA")
log("Model: 3 pętle T(7,1) przesunięte o 120 stopni (2pi/3).")

loops = []
phases = [0, 2*np.pi/3, 4*np.pi/3]

for ph in phases:
    loops.append(loop_points(7, 1, ph, N=500)) # N mniejsze dla szybkości

# E_self (energia własna pętli) - pomijamy, interesuje nas wiązanie
# E_binding = E_int(1,2) + E_int(1,3) + E_int(2,3)

log("\n[2] OBLICZENIA ENERGII WIĄZANIA")
log("-" * 60)

E_bind = 0.0
pairs = [(0,1), (0,2), (1,2)]

for i, j in pairs:
    E_pair = compute_interaction_energy(loops[i], loops[j])
    log(f"E_int(Pętla {i+1}, Pętla {j+1}) = {E_pair:.4f}")
    E_bind += E_pair

log("-" * 60)
log(f"CAŁKOWITA ENERGIA WIĄZANIA: {E_bind:.4f}")

# ANALIZA STABILNOŚCI
log("\n[3] WNIOSKI")
log("-" * 80)

if E_bind < 0:
    log("WYNIK: E_bind < 0. Splot jest STABILNY energetycznie. ✅")
    log("Interpretacja: Równoległe prądy topologiczne (wiry) przyciągają się.")
    log("Siła wiążąca 'preony' to po prostu oddziaływanie Biot-Savarta w 4D.")
else:
    log("WYNIK: E_bind > 0. Splot jest NIESTABILNY (odpychanie). ❌")

# Zapis raportu
with open(REPORT_FILE, 'w') as f:
    f.write("# QW-1208: Stabilność Splotu Elektronowego\n\n")
    f.write(f"**Data:** {datetime.now()}\n\n")
    f.write("```\n")
    f.write("\n".join(md))
    f.write("\n```\n")

print(f"Raport zapisano w {REPORT_FILE}")
