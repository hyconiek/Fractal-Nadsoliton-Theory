import numpy as np
from datetime import datetime

# ==============================================================================
# QW-1210: SPÓJNOŚĆ SPINU (SPIN CONSISTENCY CHECK)
# ==============================================================================
#
# CEL: Sprawdzić, czy splot elektronowy T(21,3) zachowuje się jak fermion
#      zgodnie z ograniczeniami Finkelsteina-Rubinsteina (FR).
#
# TEORIA:
# Twierdzenie FR wiąże topologię z statystyką (spinem):
# exp(i * 2 * pi * J) = (-1)^W
# gdzie W to indeks 'winding number' mapowania konfiguracyjnego.
# Dla Skyrmionów W = B (ładunek topologiczny).
#
# PROBLEM:
# Elektron w FIN ma Q = 24.
# Jeśli W = Q, to (-1)^24 = +1 => Spin J całkowity (Bozon!).
# To przeczy eksperymentowi (elektron to fermion, J=1/2).
#
# HIPOTEZA RATUNKOWA:
# Indeks FR (W) dla węzłów/splotów NIE JEST równy sumie Q = p+q.
# W dla węzłów to "Linking Number" z "Frame Bundle" (wstęgą).
# W = Self-Linking Number (SL).
#
# SL(K) = Wr(K) + Tw(K) (Călugăreanu Theorem)
# Ale dla węzłów torusowych SL = p*q.
#
# SPRAWDZAMY:
# 1. Oblicz SL dla T(21,3).
# 2. Sprawdź parzystość.
#
# ==============================================================================

REPORT_FILE = "RAPORT_QW1210_SPIN_CHECK.md"
md = []

def log(msg):
    print(msg)
    md.append(msg)

log("=" * 80)
log("QW-1210: TEST SPÓJNOŚCI SPINU DLA T(21,3)")
log("=" * 80)

def compute_self_linking_number(p, q):
    """
    Oblicza Self-Linking Number dla węzła/splotu torusowego T(p,q).
    Wzór (Jones, 1987): SL = p * q
    """
    return p * q

log("\n[1] OBLICZENIA DLA ELEKTRONU T(21,3)")
log("-" * 60)

p, q = 21, 3
SL = compute_self_linking_number(p, q)
Q_fin = p + q

log(f"Węzeł/Splot: T({p},{q})")
log(f"Q (FIN charge): {Q_fin}")
log(f"SL (Self-Linking): {SL}")

# Test Finkelsteina-Rubinsteina
FR_phase = (-1)**SL
log(f"Faza FR (-1)^SL: {FR_phase}")

if FR_phase == -1:
    log("Statystyka: FERMION (J = 1/2, 3/2...) ✅")
    res = "FERMION"
else:
    log("Statystyka: BOZON (J = 0, 1...) ❌")
    res = "BOZON"

log("\n[2] OBLICZENIA DLA INNYCH CZĄSTEK")
log("-" * 60)

particles = [
    ("Electron", 21, 3), # Q=24
    ("Muon", 13, 1), # Q=14 (zakładając T(13,1))
    ("Tau/Charm", 8, 1), # Q=9
    ("Up", 8, 1), # Q=9
    ("Down", 5, 2) # Q=7
]

log(f"{'Particle':<10} {'T(p,q)':<10} {'Q':<5} {'SL':<5} {'Phase':<5} {'Type'}")
log("-" * 60)

for name, p, q in particles:
    sl = p * q
    phase = (-1)**sl
    type_str = "Fermion" if phase == -1 else "Bozon"
    marker = "✅" if phase == -1 else "❌"
    log(f"{name:<10} T({p},{q})   {p+q:<5} {sl:<5} {phase:<5} {type_str} {marker}")

log("\n[3] ANALIZA KRYTYCZNA")
log("-" * 80)
log("""
WYNIK: 
Dla T(21,3): SL = 63 (nieparzyste). (-1)^63 = -1.
TO JEST FERMION!

Mechanizm:
Ograniczenie Finkelsteina-Rubinsteina zależy od Self-Linking Number (SL = pq),
a nie od ładunku Q = p+q.

Elektron T(21,3):
- Q = 24 (parzyste) -> Masa
- SL = 63 (nieparzyste) -> Spin (Fermion)

Muon T(13,1):
- Q = 14
- SL = 13 (nieparzyste) -> Fermion

To cudowne zrządzenie matematyki.
Gdyby elektron był T(13,11) (Q=24), to SL = 143 (nieparzyste) -> Fermion.
Ale gdyby był T(12,12) (Q=24), to SL = 144 (parzyste) -> Bozon.

Warunek na bycie fermionem dla T(p,q):
p * q musi być NIEPARZYSTE.
To oznacza, że p i q muszą być OBA NIEPARZYSTE.

Sprawdźmy T(21,3): 21 (niep) * 3 (niep) = Nieparzyste. OK.
Sprawdźmy T(13,1): 13 * 1 = Nieparzyste. OK.
Sprawdźmy Down Quark T(5,2): 5 * 2 = 10 (Parzyste). BOZON?!

PROBLEM: Down Quark (fermion) wychodzi jako bozon w modelu T(5,2).
Rozwiązanie: Down quark to może T(6,1)? (Q=7, SL=6 -> Bozon). T(4,3)? (SL=12 -> Bozon).
Wszystkie pary sumujące się do nieparzystego Q (7, 9) mają iloczyn parzysty!
(Suma nieparzysta = Parzysta + Nieparzysta. Iloczyn P * N = Parzysty).

WNIOSEK FUNDAMENTALNY:
Modele T(p,q) dla fermionów o NIEPARZYSTYM Q (kwarki, tau) są BOZONOWE w tym obrazie.
Fermiony o PARZYSTYM Q (elektron, mion) mogą być FERMIONOWE.

To sugeruje, że nasza identyfikacja T(p,q) dla kwarków jest BŁĘDNA lub
brakuje składnika "skręcenia" (Twist) w formule SL.
""")

# Zapis raportu
with open(REPORT_FILE, 'w') as f:
    f.write("# QW-1210: Test Spójności Spinu\n\n")
    f.write(f"**Data:** {datetime.now()}\n\n")
    f.write("```\n")
    f.write("\n".join(md))
    f.write("\n```\n")

print(f"Raport zapisano w {REPORT_FILE}")
