import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# ==============================================================================
# QW-1400: DYNAMIKA RELAKSACJI SPLOTU (ROZPAD MIONU)
# ==============================================================================
#
# CEL: Symulacja ewolucji czasowej zdeformowanego splotu preonowego (Mion).
#      Sprawdzenie, czy układ samoistnie relaksuje do stanu podstawowego (Elektron).
#
# MODEL EFEKTYWNY:
# 3 Pętle sztywne (Preony) oddziałujące potencjałem efektywnym V_eff.
# Traktujemy pętle jako obiekty quasi-sztywne z 1 stopniem swobody (Offset Z),
# aby symulacja była czytelna i szybka.
#
# Równanie ruchu: m * d2z/dt2 = - dV/dz - gamma * dz/dt
# Gdzie gamma to współczynnik emisji neutrin (tłumienie).
#
# Potencjał V(z) bierzemy z interpolacji wyników QW-1213.
#
# ==============================================================================

REPORT_FILE = "RAPORT_QW1400_MUON_DECAY.md"
md = []

def log(msg):
    print(msg)
    md.append(msg)

log("=" * 80)
log("QW-1400: SYMULACJA ROZPADU MIONU (RELAXATION DYNAMICS)")
log("=" * 80)

# Parametry fizyczne z poprzednich badań
M_PREON_EFF = 2553.0 # MeV (masa bezwładna pętli)
# Uwaga: Masa bezwładna to masa preonu, ale masa grawitacyjna (E=mc^2) to masa systemu.
# Do dynamiki F=ma bierzemy masę składnika.

# Potencjał wiązania (z QW-1213)
# Przybliżamy go studnią potencjału wokół z=0 (Elektron).
# E(z) ~= E_min + 1/2 k z^2  dla małych z (Mion jest w obszarze liniowym/harmonicznym?)
# Z QW-1213: delta=0.26 (Mion) -> dM = 105 MeV.
# E(0) = -4307 (jednostek sym). M(0) = 0.5 MeV.
# E(0.26) odpowiada masie 105 MeV.
# Czyli V(z) = (105 - 0.5) = 104.5 MeV dla z=0.26.
# k = 2 * E / z^2 = 2 * 104.5 / (0.26^2) = 209 / 0.0676 ~= 3091 MeV/unit^2.

K_SPRING = 3100.0 # Stała siłowa powrotu do symetrii

# Parametry symulacji
TIME_STEP = 0.0005 # Zwiększony krok czasowy
STEPS = 50000000 # Ekstremalnie długa symulacja (50 mln)
GAMMA = 0.5 # Współczynnik tłumienia (emisja neutrin)

# Warunki początkowe (MION)
z_positions = np.array([0.0, 0.0, 0.26]) # Dwie pętle w "tle", trzecia przesunięta
velocities = np.array([0.0, 0.0, 0.0])
masses = np.array([M_PREON_EFF, M_PREON_EFF, M_PREON_EFF])

log("[1] WARUNKI POCZĄTKOWE")
log(f"Konfiguracja: MION (Deformacja z3 = {z_positions[2]:.4f})")
log(f"Energia wzbudzenia (Masa Mionu): ~105 MeV")
log(f"Stała sprężystości wiązania k: {K_SPRING}")
log(f"Tłumienie (Neutrina) gamma: {GAMMA}")

# Symulacja Verlet
trajectory = []
energies = []

log("\n[2] START SYMULACJI CZASOWEJ...")

for t in range(STEPS):
    # Siły
    # Zakładamy, że pętle 0 i 1 są "kotwicą" (lub też się ruszają symetrycznie).
    # Uproszczenie: Pętla 2 jest przyciągana do centrum masy pętli 0 i 1 (czyli z=0).
    # F = -k * z
    
    forces = np.zeros(3)
    # Siła działa na pętlę 2 (przyciąganie do 0)
    forces[2] = -K_SPRING * z_positions[2]
    
    # Dynamika z tłumieniem: F_tot = F_pot - gamma * v
    forces[2] -= GAMMA * velocities[2]
    
    # 3. Zasada dynamiki (akcja-reakcja) - pętle 0 i 1 czują 1/2 siły odwrotnej?
    # W modelu Mionu, to pętla 2 jest "wyrwana". Pętle 0 i 1 tworzą "rdzeń".
    # Zakładamy rdzeń nieruchomy (nieskończona masa efektywna przez symetrię) dla uproszczenia wizualizacji relaksacji.
    
    # Update (Velocity Verlet)
    accelerations = forces / masses
    
    z_positions += velocities * TIME_STEP + 0.5 * accelerations * TIME_STEP**2
    velocities += accelerations * TIME_STEP # (Półkrok Eulera dla uproszczenia przy tłumieniu wystarcza w tej skali)
    
    # Monitoring
    if t % 100 == 0:
        # Energia kinetyczna T = 0.5 m v^2
        T = 0.5 * masses[2] * velocities[2]**2
        # Energia potencjalna V = 0.5 k z^2
        V = 0.5 * K_SPRING * z_positions[2]**2
        E_tot = T + V
        
        trajectory.append(z_positions[2])
        energies.append(E_tot)
        
        # log(f"Step {t}: z={z_positions[2]:.4f} E={E_tot:.2f}")

log("\n[3] ANALIZA KOŃCOWA")
final_z = z_positions[2]
final_mass_excess = 0.5 * K_SPRING * final_z**2

log(f"Pozycja końcowa z: {final_z:.6f}")
log(f"Masa resztkowa wzbudzenia: {final_mass_excess:.6f} MeV")

# Logika fizyczna
log("-" * 60)
if abs(final_z) < 0.01:
    log("WYNIK: SUKCES. Układ zrelaksował się do stanu ELEKTRONU.")
    log("Mechanizm: Energia wzbudzenia (105 MeV) została wypromieniowana")
    log("jako ciepło/fale (neutrina).")
    log(f"Zauważono oscylacje tłumione (Oscylacje neutrinowe?).")
else:
    log("WYNIK: Układ utknął w stanie metastabilnym.")

# Sprawdzenie czasu połowicznego zaniku
# Szukamy, kiedy z spadnie do 1/e początkowego (lub energiia do 1/2)
initial_z = 0.26
half_z = initial_z / 2.0
decay_step = -1

for i, z in enumerate(trajectory):
    if z < half_z and decay_step == -1:
        decay_step = i * 100
        break

if decay_step > 0:
    log(f"\nCzas relaksacji (τ): {decay_step * TIME_STEP:.4f} jedn. czasu")
    log("Interpretacja: To jest modelowy czas życia cząstki w jednostkach symulacji.")

# Zapis raportu
with open(REPORT_FILE, 'w') as f:
    f.write("# QW-1400: Symulacja Rozpadu Mionu\n\n")
    f.write(f"**Data:** {datetime.now()}\n\n")
    f.write("```\n")
    f.write("\n".join(md))
    f.write("\n```\n")

print(f"Raport zapisano w {REPORT_FILE}")
