# SUPER EXTREME LOG PHASE 2
**Weak Decay.**

## QW-1400
### S:QW-1400_Link_Relaxation_Dynamics.py
```python
REPORT_FILE = "RAPORT_QW1400_MUON_DECAY.md"
md = []
def log(msg):
    print(msg)
    md.append(msg)
log("=" * 80)
log("QW-1400: SYMULACJA ROZPADU MIONU (RELAXATION DYNAMICS)")
log("=" * 80)
M_PREON_EFF = 2553.0 
K_SPRING = 3100.0 
TIME_STEP = 0.0005 
STEPS = 50000000 
GAMMA = 0.5 
z_positions = np.array([0.0, 0.0, 0.26]) 
velocities = np.array([0.0, 0.0, 0.0])
masses = np.array([M_PREON_EFF, M_PREON_EFF, M_PREON_EFF])
log("[1] WARUNKI POCZĄTKOWE")
log(f"Konfiguracja: MION (Deformacja z3 = {z_positions[2]:.4f})")
log(f"Energia wzbudzenia (Masa Mionu): ~105 MeV")
log(f"Stała sprężystości wiązania k: {K_SPRING}")
log(f"Tłumienie (Neutrina) gamma: {GAMMA}")
trajectory = []
energies = []
log("\n[2] START SYMULACJI CZASOWEJ...")
for t in range(STEPS):
    forces = np.zeros(3)
    forces[2] = -K_SPRING * z_positions[2]
    forces[2] -= GAMMA * velocities[2]
    accelerations = forces / masses
    z_positions += velocities * TIME_STEP + 0.5 * accelerations * TIME_STEP**2
    velocities += accelerations * TIME_STEP 
    if t % 100 == 0:
        T = 0.5 * masses[2] * velocities[2]**2
        V = 0.5 * K_SPRING * z_positions[2]**2
        E_tot = T + V
        trajectory.append(z_positions[2])
        energies.append(E_tot)
log("\n[3] ANALIZA KOŃCOWA")
final_z = z_positions[2]
final_mass_excess = 0.5 * K_SPRING * final_z**2
log(f"Pozycja końcowa z: {final_z:.6f}")
log(f"Masa resztkowa wzbudzenia: {final_mass_excess:.6f} MeV")
log("-" * 60)
if abs(final_z) < 0.01:
    log("WYNIK: SUKCES. Układ zrelaksował się do stanu ELEKTRONU.")
    log("Mechanizm: Energia wzbudzenia (105 MeV) została wypromieniowana")
    log("jako ciepło/fale (neutrina).")
    log(f"Zauważono oscylacje tłumione (Oscylacje neutrinowe?).")
else:
    log("WYNIK: Układ utknął w stanie metastabilnym.")
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
    f.write("
    f.write(f"**Data:** {datetime.now()}\n\n")
    f.write("```\n")
    f.write("\n".join(md))
    f.write("\n```\n")
print(f"Raport zapisano w {REPORT_FILE}")
```
### R:RAPORT_QW1400_MUON_DECAY.md
```markdown
# QW-1400: Symulacja Rozpadu Mionu
```
--------------------
