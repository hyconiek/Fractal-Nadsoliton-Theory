#!/usr/bin/env python3
# QW-599: VACUUM TURBULENCE SPECTRUM
# Purpose: Test if vacuum fluctuations show turbulent (power-law) or thermal (white noise) spectrum
# Hypothesis H2: Próżnia = Turbulentny Eter
# Date: 2025-12-05

import numpy as np
from scipy.ndimage import convolve

print("="*80)
print("QW-599: VACUUM TURBULENCE SPECTRUM (Test H2)")
print("="*80)

# ============================================================================
# CONCEPT: Analyze power spectrum of vacuum field fluctuations
# Turbulent flow: P(k) ~ k^(-5/3) (Kolmogorov)
# Thermal noise: P(k) ~ constant (white noise)
# ============================================================================

GRID_SIZE = 64  # Larger for better FFT
DT = 0.01
EQUILIBRATION_STEPS = 300
MEASUREMENT_STEPS = 200
GAMMA = 0.3
BETA_TORS = 0.01

print(f"Grid: {GRID_SIZE}^3")
print(f"Equilibration: {EQUILIBRATION_STEPS} steps")
print(f"Measurement: {MEASUREMENT_STEPS} steps")
print("-" * 40)

# ============================================================================
# VACUUM FIELD INITIALIZATION
# ============================================================================

# Start with random vacuum fluctuations
np.random.seed(599)
psi = (np.random.randn(GRID_SIZE, GRID_SIZE, GRID_SIZE) + 
       1j * np.random.randn(GRID_SIZE, GRID_SIZE, GRID_SIZE))
psi = psi / np.sqrt(np.mean(np.abs(psi)**2))  # Normalize

print("Initialized vacuum field with random fluctuations")

# ============================================================================
# EVOLUTION
# ============================================================================

laplacian_kernel = np.zeros((3,3,3))
laplacian_kernel[1,1,1] = -6
for i in [(0,1,1), (2,1,1), (1,0,1), (1,2,1), (1,1,0), (1,1,2)]:
    laplacian_kernel[i] = 1

def laplacian(field):
    return convolve(field, laplacian_kernel, mode='wrap')

print("\nEvolving to equilibrium...")
for t in range(EQUILIBRATION_STEPS):
    rho = np.abs(psi)**2
    attractor = GAMMA * (1.0 - rho) * psi
    kin = 1j * laplacian(psi)
    back = -1j * BETA_TORS * rho * psi
    
    dpsi_dt = attractor + kin + back
    psi += DT * dpsi_dt

print("Equilibrium reached. Now measuring fluctuations...")

# Collect density snapshots
density_snapshots = []

for t in range(MEASUREMENT_STEPS):
    rho = np.abs(psi)**2
    attractor = GAMMA * (1.0 - rho) * psi
    kin = 1j * laplacian(psi)
    back = -1j * BETA_TORS * rho * psi
    
    dpsi_dt = attractor + kin + back
    psi += DT * dpsi_dt
    
    if t % 10 == 0:
        density_snapshots.append(rho.copy())

print(f"Collected {len(density_snapshots)} snapshots")

# ============================================================================
# POWER SPECTRUM ANALYSIS
# ============================================================================

print("\nComputing power spectrum...")

# Average power spectrum over all snapshots
power_spectrum_sum = np.zeros((GRID_SIZE, GRID_SIZE, GRID_SIZE))

for snapshot in density_snapshots:
    # Remove mean (fluctuations only)
    fluctuation = snapshot - np.mean(snapshot)
    
    # 3D FFT
    fft = np.fft.fftn(fluctuation)
    power = np.abs(fft)**2
    power_spectrum_sum += power

power_spectrum_avg = power_spectrum_sum / len(density_snapshots)

# Radial average (spherically symmetric)
kx = np.fft.fftfreq(GRID_SIZE) * GRID_SIZE
ky = np.fft.fftfreq(GRID_SIZE) * GRID_SIZE
kz = np.fft.fftfreq(GRID_SIZE) * GRID_SIZE
KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
K = np.sqrt(KX**2 + KY**2 + KZ**2)

# Bin by k
k_bins = np.arange(0.5, GRID_SIZE//2, 1.0)
P_k = []

for i in range(len(k_bins)-1):
    mask = (K >= k_bins[i]) & (K < k_bins[i+1])
    if np.sum(mask) > 0:
        P_k.append(np.mean(power_spectrum_avg[mask]))
    else:
        P_k.append(0)

k_values = (k_bins[:-1] + k_bins[1:]) / 2
P_k = np.array(P_k)

# Filter out zeros
valid = P_k > 0
k_values = k_values[valid]
P_k = P_k[valid]

print(f"Valid k-bins: {len(k_values)}")

# ============================================================================
# FIT POWER LAW
# ============================================================================

if len(k_values) > 5:
    # Fit P(k) = A * k^alpha in log-log space
    log_k = np.log(k_values[k_values > 1])  # Avoid k=0
    log_P = np.log(P_k[k_values > 1])
    
    if len(log_k) > 2:
        coeffs = np.polyfit(log_k, log_P, 1)
        alpha = coeffs[0]  # Power law exponent
        
        print(f"\nPower law fit: P(k) ~ k^{alpha:.2f}")
        print()
        
        # Kolmogorov turbulence: alpha = -5/3 = -1.67
        # Thermal white noise: alpha ~ 0
        
        if -2.0 < alpha < -1.2:
            print("✅ TURBULENT SPECTRUM!")
            print(f"   Exponent {alpha:.2f} consistent with Kolmogorov -5/3 = -1.67")
            turbulent = True
        elif abs(alpha) < 0.5:
            print("❌ WHITE NOISE (thermal)")
            print(f"   Exponent {alpha:.2f} ~ 0, inconsistent with turbulence")
            turbulent = False
        else:
            print("🟡 INTERMEDIATE SPECTRUM")
            print(f"   Exponent {alpha:.2f} doesn't match known patterns")
            turbulent = None
    else:
        turbulent = None
        alpha = 0
else:
    turbulent = None
    alpha = 0

print("-" * 40)

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw599_turbulence.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-599: Vacuum Turbulence Spectrum\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Hipoteza:** H2 - Próżnia = Turbulentny Eter\n\n")
    
    f.write("## 1. Metodologia\n")
    f.write("- Zainicjowano losowe fluktuacje próżni\n")
    f.write(f"- Ewolucja: {EQUILIBRATION_STEPS} kroków (równowaga) + {MEASUREMENT_STEPS} (pomiar)\n")
    f.write("- Obliczono spektrum mocy P(k) fluktuacji gęstości\n\n")
    
    f.write("## 2. Spektrum Mocy\n")
    f.write("| k (częstość) | P(k) (moc) |\n")
    f.write("|--------------|------------|\n")
    for i in range(min(10, len(k_values))):
        f.write(f"| {k_values[i]:.2f} | {P_k[i]:.2e} |\n")
    f.write("\n")
    
    f.write(f"## 3. Analiza Prawa Potęgowego\n")
    f.write(f"**Fit:** P(k) ~ k^{alpha:.2f}\n\n")
    
    if turbulent:
        f.write("### ✅ TURBULENCJA POTWIERDZONA\n")
        f.write(f"Wykładnik {alpha:.2f} jest blisko Kołmogorowskiego -5/3 = -1.67.\n\n")
        f.write("**Wniosek:** Próżnia w FIN wykazuje turbulentne właściwości!\n")
        f.write("Nie jest to zwykły szum termiczny, ale chaotyczny przepływ energii przez skale.\n")
    elif turbulent == False:
        f.write("### ❌ BRAK TURBULENCJI\n")
        f.write(f"Wykładnik {alpha:.2f} ~ 0 sugeruje szum biały (termiczny), nie turbulencję.\n")
    else:
        f.write("### 🟡 WYNIK NIEJEDNOZNACZNY\n")
        f.write("Spektrum nie pasuje do znanych modeli. Potrzebne dalsze badania.\n")

print("Report saved.")
print("="*80)
