#!/usr/bin/env python3
# QW-604: WAVE PACKET DISPERSION TEST
# Purpose: Test if FIN shows quantum (linear) or classical (diffusive) dispersion
# Prediction: Δx ∝ t^b where b=1 (quantum) or b=0.5 (classical diffusion)
# Date: 2025-12-05

import numpy as np
from scipy.ndimage import convolve

print("="*80)
print("QW-604: WAVE PACKET DISPERSION")
print("="*80)
print("Test: Quantum (Δx ∝ t) vs Classical (Δx ∝ √t)?")
print("="*80)

# ============================================================================
# PARAMETERS
# ============================================================================
GRID_SIZE = 64  # Larger for wave propagation
DT = 0.05
STEPS = 300
GAMMA = 0.1  # Lower to preserve wave packet
BETA_TORS = 0.01

print(f"Grid: {GRID_SIZE}^3")
print(f"Steps: {STEPS}")
print("-" * 40)

# ============================================================================
# INITIAL WAVE PACKET (Gaussian)
# ============================================================================

x = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
y = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
z = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# Gaussian wave packet centered at origin
# With initial momentum in +x direction
sigma_0 = 4.0  # Initial width
k_0 = 0.3       # Initial momentum (wave number)

# Amplitude: Gaussian
amplitude = np.exp(-(X**2 + Y**2 + Z**2) / (2 * sigma_0**2))

# Phase: plane wave in +x direction
phase = k_0 * X

psi = amplitude * np.exp(1j * phase)

# Normalize
psi = psi / np.sqrt(np.sum(np.abs(psi)**2))

print(f"Initial packet: σ₀={sigma_0:.1f}, k₀={k_0:.1f}")
print("Phase: plane wave in +x direction")

# ============================================================================
# EVOLUTION
# ============================================================================

laplacian_kernel = np.zeros((3,3,3))
laplacian_kernel[1,1,1] = -6
for i in [(0,1,1), (2,1,1), (1,0,1), (1,2,1), (1,1,0), (1,1,2)]:
    laplacian_kernel[i] = 1

def laplacian(field):
    return convolve(field, laplacian_kernel, mode='wrap')

def measure_width(psi):
    """Measure RMS width of wave packet"""
    prob = np.abs(psi)**2
    prob = prob / np.sum(prob)
    
    # Center of mass
    x_cm = np.sum(X * prob)
    y_cm = np.sum(Y * prob)
    z_cm = np.sum(Z * prob)
    
    # RMS width
    sigma_x = np.sqrt(np.sum(((X - x_cm)**2) * prob))
    sigma_y = np.sqrt(np.sum(((Y - y_cm)**2) * prob))
    sigma_z = np.sqrt(np.sum(((Z - z_cm)**2) * prob))
    
    # Total width (3D)
    sigma_total = np.sqrt(sigma_x**2 + sigma_y**2 + sigma_z**2)
    
    return sigma_total, x_cm

print("\nEvolving wave packet...")

history_width = []
history_center = []
history_time = []

for t in range(STEPS):
    rho = np.abs(psi)**2
    attractor = GAMMA * (1.0 - rho) * psi
    kin = 1j * laplacian(psi)
    back = -1j * BETA_TORS * rho * psi
    
    dpsi_dt = attractor + kin + back
    psi += DT * dpsi_dt
    
    # Renormalize to prevent decay
    psi = psi / np.sqrt(np.sum(np.abs(psi)**2) + 1e-10)
    
    if t % 10 == 0:
        width, center = measure_width(psi)
        history_width.append(width)
        history_center.append(center)
        history_time.append(t * DT)
        
        if t % 50 == 0:
            print(f"t={t:4d}: σ={width:.2f}, x_center={center:.2f}")

print("-" * 40)

# ============================================================================
# ANALYSIS: FIT DISPERSION LAW
# ============================================================================

widths = np.array(history_width)
times = np.array(history_time)

# Remove initial transient (first 20%)
n_skip = len(times) // 5
widths_fit = widths[n_skip:]
times_fit = times[n_skip:]

# Fit: σ(t) = σ₀ + A × t^b
# Take log: log(σ - σ₀) = log(A) + b × log(t)

# We need to determine σ₀ (initial width after relaxation)
sigma_initial = widths_fit[0]

# Dispersion: Δσ(t) = σ(t) - σ₀
delta_sigma = widths_fit - sigma_initial

# Avoid zeros/negatives
valid = delta_sigma > 0.1
delta_sigma_valid = delta_sigma[valid]
times_valid = times_fit[valid]

if len(times_valid) > 5:
    # Log-log fit
    log_t = np.log(times_valid)
    log_ds = np.log(delta_sigma_valid)
    
    # Linear fit in log space
    coeffs = np.polyfit(log_t, log_ds, 1)
    b = coeffs[0]  # Exponent!
    
    print(f"\nDispersion Law: Δσ ∝ t^{b:.3f}")
    print()
    
    # Interpret
    if 0.9 < b < 1.1:
        print("✅ QUANTUM DISPERSION!")
        print(f"   Exponent b={b:.3f} ≈ 1.0 (linear spreading)")
        print("   FIN shows Schrödinger-like dynamics")
        dispersion_type = "quantum"
    elif 0.4 < b < 0.6:
        print("❌ CLASSICAL DIFFUSION")
        print(f"   Exponent b={b:.3f} ≈ 0.5 (diffusive spreading)")
        print("   FIN shows classical random walk")
        dispersion_type = "classical"
    else:
        print("🟡 ANOMALOUS DISPERSION")
        print(f"   Exponent b={b:.3f} (neither quantum nor classical)")
        dispersion_type = "anomalous"
else:
    b = 0
    dispersion_type = "unknown"

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw604_dispersion.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-604: Wave Packet Dispersion\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Test:** Quantum vs Classical Dynamics\n\n")
    
    f.write("## 1. Metodologia\n")
    f.write(f"- Początkowy pakiet Gaussa: σ₀={sigma_0}\n")
    f.write(f"- Pęd początkowy: k₀={k_0}\n")
    f.write(f"- Ewolucja: {STEPS} kroków\n\n")
    
    f.write("## 2. Ewolucja Szerokości\n")
    f.write("| Czas | σ (width) | x_center |\n")
    f.write("|------|-----------|----------|\n")
    for i in range(min(15, len(history_time))):
        f.write(f"| {history_time[i]:.1f} | {history_width[i]:.2f} | {history_center[i]:.2f} |\n")
    f.write("\n")
    
    f.write(f"## 3. Prawo Dyspersji\n")
    f.write(f"**Fit:** Δσ ∝ t^{b:.3f}\n\n")
    
    if dispersion_type == "quantum":
        f.write("### ✅ KWANTOWA DYSPERSJA\n")
        f.write(f"Wykładnik b={b:.3f} ≈ 1.0 potwierdza że pakiety falowe\n")
        f.write("rozchodzą się **liniowo** z czasem, jak w mechanice kwantowej!\n\n")
        f.write("**Implikacja:** FIN ma kwantową dynamikę (mimo braku splątania Bell).\n")
    elif dispersion_type == "classical":
        f.write("### ❌ KLASYCZNA DYFUZJA\n")
        f.write(f"Wykładnik b={b:.3f} ≈ 0.5 pokazuje że pakiety\n")
        f.write("rozchodzą się dyfuzyjnie, jak w klasycznym błądzeniu losowym.\n")
    else:
        f.write("### 🟡 ANOMALNA DYSPERSJA\n")
        f.write(f"Wykładnik b={b:.3f} nie pasuje do standardowych modeli.\n")

print("Report saved.")
print("="*80)
