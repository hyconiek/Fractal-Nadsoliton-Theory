#!/usr/bin/env python3
# QW-604b: CLEAN WAVE DISPERSION (No Attractor, No Renorm)
# Purpose: Verify if b=2.3 is physics or numerical artifact
# Date: 2025-12-05

import numpy as np
from scipy.ndimage import convolve

print("="*80)
print("QW-604b: CLEAN WAVE DISPERSION (Validation)")
print("="*80)
print("Test: Same as QW-604 but WITHOUT gamma and renormalization")
print("="*80)

GRID_SIZE = 64
DT = 0.05
STEPS = 300
GAMMA = 0.0  # NO ATTRACTOR!
BETA_TORS = 0.01

x = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
y = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
z = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

sigma_0 = 4.0
k_0 = 0.3

amplitude = np.exp(-(X**2 + Y**2 + Z**2) / (2 * sigma_0**2))
phase = k_0 * X
psi = amplitude * np.exp(1j * phase)
psi = psi / np.sqrt(np.sum(np.abs(psi)**2))

print(f"Gamma: {GAMMA} (NO ATTRACTOR)")
print("NO RENORMALIZATION")

laplacian_kernel = np.zeros((3,3,3))
laplacian_kernel[1,1,1] = -6
for i in [(0,1,1), (2,1,1), (1,0,1), (1,2,1), (1,1,0), (1,1,2)]:
    laplacian_kernel[i] = 1

def laplacian(field):
    return convolve(field, laplacian_kernel, mode='wrap')

def measure_width(psi):
    prob = np.abs(psi)**2
    prob = prob / (np.sum(prob) + 1e-10)
    x_cm = np.sum(X * prob)
    sigma_x = np.sqrt(np.sum(((X - x_cm)**2) * prob))
    sigma_y = np.sqrt(np.sum((Y**2) * prob))
    sigma_z = np.sqrt(np.sum((Z**2) * prob))
    return np.sqrt(sigma_x**2 + sigma_y**2 + sigma_z**2), x_cm

history_width = []
history_time = []

for t in range(STEPS):
    rho = np.abs(psi)**2
    # ONLY kinetic + back-reaction (NO gamma!)
    kin = 1j * laplacian(psi)
    back = -1j * BETA_TORS * rho * psi
    
    dpsi_dt = kin + back
    psi += DT * dpsi_dt
    
    # NO RENORMALIZATION!
    
    if t % 10 == 0:
        width, center = measure_width(psi)
        history_width.append(width)
        history_time.append(t * DT)
        
        if t % 50 == 0:
            print(f"t={t:4d}: σ={width:.2f}")

widths = np.array(history_width)
times = np.array(history_time)

n_skip = len(times) // 5
widths_fit = widths[n_skip:]
times_fit = times[n_skip:]

sigma_initial = widths_fit[0]
delta_sigma = widths_fit - sigma_initial

valid = delta_sigma > 0.1
if len(delta_sigma[valid]) > 5:
    log_t = np.log(times_fit[valid])
    log_ds = np.log(delta_sigma[valid])
    coeffs = np.polyfit(log_t, log_ds, 1)
    b = coeffs[0]
    
    print(f"\nCLEAN Dispersion: Δσ ∝ t^{b:.3f}")
    print(f"Compare to QW-604: b=2.328")
    print()
    
    if 0.9 < b < 1.1:
        print("✅ QUANTUM (b≈1) - QW-604 był artefaktem!")
    elif abs(b - 2.3) < 0.5:
        print("✅ REPRODUCED (b≈2.3) - to FIZYKA, nie artefakt!")
    else:
        print(f"🟡 DIFFERENT (b={b:.2f}) - częściowy artefakt")
else:
    b = 0

with open("raport_qw604b_clean.md", "w") as f:
    f.write("# QW-604b: Clean Dispersion (Validation)\n")
    f.write("**Gamma:** 0 (no attractor)\n")
    f.write("**Renormalization:** NO\n\n")
    f.write(f"**Result:** b = {b:.3f}\n")
    f.write(f"**QW-604:** b = 2.328\n\n")
    if abs(b - 2.3) < 0.5:
        f.write("### Wniosek: FIZYKA\n")
        f.write("Super-ballistic dispersion jest REALNA, nie artefakt numeryczny!\n")
    elif 0.9 < b < 1.1:
        f.write("### Wniosek: ARTEFAKT\n")
        f.write("Gamma term był odpowiedzialny za b=2.3. Prawdziwa wartość to b≈1.\n")

print("="*80)
