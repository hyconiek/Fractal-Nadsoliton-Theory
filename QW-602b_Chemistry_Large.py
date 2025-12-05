#!/usr/bin/env python3
# QW-602b: CHEMISTRY (Larger Spacing Test)
# Quick test with 2x larger spacing to prevent collapse
# Date: 2025-12-05

import numpy as np
from scipy.ndimage import convolve

print("="*80)
print("QW-602b: CHEMISTRY (LARGER SPACING)")
print("="*80)

GRID_SIZE = 40
DT = 0.01
STEPS = 600
GAMMA = 0.3
BETA_TORS = 0.01
SEPARATION = 24.0  # 2x larger than QW-602

angle_offset = 2 * np.pi / 3

def hopfion_field(grid_size, center, R=3.0, winding=1):
    x = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    y = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    z = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    X, Y, Z = X - center[0], Y - center[1], Z - center[2]
    rho = np.sqrt(X**2 + Y**2)
    rho[rho == 0] = 1e-10
    eta = np.arctan2(Z, rho - R)
    xi = np.arctan2(Y, X)
    phase = winding * (xi + eta)
    dist_to_ring = np.sqrt((rho - R)**2 + Z**2)
    amplitude = np.tanh(dist_to_ring / 1.5)
    return amplitude * np.exp(1j * phase)

R_triangle = SEPARATION / np.sqrt(3)
center_A = (R_triangle * np.cos(0), R_triangle * np.sin(0), 0)
center_B = (R_triangle * np.cos(angle_offset), R_triangle * np.sin(angle_offset), 0)
center_C = (R_triangle * np.cos(2*angle_offset), R_triangle * np.sin(2*angle_offset), 0)

psi = 0.6 * (hopfion_field(GRID_SIZE, center_A, R=3.0, winding=+1)/3 +
             hopfion_field(GRID_SIZE, center_B, R=3.0, winding=+1)/3 +
             hopfion_field(GRID_SIZE, center_C, R=3.0, winding=-1)/3)

print(f"Spacing: {SEPARATION} (2x larger)")

laplacian_kernel = np.zeros((3,3,3))
laplacian_kernel[1,1,1] = -6
for i in [(0,1,1), (2,1,1), (1,0,1), (1,2,1), (1,1,0), (1,1,2)]:
    laplacian_kernel[i] = 1

def laplacian(field):
    return convolve(field, laplacian_kernel, mode='wrap')

history_energy = []
initial_E = np.sum(np.abs(np.gradient(psi)[0])**2 + (1.0 - np.abs(psi)**2)**2)

for t in range(STEPS):
    rho = np.abs(psi)**2
    psi += DT * (GAMMA * (1.0 - rho) * psi + 1j * laplacian(psi) - 1j * BETA_TORS * rho * psi)
    
    if t % 100 == 0:
        E = np.sum(np.abs(np.gradient(psi)[0])**2 + (1.0 - np.abs(psi)**2)**2)
        history_energy.append(E)
        print(f"t={t:4d}: E={E:.1f}")

final_E = history_energy[-1]
print(f"\nInitial E: {initial_E:.1f}")
print(f"Final E:   {final_E:.1f}")
print(f"ΔE:        {final_E - initial_E:+.1f} ({(final_E-initial_E)/initial_E*100:+.1f}%)")

if abs((final_E - initial_E) / initial_E) < 0.1:
    print("✅ STABILNY (energia zachowana ±10%)")
else:
    print("🟡 NIESTABILNY lub energetycznie aktywny")

with open("raport_qw602b_chemistry.md", "w") as f:
    f.write("# QW-602b: Larger Spacing Test\n")
    f.write(f"**Spacing:** {SEPARATION} (2x)\n\n")
    f.write(f"- Initial E: {initial_E:.1f}\n")
    f.write(f"- Final E: {final_E:.1f}\n")
    f.write(f"- Change: {(final_E-initial_E)/initial_E*100:+.1f}%\n")

print("="*80)
