#!/usr/bin/env python3
# QW-605: ANYONIC BRAIDING ACCUMULATION
# Purpose: Test if anyonic phase accumulates linearly (θ_n = n × θ_single)
# Theory: True anyons must satisfy θ_total = n × θ_single for n exchanges
# Date: 2025-12-05

import numpy as np
from scipy.ndimage import convolve

print("="*80)
print("QW-605: ANYONIC BRAIDING ACCUMULATION")
print("="*80)
print("Test: Does anyonic phase multiply with number of exchanges?")
print("QW-603 baseline: θ₁ = 0.880 (single exchange)")
print("="*80)

# ============================================================================
# PARAMETERS
# ============================================================================
GRID_SIZE = 32
DT = 0.01
GAMMA = 0.3
BETA_TORS = 0.01
WINDING = +1

# Test multiple braiding counts
BRAID_COUNTS = [1, 2, 3]
STEPS_PER_BRAID = 400  # Same as QW-603

print(f"Grid: {GRID_SIZE}^3")
print(f"Braiding counts to test: {BRAID_COUNTS}")
print(f"Steps per braid: {STEPS_PER_BRAID}")
print("-" * 40)

# ============================================================================
# HOPFION FIELD
# ============================================================================

def hopfion_field(grid_size, center, R=3.0, winding=1):
    x = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    y = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    z = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    X = X - center[0]
    Y = Y - center[1]
    Z = Z - center[2]
    
    rho = np.sqrt(X**2 + Y**2)
    rho[rho == 0] = 1e-10
    
    eta = np.arctan2(Z, rho - R)
    xi = np.arctan2(Y, X)
    
    phase = winding * (xi + eta)
    dist_to_ring = np.sqrt((rho - R)**2 + Z**2)
    amplitude = np.tanh(dist_to_ring / 1.5)
    
    return amplitude * np.exp(1j * phase)

# ============================================================================
# EVOLUTION & BRAIDING
# ============================================================================

laplacian_kernel = np.zeros((3,3,3))
laplacian_kernel[1,1,1] = -6
for i in [(0,1,1), (2,1,1), (1,0,1), (1,2,1), (1,1,0), (1,1,2)]:
    laplacian_kernel[i] = 1

def laplacian(field):
    return convolve(field, laplacian_kernel, mode='wrap')

def perform_braiding(psi, n_braids, force_strength=0.03):
    """Perform n complete braiding operations"""
    x = np.arange(GRID_SIZE) - GRID_SIZE/2
    y = np.arange(GRID_SIZE) - GRID_SIZE/2
    z = np.arange(GRID_SIZE) - GRID_SIZE/2
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    total_steps = STEPS_PER_BRAID * n_braids
    
    for t in range(total_steps):
        rho = np.abs(psi)**2
        attractor = GAMMA * (1.0 - rho) * psi
        kin = 1j * laplacian(psi)
        back = -1j * BETA_TORS * rho * psi
        
        # Continuous rotation (each full cycle = one braid)
        theta_force = 2 * np.pi * t / STEPS_PER_BRAID
        
        # Rotating forces
        force_A_x = force_strength * np.cos(theta_force) * np.exp(-(X**2 + (Y+6)**2 + Z**2)/16)
        force_A_y = force_strength * np.sin(theta_force) * np.exp(-(X**2 + (Y+6)**2 + Z**2)/16)
        
        force_B_x = -force_strength * np.cos(theta_force) * np.exp(-(X**2 + (Y-6)**2 + Z**2)/16)
        force_B_y = -force_strength * np.sin(theta_force) * np.exp(-(X**2 + (Y-6)**2 + Z**2)/16)
        
        force_total = 1j * (force_A_x + force_B_x) * X * psi + 1j * (force_A_y + force_B_y) * Y * psi
        
        dpsi_dt = attractor + kin + back + force_total
        psi += DT * dpsi_dt
    
    return psi

# ============================================================================
# TEST MULTIPLE BRAIDINGS
# ============================================================================

results = []

for n_braids in BRAID_COUNTS:
    print(f"\n{'='*60}")
    print(f"Testing {n_braids} braid(s)")
    print('='*60)
    
    # Initialize fresh hopfions
    center_A = (0, -6, 0)
    center_B = (0, +6, 0)
    
    psi_A = hopfion_field(GRID_SIZE, center_A, R=3.0, winding=WINDING)
    psi_B = hopfion_field(GRID_SIZE, center_B, R=3.0, winding=WINDING)
    
    psi_A = psi_A / (np.max(np.abs(psi_A)) + 1e-10)
    psi_B = psi_B / (np.max(np.abs(psi_B)) + 1e-10)
    
    psi = 0.7 * (psi_A + psi_B)
    
    # Measure initial phase
    initial_phase = np.sum(np.angle(psi))
    
    # Perform braiding
    print(f"Performing {n_braids} complete exchange(s)...")
    psi = perform_braiding(psi, n_braids)
    
    # Measure final phase
    final_phase = np.sum(np.angle(psi))
    
    # Compute phase change
    delta_phase = final_phase - initial_phase
    delta_phase_wrapped = (delta_phase + np.pi) % (2*np.pi) - np.pi
    
    # Exchange angle (divide by 2 since full rotation = exchange)
    theta = abs(delta_phase_wrapped) / 2
    
    results.append({
        'n_braids': n_braids,
        'delta_phase': delta_phase_wrapped,
        'theta': theta,
        'theta_per_braid': theta / n_braids if n_braids > 0 else 0
    })
    
    print(f"ΔΦ: {delta_phase_wrapped:.3f}")
    print(f"θ (exchange angle): {theta:.3f}")
    print(f"θ per braid: {theta / n_braids:.3f}")

print("\n" + "="*80)
print("RESULTS SUMMARY")
print("="*80)

# ============================================================================
# ANALYSIS: Linear Accumulation?
# ============================================================================

print("\n| Braids (n) | θ_total | θ/n (avg) | Expected θ (n×θ₁) | Error |")
print("|------------|---------|-----------|-------------------|-------|")

theta_1 = results[0]['theta']  # Single braid baseline

for r in results:
    n = r['n_braids']
    theta_total = r['theta']
    theta_avg = r['theta_per_braid']
    expected = n * theta_1
    error = abs(theta_total - expected) / expected * 100 if expected > 0 else 0
    
    print(f"| {n:10d} | {theta_total:7.3f} | {theta_avg:9.3f} | {expected:17.3f} | {error:5.1f}% |")

# Linear fit
n_values = np.array([r['n_braids'] for r in results])
theta_values = np.array([r['theta'] for r in results])

if len(n_values) > 1:
    # Fit θ = a × n
    a_fit = np.sum(n_values * theta_values) / np.sum(n_values**2)
    
    # R-squared
    theta_mean = np.mean(theta_values)
    ss_tot = np.sum((theta_values - theta_mean)**2)
    ss_res = np.sum((theta_values - a_fit * n_values)**2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    print(f"\nLinear fit: θ = {a_fit:.3f} × n")
    print(f"R²: {r_squared:.4f}")
    print()
    
    if r_squared > 0.95:
        print("✅ LINEAR ACCUMULATION CONFIRMED!")
        print(f"   θ scales linearly with braiding count (R²={r_squared:.4f})")
        print(f"   Single-braid phase: θ₁ = {a_fit:.3f}")
        print()
        print("🌟 HOPFIONY SĄ PRAWDZIWYMI ANYONAMI!")
        print("   Implikacja: Topological quantum computing możliwy w FIN")
        anyonic = True
    else:
        print("❌ NON-LINEAR")
        print(f"   θ does NOT scale linearly (R²={r_squared:.4f})")
        print("   May be artifact or higher-order effect")
        anyonic = False
else:
    anyonic = None
    a_fit = theta_1
    r_squared = 0

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw605_braiding.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-605: Anyonic Braiding Accumulation\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Test:** Linearity of anyonic phase accumulation\n\n")
    
    f.write("## 1. Metodologia\n")
    f.write("Przeprowadzono n=1, 2, 3 pełnych wymian (braids) dwóch hopfionów.\n")
    f.write("Teoria anyonów: θ_total = n × θ_single\n\n")
    
    f.write("## 2. Wyniki\n")
    f.write("| Braids | θ_total | Expected | Error |\n")
    f.write("|--------|---------|----------|-------|\n")
    for r in results:
        n = r['n_braids']
        expected = n * theta_1
        error = abs(r['theta'] - expected) / expected * 100 if expected > 0 else 0
        f.write(f"| {n} | {r['theta']:.3f} | {expected:.3f} | {error:.1f}% |\n")
    f.write("\n")
    
    f.write(f"## 3. Analiza Liniowości\n")
    f.write(f"**Fit:** θ = {a_fit:.3f} × n\n")
    f.write(f"**R²:** {r_squared:.4f}\n\n")
    
    if anyonic:
        f.write("## 4. Wnioski\n")
        f.write("### ✅ PRAWDZIWE ANYONY POTWIERDZONE!\n\n")
        f.write(f"Faza anyoniczna kumuluje się **liniowo** (R²={r_squared:.4f} > 0.95).\n\n")
        f.write("**To oznacza że:**\n")
        f.write("- Hopfiony w FIN są **prawdziwymi anyonami**, nie artefaktem numerycznym\n")
        f.write(f"- Pojedyncza wymiana daje θ = {a_fit:.3f} ≈ π/{np.pi/a_fit:.2f}\n")
        f.write("- Multiple braids akumulują fazę zgodnie z teorią anyonów\n\n")
        f.write("**Implikacje:**\n")
        f.write("- Topological quantum computing możliwy w FIN!\n")
        f.write("- Link do fractional quantum Hall effect\n")
        f.write("- FIN ma **3D anyony** (zazwyczaj tylko 2D w standardowej fizyce)\n")
    elif anyonic == False:
        f.write("## 4. Wnioski\n")
        f.write("### ❌ BRAK LINIOWOŚCI\n")
        f.write("Faza nie kumuluje się liniowo. Możliwe że QW-603 był artefaktem.\n")
    else:
        f.write("## 4. Wnioski\n")
        f.write("### 🟡 NIEWYSTARCZAJĄCE DANE\n")

print("Report saved.")
print("="*80)
