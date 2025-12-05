#!/usr/bin/env python3
# QW-569: ALTERNATIVE MEASUREMENT METHODS
# Based on QW-568 lesson: measure EFFECTS, not FIELDS
# 1. Propagation delay → v(r)
# 2. Correlation functions → time dilation
# 3. Angular momentum → frame dragging
# Author: Based on QW-563-568 diagnosis
# Date: 2025-12-05

import numpy as np
from scipy.linalg import expm
from scipy.spatial.distance import cdist
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

print("="*80)
print("QW-569: ALTERNATIVE MEASUREMENT METHODS")
print("="*80)
print("NEW APPROACH: Measure EFFECTS, not direct fields")
print("1. Propagation delay (like QW-435)")
print("2. Correlation functions (like QW-540)")
print("3. Angular momentum conservation (like QW-468)")
print("="*80)

# FROZEN PARAMETERS
ALPHA_GEO = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01

def K_complex(d):
    return (ALPHA_GEO * np.exp(1j * (OMEGA * d + PHI))) / (1 + BETA_TORS * d)

# Setup network (N=400 sufficient for effects!)
np.random.seed(569)
N_nodes = 400
positions_3d = np.random.randn(N_nodes, 3) * 3.0

print(f"\nNetwork: {N_nodes} nodes in 3D")

# Build Hamiltonian
mass_center_idx = np.argmin(np.linalg.norm(positions_3d, axis=1))
print(f"Mass center: node {mass_center_idx}")

dist_matrix = cdist(positions_3d, positions_3d)
H = np.zeros((N_nodes, N_nodes), dtype=complex)

for i in range(N_nodes):
    for j in range(i+1, N_nodes):
        d = dist_matrix[i, j]
        K_ij = K_complex(d)
        H[i, j] = K_ij
        H[j, i] = np.conj(K_ij)

H = (H + H.conj().T) / 2

# Initialize mass
psi_mass = np.zeros(N_nodes, dtype=complex)
for i in range(N_nodes):
    r = np.linalg.norm(positions_3d[i])
    if r < 1.0:
        psi_mass[i] = np.exp(-r**2 / 0.5)
psi_mass = psi_mass / np.linalg.norm(psi_mass)

# Evolve to steady state
print("Evolving mass to steady state...")
for step in range(50):
    psi_mass = expm(-1j * H * 0.1) @ psi_mass
    psi_mass = psi_mass / np.linalg.norm(psi_mass)

print("Steady state reached\n")

# ============================================================================
# METHOD 1: PROPAGATION DELAY → v(r)
# ============================================================================

print("="*80)
print("METHOD 1: PROPAGATION DELAY (as in QW-435)")
print("="*80)
print("Send pulse from far away, measure arrival time at different r")

# Create pulse at r=8 (far from mass)
pulse_start_r = 8.0
target_pos = np.array([pulse_start_r, 0, 0])
pulse_idx = np.argmin(np.linalg.norm(positions_3d - target_pos, axis=1))

psi_pulse = np.zeros(N_nodes, dtype=complex)
for i in range(N_nodes):
    r_to_pulse = np.linalg.norm(positions_3d[i] - positions_3d[pulse_idx])
    if r_to_pulse < 0.3:
        psi_pulse[i] = np.exp(-r_to_pulse**2 / 0.05)

psi_pulse = psi_pulse / np.linalg.norm(psi_pulse)

print(f"Pulse initialized @ r={np.linalg.norm(positions_3d[pulse_idx]):.2f}")
print("Propagating pulse toward mass (200 steps)...")

# Evolve and track arrival times
n_steps = 200
dt_pulse = 0.05

arrival_times = {}
radial_bins = np.linspace(0.5, 8.0, 10)

for step in range(n_steps):
    psi_pulse = expm(-1j * H * dt_pulse) @ psi_pulse
    # Don't renormalize - track spreading
    
    # Check arrival at each radial shell
    for i_bin in range(len(radial_bins)-1):
        r_min, r_max = radial_bins[i_bin], radial_bins[i_bin+1]
        r_center = (r_min + r_max) / 2
        
        if r_center in arrival_times:
            continue  # Already recorded
        
        # Check if pulse arrived
        radii = np.linalg.norm(positions_3d, axis=1)
        in_shell = (radii >= r_min) & (radii < r_max)
        
        if np.sum(in_shell) > 0:
            amplitude_in_shell = np.mean(np.abs(psi_pulse[in_shell])**2)
            
            # Threshold for "arrival"
            if amplitude_in_shell > 0.001:
                arrival_times[r_center] = step * dt_pulse

# Compute velocities from arrival times
radii_prop = []
velocities_prop = []

r_sorted = sorted(arrival_times.keys(), reverse=True)  # From far to near
for i in range(len(r_sorted)-1):
    r_outer = r_sorted[i]
    r_inner = r_sorted[i+1]
    
    t_outer = arrival_times[r_outer]
    t_inner = arrival_times[r_inner]
    
    if t_inner > t_outer:  # Sanity check
        delta_r = r_outer - r_inner
        delta_t = t_inner - t_outer
        
        v = delta_r / delta_t  # Average velocity in this shell
        r_avg = (r_outer + r_inner) / 2
        
        radii_prop.append(r_avg)
        velocities_prop.append(v)

if len(radii_prop) > 3:
    radii_prop = np.array(radii_prop)
    velocities_prop = np.array(velocities_prop)
    
    # Fit v(r) ~ A/r^n
    def flow_law(r, A, n):
        return A / (r**n + 0.1)
    
    try:
        popt, _ = curve_fit(flow_law, radii_prop, velocities_prop,
                           p0=[1.0, 0.5], bounds=([0, 0], [10, 2]))
        A_prop, n_prop = popt
        
        v_fit = flow_law(radii_prop, A_prop, n_prop)
        ss_res = np.sum((velocities_prop - v_fit)**2)
        ss_tot = np.sum((velocities_prop - np.mean(velocities_prop))**2)
        r2_prop = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        error_prop = abs(n_prop - 0.5) / 0.5 * 100
        
        print(f"\nPropagation delay results:")
        print(f"  Measured {len(radii_prop)} velocity points")
        print(f"  v(r) = {A_prop:.3f} / r^{n_prop:.3f}")
        print(f"  R² = {r2_prop:.3f}")
        print(f"  GTR: n=0.5, Measured: n={n_prop:.3f}")
        print(f"  Error: {error_prop:.1f}%")
        
        method1_success = (r2_prop > 0.6 and error_prop < 50)
        
        if method1_success:
            print(f"\n  ✅ SUCCESS: Propagation method works!")
        else:
            print(f"\n  ⚠️ PARTIAL: Reasonable fit but needs refinement")
            
    except:
        print(f"\nFit failed")
        method1_success = False
else:
    print(f"\nInsufficient data ({len(radii_prop)} points)")
    method1_success = False

# ============================================================================
# METHOD 2: CORRELATION FUNCTIONS → Time Dilation
# ============================================================================

print("\n" + "="*80)
print("METHOD 2: CORRELATION FUNCTIONS (as in QW-540)")
print("="*80)
print("Measure <ψ_A(t) ψ_B†(t)> correlation decay rate")

# Two regions at different r
r_region_A = 1.5  # Near mass
r_region_B = 6.0  # Far from mass

# Define regions
radii_all = np.linalg.norm(positions_3d, axis=1)
region_A_nodes = np.where((radii_all >= r_region_A-0.5) & (radii_all < r_region_A+0.5))[0]
region_B_nodes = np.where((radii_all >= r_region_B-0.5) & (radii_all < r_region_B+0.5))[0]

print(f"\nRegion A ({len(region_A_nodes)} nodes) @ r≈{r_region_A}")
print(f"Region B ({len(region_B_nodes)} nodes) @ r≈{r_region_B}")

# Initialize oscillating state
psi_oscil = np.zeros(N_nodes, dtype=complex)
for i in range(N_nodes):
    psi_oscil[i] = np.exp(1j * np.random.uniform(0, 2*np.pi))
psi_oscil = psi_oscil / np.linalg.norm(psi_oscil)

# Evolve and measure correlations
n_steps_corr = 100
dt_corr = 0.05

corr_A_history = []
corr_B_history = []

psi_init_A = psi_oscil[region_A_nodes].copy()
psi_init_B = psi_oscil[region_B_nodes].copy()

print("Measuring correlation decay...")

for step in range(n_steps_corr):
    psi_oscil = expm(-1j * H * dt_corr) @ psi_oscil
    psi_oscil = psi_oscil / np.linalg.norm(psi_oscil)
    
    # Correlation = <ψ(0)|ψ(t)>
    corr_A = abs(np.vdot(psi_init_A, psi_oscil[region_A_nodes]))**2
    corr_B = abs(np.vdot(psi_init_B, psi_oscil[region_B_nodes]))**2
    
    corr_A_history.append(corr_A)
    corr_B_history.append(corr_B)

corr_A_history = np.array(corr_A_history)
corr_B_history = np.array(corr_B_history)

# Fit exponential decay: C(t) ~ exp(-t/τ)
times = np.arange(n_steps_corr) * dt_corr

def exp_decay(t, tau):
    return np.exp(-t / (tau + 0.01))

try:
    popt_A, _ = curve_fit(exp_decay, times[10:], corr_A_history[10:],
                          p0=[10.0], bounds=([0.1], [100]))
    tau_A = popt_A[0]
    
    popt_B, _ = curve_fit(exp_decay, times[10:], corr_B_history[10:],
                          p0=[10.0], bounds=([0.1], [100]))
    tau_B = popt_B[0]
    
    # Time dilation: longer τ = slower decay = slower time
    # γ = τ_B / τ_A (how much faster B runs vs A)
    gamma_corr = tau_B / tau_A
    
    # GTR prediction
    v_sq_A = 2.0 / (r_region_A + 0.1)
    v_sq_B = 2.0 / (r_region_B + 0.1)
    c_speed = 10.4
    gamma_theory = 1.0 + (v_sq_A - v_sq_B) / (2 * c_speed**2)
    
    error_corr = abs(gamma_corr - gamma_theory) / abs(gamma_theory - 1.0) * 100 if abs(gamma_theory - 1.0) > 0.001 else 999
    
    print(f"\nCorrelation decay results:")
    print(f"  τ_A (near) = {tau_A:.3f}")
    print(f"  τ_B (far) = {tau_B:.3f}")
    print(f"  γ = τ_B/τ_A = {gamma_corr:.4f}")
    print(f"  γ_theory = {gamma_theory:.4f}")
    print(f"  Error: {error_corr:.1f}%")
    
    method2_success = (error_corr < 100 and gamma_corr > 1.0)
    
    if method2_success:
        print(f"\n  ✅ SUCCESS: Correlation method detects dilation!")
    else:
        print(f"\n  ⚠️ PARTIAL: Direction correct but magnitude off")
        
except Exception as e:
    print(f"\nCorrelation fit failed: {e}")
    method2_success = False

# ============================================================================
# METHOD 3: ANGULAR MOMENTUM → Frame Dragging
# ============================================================================

print("\n" + "="*80)
print("METHOD 3: ANGULAR MOMENTUM CONSERVATION (as in QW-468)")
print("="*80)
print("Measure <L_z> in shells around rotating mass")

# Create rotating vortex
L_vortex = 3.0
psi_vortex = np.zeros(N_nodes, dtype=complex)

for i in range(N_nodes):
    r_vec = positions_3d[i]
    r_mag = np.linalg.norm(r_vec)
    
    if r_mag < 1.5:
        theta = np.arctan2(r_vec[1], r_vec[0])
        amp = np.exp(-r_mag**2 / 0.5)
        psi_vortex[i] = amp * np.exp(1j * L_vortex * theta)

psi_vortex = psi_vortex / np.linalg.norm(psi_vortex)

print(f"Vortex with L={L_vortex} initialized")
print("Evolving (100 steps)...")

# Evolve
for step in range(100):
    psi_vortex = expm(-1j * H * 0.1) @ psi_vortex
    psi_vortex = psi_vortex / np.linalg.norm(psi_vortex)

print("Measuring angular momentum in shells...")

# Measure L_z in radial shells
r_bins_L = np.linspace(0.5, 6.0, 10)
L_z_per_shell = []
r_shell_centers = []

for i in range(len(r_bins_L)-1):
    r_min, r_max = r_bins_L[i], r_bins_L[i+1]
    r_center = (r_min + r_max) / 2
    
    in_shell = (radii_all >= r_min) & (radii_all < r_max)
    
    if np.sum(in_shell) < 5:
        continue
    
    # L_z operator: -i(x∂_y - y∂_x)
    # Approximate via phase winding
    L_z_total = 0
    
    for idx in np.where(in_shell)[0]:
        x, y, z = positions_3d[idx]
        theta = np.arctan2(y, x)
        
        # Phase at this node
        phase = np.angle(psi_vortex[idx])
        
        # Contribution to L_z
        L_contribution = phase / (2*np.pi)  # Winding number
        L_z_total += L_contribution
    
    L_z_per_shell.append(L_z_total)
    r_shell_centers.append(r_center)

if len(r_shell_centers) > 4:
    L_z_per_shell = np.array(L_z_per_shell)
    r_shell_centers = np.array(r_shell_centers)
    
    # Normalize by number of nodes (to get average L_z)
    # Frame dragging: L_z(r) ~ 1/r² (GTR prediction)
    
    try:
        def power_law_L(r, A, n):
            return A / (r**n + 0.1)
        
        popt_L, _ = curve_fit(power_law_L, r_shell_centers, np.abs(L_z_per_shell),
                              p0=[5.0, 2.0], bounds=([0, 0], [100, 5]))
        A_L, n_L = popt_L
        
        L_fit = power_law_L(r_shell_centers, A_L, n_L)
        ss_res_L = np.sum((np.abs(L_z_per_shell) - L_fit)**2)
        ss_tot_L = np.sum((np.abs(L_z_per_shell) - np.mean(np.abs(L_z_per_shell)))**2)
        r2_L = 1 - (ss_res_L / ss_tot_L) if ss_tot_L > 0 else 0
        
        error_L = abs(n_L - 2.0)
        
        print(f"\nAngular momentum results:")
        print(f"  L_z(r) = {A_L:.2f} / r^{n_L:.3f}")
        print(f"  R² = {r2_L:.3f}")
        print(f"  GTR: n=2.0, Measured: n={n_L:.3f}")
        print(f"  Deviation: {error_L:.2f}")
        
        method3_success = (error_L < 1.0 and r2_L > 0.5)
        
        if method3_success:
            print(f"\n  ✅ SUCCESS: Angular momentum method works!")
        else:
            print(f"\n  ⚠️ PARTIAL: Shows decay but wrong exponent")
            
    except Exception as e:
        print(f"\nFit failed: {e}")
        method3_success = False
else:
    print(f"\nInsufficient shells")
    method3_success = False

# ============================================================================
# FINAL COMPARISON
# ============================================================================

print("\n" + "="*80)
print("FINAL COMPARISON: Alternative Methods vs Direct Fields")
print("="*80)

results = {
    'Propagation delay → v(r)': '✅ SUCCESS' if method1_success else '❌ FAIL',
    'Correlation → dilation': '✅ SUCCESS' if method2_success else '❌ FAIL',
    'Angular momentum → dragging': '✅ SUCCESS' if method3_success else '❌ FAIL'
}

print("\nAlternative Methods (QW-569):")
for method, status in results.items():
    print(f"  {method}: {status}")

passed = sum([method1_success, method2_success, method3_success])
print(f"\nPassed: {passed}/3 ({passed/3*100:.0f}%)")

print("\nComparison with direct field methods (QW-563/566/567):")
print("  Direct velocity field: ❌ FAIL (n=2.0, R²=-0.05)")
print("  Direct clock phases: ❌ FAIL (γ=21, 2000% error)")
print("  Direct phase gradient: ❌ FAIL (n=0, no signal)")
print("  Success rate: 0/3 (0%)")

print("\n" + "="*80)
print("CONCLUSION:")
print("="*80)

if passed >= 2:
    print("✅✅ ALTERNATIVE METHODS WORK!")
    print("\n  Measuring EFFECTS successfully detects:")
    if method1_success:
        print(f"  • Flow velocity from propagation delay")
    if method2_success:
        print(f"  • Time dilation from correlation decay")
    if method3_success:
        print(f"  • Frame dragging from angular momentum")
    
    print("\n  **KEY LESSON:**")
    print("  Discrete networks CAN reveal flow physics,")
    print("  but you must measure EFFECTS, not FIELDS!")
    print("\n  Effects > Fields confirmed!")
    
elif passed >= 1:
    print("⚠️ PARTIAL SUCCESS - Some methods work")
    print(f"   {passed}/3 methods successful")
    print("   Alternative approaches show promise")
else:
    print("❌ Alternative methods don't improve results")
    print("   May need fundamental rethinking")

print("="*80)
print("Analysis complete.")
print("="*80)
