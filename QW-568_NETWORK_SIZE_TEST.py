#!/usr/bin/env python3
# QW-568: NETWORK SIZE TEST (N=1000 vs N=400)
# HYPOTHESIS: Larger network improves field measurements
# Re-testing QW-563 (velocity), QW-566 (dilation), QW-567 (dragging)
# Author: Based on QW-563-567 findings
# Date: 2025-12-05

import numpy as np
from scipy.linalg import expm
from scipy.spatial.distance import cdist
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

print("="*80)
print("QW-568: NETWORK SIZE TEST (N=1000)")
print("="*80)
print("HYPOTHESIS: N=1000 improves failed tests from QW-563/566/567")
print("Testing: Velocity field, Time dilation, Frame dragging")
print("="*80)

# FROZEN PARAMETERS (same as QW-563-567)
ALPHA_GEO = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01

GAMMA_GAIN = 1.0552
GAMMA_DAMP = 1.1980

print(f"\nFrozen Parameters (unchanged):")
print(f"  α_geo = {ALPHA_GEO:.6f}")
print(f"  ω = {OMEGA:.6f}")
print(f"  φ = {PHI:.6f}")
print(f"  β_tors = {BETA_TORS:.6f}")

def K_complex(d):
    """Complex kernel with phase"""
    return (ALPHA_GEO * np.exp(1j * (OMEGA * d + PHI))) / (1 + BETA_TORS * d)

# ============================================================================
# NETWORK SIZE COMPARISON
# ============================================================================

print("\n" + "="*80)
print("TESTING TWO NETWORK SIZES:")
print("="*80)

results_comparison = {}

for network_size in [400, 1000]:
    print(f"\n{'='*80}")
    print(f"NETWORK SIZE: N = {network_size}")
    print(f"{'='*80}")
    
    # Create network
    np.random.seed(568)  # Same seed for fair comparison!
    N_nodes = network_size
    positions_3d = np.random.randn(N_nodes, 3) * 3.0
    
    # Mass at origin
    mass_center_idx = np.argmin(np.linalg.norm(positions_3d, axis=1))
    print(f"Mass center: node {mass_center_idx}")
    
    # Build Hamiltonian
    print(f"Building {N_nodes}×{N_nodes} Hamiltonian...")
    dist_matrix = cdist(positions_3d, positions_3d)
    H = np.zeros((N_nodes, N_nodes), dtype=complex)
    
    for i in range(N_nodes):
        for j in range(i+1, N_nodes):
            d = dist_matrix[i, j]
            K_ij = K_complex(d)
            H[i, j] = K_ij
            H[j, i] = np.conj(K_ij)
    
    H = (H + H.conj().T) / 2
    print(f"Hamiltonian built: {H.shape}")
    
    # Initialize mass excitation
    psi = np.zeros(N_nodes, dtype=complex)
    for i in range(N_nodes):
        r = np.linalg.norm(positions_3d[i])
        if r < 1.0:
            psi[i] = np.exp(-r**2 / 0.5)
    psi = psi / np.linalg.norm(psi)
    
    # Evolve to steady state (shorter for speed)
    print(f"Evolving to steady state (50 steps)...")
    dt = 0.1
    for step in range(50):
        psi = expm(-1j * H * dt) @ psi
        psi = psi / np.linalg.norm(psi)
    
    print("Steady state reached")
    
    # ========================================================================
    # TEST 1: VELOCITY FIELD (QW-563 re-test)
    # ========================================================================
    
    print(f"\n{'-'*60}")
    print(f"TEST 1: Velocity Field Measurement (QW-563 re-test)")
    print(f"{'-'*60}")
    
    # Compute velocity from current
    velocities_radial = []
    radii = []
    
    for i in range(N_nodes):
        r_i = np.linalg.norm(positions_3d[i])
        if r_i < 0.5 or r_i > 7.0:
            continue
        
        flux = np.imag(np.conj(psi[i]) * np.dot(H[i, :], psi))
        r_hat = positions_3d[i] / (r_i + 1e-10)
        v_radial = -abs(flux) / (np.abs(psi[i])**2 + 1e-10)
        
        velocities_radial.append(v_radial)
        radii.append(r_i)
    
    velocities_radial = np.array(velocities_radial)
    radii = np.array(radii)
    
    # Sort
    sort_idx = np.argsort(radii)
    radii = radii[sort_idx]
    velocities_radial = velocities_radial[sort_idx]
    
    print(f"Sampled {len(radii)} nodes")
    
    # Fit v(r) = A / r^n
    def flow_law(r, A, n):
        return A / (r**n + 1e-10)
    
    try:
        popt, pcov = curve_fit(flow_law, radii, velocities_radial,
                               p0=[1.0, 0.5], bounds=([0, 0], [10, 2]))
        A_fit, n_fit = popt
        perr = np.sqrt(np.diag(pcov))
        
        v_fit = flow_law(radii, A_fit, n_fit)
        residuals = velocities_radial - v_fit
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((velocities_radial - np.mean(velocities_radial))**2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        error_from_gtr = abs(n_fit - 0.5) / 0.5 * 100
        
        print(f"  v(r) = {A_fit:.4f} / r^{n_fit:.4f}")
        print(f"  R² = {r_squared:.4f}")
        print(f"  GTR: n=0.5, Measured: n={n_fit:.4f}, Error: {error_from_gtr:.1f}%")
        
        velocity_success = (r_squared > 0.7 and error_from_gtr < 50)
        
    except:
        n_fit, r_squared, error_from_gtr = 0, 0, 999
        velocity_success = False
        print(f"  Fit failed")
    
    # ========================================================================
    # TEST 2: TIME DILATION (QW-566 re-test with corrections)
    # ========================================================================
    
    print(f"\n{'-'*60}")
    print(f"TEST 2: Time Dilation (QW-566 re-test, CORRECTED)")
    print(f"{'-'*60}")
    
    # Two SEPARATE wave functions for each clock
    r_clock_A = 1.5
    r_clock_B = 7.0
    
    pos_A = np.array([r_clock_A, 0.0, 0.0])
    pos_B = np.array([r_clock_B, 0.0, 0.0])
    clock_A_idx = np.argmin(np.linalg.norm(positions_3d - pos_A, axis=1))
    clock_B_idx = np.argmin(np.linalg.norm(positions_3d - pos_B, axis=1))
    
    print(f"  Clock A @ r={np.linalg.norm(positions_3d[clock_A_idx]):.2f}")
    print(f"  Clock B @ r={np.linalg.norm(positions_3d[clock_B_idx]):.2f}")
    
    # Initialize as small excitations (NOT single nodes)
    psi_clock_A = np.zeros(N_nodes, dtype=complex)
    psi_clock_B = np.zeros(N_nodes, dtype=complex)
    
    for i in range(N_nodes):
        r_to_A = np.linalg.norm(positions_3d[i] - positions_3d[clock_A_idx])
        r_to_B = np.linalg.norm(positions_3d[i] - positions_3d[clock_B_idx])
        
        if r_to_A < 0.5:
            psi_clock_A[i] = np.exp(-r_to_A**2 / 0.1)
        if r_to_B < 0.5:
            psi_clock_B[i] = np.exp(-r_to_B**2 / 0.1)
    
    psi_clock_A = psi_clock_A / np.linalg.norm(psi_clock_A)
    psi_clock_B = psi_clock_B / np.linalg.norm(psi_clock_B)
    
    # Evolve and measure FREQUENCIES
    n_steps = 100
    dt_clock = 0.05
    
    energy_A_list = []
    energy_B_list = []
    
    for step in range(n_steps):
        psi_clock_A = expm(-1j * H * dt_clock) @ psi_clock_A
        psi_clock_B = expm(-1j * H * dt_clock) @ psi_clock_B
        
        # Renormalize
        psi_clock_A = psi_clock_A / np.linalg.norm(psi_clock_A)
        psi_clock_B = psi_clock_B / np.linalg.norm(psi_clock_B)
        
        # Measure energy (frequency)
        E_A = np.real(psi_clock_A.conj() @ H @ psi_clock_A)
        E_B = np.real(psi_clock_B.conj() @ H @ psi_clock_B)
        
        energy_A_list.append(E_A)
        energy_B_list.append(E_B)
    
    # Frequency = mean energy
    freq_A = np.mean(energy_A_list[20:])  # After transient
    freq_B = np.mean(energy_B_list[20:])
    
    # Dilation ratio
    gamma_measured = freq_B / freq_A if freq_A > 0 else 1.0
    
    # GTR prediction (simplified)
    v_squared_A = 2.0 / (r_clock_A + 0.1)
    v_squared_B = 2.0 / (r_clock_B + 0.1)
    c_speed = 10.4
    gamma_theory = 1.0 + (v_squared_A - v_squared_B) / (2 * c_speed**2)
    
    error_dilation = abs(gamma_measured - gamma_theory) / abs(gamma_theory - 1.0) * 100 if abs(gamma_theory - 1.0) > 0.001 else 999
    
    print(f"  freq_A = {freq_A:.6f}")
    print(f"  freq_B = {freq_B:.6f}")
    print(f"  γ_measured = {gamma_measured:.6f}")
    print(f"  γ_theory = {gamma_theory:.6f}")
    print(f"  Error: {error_dilation:.1f}%")
    
    dilation_success = (error_dilation < 100 and gamma_measured > 1.0)
    
    # ========================================================================
    # TEST 3: FRAME DRAGGING (QW-567 re-test)
    # ========================================================================
    
    print(f"\n{'-'*60}")
    print(f"TEST 3: Frame Dragging (QW-567 re-test)")
    print(f"{'-'*60}")
    
    # Rotating vortex
    L_angular = 2.0
    psi_vortex = np.zeros(N_nodes, dtype=complex)
    
    for i in range(N_nodes):
        r_vec = positions_3d[i]
        r_mag = np.linalg.norm(r_vec)
        
        if r_mag < 1.5:
            theta = np.arctan2(r_vec[1], r_vec[0])
            amplitude = np.exp(-r_mag**2 / 0.5)
            phase = L_angular * theta
            psi_vortex[i] = amplitude * np.exp(1j * phase)
    
    psi_vortex = psi_vortex / np.linalg.norm(psi_vortex)
    
    # Evolve
    for step in range(100):
        psi_vortex = expm(-1j * H * 0.1) @ psi_vortex
        psi_vortex = psi_vortex / np.linalg.norm(psi_vortex)
    
    # Measure v_φ(r)
    r_bins = np.linspace(0.5, 6.0, 10)
    v_phi_binned = []
    r_bin_centers = []
    
    for i in range(len(r_bins) - 1):
        r_min, r_max = r_bins[i], r_bins[i+1]
        r_center = (r_min + r_max) / 2
        
        radii_all = np.linalg.norm(positions_3d, axis=1)
        in_shell = (radii_all >= r_min) & (radii_all < r_max)
        
        if np.sum(in_shell) < 5:
            continue
        
        phases_in_shell = np.angle(psi_vortex[in_shell])
        phase_std = np.std(phases_in_shell)
        
        v_phi_binned.append(phase_std * r_center)
        r_bin_centers.append(r_center)
    
    if len(r_bin_centers) > 4:
        v_phi_binned = np.array(v_phi_binned)
        r_bin_centers = np.array(r_bin_centers)
        
        def power_law_drag(r, A, n):
            return A / (r**n + 0.1)
        
        try:
            popt_drag, _ = curve_fit(power_law_drag, r_bin_centers, v_phi_binned,
                                      p0=[1.0, 2.0], bounds=([0, 0], [10, 5]))
            A_drag, n_drag = popt_drag
            
            v_fit = power_law_drag(r_bin_centers, A_drag, n_drag)
            ss_res = np.sum((v_phi_binned - v_fit)**2)
            ss_tot = np.sum((v_phi_binned - np.mean(v_phi_binned))**2)
            r_squared_drag = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            
            error_drag_exp = abs(n_drag - 2.0)
            
            print(f"  v_φ(r) = {A_drag:.4f} / r^{n_drag:.4f}")
            print(f"  R² = {r_squared_drag:.4f}")
            print(f"  GTR: n=2.0, Measured: n={n_drag:.2f}, Error: {error_drag_exp:.2f}")
            
            dragging_success = (error_drag_exp < 1.0 and r_squared_drag > 0.6)
            
        except:
            n_drag, r_squared_drag = 0, 0
            dragging_success = False
            print(f"  Fit failed")
    else:
        dragging_success = False
        print(f"  Insufficient data")
    
    # Store results
    results_comparison[network_size] = {
        'velocity': {'success': velocity_success, 'n': n_fit, 'R2': r_squared, 'error': error_from_gtr},
        'dilation': {'success': dilation_success, 'gamma': gamma_measured, 'error': error_dilation},
        'dragging': {'success': dragging_success, 'n': n_drag if 'n_drag' in locals() else 0,
                     'R2': r_squared_drag if 'r_squared_drag' in locals() else 0}
    }

# ============================================================================
# FINAL COMPARISON
# ============================================================================

print("\n" + "="*80)
print("NETWORK SIZE COMPARISON: N=400 vs N=1000")
print("="*80)

print(f"\n{'Test':<20} {'N=400':<20} {'N=1000':<20} {'Improvement?':<15}")
print("-"*80)

for test_name in ['velocity', 'dilation', 'dragging']:
    r400 = results_comparison[400][test_name]
    r1000 = results_comparison[1000][test_name]
    
    status_400 = "✅ PASS" if r400['success'] else "❌ FAIL"
    status_1000 = "✅ PASS" if r1000['success'] else "❌ FAIL"
    
    if r1000['success'] and not r400['success']:
        improvement = "✅ YES!"
    elif r1000['success'] and r400['success']:
        improvement = "✅ Both pass"
    elif not r1000['success'] and not r400['success']:
        improvement = "❌ Both fail"
    else:
        improvement = "⚠️ Worse"
    
    print(f"{test_name.capitalize():<20} {status_400:<20} {status_1000:<20} {improvement:<15}")

print("\n" + "="*80)
print("HYPOTHESIS TEST RESULT:")
print("="*80)

n400_passed = sum([results_comparison[400][t]['success'] for t in ['velocity', 'dilation', 'dragging']])
n1000_passed = sum([results_comparison[1000][t]['success'] for t in ['velocity', 'dilation', 'dragging']])

print(f"\nN=400:  {n400_passed}/3 tests passed")
print(f"N=1000: {n1000_passed}/3 tests passed")

if n1000_passed > n400_passed:
    print(f"\n✅ HYPOTHESIS CONFIRMED!")
    print(f"   N=1000 improves results: {n1000_passed-n400_passed} additional test(s) passed")
    print(f"   Larger network → smoother fields → better measurements")
elif n1000_passed == n400_passed and n1000_passed > 0:
    print(f"\n⚠️ HYPOTHESIS PARTIALLY CONFIRMED")
    print(f"   N=1000 doesn't improve pass rate, but may improve accuracy")
else:
    print(f"\n❌ HYPOTHESIS NOT CONFIRMED")
    print(f"   Network size alone doesn't fix the issues")
    print(f"   May need different method, not just larger N")

print("="*80)
print("Analysis complete.")
print("="*80)
