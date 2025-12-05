#!/usr/bin/env python3
# QW-570: FRAME DRAGGING ENHANCED (N=1500 + Active Driving)
# Replicating QW-405 success methodology on discrete network
# 1. Larger network (N=1500)
# 2. Active driving (forcing rotation in mass)
# 3. Angular momentum transfer measurement
# Author: Based on QW-405 analysis
# Date: 2025-12-05

import numpy as np
from scipy.linalg import expm
from scipy.spatial.distance import cdist
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import time

print("="*80)
print("QW-570: FRAME DRAGGING ENHANCED")
print("="*80)
print("Methodology from QW-405 (Success):")
print("1. Active driving (Omega * L_z term in Hamiltonian)")
print("2. Larger network (N=1500) for smoother gradients")
print("3. Measure angular momentum transfer to vacuum")
print("="*80)

# FROZEN PARAMETERS
ALPHA_GEO = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01

def K_complex(d):
    return (ALPHA_GEO * np.exp(1j * (OMEGA * d + PHI))) / (1 + BETA_TORS * d)

# ============================================================================
# 1. SETUP LARGE NETWORK (N=1500)
# ============================================================================

print(f"\nInitializing N=1500 network...")
start_time = time.time()

np.random.seed(570)
N_nodes = 1500
positions_3d = np.random.randn(N_nodes, 3) * 4.0  # Slightly larger volume

# Mass at origin
mass_center_idx = np.argmin(np.linalg.norm(positions_3d, axis=1))
radii = np.linalg.norm(positions_3d, axis=1)

# Define mass region (nodes to be driven)
mass_radius = 1.5
mass_indices = np.where(radii < mass_radius)[0]
vacuum_indices = np.where(radii >= mass_radius)[0]

print(f"  Nodes in mass region (r < {mass_radius}): {len(mass_indices)}")
print(f"  Nodes in vacuum region: {len(vacuum_indices)}")

# Build Base Hamiltonian (Interaction)
print("Building Hamiltonian (this may take a moment)...")
dist_matrix = cdist(positions_3d, positions_3d)
H_base = np.zeros((N_nodes, N_nodes), dtype=complex)

# Optimization: Vectorized kernel
K_matrix = K_complex(dist_matrix)
np.fill_diagonal(K_matrix, 0)  # No self-interaction
H_base = (K_matrix + K_matrix.conj().T) / 2

print(f"Hamiltonian built in {time.time() - start_time:.1f}s")

# ============================================================================
# 2. ADD ACTIVE DRIVING (ROTATION TERM)
# ============================================================================

print("\nAdding active rotation term (L_z) to mass nodes...")

# Construct L_z operator for discrete nodes
# L_z ~ -i (x ∂y - y ∂x)
# In discrete graph: Phase gradient along azimuthal direction
# We approximate L_z coupling between neighbors

Omega_rot = 5.0  # Driving strength
H_rot = np.zeros((N_nodes, N_nodes), dtype=complex)

# For each node in mass region, couple to neighbors with phase shift
for i in mass_indices:
    xi, yi, zi = positions_3d[i]
    ri = np.sqrt(xi**2 + yi**2)
    theta_i = np.arctan2(yi, xi)
    
    # Find neighbors
    neighbors = np.where(dist_matrix[i] < 0.8)[0]
    
    for j in neighbors:
        if i == j: continue
        
        xj, yj, zj = positions_3d[j]
        theta_j = np.arctan2(yj, xj)
        
        # Delta theta (shortest path)
        d_theta = theta_j - theta_i
        if d_theta > np.pi: d_theta -= 2*np.pi
        if d_theta < -np.pi: d_theta += 2*np.pi
        
        # Coupling proportional to d_theta (angular momentum transfer)
        # H_rot[i,j] adds energy if moving against rotation
        coupling = 1j * Omega_rot * d_theta
        H_rot[i, j] += coupling

# Total Hamiltonian: H = H_base + H_rot
# Only drive the mass nodes!
H_total = H_base.copy()
H_total[np.ix_(mass_indices, mass_indices)] += H_rot[np.ix_(mass_indices, mass_indices)]

# Ensure Hermiticity
H_total = (H_total + H_total.conj().T) / 2

print("Rotation term added.")

# ============================================================================
# 3. EVOLUTION
# ============================================================================

print("\nEvolving system (Active Driving)...")

# Initial state: Uniform vacuum (no rotation)
psi = np.ones(N_nodes, dtype=complex)
psi = psi / np.linalg.norm(psi)

# Measure initial angular momentum in vacuum
def measure_Lz_vacuum(psi_curr):
    L_z_total = 0
    count = 0
    
    for i in vacuum_indices:
        # Local phase gradient estimate
        neighbors = np.where(dist_matrix[i] < 0.8)[0]
        if len(neighbors) < 2: continue
        
        xi, yi, _ = positions_3d[i]
        theta_i = np.arctan2(yi, xi)
        
        phase_i = np.angle(psi_curr[i])
        
        # Average d_phase/d_theta around node
        d_phi_d_theta = 0
        w_sum = 0
        
        for j in neighbors:
            xj, yj, _ = positions_3d[j]
            theta_j = np.arctan2(yj, xj)
            
            d_theta = theta_j - theta_i
            if d_theta > np.pi: d_theta -= 2*np.pi
            if d_theta < -np.pi: d_theta += 2*np.pi
            
            if abs(d_theta) < 0.1: continue # Too close in angle
            
            phase_j = np.angle(psi_curr[j])
            d_phase = phase_j - phase_i
            if d_phase > np.pi: d_phase -= 2*np.pi
            if d_phase < -np.pi: d_phase += 2*np.pi
            
            # L_z ~ d_phase / d_theta
            est = d_phase / d_theta
            weight = np.exp(-dist_matrix[i,j])
            
            d_phi_d_theta += est * weight
            w_sum += weight
            
        if w_sum > 0:
            L_z_total += d_phi_d_theta / w_sum
            count += 1
            
    return L_z_total / count if count > 0 else 0

Lz_init = measure_Lz_vacuum(psi)
print(f"Initial Vacuum L_z: {Lz_init:.6f}")

# Evolve
n_steps = 500
dt = 0.05
history_Lz = []

for step in range(n_steps):
    psi = expm(-1j * H_total * dt) @ psi
    psi = psi / np.linalg.norm(psi)
    
    if step % 50 == 0:
        Lz_curr = measure_Lz_vacuum(psi)
        history_Lz.append(Lz_curr)
        print(f"  Step {step}: Vacuum L_z = {Lz_curr:.6f}")

Lz_final = measure_Lz_vacuum(psi)
print(f"Final Vacuum L_z: {Lz_final:.6f}")
print(f"Change ΔL_z: {Lz_final - Lz_init:.6f}")

# ============================================================================
# 4. ANALYSIS: ROTATION CURVE
# ============================================================================

print("\nAnalyzing Rotation Curve v_φ(r)...")

# Binning
r_bins = np.linspace(mass_radius, 6.0, 12)
v_phi_binned = []
r_centers = []

for i in range(len(r_bins)-1):
    r_min, r_max = r_bins[i], r_bins[i+1]
    r_c = (r_min + r_max) / 2
    
    mask = (radii >= r_min) & (radii < r_max)
    nodes_in_shell = np.where(mask)[0]
    
    if len(nodes_in_shell) < 5: continue
    
    # Calculate average azimuthal velocity in shell
    v_phi_sum = 0
    count = 0
    
    for idx in nodes_in_shell:
        # Same L_z estimator
        neighbors = np.where(dist_matrix[idx] < 0.8)[0]
        if len(neighbors) < 2: continue
        
        xi, yi, _ = positions_3d[idx]
        theta_i = np.arctan2(yi, xi)
        phase_i = np.angle(psi[idx])
        
        d_phi_d_theta = 0
        w_sum = 0
        
        for j in neighbors:
            xj, yj, _ = positions_3d[j]
            theta_j = np.arctan2(yj, xj)
            d_theta = theta_j - theta_i
            if d_theta > np.pi: d_theta -= 2*np.pi
            if d_theta < -np.pi: d_theta += 2*np.pi
            if abs(d_theta) < 0.1: continue
            
            phase_j = np.angle(psi[j])
            d_phase = phase_j - phase_i
            if d_phase > np.pi: d_phase -= 2*np.pi
            if d_phase < -np.pi: d_phase += 2*np.pi
            
            est = d_phase / d_theta
            weight = np.exp(-dist_matrix[idx,j])
            d_phi_d_theta += est * weight
            w_sum += weight
            
        if w_sum > 0:
            L_z_local = d_phi_d_theta / w_sum
            # v_phi = L_z / r
            v_phi_sum += L_z_local / r_c
            count += 1
            
    if count > 0:
        v_phi_binned.append(v_phi_sum / count)
        r_centers.append(r_c)

if len(r_centers) > 3:
    r_centers = np.array(r_centers)
    v_phi_binned = np.array(v_phi_binned)
    
    # Fit v(r) ~ 1/r^n
    def drag_law(r, A, n):
        return A / (r**n + 0.1)
    
    try:
        popt, _ = curve_fit(drag_law, r_centers, np.abs(v_phi_binned),
                           p0=[1.0, 2.0], bounds=([0, 0], [10, 5]))
        A_fit, n_fit = popt
        
        v_fit = drag_law(r_centers, A_fit, n_fit)
        ss_res = np.sum((np.abs(v_phi_binned) - v_fit)**2)
        ss_tot = np.sum((np.abs(v_phi_binned) - np.mean(np.abs(v_phi_binned)))**2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        print(f"\nRotation Curve Results:")
        print(f"  v_φ(r) = {A_fit:.4f} / r^{n_fit:.4f}")
        print(f"  R² = {r2:.4f}")
        print(f"  GTR Prediction: n ≈ 2.0 (1/r²)")
        print(f"  Measured: n = {n_fit:.4f}")
        
        if r2 > 0.6 and abs(n_fit - 2.0) < 1.0:
            print(f"\n✅ SUCCESS: Frame Dragging Detected!")
            print(f"   Vacuum acquired rotation with correct profile")
        elif r2 > 0.5:
            print(f"\n⚠️ PARTIAL: Rotation detected but profile deviates")
        else:
            print(f"\n❌ FAIL: No clear rotation profile")
            
    except Exception as e:
        print(f"Fit failed: {e}")
else:
    print("Insufficient data for curve fit")

print("="*80)
print("QW-570 ANALYSIS COMPLETE")
print("="*80)
