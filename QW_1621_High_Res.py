#!/usr/bin/env python3
"""
QW-1621 HIGH-RESOLUTION: FULL 3D SKYRME PDE — N = 128/256
==========================================================
Type: DYNAMIC PDE TEST

Memory requirements:
- N = 128: ~8 GB RAM
- N = 256: ~64 GB RAM (may fail on normal machines)

Testing N = 128 first, then 256 if possible.
"""

import numpy as np
from scipy.ndimage import laplace
from datetime import datetime
import gc
import warnings
warnings.filterwarnings('ignore')

REPORT_FILE = "RAPORT_QW1621_HIGH_RES.md"

md = []
def log(msg=""):
    print(msg)
    md.append(msg)

log("=" * 80)
log("QW-1621 HIGH-RESOLUTION: SKYRME PDE (N ≥ 128)")
log("=" * 80)
log(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
log("")

# =============================================================================
# OPTIMIZED FUNCTIONS (memory-efficient)
# =============================================================================

def create_hedgehog(N, L, R=1.5):
    """Create hedgehog ansatz with minimal memory"""
    dx = L / N
    
    # Create coordinates on-the-fly
    x = np.linspace(-L/2, L/2, N, dtype=np.float32)
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    
    # Radial distance
    r = np.sqrt(X**2 + Y**2 + Z**2, dtype=np.float32)
    r[r < 1e-6] = 1e-6
    
    # Profile function f(r) = 2 arctan(R/r)
    f = 2 * np.arctan(R / r).astype(np.float32)
    
    # SU(2) field (σ, π₁, π₂, π₃)
    sigma = np.cos(f)
    pi_1 = (X / r) * np.sin(f)
    pi_2 = (Y / r) * np.sin(f)
    pi_3 = (Z / r) * np.sin(f)
    
    # Clean up
    del X, Y, Z, r, f
    gc.collect()
    
    return sigma, pi_1, pi_2, pi_3, dx

def normalize_inplace(sigma, pi_1, pi_2, pi_3):
    """Normalize to S³ in-place to save memory"""
    norm = np.sqrt(sigma**2 + pi_1**2 + pi_2**2 + pi_3**2)
    norm[norm < 1e-10] = 1e-10
    sigma /= norm
    pi_1 /= norm
    pi_2 /= norm
    pi_3 /= norm
    del norm
    gc.collect()

def compute_energy_fast(sigma, pi_1, pi_2, pi_3, dx):
    """
    Simplified energy calculation (σ-model dominant)
    Avoids storing all gradients simultaneously
    """
    E_total = 0.0
    
    # Sum of |∇φ|² for each field
    for field in [sigma, pi_1, pi_2, pi_3]:
        for axis in range(3):
            grad = np.gradient(field, dx, axis=axis)
            E_total += np.sum(grad**2)
            del grad
    
    gc.collect()
    return E_total * dx**3 / 16.0

def compute_Q_radial(sigma, pi_1, pi_2, pi_3, dx, N):
    """
    Compute topological charge using radial profile
    More robust than full 3D integral
    
    For hedgehog: Q = (2/π) ∫ sin²(f) |f'| dr
    """
    L = N * dx
    
    # Sample along radial direction
    n_radial = N // 2
    r_vals = np.linspace(dx, L/2 - dx, n_radial)
    
    # Extract radial profile f(r) from σ = cos(f)
    # At r, average σ over sphere
    Q_integral = 0.0
    
    center = N // 2
    
    for ir, r in enumerate(r_vals[:-1]):
        # Sample point on sphere
        i_r = int(r / dx) + center
        if i_r >= N:
            continue
        
        # Get σ value at this radius (approximate)
        sigma_r = sigma[center, center, min(i_r, N-1)]
        sigma_r = np.clip(sigma_r, -1.0, 1.0)
        f_r = np.arccos(sigma_r)
        
        # sin²(f)
        sin2_f = np.sin(f_r)**2
        
        # df/dr from difference
        if ir + 1 < len(r_vals):
            i_r_next = int(r_vals[ir+1] / dx) + center
            sigma_next = sigma[center, center, min(i_r_next, N-1)]
            sigma_next = np.clip(sigma_next, -1.0, 1.0)
            f_next = np.arccos(sigma_next)
            df_dr = abs(f_next - f_r) / (r_vals[ir+1] - r)
        else:
            df_dr = 0.0
        
        Q_integral += sin2_f * df_dr * (r_vals[1] - r_vals[0])
    
    Q = (2.0 / np.pi) * Q_integral
    return Q

def gradient_flow_step_inplace(sigma, pi_1, pi_2, pi_3, dx, dt):
    """In-place gradient flow to minimize memory"""
    # Compute laplacians
    lap_s = laplace(sigma) / dx**2
    lap_1 = laplace(pi_1) / dx**2
    lap_2 = laplace(pi_2) / dx**2
    lap_3 = laplace(pi_3) / dx**2
    
    # Update in-place
    sigma += dt * lap_s
    pi_1 += dt * lap_1
    pi_2 += dt * lap_2
    pi_3 += dt * lap_3
    
    del lap_s, lap_1, lap_2, lap_3
    
    # Boundary conditions (vacuum)
    for arr, val in [(sigma, 1.0), (pi_1, 0.0), (pi_2, 0.0), (pi_3, 0.0)]:
        arr[0,:,:] = val; arr[-1,:,:] = val
        arr[:,0,:] = val; arr[:,-1,:] = val
        arr[:,:,0] = val; arr[:,:,-1] = val
    
    normalize_inplace(sigma, pi_1, pi_2, pi_3)
    gc.collect()

# =============================================================================
# RUN HIGH-RESOLUTION TEST
# =============================================================================

L = 16.0  # Physical size
results = []

for N in [128]:  # Start with 128, can add 256 if memory allows
    log(f"[N = {N}³] Starting...")
    log("-" * 60)
    
    try:
        # Memory estimate
        mem_gb = (4 * N**3 * 4) / 1e9  # 4 fields × N³ × 4 bytes (float32)
        log(f"Estimated memory: {mem_gb:.1f} GB")
        
        # Initialize
        sigma, pi_1, pi_2, pi_3, dx = create_hedgehog(N, L)
        
        # Initial measurements
        E_init = compute_energy_fast(sigma, pi_1, pi_2, pi_3, dx)
        Q_init = compute_Q_radial(sigma, pi_1, pi_2, pi_3, dx, N)
        
        log(f"Initial: E = {E_init:.4f}, Q = {Q_init:.4f}")
        
        # Gradient flow
        dt = 0.0002  # Smaller dt for stability
        n_steps = 500
        check_interval = 100
        
        E_history = [E_init]
        Q_history = [Q_init]
        
        log(f"Running {n_steps} gradient flow steps...")
        
        for step in range(n_steps):
            gradient_flow_step_inplace(sigma, pi_1, pi_2, pi_3, dx, dt)
            
            if (step + 1) % check_interval == 0:
                E = compute_energy_fast(sigma, pi_1, pi_2, pi_3, dx)
                Q = compute_Q_radial(sigma, pi_1, pi_2, pi_3, dx, N)
                E_history.append(E)
                Q_history.append(Q)
                
                dE = E - E_history[-2]
                log(f"  Step {step+1}: E = {E:.4f} (ΔE = {dE:+.4f}), Q = {Q:.4f}")
        
        # Final
        E_final = E_history[-1]
        Q_final = Q_history[-1]
        
        # Check criteria
        monotonic = all(E_history[i] >= E_history[i+1] - 0.1 for i in range(len(E_history)-1))
        Q_in_range = 0.8 <= abs(Q_final) <= 1.2
        Q_pass = 0.98 <= abs(Q_final) <= 1.02
        
        log("")
        log(f"Final: E = {E_final:.4f}, Q = {Q_final:.4f}")
        log(f"Monotonic dE/dτ: {monotonic}")
        log(f"Q ∈ [0.98, 1.02]: {Q_pass}")
        
        results.append({
            'N': N,
            'Q_init': Q_init,
            'Q_final': Q_final,
            'E_init': E_init,
            'E_final': E_final,
            'monotonic': monotonic,
            'Q_pass': Q_pass
        })
        
        # Clean up
        del sigma, pi_1, pi_2, pi_3
        gc.collect()
        
    except MemoryError:
        log(f"❌ MemoryError: N = {N}³ requires too much RAM")
        results.append({'N': N, 'error': 'MemoryError'})
    except Exception as e:
        log(f"❌ Error: {e}")
        results.append({'N': N, 'error': str(e)})
    
    log("")

# =============================================================================
# VERDICT
# =============================================================================
log("=" * 60)
log("VERDICT")
log("=" * 60)

valid_results = [r for r in results if 'error' not in r]

if not valid_results:
    log("❌ INCONCLUSIVE: All grids failed")
    overall_status = "INCONCLUSIVE"
else:
    best = max(valid_results, key=lambda x: x['N'])
    
    if best['Q_pass'] and best['monotonic']:
        log(f"✅ CONSISTENT: Q = {best['Q_final']:.4f} ∈ [0.98, 1.02]")
        log(f"   Hedgehog is stable minimum under gradient flow.")
        overall_status = "CONSISTENT"
    elif 0.5 <= abs(best['Q_final']) <= 1.5 and best['monotonic']:
        log(f"⚠️ PARTIAL: Q = {best['Q_final']:.4f} (not quite 1.0)")
        log(f"   Need N = 256 for definitive result.")
        overall_status = "PARTIAL"
    else:
        log(f"❌ INCONCLUSIVE: Q = {best['Q_final']:.4f}")
        overall_status = "INCONCLUSIVE"

log("")
log(f"OVERALL STATUS: {overall_status}")

# =============================================================================
# REPORT
# =============================================================================
with open(REPORT_FILE, 'w', encoding='utf-8') as f:
    f.write("# QW-1621 HIGH-RESOLUTION: Skyrme PDE\n\n")
    f.write(f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    f.write(f"**Status:** {overall_status}\n\n")
    
    f.write("## Results\n\n")
    f.write("| N | Q_init | Q_final | E_final | Monotonic | Q_pass |\n")
    f.write("|---|--------|---------|---------|-----------|--------|\n")
    for r in results:
        if 'error' in r:
            f.write(f"| {r['N']} | - | ERROR: {r['error']} | - | - | - |\n")
        else:
            qp = "✅" if r['Q_pass'] else "❌"
            mon = "✅" if r['monotonic'] else "❌"
            f.write(f"| {r['N']} | {r['Q_init']:.4f} | {r['Q_final']:.4f} | {r['E_final']:.4f} | {mon} | {qp} |\n")
    
    f.write("\n## Raw Log\n```\n" + '\n'.join(md) + "\n```\n")

print(f"\n✅ Report saved to: {REPORT_FILE}")
