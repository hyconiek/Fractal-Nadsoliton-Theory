#!/usr/bin/env python3
# QW-590: HOPFION STABILITY AT LAYER N=10
# Purpose: Test Hypothesis 4 (Particles as Topological Vortices) at the correct fractal scale.
# Previous attempts at N=1 failed. This simulates N=10 with hierarchical protection.
# Date: 2025-12-05

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import convolve

print("="*80)
print("QW-590: HOPFION STABILITY AT LAYER N=10")
print("="*80)

# ============================================================================
# 1. PARAMETERS
# ============================================================================
GRID_SIZE = 32
DT = 0.05
STEPS = 1000
BETA_TORS = 0.01
N_LAYER = 10

# Scale parameters for Layer N=10
# In relative units, we simulate the dynamics *at* this layer.
# The influence of other layers is modeled as an effective potential.

print(f"Grid: {GRID_SIZE}x{GRID_SIZE}x{GRID_SIZE}")
print(f"Layer: N={N_LAYER} (Proton Scale)")
print(f"Beta: {BETA_TORS}")
print("-" * 40)

# ============================================================================
# 2. FIELD INITIALIZATION (HOPFION)
# ============================================================================

def hopfion_field(grid_size, R=8.0):
    """
    Initialize a Hopfion (knot) configuration in 3D.
    Maps S3 -> S2.
    """
    x = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    y = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    z = np.linspace(-GRID_SIZE/2, GRID_SIZE/2, GRID_SIZE)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    r = np.sqrt(X**2 + Y**2 + Z**2)
    rho = np.sqrt(X**2 + Y**2)
    
    # Avoid division by zero
    rho[rho == 0] = 1e-10
    
    # Toroidal coordinates
    # eta: angle in the poloidal plane (around the torus tube)
    eta = np.arctan2(Z, rho - R)
    # xi: angle in the toroidal plane (around the z-axis)
    xi = np.arctan2(Y, X)
    
    # Profile function f(r) - goes from 0 at core to pi at infinity
    # Simple ansatz: linear core
    f = np.pi * np.tanh(r / R)
    
    # Spinor construction (simplified map)
    # Using the standard Hopf map construction via complex coordinates
    # Z0 = cos(f) + i sin(f) cos(theta) ... actually let's use the explicit form
    
    # Explicit Hopfion form with winding index = 1
    # Psi is a complex scalar field carrying the phase info? 
    # Wait, Hopfion usually requires a 2-component spinor or a vector field.
    # In FIN, we have a complex scalar field Psi.
    # A single complex scalar can encode a Vortex (winding 2D), but a Hopfion (knot) 
    # requires mapping to S2, which usually implies a spinor Psi = (u, v)^T -> n = z^dag sigma z.
    
    # HOWEVER, the FIN Master Equation is for a scalar Psi.
    # Hypothesis 4 says "Particles are Topological Vortices".
    # A scalar field in 3D can have vortex *lines* (where Psi=0).
    # A Hopfion in a scalar field would be a *closed* vortex loop (knot).
    
    print("Initializing Scalar Hopfion (Vortex Loop/Unknot)...")
    
    # We create a vortex ring (unknot) in the scalar field.
    # Phase winds around the ring.
    
    # Toroidal phase: phi = xi + eta
    # This creates a twisted phase structure.
    phase = xi + eta
    
    # Amplitude: 0 at the core (r=R ring), 1 far away
    # Core is at rho = R, Z = 0
    dist_to_ring = np.sqrt((rho - R)**2 + Z**2)
    amplitude = np.tanh(dist_to_ring / 2.0)
    
    psi = amplitude * np.exp(1j * phase)
    
    return psi

# Initialize
psi = hopfion_field(GRID_SIZE, R=GRID_SIZE/4.0)

# ============================================================================
# 3. TOPOLOGICAL CHARGE MEASUREMENT
# ============================================================================

def get_derivatives(psi):
    # Central differences
    dx = np.roll(psi, -1, axis=0) - np.roll(psi, 1, axis=0)
    dy = np.roll(psi, -1, axis=1) - np.roll(psi, 1, axis=1)
    dz = np.roll(psi, -1, axis=2) - np.roll(psi, 1, axis=2)
    return dx, dy, dz

def compute_hopf_invariant(psi):
    """
    Compute Hopf invariant Q for a scalar field vortex loop?
    For a scalar field, the topological invariant is the winding number around the loop, 
    or the 'Helicity' of the velocity field.
    
    Let v = Im(grad psi / psi) = gradient of phase
    Helicity H = Integral v . curl v dV
    
    (Note: for a pure gradient field v=grad phi, curl v = 0, except at singularities.
     The helicity captures the knotting of the vortex lines.)
    """
    # 1. Calculate phase gradient (superfluid velocity)
    # v_s = Im( (psi* . grad psi) / |psi|^2 )
    
    psi_conj = np.conj(psi)
    rho_sq = np.abs(psi)**2 + 1e-10
    
    grad_x, grad_y, grad_z = get_derivatives(psi)
    
    v_x = np.imag(psi_conj * grad_x) / rho_sq
    v_y = np.imag(psi_conj * grad_y) / rho_sq
    v_z = np.imag(psi_conj * grad_z) / rho_sq
    
    # 2. Calculate vorticity w = curl v
    # w_x = d(v_z)/dy - d(v_y)/dz
    
    dv_z_dy = np.roll(v_z, -1, axis=1) - np.roll(v_z, 1, axis=1)
    dv_y_dz = np.roll(v_y, -1, axis=2) - np.roll(v_y, 1, axis=2)
    w_x = dv_z_dy - dv_y_dz
    
    dv_x_dz = np.roll(v_x, -1, axis=2) - np.roll(v_x, 1, axis=2)
    dv_z_dx = np.roll(v_z, -1, axis=0) - np.roll(v_z, 1, axis=0)
    w_y = dv_x_dz - dv_z_dx
    
    dv_y_dx = np.roll(v_y, -1, axis=0) - np.roll(v_y, 1, axis=0)
    dv_x_dy = np.roll(v_x, -1, axis=1) - np.roll(v_x, 1, axis=1)
    w_z = dv_y_dx - dv_x_dy
    
    # 3. Calculate H = v . w
    helicity_density = v_x * w_x + v_y * w_y + v_z * w_z
    
    # Normalize by volume factor? We just want relative stability.
    H = np.sum(helicity_density)
    
    return H

initial_Q = compute_hopf_invariant(psi)
print(f"Initial Helicity (approx Q): {initial_Q:.4f}")

# ============================================================================
# 4. EVOLUTION DYNAMICS
# ============================================================================

# Laplacian convolution kernel for 3D
laplacian_kernel = np.zeros((3,3,3))
laplacian_kernel[1,1,1] = -6
for i in [(0,1,1), (2,1,1), (1,0,1), (1,2,1), (1,1,0), (1,1,2)]:
    laplacian_kernel[i] = 1

def laplacian(field):
    return convolve(field, laplacian_kernel, mode='wrap')

# Effective Potential representing Hierarchy
def V_eff(psi):
    # 1. Local Nonlinearity (Self-interaction)
    # G * |psi|^2 * psi
    # Standard NLS term
    
    # 2. Hierarchical Pressure (The new "Knowledge")
    # Deeper layers (N<10) act as a confining potential or "mass" term
    # V_deep ~ sum(Beta^(N-M)) ~ Beta^1 + Beta^2...
    # For N=10, dominant term is N=9 with weight 0.01
    
    rho = np.abs(psi)**2
    
    V_local = rho # Standard
    
    # N=10 logic:
    # "Hard" layers below push back -> positive potential term
    V_deep = 0.01 * rho**2  # Higher order correction? Or just modifying effective mass?
    
    # "Soft" layers above allow dissipation -> imaginary loss?
    # No, keep it unitary first.
    
    # Effective Mass renormalization from QW-589
    # m_eff ~ rho
    # This implies the coefficient of the Laplacian depends on rho?
    # Let's stick to the potential formulation first.
    
    return V_local + V_deep

print("Simulating evolution with Hierarchical Effective Potential...")

history_Q = []
history_max_amp = []

for t in range(STEPS):
    # Kinetic term (Laplacian)
    # i dPsi/dt = -Laplacian Psi + V Psi
    
    kin = -laplacian(psi)
    pot = V_eff(psi) * psi
    
    dpsi_dt = -1j * (kin + pot)
    
    # Correction: Dissipation needed to find ground state/attractor?
    # User said: "Nadsoliton is an Attractor".
    # Pure Hamiltonian evolution (i) conserves energy and doesn't converge to attractor.
    # We need a non-unitary term (dissipation/gain) like in QW-558.
    
    # Attractor dynamics: dPsi/dt = (Gain - Loss) Psi + i Hamiltonian
    # Loss ~ Beta_tors (viscosity)
    # Gain ~ Nonlinearity
    
    # Let's add the dissipative term:
    # dPsi/dt = (gamma - delta |psi|^2) psi - i [Hamiltonian]
    # gamma = linear gain/pumping (from vacuum energy?)
    # delta = nonlinear saturation
    
    # "Universal Attractor" term from QW-558:
    # dA/dt = A - A^3
    
    gamma = 0.1
    dissip = gamma * psi * (1.0 - np.abs(psi)**2)
    
    # Total update
    change = dpsi_dt + dissip
    
    psi += DT * change
    
    # Measure
    if t % 50 == 0:
        Q = compute_hopf_invariant(psi)
        max_amp = np.max(np.abs(psi))
        history_Q.append(Q)
        history_max_amp.append(max_amp)
        print(f"Step {t:4d}: Q={Q:.4f}, MaxAmp={max_amp:.4f}")
        
        # Check stability
        if np.abs(Q) < 0.1 * np.abs(initial_Q):
            print("❌ Topology decayed (Q -> 0)")
            break
            
print("-" * 40)
print("Simulation Complete.")

# ============================================================================
# 5. ANALYSIS
# ============================================================================

final_Q = history_Q[-1] if history_Q else 0
ratio = final_Q / initial_Q if initial_Q != 0 else 0

print(f"Initial Q: {initial_Q:.4f}")
print(f"Final Q:   {final_Q:.4f}")
print(f"Retention: {ratio*100:.1f}%")

print()
if ratio > 0.8:
    print("✅ RESULT: HOPFION IS STABLE at Layer N=10!")
    print("   Hypothesis 4 CONFIRMED via dynamic attractor.")
elif ratio > 0.3:
    print("🟡 RESULT: Partial Stability (Metastable).")
    print("   Needs fine-tuning, but better than N=1 decay.")
else:
    print("❌ RESULT: Hopfion Decayed.")
    print("   Hypothesis 4 FALSIFIED even at N=10.")

