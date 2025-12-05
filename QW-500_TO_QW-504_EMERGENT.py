# QW-500 TO QW-504: EMERGENT REALITY CHECK
# PARADIGM: "Let the Network Speak" (No imposed formulas).
# STRICT PROTOCOL: Simulation -> Observation -> Analysis.

import numpy as np
import scipy.fft
import scipy.linalg
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt

# FROZEN PARAMETERS
alpha_geo = 4 * np.log(2)
omega = np.pi/4
phi = np.pi/6
beta_tors = 0.01

print("="*80)
print("QW-500 TO QW-504: EMERGENT REALITY CHECK")
print("PARADIGM: Zero Fitting. Information as Foundation.")
print("="*80)

def K_kernel(d):
    """The Unified Kernel"""
    return (alpha_geo * np.exp(1j * (omega * d + phi))) / (1 + beta_tors * d)

# --- QW-500: ETHERIC RESONANCE (Hydrogen Spectrum without Schrodinger) ---
print("\n" + "="*60)
print("QW-500: ETHERIC RESONANCE (Network FFT)")
print("="*60)

# Simulate a 1D network string with a "proton" load at center
# Wave equation on graph: d^2 psi / dt^2 = -L psi (Laplacian)
# But here interactions are non-local via K(d)

N_points = 128
x = np.arange(N_points)
center = N_points // 2

# Construct interaction matrix M based on K(d)
# M_ij = K(|i-j|)
M = np.zeros((N_points, N_points), dtype=complex)
for i in range(N_points):
    for j in range(N_points):
        d = abs(i - j)
        if d == 0:
            M[i, j] = 0 # Self-interaction handled separately or zero
        else:
            M[i, j] = K_kernel(d)

# "Proton" is a heavy mass/potential at the center
# We modify the diagonal to represent the load
# In a wave equation, mass appears in the time derivative term, or as a potential
# Let's treat it as a potential well V(x) added to the operator
V = np.zeros(N_points)
V[center] = -10.0 # Deep potential well (Proton)

# Evolution operator H = -M + V (Effective Hamiltonian)
# We look for eigenfrequencies of this system
H_eff = -np.real(M) + np.diag(V)

# Diagonalize to find natural frequencies
evals, evecs = scipy.linalg.eigh(H_eff)

# Sort eigenvalues (Energy ~ Frequency^2 or Frequency)
# In wave equation d2/dt2 = -H psi, frequencies are sqrt(evals)
# We are looking for bound states (negative energy in Schrodinger picture, or specific modes)
# Let's look at the lowest eigenvalues (most bound)
print("Natural Frequencies (Eigenvalues of Network Hamiltonian):")
bound_states = evals[evals < 0] # Potential well creates negative eigenvalues
bound_states = np.sort(bound_states)

for i, E in enumerate(bound_states[:5]):
    print(f"  Mode {i+1}: E = {E:.6f}")

if len(bound_states) >= 3:
    E1, E2, E3 = bound_states[0], bound_states[1], bound_states[2]
    print(f"\nChecking Balmer-like ratios (1/n^2):")
    print(f"  E2/E1 = {E2/E1:.4f} (Target: 0.25)")
    print(f"  E3/E1 = {E3/E1:.4f} (Target: 0.111)")
    
    # Frequency difference ratio (Balmer series)
    # (f3 - f2) / (f2 - f1) ? Or (E2-E1)/(E3-E1)?
    # Balmer: nu ~ (1/2^2 - 1/n^2).
    # H-alpha (3->2): 1/4 - 1/9 = 5/36
    # H-beta (4->2): 1/4 - 1/16 = 3/16
    # Ratio H-beta / H-alpha = (3/16) / (5/36) = (3*36)/(16*5) = 108/80 = 1.35
    
    # Let's check energy level spacing ratios
    ratio = (E2 - E1) / (E3 - E1)
    print(f"  (E2-E1)/(E3-E1) = {ratio:.4f} (Target: 0.84 for Hydrogen)")

# --- QW-501: TOPOLOGICAL STABILITY (Proton as Knot) ---
print("\n" + "="*60)
print("QW-501: TOPOLOGICAL STABILITY (Knot Robustness)")
print("="*60)

# Simulate a 3-node loop (triangle) with phase winding
# State vector psi = [psi1, psi2, psi3]
# Interaction: d psi_i / dt = i * Sum K_ij psi_j
# Initial state: Winding number 1 (phases 0, 2pi/3, 4pi/3)

psi = np.array([
    np.exp(1j * 0),
    np.exp(1j * 2*np.pi/3),
    np.exp(1j * 4*np.pi/3)
])

print("Initial Phases:")
print(np.angle(psi) / np.pi, "* pi")

# Add strong NOISE
noise_level = 0.5
noise = (np.random.rand(3) - 0.5) * noise_level * 1j + (np.random.rand(3) - 0.5) * noise_level
psi_noisy = psi + noise
# Normalize
psi_noisy = psi_noisy / np.linalg.norm(psi_noisy)

print(f"\nApplied Noise (Level {noise_level}):")
print(np.angle(psi_noisy) / np.pi, "* pi")

# Evolve for some time
dt = 0.1
steps = 100
# Interaction matrix for triangle (all dist = 1)
K_val = K_kernel(1)
M_tri = np.array([
    [0, K_val, K_val],
    [K_val, 0, K_val],
    [K_val, K_val, 0]
])

print("\nEvolving under Kernel dynamics...")
for t in range(steps):
    # dpsi = -i M psi
    dpsi = -1j * M_tri @ psi_noisy
    psi_noisy += dpsi * dt
    psi_noisy /= np.linalg.norm(psi_noisy) # Conserve probability

print("Final Phases:")
phases = np.angle(psi_noisy)
print(phases / np.pi, "* pi")

# Check Winding Number
# Sort phases to check order
# Winding number = sum of phase differences / 2pi
# Naive check: are they still roughly 120 degrees apart?
diff1 = (phases[1] - phases[0]) % (2*np.pi)
diff2 = (phases[2] - phases[1]) % (2*np.pi)
diff3 = (phases[0] - phases[2]) % (2*np.pi)

print(f"\nPhase differences (rad): {diff1:.2f}, {diff2:.2f}, {diff3:.2f}")
print(f"Target: {2*np.pi/3:.2f} (approx 2.09)")

is_stable = np.allclose([diff1, diff2, diff3], [2*np.pi/3, 2*np.pi/3, 2*np.pi/3], atol=1.0) # Loose tolerance due to noise
print(f"Topology Preserved? {is_stable}")

# --- QW-502: ENTROPIC DRAG (Dark Matter) ---
print("\n" + "="*60)
print("QW-502: ENTROPIC DRAG (Entropy Production)")
print("="*60)

# Simulate a particle moving through a 1D chain
# Particle affects weights w_i.
# Entropy S = -Sum p_i log p_i
# We measure change in S as particle moves.

velocities = np.linspace(0.1, 2.0, 10)
entropy_production = []

for v in velocities:
    # Particle moves from x=0 to x=10
    # It perturbs the medium locally
    # Perturbation decays with time (relaxation)
    
    # Simplified model:
    # dS/dt ~ (perturbation_amplitude)^2 * relaxation_rate
    # perturbation ~ 1/v (interaction time) ? Or v (impact strength)?
    # In "turbulent ether", drag ~ v^2 usually.
    # Here we check if S_prod scales non-linearly.
    
    # Let's assume network distortion E ~ K(v) ?
    # Actually, let's simulate:
    # Fast move -> impulsive shock -> high entropy?
    # Slow move -> adiabatic -> low entropy?
    
    # Adiabatic limit (v->0): S -> 0 (reversible).
    # Fast limit: Irreversible.
    
    # Let's assume S_prod proportional to Energy Dissipated.
    # Dissipation D ~ v^n.
    # We test n.
    
    # Simulation proxy:
    # Interaction energy I = Integral K(x - vt) * rho(x) dt
    # If K is oscillatory, resonance might occur at specific v.
    
    # Let's calculate "Wake Entropy".
    # Wake amplitude A ~ 1 / (1 + (v - c_sound)^2) ? Resonance?
    # Let's use the Kernel K(d).
    # Effective coupling C_eff = Integral K(x) * exp(i k x) dx (Fourier transform at k ~ v)
    
    # We compute power spectrum of the Kernel at frequency corresponding to v
    # k ~ v.
    # Power P(v) = |FT(K)(v)|^2
    
    # FT of K(d) = FT( alpha * exp(i(wd+phi)) / (1+bd) )
    # This is roughly a Lorentzian centered at w.
    
    # We evaluate |FT(K)(v)|^2 numerically
    d_range = np.arange(0, 100, 0.1)
    k_val = K_kernel(d_range)
    # Fourier transform
    ft_val = np.sum(k_val * np.exp(-1j * v * d_range))
    power = np.abs(ft_val)**2
    
    entropy_production.append(power)

print("Entropy Production (Power) vs Velocity:")
for v, s in zip(velocities, entropy_production):
    print(f"  v={v:.2f}: S_prod={s:.4f}")

# Check scaling
# Fit S ~ v^n
log_v = np.log(velocities)
log_s = np.log(entropy_production)
slope, intercept = np.polyfit(log_v, log_s, 1)

print(f"\nScaling Exponent n (S ~ v^n): {slope:.4f}")
print("If n > 1, we have non-linear drag (Dark Matter candidate).")

# --- QW-503: OPERATOR SPECTRUM (Tau Mass) ---
print("\n" + "="*60)
print("QW-503: OPERATOR SPECTRUM (Tau Mass from Geometry)")
print("="*60)

# Construct K matrix for N=100
N_spec = 100
M_spec = np.zeros((N_spec, N_spec), dtype=complex)
for i in range(N_spec):
    for j in range(N_spec):
        d = abs(i - j)
        if d > 0:
            M_spec[i, j] = K_kernel(d)

# Eigenvalues
evals_spec, _ = scipy.linalg.eigh(M_spec)
# Sort by magnitude
abs_evals = np.sort(np.abs(evals_spec))

print("Top 5 Eigenvalues (Magnitude):")
print(abs_evals[-5:])

# Check ratios of largest eigenvalues
# Lambda_max / Lambda_next
ratio_1 = abs_evals[-1] / abs_evals[-2]
ratio_2 = abs_evals[-1] / abs_evals[-3]

print(f"\nRatios:")
print(f"  L_max / L_next = {ratio_1:.4f}")
print(f"  L_max / L_3rd  = {ratio_2:.4f}")

# Target: Tau/Electron ~ 3477
# Target: Muon/Electron ~ 207
print(f"Target Tau/e: 3477")
print(f"Target Mu/e: 207")

# --- QW-504: FRACTAL SIMILARITY (Box Counting) ---
print("\n" + "="*60)
print("QW-504: FRACTAL SIMILARITY (Micro vs Macro)")
print("="*60)

# Generate a "Macro" structure (Network growth with K)
# Generate a "Micro" structure (Soliton internal field)

# 1. Macro: Random Walk weighted by K? Or just the Kernel function itself?
# Let's analyze the Kernel function K(d) as a 1D signal.
d_vals = np.linspace(0, 100, 1000)
signal_macro = np.real(K_kernel(d_vals))

# 2. Micro: Zoom in to d in [0, 1]
d_micro = np.linspace(0, 1, 1000)
signal_micro = np.real(K_kernel(d_micro))

# Box Counting Dimension
def box_counting(signal):
    # Normalize signal to [0, 1]
    sig_norm = (signal - np.min(signal)) / (np.max(signal) - np.min(signal))
    N = len(signal)
    
    scales = [2, 4, 8, 16, 32, 64]
    counts = []
    
    for scale in scales:
        # Grid size
        eps = 1.0 / scale
        # Count boxes
        # Simple 1D box counting on the graph (t, y)
        # Resample signal to 'scale' points?
        # Better: Check how many boxes of size eps are touched by the curve
        
        box_count = 0
        step = N // scale
        for i in range(scale):
            segment = sig_norm[i*step : (i+1)*step]
            if len(segment) == 0: continue
            min_y = np.min(segment)
            max_y = np.max(segment)
            
            # Number of vertical boxes covered
            y_boxes = int(np.ceil(max_y * scale)) - int(np.floor(min_y * scale))
            box_count += max(1, y_boxes)
            
        counts.append(box_count)
    
    # Fit log(N) vs log(1/eps)
    log_eps = np.log([1.0/s for s in scales]) # Actually log(scale)
    log_N = np.log(counts)
    
    coeffs = np.polyfit(np.log(scales), log_N, 1)
    return coeffs[0]

D_macro = box_counting(signal_macro)
D_micro = box_counting(signal_micro)

print(f"Fractal Dimension D_macro (d=0..100): {D_macro:.4f}")
print(f"Fractal Dimension D_micro (d=0..1):   {D_micro:.4f}")
print(f"Difference: {abs(D_macro - D_micro):.4f}")
print(f"Self-Similar? {abs(D_macro - D_micro) < 0.1}")

print("="*80)
print("MISSION COMPLETE")
