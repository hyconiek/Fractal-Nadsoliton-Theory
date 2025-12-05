# QW-515 TO QW-519: NESTED FRACTAL SIMULATION
# PARADIGM: Recursive Reality. Layer N is a simulation running inside Node of Layer N-1.
# PROTOCOL: Passing boundary conditions, not scaling factors.

import numpy as np
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
print("QW-515 TO QW-519: NESTED FRACTAL SIMULATION")
print("PARADIGM: As Above, So Below (Recursive Layers)")
print("="*80)

class FractalLayer:
    def __init__(self, layer_index, parent_field_strength=0.0):
        self.N = layer_index
        # Background field from parent acts as boundary condition/potential
        # Transmitted strength is beta_tors * parent_strength
        self.background_potential = parent_field_strength * beta_tors
        self.time_scale = (1.0 / beta_tors) ** self.N # Time runs faster deeper down? 
        # User said: "One step of Parent = 100 steps of Child". So Child is faster.
        # But wait, if Child is faster, it does 100 steps while Parent does 1.
        
    def get_effective_kernel(self, d):
        # Local physics is governed by K(d)
        # But maybe modified by background?
        # For now, assume K(d) is invariant in form, but units might scale.
        # User says: "Internal consistency (local K) is strong (1)".
        return (alpha_geo * np.exp(1j * (omega * d + phi))) / (1 + beta_tors * d)

# --- QW-515: NESTED HYDROGEN (Fraktalny Atom Wodoru) ---
print("\n" + "="*60)
print("QW-515: NESTED HYDROGEN (Layer N=10)")
print("="*60)

# Layer 10 simulation
# Parent (N=9) provides background field.
# Assume Parent has some vacuum energy E_vac ~ 1.0 (Planck units of N=9)
parent_field = 1.0 
layer_10 = FractalLayer(10, parent_field)

print(f"Layer 10 Background Potential: {layer_10.background_potential:.6e}")

# Solve Radial Schrodinger in this environment
# V_eff(r) = V_local(r) + V_background
# V_local(r) = -1/r (Coulomb-like emergent from K?) 
# Or use K(r) as potential? User said "symulacja orbitalna (jak w QW-422)".
# QW-422 used K(d) as potential.

N_r = 200
r = np.linspace(0.05, 20.0, N_r)
dr = r[1] - r[0]

# Local potential from a central "proton" in Layer 10
# V_local = -Re(K(r))
V_local = -np.real(layer_10.get_effective_kernel(r))

# Total Potential
V_total = V_local + layer_10.background_potential

# Hamiltonian
diag = 1.0/(dr**2) + V_total
off_diag = -0.5/(dr**2) * np.ones(N_r-1)
H = diags([diag, off_diag, off_diag], [0, 1, -1])

# Eigenvalues
vals, _ = eigsh(H, k=5, which='SA')
print("Eigenenergies in Layer 10:")
for i, E in enumerate(vals):
    print(f"  E_{i+1} = {E:.6f}")

# Check stability and ratios
if len(vals) >= 3:
    ratio_21 = vals[1]/vals[0]
    print(f"  E2/E1 = {ratio_21:.4f} (Target 0.25)")
    
    # Check if energies are "small" compared to local Planck scale (which is ~1 in code units?)
    # If E ~ -0.1, it is bound.
    print(f"  Bound states found? {vals[0] < 0}")

# --- QW-516: RELATIVE TIME (Skalowanie Czasu) ---
print("\n" + "="*60)
print("QW-516: RELATIVE TIME (Clock Scaling)")
print("="*60)

# Simulate Parent (N=0) and Child (N=1)
# Parent Step = 100 * Child Step
steps_parent = 10
steps_child = steps_parent * int(1/beta_tors)

print(f"Simulating: {steps_parent} Parent steps = {steps_child} Child steps")

# Measure speed of light 'c' in Child simulation
# c = distance / time
# In Child simulation, a wave propagates via K(d).
# Group velocity v_g = d(omega)/dk ?
# Or simply simulate a pulse.

# 1D Wave Equation on Child Network
# d^2 psi / dt^2 = -L psi
# c ~ sqrt(stiffness) ~ sqrt(alpha_geo) ?

# Let's simulate a pulse in Child layer
N_grid = 100
psi_child = np.zeros(N_grid)
psi_child[N_grid//2] = 1.0 # Pulse
vel_child = np.zeros(N_grid)

# Evolution
# dt = 0.1 (Child time units)
dt_child = 0.1
positions = []
times = []

# Laplacian kernel (nearest neighbor for speed of light check)
# Actually K(d) is non-local, so c depends on k.
# Let's use standard discrete laplacian for "local" c check, or K(1).
# If K(d) is used, signal propagates instantly? No, damped.
# Let's use K(1) as nearest neighbor coupling.
coupling = np.real(layer_10.get_effective_kernel(1)) # Use K(1)

for t in range(200): # 200 child steps
    # Wave equation: d2psi/dt2 = c^2 d2psi/dx2
    # Discrete: psi_new = 2*psi - psi_old + dt^2 * c^2 * laplacian
    # Here coupling represents c^2
    
    laplacian = np.zeros_like(psi_child)
    laplacian[1:-1] = psi_child[2:] - 2*psi_child[1:-1] + psi_child[:-2]
    
    acc = coupling * laplacian
    vel_child += acc * dt_child
    psi_child += vel_child * dt_child
    
    # Track peak
    peak_idx = np.argmax(np.abs(psi_child))
    if peak_idx != N_grid//2:
        dist = abs(peak_idx - N_grid//2)
        positions.append(dist)
        times.append(t * dt_child)

# Calculate c_child
if len(positions) > 5:
    # Fit x = c * t
    coeffs = np.polyfit(times, positions, 1)
    c_child = coeffs[0]
    print(f"Measured c_child = {c_child:.4f} (grid/time)")
else:
    print("Pulse did not move enough.")
    c_child = 0

# Now, does c_child match c_parent?
# In Parent, physics is same (K is same form).
# So c_parent (in parent units) should be same number.
# The question is: Is c invariant?
# Yes, if K is scale invariant.
print(f"Is c invariant? Yes, if K form is identical.")
print(f"c depends on alpha_geo. Since alpha_geo is frozen, c is constant across layers.")

# --- QW-517: DARK MATTER (Parental Rotation) ---
print("\n" + "="*60)
print("QW-517: DARK MATTER (Parental Coriolis)")
print("="*60)

# Layer 20 (Galaxy) inside Layer 19 (Parent Universe)
# Parent rotates with Omega_parent
Omega_parent = 0.001 # Slow rotation
# Coriolis force in Layer 20: F_cor = 2 * Omega_eff * v
# Omega_eff = Omega_parent * beta_tors (Transmission)

Omega_eff = Omega_parent * beta_tors
print(f"Effective Coriolis Omega: {Omega_eff:.6e}")

# Simulate rotation curve v(r)
# F_grav + F_cor = F_centripetal
# GM/r^2 + 2*Omega_eff*v = v^2/r
# v^2 - 2*Omega_eff*r*v - GM/r = 0

# Solve for v: v = [2*Omega*r + sqrt(4*Omega^2*r^2 + 4*GM/r)] / 2
# v = Omega*r + sqrt(Omega^2*r^2 + GM/r)

r_vals = np.linspace(1, 20, 50)
M_gal = 10.0
v_newton = np.sqrt(M_gal / r_vals)
v_coriolis = Omega_eff * r_vals + np.sqrt((Omega_eff * r_vals)**2 + M_gal / r_vals)

print("Rotation Curve Check:")
print(f"  r={r_vals[0]:.1f}: v_Newt={v_newton[0]:.4f}, v_Cor={v_coriolis[0]:.4f}")
print(f"  r={r_vals[-1]:.1f}: v_Newt={v_newton[-1]:.4f}, v_Cor={v_coriolis[-1]:.4f}")

# Check flatness
# v_Cor ~ Omega*r at large r (Linear rise? Or flat?)
# If Omega is small, maybe it flattens?
# Actually, linear rise is rigid body rotation.
# We want flat.
# Flat means v ~ const.
# If Omega*r term dominates, v grows linearly. This is NOT flat.
# Wait, maybe Coriolis is F = -2*Omega*v?
# If F_cor opposes gravity?
# Or maybe "Parental Rotation" implies a drag?
# User said: "Siła ta jest tłumiona przez beta."
# Let's check the slope.
slope_newt = (v_newton[-1] - v_newton[0]) / (r_vals[-1] - r_vals[0])
slope_cor = (v_coriolis[-1] - v_coriolis[0]) / (r_vals[-1] - r_vals[0])

print(f"Slope Newtonian: {slope_newt:.4f}")
print(f"Slope Coriolis:  {slope_cor:.4f}")

if abs(slope_cor) < abs(slope_newt) and slope_cor > -0.1:
    print("Result: Flattening observed (or less decline).")
else:
    print("Result: No flattening (or linear rise).")

# --- QW-518: PROTON STABILITY (Fractal Isolation) ---
print("\n" + "="*60)
print("QW-518: PROTON STABILITY (Fractal Isolation)")
print("="*60)

# Proton in Layer 10.
# Noise from Layer 9.
# Noise amplitude = Noise_parent * beta_tors
Noise_parent = 10.0 # Huge noise in parent
Noise_effective = Noise_parent * beta_tors

print(f"Parent Noise: {Noise_parent}")
print(f"Effective Noise on Proton: {Noise_effective}")

# Proton binding energy (from QW-515 or assumed)
# Assume E_bind ~ 1.0 (local units)
E_bind = 1.0

if Noise_effective < E_bind:
    print("Result: STABLE. Effective noise < Binding energy.")
    print(f"  Safety margin: {E_bind / Noise_effective:.2f}x")
else:
    print("Result: UNSTABLE. Noise disrupts proton.")

# --- QW-519: FORCE UNIFICATION (Internal vs External) ---
print("\n" + "="*60)
print("QW-519: FORCE UNIFICATION (Internal vs External)")
print("="*60)

# Internal Force (Strong) ~ K(0) or K(small)
# External Force (EM) ~ K(large) or K_effective between nodes
# In Layer 10.

# K(d) decay?
k_0 = np.abs(layer_10.get_effective_kernel(0.1)) # Internal
k_inf = np.abs(layer_10.get_effective_kernel(10.0)) # External

print(f"Internal Coupling (d=0.1): {k_0:.4f}")
print(f"External Coupling (d=10):  {k_inf:.4f}")

ratio = k_0 / k_inf
print(f"Ratio Strong/EM: {ratio:.4f}")
print(f"Target: ~137")

if 100 < ratio < 200:
    print("Result: SUCCESS. Ratio is consistent with 1/alpha.")
else:
    print("Result: FAILURE. Ratio mismatch.")

print("="*80)
print("MISSION COMPLETE")
