
import numpy as np
import time
import json
from scipy.ndimage import laplace

print('=' * 80)
print('QW-1139 TO QW-1158: FIZYCZNA SYMULACJA WĘZŁÓW TORUSOWYCH')
print('=' * 80)

# GRID PARAMETERS
L = 32
x = np.linspace(-3, 3, L)
y = np.linspace(-3, 3, L)
z = np.linspace(-3, 3, L)
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# FIELD
# Complex scalar field psi
# A knot is a zero-locus of psi? Or a tube of phase winding?
# We use the standard construction of torus knots from accumulation of phases.
# Torus knot T(p,q) can be defined by:
# u = (X + iY) / sqrt(X^2 + Y^2) * R + ...
#
# Better: Milnor map construction for T(p,q).
# f(u,v) = u^p + v^q where u,v are complex coords in C^2 -> R^4.
# Intersection with S^3 sphere gives the knot.
# We map R^3 -> S^3 using stereographic projection.

def stereographic_projetion(X, Y, Z):
    R2 = X**2 + Y**2 + Z**2
    # Inverse stereo: R^3 -> S^3 \ {North Pole}
    # (x,y,z) -> (X,Y,Z,W) on S^3
    denom = 1 + R2
    x_s = 2*X / denom
    y_s = 2*Y / denom
    z_s = 2*Z / denom
    w_s = (1 - R2) / denom
    return x_s, y_s, z_s, w_s

def initial_knot_field(p, q):
    """Constructs complex scalar field with T(p,q) zero locus."""
    # Map grid to S^3
    x_s, y_s, z_s, w_s = stereographic_projetion(X, Y, Z)
    
    # Coordinates in C^2: (z1, z2) where |z1|^2 + |z2|^2 = 1
    # z1 = x_s + i y_s
    # z2 = z_s + i w_s
    z1 = x_s + 1j * y_s
    z2 = z_s + 1j * w_s
    
    # Milnor polynomial: f = z1^p + z2^q
    psi = z1**p + z2**q
    
    # Normalize to have magnitude 1 far away?
    # Ginzburg-Landau potential will handle magnitude.
    return psi

def energy_density(psi):
    # Kinetic: |grad psi|^2
    # Approx with finite diff
    d_psi_x = np.roll(psi, -1, axis=0) - psi
    d_psi_y = np.roll(psi, -1, axis=1) - psi
    d_psi_z = np.roll(psi, -1, axis=2) - psi
    
    grad_sq = np.abs(d_psi_x)**2 + np.abs(d_psi_y)**2 + np.abs(d_psi_z)**2
    
    # Potential: V = lambda * (1 - |psi|^2)^2
    # Standard Higgs/GL potential
    V = 1.0 * (1 - np.abs(psi)**2)**2
    
    return 0.5 * grad_sq + V

def relax_field(psi, steps=50):
    # Ginzburg-Landau relaxation (Gradient Descent on Energy)
    # dpsi/dt = - dE/dpsi* = Laplacian psi + 2*psi*(1-|psi|^2)
    dt = 0.05
    for i in range(steps):
        # Laplacian
        lap = laplace(psi) # Simple convolution
        # Non-linear term
        nl = 2 * psi * (1 - np.abs(psi)**2)
        
        # Update
        psi += dt * (lap + nl)
        
        # Dirichlet BCs (keep edges fixed? or allow float?)
        # For simplicity, torus BC (np.roll) implies periodic.
        # But we used stereographic proj which dies at infinity (0,0,0,-1).
        # We need to re-enforce boundaries?
        # Let's just evolve.
        
    return psi

def analyze_knot(p, q, name):
    print(f'\nRunning Simulation for {name} T({p},{q})...')
    
    # 1. Initialize
    psi = initial_knot_field(p, q)
    
    # 2. Relax (find ground state for this topology)
    psi = relax_field(psi, steps=50) # Short relaxation
    
    # 3. Calculate Metrics
    E_dens = energy_density(psi)
    Total_Energy = np.sum(E_dens)
    
    # Core Radius / Volume
    # Core is where |psi| < 0.5 (approx)
    core_mask = np.abs(psi) < 0.5
    Core_Volume = np.sum(core_mask)
    
    # Effective Radius ~ Volume^(1/3)
    r_core = Core_Volume**(1.0/3.0) if Core_Volume > 0 else 0
    
    # Predicted Mass from hypothesis
    # M ~ r_core^-1.52
    M_model = r_core**(-1.52) if r_core > 0 else 0
    
    print(f'  Total Energy: {Total_Energy:.1f}')
    print(f'  Core Volume: {Core_Volume}')
    print(f'  Effective r_core: {r_core:.3f}')
    print(f'  Model Mass (r^-1.52): {M_model:.4f}')
    
    return {
        'name': name,
        'p': p, 'q': q,
        'E': Total_Energy,
        'Vol': Core_Volume,
        'r_core': r_core,
        'M_model': M_model
    }

# LIST OF PARTICLES (HYPOTHETICAL KNOTS)
# Crossing numbers:
# Top: 0 (Unknot T(1,1))
# Bottom: 7? (T(2,5) c=5? No. T(3,2) c=3. T(2,7) c=?)
# Let's test a sequence of knots and see the trend.

knots_to_test = [
    (1, 1, 'Top (Unknot)'),      # c=0
    (3, 2, 'Trefoil (3)'),        # c=3
    (5, 2, 'Cinquefoil (5)'),     # c=5 (Solomon seal is link)
    (4, 3, 'T(4,3) (8)'),         # c=8 (Close to Bottom 7 / Tau 9)
    (5, 3, 'T(5,3) (10)'),
    (5, 4, 'T(5,4) (15)'),        
    (6, 5, 'T(6,5) (24)'),        # c=24 (Electron candidate?)
]

results = []
for p, q, name in knots_to_test:
    res = analyze_knot(p, q, name)
    results.append(res)

print('\n' + '='*80)
print('WYNIKI SYMULACJI')
print(f'{"Name":<20} | {"(p,q)":<8} | {"Vol":<8} | {"r_core":<8} | {"Mass (model)":<12}')
print('-' * 80)

for r in results:
    print(f'{r["name"]:<20} | ({r["p"]},{r["q"]})   | {r["Vol"]:<8} | {r["r_core"]:<8.3f} | {r["M_model"]:<12.4f}')

print('\nANALIZA TRENDU:')
# Try to match electron/muon/tau
# Assume T(1,1) is Top.
# Check if T(6,5) has small mass relative to T(1,1).
m_top = results[0]['M_model']
print(f'Ref M_top (Unknot): {m_top:.4f}')

for r in results[1:]:
    ratio = m_top / r['M_model'] if r['M_model'] > 0 else 0
    print(f'{r["name"]} Mass Ratio: 1/{ratio:.1f}')

# Electron actual ratio ~ 340,000
# Does T(6,5) give 1/340,000?
# Unlikely with such small grid, but let's see the trend.

