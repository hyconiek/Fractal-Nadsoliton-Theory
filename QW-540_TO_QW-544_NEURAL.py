# QW-540 TO QW-544: EVOLVING NEURAL UNIVERSE
# PARADIGM: Physics is Learning (H9). Reality is Maximum Resonance (H11).
# GOAL: Verify if Evolving Kernel reproduces Gravity, Fine Tuning, and Time.

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve2d

print("="*80)
print("QW-540 TO QW-544: EVOLVING NEURAL UNIVERSE")
print("Testing the 'Physics as Learning' Hypothesis")
print("="*80)

# Base Kernel Function (for reference)
def K_kernel(d, alpha, beta, omega, phi):
    denom = 1 + beta * d
    denom = np.maximum(denom, 1e-6)
    return (alpha * np.cos(omega * d + phi)) / denom

# --- QW-540: HEBBIAN GRAVITY (Wiring = Warping) ---
print("\n" + "="*60)
print("QW-540: HEBBIAN GRAVITY (Wiring = Warping)")
print("="*60)

# Hypothesis: Active nodes strengthen connections. Stronger K = Shorter Distance.
# Simulation: 1D chain. Two "masses" (clamped active nodes).
# Evolve K_ij based on Hebbian rule: dK_ij/dt = eta * |Psi_i * Psi_j| - decay * K_ij
# Measure effective force.

N_nodes = 20
# Initial K: Uniform / Distance based
coords = np.arange(N_nodes)
K_matrix = np.zeros((N_nodes, N_nodes))
for i in range(N_nodes):
    for j in range(N_nodes):
        d = abs(i - j)
        if d > 0: K_matrix[i, j] = 1.0 / d # Initial weak connectivity

# Place Masses at index 5 and 15
Psi = np.random.rand(N_nodes) * 0.1 # Background noise
Psi[5] = 1.0 # Mass 1
Psi[15] = 1.0 # Mass 2

# Hebbian Evolution
eta = 0.1
decay = 0.01
steps = 50

for t in range(steps):
    # Update K
    # Hebbian term: Correlated activity strengthens link
    # We assume "Activity" propagates.
    # Simple diffusion of activity?
    # Psi_new = K * Psi
    Psi_new = K_matrix @ Psi
    Psi_new = Psi_new / (np.max(np.abs(Psi_new)) + 1e-9) # Normalize
    Psi = 0.5 * Psi + 0.5 * Psi_new # Smooth update
    Psi[5] = 1.0 # Clamp Mass 1
    Psi[15] = 1.0 # Clamp Mass 2
    
    # Hebbian Update
    for i in range(N_nodes):
        for j in range(i+1, N_nodes):
            term = eta * Psi[i] * Psi[j]
            K_matrix[i, j] += term - decay * K_matrix[i, j]
            K_matrix[j, i] = K_matrix[i, j]

# Measure "Force" or "Distance"
# Effective Distance D_eff = 1 / K_eff
# Check K between 5 and 15 vs K between 5 and 10 (empty space)
K_mass_mass = K_matrix[5, 15]
K_mass_void = K_matrix[5, 10]

print(f"Connection Strength Mass-Mass (d=10): {K_mass_mass:.4f}")
print(f"Connection Strength Mass-Void (d=5):  {K_mass_void:.4f}")

# Gravity implies Mass-Mass connection should be anomalously strong compared to background.
if K_mass_mass > K_mass_void: # Wait, d=10 vs d=5. K should be weaker for d=10 usually.
    # But Hebbian learning might make it stronger!
    print("Result: Hebbian Bridge Formed (Action at a Distance).")
else:
    # Normalize by distance
    # K_expected ~ 1/d.
    ratio_mm = K_mass_mass * 10
    ratio_mv = K_mass_void * 5
    print(f"Distance-Normalized Strength (MM vs MV): {ratio_mm:.2f} vs {ratio_mv:.2f}")
    if ratio_mm > ratio_mv:
        print("Result: GRAVITY-LIKE (Masses pull space together).")
    else:
        print("Result: No significant warping.")

# --- QW-541: FINE TUNING (Evolutionary Constants) ---
print("\n" + "="*60)
print("QW-541: FINE TUNING (Evolutionary Constants)")
print("="*60)

# Can Evolution find alpha=2.77, beta=0.01?
# Fitness = Total Resonance Energy of the Kernel itself.
# E = Sum K(d)^2 ? Or Sum K(d) (if K is attractive).
# Let's assume the Universe maximizes "Connectivity" or "Information Flow".
# Maximize Integral |K(d)| subject to some cost (e.g. energy cost of alpha).

# Let's try to evolve alpha and beta to maximize "Resonance with Integer Lattice".
# H11 says reality is max resonance.
# Resonance R = Sum_{d=1..N} K(d) * cos(2*pi*d)? No.
# Resonance means K(d) has peaks at integer d?
# Or K(d) allows stable standing waves.

# Let's use a simple Genetic Algorithm to find params that maximize
# R = Sum_{i,j} K(|i-j|) for a fixed lattice.
# Constraint: K(0) is fixed. Energy is finite.

def fitness(alpha, beta):
    # Calculate Kernel sum on a grid
    d_vals = np.arange(1, 20)
    K_vals = K_kernel(d_vals, alpha, beta, np.pi/4, np.pi/6)
    # We want high connectivity (Sum K) but low cost (Sum K^2? or beta penalty?)
    # If beta is too small, K diverges -> infinite energy cost.
    # So Fitness = Sum(K) - Lambda * Sum(K^2)
    S1 = np.sum(np.abs(K_vals))
    S2 = np.sum(K_vals**2)
    return S1 - 0.5 * S2 # Maximize signal/noise ratio?

# Evolution
best_alpha = 0
best_beta = 0
best_fit = -1e9

for i in range(1000):
    a = np.random.rand() * 5.0
    b = np.random.rand() * 0.1
    f = fitness(a, b)
    if f > best_fit:
        best_fit = f
        best_alpha = a
        best_beta = b

print(f"Evolved Parameters:")
print(f"  Alpha: {best_alpha:.4f} (Target: ~2.77)")
print(f"  Beta:  {best_beta:.4f}  (Target: ~0.01)")

if 2.0 < best_alpha < 3.5 and best_beta < 0.05:
    print("Result: SUCCESS. Evolution finds FIN-like parameters.")
else:
    print("Result: CONVERGED to different optimum.")

# --- QW-542: ARROW OF TIME (Entropy of Learning) ---
print("\n" + "="*60)
print("QW-542: ARROW OF TIME (Entropy of Learning)")
print("="*60)

# Measure Entropy of K_matrix during Hebbian Learning (from QW-540)
# Shannon Entropy of connection weights.
# H = - Sum p log p, where p = K_ij / Sum K

# Re-run simplified Hebbian
K_mat = np.random.rand(10, 10) * 0.1
Psi = np.random.rand(10)
Entropies = []

for t in range(50):
    # Hebbian
    for i in range(10):
        for j in range(10):
            K_mat[i, j] += 0.1 * Psi[i] * Psi[j]
    # Normalize K to probability distribution
    P = np.abs(K_mat) / np.sum(np.abs(K_mat))
    H = -np.sum(P * np.log(P + 1e-9))
    Entropies.append(H)
    # Update Psi (Dynamics)
    Psi = np.tanh(K_mat @ Psi)

print(f"Initial Entropy: {Entropies[0]:.4f}")
print(f"Final Entropy:   {Entropies[-1]:.4f}")

if Entropies[-1] < Entropies[0]:
    print("Result: Entropy DECREASE (Ordering/Learning).")
    print("Interpretation: Time's Arrow is the formation of Structure.")
else:
    print("Result: Entropy INCREASE (Thermalization).")

# --- QW-543: DARK ENERGY (Neural Forgetting) ---
print("\n" + "="*60)
print("QW-543: DARK ENERGY (Neural Forgetting)")
print("="*60)

# Simulate "Forgetting" (Decay of K)
# If K decreases, Effective Distance D ~ 1/K increases.
# This looks like Expansion of Universe.

K_link = 1.0
decay_rate = 0.01
Distances = []

for t in range(20):
    K_link *= (1.0 - decay_rate)
    D_eff = 1.0 / (K_link + 1e-9)
    Distances.append(D_eff)

print(f"Distance t=0:  {Distances[0]:.2f}")
print(f"Distance t=20: {Distances[-1]:.2f}")

# Calculate Hubble Parameter H = (dD/dt) / D
H_param = (Distances[-1] - Distances[-2]) / Distances[-1]
print(f"Hubble Parameter H: {H_param:.4f}")

if H_param > 0:
    print("Result: EXPANSION CONFIRMED (Forgetting = Dark Energy).")
else:
    print("Result: Collapse.")

# --- QW-544: PARTICLE MEMORY (Stable Attractors) ---
print("\n" + "="*60)
print("QW-544: PARTICLE MEMORY (Stable Attractors)")
print("="*60)

# Can the network "remember" a particle?
# Train network on a pattern, then remove input.
# Does the pattern persist? (Hopfield Network behavior)

N_net = 50
Patterns = np.sign(np.random.randn(N_net)) # A "Particle" pattern
Weights = np.zeros((N_net, N_net))

# Train (Imprint)
for i in range(N_net):
    for j in range(N_net):
        if i != j:
            Weights[i, j] += Patterns[i] * Patterns[j]

# Test Stability
State = np.copy(Patterns)
# Corrupt state
State[:10] = -State[:10] # Flip 20% spins

print(f"Initial Overlap (Corrupted): {np.mean(State * Patterns):.2f}")

# Relax
for _ in range(5):
    State = np.sign(Weights @ State)

Final_Overlap = np.mean(State * Patterns)
print(f"Final Overlap (Recovered):   {Final_Overlap:.2f}")

if Final_Overlap > 0.9:
    print("Result: STABLE MEMORY. (Particle is an Attractor).")
else:
    print("Result: AMNESIA. (Particle decays).")

print("="*80)
print("MISSION COMPLETE")
