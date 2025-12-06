#!/usr/bin/env python3
"""
QW_Resonance_Modes.py
Purpose: Test the hypothesis that Lepton Generations (e, mu, tau) are 
         Resonance Modes (Eigenvalues) of the Information Kernel on the SAME layer.
         Hypothesis: Mass ~ Eigenvalue * "Computational Cost" (Entropy)
"""

import numpy as np
from scipy.linalg import eigh
from scipy.stats import entropy
import matplotlib.pyplot as plt

# --- Fundamental Constants (No fitting!) ---
ALPHA = 4 * np.log(2)   # Information capacity (4 bits)
OMEGA = np.pi / 4       # Octave phase
PHI = np.pi / 6         # Golden Angle
BETA = 0.01             # Torsion damping

# --- Experimental Data ---
M_E = 0.511
M_MU = 105.66
M_TAU = 1776.86

def analyze_resonance_modes(grid_size=24): # 24 for 2*12 octaves? Or just a spatial patch?
    # Let's assume the "particle" is a localized patch of the network.
    # Grid size determines the number of available modes.
    
    print(f"Analyzing Kernel Matrix ({grid_size}x{grid_size})...")
    
    # K(d) = alpha * cos(omega*d + phi) / (1 + beta*d)
    K = np.zeros((grid_size, grid_size))
    for i in range(grid_size):
        for j in range(grid_size):
            if i == j:
                K[i, j] = 1.0 # Self-interaction? Or 0? Usually diagonal is distinct.
                              # Let's say self-interaction is max intensity = ALPHA
                K[i, j] = ALPHA
            else:
                d = abs(i - j)
                K[i, j] = ALPHA * np.cos(OMEGA * d + PHI) / (1 + BETA * d)
                
    # Diagonalize
    evals, evecs = eigh(K)
    
    # Sort eigenvalues by magnitude (strength of resonance)
    # We are interested in STABLE modes, so usually large positive eigenvalues?
    # Or maybe "Zero Modes" (near zero)?
    # Let's look at absolute magnitude for now.
    idx = np.argsort(np.abs(evals))
    evals = evals[idx]
    evecs = evecs[:, idx]
    
    # We need 3 distinct stable modes corresponding to generations.
    # Hypothesis: The 3 lowest frequency modes (smoothest eigenvectors)?
    # Or the 3 highest energy modes?
    
    # Let's classify modes by "frequency" (number of sign changes / nodes in eigenvector)
    modes = []
    for i in range(grid_size):
        vec = evecs[:, i]
        # Count zero crossings
        zero_crossings = np.sum(np.diff(np.signbit(vec)))
        
        # Calculate Shannon Entropy of the squared wavefunction (probability distribution)
        prob = vec**2
        prob /= np.sum(prob)
        S = entropy(prob, base=2)
        
        modes.append({
            'id': i,
            'eigenvalue': evals[i],
            'abs_val': abs(evals[i]),
            'nodes': zero_crossings,
            'entropy': S
        })
        
    return modes

def test_mass_hypothesis(modes):
    print("\n--- Testing Mass Hypothesis ---")
    print("Hypothesis: Mass(n) ~ |Eigenvalue_n| * Cost_n")
    print("Where Cost_n ~ exp(Entropy_n)? Or Entropy_n?")
    
    # Let's filter for "Significant" modes.
    # Maybe we just take the top 3 strongest modes?
    
    sorted_by_strength = sorted(modes, key=lambda x: x['abs_val'], reverse=True)
    
    print("\nTop 5 Strongest Modes:")
    print(f"{'Rank':<5} {'EigVal':<10} {'Nodes':<5} {'Entropy':<10}")
    for i in range(5):
        m = sorted_by_strength[i]
        print(f"{i:<5} {m['abs_val']:.4f}     {m['nodes']:<5} {m['entropy']:.4f}")
        
    # Let's try to identify candidates for e, mu, tau
    # Assumption: Mass increases with Generation.
    # Does Eigenvalue decrease or increase?
    # Does Entropy increase?
    
    # Let's verify ratio mu/e ~ 206
    # and tau/mu ~ 16.8
    
    candidates = sorted_by_strength # Just take all
    
    # Brute force search for triplet with correct ratios?
    # Or look for specific "Node" counts? (Fundamental=0 nodes, 1st Excited=1 node...)
    
    # Sort by Nodes (Geometric oscillations)
    sorted_by_nodes = sorted(modes, key=lambda x: x['nodes'])
    
    print("\nModes sorted by Oscillation Nodes (Topology?):")
    print(f"{'Nodes':<5} {'EigVal':<10} {'Entropy':<10} {'Product(M)':<15}")
    
    potential_masses = []
    
    for m in sorted_by_nodes:
        # Proposed Mass Formula: M ~ |lambda| * 2^Entropy
        # Why 2^Entropy? Information cost is exponential in bits?
        # Or M ~ |lambda| / Entropy?
        
        # Let's try: M ~ |lambda| * exp(Entropy)
        mass_proxy = m['abs_val'] * np.exp(m['entropy'])
        
        # Or Inverse? Lower eigenvalue = more stable = less mass?
        # WAIT. Mass is ENERGY. Higher eigenvalue = higher binding energy? 
        # Usually E ~ frequency.
        
        potential_masses.append(mass_proxy)
        
        print(f"{m['nodes']:<5} {m['abs_val']:.4f}     {m['entropy']:.4f}     {mass_proxy:.4f}")

    # Let's check ratios for the first 3 modes (Nodes 0, 1, 2?)
    # usually Ground state has 0 nodes.
    
    m0 = potential_masses[0] # Nodes 0
    m1 = potential_masses[1] # Nodes 1
    m2 = potential_masses[2] # Nodes 2
    
    print("\n--- Ratio Check (Nodes 0,1,2) ---")
    print(f"Gen 1 (0 nodes): {m0:.4f}")
    print(f"Gen 2 (1 nodes): {m1:.4f}")
    print(f"Gen 3 (2 nodes): {m2:.4f}")
    
    if m0 > 0:
        ratio1 = m1/m0
        ratio2 = m2/m1
        print(f"Ratio 2/1 (Mu/E): {ratio1:.4f} (Exp: 206.77)")
        print(f"Ratio 3/2 (Tau/Mu): {ratio2:.4f} (Exp: 16.82)")
        
    # Alternative: Maybe generations are not 0,1,2 but specific stable points?
    pass

modes = analyze_resonance_modes(grid_size=24)
test_mass_hypothesis(modes)
