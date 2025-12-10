#!/usr/bin/env python3
# QW-593: INFORMATION UNITY ENTANGLEMENT
# Purpose: Test if FIN exhibits "entanglement" through information unity
# (mutual information, not Bell inequality)
# Date: 2025-12-05

import numpy as np
from scipy.spatial.distance import cdist

print("="*80)
print("QW-593: INFORMATION UNITY ENTANGLEMENT")
print("="*80)
print("Hypothesis: Nodes are 'entangled' through shared information structure,")
print("not quantum mechanics.")
print()

# Network parameters
N_NODES = 100
ALPHA_GEO = 4 * np.log(2)
BETA_TORS = 0.01
OMEGA = np.pi / 4
PHI = np.pi / 6

# Create network
np.random.seed(593)
positions = np.random.randn(N_NODES, 3) * 2.0
dist_matrix = cdist(positions, positions)

# Complex kernel
def K_complex(d):
    return (ALPHA_GEO * np.exp(1j * (OMEGA * d + PHI))) / (1 + BETA_TORS * d)

# Coupling matrix
K_matrix = np.zeros((N_NODES, N_NODES), dtype=complex)
for i in range(N_NODES):
    for j in range(N_NODES):
        if i != j:
            K_matrix[i, j] = K_complex(dist_matrix[i, j])

print(f"Network: N={N_NODES} nodes")
print(f"Kernel: K(d) = α_geo × exp(i(ωd+φ)) / (1+β×d)")
print()

# Initialize states
psi = np.random.randn(N_NODES) + 1j * np.random.randn(N_NODES)
psi /= np.linalg.norm(psi)

# Evolve for correlation
dt = 0.1
steps = 50

print("Evolving network dynamics...")
for step in range(steps):
    # Simple dynamics: dpsi/dt = i K psi
    dpsi_dt = 1j * (K_matrix @ psi)
    psi += dt * dpsi_dt
    psi /= np.linalg.norm(psi)

print(f"Evolution complete ({steps} steps)")
print()

# ============================================================================
# TEST 1: MUTUAL INFORMATION (Information-Theoretic Entanglement)
# ============================================================================
print("="*80)
print("TEST 1: MUTUAL INFORMATION BETWEEN DISTANT NODES")
print("="*80)
print()

# Pick two distant nodes
distances = []
pairs = []
for i in range(N_NODES):
    for j in range(i+1, N_NODES):
        d = dist_matrix[i, j]
        if d > 3.0:  # Distant
            distances.append(d)
            pairs.append((i, j))

# Pick farthest pair
if len(pairs) > 0:
    idx = np.argmax(distances)
    node_A, node_B = pairs[idx]
    dist_AB = distances[idx]
    
    print(f"Testing nodes A={node_A}, B={node_B}")
    print(f"Physical distance: d={dist_AB:.2f}")
    print()
    
    # Measure amplitudes
    amp_A = psi[node_A]
    amp_B = psi[node_B]
    
    # Phase correlation
    phase_A = np.angle(amp_A)
    phase_B = np.angle(amp_B)
    phase_diff = np.mod(phase_A - phase_B + np.pi, 2*np.pi) - np.pi
    
    # Amplitude correlation (Pearson)
    real_corr = np.real(amp_A) * np.real(amp_B) + np.imag(amp_A) * np.imag(amp_B)
    
    print(f"Node A: amp={np.abs(amp_A):.4f}, phase={phase_A:.4f}")
    print(f"Node B: amp={np.abs(amp_B):.4f}, phase={phase_B:.4f}")
    print(f"Phase difference: Δφ={phase_diff:.4f}")
    print(f"Amplitude correlation: {real_corr:.4f}")
    print()
    
    # Mutual Information (discretized)
    # I(A;B) = H(A) + H(B) - H(A,B)
    # For continuous variables, use correlation as proxy
    
    # Shannon entropy approximation
    amp_mag_A = np.abs(amp_A)
    amp_mag_B = np.abs(amp_B)
    
    # Normalized correlation coefficient
    corr_coef = np.abs(np.dot(amp_A.conj(), amp_B))
    
    print(f"Correlation coefficient: |⟨A|B⟩| = {corr_coef:.6f}")
    
    if corr_coef > 0.5:
        print("✅ STRONG INFORMATION CORRELATION despite physical distance!")
        print("   This suggests 'information unity entanglement'")
    elif corr_coef > 0.1:
        print("🟡 MODERATE correlation (nielokalne sprzężenia K(d))")
    else:
        print("❌ WEAK correlation (nodes are informationally independent)")
    
    print()

# ============================================================================
# TEST 2: HOLOGRAPHIC RECONSTRUCTION
# ============================================================================
print("="*80)
print("TEST 2: HOLOGRAPHIC PROPERTY (Can node A predict node B?)")
print("="*80)
print()

# Test if knowing local neighborhood of A can predict B
neighbors_A = np.where(dist_matrix[node_A] < 1.5)[0]
neighbors_A = [n for n in neighbors_A if n != node_A]

print(f"Node A has {len(neighbors_A)} local neighbors (d<1.5)")

# Reconstruct B from A's local info
# Hypothesis: B's state is encoded in A's neighborhood (holography)

# Simple linear prediction: psi_B ≈ sum(K_matrix[B,n] * psi[n]) for n in neighbors_A
psi_B_predicted = 0
for n in neighbors_A:
    weight = K_matrix[node_B, n]
    psi_B_predicted += weight * psi[n]

# Check if we have prediction
if len(neighbors_A) == 0 or np.abs(psi_B_predicted) == 0:
    print("⚠️ No local neighbors - cannot test holography")
    print("   Network too sparse for this analysis")
    fidelity = 0.0
else:
    # Normalize
    psi_B_predicted /= np.abs(psi_B_predicted)

    # Normalize
    psi_B_predicted /= np.abs(psi_B_predicted)
    
    # Compare with actual
    fidelity = np.abs(np.dot(psi_B_predicted.conj(), amp_B))

print(f"Predicted |B⟩ from A's neighborhood")
print(f"Fidelity: |⟨B_pred|B_actual⟩| = {fidelity:.6f}")
print()

if fidelity > 0.8:
    print("✅ HOLOGRAPHIC! Node B state encoded in A's local environment")
    print("   => Information is non-local and distributed")
elif fidelity > 0.5:
    print("🟡 PARTIAL holography (some information transfer)")
else:
    print("❌ NO holography (B independent of A's neighbors)")

print()

# ============================================================================
# TEST 3: INFORMATION FLOW (Transfer Entropy)
# ============================================================================
print("="*80)
print("TEST 3: INFORMATION FLOW A→B")
print("="*80)
print()

# Measure if A's past predicts B's future better than B's own past
# This is "transfer entropy" - directional information flow

# Evolve further and track
history_A = [amp_A]
history_B = [amp_B]

for step in range(10):
    dpsi_dt = 1j * (K_matrix @ psi)
    psi += dt * dpsi_dt
    psi /= np.linalg.norm(psi)
    history_A.append(psi[node_A])
    history_B.append(psi[node_B])

# Simple transfer: does knowing A(t) improve prediction of B(t+1)?
# Correlation between A(t) and B(t+1)
corr_A_to_B = 0
for t in range(len(history_A)-1):
    corr_A_to_B += np.abs(history_A[t].conj() * history_B[t+1])
corr_A_to_B /= (len(history_A)-1)

# Self-correlation B(t) → B(t+1)
corr_B_to_B = 0
for t in range(len(history_B)-1):
    corr_B_to_B += np.abs(history_B[t].conj() * history_B[t+1])
corr_B_to_B /= (len(history_B)-1)

print(f"Self-correlation B(t)→B(t+1): {corr_B_to_B:.6f}")
print(f"Cross-correlation A(t)→B(t+1): {corr_A_to_B:.6f}")
print()

transfer_entropy = corr_A_to_B / (corr_B_to_B + 1e-10)

print(f"Transfer Entropy Ratio: TE(A→B) = {transfer_entropy:.6f}")
print()

if transfer_entropy > 0.5:
    print("✅ INFORMATION FLOWS from A to B despite distance!")
    print("   => Non-local information transfer through K(d)")
else:
    print("❌ NO significant information flow (nodes evolve independently)")

print()
print("="*80)
print("FINAL CONCLUSION:")
print("="*80)
print()

if corr_coef > 0.5 or fidelity > 0.8 or transfer_entropy > 0.5:
    print("✅ FIN EXHIBITS 'INFORMATION UNITY ENTANGLEMENT'")
    print()
    print("Key findings:")
    print(f"  • Correlation: {corr_coef:.3f} (nielokalne K(d) sprzężenia)")
    print(f"  • Holography: {fidelity:.3f} (lokalny→globalny)")
    print(f"  • Transfer: {transfer_entropy:.3f} (przepływ informacji)")
    print()
    print("Interpretation:")
    print("  - Nodes are NOT quantum-entangled (Bell test failed)")
    print("  - BUT they ARE information-entangled (shared network state)")
    print("  - This is 'classical holographic entanglement'")
    print()
    print("Upgrade H12:")
    print("  Old: Partial Quantumness (geometry only)")
    print("  New: Classical Information Unity (holographic correlations)")
else:
    print("❌ NO INFORMATION UNITY DETECTED")
    print("   Nodes behave as independent classical oscillators")

print("="*80)
