import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# ==============================================================================
# QW-1547_Corrected: Emergent Geometric Gauge (Phase on Edges)
# ==============================================================================
# Objective: Demonstrate that the "Gauge Field" A_mu is not a fundamental field
# but simply the Phase Difference on Graph Edges.
#
# Rules (Block B):
# 1. Transform Nodes: psi_i -> e^{i theta_i} psi_i.
# 2. Action is NOT invariant unless Link Q_ij -> e^{i(theta_i - theta_j)} Q_ij.
# 3. Identify A_mu ~ arg(Q_ij).

REPORT = "RAPORT_QW1547_GEOMETRIC_GAUGE_CORRECTED.md"
md = []

def log(msg=""):
    print(msg)
    md.append(msg)

log("="*80)
log("QW-1547 CORRECTED: GAUGE FIELD AS GRAPH PHASE")
log("="*80)

# ------------------------------------------------------------------------------
# 1. Setup Graph and Action
# ------------------------------------------------------------------------------
# Link Q_ij. Node psi_i.
# Kinetic Term T = sum |psi_i - Q_ij psi_j|^2.
# This is the "Covariant Derivative on Graph".
# If Q_ij = 1 (Trivial), T ~ |grad psi|^2.
# If Q_ij = exp(i A), T ~ |(grad - iA) psi|^2.

N = 10
psi = np.ones(N, dtype=complex) # Constant field
Q = np.ones(N, dtype=complex)   # Trivial links (1D chain, link i->i+1)

def calc_action(psi, Q):
    # Sum |psi_i - Q_i * psi_{i+1}|^2
    S = 0.0
    for i in range(N-1):
        term = psi[i] - Q[i] * psi[i+1]
        S += np.abs(term)**2
    return S

S0 = calc_action(psi, Q)
log(f"Initial Action (Flat): {S0:.4f}")

# ------------------------------------------------------------------------------
# 2. Perform Local Gauge Transformation (Node Rotation)
# ------------------------------------------------------------------------------
log("\nApplying Random Local Gauge Transformation theta_i...")
theta = np.random.uniform(0, 2*np.pi, N)
psi_transformed = psi * np.exp(1j * theta)

# If Q is unchanged, Action should INCREASE (Symmetry Breaking / "Mass" term).
S_broken = calc_action(psi_transformed, Q)
log(f"Action with Transformed Psi ONLY: {S_broken:.4f} (Not Invariant!)")

# ------------------------------------------------------------------------------
# 3. Transform Connection Q_ij (Geometric Compensation)
# ------------------------------------------------------------------------------
# We must rotate Q_ij to compensate.
# psi_i -> U_i psi_i
# Term: U_i psi_i - Q_ij U_j psi_j
# We want: U_i (psi_i - Q_new psi_j)
# So Q_new = U_i Q_old U_j^dagger.
# Q_i connects i and i+1.
# Q_new[i] = exp(i theta[i]) * Q[i] * exp(-i theta[i+1])

log("\nApplying Compensating Gauge Transformation to Links Q...")
Q_transformed = np.zeros(N, dtype=complex)
for i in range(N-1):
    Q_transformed[i] = np.exp(1j * theta[i]) * Q[i] * np.exp(-1j * theta[i+1])

# Is Action invariant now?
S_restored = calc_action(psi_transformed, Q_transformed)
log(f"Action with Transformed Psi AND Links: {S_restored:.4f}")

# ------------------------------------------------------------------------------
# 4. Results
# ------------------------------------------------------------------------------
# Check invariance
diff = abs(S_restored - S0)
if diff < 1e-9:
    log("\n>> SUCCESS: Gauge Invariance confirmed.")
    log(">> The field A_mu corresponds to the phase of the Link Q_ij.")
    log(">> A_mu transforms exactly as required to preserve Action.")
else:
    log("\n>> FAILED: Action not invariant.")

# Save Report
with open(REPORT, "w", encoding="utf-8") as f:
    f.write(f"# QW-1547 Corrected: Geometric Gauge\n\n")
    f.write(f"**Date:** {datetime.now()}\n\n")
    f.write("## Hypothesis\n")
    f.write("Maxwell field A_mu is EMERGENT. It is simply the phase of the graph link $Q_{ij}$.\n")
    f.write("Local phase rotation of nodes $\\psi_i$ forces $Q_{ij}$ to rotate, creating the gauge transformation law $A \\to A + \\nabla \\theta$.\n\n")
    f.write("## Results\n")
    f.write("```\n")
    f.write("\n".join(md))
    f.write("\n```\n")

print(f"\n✅ Report saved to {REPORT}")
