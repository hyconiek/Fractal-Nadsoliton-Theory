#!/usr/bin/env python3
"""
QW-691: INTER-LAYER HORIZON (GOLDILOCKS ZONE)
==============================================
Purpose: Determine the "Thickness of the Present" (Coherence Depth).
         How many vertical layers deep can an observer maintain quantum coherence?

Hypothesis:
  There is a finite horizon N_max such that for L > N_max, S(0, L) < 2.
  This N_max defines the "Goldilocks Zone" where quantum effects are visible.
  Beyond this horizon, the fractal layers appear classical (hidden variables).

Method:
  1. Use the same Layer Chain model as QW-690 (Ising Critical Chain).
  2. Measure S(0, L) for L = 1 to 20.
  3. Find the exact crossing point where S drops below 2.0.
  
Output: RAPORT_QW691_INTERLAYER_HORIZON.md
"""

import numpy as np
from scipy.linalg import eigh
import datetime

print("="*80)
print("QW-691: INTER-LAYER HORIZON TEST")
print("="*80)
print("Testing: What is the 'Coherence Depth' of the fractal?")
print()

# --- Parameters ---
N_LAYERS = 20       # Deeper chain to find the horizon
ALPHA_LAYER = np.log(2) 

# Pauli Matrices
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
sigma_0 = np.eye(2, dtype=complex)

REPORT_FILE = "RAPORT_QW691_INTERLAYER_HORIZON.md"

def build_op(op, site, n_sites):
    # For large N, we cannot build full matrices (2^20 is too big).
    # We must use DMRG or just limit N to ~12-14 and extrapolate.
    # Given limited compute, let's do N=12 and see if we cross the horizon.
    # If horizon > 12, we extrapolate.
    if n_sites > 12:
        raise ValueError("N too large for Exact Diag")
    
    result = 1
    for i in range(n_sites):
        if i == site:
            result = np.kron(result, op)
        else:
            result = np.kron(result, sigma_0)
    return result

# Limit N for simulation
N_SIM = 12
print(f"1. Building Chain (N={N_SIM})...")
dim_total = 2**N_SIM

H = np.zeros((dim_total, dim_total), dtype=complex)

# Need to redefine build_op for N_SIM to avoid closure issues if copied
Sx = [build_op(sigma_x, i, N_SIM) for i in range(N_SIM)]
Sz = [build_op(sigma_z, i, N_SIM) for i in range(N_SIM)]

# Critial Ising Model (H = -Sum ZZ - Sum X)
for i in range(N_SIM - 1):
    H += -1.0 * Sz[i] @ Sz[i+1]
    H += -1.0 * Sx[i]
H += -1.0 * Sx[N_SIM-1]

print("   Finding Ground State... (this may take 1-2 mins)")
# Use eigh for Hermitian
evals, evecs = eigh(H)
psi_ground = evecs[:, 0]

# --- MEASUREMENT ---
print("\n2. Measuring Bell Horizon...")

theta_a, theta_ap = 0, np.pi/2
theta_b, theta_bp = np.pi/4, 3*np.pi/4

def measure_S(psi, site_A, site_B):
    def corr(thA, thB):
        opA = np.sin(thA)*Sx[site_A] + np.cos(thA)*Sz[site_A]
        opB = np.sin(thB)*Sx[site_B] + np.cos(thB)*Sz[site_B]
        obs = opA @ opB
        return np.real(psi.conj().T @ obs @ psi)
    
    E1 = corr(theta_a, theta_b)
    E2 = corr(theta_a, theta_bp)
    E3 = corr(theta_ap, theta_b)
    E4 = corr(theta_ap, theta_bp)
    return np.abs(E1 - E2 + E3 + E4)

results = []
horizon_L = None

print(f"   Observer at L=0. Scanning target L...")
for layer_L in range(1, N_SIM):
    S_val = measure_S(psi_ground, 0, layer_L)
    print(f"   L={layer_L}: S={S_val:.4f}")
    results.append((layer_L, S_val))
    
    if S_val < 2.0 and horizon_L is None:
        if len(results) >= 2:
            horizon_L = layer_L - 1 + (S_val - 2.0)/(results[-2][1] - S_val) 
        else:
            # Crossed at first step! Linear interp from S=2.82 (Theory max) at L=0?
            # Assume S(0) = 2.82
            horizon_L = (2.0 - 2.82) / (S_val - 2.82) # crude
            horizon_L = max(0.1, horizon_L) # Avoid 0 if S_val is weird
        
        print(f"   >>> HORIZON CROSSED at L={layer_L} (S < 2.0) <<<")

# --- ANALYSIS ---
print("\n3. Determining The 'Thickness of Reality'...")
if horizon_L:
    print(f"   Bell Horizon (Coherence Depth): ~{horizon_L} Layers")
else:
    print(f"   Horizon > {N_SIM} Layers (Strongly Entangled)")
    
# --- WRITE REPORT ---
with open(REPORT_FILE, "w") as f:
    f.write("# RAPORT: QW-691 INTER-LAYER HORIZON\n")
    f.write(f"**Data:** {datetime.datetime.now()}\n")
    f.write("**Cel:** Wyznaczenie 'Głębokości Koherencji' (Ile warstw jest kwantowych?)\n\n")
    
    f.write("## 1. WYNIKI\n")
    f.write("| Depth L | S (CHSH) | Status |\n")
    f.write("|---------|----------|--------|\n")
    for r in results:
        status = "Quantum" if r[1] > 2.0 else "CLASSICAL"
        f.write(f"| {r[0]} | {r[1]:.4f} | {status} |\n")
        
    f.write(f"\n**Horyzont (N_max):** {horizon_L if horizon_L else '>12'} warstw\n\n")
    
    f.write("## 2. INTERPRETACJA\n")
    if horizon_L:
        f.write(f"Dla obserwatora na powierzchni, rzeczywistość jest kwantowa tylko do głębokości **{horizon_L} warstw**.\n")
        f.write("Głębiej (L > N_max) fluktuacje są uśrednione do szumu klasycznego.\n")
        f.write("To definiuje fizyczną 'Grubość Teraźniejszości' w modelu fraktalnym.\n")

print(f"   Report saved to: {REPORT_FILE}")
