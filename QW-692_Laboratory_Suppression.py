#!/usr/bin/env python3
"""
QW-692: LABORATORY PARADOX (FRACTAL SUPPRESSION TEST)
======================================================
Purpose: Test if suppressing activity in higher fractal layers (simulating Lab Isolation/Cooling)
         restores quantum Bell violation (S > 2) from a naturally classical system.

Hypothesis: 
  "Natural State" (Active Layers 0-7) -> S < 2 (Classical Averaging)
  "Lab State" (Suppressed Layers 1-7) -> S > 2 (Quantum Fundamental Revealed)

Method:
  1. Build a Vertical Chain (N=8 Layers) engaged in thermal/fractal noise.
  2. Define "Observer Measurement" as a weighted sum over layers.
     - Natural: Weight ~ 1/L (or uniform) -> Huge averaging.
     - Lab: Weight ~ delta(L,0) -> Filters out noise, measures fundamental mode.
  3. Compare S_natural vs S_lab.
  
Output: RAPORT_QW692_LABORATORY_PARADOX.md
"""

import numpy as np
from scipy.linalg import eigh
import datetime

print("="*80)
print("QW-692: LABORATORY PARADOX TEST")
print("="*80)
print("Testing: Does suppressing fractal noise restore Bell violation?")
print()

# --- Parameters ---
N_LAYERS = 8
J_COUPLING = 1.0

# Pauli Matrices
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
sigma_0 = np.eye(2, dtype=complex)

REPORT_FILE = "RAPORT_QW692_LABORATORY_PARADOX.md"

def build_op(op, site, n_sites):
    result = 1
    for i in range(n_sites):
        if i == site:
            result = np.kron(result, op)
        else:
            result = np.kron(result, sigma_0)
    return result

# --- 1. NATURAL STATE SIMULATION ---
print(f"1. Building Natural System (N={N_LAYERS} active layers)...")
dim_total = 2**N_LAYERS

Sx = [build_op(sigma_x, i, N_LAYERS) for i in range(N_LAYERS)]
Sz = [build_op(sigma_z, i, N_LAYERS) for i in range(N_LAYERS)]

# Critical Ising Chain (Strong Entanglement between layers)
H = np.zeros((dim_total, dim_total), dtype=complex)
for i in range(N_LAYERS - 1):
    H += -1.0 * Sz[i] @ Sz[i+1] # Vertical links
    H += -1.0 * Sx[i]           # Quantum fluctuations
H += -1.0 * Sx[N_LAYERS-1]

print("   Finding Ground State...")
evals, evecs = eigh(H)
psi_ground = evecs[:, 0]

# --- 2. DEFINING OBSERVABLES ---
theta_a, theta_ap = 0, np.pi/2
theta_b, theta_bp = np.pi/4, 3*np.pi/4

def get_op_natural(theta):
    # Natural Observer averages over ALL layers weighted by 1/(L+1)
    # This simulates coupling to the whole depth.
    op_accum = np.zeros((dim_total, dim_total), dtype=complex)
    norm = 0
    for l in range(N_LAYERS):
        weight = 1.0 / (l + 1) # Decaying weight
        op_site = np.sin(theta)*Sx[l] + np.cos(theta)*Sz[l]
        op_accum += weight * op_site
        norm += weight
    return op_accum / norm

def get_op_lab(theta):
    # Lab Observer effectively measures ONLY Layer 0 (Fundamental Mode)
    # The "Isolator" (Vacuum/Cryostat) suppresses L>0 contribution to 0.
    op_site = np.sin(theta)*Sx[0] + np.cos(theta)*Sz[0]
    return op_site

def measure_S(psi, method="lab"):
    def corr(thA, thB):
        if method == "natural":
            opA = get_op_natural(thA)
            opB = get_op_natural(thB) # Self-correlation in nature? No, usually A and B are distinct.
            # Wait. Bell test needs TWO PARTS.
            # In QW-690, we measured correlation between Layer 0 and Layer L.
            # Here we need an "A" and "B" experiment.
            # Let's say we split the layers? No.
            # A and B are usually spatially separated.
            # Let's say A measures Layer 0, B measures Layer 0 (Lab).
            # Natural: A measures Avg(Layers A), B measures Avg(Layers B).
            # We need TWO CHAINS for a proper Bell test?
            pass
        else:
            pass
        return 0
    return 0
    
# CRITICAL CORRECTION:
# Bell test requires SPATIAL separation.
# Vertical layers are "internal dimensions".
# A standard Bell test uses 2 entangled particles (Particle A and Particle B).
# EACH particle has N Layers.
# So we need TWO CHAINS (Total 2*N qubits).
# 
# Lab: Measures Particle A (Layer 0) and Particle B (Layer 0).
# Natural: Measures Particle A (Avg Layers) and Particle B (Avg Layers).
pass

print("   Resizing simulation to 2 Chains of N=4 Layers (16 qubits total? Too big).")
# N=3 per chain -> 6 total.
N_PER_CHAIN = 4
N_TOTAL = 2 * N_PER_CHAIN
print(f"   Simulating 2 Entangled Particles, each with {N_PER_CHAIN} Layers (Total {N_TOTAL} qubits).")

dim_total = 2**N_TOTAL
Sx_2 = [build_op(sigma_x, i, N_TOTAL) for i in range(N_TOTAL)]
Sz_2 = [build_op(sigma_z, i, N_TOTAL) for i in range(N_TOTAL)]

# Hamiltonian: Two vertical chains, ENTANGLED at the bottom (Layer 0)
H_2 = np.zeros((dim_total, dim_total), dtype=complex)

# Chain A (qubits 0..N-1), Chain B (qubits N..2N-1)
# Chain A (qubits 0..N-1), Chain B (qubits N..2N-1)
# Vertical links within A (STRONG coupling to environment - Natural State)
J_VERTICAL = 1.0 
for i in range(N_PER_CHAIN - 1):
    H_2 += -1.0 * J_VERTICAL * Sz_2[i] @ Sz_2[i+1] 
    H_2 += -1.0 * Sx_2[i]

# Vertical links within B
offset = N_PER_CHAIN
for i in range(N_PER_CHAIN - 1):
    H_2 += -1.0 * J_VERTICAL * Sz_2[offset+i] @ Sz_2[offset+i+1]
    H_2 += -1.0 * Sx_2[offset+i]
    
# ENTANGLEMENT SOURCE: Standard interaction
J_PAIR = 1.5 
H_2 += +J_PAIR * (Sx_2[0] @ Sx_2[offset] + Sz_2[0] @ Sz_2[offset])

print("   Finding Ground State of Entangled Fractal Pair...")
evals, evecs = eigh(H_2)
psi_ground = evecs[:, 0]

# --- MEASUREMENT FUNCTION CORRECTED ---
def measure_S_corrected(psi, mode="lab"):
    
    def get_effective_op(theta, chain_offset):
        # chain_offset is 0 for A, N_PER_CHAIN for B
        if mode == "lab":
            # Lab measures only Layer 0
            idx = chain_offset + 0
            return np.sin(theta)*Sx_2[idx] + np.cos(theta)*Sz_2[idx]
        elif mode == "natural":
            # Natural measures average over all layers
            op_accum = np.zeros((dim_total, dim_total), dtype=complex)
            norm = 0
            for l in range(N_PER_CHAIN):
                # weight = 1.0 / (l+1) # Damping
                weight = 1.0 # Uniform interaction with all layers?
                # Actually "Natural" usually means we can't distinguish layers, so we sum them.
                # Or we interact with the "Bulk". Let's try uniform.
                idx = chain_offset + l
                site_op = np.sin(theta)*Sx_2[idx] + np.cos(theta)*Sz_2[idx]
                op_accum += weight * site_op
                norm += weight
            return op_accum / norm

    def corr(thA, thB):
        opA = get_effective_op(thA, 0) # Chain A
        opB = get_effective_op(thB, N_PER_CHAIN) # Chain B
        obs = opA @ opB
        return np.real(psi.conj().T @ obs @ psi)

    E1 = corr(theta_a, theta_b)
    E2 = corr(theta_a, theta_bp)
    E3 = corr(theta_ap, theta_b)
    E4 = corr(theta_ap, theta_bp)
    
    return np.abs(E1 - E2 + E3 + E4)

# --- EXECUTE ---
print("\n2. Measuring Bell Parameter S...")

S_lab = measure_S_corrected(psi_ground, mode="lab")
S_natural = measure_S_corrected(psi_ground, mode="natural")

print(f"   S_natural (Full Fractal Averaging) = {S_natural:.4f}")
print(f"   S_lab     (Layer 0 Suppression)    = {S_lab:.4f}")

# --- REPORT ---
status = "✅ SUCCESS" if (S_natural < 2.0 and S_lab > 2.0) else "❌ INCONCLUSIVE"
lab_gain = S_lab - S_natural

print(f"\n   STATUS: {status}")
print(f"   Lab restores {lab_gain:.2f} of quantumness.")

with open(REPORT_FILE, "w") as f:
    f.write("# RAPORT: QW-692 LABORATORY PARADOX\n")
    f.write(f"**Data:** {datetime.datetime.now()}\n")
    f.write("**Cel:** Czy 'wyciszenie' warstw fraktalnych (Lab) przywraca łamanie Bella?\n\n")
    f.write("## 1. WYNIKI\n")
    f.write(f"- **System:** 2 Particles x {N_PER_CHAIN} Layers\n")
    f.write(f"- **S_natural** (Average): {S_natural:.4f}\n")
    f.write(f"- **S_lab** (Layer 0):     {S_lab:.4f}\n\n")
    
    f.write("## 2. WNIOSEK\n")
    if S_lab > 2.0 and S_natural < 2.0:
        f.write("Potwierdzono: Izolacja laboratoryjna działa jako filtr modów.\n")
        f.write("W naturze (Average) kwantowość ginie w szumie warstw.\n")
        f.write("W laboratorium (Layer 0) kwantowość jest 'odzyskana'.\n")
    else:
        f.write(f"Wynik niejednoznaczny (S_lab={S_lab:.2f}). Może coupling między A-B jest za słaby?\n")

print(f"   Report saved to: {REPORT_FILE}")
