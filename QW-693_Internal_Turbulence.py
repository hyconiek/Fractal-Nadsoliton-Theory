#!/usr/bin/env python3
"""
QW-693: INTERNAL TURBULENCE VERIFICATION
========================================
Purpose: Re-investigate Hypothesis H2 ("Turbulent Ether") from the Internal Observer's perspective.
         Uses full non-local K(d) kernel instead of local Laplacian.

Method:
  1. Simulate 1D chain of Nadsolitons with K(d) coupling.
  2. Measure Internal Reynolds Number (Re_int).
  3. Measure Structure Function S2(r) to check for Kolmogorov scaling (r^2/3).

Metric for Turbulence:
  - Re_int > 1000
  - S2(r) ~ r^(2/3) (Inertial subrange)
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime

print("="*80)
print("QW-693: INTERNAL TURBULENCE VERIFICATION")
print("="*80)

# --- Parameters ---
N = 1024                # System size (large for spectral analysis)
DT = 0.01
STEPS = 2000
EQUILIBRATION = 500
ALPHA = 4 * np.log(2)   # 2.7726
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.001       # Low dissipation to allow turbulence

print(f"System Size: {N}")
print(f"Beta (Dissipation): {BETA_TORS}")
print(f"Alpha (Coupling): {ALPHA:.4f}")

# --- Initialize Field (Random Fluctuations) ---
np.random.seed(693)
psi = np.random.randn(N) + 1j * np.random.randn(N)
psi /= np.sqrt(np.mean(np.abs(psi)**2))

# --- Precompute K(d) Kernel (FFT Convolution) ---
# We use FFT for O(N log N) convolution instead of O(N^2)
d = np.arange(N)
# Distance on ring (periodic)
d_eff = np.minimum(d, N - d)
K_kernel = ALPHA * np.cos(OMEGA * d_eff + PHI) / (1 + BETA_TORS * d_eff)
K_hat = np.fft.fft(K_kernel)

print("Kernel initialized.")

# --- Dynamics ---
# dPsi/dt = -i [ V_NL * Psi ] - Dissipation
# V_NL(x) = Convolution(K, |Psi|^2) (Non-local potential)
# Plus some kinetic term? The K(d) IS the kinetic/potential mix.
# Let's model the specific Nadsoliton equation:
# i dPsi/dt = - K * Psi  (Linear non-local dispersion) + Nonlinearity?
# To get turbulence, we need non-linearity.
# FIN Theory: Self-interaction is Non-linear phase modulation?
# Let's use: i dPsi/dt = (K * Psi) + NonLinearity * Psi
# Standard NLS-like but with non-local K.
# Let's assume the "Turbulence" comes from the Advection-like term in phase space.
# Model: i dPsi/dt = Convolve(K, Psi) + |Psi|^2 Psi

def get_rhs(psi_current):
    # Non-local Kinetic/Potential term
    psi_hat = np.fft.fft(psi_current)
    # Convolution theorem: FFT(f * g) = FFT(f) * FFT(g)
    # But here interaction is likely V(x) = Int K(x-y) |Psi(y)|^2 dy ?
    # Or linear dispersion?
    # Inspecting H2 context: "Turbulent Ether" implies Flow. 
    # Let's implement the FULL interaction:
    # Hamiltonian = Sum K_ij Psi*_i Psi_j + Sum |Psi_i|^4
    # EoM: i dPsi_i/dt = Sum_j K_ij Psi_j + 2|Psi_i|^2 Psi_i
    
    # Linear Non-local term (Dispersion)
    linear_term_hat = K_hat * psi_hat
    linear_term = np.fft.ifft(linear_term_hat)
    
    # Nonlinear term (Self-interaction / pressure)
    nonlinear_term = 1.0 * (np.abs(psi_current)**2) * psi_current
    
    # Dissipation (necessary for steady state turbulence)
    # Drag term: -i * gamma * Psi (amplitude decay)
    # But we want Re number estimate.
    # Let's add small viscosity to high-k modes?
    # Or simple damping balanced by forcing? 
    # We will just observe decaying turbulence or use conserved energy.
    # Let's use conserved Hamiltonian evolution and see if energy cascades.
    
    return -1j * (linear_term + nonlinear_term)

# --- Evolution Loop ---
print(f"Evolving for {STEPS} steps...")

snapshots_v = []
timestamps = []

for t in range(STEPS):
    # RK4 Integration
    k1 = get_rhs(psi)
    k2 = get_rhs(psi + 0.5 * DT * k1)
    k3 = get_rhs(psi + 0.5 * DT * k2)
    k4 = get_rhs(psi + DT * k3)
    
    psi += (DT / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
    
    # Normalize (conserver particle number approx)
    # psi /= np.sqrt(np.mean(np.abs(psi)**2)) 
    # Actually, let's let energy flow naturally, only normalize if unstable.
    nm = np.sqrt(np.mean(np.abs(psi)**2))
    if nm > 10 or nm < 0.1: psi /= nm # Soft clamp
    
    if t > EQUILIBRATION and t % 10 == 0:
        # Calculate "Velocity" field
        # In Quantum Hydrodynamics (Madelung):
        # Psi = sqrt(rho) * exp(i S)
        # v = grad S
        phase = np.angle(psi)
        # Unwrap phase to handle 2pi jumps
        phase_unwrapped = np.unwrap(phase)
        v = np.gradient(phase_unwrapped)
        snapshots_v.append(v)
        timestamps.append(t)

snapshots_v = np.array(snapshots_v)
print(f"Collected {len(snapshots_v)} velocity snapshots.")

# --- Analysis: Internal Reynolds Number ---
# Re = L * V / nu
# L = Characteristic Length (Correlation length?)
# V = RMS Velocity
# nu = Dissipation (Beta_tors)

# 1. RMS Velocity
v_rms = np.sqrt(np.mean(snapshots_v**2))
print(f"RMS Velocity v_rms: {v_rms:.4f}")

# 2. Correlation Length L
# C(r) = <v(x)v(x+r)>
v_fluct = snapshots_v - np.mean(snapshots_v, axis=1, keepdims=True)
corr_sum = np.zeros(N)
count = 0

for v in v_fluct:
    # Autocorrelation via FFT
    ft = np.fft.fft(v)
    spec = np.abs(ft)**2
    cor = np.fft.ifft(spec).real
    cor /= cor[0] # Normalize
    corr_sum += cor
    count += 1
    
avg_corr = corr_sum / count
# Find first zero crossing or 1/e decay
L_corr = 1.0
for i in range(N // 2):
    if avg_corr[i] < 0.37: # 1/e
        L_corr = float(i)
        break

print(f"Correlation Length L_corr: {L_corr:.2f}")

# 3. Reynolds Number
# Effective viscosity ~ Beta_tors (from theory)
nu_eff = BETA_TORS
Re_int = (L_corr * v_rms) / nu_eff

print(f"Internal Reynolds Number Re_int: {Re_int:.2f}")

# --- Analysis: Structure Function S2(r) ---
# S2(r) = <|v(x+r) - v(x)|^2>
# Kolmogorov: S2(r) ~ r^(2/3)

print("Calculating Structure Function S2(r)...")
r_vals = np.arange(1, N // 10) # Investigate small scales
S2_vals = []

for r in r_vals:
    # Shift array
    diffs = []
    for v in snapshots_v:
        # Periodic difference
        d = np.abs(np.roll(v, -r) - v)
        # Only take valid range (ignore wrap-around artifact? Roll handles it)
        diffs.append(np.mean(d**2))
    S2_vals.append(np.mean(diffs))

S2_vals = np.array(S2_vals)

# Fit Power Law: S2 ~ r^zeta
log_r = np.log(r_vals)
log_S2 = np.log(S2_vals)

coeffs = np.polyfit(log_r, log_S2, 1)
zeta = coeffs[0]

print(f"Structure Function Scaling exponent zeta: {zeta:.4f}")
print("Expected Kolmogorov: 0.666 (2/3)")

# --- Verification Logic ---
scaling_error = abs(zeta - 0.666) / 0.666 * 100
is_turbulent = (Re_int > 1000) and (abs(zeta - 0.666) < 0.2)
is_partial = (Re_int > 100) and (abs(zeta - 0.666) < 0.4)

verdict_file = "raport_qw693_internal_turbulence.md"
with open(verdict_file, "w") as f:
    f.write("# RAPORT QW-693: INTERNAL TURBULENCE\n")
    f.write(f"**Date:** {datetime.datetime.now()}\n\n")
    f.write("## 1. Parameters\n")
    f.write(f"- N: {N}\n")
    f.write(f"- Beta (Viscosity): {BETA_TORS}\n")
    f.write(f"- Alpha (Coupling): {ALPHA:.4f}\n\n")
    
    f.write("## 2. Results\n")
    f.write(f"- **v_rms:** {v_rms:.4f}\n")
    f.write(f"- **L_corr:** {L_corr:.2f}\n")
    f.write(f"- **Re_int:** {Re_int:.2f}\n")
    f.write(f"- **Scaling ζ:** {zeta:.4f} (Expected 0.67)\n\n")
    
    f.write("## 3. Verdict\n")
    if is_turbulent:
        f.write("### ✅ TURBULENCE CONFIRMED\n")
        f.write("The internal degrees of freedom exhibit high Reynolds number and Kolmogorov-like scaling.\n")
        f.write("H2 (Turbulent Ether) is VALID internally.\n")
    elif is_partial:
        f.write("### 🟡 PARTIAL TURBULENCE\n")
        f.write("Signs of turbulence present, but scaling is not perfect.\n")
    else:
        f.write("### ❌ LAMINAR / NOISE\n")
        f.write("No evidence of turbulence. System is dominated by dispersion or dissipation.\n")

print("-" * 40)
print(f"Report saved: {verdict_file}")
print("VERDICT:")
if is_turbulent:
    print("✅ TURBULENCE CONFIRMED")
elif is_partial:
    print("🟡 PARTIAL TURBULENCE")
else:
    print("❌ LAMINAR / NOISE")
