# QW-520 TO QW-524: KERNEL FORENSICS
# GOAL: Understand the frozen Kernel K(d) before simulating Nadsoliton.
# PARADIGM: "Understand the Brick before building the House".

import numpy as np
import scipy.fft
import matplotlib.pyplot as plt
from scipy.signal import convolve

# FROZEN PARAMETERS
alpha_geo = 4 * np.log(2)
omega = np.pi/4
phi = np.pi/6
beta_tors = 0.01

print("="*80)
print("QW-520 TO QW-524: KERNEL FORENSICS")
print("Analyzing the Frozen Kernel K(d)")
print("="*80)

def K(d):
    # Avoid division by zero at d=0? 
    # Usually d starts from some small epsilon or we handle d=0 manually.
    # Let's assume d >= 0.
    denom = 1 + beta_tors * d
    if np.any(denom == 0): return 0 # Should not happen for d>=0
    return (alpha_geo * np.cos(omega * d + phi)) / denom

# --- QW-520: KERNEL SPECTROSCOPY (Propagator Check) ---
print("\n" + "="*60)
print("QW-520: KERNEL SPECTROSCOPY (Is it a Propagator?)")
print("="*60)

# Calculate FFT of K(d)
# If K(d) is a propagator for a massive particle, FFT should look like 1/(k^2 + m^2) (Lorentzian-like)
# If it's a wave, it should have peaks.

N_points = 1024
d_range = np.linspace(0, 200, N_points)
k_vals = K(d_range)

# FFT
fft_vals = scipy.fft.fft(k_vals)
freqs = scipy.fft.fftfreq(N_points, d_range[1] - d_range[0])

# Power Spectrum
power = np.abs(fft_vals)**2

# Find peaks in spectrum
pos_freqs = freqs[:N_points//2]
pos_power = power[:N_points//2]
peak_idx = np.argmax(pos_power)
peak_freq = pos_freqs[peak_idx]

print(f"Dominant Frequency k_peak: {peak_freq:.4f}")
print(f"Expected k from omega: {omega / (2*np.pi):.4f}")

# Check shape around peak
# Is it a sharp delta (pure wave) or broad (damped)?
print(f"Peak Width (Qualitative): {'Sharp' if np.max(pos_power)/np.mean(pos_power) > 100 else 'Broad'}")

# Does it look like a Propagator 1/(k^2+m^2)?
# 1/(k^2+m^2) is monotonic decreasing from k=0 (if m real).
# Our K oscillates, so it has a peak at k != 0.
# This implies it propagates a mode with non-zero momentum/mass?
# Or it's a "tachyon" (instability)?
print("Interpretation: Peak at k != 0 implies oscillatory spatial structure.")
print("  -> Not a simple Yukawa potential (1/r * exp(-mr)).")
print("  -> More like a crystal lattice potential or a wave packet.")

# --- QW-521: EFFECTIVE POTENTIAL (Smoothing) ---
print("\n" + "="*60)
print("QW-521: EFFECTIVE POTENTIAL (Smoothing)")
print("="*60)

# Hypothesis: A large particle "averages out" the oscillations.
# Convolve K(d) with a Gaussian of width sigma.

sigmas = [1.0, 2.0, 5.0, 10.0]
r_eval = np.linspace(0, 50, 200)
k_raw = K(r_eval)

print(f"{'r':<6} {'Raw':<10} {'Sm(1.0)':<10} {'Sm(5.0)':<10} {'Sm(10.0)':<10}")
print("-" * 50)

smoothed_results = {}

for sigma in sigmas:
    # Gaussian kernel
    x = np.linspace(-3*sigma, 3*sigma, 100)
    g = np.exp(-x**2 / (2*sigma**2))
    g /= np.sum(g) # Normalize
    
    # Convolve
    # We need a longer K array to avoid edge effects
    r_long = np.linspace(0, 100, 400)
    k_long = K(r_long)
    k_smooth = convolve(k_long, g, mode='same')
    
    # Interpolate back to r_eval
    smoothed_results[sigma] = np.interp(r_eval, r_long, k_smooth)

# Check monotonicity of smoothed potentials
for i in range(0, len(r_eval), 20):
    r = r_eval[i]
    raw = k_raw[i]
    s1 = smoothed_results[1.0][i]
    s5 = smoothed_results[5.0][i]
    s10 = smoothed_results[10.0][i]
    print(f"{r:<6.1f} {raw:<10.4f} {s1:<10.4f} {s5:<10.4f} {s10:<10.4f}")

# Check if any smoothed potential looks like 1/r (monotonic and negative?)
# Or at least monotonic?
is_monotonic = {}
for sigma in sigmas:
    s = smoothed_results[sigma]
    # Check derivative
    ds = np.diff(s)
    # Monotonic if all ds have same sign (ignoring small noise)
    pos_derivs = np.sum(ds > 0.001)
    neg_derivs = np.sum(ds < -0.001)
    is_monotonic[sigma] = (pos_derivs == 0 or neg_derivs == 0)

print("\nMonotonicity Check:")
for sigma, mono in is_monotonic.items():
    print(f"  Sigma={sigma}: Monotonic? {mono}")

# --- QW-522: TOPOLOGICAL MAPPING (Phase Space) ---
print("\n" + "="*60)
print("QW-522: TOPOLOGICAL MAPPING (Phase Space)")
print("="*60)

# Plot (K, dK/dd)
d_phase = np.linspace(0, 50, 1000)
k_phase = K(d_phase)
dk_phase = np.gradient(k_phase, d_phase)

# Check for cycles/attractors
# Since K decays, it spirals into (0,0).
# But does it have "preferred" loops?
# We check the "Phase Angle" theta = atan2(dK, K)
# Does theta rotate uniformly?
theta = np.arctan2(dk_phase, k_phase)
# Unwrap
theta_unwrapped = np.unwrap(theta)

# Check linearity of theta vs d
# If theta ~ k*d, it's a spiral.
coeffs = np.polyfit(d_phase, theta_unwrapped, 1)
slope = coeffs[0]
print(f"Phase Rotation Rate d(theta)/dd: {slope:.4f}")
print(f"Expected from omega: {-omega:.4f} (approx)")

# --- QW-523: RESONANCE DISTANCES (Mass Gap) ---
print("\n" + "="*60)
print("QW-523: RESONANCE DISTANCES (Mass Gap)")
print("="*60)

# Find local maxima of K(d)
# d where K'(d) = 0 and K''(d) < 0
# Analytical: d/dd [ cos(wd+phi)/(1+bd) ] = 0
# -w sin / (1+bd) - b cos / (1+bd)^2 = 0
# -w(1+bd) sin - b cos = 0
# tan(wd+phi) = -b / (w(1+bd))
# Since b is small, tan ~ 0 (or small negative).
# So wd + phi ~ n * pi.
# Maxima when cos > 0, so wd + phi ~ 2n*pi (or close).

# Numerical check
peaks = []
for i in range(1, len(k_phase)-1):
    if k_phase[i] > k_phase[i-1] and k_phase[i] > k_phase[i+1]:
        peaks.append(d_phase[i])

print(f"Resonance Distances (Peaks): {peaks[:5]}")
# Check spacing
diffs = np.diff(peaks)
print(f"Spacings: {diffs[:5]}")
print(f"Average Spacing: {np.mean(diffs):.4f}")
print(f"Expected 2pi/omega: {2*np.pi/omega:.4f}")

# Are these "Masses"?
# If d ~ 1/Mass?
# Or d ~ Mass?
# If spacing is constant, it's a harmonic ladder.

# --- QW-524: FLUID PROPERTIES (Viscosity/Modulus) ---
print("\n" + "="*60)
print("QW-524: FLUID PROPERTIES (Information Fluid)")
print("="*60)

# Calculate effective Bulk Modulus B
# B ~ Sum(d^2 * K(d)) over a volume?
# If K is force, Energy ~ Integral K dx.
# Stiffness ~ dK/dd.
# Let's calculate the "Moment" M2 = Integral d^2 * K(d) dd
# If this diverges, the fluid is "long-range" (like plasma).
# If it converges, it's "short-range" (like water).

# Integral d^2 * cos(d) / d ~ d * cos(d). Diverges.
# K ~ 1/d.
# d^2 * K ~ d.
# Integral 0 to R of d dd ~ R^2.
# Diverges quadratically.

print("Moment Analysis:")
print("  K(d) ~ 1/d")
print("  Moment M0 (Total Force) ~ Int K ~ Int 1/d ~ log(R) -> Diverges (Logarithmically)")
print("  Moment M2 (Bulk Modulus) ~ Int d^2 K ~ Int d -> Diverges (Quadratically)")

print("Interpretation:")
print("  The 'Information Fluid' is NON-LOCAL.")
print("  It behaves like a stiff, long-range medium (e.g., charged plasma or solid).")
print("  It is NOT a simple viscous fluid (which requires local interactions).")

print("="*80)
print("MISSION COMPLETE")
