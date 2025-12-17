#!/usr/bin/env python3
"""
QW-1518: CHIRP SIGNAL (BINARY MERGER SIMULATION)
=================================================
With PNG visualization of chirp waveform.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import fft, fftfreq
import datetime

print("="*80)
print("QW-1518: CHIRP SIGNAL (BINARY MERGER SIMULATION)")
print("="*80)

# Parameters
ALPHA_GEO = 4 * np.log(2)
PHI = np.pi / 6
BETA_TORS = 0.01

def K(d):
    if d < 0.1:
        d = 0.1
    return ALPHA_GEO * np.cos(np.pi/4 * d + PHI) / (1 + BETA_TORS * d)

c_tors = np.sqrt(abs(K(0)))

# 1D chain
N = 400
L = 100.0
dx = L / N

# Chirp parameters
f_init = 0.02
f_final = 0.2
t_merger = 150.0
chirp_amplitude = 1.0

print(f"Chirp: f_0 = {f_init} -> f_merger = {f_final}")

# Time
dt = 0.1
t_max = t_merger + 20
n_steps = int(t_max / dt)

# Fields
theta = np.zeros(N)
theta_dot = np.zeros(N)

source_pos = N // 4
detector_pos = 3 * N // 4
detector_history = []
source_history = []
freq_history = []
times = []

print(f"Simulating {n_steps} steps...")

for step in range(n_steps):
    t = step * dt
    times.append(t)
    
    # Chirp frequency
    if t < t_merger:
        tau = 1 - t / t_merger
        if tau > 0.01:
            f_t = f_init * (tau ** (-3/8))
            A_t = chirp_amplitude * (tau ** (-1/4))
        else:
            f_t = f_final
            A_t = chirp_amplitude * 10
    else:
        t_after = t - t_merger
        f_t = f_final * np.exp(-t_after / 10)
        A_t = chirp_amplitude * 10 * np.exp(-t_after / 5)
    
    freq_history.append(f_t)
    omega_gw = 2 * np.pi * f_t
    
    source = np.zeros(N)
    source[source_pos] = A_t * np.sin(omega_gw * t)
    source_history.append(A_t * np.sin(omega_gw * t))
    
    laplacian = np.zeros(N)
    for i in range(1, N-1):
        laplacian[i] = (theta[i+1] - 2*theta[i] + theta[i-1]) / dx**2
    
    theta_ddot = c_tors**2 * laplacian - BETA_TORS * theta_dot + source
    theta_dot += theta_ddot * dt
    theta += theta_dot * dt
    
    detector_history.append(theta[detector_pos])

times = np.array(times)
detector_signal = np.array(detector_history)
source_signal = np.array(source_history)
freq_signal = np.array(freq_history)

print("Simulation complete. Generating visualization...")

# Create figure with 4 subplots
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('QW-1518: Gravitational Wave Chirp Signal', fontsize=14, fontweight='bold')

# 1. Source waveform
ax1 = axes[0, 0]
ax1.plot(times, source_signal, 'b-', linewidth=0.5)
ax1.set_xlabel('Time')
ax1.set_ylabel('Source Amplitude')
ax1.set_title('Source Signal (Chirp)')
ax1.axvline(x=t_merger, color='r', linestyle='--', label='Merger')
ax1.legend()
ax1.grid(True, alpha=0.3)

# 2. Detected waveform
ax2 = axes[0, 1]
ax2.plot(times, detector_signal, 'g-', linewidth=0.5)
ax2.set_xlabel('Time')
ax2.set_ylabel('Detector Amplitude')
ax2.set_title('Detected Signal at Detector')
ax2.axvline(x=t_merger, color='r', linestyle='--', label='Merger')
ax2.legend()
ax2.grid(True, alpha=0.3)

# 3. Frequency evolution
ax3 = axes[1, 0]
ax3.plot(times, freq_signal, 'r-', linewidth=1)
ax3.set_xlabel('Time')
ax3.set_ylabel('Frequency')
ax3.set_title('Frequency Evolution (Chirp)')
ax3.axvline(x=t_merger, color='k', linestyle='--', label='Merger')
ax3.set_yscale('log')
ax3.legend()
ax3.grid(True, alpha=0.3)

# 4. Spectrogram (time-frequency)
ax4 = axes[1, 1]
# Simple spectrogram using windowed FFT
window_size = 100
n_windows = len(times) // window_size
spectrogram = np.zeros((50, n_windows))
freq_bins = np.linspace(0, 0.5, 50)

for w in range(n_windows):
    start = w * window_size
    end = start + window_size
    segment = detector_signal[start:end]
    segment_ac = segment - np.mean(segment)
    spectrum = np.abs(fft(segment_ac)[:50])**2
    spectrogram[:, w] = spectrum

extent = [0, t_max, 0, 0.5]
im = ax4.imshow(spectrogram, aspect='auto', origin='lower', extent=extent, cmap='hot')
ax4.set_xlabel('Time')
ax4.set_ylabel('Frequency')
ax4.set_title('Time-Frequency Spectrogram')
ax4.axvline(x=t_merger, color='w', linestyle='--', linewidth=2)
plt.colorbar(im, ax=ax4, label='Power')

plt.tight_layout()
plt.savefig('QW-1518_Chirp_Visualization.png', dpi=150, bbox_inches='tight')
print("[SAVED] QW-1518_Chirp_Visualization.png")

# Time-frequency analysis
window_size = 200
n_windows = len(times) // window_size
peak_freqs = []
window_times = []

for w in range(n_windows):
    start = w * window_size
    end = start + window_size
    segment = detector_signal[start:end]
    segment_ac = segment - np.mean(segment)
    spectrum = np.abs(fft(segment_ac))**2
    freqs = fftfreq(window_size, dt)
    pos_mask = freqs > 0
    if np.sum(pos_mask) > 0:
        peak_idx = np.argmax(spectrum[pos_mask])
        peak_freq = freqs[pos_mask][peak_idx]
    else:
        peak_freq = 0
    peak_freqs.append(peak_freq)
    window_times.append(times[start + window_size//2])

peak_freqs = np.array(peak_freqs)
window_times = np.array(window_times)

# Chirp detection
valid = peak_freqs > 0.001
if np.sum(valid) >= 3:
    coeffs = np.polyfit(window_times[valid], np.log(peak_freqs[valid]), 1)
    slope = coeffs[0]
else:
    slope = 0

print(f"\nFrequency evolution: d(log f)/dt = {slope:.4f}")

if slope > 0.01:
    verdict = "CHIRP DETECTED"
elif slope > 0:
    verdict = "WEAK CHIRP"
else:
    verdict = "NO CHIRP"

print(f"VERDICT: {verdict}")

# Save report
report = f"""# QW-1518: Chirp Signal

**Date:** {datetime.datetime.now()}

## Visualization
![Chirp](QW-1518_Chirp_Visualization.png)

## Results
- d(log f)/dt = {slope:.4f}
- Verdict: {verdict}
"""

with open("QW-1518_Chirp_Signal.md", "w") as f:
    f.write(report)

print("[SAVED] QW-1518_Chirp_Signal.md")
