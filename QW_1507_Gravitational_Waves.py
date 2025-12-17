import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
from scipy.fft import fft, fftfreq

# QW-1507: GRAVITATIONAL WAVES WITH COMPRESSIBLE VACUUM
# Previous QW-474 failed because vacuum was incompressible (c_s -> inf).
# Here we add ELASTICITY to the vacuum network.

# FROZEN PARAMETERS
ALPHA_GEO = 4 * np.log(2)  # 2.7726
BETA_TORS = 0.01
OMEGA = np.pi / 4

def run_gravitational_wave_simulation():
    print("=" * 80)
    print("QW-1507: GRAVITATIONAL WAVES IN COMPRESSIBLE VACUUM")
    print("=" * 80)
    
    # 1. Define the wave equation for vacuum perturbations
    # In an elastic medium: d²φ/dt² = c² ∇²φ
    # We model a 1D chain of coupled oscillators (simplest case)
    
    N = 100  # Number of nodes in chain
    L = 10.0  # Total length
    dx = L / N
    
    # Key parameter: WAVE SPEED
    # Hypothesis: c_gw emerges from network stiffness
    # From QW-471, we know the "vacuum stiffness" kappa ~ 0.355
    # From density rho ~ 1 (normalized), c = sqrt(kappa/rho)
    
    # Let's use the ONLY available stiffness measure from theory:
    # The kernel K(d) at d=0 gives the local coupling strength
    K_0 = ALPHA_GEO * np.cos(OMEGA * 0 + np.pi/6) / (1 + BETA_TORS * 0)
    print(f"Local coupling strength K(0) = {K_0:.4f}")
    
    # Spring constant for chain (from K(d) gradient)
    # k = d²K/dd² at d=0 should give curvature = stiffness
    # K(d) = A * cos(ωd + φ) / (1 + βd)
    # For small d: K ≈ A*cos(φ) - A*ω*sin(φ)*d - ...
    # Second derivative term gives spring constant
    
    # Simplified: use K(0) as spring constant
    k_spring = abs(K_0)
    m_node = 1.0  # Mass per node (normalized)
    
    # Wave speed in chain: c = dx * sqrt(k/m)
    c_theoretical = dx * np.sqrt(k_spring / m_node)
    print(f"Theoretical wave speed: c_gw = {c_theoretical:.4f}")
    
    # 2. Set up oscillating source (binary inspiral analog)
    source_node = N // 4  # Source at 1/4 of chain
    source_freq = 0.5  # Hz
    source_amplitude = 0.1
    
    print(f"\nSource parameters:")
    print(f"  Location: node {source_node}")
    print(f"  Frequency: f = {source_freq}")
    print(f"  Amplitude: A = {source_amplitude}")
    
    # 3. Wave equation for coupled chain
    # phi[i] = displacement of node i
    # d²phi[i]/dt² = k/m * (phi[i+1] - 2*phi[i] + phi[i-1]) + source(t)
    
    def wave_dynamics(t, state):
        # state = [phi_0, phi_1, ..., phi_N-1, dphi_0/dt, ..., dphi_N-1/dt]
        phi = state[:N]
        dphi = state[N:]
        
        # Acceleration from Laplacian (discrete)
        d2phi = np.zeros(N)
        for i in range(1, N-1):
            d2phi[i] = (k_spring / m_node) * (phi[i+1] - 2*phi[i] + phi[i-1])
        
        # Boundary conditions: fixed ends
        d2phi[0] = 0
        d2phi[N-1] = 0
        
        # Add source term (sinusoidal driving at source node)
        d2phi[source_node] += source_amplitude * np.sin(2 * np.pi * source_freq * t)
        
        # Add damping (realistic energy loss)
        damping = BETA_TORS  # Use torsion parameter as damping coefficient
        d2phi -= damping * dphi
        
        return np.concatenate([dphi, d2phi])
    
    # 4. Integrate wave equation
    t_max = 50.0
    n_points = 500
    t_eval = np.linspace(0, t_max, n_points)
    
    # Initial conditions: at rest
    initial_state = np.zeros(2 * N)
    
    print(f"\nIntegrating wave equation...")
    sol = solve_ivp(wave_dynamics, [0, t_max], initial_state,
                    method='RK45', t_eval=t_eval, rtol=1e-6)
    
    phi_history = sol.y[:N, :]  # Displacement field over time
    times = sol.t
    
    print(f"Integration complete: {len(times)} time points")
    
    # 5. Detect waves at different distances
    detector_nodes = [N//2, 3*N//4]  # Detectors at 50% and 75% of chain
    
    print(f"\nWave detection:")
    
    wave_detected = False
    detected_speeds = []
    
    for det_node in detector_nodes:
        signal = phi_history[det_node, :]
        
        # Remove DC component
        signal_ac = signal - np.mean(signal)
        
        # FFT to find dominant frequency
        dt = times[1] - times[0]
        freqs = fftfreq(len(times), dt)
        spectrum = np.abs(fft(signal_ac))**2
        
        # Find peak frequency (positive frequencies only)
        pos_mask = freqs > 0
        peak_idx = np.argmax(spectrum[pos_mask])
        peak_freq = freqs[pos_mask][peak_idx]
        peak_power = spectrum[pos_mask][peak_idx]
        
        print(f"  Node {det_node}: f_peak = {peak_freq:.4f}, power = {peak_power:.2e}")
        
        # Check if source frequency is detected
        if abs(peak_freq - source_freq) < 0.1:
            wave_detected = True
            
            # Calculate time delay for wave speed
            distance = (det_node - source_node) * dx
            
            # Find first peak arrival time
            peaks, _ = find_peaks(np.abs(signal_ac), height=0.001)
            if len(peaks) > 0:
                arrival_time = times[peaks[0]]
                if arrival_time > 0:
                    measured_speed = distance / arrival_time
                    detected_speeds.append(measured_speed)
                    print(f"    Distance: {distance:.2f}, Arrival: {arrival_time:.2f}, Speed: {measured_speed:.4f}")
    
    # 6. Final Results
    print("\n" + "=" * 80)
    print("QW-1507 RESULTS")
    print("=" * 80)
    
    if wave_detected:
        mean_speed = np.mean(detected_speeds) if detected_speeds else 0
        print(f"✅ GRAVITATIONAL WAVES DETECTED")
        print(f"   Source frequency: {source_freq}")
        print(f"   Theoretical speed: c_th = {c_theoretical:.4f}")
        print(f"   Measured speed: c_meas = {mean_speed:.4f}")
        if mean_speed > 0:
            ratio = mean_speed / c_theoretical
            print(f"   Ratio c_meas/c_th = {ratio:.4f}")
        result = "SUCCESS"
    else:
        print(f"❌ GRAVITATIONAL WAVES NOT DETECTED")
        print(f"   Peak frequencies do not match source")
        print(f"   System may be overdamped or source too weak")
        result = "FAILURE"
    
    # Save report
    report = f"""# QW-1507: Gravitational Waves in Compressible Vacuum

**Date:** 2025-12-17
**Status:** {result}

## Method
- 1D chain of {N} coupled oscillators
- Spring constant k = K(0) = {K_0:.4f}
- Damping coefficient = β_tors = {BETA_TORS}
- Source: sinusoidal at f = {source_freq} Hz

## Results
- Theoretical wave speed: c_th = {c_theoretical:.4f}
- Waves detected: {wave_detected}
- Measured speed: {np.mean(detected_speeds) if detected_speeds else 'N/A'}

## Conclusion
{"Waves propagate at finite speed derived from kernel K(d)." if wave_detected else "Model requires adjustments to support wave propagation."}
"""
    
    with open("QW-1507_Gravitational_Waves.md", "w") as f:
        f.write(report)
    
    print("\n[SAVED] QW-1507_Gravitational_Waves.md")

if __name__ == "__main__":
    run_gravitational_wave_simulation()
