import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
from scipy.fft import fft, fftfreq

# QW-1508: GRAVITATIONAL WAVES - PARAMETER SCAN
# QW-1507 failed due to overdamping. Here we systematically scan
# the damping coefficient to find the propagation regime.

ALPHA_GEO = 4 * np.log(2)
OMEGA = np.pi / 4

def run_wave_parameter_scan():
    print("=" * 80)
    print("QW-1508: GRAVITATIONAL WAVES - PARAMETER SCAN")
    print("=" * 80)
    
    N = 100
    L = 10.0
    dx = L / N
    
    K_0 = ALPHA_GEO * np.cos(OMEGA * 0 + np.pi/6)
    k_spring = abs(K_0)
    m_node = 1.0
    c_theoretical = dx * np.sqrt(k_spring / m_node)
    
    source_node = N // 4
    source_freq = 0.5
    source_amplitude = 1.0  # INCREASED from 0.1
    
    # SCAN damping values
    damping_values = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1]
    
    print(f"Theoretical wave speed: c_th = {c_theoretical:.4f}")
    print(f"Source frequency: f = {source_freq}")
    print(f"Source amplitude: A = {source_amplitude}")
    print(f"\nScanning damping coefficients: {damping_values}")
    print("-" * 80)
    
    results = []
    
    for damping in damping_values:
        
        def wave_dynamics(t, state):
            phi = state[:N]
            dphi = state[N:]
            
            d2phi = np.zeros(N)
            for i in range(1, N-1):
                d2phi[i] = (k_spring / m_node) * (phi[i+1] - 2*phi[i] + phi[i-1])
            
            d2phi[0] = 0
            d2phi[N-1] = 0
            d2phi[source_node] += source_amplitude * np.sin(2 * np.pi * source_freq * t)
            d2phi -= damping * dphi
            
            return np.concatenate([dphi, d2phi])
        
        t_max = 100.0
        n_points = 1000
        t_eval = np.linspace(0, t_max, n_points)
        initial_state = np.zeros(2 * N)
        
        sol = solve_ivp(wave_dynamics, [0, t_max], initial_state,
                        method='RK45', t_eval=t_eval, rtol=1e-6)
        
        phi_history = sol.y[:N, :]
        times = sol.t
        dt = times[1] - times[0]
        
        # Detect at far node
        det_node = 3 * N // 4
        signal = phi_history[det_node, :]
        signal_ac = signal - np.mean(signal)
        
        freqs = fftfreq(len(times), dt)
        spectrum = np.abs(fft(signal_ac))**2
        
        pos_mask = freqs > 0
        peak_idx = np.argmax(spectrum[pos_mask])
        peak_freq = freqs[pos_mask][peak_idx]
        peak_power = spectrum[pos_mask][peak_idx]
        
        # Check frequency match
        freq_match = abs(peak_freq - source_freq) < 0.15
        
        # Measure amplitude at detector
        amplitude = np.std(signal_ac)
        
        status = "✅ DETECTED" if freq_match else "❌ NOT DETECTED"
        print(f"β = {damping:.3f}: f_peak = {peak_freq:.3f}, A = {amplitude:.4f} -> {status}")
        
        results.append({
            "damping": damping,
            "peak_freq": peak_freq,
            "amplitude": amplitude,
            "detected": freq_match
        })
    
    # Analysis
    print("\n" + "=" * 80)
    print("QW-1508 SUMMARY")
    print("=" * 80)
    
    detected_runs = [r for r in results if r["detected"]]
    
    if detected_runs:
        print(f"✅ GRAVITATIONAL WAVES DETECTED in {len(detected_runs)}/{len(results)} runs")
        for r in detected_runs:
            print(f"   β = {r['damping']:.3f}: f = {r['peak_freq']:.3f}, A = {r['amplitude']:.4f}")
        
        # Find optimal damping (lowest damping where detected)
        optimal = min(detected_runs, key=lambda x: x['damping'])
        conclusion = f"Waves propagate when β ≤ {optimal['damping']:.3f}"
    else:
        print("❌ GRAVITATIONAL WAVES NOT DETECTED in any configuration")
        conclusion = "Model fundamentally does not support wave propagation"
    
    # Critical damping analysis
    # For wave propagation: damping < 2*sqrt(k*m) per node
    critical_damping = 2 * np.sqrt(k_spring * m_node)
    print(f"\nCritical damping (per node): γ_c = {critical_damping:.4f}")
    print(f"β_tors (FIN parameter): {0.01:.4f}")
    print(f"Ratio β/γ_c = {0.01 / critical_damping:.4f}")
    
    if 0.01 / critical_damping < 1.0:
        wave_support = "UNDERDAMPED (Waves should propagate)"
    else:
        wave_support = "OVERDAMPED (Waves suppressed)"
    print(f"Regime: {wave_support}")
    
    # Save report
    report = f"""# QW-1508: Gravitational Waves - Parameter Scan

**Date:** 2025-12-17

## Scan Results
| Damping β | Peak Freq | Amplitude | Detected? |
|-----------|-----------|-----------|-----------|
"""
    for r in results:
        status = "Yes" if r["detected"] else "No"
        report += f"| {r['damping']:.3f} | {r['peak_freq']:.3f} | {r['amplitude']:.4f} | {status} |\n"
    
    report += f"""
## Analysis
- Critical damping: γ_c = {critical_damping:.4f}
- FIN β_tors = 0.01
- Regime: {wave_support}

## Conclusion
{conclusion}

## Physical Interpretation
{"The Nadsoliton vacuum supports gravitational wave propagation when torsion damping is below critical value." if detected_runs else "The current model requires modification to support wave propagation. Possible issues: wrong coupling form, missing compressibility term, or fundamental incompatibility."}
"""
    
    with open("QW-1508_Wave_Parameter_Scan.md", "w") as f:
        f.write(report)
    
    print(f"\n[SAVED] QW-1508_Wave_Parameter_Scan.md")

if __name__ == "__main__":
    run_wave_parameter_scan()
