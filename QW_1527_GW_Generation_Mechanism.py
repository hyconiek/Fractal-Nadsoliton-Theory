import numpy as np
import matplotlib.pyplot as plt

def correlation_tensor_Qij(masses, positions):
    """
    Calculates the Correlation Tensor Q_ij for a set of point-like nadsolitons.
    Q_ij = Sum_a M_a (x_ai x_aj - 1/3 delta_ij r_a^2)
    This is formally identical to the quadrupole moment but interpreted as info density.
    """
    Q = np.zeros((3, 3))
    for m, r_vec in zip(masses, positions):
        r_sq = np.dot(r_vec, r_vec)
        for i in range(3):
            for j in range(3):
                delta = 1.0 if i == j else 0.0
                Q[i, j] += m * (r_vec[i] * r_vec[j] - (1.0/3.0) * delta * r_sq)
    return Q

def simulate_binary_evolution(m1, m2, initial_separation, orbits=5.0, steps=1000):
    """
    Simulates a binary system using effective Keplerian dynamics (proven in QW-470).
    Returns the time series of the correlation tensor.
    """
    G_eff = 1.0  # Normalized effective gravity
    mu = (m1 * m2) / (m1 + m2)
    M_tot = m1 + m2
    
    # Keplerian Orbital Frequency (Circular)
    Omega = np.sqrt(G_eff * M_tot / initial_separation**3)
    period = 2 * np.pi / Omega
    
    times = np.linspace(0, orbits * period, steps)
    dt = times[1] - times[0]
    
    Q_series = []
    
    # Center of Mass Frame
    r1_dist = (m2 / M_tot) * initial_separation
    r2_dist = (m1 / M_tot) * initial_separation
    
    for t in times:
        phase = Omega * t
        
        # Positions in xy plane
        pos1 = np.array([r1_dist * np.cos(phase), r1_dist * np.sin(phase), 0.0])
        pos2 = np.array([r2_dist * np.cos(phase + np.pi), r2_dist * np.sin(phase + np.pi), 0.0])
        
        Q = correlation_tensor_Qij([m1, m2], [pos1, pos2])
        Q_series.append(Q)
        
    return np.array(Q_series), times, Omega

# --- VERIFICATION 1: Is Q_ij oscillatory? ---
print("Running QW-1527: FIN GW Generation Verification")
print("="*60)

m1 = 10.0
m2 = 10.0
r_sep = 5.0

Q_t, times, Omega_true = simulate_binary_evolution(m1, m2, r_sep)
print(f"[Binary] M1={m1}, M2={m2}, R={r_sep}")
print(f"[Orbit] Omega={Omega_true:.4f}")

# Calculate second derivative d2Q/dt2 (numerical)
dt = times[1] - times[0]
dQ_dt = np.gradient(Q_t, dt, axis=0)
d2Q_dt2 = np.gradient(dQ_dt, dt, axis=0)

# Amplitude of Q_xx component oscillation
Qxx_amp = (np.max(Q_t[:, 0, 0]) - np.min(Q_t[:, 0, 0])) / 2
d2Qxx_amp = (np.max(d2Q_dt2[:, 0, 0]) - np.min(d2Q_dt2[:, 0, 0])) / 2

print(f"[Result] Q_xx Amplitude: {Qxx_amp:.4f}")
print(f"[Result] d2Q_xx/dt2 Amplitude: {d2Qxx_amp:.4f}")

# Check scaling with Omega^2
ratio = d2Qxx_amp / Qxx_amp
expected_ratio = (2*Omega_true)**2 # GW is at 2*Omega
print(f"[Check] Ratio d2Q/Q = {ratio:.4f}")
print(f"[Check] Expected (2*Omega)^2 = {expected_ratio:.4f}")

if abs(ratio - expected_ratio) / expected_ratio < 0.01:
    print("✅ FREQUENCY CHECK PASSED: Emits at 2*Omega")
else:
    print("❌ FREQUENCY CHECK FAILED")
    
# --- VERIFICATION 2: Chirp Mass Scaling ---
print("\n--- Chirp Mass Scaling Verification ---")
mc_values = []
amplitude_values = []

test_pairs = [
    (10.0, 10.0),
    (20.0, 20.0),
    (30.0, 30.0),
    (10.0, 30.0),
    (5.0, 40.0)
]

print("Simulating different mass pairs...")
for m1_test, m2_test in test_pairs:
    # Use fixed frequency, vary masses (to isolate mass dependence)
    # OR vary masses and let dynamics dictate (more physical) -> Let's use physical dynamics at fixed separation
    Q_t_test, t_test, Om_test = simulate_binary_evolution(m1_test, m2_test, r_sep)
    
    # Calculate Chirp Mass
    mc = (m1_test * m2_test)**(3/5) / (m1_test + m2_test)**(1/5)
    
    # Measure d2Q amplitude (source strength)
    dQ_dt_test = np.gradient(Q_t_test, t_test[1]-t_test[0], axis=0)
    d2Q_dt2_test = np.gradient(dQ_dt_test, t_test[1]-t_test[0], axis=0)
    amp = (np.max(d2Q_dt2_test[:, 0, 0]) - np.min(d2Q_dt2_test[:, 0, 0])) / 2
    
    # We want to check relation: Amp ~ Mc^(5/3) * Omega^(2/3) roughly?
    # Analytical derivation says: d2Q ~ G^(2/3) * Mc^(5/3) * Omega^(2/3)
    # Let's normalize by Omega^(2/3) to isolate Mc^(5/3)
    
    metric = amp / (Om_test**(2/3))
    
    mc_values.append(mc)
    amplitude_values.append(metric)
    print(f"  M1={m1_test:4.1f} M2={m2_test:4.1f} -> Mc={mc:6.3f} | Normalized Amp={metric:8.3f}")

# Fit Power Law: Amp_norm = k * Mc^p
mc_log = np.log(mc_values)
amp_log = np.log(amplitude_values)
slope, intercept = np.polyfit(mc_log, amp_log, 1)

print(f"\n[Scaling Result] Power Law Exponent: {slope:.4f}")
print(f"[Expected]       Exponent: 5/3 = {5/3:.4f}")

if abs(slope - 5/3) < 0.05:
    print("✅ CHIRP MASS SCALING CONFIRMED (Correct 5/3 exponent!)")
    scaling_status = "✅ CONFIRMED"
else:
    print("❌ SCALING MISMATCH")
    scaling_status = "❌ FAILED"

# --- Save Results to Markdown ---
with open("QW_1527_Simulation_Output.md", "w") as f:
    f.write("# QW-1527: Simulation Verification Output\n\n")
    f.write("## 1. Frequency Analysis\n")
    f.write(f"- **Omega (Orbital):** {Omega_true:.4f} rad/s\n")
    f.write(f"- **Q_xx Amplitude:** {Qxx_amp:.4f}\n")
    f.write(f"- **d2Q/dt2 Amplitude:** {d2Qxx_amp:.4f}\n")
    f.write(f"- **Frequency Ratio (d2Q/Q):** {ratio:.4f}\n")
    f.write(f"- **Expected Ratio (2*Omega)^2:** {expected_ratio:.4f}\n")
    if abs(ratio - expected_ratio) / expected_ratio < 0.01:
        f.write("- **Status:** ✅ FREQUENCY CHECK PASSED (Emission at 2*Omega)\n")
    else:
        f.write("- **Status:** ❌ FREQUENCY CHECK FAILED\n")

    f.write("\n## 2. Chirp Mass Scaling Analysis\n")
    f.write("| Mass 1 | Mass 2 | Chirp Mass (Mc) | Normalized Amp |\n")
    f.write("|---|---|---|---|\n")
    for i in range(len(mc_values)):
        f.write(f"| {test_pairs[i][0]:.1f} | {test_pairs[i][1]:.1f} | {mc_values[i]:.3f} | {amplitude_values[i]:.3f} |\n")
    
    f.write(f"\n- **Power Law Exponent (Fitted):** {slope:.4f}\n")
    f.write(f"- **Theoretical Exponent (GR/FIN):** {5/3:.4f}\n")
    f.write(f"- **Status:** {scaling_status}\n")
    
print("\n[Output Saved] Results written to QW_1527_Simulation_Output.md")

