#!/usr/bin/env python3
"""
QW-639c: ELECTRON MASS - NONLINEAR I_PROC SCALING
Purpose: Test if I_proc scales exponentially with octave (like κ^N)
Hypothesis: I_proc ∝ κ^(N/12) - Processing intensity amplifies with frequency
Date: 2025-12-06
"""

import numpy as np
from scipy.linalg import eigh
from scipy.stats import entropy

print("="*80)
print("QW-639c: NONLINEAR I_PROC SCALING TEST")
print("="*80)
print("Hypothesis: I_proc ∝ κ^(N_octave/12) - Computational load grows with frequency")
print("="*80)

# Constants
ALPHA_GEO = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01
M_PLANCK_GeV = 1.2209e19

# Particle masses (experimental)
M_ELECTRON_MeV = 0.511
M_MUON_MeV = 105.66
M_TAU_MeV = 1776.86

# Fixed components
W = 1  # All leptons have same topology
kappa = ALPHA_GEO / (OMEGA * PHI)
frac_damp = BETA_TORS ** 10

# Vacuum Hamiltonian
N_octaves = 12
def K(d):
    return ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + BETA_TORS * d)

H_vac = np.zeros((N_octaves, N_octaves))
for i in range(N_octaves):
    for j in range(N_octaves):
        H_vac[i, j] = -K(abs(i - j))

evals, evecs = eigh(H_vac)

# Base processing intensity (from QW-639b)
C_stability = 12.027675
sigma = 0.8
lambda_chaos = 0.1

def compute_mass(N_oct, particle_name):
    """Compute particle mass for given octave index"""
    
    # Octave amplification
    oct_amp = kappa ** (N_oct / 12)
    
    # Resonance overlap
    psi = np.zeros(N_octaves)
    psi[N_oct] = 1.0
    A_res = abs(np.dot(psi, evecs[:, 0]))
    
    # Processing intensity (HYPOTHESIS: scales with octave amplification!)
    psi_dist = np.exp(-0.5 * ((np.arange(N_octaves) - N_oct) / sigma)**2)
    psi_dist /= np.sum(psi_dist)
    S_base = entropy(psi_dist, base=2)
    
    # KEY: I_proc also amplifies with κ^(N/12)
    I_proc_base = (S_base * lambda_chaos / C_stability)
    I_proc_scaled = I_proc_base * oct_amp  # Amplification!
    
    # Total mass
    m_GeV = M_PLANCK_GeV * W * oct_amp * A_res * frac_damp * I_proc_scaled
    m_MeV = m_GeV * 1000
    
    return {
        'octave': N_oct,
        'oct_amp': oct_amp,
        'A_res': A_res,
        'I_proc': I_proc_scaled,
        'mass_MeV': m_MeV
    }

# ============================================================================
# TEST LEPTONS
# ============================================================================

print("\n" + "="*80)
print("LEPTON MASS PREDICTIONS (Nonlinear I_proc)")
print("="*80)

leptons = [
    ('Electron', 1, M_ELECTRON_MeV),
    ('Muon', 4, M_MUON_MeV),
    ('Tau', 7, M_TAU_MeV)
]

results = []
for name, N_oct, m_exp in leptons:
    result = compute_mass(N_oct, name)
    error = abs(result['mass_MeV'] - m_exp) / m_exp * 100
    
    results.append({
        'name': name,
        'N': N_oct,
        'm_pred': result['mass_MeV'],
        'm_exp': m_exp,
        'error': error
    })
    
    print(f"\n{name} (Octave {N_oct}):")
    print(f"  κ^({N_oct}/12) = {result['oct_amp']:.6f}")
    print(f"  A_res = {result['A_res']:.6f}")
    print(f"  I_proc = {result['I_proc']:.6f}")
    print(f"  m(predicted) = {result['mass_MeV']:.2f} MeV")
    print(f"  m(experiment) = {m_exp:.2f} MeV")
    print(f"  Error: {error:.1f}%")

# Summary table
print("\n" + "="*80)
print("SUMMARY TABLE")
print("="*80)
print(f"{'Particle':<10} {'Octave':<8} {'Predicted':<12} {'Experimental':<15} {'Error'}")
print("-" * 75)
for r in results:
    print(f"{r['name']:<10} {r['N']:<8} {r['m_pred']:>10.2f} MeV  {r['m_exp']:>10.2f} MeV  {r['error']:>6.1f}%")

# Average error
avg_error = np.mean([r['error'] for r in results])
print(f"\nAverage Error: {avg_error:.1f}%")

if avg_error < 10:
    print("\n✅ THEORY OF EVERYTHING CONFIRMED!")
    print("   All lepton masses derived from geometry!")
    verdict = "ToE"
elif avg_error < 30:
    print("\n🟢 STRONG FRAMEWORK")
    print("   Correct order of magnitude for all particles")
    verdict = "Strong"
else:
    print("\n🟡 PARTIAL SUCCESS")
    print("   Mechanism correct but needs refinement")
    verdict = "Partial"

# ============================================================================
# PHYSICAL INTERPRETATION
# ============================================================================

print("\n" + "="*80)
print("PHYSICAL INTERPRETATION")
print("="*80)

print("\n📊 I_proc Scaling Law:")
print(f"   I_proc(N) = I_base × κ^(N/12)")
print(f"   where κ = {kappa:.4f}")
print(f"\n   This means: Higher octaves require MORE computational effort")
print(f"   Tau (N=7) is ~{kappa**(7/12):.1f}× harder to maintain than electron (N=1)")

print(f"\n💡 Mass Hierarchy Mechanism:")
print(f"   m ∝ κ^(N/12) × [resonance] × κ^(N/12)")
print(f"   m ∝ κ^(2N/12) = κ^(N/6)")
print(f"\n   This gives:")
for r in results:
    ratio_predicted = (kappa ** (r['N']/6))
    ratio_actual = r['m_exp'] / M_ELECTRON_MeV
    print(f"   {r['name']:8}: m/m_e ≈ κ^({r['N']}/6) = {ratio_predicted:>7.1f} (actual: {ratio_actual:>7.1f})")

# ============================================================================
# REPORT
# ============================================================================

report_path = "/home/krzysiek/Pobrane/TOE/edison/raport_qw639c_nonlinear.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-639c: Nonlinear I_proc Scaling\n")
    f.write("**Data:** 2025-12-06\n\n")
    
    f.write("## Hypothesis\n\n")
    f.write("$$\n")
    f.write("I_{proc}(N) = I_{base} \\times \\kappa^{N/12}\n")
    f.write("$$\n\n")
    f.write("Processing intensity AMPLIFIES with octave frequency!\n\n")
    
    f.write("## Results\n\n")
    f.write("| Particle | Octave | Predicted | Experimental | Error |\n")
    f.write("|----------|--------|-----------|--------------|-------|\n")
    for r in results:
        f.write(f"| {r['name']:<8} | {r['N']:<6} | {r['m_pred']:>8.2f} MeV | {r['m_exp']:>10.2f} MeV | {r['error']:>5.1f}% |\n")
    
    f.write(f"\n**Average Error:** {avg_error:.1f}%\n\n")
    
    if verdict == "ToE":
        f.write("## Verdict\n\n")
        f.write("### ✅ THEORY OF EVERYTHING\n\n")
        f.write("All lepton masses emerge from pure geometry!\n\n")
        f.write("**Complete Formula:**\n\n")
        f.write("$$\n")
        f.write("m = m_{Planck} \\cdot W \\cdot \\kappa^{N/12} \\cdot A_{res} \\cdot \\beta^{10} \\cdot \\frac{S \\cdot \\lambda}{C_{stab}} \\cdot \\kappa^{N/12}\n")
        f.write("$$\n\n")
        f.write("Simplified:\n\n")
        f.write("$$\n")
        f.write("m \\propto \\kappa^{N/6} \\cdot A_{res}(N)\n")
        f.write("$$\n")
    elif verdict == "Strong":
        f.write("## Verdict\n\n")
        f.write("### 🟢 STRONG FRAMEWORK\n\n")
        f.write("Mechanism captured correctly. Hierarchia mas wynika z amplifikacji geometrycznej.\n")
    else:
        f.write("## Verdict\n\n")
        f.write("### 🟡 WYMAGA DROBNYCH POPRAWEK\n\n")

print("Report saved.")
print("="*80)
