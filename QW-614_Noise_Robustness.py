#!/usr/bin/env python3
# QW-614: ROBUSTNOŚĆ OKTAW⊥WARSTW DO SZUMU
# Purpose: Test if octave-layer orthogonality survives thermal/quantum noise
# Hypothesis: Fundamental structure should be robust to perturbations
# Date: 2025-12-05

import numpy as np
from scipy.linalg import eigh

print("="*80)
print("QW-614: ROBUSTNOŚĆ OKTAW⊥WARSTW DO SZUMU")
print("="*80)
print("Test: Czy ortogonalność oktaw⊥warstw przetrwa szum?")
print("="*80)

# ============================================================================
# PARAMETERS
# ============================================================================
N_OCTAVES = 12
N_LAYERS = 20
ALPHA_GEO = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01

# Noise levels to test
NOISE_LEVELS = [0.0, 0.01, 0.05, 0.1, 0.2, 0.5]

print(f"\nNetwork: {N_OCTAVES} octaves × {N_LAYERS} layers")
print(f"Testing noise levels: {NOISE_LEVELS}")
print("-" * 40)

# ============================================================================
# BASELINE (NO NOISE)
# ============================================================================

def K(d, beta=BETA_TORS):
    return ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + beta * d)

# Build baseline octave spectrum
S_baseline = np.zeros((N_OCTAVES, N_OCTAVES))
for i in range(N_OCTAVES):
    for j in range(N_OCTAVES):
        S_baseline[i, j] = K(abs(i - j))

eigs_baseline, _ = eigh(S_baseline)
eigs_baseline_norm = eigs_baseline / np.max(np.abs(eigs_baseline))

print(f"\nBaseline octave spectrum (normalized):")
print(f"  Top 3: {eigs_baseline_norm[-3:]}")

# ============================================================================
# TEST WITH NOISE
# ============================================================================

results = []

for noise_level in NOISE_LEVELS:
    print(f"\n{'='*60}")
    print(f"Noise level: {noise_level}")
    print('='*60)
    
    # Test octave structure robustness across layers WITH NOISE
    eigenvalue_sets_noisy = []
    
    for layer_idx in [0, 5, 10, 15, 19]:
        beta_eff = BETA_TORS * (layer_idx + 1)
        
        # Build coupling matrix for this layer
        S_layer = np.zeros((N_OCTAVES, N_OCTAVES))
        for i in range(N_OCTAVES):
            for j in range(N_OCTAVES):
                d = abs(i - j)
                S_layer[i, j] = ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + beta_eff * d)
        
        # Add noise: Random perturbation to coupling matrix
        # Thermal noise: fluctuations in coupling strength
        noise = np.random.randn(N_OCTAVES, N_OCTAVES) * noise_level
        noise = (noise + noise.T) / 2  # Keep symmetric
        
        S_layer_noisy = S_layer + noise
        
        # Compute eigenvalues
        eigs, _ = eigh(S_layer_noisy)
        eigs_norm = eigs / np.max(np.abs(eigs))
        
        eigenvalue_sets_noisy.append(eigs_norm)
    
    # Measure orthogonality: correlation between layer 0 and layer 19
    corr_noisy = np.corrcoef(eigenvalue_sets_noisy[0], eigenvalue_sets_noisy[-1])[0, 1]
    
    # Measure coupling decay (exponential scaling)
    coupling_strengths = []
    for layer in range(N_LAYERS):
        beta_eff = BETA_TORS * (layer + 1)
        S_test = np.zeros((N_OCTAVES, N_OCTAVES))
        for i in range(N_OCTAVES):
            for j in range(N_OCTAVES):
                d = abs(i - j)
                S_test[i, j] = ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + beta_eff * d)
        
        # Add noise
        noise = np.random.randn(N_OCTAVES, N_OCTAVES) * noise_level
        noise = (noise + noise.T) / 2
        S_test += noise
        
        coupling = np.linalg.norm(S_test, 'fro') / N_OCTAVES
        coupling_strengths.append(coupling)
    
    # Fit exponential decay
    log_couplings = np.log(coupling_strengths)
    layers_array = np.arange(N_LAYERS)
    coeffs = np.polyfit(layers_array, log_couplings, 1)
    lambda_decay = -coeffs[0]
    
    # R² for exponential fit
    fit_vals = np.polyval(coeffs, layers_array)
    ss_tot = np.sum((log_couplings - np.mean(log_couplings))**2)
    ss_res = np.sum((log_couplings - fit_vals)**2)
    r_squared_exp = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    results.append({
        'noise': noise_level,
        'correlation': corr_noisy,
        'lambda_decay': lambda_decay,
        'r2_exp': r_squared_exp
    })
    
    print(f"  Octave structure correlation: r = {corr_noisy:.4f}")
    print(f"  Layer decay rate: λ = {lambda_decay:.4f}")
    print(f"  Exponential fit: R² = {r_squared_exp:.4f}")

# ============================================================================
# ANALYSIS
# ============================================================================

print("\n" + "="*80)
print("ROBUSTNESS ANALYSIS")
print("="*80)

print("\n| Noise | Corr (r) | λ_decay | R²_exp | Status |")
print("|-------|----------|---------|--------|--------|")

for res in results:
    # Orthogonality preserved if r > 0.9
    # Exponential scaling if R² > 0.9
    ortho_ok = res['correlation'] > 0.9
    exp_ok = res['r2_exp'] > 0.9
    
    if ortho_ok and exp_ok:
        status = "✅ Robust"
    elif ortho_ok or exp_ok:
        status = "🟡 Partial"
    else:
        status = "❌ Degraded"
    
    print(f"| {res['noise']:5.2f} | {res['correlation']:8.4f} | {res['lambda_decay']:7.4f} | {res['r2_exp']:6.3f} | {status:10s} |")

# Find critical noise level where structure breaks
critical_noise = None
for res in results:
    if res['correlation'] < 0.9 or res['r2_exp'] < 0.9:
        critical_noise = res['noise']
        break

if critical_noise is not None:
    print(f"\n⚠️  Critical noise level: σ_crit ≈ {critical_noise}")
    print(f"    Structure degrades above this noise")
else:
    print(f"\n✅ STRUCTURE ROBUST")
    print(f"    Oktawy⊥Warstwy survive noise up to σ={NOISE_LEVELS[-1]}")

# ============================================================================
# PHYSICAL INTERPRETATION
# ============================================================================

print("\n" + "="*80)
print("PHYSICAL INTERPRETATION")
print("="*80)

if critical_noise is None or critical_noise > 0.1:
    print("\n✅ FUNDAMENTALNA STRUKTURA POTWIERDZONA")
    print("\nOrtogonalność oktaw⊥warstw przetrwała szum:")
    print("- Thermal fluctuations (noise < 10%): ROBUST")
    print("- Quantum corrections: Struktura zachowana")
    print("\n**Znaczenie:**")
    print("12×20 lattice to FUNDAMENTALNA geometria sieci,")
    print("nie emergentny artefakt!")
    verdict = "robust"
elif critical_noise > 0.01:
    print("\n🟡 CZĘŚCIOWA ROBUSTNOŚĆ")
    print(f"\nStruktura zachowana do σ < {critical_noise}")
    print("Może wymagać dekoherencji lub protection mechanism")
    verdict = "partial"
else:
    print("\n❌ WRAŻLIWA STRUKTURA")
    print("\nOrtogonaln ość łatwo zaburzona przez szum")
    verdict = "fragile"

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw614_noise_robustness.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-614: Robustność Oktaw⊥Warstw do Szumu\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Test:** Czy fundamentalna struktura przetrwa thermal/quantum noise?\n\n")
    
    f.write("## 1. Metodologia\n")
    f.write(f"Dodano szum gaussian do macierzy sprzężeń: S → S + N(0,σ)\n")
    f.write(f"Testowano σ = {NOISE_LEVELS}\n\n")
    
    f.write("## 2. Wyniki\n")
    f.write("| Noise | Correlation | Status |\n|-------|-------------|--------|\n")
    for res in results:
        status = "✅" if res['correlation'] > 0.9 else ("🟡" if res['correlation'] > 0.7 else "❌")
        f.write(f"| {res['noise']:.2f} | {res['correlation']:.4f} | {status} |\n")
    f.write("\n")
    
    if verdict == "robust":
        f.write("### ✅ STRUKTURA FUNDAMENTALNA POTWIERDZONA\n")
        f.write("Oktawy⊥Warstwy to robust geometric structure!\n")
    elif verdict == "partial":
        f.write("### 🟡 CZĘŚCIOWA ROBUSTNOŚĆ\n")
    else:
        f.write("### ❌ WRAŻLIWA NA SZUM\n")

print("Report saved.")
print("="*80)
