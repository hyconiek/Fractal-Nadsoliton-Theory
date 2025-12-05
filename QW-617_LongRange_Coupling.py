#!/usr/bin/env python3
# QW-617: LONG-RANGE COUPLING EFFECTS
# Purpose: Test how β_tors (long-range coupling strength) affects emergent structure
# Question: Czy silniejsze/słabsze coupling dla oddalonych oktaw zmienia geometrię 12×20?
# Date: 2025-12-05

import numpy as np
from scipy.linalg import eigh

print("="*80)
print("QW-617: WPŁYW LONG-RANGE COUPLING NA EMERGENCJĘ")
print("="*80)
print("Test: Jak β_tors wpływa na strukturę oktaw×warstw?")
print("="*80)

# ============================================================================
# PARAMETERS
# ============================================================================
N_OCTAVES = 12
N_LAYERS = 20
ALPHA_GEO = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6

# Test different β_tors values
# Baseline: 0.01
# Weaker damping → stronger long-range: 0.001, 0.005
# Stronger damping → weaker long-range: 0.02, 0.05, 0.1
BETA_VALUES = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1]

print(f"\nNetwork: {N_OCTAVES} oktaw × {N_LAYERS} warstw")
print(f"Testing β_tors: {BETA_VALUES}")
print(f"Baseline (frozen): 0.01")
print("-" * 40)

# ============================================================================
# TEST FOR EACH β_tors
# ============================================================================

results = []

for beta_tors in BETA_VALUES:
    print(f"\n{'='*60}")
    print(f"β_tors = {beta_tors}")
    print('='*60)
    
    # Build coupling kernel
    def K(d):
        return ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + beta_tors * d)
    
    # 1) OCTAVE STRUCTURE across layers
    eigenvalue_sets = []
    
    for layer_idx in [0, 5, 10, 15, 19]:
        beta_eff = beta_tors * (layer_idx + 1)
        
        S_layer = np.zeros((N_OCTAVES, N_OCTAVES))
        for i in range(N_OCTAVES):
            for j in range(N_OCTAVES):
                d = abs(i - j)
                S_layer[i, j] = ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + beta_eff * d)
        
        eigs, _ = eigh(S_layer)
        eigs_norm = eigs / np.max(np.abs(eigs))
        eigenvalue_sets.append(eigs_norm)
    
    # Measure orthogonality (oktawy⊥warstwy)
    corr = np.corrcoef(eigenvalue_sets[0], eigenvalue_sets[-1])[0, 1]
    
    # 2) LAYER COUPLING STRENGTH decay
    coupling_strengths = []
    for layer in range(N_LAYERS):
        beta_eff = beta_tors * (layer + 1)
        S = np.zeros((N_OCTAVES, N_OCTAVES))
        for i in range(N_OCTAVES):
            for j in range(N_OCTAVES):
                d = abs(i - j)
                S[i, j] = ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + beta_eff * d)
        coupling = np.linalg.norm(S, 'fro') / N_OCTAVES
        coupling_strengths.append(coupling)
    
    # Exponential fit
    log_couplings = np.log(coupling_strengths)
    layers_array = np.arange(N_LAYERS)
    coeffs = np.polyfit(layers_array, log_couplings, 1)
    lambda_decay = -coeffs[0]
    
    # R²
    fit_vals = np.polyval(coeffs, layers_array)
    ss_tot = np.sum((log_couplings - np.mean(log_couplings))**2)
    ss_res = np.sum((log_couplings - fit_vals)**2)
    r2_exp = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    # 3) LONG-RANGE STRENGTH (coupling at max distance)
    S_full = np.zeros((N_OCTAVES, N_OCTAVES))
    for i in range(N_OCTAVES):
        for j in range(N_OCTAVES):
            S_full[i, j] = K(abs(i - j))
    
    # Ratio: long-range (d=11) vs short-range (d=1)
    K_short = K(1)
    K_long = K(11)  # Max distance in 12 octaves
    long_range_ratio = abs(K_long / K_short) if K_short != 0 else 0
    
    results.append({
        'beta': beta_tors,
        'orthogonality': corr,
        'lambda_decay': lambda_decay,
        'r2_exp': r2_exp,
        'long_range_ratio': long_range_ratio,
        'K_short': K_short,
        'K_long': K_long
    })
    
    print(f"  Orthogonality (oktawy⊥warstwy): r = {corr:.4f}")
    print(f"  Layer decay: λ = {lambda_decay:.4f}")
    print(f"  Long-range ratio K(11)/K(1): {long_range_ratio:.4f}")

# ============================================================================
# ANALYSIS
# ============================================================================

print("\n" + "="*80)
print("ANALIZA WPŁYWU β_tors")
print("="*80)

print("\n| β_tors | Ortho (r) | λ_decay | K_long/K_short | Status |")
print("|--------|-----------|---------|----------------|--------|")

# Baseline
baseline = next((r for r in results if abs(r['beta'] - 0.01) < 1e-6), results[0])

for res in results:
    # Compare to baseline
    ortho_ok = res['orthogonality'] > 0.9
    exp_ok = res['r2_exp'] > 0.9
    
    # Deviation from baseline
    ortho_change = abs(res['orthogonality'] - baseline['orthogonality'])
    
    if ortho_ok and ortho_change < 0.05:
        status = "✅ Stable"
    elif ortho_ok:
        status = "🟡 Changed"
    else:
        status = "❌ Broken"
    
    is_baseline = abs(res['beta'] - 0.01) < 1e-6
    marker = " *" if is_baseline else ""
    
    print(f"| {res['beta']:6.3f}{marker} | {res['orthogonality']:9.4f} | {res['lambda_decay']:7.4f} | {res['long_range_ratio']:14.4f} | {status:10s} |")

# ============================================================================
# CRITICAL INSIGHTS
# ============================================================================

print("\n" + "="*80)
print("KLUCZOWE WNIOSKI")
print("="*80)

# How orthogonality changes with beta
ortho_values = [r['orthogonality'] for r in results]
beta_values_array = np.array([r['beta'] for r in results])

ortho_std = np.std(ortho_values)
ortho_range = max(ortho_values) - min(ortho_values)

print(f"\nOrthogonality variation:")
print(f"  Range: {min(ortho_values):.4f} → {max(ortho_values):.4f}")
print(f"  Std: {ortho_std:.4f}")
print(f"  Relative change: {ortho_range/np.mean(ortho_values) * 100:.1f}%")

# Long-range coupling effect
lr_ratios = [r['long_range_ratio'] for r in results]
print(f"\nLong-range coupling K(11)/K(1):")
print(f"  Range: {min(lr_ratios):.4f} → {max(lr_ratios):.4f}")
print(f"  β=0.001 (strong LR): {results[0]['long_range_ratio']:.4f}")
print(f"  β=0.01 (baseline): {baseline['long_range_ratio']:.4f}")
print(f"  β=0.1 (weak LR): {results[-1]['long_range_ratio']:.4f}")

if ortho_std < 0.01:
    print("\n✅ STRUKTURA NIEZALEŻNA OD β_tors!")
    print("   Orthogonalność oktaw⊥warstw jest ROBUST")
    print("   Long-range coupling NIE zmienia fundamentalnej geometrii")
    verdict = "independent"
elif ortho_range < 0.05:
    print("\n🟡 SŁABA ZALEŻNOŚĆ")
    print(f"   Orthogonalność zmienia się o {ortho_range*100:.1f}%")
    verdict = "weak_dep"
else:
    print("\n⚠️  SILNA ZALEŻNOŚĆ")
    print(f"   β_tors znacząco wpływa na strukturę ({ortho_range*100:.1f}% change)")
    verdict = "strong_dep"

# ============================================================================
# PHYSICAL INTERPRETATION
# ============================================================================

print("\n" + "="*80)
print("INTERPRETACJA FIZYCZNA")
print("="*80)

if verdict == "independent":
    print("\n✅ FUNDAMENTALNA GEOMETRIA POTWIERDZONA!")
    print("\n**Odkrycie:**")
    print("  Struktura 12×20 (oktawy⊥warstwy) jest NIEZALEŻNA")
    print("  od siły long-range coupling!")
    print("\n**Znaczenie:**")
    print("  - β_tors = fenomenologiczny parametr skali")
    print("  - Orthogonalność = FUNDAMENTALNA własność")
    print("  - Geometria wynika z TOPOLOGII sieci, nie z coupling strength")
    print("\n**Implikacja:**")
    print("  12 oktaw i 20 warstw to NIEZMIENNIK topologiczny,")
    print("  nawet jeśli long-range coupling się zmienia!")
elif verdict == "weak_dep":
    print("\n🟡 CZĘŚCIOWA NIEZALEŻNOŚĆ")
    print("β_tors wpływa słabo - struktura w większości robust")
else:
    print("\n⚠️  SILNA ZALEŻNOŚĆ")
    print("Long-range coupling kształtuje emergentną strukturę")

# Test if stronger LR → better dimensional amplification?
print("\n" + "="*80)
print("BONUS: Long-Range vs Dimensional Amplification")
print("="*80)

print("\nHipoteza: Silniejszy LR coupling → wyższy wymiar emergentny?")
print("(Stronger distant interactions → more dimensions emerge)")

# Approximate: stronger LR (lower beta) → flatter coupling → more uniform network
print(f"\n| β_tors | K_long/K_short | Prediction |")
print(f"|--------|----------------|------------|")
for res in results[:3]:  # Show a few
    if res['long_range_ratio'] > 0.5:
        pred = "Higher d_eff"
    elif res['long_range_ratio'] > 0.1:
        pred = "Balanced"
    else:
        pred = "Lower d_eff"
    print(f"| {res['beta']:.3f}  | {res['long_range_ratio']:14.4f} | {pred:11s} |")

print("\n(Requires QW-618 to confirm dimensional effects)")

# ============================================================================
# REPORT
# ============================================================================

report_path = "raport_qw617_longrange.md"
print(f"\nGenerating report: {report_path}...")

with open(report_path, "w") as f:
    f.write("# Raport QW-617: Wpływ Long-Range Coupling na Emergencję\n")
    f.write("**Data:** 2025-12-05\n")
    f.write("**Test:** Jak β_tors wpływa na strukturę 12×20?\n\n")
    
    f.write("## 1. Wyniki\n")
    f.write("| β_tors | Orthogonality | Status |\n|--------|---------------|--------|\n")
    for res in results:
        status = "✅" if res['orthogonality'] > 0.9 else "🟡"
        f.write(f"| {res['beta']:.3f} | {res['orthogonality']:.4f} | {status} |\n")
    f.write("\n")
    
    f.write(f"**Orthogonality range:** {ortho_range:.4f} ({ortho_range/np.mean(ortho_values)*100:.1f}%)\n\n")
    
    if verdict == "independent":
        f.write("### ✅ STRUKTURA NIEZALEŻNA!\n")
        f.write("12×20 lattice geometria jest ROBUST do zmian β_tors.\n")
        f.write("Fundamentalna topologia, nie fenomenologia!\n")
    elif verdict == "weak_dep":
        f.write("### 🟡 SŁABA ZALEŻNOŚĆ\n")
    else:
        f.write("### ⚠️ SILNA ZALEŻNOŚĆ\n")

print("Report saved.")
print("="*80)
