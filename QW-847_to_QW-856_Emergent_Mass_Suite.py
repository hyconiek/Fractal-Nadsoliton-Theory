#!/usr/bin/env python3
"""
QW-847 to QW-856: TRUE EMERGENT MASS HIERARCHY
===============================================
"What emerges from K(d)?" - NOT "Does it match Standard Model?"

PHILOSOPHY:
- ONLY INPUT: K(d) = α·cos(ωd+φ)/(1+βd) with frozen parameters
- ZERO calibration to known masses
- ZERO external constants
- Analyze INTERNAL ratios and structure

Author: Krzysztof Żuchowski / Antigravity System
Date: 2025-12-09
"""

import numpy as np
from scipy.linalg import eigh
from scipy.cluster.hierarchy import linkage, fcluster
import time
import json

# ============================================================================
# FROZEN KERNEL PARAMETERS - THE ONLY INPUT
# ============================================================================
ALPHA_GEO = 4 * np.log(2)  # 2.7726
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01

def K_complex(d):
    """The Unified Kernel - FROZEN. This is the ONLY foundation."""
    d = np.asarray(d)
    denom = 1 + BETA_TORS * np.abs(d)
    return (ALPHA_GEO * np.exp(1j * (OMEGA * d + PHI))) / denom

def K_real(d):
    return np.real(K_complex(d))

def K_magnitude(d):
    return np.abs(K_complex(d))

# ============================================================================
# BUILD KERNEL MATRIX
# ============================================================================
def build_K_matrix(N, d_max):
    """Build NxN matrix where M[i,j] = K(|i-j| * d_max/N)"""
    M = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            d = abs(i - j) * d_max / N
            M[i, j] = K_complex(d)
    # Make Hermitian
    H = (M + M.conj().T) / 2
    return H

# ============================================================================
# QW-847: EIGENSPECTRUM OF K(d)
# ============================================================================
def qw847_eigenspectrum_K():
    """
    What is the raw eigenspectrum of K(d)?
    No comparison to physics - just the numbers.
    """
    print("[QW-847] Raw Eigenspectrum of K(d)")
    
    N = 200
    d_max = 20.0
    H = build_K_matrix(N, d_max)
    
    eigenvalues, eigenvectors = eigh(H)
    eigenvalues = np.sort(eigenvalues)[::-1]  # Descending
    
    return {
        "N": N,
        "d_max": d_max,
        "N_eigenvalues": len(eigenvalues),
        "lambda_max": float(eigenvalues[0]),
        "lambda_min": float(eigenvalues[-1]),
        "lambda_range": float(eigenvalues[0] - eigenvalues[-1]),
        "top_10": [float(e) for e in eigenvalues[:10]],
        "bottom_10": [float(e) for e in eigenvalues[-10:]],
        "eigenvalues": eigenvalues
    }

# ============================================================================
# QW-848: INTERNAL RATIOS
# ============================================================================
def qw848_internal_ratios(eigenvalues):
    """
    What are the internal ratios λᵢ/λⱼ?
    Do they cluster around specific values?
    """
    print("[QW-848] Internal Ratios")
    
    # Use only positive eigenvalues for ratios
    pos = eigenvalues[eigenvalues > 0.1]
    
    if len(pos) < 5:
        return {"ERROR": "Insufficient positive eigenvalues"}
    
    # Compute all ratios
    ratios = []
    for i in range(min(20, len(pos))):
        for j in range(i+1, min(20, len(pos))):
            if pos[j] > 0.01:
                ratio = pos[i] / pos[j]
                ratios.append(ratio)
    
    ratios = np.array(ratios)
    
    # Statistics
    unique_ratios = np.unique(np.round(ratios, 2))
    
    # Check for powers of 2 or 4
    log2_ratios = np.log2(ratios[ratios > 0])
    log4_ratios = np.log(ratios[ratios > 0]) / np.log(4)
    
    # Are any close to integers?
    near_int_log2 = [r for r in log2_ratios if abs(r - round(r)) < 0.1]
    near_int_log4 = [r for r in log4_ratios if abs(r - round(r)) < 0.1]
    
    return {
        "N_ratios": len(ratios),
        "ratio_min": float(ratios.min()),
        "ratio_max": float(ratios.max()),
        "ratio_mean": float(ratios.mean()),
        "ratio_std": float(ratios.std()),
        "N_unique": len(unique_ratios),
        "near_power_of_2": len(near_int_log2),
        "near_power_of_4": len(near_int_log4),
        "sample_ratios": [float(r) for r in ratios[:10]]
    }

# ============================================================================
# QW-849: GAP STRUCTURE
# ============================================================================
def qw849_gap_structure(eigenvalues):
    """
    What is the gap structure Δλᵢ = λᵢ₊₁ - λᵢ?
    Is it regular? Does it reveal quantization?
    """
    print("[QW-849] Gap Structure")
    
    sorted_eig = np.sort(eigenvalues)[::-1]  # Descending
    gaps = -np.diff(sorted_eig)  # Positive gaps
    
    # Remove tiny gaps (numerical noise)
    significant_gaps = gaps[gaps > 0.01]
    
    if len(significant_gaps) < 5:
        return {"ERROR": "Insufficient significant gaps"}
    
    # Statistics
    gap_mean = np.mean(significant_gaps)
    gap_std = np.std(significant_gaps)
    gap_cv = gap_std / gap_mean  # Coefficient of variation
    
    # Is it regular? (CV < 0.3 = regular)
    is_regular = gap_cv < 0.3
    
    # Check for quantization: are gaps multiples of a base value?
    gap_min = significant_gaps.min()
    gap_ratios = significant_gaps / gap_min
    near_integer = np.sum(np.abs(gap_ratios - np.round(gap_ratios)) < 0.1)
    
    return {
        "N_gaps": len(significant_gaps),
        "gap_mean": float(gap_mean),
        "gap_std": float(gap_std),
        "gap_cv": float(gap_cv),
        "is_regular": bool(is_regular),
        "gap_min": float(gap_min),
        "gap_max": float(significant_gaps.max()),
        "quantization_fraction": float(near_integer / len(significant_gaps)),
        "sample_gaps": [float(g) for g in significant_gaps[:10]]
    }

# ============================================================================
# QW-850: RATIO HISTOGRAM
# ============================================================================
def qw850_ratio_histogram(eigenvalues):
    """
    Histogram of log(λᵢ/λⱼ) - do peaks appear at log(4^n)?
    """
    print("[QW-850] Ratio Histogram")
    
    pos = eigenvalues[eigenvalues > 0.1]
    
    if len(pos) < 5:
        return {"ERROR": "Insufficient data"}
    
    # All ratios
    log_ratios = []
    for i in range(min(30, len(pos))):
        for j in range(i+1, min(30, len(pos))):
            if pos[j] > 0.01:
                log_ratios.append(np.log(pos[i] / pos[j]))
    
    log_ratios = np.array(log_ratios)
    
    # Histogram
    hist, bin_edges = np.histogram(log_ratios, bins=20)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Find peaks
    peak_indices = []
    for i in range(1, len(hist)-1):
        if hist[i] > hist[i-1] and hist[i] > hist[i+1]:
            peak_indices.append(i)
    
    peak_values = [bin_centers[i] for i in peak_indices]
    
    # Check if peaks are at ln(4), 2*ln(4), etc.
    ln4 = np.log(4)
    near_ln4_multiples = []
    for p in peak_values:
        n = p / ln4
        if abs(n - round(n)) < 0.2:
            near_ln4_multiples.append((p, round(n)))
    
    return {
        "N_log_ratios": len(log_ratios),
        "log_ratio_mean": float(log_ratios.mean()),
        "log_ratio_std": float(log_ratios.std()),
        "N_peaks": len(peak_indices),
        "peak_positions": peak_values,
        "peaks_at_ln4_multiples": near_ln4_multiples,
        "ln4": float(ln4),
        "NATURAL_BASE_DETECTED": len(near_ln4_multiples) > 0
    }

# ============================================================================
# QW-851: RESONANCE SEARCH
# ============================================================================
def qw851_resonance_search(eigenvalues):
    """
    Search for: log(λ) = a + b*n (exponential quantization)
    """
    print("[QW-851] Resonance Search")
    
    pos = eigenvalues[eigenvalues > 0.1]
    pos = np.sort(pos)[::-1]  # Descending
    
    if len(pos) < 10:
        return {"ERROR": "Insufficient data"}
    
    # Fit: log(λ_n) = a - b*n
    n = np.arange(len(pos[:20]))
    log_lambda = np.log(pos[:20])
    
    # Linear regression
    coeffs = np.polyfit(n, log_lambda, 1)
    b = -coeffs[0]  # Decay rate
    a = coeffs[1]   # Intercept
    
    # Residuals
    fitted = a - b * n
    residuals = log_lambda - fitted
    r_squared = 1 - np.var(residuals) / np.var(log_lambda)
    
    # What is the natural base?
    natural_base = np.exp(b)
    
    # Compare to known bases
    bases = {
        "2": 2,
        "e": np.e,
        "4": 4,
        "π": np.pi,
        "4ln2": 4*np.log(2)
    }
    
    closest_base = min(bases.items(), key=lambda x: abs(x[1] - natural_base))
    
    return {
        "decay_rate_b": float(b),
        "intercept_a": float(a),
        "r_squared": float(r_squared),
        "natural_base": float(natural_base),
        "closest_known_base": closest_base[0],
        "closest_base_value": float(closest_base[1]),
        "base_error_percent": float(abs(natural_base - closest_base[1]) / closest_base[1] * 100),
        "EXPONENTIAL_QUANTIZATION": r_squared > 0.9
    }

# ============================================================================
# QW-852: DEGENERACY COUNT
# ============================================================================
def qw852_degeneracy_count(eigenvalues):
    """
    How many eigenvalues are nearly degenerate?
    Could indicate "generations" or "families".
    """
    print("[QW-852] Degeneracy Count")
    
    sorted_eig = np.sort(eigenvalues)[::-1]
    
    # Count quasi-degenerate pairs (within 1%)
    threshold = 0.01
    degenerate_groups = []
    current_group = [sorted_eig[0]]
    
    for i in range(1, len(sorted_eig)):
        if abs(sorted_eig[i] - current_group[-1]) / abs(current_group[-1] + 1e-10) < threshold:
            current_group.append(sorted_eig[i])
        else:
            if len(current_group) > 1:
                degenerate_groups.append(len(current_group))
            current_group = [sorted_eig[i]]
    
    if len(current_group) > 1:
        degenerate_groups.append(len(current_group))
    
    # Count isolated eigenvalues
    n_isolated = len(sorted_eig) - sum(degenerate_groups)
    
    return {
        "N_degenerate_groups": len(degenerate_groups),
        "group_sizes": degenerate_groups[:10],
        "N_isolated": n_isolated,
        "degeneracy_fraction": float(sum(degenerate_groups) / len(sorted_eig)),
        "FAMILIES_DETECTED": len(set(degenerate_groups)) <= 3
    }

# ============================================================================
# QW-853: LOCALIZATION (IPR)
# ============================================================================
def qw853_localization(N=200, d_max=20.0):
    """
    Inverse Participation Ratio IPR = Σ|ψ|⁴
    High IPR = localized (particle-like)
    Low IPR = delocalized (wave-like)
    """
    print("[QW-853] Localization (IPR)")
    
    H = build_K_matrix(N, d_max)
    eigenvalues, eigenvectors = eigh(H)
    
    # IPR for each eigenstate
    iprs = []
    for i in range(len(eigenvalues)):
        psi = eigenvectors[:, i]
        psi_norm = psi / np.linalg.norm(psi)
        ipr = np.sum(np.abs(psi_norm)**4)
        iprs.append(ipr)
    
    iprs = np.array(iprs)
    
    # Sort by eigenvalue
    idx_sorted = np.argsort(eigenvalues)[::-1]
    iprs_sorted = iprs[idx_sorted]
    eigenvalues_sorted = eigenvalues[idx_sorted]
    
    # Find most localized states
    most_localized_idx = np.argmax(iprs)
    
    # Threshold for "particle-like" (IPR > 3/N)
    particle_threshold = 3.0 / N
    n_particle_like = np.sum(iprs > particle_threshold)
    
    return {
        "ipr_mean": float(iprs.mean()),
        "ipr_max": float(iprs.max()),
        "ipr_min": float(iprs.min()),
        "most_localized_eigenvalue": float(eigenvalues[most_localized_idx]),
        "particle_threshold": float(particle_threshold),
        "N_particle_like": int(n_particle_like),
        "STABLE_PARTICLES_FOUND": n_particle_like > 0,
        "top_5_localized_eigenvalues": [float(eigenvalues[i]) for i in np.argsort(iprs)[::-1][:5]]
    }

# ============================================================================
# QW-854: STABILITY VS λ
# ============================================================================
def qw854_stability_vs_lambda(N=100, d_max=15.0):
    """
    Which eigenvalues survive under time evolution?
    """
    print("[QW-854] Stability vs λ")
    
    H = build_K_matrix(N, d_max)
    eigenvalues, eigenvectors = eigh(H)
    
    # Start with superposition of all eigenstates
    psi_0 = np.sum(eigenvectors, axis=1)
    psi_0 = psi_0 / np.linalg.norm(psi_0)
    
    dt = 0.01
    n_steps = 500
    
    # Track projection onto each eigenstate
    projections = []
    for i in range(len(eigenvalues)):
        proj = abs(np.vdot(eigenvectors[:, i], psi_0))**2
        projections.append(proj)
    
    projections = np.array(projections)
    
    # Evolve and check stability
    psi = psi_0.copy()
    for step in range(n_steps):
        # Unitary evolution under H
        psi = psi - 1j * dt * (H @ psi)
        psi = psi / np.linalg.norm(psi)
    
    # Final projections
    final_projections = []
    for i in range(len(eigenvalues)):
        proj = abs(np.vdot(eigenvectors[:, i], psi))**2
        final_projections.append(proj)
    
    final_projections = np.array(final_projections)
    
    # Which survived?
    survived = np.where(final_projections > 0.01)[0]
    surviving_eigenvalues = eigenvalues[survived]
    
    return {
        "N_initial_projections": int(np.sum(projections > 0.01)),
        "N_surviving": len(survived),
        "surviving_eigenvalues": [float(e) for e in surviving_eigenvalues[:10]],
        "dominant_eigenvalue": float(eigenvalues[np.argmax(final_projections)]),
        "STABLE_MODES_EXIST": len(survived) > 0 and len(survived) < len(eigenvalues)
    }

# ============================================================================
# QW-855: HIERARCHY DEPTH
# ============================================================================
def qw855_hierarchy_depth(eigenvalues):
    """
    How many orders of magnitude does the spectrum span?
    """
    print("[QW-855] Hierarchy Depth")
    
    pos = eigenvalues[eigenvalues > 1e-10]
    
    if len(pos) < 2:
        return {"ERROR": "Insufficient data"}
    
    lambda_max = np.max(np.abs(pos))
    lambda_min = np.min(np.abs(pos))
    
    hierarchy_ratio = lambda_max / lambda_min
    orders_of_magnitude = np.log10(hierarchy_ratio)
    
    # Compare to SM hierarchy (top/electron ~ 3.4×10^5)
    sm_hierarchy = 3.4e5
    
    return {
        "lambda_max": float(lambda_max),
        "lambda_min": float(lambda_min),
        "hierarchy_ratio": float(hierarchy_ratio),
        "orders_of_magnitude": float(orders_of_magnitude),
        "SM_comparison": float(np.log10(sm_hierarchy)),
        "SUFFICIENT_HIERARCHY": orders_of_magnitude > 2
    }

# ============================================================================
# QW-856: EMERGENCE REPORT
# ============================================================================
def qw856_emergence_report(all_results):
    """
    Synthesize: What EMERGED from K(d)?
    """
    print("[QW-856] Emergence Report")
    
    report = {
        "KERNEL_PARAMETERS": {
            "alpha": ALPHA_GEO,
            "omega": OMEGA,
            "phi": PHI,
            "beta": BETA_TORS
        },
        "EMERGENT_PROPERTIES": {}
    }
    
    # Check each finding
    if "qw850" in all_results and all_results["qw850"].get("NATURAL_BASE_DETECTED"):
        report["EMERGENT_PROPERTIES"]["Natural_Base"] = all_results["qw850"]["peaks_at_ln4_multiples"]
    
    if "qw851" in all_results and all_results["qw851"].get("EXPONENTIAL_QUANTIZATION"):
        report["EMERGENT_PROPERTIES"]["Quantization_Base"] = all_results["qw851"]["natural_base"]
        report["EMERGENT_PROPERTIES"]["Closest_Known_Base"] = all_results["qw851"]["closest_known_base"]
    
    if "qw852" in all_results and all_results["qw852"].get("FAMILIES_DETECTED"):
        report["EMERGENT_PROPERTIES"]["N_Families"] = len(set(all_results["qw852"]["group_sizes"]))
    
    if "qw853" in all_results and all_results["qw853"].get("STABLE_PARTICLES_FOUND"):
        report["EMERGENT_PROPERTIES"]["N_Particle_Like"] = all_results["qw853"]["N_particle_like"]
    
    if "qw855" in all_results and all_results["qw855"].get("SUFFICIENT_HIERARCHY"):
        report["EMERGENT_PROPERTIES"]["Hierarchy_Orders"] = all_results["qw855"]["orders_of_magnitude"]
    
    # Summary
    n_emergent = len(report["EMERGENT_PROPERTIES"])
    report["SUMMARY"] = f"{n_emergent}/5 key properties emerged from K(d)"
    report["SUCCESS"] = n_emergent >= 3
    
    return report

# ============================================================================
# MAIN RUNNER
# ============================================================================
def run_emergence_suite():
    print("=" * 80)
    print("QW-847 TO QW-856: TRUE EMERGENT MASS HIERARCHY")
    print("=" * 80)
    print(f"ONLY INPUT: K(d) = {ALPHA_GEO:.4f} * cos({OMEGA:.4f}*d + {PHI:.4f}) / (1 + {BETA_TORS}*d)")
    print("ZERO external constants. ZERO calibration.")
    print("=" * 80)
    
    all_results = {}
    
    # QW-847
    result_847 = qw847_eigenspectrum_K()
    eigenvalues = result_847.pop("eigenvalues")
    all_results["qw847"] = result_847
    print(f"  λ_max = {result_847['lambda_max']:.4f}, λ_min = {result_847['lambda_min']:.4f}")
    
    # QW-848
    result_848 = qw848_internal_ratios(eigenvalues)
    all_results["qw848"] = result_848
    print(f"  {result_848['N_ratios']} ratios, {result_848['near_power_of_4']} near power-of-4")
    
    # QW-849
    result_849 = qw849_gap_structure(eigenvalues)
    all_results["qw849"] = result_849
    print(f"  Gap CV = {result_849['gap_cv']:.3f}, Regular: {result_849['is_regular']}")
    
    # QW-850
    result_850 = qw850_ratio_histogram(eigenvalues)
    all_results["qw850"] = result_850
    print(f"  {result_850['N_peaks']} peaks, Base-4 detected: {result_850['NATURAL_BASE_DETECTED']}")
    
    # QW-851
    result_851 = qw851_resonance_search(eigenvalues)
    all_results["qw851"] = result_851
    print(f"  Natural base = {result_851['natural_base']:.4f} (closest: {result_851['closest_known_base']})")
    print(f"  Exponential quantization: {result_851['EXPONENTIAL_QUANTIZATION']}")
    
    # QW-852
    result_852 = qw852_degeneracy_count(eigenvalues)
    all_results["qw852"] = result_852
    print(f"  {result_852['N_degenerate_groups']} degenerate groups, Families: {result_852['FAMILIES_DETECTED']}")
    
    # QW-853
    result_853 = qw853_localization()
    all_results["qw853"] = result_853
    print(f"  {result_853['N_particle_like']} particle-like modes")
    
    # QW-854
    result_854 = qw854_stability_vs_lambda()
    all_results["qw854"] = result_854
    print(f"  {result_854['N_surviving']} stable modes")
    
    # QW-855
    result_855 = qw855_hierarchy_depth(eigenvalues)
    all_results["qw855"] = result_855
    print(f"  Hierarchy: {result_855['orders_of_magnitude']:.2f} orders of magnitude")
    
    # QW-856
    report = qw856_emergence_report(all_results)
    all_results["qw856"] = report
    print(f"\n{report['SUMMARY']}")
    
    # Save
    with open("RAPORT_QW847_QW856_EMERGENT_MASS.md", "w") as f:
        f.write("# RAPORT: PRAWDZIWIE EMERGENTNA HIERARCHIA MAS (QW-847 – QW-856)\n\n")
        f.write(f"**Data:** {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write("**Metodologia:** Zero zewnętrznych stałych, zero kalibracji\n\n")
        
        f.write("## Parametry Jądra (JEDYNE WEJŚCIE)\n")
        f.write(f"- α = {ALPHA_GEO:.6f}\n")
        f.write(f"- ω = {OMEGA:.6f}\n")
        f.write(f"- φ = {PHI:.6f}\n")
        f.write(f"- β = {BETA_TORS:.6f}\n\n")
        
        f.write("## Co Wyłoniło Się z K(d)?\n\n")
        for prop, value in report["EMERGENT_PROPERTIES"].items():
            f.write(f"- **{prop}:** {value}\n")
        
        f.write(f"\n## Podsumowanie\n{report['SUMMARY']}\n")
    
    with open("RAPORT_QW847_QW856_DATA.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    
    print("\n" + "=" * 80)
    print("EMERGENCE SUITE COMPLETE")
    print("=" * 80)
    
    return all_results

if __name__ == "__main__":
    run_emergence_suite()
