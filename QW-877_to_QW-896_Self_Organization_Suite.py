#!/usr/bin/env python3
"""
QW-877 to QW-896: DOES K(d) SELF-ORGANIZE INTO LAYERS?
======================================================
"The only honest question: Does K(d) spontaneously create hierarchy?"

CRITICAL DIFFERENCE FROM QW-857-876:
- NO manual β^N scaling
- NO assumed layer structure
- ONLY K(d) as input
- Look for SPONTANEOUS emergence of:
  - Periodic/quasi-periodic structure
  - Eigenvalue clustering
  - Self-similar scales
  - Localization at different d ranges

Author: Krzysztof Żuchowski / Antigravity System
Date: 2025-12-09
"""

import numpy as np
from scipy.linalg import eigh
from scipy.signal import find_peaks, periodogram
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
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
    d = np.asarray(d)
    denom = 1 + BETA_TORS * np.abs(d)
    return (ALPHA_GEO * np.exp(1j * (OMEGA * d + PHI))) / denom

def K_real(d):
    return np.real(K_complex(d))

def K_magnitude(d):
    return np.abs(K_complex(d))

# ============================================================================
# QW-877: INTRINSIC PERIODICITY OF K(d)
# ============================================================================
def qw877_intrinsic_periodicity():
    """
    Does K(d) have intrinsic periodicity?
    cos(ωd) has period = 2π/ω = 8
    But what about K(d) as a whole?
    """
    print("[QW-877] Intrinsic Periodicity of K(d)")
    
    d_vals = np.linspace(0, 100, 1000)
    K_vals = np.real(K_complex(d_vals))
    
    # Find peaks (local maxima)
    peaks, properties = find_peaks(K_vals, height=0)
    
    if len(peaks) >= 2:
        peak_positions = d_vals[peaks]
        peak_spacing = np.diff(peak_positions)
        avg_period = np.mean(peak_spacing)
        period_std = np.std(peak_spacing)
        
        # Theoretical period of cos(ωd) = 2π/ω
        theoretical_period = 2 * np.pi / OMEGA
    else:
        avg_period = 0
        period_std = 0
        theoretical_period = 2 * np.pi / OMEGA
    
    # FFT analysis
    freqs, psd = periodogram(K_vals, fs=1000/100)  # fs = samples per unit d
    dominant_idx = np.argmax(psd[1:]) + 1  # Skip DC
    dominant_freq = freqs[dominant_idx]
    dominant_period = 1 / dominant_freq if dominant_freq > 0 else 0
    
    return {
        "theoretical_period": float(theoretical_period),
        "observed_avg_period": float(avg_period),
        "period_std": float(period_std),
        "fft_dominant_period": float(dominant_period),
        "N_peaks_found": len(peaks),
        "IS_PERIODIC": period_std / avg_period < 0.1 if avg_period > 0 else False,
        "NATURAL_LAYER_SIZE": float(avg_period) if avg_period > 0 else float(theoretical_period)
    }

# ============================================================================
# QW-878: ZEROS OF K(d) - NATURAL BOUNDARIES
# ============================================================================
def qw878_zeros_of_K():
    """
    Where does K(d) = 0? These are natural layer boundaries.
    """
    print("[QW-878] Zeros of K(d)")
    
    d_vals = np.linspace(0.1, 100, 10000)
    K_vals = np.real(K_complex(d_vals))
    
    # Find zero crossings
    zero_crossings = []
    for i in range(len(K_vals) - 1):
        if K_vals[i] * K_vals[i+1] < 0:
            # Linear interpolation
            d_zero = d_vals[i] - K_vals[i] * (d_vals[i+1] - d_vals[i]) / (K_vals[i+1] - K_vals[i])
            zero_crossings.append(d_zero)
    
    zero_crossings = np.array(zero_crossings)
    
    if len(zero_crossings) >= 2:
        spacing = np.diff(zero_crossings)
        avg_spacing = np.mean(spacing)
        
        # Check if spacing is regular (defines layers)
        regularity = np.std(spacing) / avg_spacing if avg_spacing > 0 else np.inf
    else:
        avg_spacing = 0
        regularity = np.inf
    
    # Theoretical: cos(ωd + φ) = 0 when ωd + φ = π/2 + nπ
    # d_n = (π/2 - φ + nπ) / ω
    theoretical_zeros = [(np.pi/2 - PHI + n*np.pi) / OMEGA for n in range(20)]
    theoretical_zeros = [z for z in theoretical_zeros if z > 0]
    
    return {
        "N_zeros_found": len(zero_crossings),
        "zero_positions": [float(z) for z in zero_crossings[:10]],
        "avg_spacing": float(avg_spacing),
        "spacing_regularity": float(regularity),
        "theoretical_zeros": [float(z) for z in theoretical_zeros[:10]],
        "REGULAR_LAYERS": regularity < 0.2,
        "LAYER_WIDTH": float(avg_spacing) if avg_spacing > 0 else 4.0
    }

# ============================================================================
# QW-879: EIGENVALUE CLUSTERING (NO MANUAL LAYERS)
# ============================================================================
def qw879_eigenvalue_clustering():
    """
    Do eigenvalues of K matrix naturally cluster?
    Use hierarchical clustering on eigenvalues.
    """
    print("[QW-879] Eigenvalue Clustering")
    
    N = 200
    d_max = 50.0  # Large enough for multiple "periods"
    
    # Build K matrix
    M = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            d = abs(i - j) * d_max / N
            M[i, j] = K_complex(d)
    H = (M + M.conj().T) / 2
    
    eigenvalues = np.linalg.eigvalsh(H)
    eigenvalues = np.sort(eigenvalues)[::-1]
    
    # Hierarchical clustering on eigenvalues
    eig_2d = eigenvalues.reshape(-1, 1)
    
    if len(eigenvalues) > 2:
        linkage_matrix = linkage(eig_2d, method='ward')
        
        # Try different numbers of clusters
        cluster_results = []
        for n_clusters in [2, 3, 4, 5, 6]:
            clusters = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
            
            # Compute within-cluster variance
            wc_var = 0
            for c in range(1, n_clusters + 1):
                cluster_eig = eigenvalues[clusters == c]
                if len(cluster_eig) > 1:
                    wc_var += np.var(cluster_eig) * len(cluster_eig)
            wc_var /= len(eigenvalues)
            
            cluster_results.append({
                "n_clusters": n_clusters,
                "within_cluster_variance": float(wc_var),
                "sizes": [int(np.sum(clusters == c)) for c in range(1, n_clusters + 1)]
            })
        
        # Find optimal (elbow method: look for biggest drop in variance)
        variances = [r["within_cluster_variance"] for r in cluster_results]
        var_drops = [variances[i] - variances[i+1] for i in range(len(variances)-1)]
        optimal_n = cluster_results[np.argmax(var_drops)]["n_clusters"] + 1
    else:
        cluster_results = []
        optimal_n = 1
    
    return {
        "N_eigenvalues": len(eigenvalues),
        "lambda_range": [float(eigenvalues[-1]), float(eigenvalues[0])],
        "cluster_analysis": cluster_results,
        "optimal_n_clusters": optimal_n,
        "NATURAL_CLUSTERS_FOUND": optimal_n >= 3
    }

# ============================================================================
# QW-880: EIGENVECTOR LOCALIZATION (SPONTANEOUS)
# ============================================================================
def qw880_eigenvector_localization():
    """
    Do eigenvectors localize at different d ranges spontaneously?
    """
    print("[QW-880] Eigenvector Localization")
    
    N = 150
    d_max = 40.0
    
    M = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            d = abs(i - j) * d_max / N
            M[i, j] = K_complex(d)
    H = (M + M.conj().T) / 2
    
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    
    # For each eigenvector, find center of mass (localization point)
    centers = []
    widths = []
    for i in range(len(eigenvalues)):
        psi = np.abs(eigenvectors[:, i])**2
        psi = psi / np.sum(psi)
        
        positions = np.arange(N) * d_max / N
        center = np.sum(positions * psi)
        width = np.sqrt(np.sum((positions - center)**2 * psi))
        
        centers.append(center)
        widths.append(width)
    
    centers = np.array(centers)
    widths = np.array(widths)
    
    # Are centers spread across d, or clustered?
    center_spread = np.std(centers)
    avg_width = np.mean(widths)
    
    # Localization ratio: if width << d_max, states are localized
    localization_ratio = avg_width / d_max
    
    # Check for clustering of centers (potential layers)
    center_hist, bin_edges = np.histogram(centers, bins=10)
    n_empty_bins = np.sum(center_hist == 0)
    
    return {
        "center_spread": float(center_spread),
        "avg_width": float(avg_width),
        "localization_ratio": float(localization_ratio),
        "n_empty_bins_in_d": n_empty_bins,
        "LOCALIZED_STATES": localization_ratio < 0.3,
        "LAYER_STRUCTURE": n_empty_bins >= 3
    }

# ============================================================================
# QW-881: SELF-SIMILARITY OF EIGENSPECTRUM
# ============================================================================
def qw881_self_similarity():
    """
    Is the eigenspectrum self-similar at different scales?
    """
    print("[QW-881] Self-Similarity of Eigenspectrum")
    
    # Compare eigenspectra at different d_max
    results = []
    for d_max in [10, 20, 40, 80]:
        N = 100
        M = np.zeros((N, N), dtype=complex)
        for i in range(N):
            for j in range(N):
                d = abs(i - j) * d_max / N
                M[i, j] = K_complex(d)
        H = (M + M.conj().T) / 2
        
        eigenvalues = np.linalg.eigvalsh(H)
        eigenvalues = np.sort(eigenvalues)[::-1][:20]  # Top 20
        
        # Normalize to [0, 1]
        if eigenvalues.max() != eigenvalues.min():
            normalized = (eigenvalues - eigenvalues.min()) / (eigenvalues.max() - eigenvalues.min())
        else:
            normalized = np.zeros_like(eigenvalues)
        
        results.append({
            "d_max": d_max,
            "normalized_spectrum": [float(x) for x in normalized]
        })
    
    # Compare spectra (correlation)
    correlations = []
    for i in range(len(results)):
        for j in range(i+1, len(results)):
            spec1 = np.array(results[i]["normalized_spectrum"])
            spec2 = np.array(results[j]["normalized_spectrum"])
            corr = np.corrcoef(spec1, spec2)[0, 1]
            correlations.append({
                "scales": (results[i]["d_max"], results[j]["d_max"]),
                "correlation": float(corr)
            })
    
    avg_corr = np.mean([c["correlation"] for c in correlations])
    
    return {
        "scale_comparisons": correlations,
        "avg_cross_scale_correlation": float(avg_corr),
        "SELF_SIMILAR": avg_corr > 0.9
    }

# ============================================================================
# QW-882: NATURAL SCALE RATIOS
# ============================================================================
def qw882_natural_scale_ratios():
    """
    What scale ratios appear naturally in K(d)?
    """
    print("[QW-882] Natural Scale Ratios")
    
    # K(d) has two characteristic scales:
    # 1. Period of cos: 2π/ω = 8
    # 2. Decay scale: 1/β = 100
    
    period = 2 * np.pi / OMEGA
    decay_scale = 1 / BETA_TORS
    
    scale_ratio = decay_scale / period
    
    # Check if this ratio appears in eigenspectrum
    N = 150
    d_max = 50.0
    M = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            d = abs(i - j) * d_max / N
            M[i, j] = K_complex(d)
    H = (M + M.conj().T) / 2
    
    eigenvalues = np.linalg.eigvalsh(H)
    pos_eig = eigenvalues[eigenvalues > 0.1]
    pos_eig = np.sort(pos_eig)[::-1]
    
    if len(pos_eig) >= 2:
        # All ratios
        ratios = []
        for i in range(min(10, len(pos_eig))):
            for j in range(i+1, min(10, len(pos_eig))):
                ratios.append(pos_eig[i] / pos_eig[j])
        
        ratios = np.array(ratios)
        
        # Check if scale_ratio appears
        near_scale_ratio = np.sum(np.abs(ratios - scale_ratio) / scale_ratio < 0.1)
        near_period_ratio = np.sum(np.abs(ratios - period) / period < 0.2)
    else:
        near_scale_ratio = 0
        near_period_ratio = 0
    
    return {
        "oscillation_period": float(period),
        "decay_scale": float(decay_scale),
        "natural_scale_ratio": float(scale_ratio),
        "eigenvalue_ratios_near_scale_ratio": int(near_scale_ratio),
        "SCALE_RATIO_IN_SPECTRUM": near_scale_ratio > 0
    }

# ============================================================================
# QW-883: SPECTRAL GAP STRUCTURE (WITHOUT LAYERS)
# ============================================================================
def qw883_spectral_gaps():
    """
    Are there natural gaps in the spectrum?
    Gaps would indicate layer boundaries.
    """
    print("[QW-883] Spectral Gap Structure")
    
    N = 200
    d_max = 60.0
    
    M = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            d = abs(i - j) * d_max / N
            M[i, j] = K_complex(d)
    H = (M + M.conj().T) / 2
    
    eigenvalues = np.linalg.eigvalsh(H)
    eigenvalues = np.sort(eigenvalues)[::-1]
    
    # Compute gaps
    gaps = -np.diff(eigenvalues)
    
    # Find significant gaps (> 2× average)
    avg_gap = np.mean(gaps)
    significant_gaps = []
    for i, gap in enumerate(gaps):
        if gap > 2 * avg_gap:
            significant_gaps.append({
                "position": i,
                "eigenvalue": float(eigenvalues[i]),
                "gap_size": float(gap),
                "relative_size": float(gap / avg_gap)
            })
    
    # Does gap structure reveal layers?
    if len(significant_gaps) >= 2:
        gap_positions = [g["position"] for g in significant_gaps]
        gap_spacing = np.diff(gap_positions)
        regular_gaps = np.std(gap_spacing) / np.mean(gap_spacing) < 0.3 if len(gap_spacing) > 0 else False
    else:
        regular_gaps = False
    
    return {
        "N_significant_gaps": len(significant_gaps),
        "significant_gaps": significant_gaps[:5],
        "avg_gap": float(avg_gap),
        "REGULAR_GAP_STRUCTURE": regular_gaps,
        "NATURAL_LAYERS_FROM_GAPS": len(significant_gaps) >= 3
    }

# ============================================================================
# QW-884: RESONANCE STRUCTURE
# ============================================================================
def qw884_resonance_structure():
    """
    Does K(d) create resonance peaks that define layers?
    """
    print("[QW-884] Resonance Structure")
    
    d_vals = np.linspace(0.1, 100, 2000)
    K_mag = K_magnitude(d_vals)
    
    # Find resonance peaks
    peaks, properties = find_peaks(K_mag, height=0.5, distance=10)
    
    if len(peaks) >= 2:
        peak_d = d_vals[peaks]
        peak_heights = K_mag[peaks]
        
        # Are peaks evenly spaced?
        spacing = np.diff(peak_d)
        regularity = np.std(spacing) / np.mean(spacing) if len(spacing) > 0 else np.inf
        
        # Peak height decay
        if len(peak_heights) >= 2:
            log_heights = np.log(peak_heights)
            decay_fit = np.polyfit(peak_d, log_heights, 1)
            decay_rate = -decay_fit[0]
        else:
            decay_rate = 0
    else:
        regularity = np.inf
        decay_rate = 0
        peak_d = []
    
    return {
        "N_resonance_peaks": len(peaks) if 'peaks' in dir() else 0,
        "peak_positions": [float(p) for p in peak_d[:10]] if len(peak_d) > 0 else [],
        "spacing_regularity": float(regularity),
        "peak_decay_rate": float(decay_rate),
        "REGULAR_RESONANCES": regularity < 0.2 if regularity < np.inf else False
    }

# ============================================================================
# QW-885: EFFECTIVE DIMENSION VS d
# ============================================================================
def qw885_effective_dimension():
    """
    Does the effective dimension change with d?
    This would indicate layer transitions.
    """
    print("[QW-885] Effective Dimension vs d")
    
    results = []
    for d_center in [5, 10, 20, 40, 80]:
        N = 50
        d_range = 5  # Window around center
        
        M = np.zeros((N, N), dtype=complex)
        for i in range(N):
            for j in range(N):
                d = d_center - d_range/2 + abs(i - j) * d_range / N
                M[i, j] = K_complex(d)
        H = (M + M.conj().T) / 2
        
        # Spectral dimension from Laplacian
        A = np.abs(H)
        D = np.diag(np.sum(A, axis=1))
        L = D - A
        
        eigenvalues = np.linalg.eigvalsh(L)
        eigenvalues = eigenvalues[eigenvalues > 1e-10]
        
        if len(eigenvalues) >= 5:
            # Heat trace approximation
            ts = np.logspace(-1, 1, 10)
            Ps = [np.sum(np.exp(-eigenvalues * t)) for t in ts]
            
            log_t = np.log(ts)
            log_P = np.log(np.array(Ps) + 1e-10)
            
            slopes = -2 * np.gradient(log_P, log_t)
            d_s = np.mean(slopes[3:7])
        else:
            d_s = 0
        
        results.append({
            "d_center": d_center,
            "d_spectral": float(d_s)
        })
    
    # Does dimension vary significantly?
    dimensions = [r["d_spectral"] for r in results]
    dim_variation = np.std(dimensions) / np.mean(dimensions) if np.mean(dimensions) > 0 else 0
    
    return {
        "dimension_vs_d": results,
        "dimension_variation": float(dim_variation),
        "DIMENSION_VARIES": dim_variation > 0.1,
        "LAYER_TRANSITIONS": dim_variation > 0.3
    }

# ============================================================================
# QW-886: SPONTANEOUS HIERARCHY MEASURE
# ============================================================================
def qw886_spontaneous_hierarchy():
    """
    Without any manual scaling, what hierarchy emerges purely from K(d)?
    """
    print("[QW-886] Spontaneous Hierarchy Measure")
    
    N = 300
    d_max = 100.0  # Large d_max to capture all structure
    
    M = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            d = abs(i - j) * d_max / N
            M[i, j] = K_complex(d)
    H = (M + M.conj().T) / 2
    
    eigenvalues = np.linalg.eigvalsh(H)
    pos_eig = eigenvalues[eigenvalues > 1e-10]
    
    if len(pos_eig) >= 2:
        lambda_max = np.max(pos_eig)
        lambda_min = np.min(pos_eig)
        hierarchy = lambda_max / lambda_min
        orders = np.log10(hierarchy)
    else:
        orders = 0
    
    # For comparison: SM needs 5.5 orders
    target = 5.5
    
    return {
        "lambda_max": float(lambda_max) if 'lambda_max' in dir() else 0,
        "lambda_min": float(lambda_min) if 'lambda_min' in dir() else 0,
        "SPONTANEOUS_HIERARCHY_ORDERS": float(orders),
        "SM_TARGET": target,
        "GAP_TO_SM": float(target - orders),
        "SUFFICIENT_FOR_SM": orders >= target
    }

# ============================================================================
# QW-887: β_TORS DERIVED FROM K(d)?
# ============================================================================
def qw887_beta_from_Kd():
    """
    Can β_tors be derived from K(d) structure rather than assumed?
    """
    print("[QW-887] β_tors Derived from K(d)?")
    
    # K(d) = α cos(ωd+φ) / (1 + βd)
    # At d = 1/β, the denominator doubles
    # This defines a natural transition scale
    
    d_transition = 1 / BETA_TORS
    
    # What is K(d) ratio at this scale?
    K_at_0 = K_magnitude(0.1)
    K_at_transition = K_magnitude(d_transition)
    K_ratio = K_at_0 / K_at_transition
    
    # From eigenspectrum, can we recover this ratio?
    N = 100
    
    # Small d spectrum
    M_small = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            d = abs(i - j) * 5.0 / N  # d in [0, 5]
            M_small[i, j] = K_complex(d)
    H_small = (M_small + M_small.conj().T) / 2
    eig_small = np.linalg.eigvalsh(H_small)
    
    # Large d spectrum
    M_large = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            d = d_transition + abs(i - j) * 5.0 / N  # d in [100, 105]
            M_large[i, j] = K_complex(d)
    H_large = (M_large + M_large.conj().T) / 2
    eig_large = np.linalg.eigvalsh(H_large)
    
    # Ratio of max eigenvalues
    eig_ratio = np.max(np.abs(eig_small)) / np.max(np.abs(eig_large))
    
    # Recovered β
    # If K(d) ~ 1/(1+βd), then ratio ~ (1 + βd_large) ~ βd_large for large d
    # So β ~ eig_ratio / d_transition
    recovered_beta = eig_ratio / d_transition / 10  # Rough estimate
    
    return {
        "assumed_beta": BETA_TORS,
        "d_transition": float(d_transition),
        "K_ratio_at_transition": float(K_ratio),
        "eigenvalue_ratio": float(eig_ratio),
        "recovered_beta_estimate": float(recovered_beta),
        "BETA_DERIVABLE": abs(recovered_beta - BETA_TORS) / BETA_TORS < 0.5
    }

# ============================================================================
# QW-888: LAYER COUNT FROM K(d)
# ============================================================================
def qw888_layer_count():
    """
    How many "natural layers" does K(d) structure suggest?
    """
    print("[QW-888] Layer Count from K(d)")
    
    # Count oscillations until amplitude drops to 1%
    d_test = np.linspace(0, 500, 5000)
    K_test = K_magnitude(d_test)
    
    threshold = K_test[0] * 0.01
    
    # Find where K drops below threshold
    d_cutoff_idx = np.where(K_test < threshold)[0]
    d_cutoff = d_test[d_cutoff_idx[0]] if len(d_cutoff_idx) > 0 else 500
    
    # Count oscillations up to cutoff
    period = 2 * np.pi / OMEGA
    n_oscillations = d_cutoff / period
    
    # Each oscillation could be a "layer"
    # Or layers could be defined by zeros
    zeros = qw878_zeros_of_K()
    n_zeros = zeros["N_zeros_found"]
    
    return {
        "d_cutoff_1percent": float(d_cutoff),
        "period": float(period),
        "n_oscillations": float(n_oscillations),
        "n_zeros_found": n_zeros,
        "NATURAL_LAYER_COUNT": int(n_oscillations),
        "LAYERS_FROM_ZEROS": n_zeros // 2  # Each zero pair bounds a layer
    }

# ============================================================================
# QW-889: HIERARCHY AMPLIFICATION MECHANISM
# ============================================================================
def qw889_hierarchy_amplification():
    """
    Is there any mechanism in K(d) that amplifies hierarchy?
    """
    print("[QW-889] Hierarchy Amplification Mechanism")
    
    # Test: does coupling between "layers" (d ranges) amplify hierarchy?
    N = 50
    n_regions = 4
    region_width = 10.0
    
    # Build coupled system (regions at different d)
    total_dim = N * n_regions
    H = np.zeros((total_dim, total_dim), dtype=complex)
    
    for region in range(n_regions):
        d_start = region * region_width
        i_start = region * N
        
        # Intra-region
        for i in range(N):
            for j in range(N):
                d = d_start + abs(i - j) * region_width / N
                H[i_start + i, i_start + j] = K_complex(d)
        
        # Inter-region coupling (use K at boundary)
        if region < n_regions - 1:
            d_boundary = (region + 1) * region_width
            coupling_strength = K_magnitude(d_boundary)
            
            j_start = (region + 1) * N
            for i in range(N):
                H[i_start + i, j_start + i] = coupling_strength
                H[j_start + i, i_start + i] = coupling_strength
    
    H = (H + H.conj().T) / 2
    eigenvalues = np.linalg.eigvalsh(H)
    pos_eig = eigenvalues[eigenvalues > 1e-10]
    
    if len(pos_eig) >= 2:
        coupled_hierarchy = np.max(pos_eig) / np.min(pos_eig)
        coupled_orders = np.log10(coupled_hierarchy)
    else:
        coupled_orders = 0
    
    # Compare to single region
    single_res = qw886_spontaneous_hierarchy()
    single_orders = single_res["SPONTANEOUS_HIERARCHY_ORDERS"]
    
    amplification = coupled_orders / single_orders if single_orders > 0 else 0
    
    return {
        "single_region_orders": float(single_orders),
        "coupled_regions_orders": float(coupled_orders),
        "amplification_factor": float(amplification),
        "NATURAL_AMPLIFICATION": amplification > 1.2
    }

# ============================================================================
# QW-890: CUMULATIVE HIERARCHY
# ============================================================================
def qw890_cumulative_hierarchy():
    """
    Does hierarchy accumulate naturally as we include more d range?
    """
    print("[QW-890] Cumulative Hierarchy")
    
    results = []
    for d_max in [10, 20, 40, 80, 160, 320]:
        N = min(200, d_max * 5)
        
        M = np.zeros((N, N), dtype=complex)
        for i in range(N):
            for j in range(N):
                d = abs(i - j) * d_max / N
                M[i, j] = K_complex(d)
        H = (M + M.conj().T) / 2
        
        eigenvalues = np.linalg.eigvalsh(H)
        pos_eig = eigenvalues[eigenvalues > 1e-10]
        
        if len(pos_eig) >= 2:
            orders = np.log10(np.max(pos_eig) / np.min(pos_eig))
        else:
            orders = 0
        
        results.append({
            "d_max": d_max,
            "orders": float(orders)
        })
    
    # Does hierarchy grow with d_max?
    d_vals = [r["d_max"] for r in results]
    order_vals = [r["orders"] for r in results]
    
    growth = (order_vals[-1] - order_vals[0]) / (d_vals[-1] - d_vals[0]) if len(results) > 1 else 0
    
    return {
        "results": results,
        "hierarchy_growth_per_d": float(growth),
        "CUMULATIVE_GROWTH": order_vals[-1] > order_vals[0] * 1.5,
        "final_orders": float(order_vals[-1])
    }

# ============================================================================
# QW-891 - QW-896: SYNTHESIS AND DIAGNOSTICS
# ============================================================================

def qw891_oscillation_contribution():
    """Does the oscillatory part (cos) contribute to hierarchy?"""
    print("[QW-891] Oscillation Contribution")
    
    N = 150
    d_max = 50.0
    
    # With oscillation
    M_osc = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            d = abs(i - j) * d_max / N
            M_osc[i, j] = K_complex(d)
    H_osc = (M_osc + M_osc.conj().T) / 2
    eig_osc = np.linalg.eigvalsh(H_osc)
    pos_osc = eig_osc[eig_osc > 1e-10]
    orders_osc = np.log10(np.max(pos_osc) / np.min(pos_osc)) if len(pos_osc) >= 2 else 0
    
    # Without oscillation (pure decay)
    M_decay = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            d = abs(i - j) * d_max / N
            M_decay[i, j] = ALPHA_GEO / (1 + BETA_TORS * d)
    H_decay = (M_decay + M_decay.conj().T) / 2
    eig_decay = np.linalg.eigvalsh(H_decay)
    pos_decay = eig_decay[eig_decay > 1e-10]
    orders_decay = np.log10(np.max(pos_decay) / np.min(pos_decay)) if len(pos_decay) >= 2 else 0
    
    return {
        "orders_with_oscillation": float(orders_osc),
        "orders_pure_decay": float(orders_decay),
        "oscillation_contribution": float(orders_osc - orders_decay),
        "OSCILLATION_HELPS": orders_osc > orders_decay
    }

def qw892_decay_contribution():
    """Does the decay part (1+βd) contribute to hierarchy?"""
    print("[QW-892] Decay Contribution")
    
    N = 150
    d_max = 50.0
    
    # With decay
    orders_with_decay = qw886_spontaneous_hierarchy()["SPONTANEOUS_HIERARCHY_ORDERS"]
    
    # Without decay (pure oscillation, manually set β=0)
    M_pure = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            d = abs(i - j) * d_max / N
            M_pure[i, j] = ALPHA_GEO * np.exp(1j * (OMEGA * d + PHI))
    H_pure = (M_pure + M_pure.conj().T) / 2
    eig_pure = np.linalg.eigvalsh(H_pure)
    pos_pure = eig_pure[eig_pure > 1e-10]
    orders_pure = np.log10(np.max(pos_pure) / np.min(pos_pure)) if len(pos_pure) >= 2 else 0
    
    return {
        "orders_with_decay": float(orders_with_decay),
        "orders_pure_oscillation": float(orders_pure),
        "decay_contribution": float(orders_with_decay - orders_pure),
        "DECAY_DOMINATES": orders_with_decay > orders_pure * 1.5
    }

def qw893_alpha_contribution():
    """Does α = 4ln(2) matter for hierarchy?"""
    print("[QW-893] Alpha Contribution")
    
    N = 100
    d_max = 40.0
    
    results = []
    for alpha in [1.0, 2.0, ALPHA_GEO, 4.0, 8.0]:
        M = np.zeros((N, N), dtype=complex)
        for i in range(N):
            for j in range(N):
                d = abs(i - j) * d_max / N
                denom = 1 + BETA_TORS * d
                M[i, j] = alpha * np.exp(1j * (OMEGA * d + PHI)) / denom
        H = (M + M.conj().T) / 2
        eig = np.linalg.eigvalsh(H)
        pos = eig[eig > 1e-10]
        orders = np.log10(np.max(pos) / np.min(pos)) if len(pos) >= 2 else 0
        
        results.append({"alpha": alpha, "orders": float(orders)})
    
    # Alpha doesn't change orders (only scales eigenvalues uniformly)
    order_values = [r["orders"] for r in results]
    variation = np.std(order_values) / np.mean(order_values) if np.mean(order_values) > 0 else 0
    
    return {
        "results": results,
        "order_variation": float(variation),
        "ALPHA_AFFECTS_HIERARCHY": variation > 0.1
    }

def qw894_omega_contribution():
    """Does ω (oscillation frequency) affect hierarchy?"""
    print("[QW-894] Omega Contribution")
    
    N = 100
    d_max = 40.0
    
    results = []
    for omega in [np.pi/8, np.pi/4, np.pi/2, np.pi]:
        M = np.zeros((N, N), dtype=complex)
        for i in range(N):
            for j in range(N):
                d = abs(i - j) * d_max / N
                denom = 1 + BETA_TORS * d
                M[i, j] = ALPHA_GEO * np.exp(1j * (omega * d + PHI)) / denom
        H = (M + M.conj().T) / 2
        eig = np.linalg.eigvalsh(H)
        pos = eig[eig > 1e-10]
        orders = np.log10(np.max(pos) / np.min(pos)) if len(pos) >= 2 else 0
        
        results.append({"omega": omega, "period": float(2*np.pi/omega), "orders": float(orders)})
    
    return {
        "results": results,
        "OMEGA_AFFECTS_HIERARCHY": np.std([r["orders"] for r in results]) > 0.2
    }

def qw895_beta_contribution():
    """Does β (decay rate) affect hierarchy?"""
    print("[QW-895] Beta Contribution")
    
    N = 100
    d_max = 40.0
    
    results = []
    for beta in [0.001, 0.01, 0.05, 0.1, 0.5]:
        M = np.zeros((N, N), dtype=complex)
        for i in range(N):
            for j in range(N):
                d = abs(i - j) * d_max / N
                denom = 1 + beta * d
                M[i, j] = ALPHA_GEO * np.exp(1j * (OMEGA * d + PHI)) / denom
        H = (M + M.conj().T) / 2
        eig = np.linalg.eigvalsh(H)
        pos = eig[eig > 1e-10]
        orders = np.log10(np.max(pos) / np.min(pos)) if len(pos) >= 2 else 0
        
        results.append({"beta": beta, "decay_scale": float(1/beta), "orders": float(orders)})
    
    return {
        "results": results,
        "BETA_AFFECTS_HIERARCHY": np.std([r["orders"] for r in results]) > 0.5
    }

def qw896_final_verdict(all_results):
    """Final synthesis: Does K(d) self-organize into hierarchy?"""
    print("[QW-896] Final Verdict")
    
    verdict = {
        "SPONTANEOUS_ORDERS": 0,
        "KEY_FINDINGS": [],
        "WHAT_WORKS": [],
        "WHAT_FAILS": [],
        "VERDICT": ""
    }
    
    # Check spontaneous hierarchy
    if "qw886" in all_results:
        verdict["SPONTANEOUS_ORDERS"] = all_results["qw886"].get("SPONTANEOUS_HIERARCHY_ORDERS", 0)
        if verdict["SPONTANEOUS_ORDERS"] >= 5.5:
            verdict["KEY_FINDINGS"].append("K(d) alone gives SM hierarchy!")
        else:
            verdict["KEY_FINDINGS"].append(f"K(d) gives {verdict['SPONTANEOUS_ORDERS']:.2f} orders (need 5.5)")
    
    # Check what creates layers
    if "qw877" in all_results and all_results["qw877"].get("IS_PERIODIC"):
        verdict["WHAT_WORKS"].append("Intrinsic periodicity detected")
    
    if "qw879" in all_results and all_results["qw879"].get("NATURAL_CLUSTERS_FOUND"):
        verdict["WHAT_WORKS"].append("Eigenvalue clusters form naturally")
    
    if "qw883" in all_results and all_results["qw883"].get("NATURAL_LAYERS_FROM_GAPS"):
        verdict["WHAT_WORKS"].append("Spectral gaps define natural layers")
    
    if "qw889" in all_results and all_results["qw889"].get("NATURAL_AMPLIFICATION"):
        verdict["WHAT_WORKS"].append("Natural amplification mechanism exists")
    
    # What fails
    if verdict["SPONTANEOUS_ORDERS"] < 5.5:
        verdict["WHAT_FAILS"].append(f"Gap of {5.5 - verdict['SPONTANEOUS_ORDERS']:.2f} orders to SM")
    
    if "qw895" in all_results and all_results["qw895"].get("BETA_AFFECTS_HIERARCHY"):
        verdict["KEY_FINDINGS"].append("β is critical for hierarchy")
    
    # Final verdict
    if len(verdict["WHAT_WORKS"]) >= 3:
        verdict["VERDICT"] = "K(d) has intrinsic self-organizing structure"
    else:
        verdict["VERDICT"] = "K(d) structure is insufficient for full hierarchy"
    
    return verdict

# ============================================================================
# MAIN RUNNER
# ============================================================================
def run_self_organization_suite():
    print("=" * 80)
    print("QW-877 TO QW-896: DOES K(d) SELF-ORGANIZE INTO LAYERS?")
    print("=" * 80)
    print("NO manual β^N scaling. ONLY K(d) structure analyzed.")
    print("=" * 80)
    
    all_results = {}
    
    tests = [
        ("qw877", qw877_intrinsic_periodicity),
        ("qw878", qw878_zeros_of_K),
        ("qw879", qw879_eigenvalue_clustering),
        ("qw880", qw880_eigenvector_localization),
        ("qw881", qw881_self_similarity),
        ("qw882", qw882_natural_scale_ratios),
        ("qw883", qw883_spectral_gaps),
        ("qw884", qw884_resonance_structure),
        ("qw885", qw885_effective_dimension),
        ("qw886", qw886_spontaneous_hierarchy),
        ("qw887", qw887_beta_from_Kd),
        ("qw888", qw888_layer_count),
        ("qw889", qw889_hierarchy_amplification),
        ("qw890", qw890_cumulative_hierarchy),
        ("qw891", qw891_oscillation_contribution),
        ("qw892", qw892_decay_contribution),
        ("qw893", qw893_alpha_contribution),
        ("qw894", qw894_omega_contribution),
        ("qw895", qw895_beta_contribution),
    ]
    
    for test_id, test_func in tests:
        try:
            result = test_func()
            all_results[test_id] = result
            
            # Print key result
            for key in ["SPONTANEOUS_HIERARCHY_ORDERS", "IS_PERIODIC", "NATURAL_CLUSTERS_FOUND", 
                       "NATURAL_LAYERS_FROM_GAPS", "NATURAL_AMPLIFICATION", "orders"]:
                if key in result:
                    print(f"  → {key}: {result[key]}")
                    break
        except Exception as e:
            print(f"  ❌ ERROR: {e}")
            all_results[test_id] = {"ERROR": str(e)}
    
    # Final verdict
    verdict = qw896_final_verdict(all_results)
    all_results["qw896"] = verdict
    print(f"\n[QW-896] Verdict: {verdict['VERDICT']}")
    print(f"  Spontaneous orders: {verdict['SPONTANEOUS_ORDERS']:.2f}")
    
    # Save
    with open("RAPORT_QW877_QW896_SELF_ORGANIZATION.md", "w") as f:
        f.write("# RAPORT: CZY K(d) SAMOORGANIZUJE SIĘ W WARSTWY? (QW-877 – QW-896)\n\n")
        f.write(f"**Data:** {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write("**Metodologia:** Zero ręcznego skalowania β^N\n\n")
        f.write("## Kluczowe Wyniki\n\n")
        f.write(f"**Spontaniczne rzędy hierarchii:** {verdict['SPONTANEOUS_ORDERS']:.2f}\n")
        f.write(f"**Cel SM:** 5.5\n\n")
        f.write("## Co Działa\n")
        for item in verdict["WHAT_WORKS"]:
            f.write(f"- ✅ {item}\n")
        f.write("\n## Co Nie Działa\n")
        for item in verdict["WHAT_FAILS"]:
            f.write(f"- ❌ {item}\n")
        f.write(f"\n## Werdykt\n**{verdict['VERDICT']}**\n")
    
    with open("RAPORT_QW877_QW896_DATA.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    
    print("\n" + "=" * 80)
    print("SELF-ORGANIZATION SUITE COMPLETE")
    print("=" * 80)
    
    return all_results

if __name__ == "__main__":
    run_self_organization_suite()
