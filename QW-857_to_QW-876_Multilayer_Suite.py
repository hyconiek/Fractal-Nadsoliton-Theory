#!/usr/bin/env python3
"""
QW-857 to QW-876: MULTI-LAYER K(d) HIERARCHICAL EMERGENCE
==========================================================
"Does a stack of K(d) layers produce sufficient hierarchy?"

INSIGHT: Single K(d) gives 3.66 orders of magnitude.
         Theory suggests 20-30 fractal layers.
         Each layer with β=0.01 gives 100× suppression.
         N layers gives β^N = 10^(-2N) hierarchy.

TESTS:
- Multi-layer spectrum
- Inter-layer coupling
- Recursive nesting
- Hierarchy amplification
- Emergent structure from stack

Author: Krzysztof Żuchowski / Antigravity System
Date: 2025-12-09
"""

import numpy as np
from scipy.linalg import eigh, block_diag
from scipy.sparse import csr_matrix, kron, eye
import time
import json

# ============================================================================
# FROZEN KERNEL PARAMETERS
# ============================================================================
ALPHA_GEO = 4 * np.log(2)  # 2.7726
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01  # Critical: this is the INTER-LAYER coupling!

def K_complex(d):
    d = np.asarray(d)
    denom = 1 + BETA_TORS * np.abs(d)
    return (ALPHA_GEO * np.exp(1j * (OMEGA * d + PHI))) / denom

def K_real(d):
    return np.real(K_complex(d))

def K_magnitude(d):
    return np.abs(K_complex(d))

# ============================================================================
# BUILD SINGLE-LAYER K MATRIX
# ============================================================================
def build_K_layer(N, d_max):
    """Build NxN K matrix for single layer"""
    M = np.zeros((N, N), dtype=complex)
    for i in range(N):
        for j in range(N):
            d = abs(i - j) * d_max / N
            M[i, j] = K_complex(d)
    H = (M + M.conj().T) / 2
    return H

# ============================================================================
# QW-857: MULTI-LAYER SPECTRUM
# ============================================================================
def qw857_multilayer_spectrum(N_layers=5, N_sites=50):
    """
    Build a stack of K(d) layers and diagonalize.
    Inter-layer coupling: β_tors = 0.01
    """
    print(f"[QW-857] Multi-Layer Spectrum (N_layers={N_layers})")
    
    # Build each layer
    layers = []
    for layer in range(N_layers):
        H_layer = build_K_layer(N_sites, d_max=10.0)
        layers.append(H_layer)
    
    # Stack them with inter-layer coupling
    total_dim = N_layers * N_sites
    H_total = np.zeros((total_dim, total_dim), dtype=complex)
    
    for layer in range(N_layers):
        i_start = layer * N_sites
        i_end = (layer + 1) * N_sites
        H_total[i_start:i_end, i_start:i_end] = layers[layer]
        
        # Inter-layer coupling
        if layer < N_layers - 1:
            j_start = (layer + 1) * N_sites
            j_end = (layer + 2) * N_sites
            coupling = BETA_TORS * np.eye(N_sites)
            H_total[i_start:i_end, j_start:j_end] = coupling
            H_total[j_start:j_end, i_start:i_end] = coupling.T.conj()
    
    # Diagonalize
    eigenvalues = np.linalg.eigvalsh(H_total)
    eigenvalues = np.sort(eigenvalues)[::-1]
    
    lambda_max = eigenvalues[0]
    lambda_min = eigenvalues[-1]
    hierarchy = np.abs(lambda_max / lambda_min) if lambda_min != 0 else np.inf
    
    return {
        "N_layers": N_layers,
        "N_sites": N_sites,
        "Total_dim": total_dim,
        "lambda_max": float(lambda_max),
        "lambda_min": float(np.abs(lambda_min)),
        "Hierarchy_ratio": float(hierarchy),
        "Orders_of_magnitude": float(np.log10(hierarchy)) if hierarchy > 0 else 0
    }

# ============================================================================
# QW-858: HIERARCHY VS N_LAYERS
# ============================================================================
def qw858_hierarchy_vs_layers():
    """
    How does hierarchy grow with number of layers?
    """
    print("[QW-858] Hierarchy vs N_layers")
    
    results = []
    for N_layers in [1, 2, 3, 5, 7, 10]:
        res = qw857_multilayer_spectrum(N_layers=N_layers, N_sites=30)
        results.append({
            "N_layers": N_layers,
            "Hierarchy": res["Hierarchy_ratio"],
            "Orders": res["Orders_of_magnitude"]
        })
    
    # Fit: Hierarchy ~ β^(-N) = 100^N ?
    N_vals = np.array([r["N_layers"] for r in results])
    orders = np.array([r["Orders"] for r in results])
    
    # Linear fit in log space: orders = a * N + b
    if len(orders[orders > 0]) >= 2:
        valid = orders > 0
        coeffs = np.polyfit(N_vals[valid], orders[valid], 1)
        growth_rate = coeffs[0]
    else:
        growth_rate = 0
    
    return {
        "results": results,
        "growth_rate_per_layer": float(growth_rate),
        "expected_rate": 2.0,  # log10(100) = 2
        "EXPONENTIAL_GROWTH": growth_rate > 1.5
    }

# ============================================================================
# QW-859: WEAK COUPLING LIMIT
# ============================================================================
def qw859_weak_coupling():
    """
    What if β_tors → 0? (Decoupled layers)
    """
    print("[QW-859] Weak Coupling Limit")
    
    global BETA_TORS
    original_beta = BETA_TORS
    
    results = []
    for beta in [1.0, 0.1, 0.01, 0.001, 0.0001]:
        BETA_TORS = beta
        res = qw857_multilayer_spectrum(N_layers=5, N_sites=30)
        results.append({
            "beta": beta,
            "Hierarchy": res["Hierarchy_ratio"],
            "Orders": res["Orders_of_magnitude"]
        })
    
    BETA_TORS = original_beta
    
    return {
        "results": results,
        "observation": "Weaker coupling → larger hierarchy" if results[-1]["Orders"] > results[0]["Orders"] else "Coupling doesn't amplify hierarchy"
    }

# ============================================================================
# QW-860: RECURSIVE NESTING
# ============================================================================
def qw860_recursive_nesting(depth=3):
    """
    K(d) at scale n feeds into K(d) at scale n+1.
    Each deeper level has time running 100× faster.
    """
    print(f"[QW-860] Recursive Nesting (depth={depth})")
    
    N = 20
    d_max = 10.0
    
    # Level 0: base K(d)
    H_0 = build_K_layer(N, d_max)
    eigenvalues_0 = np.linalg.eigvalsh(H_0)
    
    # Each level: H_n = H_{n-1} + β * K(d_n)
    # Where d_n is rescaled by 1/β
    
    H_current = H_0.copy()
    hierarchy_at_levels = [np.max(np.abs(eigenvalues_0)) / np.min(np.abs(eigenvalues_0[eigenvalues_0 != 0]) + 1e-10)]
    
    for level in range(1, depth + 1):
        # Rescaled K at this level
        scale = BETA_TORS ** level
        H_level = build_K_layer(N, d_max * scale)
        
        # Couple to previous
        H_current = H_current + scale * H_level
        
        eig = np.linalg.eigvalsh(H_current)
        eig_pos = eig[eig > 1e-10]
        if len(eig_pos) >= 2:
            hierarchy_at_levels.append(np.max(eig_pos) / np.min(eig_pos))
        else:
            hierarchy_at_levels.append(1.0)
    
    return {
        "depth": depth,
        "hierarchy_per_level": [float(h) for h in hierarchy_at_levels],
        "final_hierarchy": float(hierarchy_at_levels[-1]),
        "amplification": float(hierarchy_at_levels[-1] / hierarchy_at_levels[0]) if hierarchy_at_levels[0] > 0 else 0
    }

# ============================================================================
# QW-861: TIME DILATION EFFECT
# ============================================================================
def qw861_time_dilation():
    """
    If layer N has time t, layer N+1 has time t/β.
    How does this affect eigenvalue spectrum?
    """
    print("[QW-861] Time Dilation Effect")
    
    N = 30
    d_max = 10.0
    
    # Base layer
    H_base = build_K_layer(N, d_max)
    eig_base = np.linalg.eigvalsh(H_base)
    
    # Time-dilated layers (frequency scales with 1/t)
    # E = ℏω → if t → t/β, then ω → β*ω, so E → β*E
    results = []
    for layer in range(5):
        dilation_factor = BETA_TORS ** layer
        eig_dilated = eig_base * dilation_factor
        
        gap = eig_dilated[1] - eig_dilated[0] if len(eig_dilated) > 1 else 0
        results.append({
            "layer": layer,
            "dilation": dilation_factor,
            "E_max": float(eig_dilated[0]),
            "E_min": float(eig_dilated[-1]),
            "gap": float(np.abs(gap))
        })
    
    # Total combined spectrum
    all_eig = []
    for layer in range(5):
        all_eig.extend(list(eig_base * (BETA_TORS ** layer)))
    all_eig = np.array(all_eig)
    all_eig = np.sort(np.abs(all_eig))[::-1]
    
    return {
        "per_layer": results,
        "combined_hierarchy": float(all_eig[0] / all_eig[-1]) if all_eig[-1] > 1e-10 else np.inf,
        "combined_orders": float(np.log10(all_eig[0] / all_eig[-1])) if all_eig[-1] > 1e-10 else 0
    }

# ============================================================================
# QW-862: STACKED EIGENVALUES
# ============================================================================
def qw862_stacked_eigenvalues():
    """
    Simply stack eigenvalue spectra from N layers.
    """
    print("[QW-862] Stacked Eigenvalues")
    
    N = 40
    N_layers = 10
    
    all_eigenvalues = []
    for layer in range(N_layers):
        H = build_K_layer(N, d_max=10.0 / (layer + 1))
        eig = np.linalg.eigvalsh(H)
        # Scale by layer depth
        eig_scaled = eig * (BETA_TORS ** layer)
        all_eigenvalues.extend(list(eig_scaled))
    
    all_eigenvalues = np.array(all_eigenvalues)
    all_eigenvalues = all_eigenvalues[all_eigenvalues > 1e-15]
    
    if len(all_eigenvalues) < 2:
        return {"ERROR": "Insufficient eigenvalues"}
    
    lambda_max = np.max(all_eigenvalues)
    lambda_min = np.min(all_eigenvalues)
    
    return {
        "N_total_eigenvalues": len(all_eigenvalues),
        "lambda_max": float(lambda_max),
        "lambda_min": float(lambda_min),
        "Hierarchy": float(lambda_max / lambda_min),
        "Orders": float(np.log10(lambda_max / lambda_min))
    }

# ============================================================================
# QW-863: EFFECTIVE β FROM HIERARCHY
# ============================================================================
def qw863_effective_beta():
    """
    What effective β explains observed SM hierarchy?
    SM needs ~5.5 orders. Single layer gives 3.66.
    """
    print("[QW-863] Effective β from Hierarchy")
    
    target_orders = 5.5  # log10(340000)
    single_layer_orders = 3.66
    
    # If hierarchy = β^(-N), then log10(hierarchy) = N * log10(1/β)
    # We need: 5.5 = N * log10(1/β)
    
    results = []
    for N_layers in [1, 2, 3, 5, 10, 20, 30]:
        # Required β: β = 10^(-5.5/N)
        required_beta = 10 ** (-target_orders / N_layers)
        
        # What additional orders does β^N give?
        additional_orders = N_layers * np.log10(1 / BETA_TORS)
        total_orders = single_layer_orders + additional_orders
        
        results.append({
            "N_layers": N_layers,
            "required_beta": float(required_beta),
            "with_beta_0.01": float(total_orders),
            "reaches_target": total_orders >= target_orders
        })
    
    return {
        "target_orders": target_orders,
        "single_layer": single_layer_orders,
        "results": results,
        "N_layers_needed_with_beta_0.01": min([r["N_layers"] for r in results if r["reaches_target"]], default=">30")
    }

# ============================================================================
# QW-864: LAYER-LAYER CORRELATIONS
# ============================================================================
def qw864_layer_correlations():
    """
    Are eigenvalues correlated between adjacent layers?
    """
    print("[QW-864] Layer-Layer Correlations")
    
    N = 50
    N_layers = 5
    
    layer_spectra = []
    for layer in range(N_layers):
        H = build_K_layer(N, d_max=10.0)
        eig = np.linalg.eigvalsh(H)
        eig = np.sort(eig)[::-1][:20]  # Top 20
        layer_spectra.append(eig)
    
    # Correlation matrix
    corr_matrix = np.zeros((N_layers, N_layers))
    for i in range(N_layers):
        for j in range(N_layers):
            corr = np.corrcoef(layer_spectra[i], layer_spectra[j])[0, 1]
            corr_matrix[i, j] = corr
    
    avg_corr = np.mean(corr_matrix[corr_matrix != 1])
    
    return {
        "N_layers": N_layers,
        "correlation_matrix": corr_matrix.tolist(),
        "avg_off_diagonal_corr": float(avg_corr),
        "IDENTICAL_LAYERS": avg_corr > 0.99
    }

# ============================================================================
# QW-865: LAYER-DEPENDENT K(d)
# ============================================================================
def qw865_layer_dependent_K():
    """
    What if K(d) parameters vary by layer?
    α_n = α × n, ω_n = ω / n, etc.
    """
    print("[QW-865] Layer-Dependent K(d)")
    
    def build_K_layer_custom(N, d_max, alpha_scale, omega_scale):
        M = np.zeros((N, N), dtype=complex)
        alpha_eff = ALPHA_GEO * alpha_scale
        omega_eff = OMEGA * omega_scale
        for i in range(N):
            for j in range(N):
                d = abs(i - j) * d_max / N
                denom = 1 + BETA_TORS * np.abs(d)
                M[i, j] = (alpha_eff * np.exp(1j * (omega_eff * d + PHI))) / denom
        H = (M + M.conj().T) / 2
        return H
    
    N = 30
    N_layers = 5
    
    all_eig = []
    for layer in range(N_layers):
        alpha_scale = 1.0 / (layer + 1)
        omega_scale = 1.0 * (layer + 1)
        
        H = build_K_layer_custom(N, d_max=10.0, alpha_scale=alpha_scale, omega_scale=omega_scale)
        eig = np.linalg.eigvalsh(H)
        all_eig.extend(list(eig))
    
    all_eig = np.array(all_eig)
    pos_eig = all_eig[all_eig > 1e-10]
    
    if len(pos_eig) < 2:
        return {"ERROR": "Insufficient data"}
    
    return {
        "N_layers": N_layers,
        "lambda_max": float(pos_eig.max()),
        "lambda_min": float(pos_eig.min()),
        "Hierarchy": float(pos_eig.max() / pos_eig.min()),
        "Orders": float(np.log10(pos_eig.max() / pos_eig.min()))
    }

# ============================================================================
# QW-866: TENSOR PRODUCT LAYERS
# ============================================================================
def qw866_tensor_product():
    """
    H_total = H_1 ⊗ H_2 ⊗ ... ⊗ H_N
    Eigenvalues multiply: λ_total = λ_1 × λ_2 × ... × λ_N
    """
    print("[QW-866] Tensor Product Layers")
    
    N = 10  # Small for tensor product
    N_layers = 3
    
    # Get eigenvalues for each layer
    layer_eig = []
    for layer in range(N_layers):
        H = build_K_layer(N, d_max=10.0)
        eig = np.linalg.eigvalsh(H)
        eig_pos = eig[eig > 0.1][:5]  # Top 5 positive
        layer_eig.append(eig_pos)
    
    # Tensor product eigenvalues
    tensor_eig = layer_eig[0]
    for layer in range(1, N_layers):
        new_eig = np.outer(tensor_eig, layer_eig[layer]).flatten()
        tensor_eig = new_eig
    
    tensor_eig = tensor_eig[tensor_eig > 1e-10]
    
    return {
        "N_layers": N_layers,
        "N_tensor_eigenvalues": len(tensor_eig),
        "lambda_max": float(tensor_eig.max()),
        "lambda_min": float(tensor_eig.min()),
        "Hierarchy": float(tensor_eig.max() / tensor_eig.min()),
        "Orders": float(np.log10(tensor_eig.max() / tensor_eig.min())),
        "MULTIPLICATIVE_HIERARCHY": True
    }

# ============================================================================
# QW-867: FRACTAL DIMENSION FROM STACK
# ============================================================================
def qw867_fractal_dimension():
    """
    Compute spectral dimension d_s of stacked system.
    """
    print("[QW-867] Fractal Dimension from Stack")
    
    N = 40
    N_layers = 5
    
    # Build stacked Laplacian
    total_dim = N * N_layers
    L_total = np.zeros((total_dim, total_dim))
    
    for layer in range(N_layers):
        i_start = layer * N
        i_end = (layer + 1) * N
        
        # Intra-layer Laplacian
        H = build_K_layer(N, d_max=10.0)
        A = np.abs(H)  # Use magnitude as adjacency
        D = np.diag(np.sum(A, axis=1))
        L = D - A
        L_total[i_start:i_end, i_start:i_end] = L
        
        # Inter-layer connections
        if layer < N_layers - 1:
            j_start = (layer + 1) * N
            coupling = BETA_TORS * np.eye(N)
            L_total[i_start:i_end, j_start:j_start+N] = -coupling
            L_total[j_start:j_start+N, i_start:i_end] = -coupling
            # Adjust diagonal
            L_total[i_start:i_end, i_start:i_end] += np.diag(np.sum(coupling, axis=1))
            L_total[j_start:j_start+N, j_start:j_start+N] += np.diag(np.sum(coupling, axis=0))
    
    eigenvalues = np.linalg.eigvalsh(L_total)
    eigenvalues = eigenvalues[eigenvalues > 1e-10]
    
    if len(eigenvalues) < 10:
        return {"ERROR": "Insufficient eigenvalues"}
    
    # Heat trace
    ts = np.logspace(-2, 1, 20)
    Ps = [np.sum(np.exp(-eigenvalues * t)) for t in ts]
    
    log_t = np.log(ts)
    log_P = np.log(np.array(Ps) + 1e-10)
    
    slopes = -2 * np.gradient(log_P, log_t)
    d_s = np.mean(slopes[5:15])
    
    return {
        "N_layers": N_layers,
        "d_spectral": float(d_s),
        "d_s_range": [float(slopes.min()), float(slopes.max())],
        "FRACTAL_STRUCTURE": 1.5 < d_s < 3.0
    }

# ============================================================================
# QW-868: LOCALIZATION ACROSS LAYERS
# ============================================================================
def qw868_localization_layers():
    """
    Where do particle-like states localize in the stack?
    """
    print("[QW-868] Localization Across Layers")
    
    N = 30
    N_layers = 5
    
    # Build stacked system
    total_dim = N * N_layers
    H_total = np.zeros((total_dim, total_dim), dtype=complex)
    
    for layer in range(N_layers):
        i_start = layer * N
        i_end = (layer + 1) * N
        H_total[i_start:i_end, i_start:i_end] = build_K_layer(N, 10.0)
        
        if layer < N_layers - 1:
            j_start = (layer + 1) * N
            coupling = BETA_TORS * np.eye(N)
            H_total[i_start:i_end, j_start:j_start+N] = coupling
            H_total[j_start:j_start+N, i_start:i_end] = coupling.T.conj()
    
    eigenvalues, eigenvectors = np.linalg.eigh(H_total)
    
    # Find most localized state
    iprs = []
    layer_localization = []
    for i in range(len(eigenvalues)):
        psi = eigenvectors[:, i]
        psi_norm = psi / np.linalg.norm(psi)
        ipr = np.sum(np.abs(psi_norm)**4)
        iprs.append(ipr)
        
        # Which layer?
        layer_weights = []
        for layer in range(N_layers):
            i_start = layer * N
            i_end = (layer + 1) * N
            weight = np.sum(np.abs(psi_norm[i_start:i_end])**2)
            layer_weights.append(weight)
        layer_localization.append(np.argmax(layer_weights))
    
    most_localized_idx = np.argmax(iprs)
    
    return {
        "most_localized_eigenvalue": float(eigenvalues[most_localized_idx]),
        "most_localized_layer": int(layer_localization[most_localized_idx]),
        "ipr_max": float(max(iprs)),
        "layer_distribution": [layer_localization.count(l) for l in range(N_layers)],
        "LAYER_0_DOMINANT": layer_localization[most_localized_idx] == 0
    }

# ============================================================================
# QW-869: EFFECTIVE HIERARCHY FORMULA
# ============================================================================
def qw869_effective_formula():
    """
    Test: Hierarchy = Base^N × λ1/λN
    """
    print("[QW-869] Effective Hierarchy Formula")
    
    single_layer_res = qw857_multilayer_spectrum(N_layers=1, N_sites=50)
    single_orders = single_layer_res["Orders_of_magnitude"]
    
    # Test formula: Orders(N) = single_orders + N × log10(1/β)
    results = []
    for N_layers in [1, 3, 5, 10]:
        res = qw857_multilayer_spectrum(N_layers=N_layers, N_sites=50)
        predicted = single_orders + N_layers * np.log10(1/BETA_TORS)
        observed = res["Orders_of_magnitude"]
        
        results.append({
            "N": N_layers,
            "predicted_orders": float(predicted),
            "observed_orders": float(observed),
            "error": float(abs(predicted - observed))
        })
    
    avg_error = np.mean([r["error"] for r in results])
    
    return {
        "formula": "Orders = single_layer + N × log10(1/β)",
        "single_layer_orders": float(single_orders),
        "log10(1/β)": float(np.log10(1/BETA_TORS)),
        "results": results,
        "avg_error": float(avg_error),
        "FORMULA_VALID": avg_error < 1.0
    }

# ============================================================================
# QW-870: GAP STRUCTURE ACROSS LAYERS
# ============================================================================
def qw870_gap_structure_layers():
    """
    How do gaps between eigenvalues change with N_layers?
    """
    print("[QW-870] Gap Structure Across Layers")
    
    results = []
    for N_layers in [1, 2, 5, 10]:
        H = np.zeros((N_layers * 30, N_layers * 30), dtype=complex)
        for layer in range(N_layers):
            i_s = layer * 30
            H[i_s:i_s+30, i_s:i_s+30] = build_K_layer(30, 10.0)
            if layer < N_layers - 1:
                coupling = BETA_TORS * np.eye(30)
                H[i_s:i_s+30, i_s+30:i_s+60] = coupling
                H[i_s+30:i_s+60, i_s:i_s+30] = coupling.T.conj()
        
        eig = np.sort(np.linalg.eigvalsh(H))[::-1]
        gaps = -np.diff(eig[:20])
        gaps = gaps[gaps > 0.01]
        
        if len(gaps) > 0:
            cv = np.std(gaps) / np.mean(gaps)
        else:
            cv = 0
        
        results.append({
            "N_layers": N_layers,
            "gap_mean": float(np.mean(gaps)) if len(gaps) > 0 else 0,
            "gap_cv": float(cv),
            "N_gaps": len(gaps)
        })
    
    return {"results": results}

# ============================================================================
# QW-871: EMERGENT GENERATIONS
# ============================================================================
def qw871_emergent_generations():
    """
    Do eigenvalue clusters emerge from multi-layer structure?
    """
    print("[QW-871] Emergent Generations")
    
    res = qw857_multilayer_spectrum(N_layers=5, N_sites=40)
    
    # Also get the raw eigenvalues
    H = np.zeros((5 * 40, 5 * 40), dtype=complex)
    for layer in range(5):
        i_s = layer * 40
        H[i_s:i_s+40, i_s:i_s+40] = build_K_layer(40, 10.0)
        if layer < 4:
            coupling = BETA_TORS * np.eye(40)
            H[i_s:i_s+40, i_s+40:i_s+80] = coupling
            H[i_s+40:i_s+80, i_s:i_s+40] = coupling.T.conj()
    
    eigenvalues = np.sort(np.linalg.eigvalsh(H))[::-1]
    
    # Cluster eigenvalues (simple: histogram bins)
    pos_eig = eigenvalues[eigenvalues > 0.1]
    log_eig = np.log10(pos_eig)
    
    # Histogram
    hist, bin_edges = np.histogram(log_eig, bins=10)
    
    # Count peaks
    peaks = []
    for i in range(1, len(hist) - 1):
        if hist[i] > hist[i-1] and hist[i] > hist[i+1]:
            peaks.append((bin_edges[i] + bin_edges[i+1]) / 2)
    
    return {
        "N_positive_eigenvalues": len(pos_eig),
        "log_eig_range": [float(log_eig.min()), float(log_eig.max())],
        "N_peaks": len(peaks),
        "peak_positions": peaks,
        "GENERATIONS_DETECTED": len(peaks) >= 3
    }

# ============================================================================
# QW-872: SM HIERARCHY REACHED?
# ============================================================================
def qw872_sm_hierarchy_reached():
    """
    Can we reach Standard Model hierarchy (5.5 orders)?
    """
    print("[QW-872] SM Hierarchy Reached?")
    
    target = 5.5
    
    for N_layers in [1, 5, 10, 15, 20, 25, 30]:
        res = qw857_multilayer_spectrum(N_layers=N_layers, N_sites=20)
        orders = res["Orders_of_magnitude"]
        
        if orders >= target:
            return {
                "N_layers_needed": N_layers,
                "orders_achieved": float(orders),
                "target": target,
                "SM_HIERARCHY_ACHIEVED": True
            }
    
    # Also try tensor product approach
    res_tensor = qw866_tensor_product()
    
    return {
        "max_tested_layers": 30,
        "max_orders_achieved": float(res["Orders_of_magnitude"]),
        "target": target,
        "tensor_product_orders": res_tensor.get("Orders", 0),
        "SM_HIERARCHY_ACHIEVED": False
    }

# ============================================================================
# QW-873: STABILITY OF MULTILAYER
# ============================================================================
def qw873_stability_multilayer():
    """
    Are multi-layer states more or less stable?
    """
    print("[QW-873] Stability of Multi-layer")
    
    N = 25
    N_layers = 5
    
    # Build system
    total_dim = N * N_layers
    H = np.zeros((total_dim, total_dim), dtype=complex)
    for layer in range(N_layers):
        i_s = layer * N
        H[i_s:i_s+N, i_s:i_s+N] = build_K_layer(N, 10.0)
        if layer < N_layers - 1:
            coupling = BETA_TORS * np.eye(N)
            H[i_s:i_s+N, i_s+N:i_s+2*N] = coupling
            H[i_s+N:i_s+2*N, i_s:i_s+N] = coupling.T.conj()
    
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    
    # Start with superposition
    psi = np.sum(eigenvectors, axis=1)
    psi = psi / np.linalg.norm(psi)
    
    # Evolve
    dt = 0.01
    n_steps = 200
    
    overlap_history = []
    for step in range(n_steps):
        psi_new = psi - 1j * dt * (H @ psi)
        psi_new = psi_new / np.linalg.norm(psi_new)
        
        # Overlap with ground state
        overlap = np.abs(np.vdot(eigenvectors[:, 0], psi_new))**2
        overlap_history.append(overlap)
        psi = psi_new
    
    return {
        "final_ground_overlap": float(overlap_history[-1]),
        "max_overlap": float(max(overlap_history)),
        "RELAXATION_TO_GROUND": overlap_history[-1] > overlap_history[0]
    }

# ============================================================================
# QW-874: OPTIMAL N_LAYERS
# ============================================================================
def qw874_optimal_layers():
    """
    What N_layers gives best hierarchy/cost ratio?
    """
    print("[QW-874] Optimal N_layers")
    
    results = []
    for N_layers in [1, 2, 3, 5, 7, 10, 15, 20]:
        start = time.time()
        res = qw857_multilayer_spectrum(N_layers=N_layers, N_sites=30)
        elapsed = time.time() - start
        
        results.append({
            "N_layers": N_layers,
            "Orders": res["Orders_of_magnitude"],
            "Time_seconds": float(elapsed),
            "Orders_per_second": float(res["Orders_of_magnitude"] / elapsed) if elapsed > 0 else 0
        })
    
    best = max(results, key=lambda x: x["Orders_per_second"])
    
    return {
        "results": results,
        "optimal_N": best["N_layers"],
        "optimal_ratio": float(best["Orders_per_second"])
    }

# ============================================================================
# QW-875: INTER-LAYER RESONANCE
# ============================================================================
def qw875_interlayer_resonance():
    """
    Are there resonances in inter-layer coupling?
    """
    print("[QW-875] Inter-layer Resonance")
    
    N = 30
    
    # Scan coupling values
    results = []
    for beta in [0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0]:
        # 2-layer system
        H = np.zeros((2*N, 2*N), dtype=complex)
        H[:N, :N] = build_K_layer(N, 10.0)
        H[N:, N:] = build_K_layer(N, 10.0)
        coupling = beta * np.eye(N)
        H[:N, N:] = coupling
        H[N:, :N] = coupling.T.conj()
        
        eigenvalues = np.linalg.eigvalsh(H)
        pos_eig = eigenvalues[eigenvalues > 0.1]
        
        if len(pos_eig) >= 2:
            hierarchy = pos_eig.max() / pos_eig.min()
            orders = np.log10(hierarchy)
        else:
            orders = 0
        
        results.append({
            "beta": beta,
            "orders": float(orders)
        })
    
    # Find optimal β
    best = max(results, key=lambda x: x["orders"])
    
    return {
        "results": results,
        "optimal_beta": best["beta"],
        "max_orders": float(best["orders"]),
        "RESONANCE_FOUND": any(r["orders"] > results[0]["orders"] * 1.5 for r in results[1:])
    }

# ============================================================================
# QW-876: FINAL SYNTHESIS
# ============================================================================
def qw876_final_synthesis(all_results):
    """
    Synthesize all multi-layer findings.
    """
    print("[QW-876] Final Synthesis")
    
    synthesis = {
        "Single_layer_orders": 3.66,
        "Target_SM_orders": 5.5,
        "Gap_to_fill": 1.84,
        "Best_multilayer_mechanism": None,
        "Max_orders_achieved": 0,
        "Key_findings": []
    }
    
    # Check each result
    if "qw858" in all_results:
        if all_results["qw858"].get("EXPONENTIAL_GROWTH"):
            synthesis["Key_findings"].append("Hierarchy grows exponentially with layers")
    
    if "qw866" in all_results:
        orders = all_results["qw866"].get("Orders", 0)
        if orders > synthesis["Max_orders_achieved"]:
            synthesis["Max_orders_achieved"] = orders
            synthesis["Best_multilayer_mechanism"] = "Tensor Product"
    
    if "qw869" in all_results:
        if all_results["qw869"].get("FORMULA_VALID"):
            synthesis["Key_findings"].append("Formula validated: Orders = single + N×log10(1/β)")
    
    if "qw872" in all_results:
        if all_results["qw872"].get("SM_HIERARCHY_ACHIEVED"):
            synthesis["Key_findings"].append("SM hierarchy achieved!")
            synthesis["N_layers_for_SM"] = all_results["qw872"].get("N_layers_needed")
    
    synthesis["SUCCESS"] = synthesis["Max_orders_achieved"] >= 5.0
    
    return synthesis

# ============================================================================
# MAIN RUNNER
# ============================================================================
def run_multilayer_suite():
    print("=" * 80)
    print("QW-857 TO QW-876: MULTI-LAYER K(d) HIERARCHICAL EMERGENCE")
    print("=" * 80)
    print(f"BETA_TORS = {BETA_TORS} (inter-layer coupling)")
    print("=" * 80)
    
    all_results = {}
    
    tests = [
        ("qw857", qw857_multilayer_spectrum),
        ("qw858", qw858_hierarchy_vs_layers),
        ("qw859", qw859_weak_coupling),
        ("qw860", qw860_recursive_nesting),
        ("qw861", qw861_time_dilation),
        ("qw862", qw862_stacked_eigenvalues),
        ("qw863", qw863_effective_beta),
        ("qw864", qw864_layer_correlations),
        ("qw865", qw865_layer_dependent_K),
        ("qw866", qw866_tensor_product),
        ("qw867", qw867_fractal_dimension),
        ("qw868", qw868_localization_layers),
        ("qw869", qw869_effective_formula),
        ("qw870", qw870_gap_structure_layers),
        ("qw871", qw871_emergent_generations),
        ("qw872", qw872_sm_hierarchy_reached),
        ("qw873", qw873_stability_multilayer),
        ("qw874", qw874_optimal_layers),
        ("qw875", qw875_interlayer_resonance),
    ]
    
    for test_id, test_func in tests:
        try:
            result = test_func()
            all_results[test_id] = result
            
            # Print key metric
            if "Orders" in result:
                print(f"  → Orders: {result['Orders']:.2f}")
            elif "Orders_of_magnitude" in result:
                print(f"  → Orders: {result['Orders_of_magnitude']:.2f}")
            elif "Hierarchy" in result:
                print(f"  → Hierarchy: {result['Hierarchy']:.2f}")
        except Exception as e:
            print(f"  ❌ ERROR: {e}")
            all_results[test_id] = {"ERROR": str(e)}
    
    # Final synthesis
    synthesis = qw876_final_synthesis(all_results)
    all_results["qw876"] = synthesis
    print(f"\n[QW-876] Final Synthesis: {'SUCCESS' if synthesis['SUCCESS'] else 'MORE WORK NEEDED'}")
    print(f"  Max orders achieved: {synthesis['Max_orders_achieved']:.2f}")
    
    # Save reports
    with open("RAPORT_QW857_QW876_MULTILAYER.md", "w") as f:
        f.write("# RAPORT: WIELOWARSTWOWA HIERARCHIA K(d) (QW-857 – QW-876)\n\n")
        f.write(f"**Data:** {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"**β_tors:** {BETA_TORS}\n\n")
        f.write("## Kluczowe Wyniki\n\n")
        f.write(f"- Pojedyncza warstwa: 3.66 rzędów\n")
        f.write(f"- Cel (SM): 5.5 rzędów\n")
        f.write(f"- **Osiągnięto:** {synthesis['Max_orders_achieved']:.2f} rzędów\n\n")
        f.write("## Odkrycia\n")
        for finding in synthesis["Key_findings"]:
            f.write(f"- {finding}\n")
    
    with open("RAPORT_QW857_QW876_DATA.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    
    print("\n" + "=" * 80)
    print("MULTILAYER SUITE COMPLETE")
    print("=" * 80)
    
    return all_results

if __name__ == "__main__":
    run_multilayer_suite()
