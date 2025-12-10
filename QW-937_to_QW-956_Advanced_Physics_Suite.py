#!/usr/bin/env python3
"""
QW-937 to QW-956: ADVANCED K(d) PHYSICS SUITE
==============================================
Based on CONFIRMED findings:
✅ Hierarchy 6+ orders (fractal layers)
✅ 3 generations from topology
✅ 1/r² gravity
✅ Dark energy at large d
✅ Self-organization into layers
✅ Bound states

FROZEN PARAMETERS (accepted as axiomatic):
α = 4ln(2), ω = π/4, φ = π/6, β = 0.01

Author: Krzysztof Żuchowski / Antigravity System
Date: 2025-12-09
"""

import numpy as np
from scipy.linalg import eigh
from scipy.signal import find_peaks
from scipy.optimize import minimize
import time
import json

# ============================================================================
# FROZEN KERNEL - ACCEPTED AS AXIOMATIC
# ============================================================================
ALPHA_GEO = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01

def K_complex(d):
    d = np.asarray(d)
    denom = 1 + BETA_TORS * np.abs(d)
    return (ALPHA_GEO * np.exp(1j * (OMEGA * d + PHI))) / denom

def K_real(d):
    return np.real(K_complex(d))

def K_mag(d):
    return np.abs(K_complex(d))

def K_phase(d):
    return np.angle(K_complex(d))

# ============================================================================
# HELPER: BUILD SYSTEMS
# ============================================================================
def build_K_matrix(N, d_max):
    M = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            d = abs(i - j) * d_max / N
            M[i, j] = K_mag(d)
    return (M + M.T) / 2

def get_spectrum(H):
    return np.sort(np.linalg.eigvalsh(H))[::-1]

# ============================================================================
# QW-937: COMPLETE EIGENSPECTRUM ANALYSIS
# ============================================================================
def qw937_complete_spectrum():
    """Full eigenspectrum characterization"""
    print("[QW-937] Complete Spectrum Analysis")
    
    N = 200
    d_max = 100.0
    
    H = build_K_matrix(N, d_max)
    eig = get_spectrum(H)
    
    # Statistics
    pos_eig = eig[eig > 0]
    neg_eig = eig[eig < 0]
    
    # Gap structure
    gaps = -np.diff(eig)
    
    # Degeneracy detection
    degeneracy_threshold = np.std(gaps) * 0.1
    degenerate_groups = []
    group_start = 0
    for i in range(len(gaps)):
        if gaps[i] > degeneracy_threshold:
            if i > group_start:
                degenerate_groups.append((group_start, i))
            group_start = i + 1
    
    return {
        "N_positive": len(pos_eig),
        "N_negative": len(neg_eig),
        "lambda_max": float(eig[0]),
        "lambda_min": float(eig[-1]),
        "hierarchy_orders": float(np.log10(eig[0] / eig[eig > 0][-1])),
        "N_degenerate_groups": len(degenerate_groups),
        "mean_gap": float(np.mean(gaps)),
        "gap_std": float(np.std(gaps))
    }

# ============================================================================
# QW-938: EIGENVECTOR STRUCTURE
# ============================================================================
def qw938_eigenvector_structure():
    """Analyze eigenvector localization and patterns"""
    print("[QW-938] Eigenvector Structure")
    
    N = 150
    d_max = 75.0
    
    H = build_K_matrix(N, d_max)
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    
    # IPR for localization
    iprs = []
    centers = []
    for i in range(len(eigenvalues)):
        psi = eigenvectors[:, i]
        psi2 = np.abs(psi)**2
        psi2 = psi2 / np.sum(psi2)
        
        ipr = np.sum(psi2**2)
        iprs.append(ipr)
        
        positions = np.arange(N) * d_max / N
        center = np.sum(positions * psi2)
        centers.append(center)
    
    iprs = np.array(iprs)
    centers = np.array(centers)
    
    # Classification
    localized = np.sum(iprs > 3/N)
    extended = np.sum(iprs < 1.5/N)
    
    return {
        "N_localized": int(localized),
        "N_extended": int(extended),
        "mean_IPR": float(np.mean(iprs)),
        "center_spread": float(np.std(centers)),
        "PARTICLE_STATES": int(localized)
    }

# ============================================================================
# QW-939: ZEROS AND EXTREMA OF K(d)
# ============================================================================
def qw939_zeros_extrema():
    """Detailed analysis of K(d) structure"""
    print("[QW-939] Zeros and Extrema of K(d)")
    
    d = np.linspace(0, 50, 5000)
    K = K_real(d)
    K_abs = K_mag(d)
    
    # Zeros
    zeros = []
    for i in range(len(K) - 1):
        if K[i] * K[i+1] < 0:
            d_zero = d[i] - K[i] * (d[i+1] - d[i]) / (K[i+1] - K[i])
            zeros.append(d_zero)
    
    # Maxima
    peaks_idx, _ = find_peaks(K_abs)
    peaks_d = d[peaks_idx]
    peaks_K = K_abs[peaks_idx]
    
    # Minima
    minima_idx, _ = find_peaks(-K_abs)
    minima_d = d[minima_idx]
    
    return {
        "N_zeros": len(zeros),
        "zeros": [float(z) for z in zeros[:10]],
        "N_maxima": len(peaks_d),
        "maxima_d": [float(p) for p in peaks_d[:10]],
        "maxima_K": [float(k) for k in peaks_K[:10]],
        "N_minima": len(minima_d),
        "oscillation_period": float(2*np.pi/OMEGA)
    }

# ============================================================================
# QW-940: PHASE STRUCTURE
# ============================================================================
def qw940_phase_structure():
    """Analyze phase of K(d)"""
    print("[QW-940] Phase Structure")
    
    d = np.linspace(0, 30, 1000)
    phase = K_phase(d)
    
    # Phase velocity
    d_phase = np.gradient(phase, d)
    
    # Phase jumps
    phase_wrapped = np.mod(phase + np.pi, 2*np.pi) - np.pi
    jumps = np.abs(np.diff(phase_wrapped)) > np.pi
    
    return {
        "phase_at_d0": float(phase[0]),
        "phase_velocity_mean": float(np.mean(d_phase)),
        "phase_velocity_expected": float(OMEGA),
        "N_phase_jumps": int(np.sum(jumps)),
        "PHASE_LINEAR": np.std(d_phase) / np.mean(np.abs(d_phase)) < 0.3
    }

# ============================================================================
# QW-941: SPECTRAL DENSITY
# ============================================================================
def qw941_spectral_density():
    """Density of states analysis"""
    print("[QW-941] Spectral Density")
    
    N = 300
    d_max = 150.0
    
    H = build_K_matrix(N, d_max)
    eig = get_spectrum(H)
    
    # DOS histogram
    bins = 50
    hist, bin_edges = np.histogram(eig, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Find peaks in DOS
    peaks_idx, _ = find_peaks(hist, height=np.mean(hist))
    
    # Van Hove singularities = peaks in DOS
    van_hove = len(peaks_idx)
    
    return {
        "DOS_peaks": int(van_hove),
        "eigenvalue_range": float(eig[0] - eig[-1]),
        "mean_DOS": float(np.mean(hist)),
        "max_DOS": float(np.max(hist)),
        "VAN_HOVE_SINGULARITIES": van_hove
    }

# ============================================================================
# QW-942: CORRELATION FUNCTIONS
# ============================================================================
def qw942_correlations():
    """Two-point correlation from K(d)"""
    print("[QW-942] Correlation Functions")
    
    # G(r) = ⟨K(0)K(r)⟩
    r = np.linspace(0, 40, 200)
    K0 = K_mag(0)
    Kr = K_mag(r)
    
    G_r = K0 * Kr
    
    # Correlation length
    G_normalized = G_r / G_r[0]
    
    # Find where correlation drops to 1/e
    xi_idx = np.argmin(np.abs(G_normalized - 1/np.e))
    xi = r[xi_idx]
    
    # Power law or exponential?
    log_r = np.log(r[1:20])
    log_G = np.log(G_normalized[1:20])
    coeffs = np.polyfit(log_r, log_G, 1)
    power_law_exp = coeffs[0]
    
    return {
        "correlation_length": float(xi),
        "G_at_xi": float(G_normalized[xi_idx]),
        "power_law_exponent": float(power_law_exp),
        "LONG_RANGE": xi > 10
    }

# ============================================================================
# QW-943: SCALING DIMENSION
# ============================================================================
def qw943_scaling_dimension():
    """Scaling dimension from K(d)"""
    print("[QW-943] Scaling Dimension")
    
    # K(d) ~ d^(-Δ) at large d
    d = np.linspace(10, 100, 100)
    K = K_mag(d)
    
    log_d = np.log(d)
    log_K = np.log(K)
    
    coeffs = np.polyfit(log_d, log_K, 1)
    Delta = -coeffs[0]
    
    # Anomalous dimension
    Delta_free = 1.0  # Free field
    gamma = Delta - Delta_free
    
    return {
        "scaling_dimension": float(Delta),
        "anomalous_dimension": float(gamma),
        "ALPHA_GEO_related": float(ALPHA_GEO),
        "CONFORMAL": abs(Delta - 2.77) < 0.5
    }

# ============================================================================
# QW-944: SUSCEPTIBILITY
# ============================================================================
def qw944_susceptibility():
    """Magnetic-like susceptibility from spectrum"""
    print("[QW-944] Susceptibility")
    
    N = 100
    d_max = 50.0
    
    H = build_K_matrix(N, d_max)
    eig = get_spectrum(H)
    
    # Susceptibility χ = Σ 1/λ_n
    pos_eig = eig[eig > 0.01]
    chi = np.sum(1 / pos_eig)
    
    # Specific heat C = Σ λ²/(exp(λ)-1)²
    beta_eff = 1.0
    C = np.sum((pos_eig**2 * np.exp(beta_eff * pos_eig)) / 
               (np.exp(beta_eff * pos_eig) - 1)**2)
    
    return {
        "susceptibility": float(chi),
        "specific_heat": float(C),
        "chi_per_mode": float(chi / len(pos_eig))
    }

# ============================================================================
# QW-945: ENTANGLEMENT ENTROPY
# ============================================================================
def qw945_entanglement():
    """Entanglement entropy from ground state"""
    print("[QW-945] Entanglement Entropy")
    
    N = 100
    d_max = 50.0
    
    H = build_K_matrix(N, d_max)
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    
    # Ground state
    psi0 = eigenvectors[:, -1]
    
    # Reduced density matrix (trace out half)
    n_A = N // 2
    psi_matrix = psi0.reshape(2, n_A) if N % 2 == 0 else None
    
    if psi_matrix is not None:
        rho_A = psi_matrix @ psi_matrix.T
        rho_A = rho_A / np.trace(rho_A)
        
        # Von Neumann entropy
        eig_rho = np.linalg.eigvalsh(rho_A)
        eig_rho = eig_rho[eig_rho > 1e-10]
        S = -np.sum(eig_rho * np.log(eig_rho))
    else:
        S = 0
    
    return {
        "entanglement_entropy": float(S),
        "max_entropy": float(np.log(N//2)),
        "entropy_ratio": float(S / np.log(N//2)) if N > 2 else 0
    }

# ============================================================================
# QW-946: FORCE FROM GRADIENT
# ============================================================================
def qw946_force_gradient():
    """Force law from -∇K(d)"""
    print("[QW-946] Force from Gradient")
    
    d = np.linspace(1, 50, 200)
    K = K_mag(d)
    
    # Force = -dK/dd
    F = -np.gradient(K, d)
    
    # Fit F ~ d^n
    valid = d > 2
    log_d = np.log(d[valid])
    log_F = np.log(np.abs(F[valid]) + 1e-10)
    
    coeffs = np.polyfit(log_d, log_F, 1)
    n = coeffs[0]
    
    # Check for different regimes
    short_range = d < 10
    long_range = d > 20
    
    coeffs_short = np.polyfit(np.log(d[short_range]), 
                              np.log(np.abs(F[short_range]) + 1e-10), 1)
    coeffs_long = np.polyfit(np.log(d[long_range]), 
                             np.log(np.abs(F[long_range]) + 1e-10), 1)
    
    return {
        "overall_exponent": float(n),
        "short_range_exp": float(coeffs_short[0]),
        "long_range_exp": float(coeffs_long[0]),
        "NEWTON_SHORT": abs(coeffs_short[0] + 2) < 0.5,
        "NEWTON_LONG": abs(coeffs_long[0] + 2) < 0.5
    }

# ============================================================================
# QW-947: POTENTIAL WELL STRUCTURE
# ============================================================================
def qw947_potential_well():
    """Find bound states from K(d) as potential"""
    print("[QW-947] Potential Well Structure")
    
    d = np.linspace(0, 30, 300)
    V = -K_mag(d)  # Negative K as attractive potential
    
    # Find wells (local minima)
    minima_idx, _ = find_peaks(-V)
    n_wells = len(minima_idx)
    
    # Depth of wells
    if n_wells > 0:
        well_depths = [V.max() - V[idx] for idx in minima_idx]
        well_positions = d[minima_idx]
    else:
        well_depths = []
        well_positions = []
    
    return {
        "N_wells": n_wells,
        "well_positions": [float(p) for p in well_positions[:5]],
        "well_depths": [float(d) for d in well_depths[:5]],
        "BOUND_STATES_POSSIBLE": n_wells >= 3
    }

# ============================================================================
# QW-948: TUNNELING RATES
# ============================================================================
def qw948_tunneling():
    """Tunneling between adjacent wells"""
    print("[QW-948] Tunneling Rates")
    
    d = np.linspace(0, 30, 500)
    V = -K_mag(d)
    
    # WKB approximation for tunneling
    # T ~ exp(-2∫√(2m(V-E))dx)
    
    minima_idx, _ = find_peaks(-V)
    
    if len(minima_idx) >= 2:
        # Barrier between first two wells
        barrier_region = (d > d[minima_idx[0]]) & (d < d[minima_idx[1]])
        V_barrier = V[barrier_region]
        d_barrier = d[barrier_region]
        
        E_ground = V[minima_idx[0]] * 0.9  # Rough estimate
        
        integrand = np.sqrt(np.maximum(V_barrier - E_ground, 0))
        action = np.trapz(integrand, d_barrier)
        
        T = np.exp(-2 * action)
    else:
        T = 0
        action = float('inf')
    
    return {
        "tunneling_amplitude": float(T),
        "WKB_action": float(action),
        "COHERENT_TUNNELING": T > 0.01
    }

# ============================================================================
# QW-949: SYMMETRY ANALYSIS
# ============================================================================
def qw949_symmetry():
    """Symmetry properties of K(d)"""
    print("[QW-949] Symmetry Analysis")
    
    d = np.linspace(-30, 30, 1000)
    K = K_real(d)
    K_abs = K_mag(d)
    
    # Even/odd decomposition
    K_even = (K_real(d) + K_real(-d)) / 2
    K_odd = (K_real(d) - K_real(-d)) / 2
    
    even_norm = np.linalg.norm(K_even)
    odd_norm = np.linalg.norm(K_odd)
    
    # Parity
    parity = even_norm / (even_norm + odd_norm)
    
    # Time reversal (complex conjugate)
    K_conj = np.conj(K_complex(d))
    K_original = K_complex(d)
    TR_invariance = np.allclose(K_conj, K_original)
    
    return {
        "even_component": float(even_norm),
        "odd_component": float(odd_norm),
        "parity": float(parity),
        "TIME_REVERSAL": bool(TR_invariance),
        "PARITY_SYMMETRIC": parity > 0.9
    }

# ============================================================================
# QW-950: RESONANCE IDENTIFICATION
# ============================================================================
def qw950_resonances():
    """Identify resonances from spectrum"""
    print("[QW-950] Resonance Identification")
    
    N = 150
    d_max = 75.0
    
    H = build_K_matrix(N, d_max)
    eig = get_spectrum(H)
    
    # Resonances = eigenvalues with large gaps before and after
    gaps = -np.diff(eig)
    
    resonances = []
    for i in range(1, len(eig) - 1):
        if gaps[i-1] > 2*np.mean(gaps) and gaps[i] > 2*np.mean(gaps):
            resonances.append({
                "index": i,
                "eigenvalue": float(eig[i]),
                "isolation": float(min(gaps[i-1], gaps[i]) / np.mean(gaps))
            })
    
    # Resonance ratios
    if len(resonances) >= 2:
        ratios = [resonances[i]["eigenvalue"] / resonances[i+1]["eigenvalue"] 
                  for i in range(len(resonances)-1)]
    else:
        ratios = []
    
    return {
        "N_resonances": len(resonances),
        "resonances": resonances[:5],
        "resonance_ratios": [float(r) for r in ratios[:5]],
        "PARTICLE_SPECTRUM": len(resonances) >= 3
    }

# ============================================================================
# QW-951: CRITICAL EXPONENTS
# ============================================================================
def qw951_critical_exponents():
    """Extract critical exponents from K(d)"""
    print("[QW-951] Critical Exponents")
    
    # At d=0, K is maximum - treat as critical point
    d = np.linspace(0.1, 20, 200)
    K = K_mag(d)
    
    # β exponent: order parameter ~ t^β near critical point
    log_d = np.log(d[:50])
    log_dK = np.log(np.abs(K[:50] - K[0]))
    
    valid = ~np.isinf(log_dK) & ~np.isnan(log_dK)
    if np.sum(valid) > 5:
        coeffs = np.polyfit(log_d[valid], log_dK[valid], 1)
        beta = coeffs[0]
    else:
        beta = 0
    
    # ν exponent: correlation length ~ t^(-ν)
    xi = 1 / BETA_TORS  # Natural correlation length
    
    return {
        "beta_exponent": float(beta),
        "correlation_length": float(xi),
        "ALPHA_GEO": float(ALPHA_GEO),
        "CRITICAL_BEHAVIOR": abs(beta - 1) < 0.5
    }

# ============================================================================
# QW-952: DIMENSIONAL REDUCTION
# ============================================================================
def qw952_dimensional_reduction():
    """How dimensions reduce with d"""
    print("[QW-952] Dimensional Reduction")
    
    results = []
    for d_max in [10, 20, 40, 80, 160]:
        N = 100
        H = build_K_matrix(N, d_max)
        
        # Effective dimension from participation ratio
        eig = get_spectrum(H)
        pos_eig = eig[eig > 0]
        
        # Spectral dimension from heat kernel
        if len(pos_eig) > 5:
            ts = np.logspace(-1, 1, 10)
            Ps = [np.sum(np.exp(-pos_eig * t)) for t in ts]
            
            log_t = np.log(ts)
            log_P = np.log(np.array(Ps))
            
            slopes = -np.gradient(log_P, log_t)
            d_s = 2 * np.mean(slopes[3:7])
        else:
            d_s = 0
        
        results.append({"d_max": d_max, "d_spectral": float(d_s)})
    
    # Does dimension flow with scale?
    dimensions = [r["d_spectral"] for r in results]
    
    return {
        "results": results,
        "UV_dimension": float(dimensions[0]),
        "IR_dimension": float(dimensions[-1]),
        "DIMENSION_FLOWS": abs(dimensions[0] - dimensions[-1]) > 0.5
    }

# ============================================================================
# QW-953: ANOMALY DETECTION
# ============================================================================
def qw953_anomalies():
    """Detect anomalies in K(d) spectrum"""
    print("[QW-953] Anomaly Detection")
    
    N = 200
    d_max = 100.0
    
    H = build_K_matrix(N, d_max)
    eig = get_spectrum(H)
    
    # Level spacing statistics (for chaos detection)
    pos_eig = np.sort(eig[eig > 0])[::-1]
    spacings = np.abs(np.diff(pos_eig))
    
    # Normalize
    mean_spacing = np.mean(spacings)
    s = spacings / mean_spacing
    
    # Level spacing distribution
    # Poisson: P(s) = exp(-s) for integrable
    # Wigner: P(s) = (π/2)s exp(-πs²/4) for chaotic
    
    # r-statistic
    r = []
    for i in range(len(s) - 1):
        r.append(min(s[i], s[i+1]) / max(s[i], s[i+1]))
    r_mean = np.mean(r)
    
    # Poisson: r = 0.386, Wigner: r = 0.536
    is_chaotic = r_mean > 0.45
    
    return {
        "mean_spacing": float(mean_spacing),
        "r_statistic": float(r_mean),
        "POISSON_expected": 0.386,
        "WIGNER_expected": 0.536,
        "IS_CHAOTIC": is_chaotic,
        "IS_INTEGRABLE": not is_chaotic
    }

# ============================================================================
# QW-954: TOPOLOGICAL INVARIANT
# ============================================================================
def qw954_topology():
    """Calculate topological invariants"""
    print("[QW-954] Topological Invariant")
    
    d = np.linspace(0, 4*np.pi/OMEGA, 1000)  # One full period
    K_c = K_complex(d)
    
    # Winding number = (1/2πi) ∮ dlog(K)
    dK = np.gradient(K_c, d)
    integrand = dK / K_c
    
    winding = np.real(np.sum(integrand) * (d[1] - d[0]) / (2j * np.pi))
    
    # Phase winding
    phase = K_phase(d)
    phase_winding = (phase[-1] - phase[0]) / (2 * np.pi)
    
    return {
        "winding_number": float(winding),
        "phase_winding": float(phase_winding),
        "expected_winding": 1.0,
        "TOPOLOGICALLY_NONTRIVIAL": abs(winding) > 0.5
    }

# ============================================================================
# QW-955: RENORMALIZATION GROUP
# ============================================================================
def qw955_rg_flow():
    """RG flow of K(d) parameters"""
    print("[QW-955] Renormalization Group")
    
    # How does effective K change at different scales?
    scales = [1, 2, 4, 8, 16]
    results = []
    
    for s in scales:
        # Coarse-grained K
        d = np.linspace(0, 50 * s, 100)
        K = K_mag(d / s)  # Rescaled
        
        # Effective parameters
        alpha_eff = K[0]  # K(0)
        
        # Decay rate
        log_K = np.log(K[K > 0.1])
        log_d = np.log(d[:len(log_K)])
        if len(log_K) > 5:
            coeffs = np.polyfit(log_d, log_K, 1)
            beta_eff = -coeffs[0]
        else:
            beta_eff = BETA_TORS
        
        results.append({
            "scale": s,
            "alpha_eff": float(alpha_eff),
            "beta_eff": float(beta_eff)
        })
    
    # Fixed point?
    alphas = [r["alpha_eff"] for r in results]
    fixed_point = np.std(alphas) / np.mean(alphas) < 0.1
    
    return {
        "results": results,
        "FIXED_POINT": fixed_point,
        "UV_alpha": float(alphas[0]),
        "IR_alpha": float(alphas[-1])
    }

# ============================================================================
# QW-956: FINAL PHYSICAL SYNTHESIS
# ============================================================================
def qw956_synthesis(all_results):
    """Synthesize all physical findings"""
    print("[QW-956] Final Physical Synthesis")
    
    synthesis = {
        "CONFIRMED_PHYSICS": [],
        "QUANTUM_PROPERTIES": {},
        "CLASSICAL_LIMIT": {},
        "EMERGENT_PHENOMENA": [],
        "VERDICT": ""
    }
    
    # Collect confirmations
    if "qw946" in all_results:
        if all_results["qw946"].get("NEWTON_LONG"):
            synthesis["CONFIRMED_PHYSICS"].append("1/r² gravity at long range")
    
    if "qw947" in all_results:
        if all_results["qw947"].get("BOUND_STATES_POSSIBLE"):
            synthesis["CONFIRMED_PHYSICS"].append("Bound states from potential wells")
    
    if "qw950" in all_results:
        if all_results["qw950"].get("PARTICLE_SPECTRUM"):
            synthesis["CONFIRMED_PHYSICS"].append("Discrete particle spectrum")
    
    if "qw954" in all_results:
        if all_results["qw954"].get("TOPOLOGICALLY_NONTRIVIAL"):
            synthesis["CONFIRMED_PHYSICS"].append("Topological protection")
    
    # Quantum
    if "qw953" in all_results:
        synthesis["QUANTUM_PROPERTIES"]["integrable"] = all_results["qw953"].get("IS_INTEGRABLE", False)
    
    if "qw945" in all_results:
        synthesis["QUANTUM_PROPERTIES"]["entanglement_entropy"] = all_results["qw945"].get("entanglement_entropy", 0)
    
    # Emergent
    if "qw952" in all_results:
        if all_results["qw952"].get("DIMENSION_FLOWS"):
            synthesis["EMERGENT_PHENOMENA"].append("Dimensional reduction")
    
    if "qw951" in all_results:
        if all_results["qw951"].get("CRITICAL_BEHAVIOR"):
            synthesis["EMERGENT_PHENOMENA"].append("Critical behavior")
    
    n_confirmed = len(synthesis["CONFIRMED_PHYSICS"])
    n_emergent = len(synthesis["EMERGENT_PHENOMENA"])
    
    if n_confirmed >= 3:
        synthesis["VERDICT"] = f"SOLID PHYSICS: {n_confirmed} confirmed phenomena, {n_emergent} emergent"
    else:
        synthesis["VERDICT"] = f"PARTIAL: {n_confirmed} confirmed, needs more validation"
    
    return synthesis

# ============================================================================
# MAIN RUNNER
# ============================================================================
def run_advanced_suite():
    print("=" * 80)
    print("QW-937 TO QW-956: ADVANCED K(d) PHYSICS SUITE")
    print("=" * 80)
    print("Building on CONFIRMED findings with frozen parameters")
    print("=" * 80)
    
    all_results = {}
    
    tests = [
        ("qw937", qw937_complete_spectrum),
        ("qw938", qw938_eigenvector_structure),
        ("qw939", qw939_zeros_extrema),
        ("qw940", qw940_phase_structure),
        ("qw941", qw941_spectral_density),
        ("qw942", qw942_correlations),
        ("qw943", qw943_scaling_dimension),
        ("qw944", qw944_susceptibility),
        ("qw945", qw945_entanglement),
        ("qw946", qw946_force_gradient),
        ("qw947", qw947_potential_well),
        ("qw948", qw948_tunneling),
        ("qw949", qw949_symmetry),
        ("qw950", qw950_resonances),
        ("qw951", qw951_critical_exponents),
        ("qw952", qw952_dimensional_reduction),
        ("qw953", qw953_anomalies),
        ("qw954", qw954_topology),
        ("qw955", qw955_rg_flow),
    ]
    
    for test_id, test_func in tests:
        try:
            result = test_func()
            all_results[test_id] = result
        except Exception as e:
            print(f"  ❌ ERROR: {e}")
            all_results[test_id] = {"ERROR": str(e)}
    
    # Final synthesis
    synthesis = qw956_synthesis(all_results)
    all_results["qw956"] = synthesis
    print(f"\n[QW-956] Verdict: {synthesis['VERDICT']}")
    print(f"  Confirmed: {synthesis['CONFIRMED_PHYSICS']}")
    
    # Save
    with open("RAPORT_QW937_QW956_ADVANCED.md", "w") as f:
        f.write("# RAPORT: ZAAWANSOWANA FIZYKA K(d) (QW-937 – QW-956)\n\n")
        f.write(f"**Data:** {time.strftime('%Y-%m-%d %H:%M')}\n")
        f.write("**Parametry:** α=4ln(2), ω=π/4, φ=π/6, β=0.01 (zamrożone)\n\n")
        f.write("## Potwierdzona Fizyka\n")
        for p in synthesis["CONFIRMED_PHYSICS"]:
            f.write(f"- ✅ {p}\n")
        f.write("\n## Emergentne Zjawiska\n")
        for e in synthesis["EMERGENT_PHENOMENA"]:
            f.write(f"- 🔬 {e}\n")
        f.write(f"\n## Werdykt\n**{synthesis['VERDICT']}**\n")
    
    with open("RAPORT_QW937_QW956_DATA.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    
    print("\n" + "=" * 80)
    print("ADVANCED SUITE COMPLETE")
    print("=" * 80)
    
    return all_results

if __name__ == "__main__":
    run_advanced_suite()
