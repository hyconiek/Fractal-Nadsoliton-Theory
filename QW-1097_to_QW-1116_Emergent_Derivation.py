#!/usr/bin/env python3
"""
QW-1097 to QW-1116: EMERGENTNE WYPROWADZENIE POZYCJI CZĄSTEK
============================================================
CEL: Pokazać DLACZEGO elektron ma d=6, mion d=3.5, top d=0 
     wynikające Z DYNAMIKI K(d), nie z fittingu.

HIPOTEZY DO TESTOWANIA:
1. Pozycje = stany stacjonarne Hamiltonianu K(d)
2. Pozycje = stabilne topologie (winding numbers)
3. Pozycje = minimα energii swobodnej
4. Pozycje = resonanse 4-bitowego oscylatora
5. Pozycje = kwantowane orbity w potencjale K(d)

Author: Krzysztof Żuchowski / Antigravity System
Date: 2025-12-10
"""

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize, brentq
from scipy.integrate import solve_ivp, quad
import time
import json

# ============================================================================
# FROZEN KERNEL
# ============================================================================
ALPHA_GEO = 4 * np.log(2)
OMEGA = np.pi / 4
PHI = np.pi / 6
BETA_TORS = 0.01
N_OCTAVES = 12

def K(d):
    d = np.asarray(d)
    return ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + BETA_TORS * np.abs(d))

def K_complex(d):
    d = np.asarray(d)
    return ALPHA_GEO * np.exp(1j * (OMEGA * d + PHI)) / (1 + BETA_TORS * np.abs(d))

# Znane masy (MeV) i pozycje
PARTICLES = {
    "top": {"mass": 172760, "d": 0.00, "winding": 0},
    "bottom": {"mass": 4180, "d": 1.75, "winding": 1},
    "tau": {"mass": 1777, "d": 2.25, "winding": 2},
    "muon": {"mass": 105.7, "d": 3.50, "winding": 3},
    "electron": {"mass": 0.511, "d": 6.00, "winding": 4},
}

# ============================================================================
# QW-1097: BOHR-LIKE QUANTIZATION IN K(d) POTENTIAL
# ============================================================================
def qw1097_bohr_quantization():
    """Bohr quantization: ∫p dd = 2πn, derive allowed d values"""
    print("[QW-1097] Bohr Quantization")
    
    # "Momentum" from K(d): p(d) ~ sqrt(K(d) - E)
    # For bound state E < min(K), orbit where K(d) > E
    
    def action_integral(E, d_max=20):
        """Compute ∫sqrt(|K(d) - E|) dd"""
        d_vals = np.linspace(0, d_max, 500)
        K_vals = K(d_vals)
        integrand = np.sqrt(np.abs(K_vals - E) + 1e-10)
        return np.trapz(integrand, d_vals)
    
    # Find quantized energies where action = n × 2π
    quantized_positions = []
    
    for n in range(1, 20):
        target_action = n * 2 * np.pi
        
        # Find E such that action = target
        def objective(E):
            return action_integral(E) - target_action
        
        try:
            E_n = brentq(objective, -3, 3)
            
            # Find classical turning point (where K(d) = E_n)
            d_vals = np.linspace(0, 15, 300)
            K_vals = K(d_vals)
            crossings = np.where(np.diff(np.sign(K_vals - E_n)))[0]
            
            if len(crossings) > 0:
                d_turn = d_vals[crossings[0]]
                quantized_positions.append({
                    "n": n,
                    "E": float(E_n),
                    "d_turning": float(d_turn)
                })
        except:
            pass
    
    # Extract d values
    d_values = [p["d_turning"] for p in quantized_positions]
    
    # Match to known positions
    matches = []
    for name, data in PARTICLES.items():
        d_exp = data["d"]
        closest = min(d_values, key=lambda x: abs(x - d_exp), default=None)
        if closest is not None:
            error = abs(closest - d_exp)
            matches.append((name, d_exp, closest, error))
    
    return {
        "quantized_levels": quantized_positions[:10],
        "matches": matches,
        "INSIGHT": "Bohr quantization gives allowed d values"
    }

# ============================================================================
# QW-1098: WINDING NUMBER CONSTRAINT
# ============================================================================
def qw1098_winding_constraint():
    """Winding number w = 0, 1, 2, 3, 4 → d = f(w)"""
    print("[QW-1098] Winding Number → Position")
    
    # Hypothesis: d = A × w^2 or d = B × w × log(w+1) 
    # to match d = 0, 1.75, 2.25, 3.5, 6.0
    
    winding_numbers = [0, 1, 2, 3, 4]
    expected_d = [0, 1.75, 2.25, 3.5, 6.0]
    
    # Try different formulae
    results = {}
    
    # Formula 1: d = A × w²
    def f1(A, w): return A * w**2
    A_opt = np.mean([d/w**2 if w > 0 else 0 for d, w in zip(expected_d[1:], winding_numbers[1:])])
    pred_1 = [f1(A_opt, w) for w in winding_numbers]
    err_1 = np.mean([abs(p - e) for p, e in zip(pred_1, expected_d)])
    results["d = A×w²"] = {"A": float(A_opt), "predictions": [float(p) for p in pred_1], "error": float(err_1)}
    
    # Formula 2: d = A × w × (1 + B×w)
    def objective2(params):
        A, B = params
        pred = [A * w * (1 + B * w) if w > 0 else 0 for w in winding_numbers]
        return sum((p - e)**2 for p, e in zip(pred, expected_d))
    
    from scipy.optimize import minimize
    res = minimize(objective2, [1, 0.5])
    A2, B2 = res.x
    pred_2 = [A2 * w * (1 + B2 * w) if w > 0 else 0 for w in winding_numbers]
    err_2 = np.mean([abs(p - e) for p, e in zip(pred_2, expected_d)])
    results["d = A×w×(1+B×w)"] = {"A": float(A2), "B": float(B2), "predictions": [float(p) for p in pred_2], "error": float(err_2)}
    
    # Formula 3: d from K(d) winding phase
    # Accumulated phase = winding × some unit
    phase_unit = np.pi / 4  # From OMEGA
    d_from_phase = [w * phase_unit * 2 for w in winding_numbers]  # 2 = period factor
    err_3 = np.mean([abs(p - e) for p, e in zip(d_from_phase, expected_d)])
    results["d = 2×ω×w"] = {"predictions": [float(p) for p in d_from_phase], "error": float(err_3)}
    
    # Best formula
    best = min(results.items(), key=lambda x: x[1]["error"])
    
    return {
        "formulas_tested": results,
        "best_formula": best[0],
        "best_error": float(best[1]["error"]),
        "INSIGHT": f"Best: {best[0]} with error {best[1]['error']:.2f}"
    }

# ============================================================================
# QW-1099: STABLE TOPOLOGICAL CONFIGURATIONS
# ============================================================================
def qw1099_topological_stability():
    """Find d values where topology is stable"""
    print("[QW-1099] Topological Stability")
    
    # Topological energy: E_topo = ∫|∇θ|² dx where θ = arg(K_complex)
    
    def topological_energy(d_center, width=1.0):
        d_vals = np.linspace(d_center - width, d_center + width, 50)
        phases = np.angle(K_complex(d_vals))
        # Gradient of phase
        grad_phase = np.gradient(phases, d_vals)
        return np.sum(grad_phase**2)
    
    # Find local minima of topological energy
    d_range = np.linspace(0.5, 15, 100)
    E_topo = [topological_energy(d) for d in d_range]
    
    # Find minima
    minima_idx = []
    for i in range(1, len(E_topo) - 1):
        if E_topo[i] < E_topo[i-1] and E_topo[i] < E_topo[i+1]:
            minima_idx.append(i)
    
    stable_d = [d_range[i] for i in minima_idx]
    
    # Match to particles
    matches = []
    for name, data in PARTICLES.items():
        d_exp = data["d"]
        if len(stable_d) > 0:
            closest = min(stable_d, key=lambda x: abs(x - d_exp))
            matches.append((name, d_exp, closest, abs(closest - d_exp)))
    
    return {
        "stable_positions": [float(d) for d in stable_d[:10]],
        "matches": matches,
        "INSIGHT": f"Found {len(stable_d)} topologically stable positions"
    }

# ============================================================================
# QW-1100: 4-BIT STATE ENUMERATION
# ============================================================================
def qw1100_4bit_states():
    """Enumerate 4-bit states and find energy ordering"""
    print("[QW-1100] 4-bit State Enumeration")
    
    # 4 bits = 16 states. Each state has energy from K(d)
    # Hypothesis: state energy = sum of K at bit positions
    
    states = []
    for s in range(16):
        bits = [(s >> i) & 1 for i in range(4)]
        n_on = sum(bits)  # Number of ON bits
        
        # Energy = sum of (-1)^bit[i] × K(i) 
        energy = 0
        for i, b in enumerate(bits):
            energy += (2*b - 1) * K(i * 0.25)  # Bit contributes ±K
        
        # Position: determined by bit pattern
        # Try: d = sum(bit[i] × (i+1) × 0.25) / (n_on + 1)
        d = sum(b * (i + 1) * 0.25 for i, b in enumerate(bits)) 
        
        states.append({
            "state": format(s, '04b'),
            "n_on": n_on,
            "energy": float(energy),
            "d": float(d)
        })
    
    # Sort by energy (lowest = most stable = highest mass = smallest d)
    states_sorted = sorted(states, key=lambda x: x["energy"])
    
    # Take lowest 5 as "particle states"
    particle_states = states_sorted[:5]
    
    return {
        "all_states": states,
        "lowest_energy_states": particle_states,
        "d_values": [s["d"] for s in particle_states],
        "INSIGHT": "4-bit states give discrete d values"
    }

# ============================================================================
# QW-1101: CHANNEL DECOUPLING SEQUENCE
# ============================================================================
def qw1101_channel_decoupling():
    """d increases as channels decouple one by one"""
    print("[QW-1101] Channel Decoupling")
    
    # Start with 4 coupled channels (top quark, d=0)
    # Decouple 1 channel → d increases by Δ
    # Energy drops by factor 4 per channel
    
    # If channels couple with strength K(Δ), then:
    # Total coupling = ∑ K(k×Δ) for k coupled channels
    
    def total_coupling(n_channels, delta):
        return sum(K(k * delta) for k in range(n_channels))
    
    # Find delta such that decoupling matches known d values
    def objective(delta):
        d_pred = [0]  # Top at d=0
        for n in [3, 2, 1, 0]:  # Progressively fewer channels
            d_next = d_pred[-1] + delta * (4 - n) / 4
            d_pred.append(d_next)
        
        expected = [0, 1.75, 3.5, 5.0, 6.0]
        return sum((p - e)**2 for p, e in zip(d_pred, expected))
    
    from scipy.optimize import minimize_scalar
    res = minimize_scalar(objective, bounds=(0.5, 3), method='bounded')
    delta_opt = res.x
    
    # Predict positions
    d_pred = [0]
    for n in [3, 2, 1, 0]:
        d_next = d_pred[-1] + delta_opt * (4 - n) / 4
        d_pred.append(d_next)
    
    return {
        "delta_optimal": float(delta_opt),
        "predicted_d": [float(d) for d in d_pred],
        "expected_d": [0, 1.75, 3.5, 5.0, 6.0],
        "INSIGHT": f"Channel decoupling with Δ={delta_opt:.2f} gives d sequence"
    }

# ============================================================================
# QW-1102: INFORMATION PROCESSING DEPTH
# ============================================================================
def qw1102_processing_depth():
    """d = number of bits processed by Nadsoliton"""
    print("[QW-1102] Information Processing Depth")
    
    # Hypothesis: each particle corresponds to processing N bits
    # Mass decreases as more bits are processed (entropy increases)
    
    # m / m_Planck = 4^(-N_bits)
    # From formula: d = N_bits (in units of 0.25 per sub-bit)
    
    # Top: 0 bits processed → d=0
    # Electron: 24 bits processed → d = 24 × 0.25 = 6
    
    bits_per_step = 0.25  # 1 sub-bit
    
    particles_pred = {}
    for name, data in PARTICLES.items():
        d_exp = data["d"]
        n_bits = d_exp / bits_per_step
        # "Bit address" in 4-bit word
        word = int(n_bits // 4)
        sub_bit = int(n_bits % 4)
        particles_pred[name] = {
            "d_expected": d_exp,
            "n_bits_total": float(n_bits),
            "n_words": word,
            "sub_bit": sub_bit,
            "address": f"({word}, {sub_bit})"
        }
    
    return {
        "particles": particles_pred,
        "INSIGHT": "Each particle = unique information processing depth (Word, Sub-bit)"
    }

# ============================================================================
# QW-1103: EIGENSTATE SELECTION
# ============================================================================
def qw1103_eigenstate_selection():
    """Which eigenstates of K-Hamiltonian correspond to particles?"""
    print("[QW-1103] Eigenstate Selection")
    
    N = 50
    sites = np.linspace(0, 12, N)
    
    # Build Hamiltonian
    H = np.zeros((N, N))
    for i in range(N):
        H[i, i] = -K(sites[i])  # On-site potential
        if i > 0:
            H[i, i-1] = -1  # Hopping
        if i < N - 1:
            H[i, i+1] = -1
    H = (H + H.T) / 2
    
    eigenvalues, eigenvectors = np.linalg.eigh(H)
    
    # Select eigenstates by criteria
    selected = []
    
    for k in range(min(10, N)):
        vec = eigenvectors[:, k]
        
        # Peak position
        peak_idx = np.argmax(np.abs(vec))
        d_peak = sites[peak_idx]
        
        # Winding (number of nodes)
        nodes = np.sum(np.diff(np.sign(vec)) != 0)
        
        # Localization (IPR)
        ipr = np.sum(vec**4)
        
        selected.append({
            "eigenvalue": float(eigenvalues[k]),
            "d_peak": float(d_peak),
            "nodes": int(nodes),
            "IPR": float(ipr)
        })
    
    # Match to particles
    d_peaks = [s["d_peak"] for s in selected]
    matches = []
    for name, data in PARTICLES.items():
        d_exp = data["d"]
        closest = min(d_peaks, key=lambda x: abs(x - d_exp))
        matches.append((name, d_exp, closest, abs(closest - d_exp)))
    
    return {
        "selected_states": selected,
        "matches": matches,
        "INSIGHT": "Eigenstates peak at specific d values"
    }

# ============================================================================
# QW-1104: FREE ENERGY MINIMIZATION
# ============================================================================
def qw1104_free_energy():
    """F = E - TS, find d that minimizes free energy"""
    print("[QW-1104] Free Energy Minimization")
    
    def free_energy(d, T=1.0):
        E = -K(d)  # Energy (potential)
        S = ALPHA_GEO * d  # Entropy increases with d (more disorder)
        F = E - T * S
        return F
    
    # Find minima at different temperatures
    minima_by_T = {}
    
    for T in [0.1, 0.5, 1.0, 2.0]:
        local_minima = []
        for d_init in np.arange(0.5, 12, 1):
            res = minimize(lambda x: free_energy(x[0], T), [d_init], bounds=[(0, 12)])
            if res.success:
                d_min = res.x[0]
                if not any(abs(d_min - m) < 0.2 for m in local_minima):
                    local_minima.append(d_min)
        minima_by_T[f"T={T}"] = sorted(local_minima)
    
    # Best temperature
    best_T = None
    best_match = 0
    for T_key, minima in minima_by_T.items():
        matches = sum(1 for d in [0, 1.75, 3.5, 6.0] if any(abs(d - m) < 0.3 for m in minima))
        if matches > best_match:
            best_match = matches
            best_T = T_key
    
    return {
        "minima_by_temperature": minima_by_T,
        "best_temperature": best_T,
        "INSIGHT": f"Free energy minima at {best_T} match {best_match}/4 particles"
    }

# ============================================================================
# QW-1105: QUANTIZED ORBIT RADII
# ============================================================================
def qw1105_orbit_radii():
    """Bohr-like: r_n ~ n² → d_n ~ n²"""
    print("[QW-1105] Quantized Orbit Radii")
    
    # If d ~ n² where n is principal quantum number
    # d = A × (n + δ)² for some offset δ
    
    n_values = [0, 1, 2, 3, 4]
    d_expected = [0, 1.75, 2.25, 3.5, 6.0]
    
    # Fit: d = A × n²
    # But d(1) = 1.75, d(2) = 2.25 doesn't follow n²
    
    # Alternative: d ~ cumulative sum
    # d(n) = Σ k from k=0 to n-1
    
    # Or: d(n) = Σ (α_k) where α_k = contribution of k-th quantum number
    
    # From data:
    # Δd(0→1) = 1.75
    # Δd(1→2) = 0.50
    # Δd(2→3) = 1.25
    # Δd(3→4) = 2.50
    
    deltas = [d_expected[i+1] - d_expected[i] for i in range(len(d_expected)-1)]
    
    # Pattern in deltas?
    # 1.75, 0.5, 1.25, 2.5 → not obvious
    # Ratio: 0.5/1.75 = 0.29, 1.25/0.5 = 2.5, 2.5/1.25 = 2
    
    # Alternative: 0.25 × [7, 2, 5, 10] steps
    steps_in_units = [d / 0.25 for d in deltas]  # [7, 2, 5, 10]
    
    return {
        "n_values": n_values,
        "d_expected": d_expected,
        "deltas": [float(d) for d in deltas],
        "steps_in_0.25_units": [float(s) for s in steps_in_units],
        "INSIGHT": "Steps are 7, 2, 5, 10 in units of 0.25 — not simple pattern"
    }

# ============================================================================
# QW-1106: RG FLOW FIXED POINTS
# ============================================================================
def qw1106_rg_fixed_points():
    """RG flow: d(s) evolves with scale s, find fixed points"""
    print("[QW-1106] RG Fixed Points")
    
    def beta_function(d):
        """β(d) = dd/ds from K(d) structure"""
        dd = 0.01
        dK = (K(d + dd) - K(d - dd)) / (2 * dd)
        return -dK / (ALPHA_GEO + 0.01)  # Normalize
    
    # RG flow: dd/ds = β(d)
    # Fixed points: β(d*) = 0
    
    d_vals = np.linspace(0.1, 12, 200)
    beta_vals = [beta_function(d) for d in d_vals]
    
    # Find zeros (fixed points)
    fixed_points = []
    for i in range(1, len(beta_vals)):
        if beta_vals[i] * beta_vals[i-1] < 0:
            d_fixed = d_vals[i]
            fixed_points.append(d_fixed)
    
    # Match to particles
    matches = []
    for name, data in PARTICLES.items():
        d_exp = data["d"]
        if len(fixed_points) > 0:
            closest = min(fixed_points, key=lambda x: abs(x - d_exp))
            matches.append((name, d_exp, closest, abs(closest - d_exp)))
    
    return {
        "fixed_points": [float(f) for f in fixed_points[:10]],
        "matches": matches,
        "INSIGHT": f"Found {len(fixed_points)} RG fixed points"
    }

# ============================================================================
# QW-1107: RESONANCE CONDITION
# ============================================================================
def qw1107_resonance():
    """Particles at resonance: dK/dd = 0 or K(d) = integer × α"""
    print("[QW-1107] Resonance Condition")
    
    d_vals = np.linspace(0, 12, 300)
    K_vals = K(d_vals)
    
    # Condition 1: K(d) = n × α/4 (quantized)
    resonance_d = []
    for n in range(-10, 10):
        target_K = n * ALPHA_GEO / 4
        idx = np.argmin(np.abs(K_vals - target_K))
        if np.abs(K_vals[idx] - target_K) < 0.1:
            resonance_d.append(d_vals[idx])
    
    # Remove duplicates
    resonance_d = sorted(set([round(d, 1) for d in resonance_d]))
    
    # Match
    matches = []
    for name, data in PARTICLES.items():
        d_exp = data["d"]
        if len(resonance_d) > 0:
            closest = min(resonance_d, key=lambda x: abs(x - d_exp))
            matches.append((name, d_exp, closest, abs(closest - d_exp)))
    
    return {
        "resonance_positions": [float(r) for r in resonance_d[:15]],
        "matches": matches,
        "INSIGHT": "Resonances where K(d) = n × α/4"
    }

# ============================================================================
# QW-1108: INFORMATION ENTROPY QUANTIZATION
# ============================================================================
def qw1108_entropy_quantization():
    """S(d) = ∫_0^d ln|K| dd', find where S = n × ln(2)"""
    print("[QW-1108] Entropy Quantization")
    
    d_vals = np.linspace(0.1, 12, 300)
    
    # Cumulative "entropy"
    entropy = []
    for d in d_vals:
        d_int = np.linspace(0.1, d, 50)
        K_int = np.abs(K(d_int)) + 1e-10
        S = np.trapz(np.log(K_int), d_int)
        entropy.append(S)
    
    entropy = np.array(entropy)
    
    # Normalize by ln(2)
    entropy_bits = entropy / np.log(2)
    
    # Find where entropy = 0, 4, 8, 16, 24 bits
    target_bits = [0, 4, 8, 12, 16, 20, 24]
    quantized_d = []
    for n in target_bits:
        idx = np.argmin(np.abs(entropy_bits - n))
        quantized_d.append(d_vals[idx])
    
    # Match (map 4, 8, 12, 16, 24 to particles)
    # Top = 0 bits, Electron = 24 bits?
    
    return {
        "target_bits": target_bits,
        "quantized_d": [float(d) for d in quantized_d],
        "entropy_at_6": float(np.interp(6.0, d_vals, entropy_bits)),
        "INSIGHT": "Entropy-quantized positions"
    }

# ============================================================================
# QW-1109: ACTION QUANTIZATION (IMPROVED)
# ============================================================================
def qw1109_action_improved():
    """Improved: S = ∫ K dd = n × (2πα), find d_n"""
    print("[QW-1109] Action Quantization (Improved)")
    
    d_vals = np.linspace(0, 12, 300)
    K_vals = K(d_vals)
    
    # Cumulative action
    action = np.cumsum(K_vals) * (d_vals[1] - d_vals[0])
    
    # Quantum of action = 2π × α or just π
    action_quantum = np.pi
    
    # Find d where action = n × quantum
    quantized_d = []
    for n in range(1, 30):
        target = n * action_quantum
        idx = np.argmin(np.abs(action - target))
        quantized_d.append(d_vals[idx])
    
    # Round to 0.25 and remove duplicates
    quantized_d = sorted(set([round(d * 4) / 4 for d in quantized_d]))
    
    # Match
    matches = []
    for name, data in PARTICLES.items():
        d_exp = data["d"]
        if len(quantized_d) > 0:
            closest = min(quantized_d, key=lambda x: abs(x - d_exp))
            if abs(closest - d_exp) < 0.3:
                matches.append((name, d_exp, closest, abs(closest - d_exp)))
    
    return {
        "action_quantum": float(action_quantum),
        "quantized_positions": [float(d) for d in quantized_d[:12]],
        "matches": matches,
        "match_count": len(matches),
        "INSIGHT": f"Action quantization matches {len(matches)}/5 particles"
    }

# ============================================================================
# QW-1110: STABILITY ANALYSIS
# ============================================================================
def qw1110_stability():
    """d²V/dd² > 0 at particle positions (stable minima)"""
    print("[QW-1110] Stability Analysis")
    
    d_vals = np.linspace(0.1, 12, 200)
    
    # Second derivative of potential V = -K
    d2V = []
    dd = 0.01
    for d in d_vals:
        K_plus = K(d + dd)
        K_zero = K(d)
        K_minus = K(d - dd)
        d2K = (K_plus - 2*K_zero + K_minus) / dd**2
        d2V.append(-d2K)  # V = -K, so d²V/dd² = -d²K/dd²
    
    d2V = np.array(d2V)
    
    # Stable: d2V > 0
    stable_mask = d2V > 0
    stable_regions = d_vals[stable_mask]
    
    # Check stability at particle positions
    stability_check = {}
    for name, data in PARTICLES.items():
        d_exp = data["d"]
        d2V_at_d = np.interp(d_exp, d_vals, d2V)
        stability_check[name] = {
            "d": float(d_exp),
            "d2V": float(d2V_at_d),
            "stable": bool(d2V_at_d > 0)
        }
    
    n_stable = sum(1 for v in stability_check.values() if v["stable"])
    
    return {
        "stability_at_particles": stability_check,
        "n_stable": n_stable,
        "INSIGHT": f"{n_stable}/5 particle positions are stable"
    }

# ============================================================================
# SYNTHESIS
# ============================================================================
def qw1116_synthesis(all_results):
    """Final synthesis: which mechanism derives positions?"""
    print("[QW-1116] Synthesis")
    
    # Count matches for each mechanism
    mechanism_scores = {}
    
    for key, result in all_results.items():
        if "matches" in result:
            matches = result["matches"]
            if isinstance(matches, list) and len(matches) > 0:
                if isinstance(matches[0], tuple):
                    good_matches = sum(1 for m in matches if len(m) > 3 and m[3] < 0.3)
                    mechanism_scores[key] = good_matches
        if "match_count" in result:
            mechanism_scores[key] = result["match_count"]
    
    sorted_scores = sorted(mechanism_scores.items(), key=lambda x: x[1], reverse=True)
    
    synthesis = {
        "mechanism_scores": mechanism_scores,
        "best_mechanisms": sorted_scores[:5],
        "BEST": sorted_scores[0] if sorted_scores else ("none", 0),
        "FULL_EMERGENCE": sorted_scores[0][1] >= 4 if sorted_scores else False,
        "VERDICT": ""
    }
    
    if synthesis["FULL_EMERGENCE"]:
        synthesis["VERDICT"] = f"SUKCES: {synthesis['BEST'][0]} wyprowadza {synthesis['BEST'][1]}/5 pozycji!"
    else:
        synthesis["VERDICT"] = f"CZĘŚCIOWY: Najlepszy = {synthesis['BEST'][0]} ({synthesis['BEST'][1]}/5)"
    
    return synthesis

# ============================================================================
# MAIN
# ============================================================================
def run_derivation_suite():
    print("=" * 80)
    print("QW-1097 TO QW-1116: EMERGENTNE WYPROWADZENIE POZYCJI")
    print("=" * 80)
    
    all_results = {}
    
    tests = [
        ("qw1097", qw1097_bohr_quantization),
        ("qw1098", qw1098_winding_constraint),
        ("qw1099", qw1099_topological_stability),
        ("qw1100", qw1100_4bit_states),
        ("qw1101", qw1101_channel_decoupling),
        ("qw1102", qw1102_processing_depth),
        ("qw1103", qw1103_eigenstate_selection),
        ("qw1104", qw1104_free_energy),
        ("qw1105", qw1105_orbit_radii),
        ("qw1106", qw1106_rg_fixed_points),
        ("qw1107", qw1107_resonance),
        ("qw1108", qw1108_entropy_quantization),
        ("qw1109", qw1109_action_improved),
        ("qw1110", qw1110_stability),
    ]
    
    for test_id, test_func in tests:
        try:
            result = test_func()
            all_results[test_id] = result
            print(f"  → {result.get('INSIGHT', '')}")
        except Exception as e:
            print(f"  ❌ {e}")
            all_results[test_id] = {"ERROR": str(e)}
    
    # Synthesis
    synthesis = qw1116_synthesis(all_results)
    all_results["qw1116"] = synthesis
    print(f"\n[QW-1116] {synthesis['VERDICT']}")
    
    # Save
    with open("RAPORT_QW1097_QW1116_DERIVATION.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    
    print("\n" + "=" * 80)
    
    return all_results

if __name__ == "__main__":
    run_derivation_suite()
