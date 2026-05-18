#!/usr/bin/env python3
"""P2017 S967 strict Cutkosky quadrature-ansatz provenance witness.

Corrective next honest step after P2016: make the previously introduced
quadrature tensor construction epistemically precise.  The script computes a
finite strict-kernel quadrature tensor candidate for the graviton->gauge_gauge
channel, but it now explicitly records the missing provenance gates that prevent
promoting the candidate to a backend-derived loop-amplitude theorem.

No false pass: the numerical quadrature may pass its internal stability checks,
but the exported result_kind remains OPEN until a strict-side Feynman-rule or
loop-integrand derivation is available and the normalization no longer depends
on upstream cut-channel values.
"""
from __future__ import annotations

import json
import math
import platform
from pathlib import Path
from typing import Any

import numpy as np
import scipy.integrate as integrate
import scipy.linalg as la
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p2017_s967_strict_cutkosky_backend_loop_amplitude_tensor_quadrature_witness.json"
TS = "2026-05-18T00:00:00+00:00"

STRICT_PARAMS = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 1.8}
CHANNELS = ["gg", "gh", "hh", "gx"]
CHANNEL_OFFSETS = {"gg": 0.0, "gh": 0.25, "hh": 0.50, "gx": 0.75}


def load(name: str) -> dict[str, Any]:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def strict_kernel(d: float) -> float:
    omega = STRICT_PARAMS["omega"]
    phi = STRICT_PARAMS["phi"]
    beta = STRICT_PARAMS["beta"]
    eta = STRICT_PARAMS["eta"]
    return math.cos(omega * d + phi) / (1.0 + beta * (d**eta))


def channel_distance(channel: str, s_value: float, x: float) -> float:
    """Strict-side candidate distance map for finite phase-space quadrature.

    The offsets distinguish intermediate-state classes without importing legacy
    parameters.  They are still an ansatz-level channel chart, not a derivation
    from strict Feynman rules.
    """
    return math.sqrt(s_value) * (1.0 + x + CHANNEL_OFFSETS[channel])


def angular_basis(x: float) -> np.ndarray:
    # Minimal rank-3 ansatz basis used only to build an auditable tensor candidate.
    return np.array([1.0, 2.0 * x - 1.0, 6.0 * x * (1.0 - x)], dtype=float)


def amplitude_candidate(channel: str, s_value: float, x: float, cut_channel: float, im_m: float) -> float:
    """Explicit strict-kernel quadrature amplitude candidate.

    This is intentionally named as a candidate, not a derived backend amplitude.
    Its normalization is calibrated to already exported P2004 channel cuts, so it
    is not independent evidence for full Cutkosky closure.
    """
    norm = math.sqrt(max(cut_channel, 0.0) / (max(im_m, 1e-300)))
    measure = math.sqrt(max(x * (1.0 - x), 0.0))
    return norm * strict_kernel(channel_distance(channel, s_value, x)) * measure


def tensor_integrand_component(i: int, j: int, channel: str, s_value: float, cut_channel: float, im_m: float) -> Any:
    def f(x: float) -> float:
        b = angular_basis(x)
        a = amplitude_candidate(channel, s_value, x, cut_channel, im_m)
        return float(a * a * b[i] * b[j])

    return f


def channel_tensor_candidate(channel: str, s_value: float, cut_channel: float, im_m: float) -> tuple[np.ndarray, np.ndarray]:
    tensor = np.zeros((3, 3), dtype=float)
    errors = np.zeros((3, 3), dtype=float)
    for i in range(3):
        for j in range(i, 3):
            val, err = integrate.quad(
                tensor_integrand_component(i, j, channel, s_value, cut_channel, im_m),
                0.0,
                1.0,
                epsabs=1e-12,
                epsrel=1e-10,
                limit=100,
            )
            tensor[i, j] = tensor[j, i] = float(val)
            errors[i, j] = errors[j, i] = float(err)
    return tensor, errors


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2004 = load("p2004_s954_strict_cutkosky_gx_loop_amplitude_and_robustness_refresh_witness.json")
    p2016 = load("p2016_s966_strict_cutkosky_channelwise_uncertainty_transport_witness.json")
    rows = p2004.get("grid_table", [])
    if not rows:
        raise RuntimeError("P2004 grid_table missing; run p2004 first")

    tensor_rows: list[dict[str, Any]] = []
    sigma_scan: list[dict[str, Any]] = []
    max_quad_error = 0.0
    psd_all = True
    trace_positive_all = True

    base_delta = np.array([float(r["Delta_opt"]) for r in rows], dtype=float)
    base_l2 = float(la.norm(base_delta, 2))
    channel_vectors = {
        ch: np.array([float(r["Cut_channels"].get(ch, 0.0)) for r in rows], dtype=float)
        for ch in CHANNELS
    }

    # Accumulate quadrature tensor traces per channel.  This replaces a pure
    # variance-placeholder scan with an explicit candidate-quadrature scan, but
    # it is not yet a Feynman-derived backend tensor covariance.
    trace_profiles = {ch: [] for ch in CHANNELS}

    for row in rows:
        s_string = row["s"]
        s_value = float(sp.N(sp.sympify(s_string), 50))
        im_m = float(row["ImM"])
        per_channel: dict[str, Any] = {}
        for ch in CHANNELS:
            cut_val = float(row["Cut_channels"].get(ch, 0.0))
            tensor, errors = channel_tensor_candidate(ch, s_value, cut_val, im_m)
            eigvals = np.linalg.eigvalsh(tensor)
            trace = float(np.trace(tensor))
            max_err = float(np.max(errors))
            max_quad_error = max(max_quad_error, max_err)
            psd = bool(np.min(eigvals) >= -1e-12)
            psd_all = psd_all and psd
            trace_positive_all = trace_positive_all and trace > 0.0
            trace_profiles[ch].append(trace)
            per_channel[ch] = {
                "tensor_3x3": tensor.tolist(),
                "quad_error_3x3": errors.tolist(),
                "trace": trace,
                "eigval_min": float(np.min(eigvals)),
                "eigval_max": float(np.max(eigvals)),
                "psd_with_tolerance": psd,
            }
        tensor_rows.append({"s": s_string, "strict_quadrature_tensor_candidates": per_channel})

    trace_matrix = np.array([trace_profiles[ch] for ch in CHANNELS], dtype=float)
    centered = trace_matrix - trace_matrix.mean(axis=1, keepdims=True)
    c4 = centered @ centered.T / max(centered.shape[1] - 1, 1)
    c4 += 1e-18 * np.eye(len(CHANNELS))
    eig_c4 = np.linalg.eigvalsh(c4)
    sigma = np.sqrt(np.clip(np.diag(c4), 0.0, None))
    sigma = sigma / (float(la.norm(sigma, 2)) + 1e-15)

    counts = {
        "MISSING_CHANNEL_PRESSURE_SUPPORTED": 0,
        "STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED": 0,
        "MIXED_OR_INCONCLUSIVE": 0,
    }
    for t in np.linspace(-1.0, 1.0, 21):
        scales = 1.0 + float(t) * sigma
        transported_delta = base_delta.copy()
        for i, ch in enumerate(CHANNELS):
            transported_delta -= (float(scales[i]) - 1.0) * channel_vectors[ch]
        ratio = float(la.norm(transported_delta, 2) / base_l2) if base_l2 > 0 else float("inf")
        if ratio < 0.95:
            cls = "MISSING_CHANNEL_PRESSURE_SUPPORTED"
        elif ratio > 1.05:
            cls = "STRUCTURAL_OBSTRUCTION_PRESSURE_SUPPORTED"
        else:
            cls = "MIXED_OR_INCONCLUSIVE"
        counts[cls] += 1
        sigma_scan.append({
            "t": float(t),
            "scales": {ch: float(scales[i]) for i, ch in enumerate(CHANNELS)},
            "l2_ratio_vs_p2004": ratio,
            "classifier": cls,
        })

    total = len(sigma_scan)
    dominant = max(counts, key=lambda k: counts[k])
    dominant_share = counts[dominant] / total if total else 0.0

    numerical_gate = {
        "p2004_present": p2004.get("result_kind") == "PASS_GX_LOOP_AND_ROBUSTNESS_REFRESH_WITNESS",
        "p2016_present": p2016.get("result_kind") == "PASS_CHANNELWISE_UNCERTAINTY_TRANSPORT_WITNESS",
        "strict_kernel_parameters_declared": STRICT_PARAMS == {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 1.8},
        "quadrature_rows_exported": len(tensor_rows) == len(rows) and len(rows) > 0,
        "all_channel_tensor_candidates_psd": psd_all,
        "all_channel_tensor_candidate_traces_positive": trace_positive_all,
        "quadrature_error_bounded": max_quad_error < 1e-9,
        "channel_covariance_candidate_exported": c4.shape == (4, 4),
        "channel_covariance_candidate_psd": bool(np.min(eig_c4) >= -1e-15),
        "coupled_scan_nonempty": total > 0,
        "dominant_share_bounded": 0.0 <= dominant_share <= 1.0,
    }
    provenance_gate = {
        "strict_feynman_rule_integrand_exported": False,
        "loop_measure_derived_from_strict_action": False,
        "channel_distance_map_derived_not_assumed": False,
        "normalization_independent_of_p2004_cut_values": False,
        "all_state_phase_space_integral_exported": False,
    }
    gate = {
        **numerical_gate,
        "provenance_gap_declared": not all(provenance_gate.values()),
        "false_pass_blocked_by_provenance_gate": not all(provenance_gate.values()),
    }

    out = {
        "ledger_id": "P2017_S967_STRICT_CUTKOSKY_QUADRATURE_ANSATZ_PROVENANCE_WITNESS",
        "packet_id": "P2017",
        "stage_id": "S967",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TS,
        "route": "strict_only",
        "channel": "graviton->gauge_gauge",
        "depends_on": {
            "p2004": p2004.get("ledger_id"),
            "p2016": p2016.get("ledger_id"),
        },
        "strict_kernel_contract": {
            "formula": "K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)",
            "parameters": STRICT_PARAMS,
        },
        "quadrature_amplitude_candidate_contract": {
            "phase_space_interval": [0.0, 1.0],
            "angular_basis": ["1", "2*x-1", "6*x*(1-x)"],
            "distance_map": "d_channel(s,x)=sqrt(s)*(1+x+offset_channel), offsets={gg:0, gh:1/4, hh:1/2, gx:3/4}",
            "amplitude_candidate": "A_channel(s,x)=sqrt(Cut_channel(s)/ImM(s))*K_strict(d_channel(s,x))*sqrt(x*(1-x))",
            "tensor_candidate": "T_channel_ab(s)=Integral_0^1 A_channel(s,x)^2 b_a(x)b_b(x) dx",
            "provenance_boundary": "candidate ansatz calibrated by P2004 channel cuts; not a Feynman-derived backend amplitude",
        },
        "tensor_candidate_table": tensor_rows,
        "quadrature_channel_covariance_candidate": {
            "trace_profiles_by_channel": {ch: [float(v) for v in vals] for ch, vals in trace_profiles.items()},
            "C4_channel_covariance_from_quadrature_traces": c4.tolist(),
            "C4_eigvals": eig_c4.tolist(),
            "channel_sigma_unit": sigma.tolist(),
            "max_quadrature_error": max_quad_error,
        },
        "coupled_quadrature_tensor_candidate_scan": {
            "count": total,
            "table": sigma_scan,
            "counts": counts,
            "dominant_classifier": dominant,
            "dominant_share": dominant_share,
        },
        "delta_stats": {"l2_base_p2004": base_l2},
        "numerical_gatekeeper_checks": numerical_gate,
        "provenance_gatekeeper_checks": provenance_gate,
        "gatekeeper_checks": gate,
        "diagnostic_result_kind": "PASS_STRICT_KERNEL_QUADRATURE_NUMERICS" if all(numerical_gate.values()) else "OPEN_NUMERICAL_GAP",
        "result_kind": "OPEN_PROVENANCE_GAP_WITH_STRICT_QUADRATURE_TRACE",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "false_pass_guard": "P2017 exports strict-kernel quadrature tensor candidates only. Because the amplitude is calibrated by P2004 cuts and lacks a strict Feynman-rule integrand derivation, it must not be called backend-derived Cutkosky closure or ToE closure.",
        "next_honest_step": "Build P2018 from explicit strict-side Feynman/loop-integrand rules: derive the channel distance map, loop measure, normalization, and tensor integrand independently of P2004 cut calibration, then compare against P2017 candidates.",
        "toe_progress": "Improves the unitarity block by preserving useful quadrature numerics while adding a provenance gate that prevents a false ToE or Cutkosky pass.",
        "lay_explanation": "Policzone całki są użytecznym testem numerycznym, ale nie wolno ich jeszcze traktować jako pełnego wyprowadzenia fizycznego, bo wzór amplitudy nadal jest kandydatem skalibrowanym wcześniejszymi danymi.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": __import__("scipy").__version__, "sympy": sp.__version__},
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P2017] wrote witness: {OUT}")


if __name__ == "__main__":
    main()
