#!/usr/bin/env python3
"""
QW-2090: H0-Lambda decoupling gate (deterministic, no-retune, no-scan).

Goal:
- attempt an internally consistent decoupling of (h0, lambda_cosmological)
  from an independent H(z) node system,
- block false closure claims when only coupled/anchored relations are available.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "h0_lambda_decoupling_input_qw2090.json"
OUT_JSON = ROOT / "report_qw2090_h0_lambda_decoupling_gate.json"
OUT_MD = ROOT / "RAPORT_QW2090_H0_LAMBDA_DECOUPLING_GATE.md"


def load_json(path: Path) -> Dict:
    return json.loads(path.read_text(encoding="utf-8"))


def get_ref(groups: Dict, pid: str) -> Dict:
    for _, items in groups.items():
        for item in items:
            if item.get("id") == pid:
                return item
    raise KeyError(f"Missing registry parameter: {pid}")


def rel_err_pct(pred: float, ref: float) -> float:
    if ref == 0.0:
        return float("inf") if pred != 0.0 else 0.0
    return abs(pred - ref) / abs(ref) * 100.0


def source_metadata_complete(in_data: Dict) -> bool:
    src = str(in_data.get("source", "")).strip()
    citation = str(in_data.get("citation", "")).strip()
    ref_url = str(in_data.get("reference_url", "")).strip()
    src_sha = str(in_data.get("source_sha256", "")).strip()
    src_ver = str(in_data.get("source_version", "")).strip()
    low = f"{src} {citation} {ref_url}".lower()
    not_placeholder = bool(src) and ("placeholder" not in low)
    has_reference = bool(citation or ref_url)
    has_integrity = bool(src_sha or src_ver)
    return bool(not_placeholder and has_reference and has_integrity)


def solve_linear_decoupling(
    nodes: List[Dict], omega_m: float, omega_r: float
) -> Tuple[float, float, float]:
    """
    Solve weighted linear system:
      H(z)^2 = x * E(z) + y,
      E(z) = Omega_m(1+z)^3 + Omega_r(1+z)^4
    Then:
      H0 = sqrt(x), Omega_Lambda = y/x.
    Returns: (h0_km_s_mpc, omega_lambda, max_rel_residual).
    """
    a_rows = []
    b_rows = []
    w_rows = []
    for row in nodes:
        z = float(row["z"])
        h = float(row["h_km_s_mpc"])
        s = abs(float(row["sigma_total"]))
        if h <= 0.0 or s <= 0.0:
            continue
        e = omega_m * (1.0 + z) ** 3 + omega_r * (1.0 + z) ** 4
        # Var(H^2) ~ (2 H sigma)^2
        var_h2 = max((2.0 * h * s) ** 2, 1e-12)
        w = 1.0 / var_h2
        a_rows.append([e, 1.0])
        b_rows.append(h * h)
        w_rows.append(w)

    if len(a_rows) < 2:
        raise RuntimeError("Need at least 2 valid H(z) nodes for decoupling.")

    a = np.array(a_rows, dtype=float)
    b = np.array(b_rows, dtype=float)
    w = np.sqrt(np.array(w_rows, dtype=float))
    aw = a * w[:, None]
    bw = b * w
    x, y = np.linalg.lstsq(aw, bw, rcond=None)[0]
    rank = int(np.linalg.matrix_rank(aw))
    if rank < 2:
        raise RuntimeError("Decoupling regression is rank-deficient.")
    if x <= 0.0:
        raise RuntimeError("Unphysical decoupling solution: H0^2 <= 0.")

    h0 = math.sqrt(float(x))
    omega_lambda = float(y / x)

    # Relative residual on H(z).
    h_pred = np.sqrt(np.clip(x * a[:, 0] + y, 0.0, None))
    rel = np.abs(h_pred - np.sqrt(np.clip(b, 0.0, None))) / np.clip(
        np.sqrt(np.clip(b, 0.0, None)), 1e-12, None
    )
    max_rel_residual = float(np.max(rel))
    return h0, omega_lambda, max_rel_residual


def solve_flatness_projection(
    nodes: List[Dict], omega_m: float, omega_r: float, c_si: float, mpc_m: float
) -> Tuple[float, float, float, float]:
    """
    Deterministic diagnostic projection under flatness prior:
      Omega_Lambda = 1 - Omega_m - Omega_r.
    Returns:
      (h0_km_s_mpc, lambda_m_inv2, omega_lambda_flat, max_rel_residual_hz)
    """
    omega_lambda_flat = float(1.0 - omega_m - omega_r)
    if omega_lambda_flat <= 0.0:
        raise RuntimeError("Unphysical flatness projection: Omega_Lambda <= 0.")

    h0_vals = []
    h0_w = []
    hz_rel = []
    for row in nodes:
        z = float(row["z"])
        h = float(row["h_km_s_mpc"])
        s = abs(float(row["sigma_total"]))
        if h <= 0.0 or s <= 0.0:
            continue
        e = omega_m * (1.0 + z) ** 3 + omega_r * (1.0 + z) ** 4 + omega_lambda_flat
        if e <= 0.0:
            continue
        h0_i = h / math.sqrt(e)
        sig_h0_i = s / math.sqrt(e)
        w_i = 1.0 / max(sig_h0_i * sig_h0_i, 1e-12)
        h0_vals.append(h0_i)
        h0_w.append(w_i)

    if len(h0_vals) < 2:
        raise RuntimeError("Need at least 2 valid H(z) nodes for flatness projection.")

    h0 = float(np.average(np.array(h0_vals, dtype=float), weights=np.array(h0_w, dtype=float)))
    h0_si = h0 * 1000.0 / mpc_m
    lam = float(3.0 * omega_lambda_flat * h0_si * h0_si / (c_si * c_si))

    for row in nodes:
        z = float(row["z"])
        h = float(row["h_km_s_mpc"])
        if h <= 0.0:
            continue
        e = omega_m * (1.0 + z) ** 3 + omega_r * (1.0 + z) ** 4 + omega_lambda_flat
        h_pred = h0 * math.sqrt(max(e, 0.0))
        hz_rel.append(abs(h_pred - h) / max(abs(h), 1e-12))
    max_rel = float(max(hz_rel)) if hz_rel else float("inf")

    return h0, lam, omega_lambda_flat, max_rel


def main() -> None:
    p = argparse.ArgumentParser(description="QW-2090 H0-Lambda decoupling gate")
    p.add_argument("--input", default=str(DEFAULT_INPUT), help="Decoupling input JSON.")
    args = p.parse_args()

    reg = load_json(ROOT / "report_qw2068_sm_gr_parameter_registry.json")
    r2073 = load_json(ROOT / "report_qw2073_radiative_channels_closure_upgrade.json")
    groups = reg["groups"]

    h0_ref_item = get_ref(groups, "h0")
    lam_ref_item = get_ref(groups, "lambda_cosmological")
    h0_ref = float(h0_ref_item["value"])
    h0_tol = float(h0_ref_item["tolerance_rel_pct"])
    lam_ref = float(lam_ref_item["value"])
    lam_tol = float(lam_ref_item["tolerance_rel_pct"])

    c_si = float(get_ref(groups, "c_light")["value"])
    mpc_m = 3.085677581e22

    in_path = Path(args.input).resolve()
    in_data = {}
    input_present = False
    if in_path.exists():
        in_data = load_json(in_path)
        input_present = True

    source_label = "no_external_input"
    strict_candidate_available = False
    fit_failed_reason = ""
    h0_candidate = None
    lam_candidate = None
    omega_lambda_candidate = None
    max_rel_residual = None
    identifiability_cond_number = None
    identifiability_e_span = None

    flat_h0_candidate = None
    flat_lam_candidate = None
    flat_omega_lambda = None
    flat_max_rel_residual = None
    flat_projection_failed_reason = ""

    if input_present:
        nodes = in_data.get("h_of_z_nodes", [])
        omega_m = float(in_data.get("omega_m", 0.315))
        omega_r = float(in_data.get("omega_r", 9.0e-5))
        source_label = str(in_data.get("source", "input_json"))

        # Identifiability diagnostics for two-parameter decoupling matrix [E(z), 1].
        valid_nodes = []
        for row in nodes:
            z = float(row["z"])
            h = float(row["h_km_s_mpc"])
            s = abs(float(row["sigma_total"]))
            if h <= 0.0 or s <= 0.0:
                continue
            e = omega_m * (1.0 + z) ** 3 + omega_r * (1.0 + z) ** 4
            valid_nodes.append((e, h, s))
        if len(valid_nodes) >= 2:
            e_arr = np.array([r[0] for r in valid_nodes], dtype=float)
            a = np.array([[r[0], 1.0] for r in valid_nodes], dtype=float)
            identifiability_cond_number = float(np.linalg.cond(a))
            identifiability_e_span = float(np.max(e_arr) - np.min(e_arr))

        try:
            h0_fit, omega_lambda_fit, max_rel = solve_linear_decoupling(nodes, omega_m, omega_r)
            h0_candidate = float(h0_fit)
            omega_lambda_candidate = float(omega_lambda_fit)
            h0_si = h0_candidate * 1000.0 / mpc_m
            lam_candidate = float(3.0 * omega_lambda_candidate * h0_si * h0_si / (c_si * c_si))
            max_rel_residual = float(max_rel)
            strict_candidate_available = True
        except Exception as exc:  # noqa: BLE001
            fit_failed_reason = str(exc)

        try:
            flat_h0, flat_lam, flat_ol, flat_res = solve_flatness_projection(
                nodes, omega_m, omega_r, c_si, mpc_m
            )
            flat_h0_candidate = float(flat_h0)
            flat_lam_candidate = float(flat_lam)
            flat_omega_lambda = float(flat_ol)
            flat_max_rel_residual = float(flat_res)
        except Exception as exc:  # noqa: BLE001
            flat_projection_failed_reason = str(exc)

    # Deterministic fallback (explicitly coupled/anchor-dependent).
    if not strict_candidate_available:
        omega_lambda_flat = float(r2073["diagnostics"]["omega_lambda_flat"])
        h0_si = math.sqrt(lam_ref * c_si * c_si / (3.0 * omega_lambda_flat))
        h0_candidate = float(h0_si * (mpc_m / 1000.0))
        lam_candidate = float(lam_ref)
        omega_lambda_candidate = float(omega_lambda_flat)
        max_rel_residual = None
        if not fit_failed_reason:
            fit_failed_reason = "independent_hz_nodes_not_available"

    h0_rel = rel_err_pct(h0_candidate, h0_ref)
    lam_rel = rel_err_pct(lam_candidate, lam_ref)
    flat_h0_rel = (
        rel_err_pct(float(flat_h0_candidate), h0_ref) if flat_h0_candidate is not None else None
    )
    flat_lam_rel = (
        rel_err_pct(float(flat_lam_candidate), lam_ref) if flat_lam_candidate is not None else None
    )

    flags = {
        "independent_input_present": bool(input_present),
        "source_metadata_complete": bool(source_metadata_complete(in_data) if input_present else False),
        "strict_candidate_available": bool(strict_candidate_available),
        "solution_finite": bool(
            np.isfinite(h0_candidate) and np.isfinite(lam_candidate) and np.isfinite(omega_lambda_candidate)
        ),
        "omega_lambda_physical_window": bool(0.0 < omega_lambda_candidate < 1.5),
        "fit_residual_lt_0p25": bool(
            max_rel_residual is None or float(max_rel_residual) <= 0.25
        ),
        "h0_within_registry_tolerance": bool(h0_rel <= h0_tol),
        "lambda_within_registry_tolerance": bool(lam_rel <= lam_tol),
        "no_anchor_feedback_loop": bool(in_data.get("provenance_anchor_free", False)),
    }
    pass_count = int(sum(1 for v in flags.values() if bool(v)))
    total_flags = int(len(flags))

    strict_candidate_governance_ready = bool(
        flags["independent_input_present"]
        and flags["source_metadata_complete"]
        and flags["strict_candidate_available"]
        and flags["solution_finite"]
        and flags["omega_lambda_physical_window"]
        and flags["fit_residual_lt_0p25"]
        and flags["no_anchor_feedback_loop"]
    )
    strict_targets_within_tol = bool(
        flags["h0_within_registry_tolerance"] and flags["lambda_within_registry_tolerance"]
    )
    strict_pass = bool(strict_candidate_governance_ready and strict_targets_within_tol)
    strict_target_miss = bool(strict_candidate_governance_ready and not strict_targets_within_tol)

    if strict_pass:
        verdict = "H0_LAMBDA_DECOUPLING_GATE_PASS_STRICT"
        updates = [
            {
                "id": "h0",
                "predicted_value": float(h0_candidate),
                "status": "derived",
                "strict_level": "strict_internal_gate",
                "method": "qw2090_weighted_hz_decoupling_two_parameter_fit_nonanchor",
                "notes": "Decoupled from independent H(z) node system under strict no-retune constraints.",
            },
            {
                "id": "lambda_cosmological",
                "predicted_value": float(lam_candidate),
                "status": "derived",
                "strict_level": "strict_internal_gate",
                "method": "qw2090_weighted_hz_decoupling_two_parameter_fit_nonanchor",
                "notes": "Derived jointly with H0 from independent H(z) node system.",
            },
        ]
    elif strict_target_miss:
        verdict = "H0_LAMBDA_DECOUPLING_GATE_TARGET_MISS"
        h0_status = "derived" if flags["h0_within_registry_tolerance"] else "derived_strict_target_miss"
        lam_status = (
            "derived" if flags["lambda_within_registry_tolerance"] else "derived_strict_target_miss"
        )
        updates = [
            {
                "id": "h0",
                "predicted_value": float(h0_candidate),
                "status": h0_status,
                "strict_level": "strict_internal_gate",
                "method": "qw2090_weighted_hz_decoupling_two_parameter_fit_nonanchor",
                "notes": "Strict decoupling candidate exists; target miss vs registry tolerance.",
            },
            {
                "id": "lambda_cosmological",
                "predicted_value": float(lam_candidate),
                "status": lam_status,
                "strict_level": "strict_internal_gate",
                "method": "qw2090_weighted_hz_decoupling_two_parameter_fit_nonanchor",
                "notes": "Strict decoupling candidate exists; target miss vs registry tolerance.",
            },
        ]
    else:
        verdict = "H0_LAMBDA_DECOUPLING_GATE_PENDING_NONCLOSING"
        updates = [
            {
                "id": "h0",
                "predicted_value": float(h0_candidate),
                "status": "derived_coupled_anchor_dependent",
                "strict_level": "coupled_anchor_dependent",
                "method": "qw2090_coupled_fallback_from_lambda_omega_relation",
                "notes": "Decoupling not strict: independent anchor-free H(z) constraints are missing or insufficient.",
            },
            {
                "id": "lambda_cosmological",
                "predicted_value": float(lam_candidate),
                "status": "derived_coupled_anchor_dependent",
                "strict_level": "coupled_anchor_dependent",
                "method": "qw2090_coupled_fallback_from_lambda_omega_relation",
                "notes": "Absolute closure remains coupled without independent decoupling observables.",
            },
        ]

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            "registry": "report_qw2068_sm_gr_parameter_registry.json",
            "radiative_upgrade": "report_qw2073_radiative_channels_closure_upgrade.json",
            "input_json": str(in_path) if input_present else None,
        },
        "inputs": {
            "input_present": input_present,
            "source_label": source_label,
        },
        "fit": {
            "fit_failed_reason": fit_failed_reason,
            "h0_candidate_km_s_mpc": h0_candidate,
            "lambda_candidate_m_inv2": lam_candidate,
            "omega_lambda_candidate": omega_lambda_candidate,
            "max_rel_residual": max_rel_residual,
        },
        "reference": {
            "h0_ref": h0_ref,
            "h0_tol_rel_pct": h0_tol,
            "lambda_ref": lam_ref,
            "lambda_tol_rel_pct": lam_tol,
            "h0_rel_err_pct": h0_rel,
            "lambda_rel_err_pct": lam_rel,
        },
        "identifiability": {
            "design_condition_number": identifiability_cond_number,
            "e_column_span": identifiability_e_span,
            "note": "For [E(z),1] design; small E-span implies weak separation of two-parameter decoupling.",
        },
        "flatness_projection_diagnostic": {
            "flat_projection_failed_reason": flat_projection_failed_reason,
            "h0_candidate_km_s_mpc": flat_h0_candidate,
            "lambda_candidate_m_inv2": flat_lam_candidate,
            "omega_lambda_flat": flat_omega_lambda,
            "max_rel_residual": flat_max_rel_residual,
            "h0_rel_err_pct": flat_h0_rel,
            "lambda_rel_err_pct": flat_lam_rel,
            "h0_within_registry_tolerance": (
                bool(flat_h0_rel <= h0_tol) if flat_h0_rel is not None else False
            ),
            "lambda_within_registry_tolerance": (
                bool(flat_lam_rel <= lam_tol) if flat_lam_rel is not None else False
            ),
        },
        "flags": flags,
        "pass_count": pass_count,
        "total_flags": total_flags,
        "updates": updates,
        "verdict": verdict,
        "required_next_step": (
            "INTEGRATE_QW2090_IN_QW2069_AND_RERUN_FULL_CLOSURE_CHAIN"
            if strict_pass
            else (
                "REFINE_NONANCHOR_HZ_NODE_SYSTEM_OR_EXTEND_MODEL_BEFORE_STRICT_RETEST"
                if strict_target_miss
                else "PROVIDE_INDEPENDENT_ANCHOR_FREE_HZ_NODE_INPUT_AND_RERUN_QW2090"
            )
        ),
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-2090: H0 LAMBDA DECOUPLING GATE",
        "",
        f"- Date UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- pass_count: {pass_count}/{total_flags}",
        f"- input_present: {input_present}",
        f"- fit_failed_reason: {fit_failed_reason}",
        "",
        "## Candidate",
        f"- H0: {h0_candidate:.9f} km s^-1 Mpc^-1",
        f"- Lambda: {lam_candidate:.12e} m^-2",
        f"- Omega_lambda: {omega_lambda_candidate:.9f}",
        f"- max_rel_residual: {max_rel_residual}",
        "",
        "## Registry Error",
        f"- H0 rel_err_pct: {h0_rel:.6f} (tol={h0_tol:.3f})",
        f"- Lambda rel_err_pct: {lam_rel:.6f} (tol={lam_tol:.3f})",
        "",
        "## Identifiability + Flatness Diagnostic",
        f"- design_condition_number([E,1]): {identifiability_cond_number}",
        f"- E-column span: {identifiability_e_span}",
        f"- flat_h0_candidate: {flat_h0_candidate}",
        f"- flat_lambda_candidate: {flat_lam_candidate}",
        f"- flat_h0_rel_err_pct: {flat_h0_rel}",
        f"- flat_lambda_rel_err_pct: {flat_lam_rel}",
        "",
        "## Artifact",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-2090] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-2090] Saved MD:   {OUT_MD.name}")
    print(f"[QW-2090] verdict={verdict} pass_count={pass_count}/{total_flags}")


if __name__ == "__main__":
    main()
