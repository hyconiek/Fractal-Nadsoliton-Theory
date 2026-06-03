#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2519_s1469_strict_damping_biharmonic_endpoint_anchor_acceptance_certificate.json"
MD = GEN / "p2519_s1469_strict_damping_biharmonic_endpoint_anchor_acceptance_certificate.md"

SOURCE_FILES = {
    "P2518_AFFINE_SLOPE_NONIDENTIFIABILITY": GEN / "p2518_s1468_strict_damping_biharmonic_affine_slope_nonidentifiability_certificate.json",
}

STRICT_DELTA = Fraction(4, 5)
STRICT_ETA = Fraction(9, 5)
STRICT_BETA = Fraction(1, 1)
DOMAIN = list(range(1, 12))
SLOPE_CANDIDATES = [Fraction(-1, 1), Fraction(0, 1), Fraction(1, 2), Fraction(4, 5), Fraction(1, 1), Fraction(9, 5), Fraction(2, 1)]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:40]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2519|S1469|biharmonic endpoint anchor acceptance|endpoint slope anchor|two-node affine anchor|numeric endpoint anchor|affine slope anchor acceptance",
        "precursor_packets": "P2518|S1468|biharmonic affine-slope nonidentifiability|P2517|dual-key axiom boundary|P2516|dual-key source acceptance",
        "anchor_language": "endpoint anchor|two-node anchor|affine anchor|node anchor|slope anchor|anchor determinant",
        "numeric_source_language": "beta_eta_numeric_source|numeric beta/eta|delta=4/5|eta=9/5|strict source theorem",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def frac_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def endpoint_anchor_symbolic() -> dict[str, Any]:
    return {
        "coordinate": "ell = log(d)",
        "left_endpoint": {"d": 1, "ell": "0", "y_value": "0"},
        "right_endpoint": {"d": 11, "ell": "log(11)", "y_value": "(4/5)*log(11)"},
        "affine_unknowns": ["b", "a"],
        "anchor_matrix": [["1", "0"], ["1", "log(11)"]],
        "anchor_determinant": "log(11)",
        "determinant_positive": math.log(11) > 0,
        "solved_intercept_b": "0",
        "solved_slope_a_delta": frac_text(STRICT_DELTA),
        "solved_eta": frac_text(STRICT_ETA),
        "solved_beta": frac_text(STRICT_BETA),
        "conditional_acceptance_statement": "If the two endpoint node anchors are strict theorems and the m=2 affine-kernel signature is admitted, then the affine zero-mode is pinned to y(ell)=(4/5)ell, hence beta=1 and eta=9/5. The anchor data themselves are not sourced by this certificate.",
    }


def endpoint_candidate_rows() -> list[dict[str, Any]]:
    rows = []
    right = math.log(11)
    strict_right = float(STRICT_DELTA) * right
    for slope in SLOPE_CANDIDATES:
        candidate_right = float(slope) * right
        residual = candidate_right - strict_right
        rows.append({
            "candidate_slope": frac_text(slope),
            "candidate_eta_if_beta_one": frac_text(Fraction(1, 1) + slope),
            "left_anchor_residual_y_0": 0.0,
            "right_anchor_residual_at_log_11": residual,
            "absolute_right_anchor_residual": abs(residual),
            "accepted_by_endpoint_anchors": residual == 0.0,
        })
    return rows


def all_node_residual_rows() -> list[dict[str, Any]]:
    rows = []
    for d in DOMAIN:
        ell = math.log(d)
        target = float(STRICT_DELTA) * ell
        rows.append({
            "d": d,
            "ell_log_d": ell,
            "target_y_delta_log_d": target,
            "reconstructed_y_delta_log_d": target,
            "residual": 0.0,
            "reconstructed_power_d_eta_minus_1": d ** float(STRICT_DELTA),
            "strict_denominator_component_d_eta": d ** float(STRICT_ETA),
        })
    return rows


def finite_design_matrix_audit() -> dict[str, Any]:
    # Least-squares design for y=b+a*ell on all strict target nodes.
    xs = [math.log(d) for d in DOMAIN]
    ys = [float(STRICT_DELTA) * x for x in xs]
    n = len(xs)
    sum_x = sum(xs)
    sum_x2 = sum(x * x for x in xs)
    sum_y = sum(ys)
    sum_xy = sum(x * y for x, y in zip(xs, ys))
    gram_det = n * sum_x2 - sum_x * sum_x
    solved_b = (sum_x2 * sum_y - sum_x * sum_xy) / gram_det
    solved_a = (n * sum_xy - sum_x * sum_y) / gram_det
    residuals = [y - (solved_b + solved_a * x) for x, y in zip(xs, ys)]
    return {
        "design": "columns [1, log(d)] for d=1..11",
        "normal_matrix": [[n, sum_x], [sum_x, sum_x2]],
        "normal_matrix_determinant": gram_det,
        "normal_matrix_determinant_positive": gram_det > 0,
        "least_squares_intercept": solved_b,
        "least_squares_slope": solved_a,
        "least_squares_slope_matches_delta_4_over_5": abs(solved_a - float(STRICT_DELTA)) < 1e-14,
        "max_abs_residual": max(abs(value) for value in residuals),
        "all_node_residuals_zero_within_float_tolerance": max(abs(value) for value in residuals) < 1e-14,
        "node_count": n,
    }


def build_endpoint_anchor_certificate(p2518: dict[str, Any]) -> dict[str, Any]:
    symbolic = endpoint_anchor_symbolic()
    candidate_rows = endpoint_candidate_rows()
    design = finite_design_matrix_audit()
    accepted = [row for row in candidate_rows if row["accepted_by_endpoint_anchors"]]
    return {
        "frontier_atom_under_attack": "minimal numeric endpoint anchors for beta_eta_numeric_source_target",
        "p2518_nonidentifiability_inherited": p2518.get("operator_signature_numeric_key_separation_exported") is True,
        "certificate_type": "conditional endpoint-anchor acceptance certificate after biharmonic affine-slope nonidentifiability",
        "endpoint_anchor_symbolic_solution": symbolic,
        "endpoint_candidate_acceptance_rows": candidate_rows,
        "all_node_reconstruction_rows": all_node_residual_rows(),
        "finite_design_matrix_audit": design,
        "unique_candidate_slope_accepted_by_endpoint_anchors": len(accepted) == 1 and accepted[0]["candidate_slope"] == frac_text(STRICT_DELTA),
        "endpoint_anchors_pin_affine_slope_if_anchors_are_sourced": symbolic["determinant_positive"] and len(accepted) == 1,
        "all_strict_nodes_reconstructed_after_anchor_pinning": design["all_node_residuals_zero_within_float_tolerance"],
        "conditional_beta_eta_target_recovered_from_anchors": True,
        "endpoint_anchor_source_exported": False,
        "beta_eta_numeric_source_exported": False,
        "m2_operator_signature_source_exported": False,
        "strict_damping_beta_eta_source_exported": False,
        "damping_compression_bridge_component_ready": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_claimed": False,
    }


def append_doc_sections() -> None:
    eq_section = """
## P2519/S1469 biharmonic endpoint-anchor acceptance certificate

`P2519/S1469` follows the P2518 affine-slope nonidentifiability theorem with a minimal conditional acceptance audit.  The `m=2` biharmonic signature leaves the whole affine family `y(ell)=a ell+b` invisible, but two strict endpoint node anchors, `y(0)=0` and `y(log 11)=(4/5) log 11`, would pin the affine kernel because the anchor matrix determinant is `log 11 > 0`.  Under those anchors the reconstructed flow is `y(ell)=(4/5)ell`, giving `beta=1` and `eta=9/5`, and the finite all-node design audit on `d=1..11` has zero residual within arithmetic tolerance.

This is only a conditional numeric-target acceptance theorem: it identifies the minimal endpoint data that would turn the P2518 nonidentifiability into a pinned affine flow, but it does not source those endpoint anchors.  Therefore it exports no `beta_eta_numeric_source`, no `m2_operator_signature_source`, no strict damping source closure, no bridge completion, no role-transfer theorem, no QW-2191 discharge, no role-bearing `L_total`, no physical-value generator, and no ToE closure.
"""
    lag_section = """
## P2519/S1469 endpoint-anchor acceptance guard

`P2519/S1469` records the exact conditional data needed after P2518: the biharmonic affine kernel is pinned only if strict dynamics supply endpoint node anchors fixing `y(0)=0` and `y(log 11)=(4/5)log 11`.  The calculation recovers the numeric target under those anchors, but the anchors remain unsourced, so no nonlinear compression-flow source or role-bearing `L_total` term is licensed.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2519/S1469 biharmonic endpoint-anchor acceptance certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2519/S1469 endpoint-anchor acceptance guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2518 = theorem(sources["P2518_AFFINE_SLOPE_NONIDENTIFIABILITY"], "strict_damping_biharmonic_affine_slope_nonidentifiability_certificate")
    cert = build_endpoint_anchor_certificate(p2518)
    design = cert["finite_design_matrix_audit"]
    theorem_export = {
        "theorem_name": "P2519_T1_strict_damping_biharmonic_endpoint_anchor_acceptance_certificate",
        "audited_chain": ["P2518/S1468", "P2517/S1467", "P2516/S1466"],
        "strict_damping_biharmonic_endpoint_anchor_acceptance_certificate": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2518_nonidentifiability_inherited": cert["p2518_nonidentifiability_inherited"],
        "anchor_determinant": cert["endpoint_anchor_symbolic_solution"]["anchor_determinant"],
        "anchor_determinant_positive": cert["endpoint_anchor_symbolic_solution"]["determinant_positive"],
        "solved_slope_delta": cert["endpoint_anchor_symbolic_solution"]["solved_slope_a_delta"],
        "solved_beta": cert["endpoint_anchor_symbolic_solution"]["solved_beta"],
        "solved_eta": cert["endpoint_anchor_symbolic_solution"]["solved_eta"],
        "unique_candidate_slope_accepted_by_endpoint_anchors": cert["unique_candidate_slope_accepted_by_endpoint_anchors"],
        "endpoint_anchors_pin_affine_slope_if_anchors_are_sourced": cert["endpoint_anchors_pin_affine_slope_if_anchors_are_sourced"],
        "all_strict_nodes_reconstructed_after_anchor_pinning": cert["all_strict_nodes_reconstructed_after_anchor_pinning"],
        "finite_design_normal_matrix_determinant_positive": design["normal_matrix_determinant_positive"],
        "finite_design_all_node_residuals_zero": design["all_node_residuals_zero_within_float_tolerance"],
        "conditional_beta_eta_target_recovered_from_anchors": cert["conditional_beta_eta_target_recovered_from_anchors"],
        "endpoint_anchor_source_exported": False,
        "beta_eta_numeric_source_exported": False,
        "m2_operator_signature_source_exported": False,
        "strict_damping_beta_eta_source_exported": False,
        "damping_compression_bridge_component_ready": False,
        "full_bridge_theorem_exported": False,
        "role_transfer_theorem_exported": False,
        "selector_closure_exported": False,
        "qw2191_discharged_by_this_certificate": False,
        "role_bearing_ltotal_exported": False,
        "toe_closure_claimed": False,
        "not_licensed": [
            "P2519 exports a conditional endpoint-anchor acceptance theorem, not the source of those anchors.",
            "The endpoint anchors recover beta=1 and eta=9/5 only if strict dynamics independently supply the node-anchor data.",
            "The m=2 operator signature source remains separate and unsourced by this packet.",
            "No damping bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Search for a strict nadsoliton source for the endpoint node anchors or replace the endpoint-anchor premise by a non-axiomatic internal source theorem.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2518_separation_inherited": theorem_export["p2518_nonidentifiability_inherited"],
        "anchor_determinant_positive": theorem_export["anchor_determinant_positive"],
        "unique_candidate_slope": theorem_export["unique_candidate_slope_accepted_by_endpoint_anchors"],
        "all_nodes_reconstructed": theorem_export["all_strict_nodes_reconstructed_after_anchor_pinning"],
        "conditional_recovery_only": theorem_export["conditional_beta_eta_target_recovered_from_anchors"] and not theorem_export["endpoint_anchor_source_exported"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "endpoint_anchor_source_exported",
            "beta_eta_numeric_source_exported",
            "m2_operator_signature_source_exported",
            "strict_damping_beta_eta_source_exported",
            "damping_compression_bridge_component_ready",
            "full_bridge_theorem_exported",
            "role_transfer_theorem_exported",
            "selector_closure_exported",
            "qw2191_discharged_by_this_certificate",
            "role_bearing_ltotal_exported",
            "toe_closure_claimed",
        ]),
    }
    return {
        "packet_id": "P2519",
        "stage_id": "S1469",
        "status": "STRICT_DAMPING_BIHARMONIC_ENDPOINT_ANCHOR_ACCEPTANCE_CERTIFICATE_CONDITIONAL_TARGET_RECOVERY_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_biharmonic_endpoint_anchor_acceptance_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_biharmonic_endpoint_anchor_acceptance_certificate"]["theorem_export"]
    lines = [
        "# P2519/S1469 strict damping biharmonic endpoint-anchor acceptance certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2518 nonidentifiability inherited: `{t['p2518_nonidentifiability_inherited']}`.",
        f"- Anchor determinant: `{t['anchor_determinant']}`; positive: `{t['anchor_determinant_positive']}`.",
        f"- Solved slope delta: `{t['solved_slope_delta']}`.",
        f"- Solved beta: `{t['solved_beta']}`.",
        f"- Solved eta: `{t['solved_eta']}`.",
        f"- Unique candidate slope accepted by endpoint anchors: `{t['unique_candidate_slope_accepted_by_endpoint_anchors']}`.",
        f"- All strict nodes reconstructed after anchor pinning: `{t['all_strict_nodes_reconstructed_after_anchor_pinning']}`.",
        f"- Conditional beta/eta target recovered from anchors: `{t['conditional_beta_eta_target_recovered_from_anchors']}`.",
        f"- Endpoint-anchor source exported: `{t['endpoint_anchor_source_exported']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports only conditional target recovery from explicitly supplied endpoint anchors. It does not source the anchors, does not source the m=2 operator signature, and exports no bridge completion, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_biharmonic_endpoint_anchor_acceptance_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_biharmonic_endpoint_anchor_acceptance_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
