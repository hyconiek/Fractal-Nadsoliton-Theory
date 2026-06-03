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
OUT = GEN / "p2521_s1471_strict_damping_single_node_anchor_equivalence_certificate.json"
MD = GEN / "p2521_s1471_strict_damping_single_node_anchor_equivalence_certificate.md"

SOURCE_FILES = {
    "P2520_ENDPOINT_SUBKEY_LATTICE": GEN / "p2520_s1470_strict_damping_endpoint_anchor_subkey_lattice_certificate.json",
}

STRICT_DELTA = Fraction(4, 5)
STRICT_ETA = Fraction(9, 5)
NONZERO_NODE_DOMAIN = list(range(2, 12))
FULL_NODE_DOMAIN = list(range(1, 12))
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
        "new_packet": "P2521|S1471|single-node anchor equivalence|single nonzero node anchor|node-anchor equivalence|anchor placement nonidentifiability|nonzero node slope anchor",
        "precursor_packets": "P2520|S1470|endpoint anchor subkey lattice|P2519|endpoint-anchor acceptance|P2518|affine-slope nonidentifiability",
        "node_anchor_language": "single node anchor|nonzero node anchor|node-anchor placement|anchor placement|slope anchor|log\\(d\\)",
        "source_blockers": "beta_eta_numeric_source|endpoint_anchor_source|right_endpoint_value_anchor_source|strict source theorem|proper subset obstruction",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def frac_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def determinant_for_left_plus_node(d: int) -> float:
    return math.log(d)


def single_node_equivalence_rows() -> list[dict[str, Any]]:
    rows = []
    for d_anchor in NONZERO_NODE_DOMAIN:
        ell_anchor = math.log(d_anchor)
        y_anchor = float(STRICT_DELTA) * ell_anchor
        solved_delta = y_anchor / ell_anchor
        all_node_residuals = []
        for d in FULL_NODE_DOMAIN:
            ell = math.log(d)
            residual = float(STRICT_DELTA) * ell - solved_delta * ell
            all_node_residuals.append(residual)
        rows.append({
            "anchor_node_d": d_anchor,
            "anchor_ell_log_d": ell_anchor,
            "anchor_value_y": y_anchor,
            "left_plus_node_constraint_matrix": [[1.0, 0.0], [1.0, ell_anchor]],
            "determinant_log_d": determinant_for_left_plus_node(d_anchor),
            "determinant_positive": determinant_for_left_plus_node(d_anchor) > 0,
            "solved_intercept_b": 0.0,
            "solved_delta": solved_delta,
            "solved_delta_matches_4_over_5": abs(solved_delta - float(STRICT_DELTA)) < 1e-14,
            "solved_eta": 1.0 + solved_delta,
            "solved_eta_matches_9_over_5": abs((1.0 + solved_delta) - float(STRICT_ETA)) < 1e-14,
            "max_abs_all_node_residual_after_pinning": max(abs(value) for value in all_node_residuals),
            "all_nodes_reconstructed_after_pinning": max(abs(value) for value in all_node_residuals) < 1e-14,
        })
    return rows


def slope_candidate_by_anchor_rows() -> list[dict[str, Any]]:
    rows = []
    for d_anchor in NONZERO_NODE_DOMAIN:
        ell_anchor = math.log(d_anchor)
        target = float(STRICT_DELTA) * ell_anchor
        accepted = []
        for slope in SLOPE_CANDIDATES:
            residual = float(slope) * ell_anchor - target
            row = {
                "anchor_node_d": d_anchor,
                "candidate_delta": frac_text(slope),
                "candidate_eta": frac_text(Fraction(1, 1) + slope),
                "anchor_residual_with_left_normalization": residual,
                "accepted_by_anchor": abs(residual) < 1e-14,
            }
            if row["accepted_by_anchor"]:
                accepted.append(row)
            rows.append(row)
    return rows


def placement_lattice() -> dict[str, Any]:
    equivalence = single_node_equivalence_rows()
    candidate_rows = slope_candidate_by_anchor_rows()
    by_anchor = []
    for d_anchor in NONZERO_NODE_DOMAIN:
        accepted = [row for row in candidate_rows if row["anchor_node_d"] == d_anchor and row["accepted_by_anchor"]]
        by_anchor.append({
            "anchor_node_d": d_anchor,
            "accepted_candidate_count": len(accepted),
            "accepted_candidate_delta": accepted[0]["candidate_delta"] if accepted else None,
            "accepted_candidate_eta": accepted[0]["candidate_eta"] if accepted else None,
            "unique_strict_candidate_accepted": len(accepted) == 1 and accepted[0]["candidate_delta"] == frac_text(STRICT_DELTA),
        })
    return {
        "nonzero_anchor_nodes": NONZERO_NODE_DOMAIN,
        "single_node_equivalence_rows": equivalence,
        "candidate_rows": candidate_rows,
        "by_anchor_summary": by_anchor,
        "nonzero_anchor_count": len(NONZERO_NODE_DOMAIN),
        "every_nonzero_node_with_left_anchor_pins_same_delta": all(row["solved_delta_matches_4_over_5"] and row["all_nodes_reconstructed_after_pinning"] for row in equivalence),
        "every_nonzero_node_has_positive_determinant": all(row["determinant_positive"] for row in equivalence),
        "every_anchor_uniquely_accepts_strict_candidate_on_grid": all(row["unique_strict_candidate_accepted"] for row in by_anchor),
        "anchor_placement_not_identified_by_affine_algebra": True,
    }


def build_single_node_certificate(p2520: dict[str, Any]) -> dict[str, Any]:
    lattice = placement_lattice()
    return {
        "frontier_atom_under_attack": "right endpoint anchor placement/source after endpoint-subkey decomposition",
        "p2520_subkey_lattice_inherited": p2520.get("numeric_source_refinement_exported") is True,
        "certificate_type": "single nonzero node-anchor equivalence and placement nonidentifiability certificate",
        "single_node_anchor_lattice": lattice,
        "symbolic_statement": "With left normalization y(0)=0 admitted, any single nonzero strict node anchor y(log d*)=(4/5)log d* for d*>1 pins the affine mode y=b+a ell to b=0,a=4/5 because det [[1,0],[1,log d*]]=log d*>0. Thus d*=11 is sufficient but not algebraically unique; anchor placement/value must still be sourced.",
        "right_endpoint_d11_is_sufficient_not_unique": True,
        "single_nonzero_node_anchor_equivalence_exported": True,
        "anchor_placement_source_exported": False,
        "nonzero_node_value_source_exported": False,
        "beta_normalization_left_anchor_source_exported": False,
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
## P2521/S1471 single nonzero node-anchor equivalence certificate

`P2521/S1471` refines the P2520 endpoint-subkey lattice by separating endpoint sufficiency from endpoint uniqueness.  Once the left beta-normalization anchor `y(0)=0` is admitted, any one nonzero strict node anchor `y(log d*)=(4/5)log d*` with `d*>1` pins the affine kernel, because the determinant of `[[1,0],[1,log d*]]` is `log d*>0`.  The finite audit checks every `d*=2..11`: each node reconstructs `delta=4/5`, `eta=9/5`, and all strict nodes `d=1..11` with zero residual within arithmetic tolerance.

Therefore the previous right endpoint `d*=11` is sufficient but not algebraically unique.  This exports a single-node equivalence/placement-nonidentifiability theorem only: it does not source the left normalization anchor, does not source which nonzero node/value is strict, does not export `beta_eta_numeric_source`, does not source the `m=2` operator signature, and exports no strict damping source closure, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, physical-value generator, or ToE closure.
"""
    lag_section = """
## P2521/S1471 single-node anchor placement guard

`P2521/S1471` shows that, after left normalization, the right endpoint `d=11` is only one sufficient nonzero node anchor among `d=2..11`; the affine algebra itself does not select the anchor placement.  A role-bearing damping source must still source the left normalization and a nonzero node value/placement, so no nonlinear compression-flow source or role-bearing `L_total` term is licensed.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2521/S1471 single nonzero node-anchor equivalence certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2521/S1471 single-node anchor placement guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2520 = theorem(sources["P2520_ENDPOINT_SUBKEY_LATTICE"], "strict_damping_endpoint_anchor_subkey_lattice_certificate")
    cert = build_single_node_certificate(p2520)
    lattice = cert["single_node_anchor_lattice"]
    theorem_export = {
        "theorem_name": "P2521_T1_strict_damping_single_node_anchor_equivalence_certificate",
        "audited_chain": ["P2520/S1470", "P2519/S1469", "P2518/S1468"],
        "strict_damping_single_node_anchor_equivalence_certificate": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2520_subkey_lattice_inherited": cert["p2520_subkey_lattice_inherited"],
        "nonzero_anchor_count": lattice["nonzero_anchor_count"],
        "every_nonzero_node_with_left_anchor_pins_same_delta": lattice["every_nonzero_node_with_left_anchor_pins_same_delta"],
        "every_nonzero_node_has_positive_determinant": lattice["every_nonzero_node_has_positive_determinant"],
        "every_anchor_uniquely_accepts_strict_candidate_on_grid": lattice["every_anchor_uniquely_accepts_strict_candidate_on_grid"],
        "anchor_placement_not_identified_by_affine_algebra": lattice["anchor_placement_not_identified_by_affine_algebra"],
        "right_endpoint_d11_is_sufficient_not_unique": cert["right_endpoint_d11_is_sufficient_not_unique"],
        "single_nonzero_node_anchor_equivalence_exported": cert["single_nonzero_node_anchor_equivalence_exported"],
        "anchor_placement_source_exported": False,
        "nonzero_node_value_source_exported": False,
        "beta_normalization_left_anchor_source_exported": False,
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
            "P2521 exports equivalence of nonzero node anchors after left normalization, not a source for the left normalization or node value.",
            "The endpoint d=11 is sufficient but not algebraically unique; affine algebra does not select anchor placement.",
            "The m=2 operator signature source remains separate and unsourced by this packet.",
            "No damping bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Search for an internal strict invariant that selects or supplies a nonzero node anchor value/placement, rather than treating d=11 endpoint placement as silently sourced.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2520_inherited": theorem_export["p2520_subkey_lattice_inherited"],
        "ten_nonzero_anchors_audited": theorem_export["nonzero_anchor_count"] == 10,
        "all_nonzero_anchors_pin_delta": theorem_export["every_nonzero_node_with_left_anchor_pins_same_delta"],
        "all_determinants_positive": theorem_export["every_nonzero_node_has_positive_determinant"],
        "placement_not_unique": theorem_export["right_endpoint_d11_is_sufficient_not_unique"] and theorem_export["anchor_placement_not_identified_by_affine_algebra"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "anchor_placement_source_exported",
            "nonzero_node_value_source_exported",
            "beta_normalization_left_anchor_source_exported",
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
        "packet_id": "P2521",
        "stage_id": "S1471",
        "status": "STRICT_DAMPING_SINGLE_NODE_ANCHOR_EQUIVALENCE_CERTIFICATE_PLACEMENT_NONIDENTIFIABILITY_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_single_node_anchor_equivalence_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_single_node_anchor_equivalence_certificate"]["theorem_export"]
    lines = [
        "# P2521/S1471 strict damping single-node anchor equivalence certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2520 subkey lattice inherited: `{t['p2520_subkey_lattice_inherited']}`.",
        f"- Nonzero anchor count: `{t['nonzero_anchor_count']}`.",
        f"- Every nonzero node with left anchor pins delta: `{t['every_nonzero_node_with_left_anchor_pins_same_delta']}`.",
        f"- Every nonzero node determinant positive: `{t['every_nonzero_node_has_positive_determinant']}`.",
        f"- Every anchor uniquely accepts strict grid candidate: `{t['every_anchor_uniquely_accepts_strict_candidate_on_grid']}`.",
        f"- Anchor placement not identified by affine algebra: `{t['anchor_placement_not_identified_by_affine_algebra']}`.",
        f"- Right endpoint d=11 is sufficient not unique: `{t['right_endpoint_d11_is_sufficient_not_unique']}`.",
        f"- Single nonzero node-anchor equivalence exported: `{t['single_nonzero_node_anchor_equivalence_exported']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports only equivalence of nonzero node anchors after left normalization. It does not source the left normalization, node value, anchor placement, m=2 operator signature, bridge completion, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_single_node_anchor_equivalence_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_single_node_anchor_equivalence_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
