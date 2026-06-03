#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from fractions import Fraction
from itertools import combinations
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
OUT = GEN / "p2522_s1472_strict_damping_two_node_anchor_basis_equivalence_certificate.json"
MD = GEN / "p2522_s1472_strict_damping_two_node_anchor_basis_equivalence_certificate.md"

SOURCE_FILES = {
    "P2521_SINGLE_NODE_ANCHOR_EQUIVALENCE": GEN / "p2521_s1471_strict_damping_single_node_anchor_equivalence_certificate.json",
}

STRICT_DELTA = Fraction(4, 5)
STRICT_ETA = Fraction(9, 5)
NODE_DOMAIN = list(range(1, 12))
REPRESENTATIVE_PAIRS = [(1, 11), (2, 3), (2, 11), (10, 11)]
INTERCEPT_CANDIDATES = [Fraction(-1, 1), Fraction(-1, 2), Fraction(0, 1), Fraction(1, 2), Fraction(1, 1)]
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
        "new_packet": "P2522|S1472|two-node anchor basis equivalence|two nonzero node anchors|node-pair anchor basis|anchor basis nonuniqueness|two-node source basis",
        "precursor_packets": "P2521|S1471|single-node anchor equivalence|P2520|endpoint anchor subkey lattice|P2519|endpoint-anchor acceptance",
        "basis_language": "two-node anchor|node-pair anchor|anchor basis|basis equivalence|distinct node anchors|log\\(d_j/d_i\\)",
        "source_blockers": "beta_eta_numeric_source|anchor placement source|node value source|strict source theorem|proper subset obstruction",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def frac_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def strict_y(d: int) -> float:
    return float(STRICT_DELTA) * math.log(d)


def solve_pair(d_i: int, d_j: int) -> dict[str, float]:
    ell_i = math.log(d_i)
    ell_j = math.log(d_j)
    y_i = strict_y(d_i)
    y_j = strict_y(d_j)
    determinant = ell_j - ell_i
    slope = (y_j - y_i) / determinant
    intercept = y_i - slope * ell_i
    return {
        "ell_i": ell_i,
        "ell_j": ell_j,
        "y_i": y_i,
        "y_j": y_j,
        "determinant": determinant,
        "intercept": intercept,
        "slope": slope,
    }


def two_node_basis_rows() -> list[dict[str, Any]]:
    rows = []
    for d_i, d_j in combinations(NODE_DOMAIN, 2):
        solved = solve_pair(d_i, d_j)
        residuals = []
        for d in NODE_DOMAIN:
            ell = math.log(d)
            residuals.append(strict_y(d) - (solved["intercept"] + solved["slope"] * ell))
        rows.append({
            "node_pair": [d_i, d_j],
            "constraint_matrix": [[1.0, solved["ell_i"]], [1.0, solved["ell_j"]]],
            "constraint_rhs": [solved["y_i"], solved["y_j"]],
            "determinant_log_ratio": solved["determinant"],
            "determinant_positive": solved["determinant"] > 0,
            "solved_intercept_log_beta": solved["intercept"],
            "solved_intercept_matches_zero": abs(solved["intercept"]) < 1e-14,
            "solved_delta": solved["slope"],
            "solved_delta_matches_4_over_5": abs(solved["slope"] - float(STRICT_DELTA)) < 1e-14,
            "solved_eta": 1.0 + solved["slope"],
            "solved_eta_matches_9_over_5": abs((1.0 + solved["slope"]) - float(STRICT_ETA)) < 1e-14,
            "max_abs_all_node_residual_after_pair_pinning": max(abs(value) for value in residuals),
            "all_nodes_reconstructed_after_pair_pinning": max(abs(value) for value in residuals) < 1e-14,
            "derives_left_normalization_y0": abs(solved["intercept"]) < 1e-14,
        })
    return rows


def representative_candidate_grid() -> dict[str, Any]:
    rows = []
    for d_i, d_j in REPRESENTATIVE_PAIRS:
        ell_i = math.log(d_i)
        ell_j = math.log(d_j)
        target_i = strict_y(d_i)
        target_j = strict_y(d_j)
        for intercept in INTERCEPT_CANDIDATES:
            for slope in SLOPE_CANDIDATES:
                residual_i = float(intercept) + float(slope) * ell_i - target_i
                residual_j = float(intercept) + float(slope) * ell_j - target_j
                rows.append({
                    "node_pair": [d_i, d_j],
                    "intercept_log_beta_candidate": frac_text(intercept),
                    "slope_delta_candidate": frac_text(slope),
                    "eta_candidate_if_slope_delta": frac_text(Fraction(1, 1) + slope),
                    "residual_i": residual_i,
                    "residual_j": residual_j,
                    "accepted_by_pair": abs(residual_i) < 1e-14 and abs(residual_j) < 1e-14,
                })
    by_pair = []
    for pair in REPRESENTATIVE_PAIRS:
        accepted = [row for row in rows if tuple(row["node_pair"]) == pair and row["accepted_by_pair"]]
        by_pair.append({
            "node_pair": list(pair),
            "accepted_count": len(accepted),
            "accepted_pair": accepted[0] if len(accepted) == 1 else {},
            "unique_strict_pair_accepted": len(accepted) == 1 and accepted[0]["intercept_log_beta_candidate"] == "0" and accepted[0]["slope_delta_candidate"] == "4/5",
        })
    return {
        "representative_pairs": [list(pair) for pair in REPRESENTATIVE_PAIRS],
        "candidate_intercepts": [frac_text(value) for value in INTERCEPT_CANDIDATES],
        "candidate_slopes": [frac_text(value) for value in SLOPE_CANDIDATES],
        "rows": rows,
        "row_count": len(rows),
        "by_pair_summary": by_pair,
        "every_representative_pair_uniquely_accepts_strict_pair": all(row["unique_strict_pair_accepted"] for row in by_pair),
    }


def build_two_node_basis_certificate(p2521: dict[str, Any]) -> dict[str, Any]:
    basis_rows = two_node_basis_rows()
    candidate_grid = representative_candidate_grid()
    return {
        "frontier_atom_under_attack": "node-anchor basis/source nonuniqueness for beta_eta_numeric_target",
        "p2521_single_node_equivalence_inherited": p2521.get("single_nonzero_node_anchor_equivalence_exported") is True,
        "certificate_type": "two-node anchor basis equivalence and basis-nonuniqueness certificate",
        "symbolic_statement": "For affine y=b+a ell, any two distinct strict node anchors y(log d_i)=(4/5)log d_i and y(log d_j)=(4/5)log d_j have determinant log(d_j/d_i) != 0 and solve b=0,a=4/5. Therefore the left-normalization-plus-one-node basis is sufficient but not unique; any distinct node pair is an equivalent conditional basis if its node values are strictly sourced.",
        "two_node_basis_rows": basis_rows,
        "representative_candidate_grid": candidate_grid,
        "node_pair_count": len(basis_rows),
        "every_distinct_node_pair_has_positive_determinant_in_ordered_domain": all(row["determinant_positive"] for row in basis_rows),
        "every_distinct_node_pair_derives_beta_normalization": all(row["derives_left_normalization_y0"] for row in basis_rows),
        "every_distinct_node_pair_pins_delta_eta": all(row["solved_delta_matches_4_over_5"] and row["solved_eta_matches_9_over_5"] for row in basis_rows),
        "every_distinct_node_pair_reconstructs_all_nodes": all(row["all_nodes_reconstructed_after_pair_pinning"] for row in basis_rows),
        "representative_grid_unique_strict_pair": candidate_grid["every_representative_pair_uniquely_accepts_strict_pair"],
        "anchor_basis_equivalence_exported": True,
        "anchor_basis_source_exported": False,
        "node_pair_value_source_exported": False,
        "anchor_placement_source_exported": False,
        "beta_normalization_left_anchor_source_exported": False,
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
## P2522/S1472 two-node anchor basis equivalence certificate

`P2522/S1472` continues P2521 by showing that the left-normalization-plus-one-node basis is itself not unique.  For the affine running exponent `y=b+a ell`, any two distinct strict node anchors `y(log d_i)=(4/5)log d_i` and `y(log d_j)=(4/5)log d_j` have determinant `log(d_j/d_i) != 0` and solve `b=0`, `a=4/5`, hence `beta=1`, `eta=9/5`.  The finite audit checks all 55 distinct node pairs in `d=1..11`: every pair derives the left normalization, pins the same slope, and reconstructs all audited strict nodes with zero residual within arithmetic tolerance.

This exports an anchor-basis equivalence theorem only.  It does not source which node pair/value basis strict dynamics supplies, does not export `beta_eta_numeric_source`, does not source the `m=2` operator signature, and exports no strict damping source closure, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, physical-value generator, or ToE closure.
"""
    lag_section = """
## P2522/S1472 two-node anchor basis guard

`P2522/S1472` shows that two distinct strict node values can conditionally replace the left-normalization-plus-one-node anchor basis and derive `beta=1, eta=9/5`.  Because the node-pair basis and its values remain unsourced, this does not license a nonlinear compression-flow source or role-bearing `L_total` term.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2522/S1472 two-node anchor basis equivalence certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2522/S1472 two-node anchor basis guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2521 = theorem(sources["P2521_SINGLE_NODE_ANCHOR_EQUIVALENCE"], "strict_damping_single_node_anchor_equivalence_certificate")
    cert = build_two_node_basis_certificate(p2521)
    theorem_export = {
        "theorem_name": "P2522_T1_strict_damping_two_node_anchor_basis_equivalence_certificate",
        "audited_chain": ["P2521/S1471", "P2520/S1470", "P2519/S1469"],
        "strict_damping_two_node_anchor_basis_equivalence_certificate": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2521_single_node_equivalence_inherited": cert["p2521_single_node_equivalence_inherited"],
        "node_pair_count": cert["node_pair_count"],
        "every_distinct_node_pair_has_positive_determinant_in_ordered_domain": cert["every_distinct_node_pair_has_positive_determinant_in_ordered_domain"],
        "every_distinct_node_pair_derives_beta_normalization": cert["every_distinct_node_pair_derives_beta_normalization"],
        "every_distinct_node_pair_pins_delta_eta": cert["every_distinct_node_pair_pins_delta_eta"],
        "every_distinct_node_pair_reconstructs_all_nodes": cert["every_distinct_node_pair_reconstructs_all_nodes"],
        "representative_grid_unique_strict_pair": cert["representative_grid_unique_strict_pair"],
        "anchor_basis_equivalence_exported": cert["anchor_basis_equivalence_exported"],
        "anchor_basis_source_exported": False,
        "node_pair_value_source_exported": False,
        "anchor_placement_source_exported": False,
        "beta_normalization_left_anchor_source_exported": False,
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
            "P2522 exports equivalence of two-node anchor bases, not a source for any node-pair basis or node values.",
            "The left-normalization-plus-one-node basis is sufficient but not unique; affine algebra does not select the source basis.",
            "The m=2 operator signature source remains separate and unsourced by this packet.",
            "No damping bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Search for strict dynamics that exports a node-pair/value basis or a basis-independent invariant producing the same affine numeric target without axiomatically choosing anchors.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2521_inherited": theorem_export["p2521_single_node_equivalence_inherited"],
        "all_55_pairs_audited": theorem_export["node_pair_count"] == 55,
        "all_pairs_pin_target": theorem_export["every_distinct_node_pair_derives_beta_normalization"] and theorem_export["every_distinct_node_pair_pins_delta_eta"],
        "all_pairs_reconstruct_nodes": theorem_export["every_distinct_node_pair_reconstructs_all_nodes"],
        "representative_grid_unique": theorem_export["representative_grid_unique_strict_pair"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "anchor_basis_source_exported",
            "node_pair_value_source_exported",
            "anchor_placement_source_exported",
            "beta_normalization_left_anchor_source_exported",
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
        "packet_id": "P2522",
        "stage_id": "S1472",
        "status": "STRICT_DAMPING_TWO_NODE_ANCHOR_BASIS_EQUIVALENCE_CERTIFICATE_BASIS_NONUNIQUENESS_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_two_node_anchor_basis_equivalence_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_two_node_anchor_basis_equivalence_certificate"]["theorem_export"]
    lines = [
        "# P2522/S1472 strict damping two-node anchor basis equivalence certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2521 single-node equivalence inherited: `{t['p2521_single_node_equivalence_inherited']}`.",
        f"- Node-pair count: `{t['node_pair_count']}`.",
        f"- Every ordered-domain pair determinant positive: `{t['every_distinct_node_pair_has_positive_determinant_in_ordered_domain']}`.",
        f"- Every distinct node pair derives beta normalization: `{t['every_distinct_node_pair_derives_beta_normalization']}`.",
        f"- Every distinct node pair pins delta/eta: `{t['every_distinct_node_pair_pins_delta_eta']}`.",
        f"- Every distinct node pair reconstructs all nodes: `{t['every_distinct_node_pair_reconstructs_all_nodes']}`.",
        f"- Representative grid unique strict pair: `{t['representative_grid_unique_strict_pair']}`.",
        f"- Anchor-basis equivalence exported: `{t['anchor_basis_equivalence_exported']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports only equivalence of two-node anchor bases. It does not source the node-pair basis, node values, m=2 operator signature, bridge completion, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_two_node_anchor_basis_equivalence_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_two_node_anchor_basis_equivalence_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
