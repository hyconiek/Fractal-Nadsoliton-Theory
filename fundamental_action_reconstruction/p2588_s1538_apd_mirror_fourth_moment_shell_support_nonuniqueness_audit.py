#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from typing import Any

import numpy as np

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2569_s1519_apd_value_bridge_interpolation_dynamic_nonuniqueness_certificate import apd_rows
from p2578_s1528_apd_augmented_boundary_basis_metric_dependence_audit import BOUNDARY_TARGETS, base_and_vanish
from p2580_s1530_apd_inner_product_basis_covariance_requirement_certificate import boundary_matrix
from p2581_s1531_apd_gram_measure_moment_dependence_audit import reference_basis
from p2584_s1534_apd_full_moment_finite_support_conditional_uniqueness_audit import target_row

GEN = ROOT / "generated"
OUT = GEN / "p2588_s1538_apd_mirror_fourth_moment_shell_support_nonuniqueness_audit.json"
MD = GEN / "p2588_s1538_apd_mirror_fourth_moment_shell_support_nonuniqueness_audit.md"

SOURCE_FILES = {
    "P2416_APD_VALUE_BRIDGE": GEN / "p2416_s1366_apd_multiplicative_bridge_assembly_necessity_certificate.json",
    "P2584_FULL_MOMENT_FIXED_SUPPORT": GEN / "p2584_s1534_apd_full_moment_finite_support_conditional_uniqueness_audit.json",
    "P2587_SECOND_MOMENT_SHELL_NONUNIQUENESS": GEN / "p2587_s1537_apd_mirror_second_moment_shell_support_nonuniqueness_audit.json",
}
MIRROR_CENTER = 5.5
MIRROR_PAIR_SUM = 2.0 * MIRROR_CENTER
ENDPOINTS = [0.25, 10.75]
ENDPOINT_OFFSET = 5.25
EXPECTED_CARDINALITY = 8
INTERNAL_SECOND_SHELL = 14.0
INTERNAL_FOURTH_SHELL = 98.0
FULL_SUPPORT_MOMENT_ORDERS = list(range(EXPECTED_CARDINALITY))
SHELL_ANGLES = [0.2, 1.0, 2.0]
NEGATIVE_EXPORT_FLAGS = [
    "apd_mirror_fourth_moment_shell_source_exported",
    "apd_support_kurtosis_shell_source_exported",
    "apd_support_variance_shell_source_exported",
    "apd_support_selection_source_exported",
    "apd_finite_support_source_exported",
    "apd_positive_measure_source_exported",
    "apd_l2_inner_product_source_exported",
    "strict_dynamical_source_for_A_P_D_exported",
    "strict_phase_frequency_source_exported",
    "strict_damping_beta_eta_source_exported",
    "bridge_theorem_exported",
    "legacy_to_strict_completion_bridge_exported",
    "role_transfer_theorem_exported",
    "role_bearing_ltotal_exported",
    "qw2191_discharged_by_this_certificate",
    "toe_closure_claimed",
]


def shell_squared_offsets(angle: float) -> list[float]:
    mean = INTERNAL_SECOND_SHELL / 3.0
    radius = math.sqrt(INTERNAL_FOURTH_SHELL - INTERNAL_SECOND_SHELL ** 2 / 3.0)
    e1 = np.array([1.0, -1.0, 0.0], dtype=float) / math.sqrt(2.0)
    e2 = np.array([1.0, 1.0, -2.0], dtype=float) / math.sqrt(6.0)
    squared = mean + radius * (math.cos(angle) * e1 + math.sin(angle) * e2)
    return [float(value) for value in squared]


def support_from_angle(angle: float) -> dict[str, Any]:
    offsets = [math.sqrt(value) for value in shell_squared_offsets(angle)]
    internal_points = [MIRROR_CENTER - offset for offset in offsets] + [MIRROR_CENTER + offset for offset in offsets]
    points = sorted([ENDPOINTS[0], *internal_points, ENDPOINTS[1]])
    return {
        "name": f"eight_point_second_fourth_shell_angle_{str(angle).replace('.', 'p')}",
        "angle_radians": angle,
        "squared_offsets": shell_squared_offsets(angle),
        "offsets": offsets,
        "points": points,
        "weights": [1.0] * EXPECTED_CARDINALITY,
    }


MIRROR_FOURTH_MOMENT_SUPPORTS = [support_from_angle(angle) for angle in SHELL_ANGLES]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
    ], cwd=REPO, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2588|S1538|APD mirror fourth moment|mirror fourth moment.*APD",
        "intended_research_nonduplication": "APD.*fourth moment shell|fourth moment shell.*APD|APD.*support kurtosis shell|support kurtosis shell.*APD|APD.*central fourth support moment|central fourth support moment.*APD|APD.*mirror kurtosis support|mirror kurtosis support.*APD",
        "apd_precursors": "P2416|S1366|P2584|S1534|P2587|S1537|APD.*second moment shell|APD.*full moment|strict_dynamical_source_for_A_P_D",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer theorem|QW-2191|ToE closure|source theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def vandermonde_moment_matrix(points: np.ndarray) -> np.ndarray:
    return np.vstack([points ** order for order in FULL_SUPPORT_MOMENT_ORDERS])


def full_support_moments(points: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return vandermonde_moment_matrix(points) @ weights


def shell_support_witness(support: dict[str, Any]) -> dict[str, Any]:
    points = np.array(sorted(support["points"]), dtype=float)
    weights = np.array(support["weights"], dtype=float)
    matrix = vandermonde_moment_matrix(points)
    moments = full_support_moments(points, weights)
    recovered = np.linalg.solve(matrix, moments)
    pair_sums = [float(points[index] + points[-1 - index]) for index in range(len(points) // 2)]
    internal_squared = np.array(support["squared_offsets"], dtype=float)
    internal_second = float(np.sum(internal_squared))
    internal_fourth = float(np.sum(internal_squared ** 2))
    central_second = float(np.sum((points - MIRROR_CENTER) ** 2))
    central_fourth = float(np.sum((points - MIRROR_CENTER) ** 4))
    return {
        "support_name": support["name"],
        "angle_radians": float(support["angle_radians"]),
        "points": [float(value) for value in points],
        "input_weights": [float(value) for value in weights],
        "squared_offsets": [float(value) for value in support["squared_offsets"]],
        "offsets": [float(value) for value in support["offsets"]],
        "support_cardinality": int(len(points)),
        "endpoints": [float(np.min(points)), float(np.max(points))],
        "mirror_center": MIRROR_CENTER,
        "pair_sums": pair_sums,
        "is_pair_symmetric_about_center": bool(all(abs(value - MIRROR_PAIR_SUM) <= 1.0e-12 for value in pair_sums)),
        "internal_second_shell": internal_second,
        "internal_fourth_shell": internal_fourth,
        "has_expected_internal_second_fourth_shell": bool(abs(internal_second - INTERNAL_SECOND_SHELL) <= 1.0e-10 and abs(internal_fourth - INTERNAL_FOURTH_SHELL) <= 1.0e-10),
        "central_second_moment_unweighted": central_second,
        "central_fourth_moment_unweighted": central_fourth,
        "full_moment_orders": FULL_SUPPORT_MOMENT_ORDERS,
        "full_moments": [float(value) for value in moments],
        "vandermonde_determinant": float(np.linalg.det(matrix)),
        "vandermonde_rank": int(np.linalg.matrix_rank(matrix)),
        "recovered_weights_from_full_moments": [float(value) for value in recovered],
        "max_abs_recovered_weight_error": float(np.max(np.abs(recovered - weights))),
        "full_moments_conditionally_select_weights_on_fixed_support": bool(np.max(np.abs(recovered - weights)) <= 1.0e-5),
        "all_recovered_weights_positive": bool(np.min(recovered) > 0.0),
    }


def shell_support_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base, vanish = base_and_vanish(rows)
    basis = reference_basis(vanish)
    matrix = boundary_matrix(basis)
    witnesses = [shell_support_witness(support) for support in MIRROR_FOURTH_MOMENT_SUPPORTS]
    target_rows = [target_row(target, base, basis, matrix, rows, witnesses) for target in BOUNDARY_TARGETS]
    central_seconds = {round(witness["central_second_moment_unweighted"], 10) for witness in witnesses}
    central_fourths = {round(witness["central_fourth_moment_unweighted"], 10) for witness in witnesses}
    return {
        "numpy_version": np.__version__,
        "mirror_center": MIRROR_CENTER,
        "mirror_pair_sum": MIRROR_PAIR_SUM,
        "expected_cardinality": EXPECTED_CARDINALITY,
        "endpoints": ENDPOINTS,
        "internal_second_shell": INTERNAL_SECOND_SHELL,
        "internal_fourth_shell": INTERNAL_FOURTH_SHELL,
        "full_support_moment_orders": FULL_SUPPORT_MOMENT_ORDERS,
        "support_witnesses": witnesses,
        "target_rows": target_rows,
        "support_count": len(witnesses),
        "target_count": len(target_rows),
        "all_supports_have_expected_cardinality": all(witness["support_cardinality"] == EXPECTED_CARDINALITY for witness in witnesses),
        "all_supports_share_endpoints": all(witness["endpoints"] == ENDPOINTS for witness in witnesses),
        "all_supports_pair_symmetric_about_center": all(witness["is_pair_symmetric_about_center"] for witness in witnesses),
        "all_supports_share_second_fourth_shell": all(witness["has_expected_internal_second_fourth_shell"] for witness in witnesses),
        "all_supports_share_unweighted_central_second_moment": len(central_seconds) == 1,
        "all_supports_share_unweighted_central_fourth_moment": len(central_fourths) == 1,
        "all_fixed_support_weights_conditionally_unique": all(witness["full_moments_conditionally_select_weights_on_fixed_support"] for witness in witnesses),
        "all_fixed_supports_have_vandermonde_rank_eight": all(witness["vandermonde_rank"] == EXPECTED_CARDINALITY for witness in witnesses),
        "all_targets_shell_support_nonunique": all(row["fixed_support_full_moments_do_not_select_support"] for row in target_rows),
        "all_gram_metrics_positive_definite": all(row["all_gram_metrics_positive_definite"] for row in target_rows),
        "all_minimizers_preserve_nodes_and_boundaries": all(row["all_minimizers_preserve_apd_nodes"] and row["all_minimizers_satisfy_boundary_targets"] for row in target_rows),
        "second_fourth_moment_shell_is_unsourced_support_selector": True,
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2416_payload = load_json(SOURCE_FILES["P2416_APD_VALUE_BRIDGE"])
    p2584_payload = load_json(SOURCE_FILES["P2584_FULL_MOMENT_FIXED_SUPPORT"])
    p2587_payload = load_json(SOURCE_FILES["P2587_SECOND_MOMENT_SHELL_NONUNIQUENESS"])
    p2584 = theorem(p2584_payload, "apd_full_moment_finite_support_conditional_uniqueness_audit")
    p2587 = theorem(p2587_payload, "apd_mirror_second_moment_shell_support_nonuniqueness_audit")
    rows = apd_rows(p2416_payload)
    audit = shell_support_audit(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2588_T1_apd_mirror_fourth_moment_shell_support_nonuniqueness_audit",
        "audited_chain": ["P2416/S1366", "P2584/S1534", "P2587/S1537"],
        "frontier_atom_under_attack": "strict_dynamical_source_for_A_P_D",
        "candidate_principle_under_test": "select APD finite support by endpoints, mirror symmetry, cardinality eight, fixed second/fourth support shells, and fixed-support full moments",
        "p2416_apd_value_bridge_inherited": p2416_payload.get("apd_multiplicative_bridge_assembly_necessity_certificate", {}).get("finite_witness_certificate", {}).get("apd_value_assembly_witness_ready") is True,
        "p2584_fixed_support_conditional_uniqueness_inherited": p2584.get("full_moments_select_weights_only_after_support_is_fixed") is True,
        "p2587_second_moment_shell_nonuniqueness_inherited": p2587.get("mirror_second_moment_shell_does_not_select_apd_support") is True,
        "apd_node_rows": rows,
        "apd_mirror_fourth_moment_shell_support_nonuniqueness_audit": audit,
        "mirror_second_fourth_moment_shell_does_not_select_apd_support": audit["all_targets_shell_support_nonunique"],
        "fixed_support_full_moments_remain_conditional": True,
        "apd_dynamic_source_remains_unsourced": True,
        "recommended_next_honest_step": (
            "Do not promote endpoints, mirror symmetry, cardinality, or fixed second/fourth support shells into an APD support source. P2588 keeps all of those constraints and still obtains distinct supports with conditionally unique weights and support-dependent APD dynamics. The next honest step is a strict nadsoliton-derived support/density law or a genuine internal APD dynamics theorem."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2416_apd_value_bridge_inherited": theorem_export["p2416_apd_value_bridge_inherited"],
        "p2584_fixed_support_conditional_uniqueness_inherited": theorem_export["p2584_fixed_support_conditional_uniqueness_inherited"],
        "p2587_second_moment_shell_nonuniqueness_inherited": theorem_export["p2587_second_moment_shell_nonuniqueness_inherited"],
        "three_supports_audited": audit["support_count"] == 3,
        "three_targets_audited": audit["target_count"] == 3,
        "all_supports_have_expected_cardinality": audit["all_supports_have_expected_cardinality"],
        "all_supports_share_endpoints": audit["all_supports_share_endpoints"],
        "all_supports_pair_symmetric_about_center": audit["all_supports_pair_symmetric_about_center"],
        "all_supports_share_second_fourth_shell": audit["all_supports_share_second_fourth_shell"],
        "all_fixed_support_weights_conditionally_unique": audit["all_fixed_support_weights_conditionally_unique"],
        "all_targets_shell_support_nonunique": audit["all_targets_shell_support_nonunique"],
        "all_gram_metrics_positive_definite": audit["all_gram_metrics_positive_definite"],
        "all_minimizers_preserve_nodes_and_boundaries": audit["all_minimizers_preserve_nodes_and_boundaries"],
        "apd_dynamic_source_not_exported": theorem_export["strict_dynamical_source_for_A_P_D_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2588",
        "stage_id": "S1538",
        "status": "P2588_APD_MIRROR_FOURTH_MOMENT_SHELL_SUPPORT_NONUNIQUENESS_AUDIT_NO_APD_DYNAMIC_SOURCE_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "apd_mirror_fourth_moment_shell_support_nonuniqueness_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2416_APD_VALUE_BRIDGE": sha256_json(p2416_payload),
                "P2584_FULL_MOMENT_FIXED_SUPPORT": sha256_json(p2584_payload),
                "P2587_SECOND_MOMENT_SHELL_NONUNIQUENESS": sha256_json(p2587_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["apd_mirror_fourth_moment_shell_support_nonuniqueness_audit"]["theorem_export"]
    audit = t["apd_mirror_fourth_moment_shell_support_nonuniqueness_audit"]
    lines = [
        "# P2588/S1538 APD mirror fourth-moment shell support nonuniqueness audit", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- Supports audited: `{audit['support_count']}`.",
        f"- Cardinality: `{audit['expected_cardinality']}`.",
        f"- Mirror center: `{audit['mirror_center']}`.",
        f"- Internal second/fourth shells: `{audit['internal_second_shell']}` / `{audit['internal_fourth_shell']}`.",
        f"- Shell support nonunique: `{audit['all_targets_shell_support_nonunique']}`.",
        f"- Strict APD dynamic source exported: `{t['strict_dynamical_source_for_A_P_D_exported']}`.", "",
        "## Interpretation", "",
        "P2588 strengthens P2587 by adding a fixed fourth-moment shell on eight-point mirror supports.  The supports share endpoints, cardinality, reflection symmetry, and unweighted central second/fourth moments; fixed-support full moments recover positive weights, but selected off-node APD dynamics still changes with the support.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No fourth-moment shell support source, kurtosis-shell source, variance-shell source, support-selection source, finite-support source, positive-measure source, strict A/P/D dynamic source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['apd_mirror_fourth_moment_shell_support_nonuniqueness_audit']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2588/S1538 APD mirror fourth-moment shell support nonuniqueness guard

`P2588/S1538` strengthens P2587 by imposing eight-point mirror supports with shared endpoints, cardinality `8`, and common internal second/fourth moment shells.  Fixed-support full moments still recover weights only after the support is chosen, and the off-node APD dynamics remains support-dependent.  Therefore a mirror kurtosis shell is not a strict APD support source.
""".strip()
    lag_section = """
## P2588/S1538 APD mirror fourth-moment shell support nonuniqueness Ltotal guard

`P2588/S1538` blocks a role-bearing APD Gram term in `L_total` from being justified by finite support variance/kurtosis shells.  Moment-shell geometry is still selector data; strict nadsoliton dynamics must source the actual support/density law.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2588/S1538 APD mirror fourth-moment shell support nonuniqueness guard", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2588/S1538 APD mirror fourth-moment shell support nonuniqueness Ltotal guard", lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
