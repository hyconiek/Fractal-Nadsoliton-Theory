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
OUT = GEN / "p2520_s1470_strict_damping_endpoint_anchor_subkey_lattice_certificate.json"
MD = GEN / "p2520_s1470_strict_damping_endpoint_anchor_subkey_lattice_certificate.md"

SOURCE_FILES = {
    "P2519_ENDPOINT_ANCHOR_ACCEPTANCE": GEN / "p2519_s1469_strict_damping_biharmonic_endpoint_anchor_acceptance_certificate.json",
}

STRICT_DELTA = Fraction(4, 5)
STRICT_BETA_LOG = Fraction(0, 1)
STRICT_ETA = Fraction(9, 5)
INTERCEPT_CANDIDATES = [Fraction(-1, 1), Fraction(-1, 2), Fraction(0, 1), Fraction(1, 2), Fraction(1, 1)]
SLOPE_CANDIDATES = [Fraction(-1, 1), Fraction(0, 1), Fraction(1, 2), Fraction(4, 5), Fraction(1, 1), Fraction(9, 5), Fraction(2, 1)]
ANCHORS = ["beta_normalization_left_anchor", "right_endpoint_value_anchor"]


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
        "new_packet": "P2520|S1470|endpoint anchor subkey lattice|anchor subkey lattice|beta normalization anchor|right endpoint value anchor|numeric-source subkey lattice",
        "precursor_packets": "P2519|S1469|biharmonic endpoint-anchor acceptance|P2518|affine-slope nonidentifiability|P2516|dual-key source acceptance",
        "subkey_language": "beta normalization anchor|eta slope anchor|right endpoint value|anchor rank lattice|affine anchor rank|two endpoint subkeys",
        "source_blockers": "beta_eta_numeric_source|endpoint_anchor_source|strict source theorem|proper subset obstruction|non-strict",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def frac_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def rank_2col(rows: list[list[float]]) -> int:
    if not rows:
        return 0
    nonzero = [row for row in rows if abs(row[0]) > 1e-14 or abs(row[1]) > 1e-14]
    if not nonzero:
        return 0
    first = nonzero[0]
    for row in nonzero[1:]:
        det = first[0] * row[1] - first[1] * row[0]
        if abs(det) > 1e-14:
            return 2
    return 1


def subkey_rank_lattice() -> dict[str, Any]:
    log11 = math.log(11)
    rows = []
    for mask in range(1 << len(ANCHORS)):
        active = [ANCHORS[index] for index in range(len(ANCHORS)) if mask & (1 << index)]
        matrix = []
        rhs = []
        if "beta_normalization_left_anchor" in active:
            matrix.append([1.0, 0.0])
            rhs.append(0.0)
        if "right_endpoint_value_anchor" in active:
            matrix.append([1.0, log11])
            rhs.append(float(STRICT_DELTA) * log11)
        rank = rank_2col(matrix)
        nullity = 2 - rank
        beta_fixed = "beta_normalization_left_anchor" in active and rank >= 1
        delta_fixed = set(active) == set(ANCHORS) and rank == 2
        if not active:
            solution_form = "b free, a free"
        elif active == ["beta_normalization_left_anchor"]:
            solution_form = "b=0, a free"
        elif active == ["right_endpoint_value_anchor"]:
            solution_form = "b=(4/5-a)*log(11), a free"
        else:
            solution_form = "b=0, a=4/5"
        rows.append({
            "mask": format(mask, "02b"),
            "active_anchors": active,
            "constraint_matrix": matrix,
            "constraint_rhs": rhs,
            "rank": rank,
            "nullity": nullity,
            "beta_log_fixed_to_zero": beta_fixed,
            "delta_fixed_to_4_over_5": delta_fixed,
            "eta_fixed_to_9_over_5": delta_fixed,
            "unique_numeric_beta_eta_target": beta_fixed and delta_fixed,
            "solution_form": solution_form,
            "classification": "accepted_full_numeric_anchor_pair" if beta_fixed and delta_fixed else "proper_subkey_obstruction",
        })
    return {
        "anchors": ANCHORS,
        "rows": rows,
        "row_count": len(rows),
        "accepted_row_count": sum(1 for row in rows if row["unique_numeric_beta_eta_target"]),
        "proper_subkey_obstruction_count": sum(1 for row in rows if not row["unique_numeric_beta_eta_target"]),
        "only_both_subkeys_accept": all((set(row["active_anchors"]) == set(ANCHORS)) == row["unique_numeric_beta_eta_target"] for row in rows),
        "left_only_fixes_beta_not_eta": any(row["active_anchors"] == ["beta_normalization_left_anchor"] and row["beta_log_fixed_to_zero"] and not row["eta_fixed_to_9_over_5"] for row in rows),
        "right_only_fixes_neither_beta_nor_eta": any(row["active_anchors"] == ["right_endpoint_value_anchor"] and not row["beta_log_fixed_to_zero"] and not row["eta_fixed_to_9_over_5"] for row in rows),
    }


def candidate_pair_audit() -> dict[str, Any]:
    log11 = math.log(11)
    rows = []
    for intercept in INTERCEPT_CANDIDATES:
        for slope in SLOPE_CANDIDATES:
            left_residual = float(intercept)
            right_residual = float(intercept) + float(slope) * log11 - float(STRICT_DELTA) * log11
            rows.append({
                "intercept_log_beta_candidate": frac_text(intercept),
                "slope_delta_candidate": frac_text(slope),
                "eta_candidate_if_slope_delta": frac_text(Fraction(1, 1) + slope),
                "left_anchor_residual": left_residual,
                "right_anchor_residual": right_residual,
                "accepted_by_both_anchors": abs(left_residual) < 1e-14 and abs(right_residual) < 1e-14,
                "accepted_by_left_only": abs(left_residual) < 1e-14,
                "accepted_by_right_only": abs(right_residual) < 1e-14,
            })
    accepted_both = [row for row in rows if row["accepted_by_both_anchors"]]
    return {
        "candidate_intercepts": [frac_text(value) for value in INTERCEPT_CANDIDATES],
        "candidate_slopes": [frac_text(value) for value in SLOPE_CANDIDATES],
        "rows": rows,
        "grid_row_count": len(rows),
        "accepted_by_both_count": len(accepted_both),
        "accepted_by_both_unique_pair": accepted_both[0] if len(accepted_both) == 1 else {},
        "accepted_by_left_only_count": sum(1 for row in rows if row["accepted_by_left_only"]),
        "accepted_by_right_only_count": sum(1 for row in rows if row["accepted_by_right_only"]),
        "only_strict_pair_accepted_by_both": len(accepted_both) == 1 and accepted_both[0]["intercept_log_beta_candidate"] == "0" and accepted_both[0]["slope_delta_candidate"] == "4/5",
    }


def source_obligation_normal_form() -> dict[str, Any]:
    return {
        "beta_eta_numeric_source_refinement": "beta_eta_numeric_source = beta_normalization_left_anchor_source AND right_endpoint_value_anchor_source",
        "strict_damping_source_refinement": "strict_damping_beta_eta_source = m2_operator_signature_source AND beta_normalization_left_anchor_source AND right_endpoint_value_anchor_source",
        "proper_subkeys_rejected": True,
        "subkeys_are_targets_not_sources_here": True,
    }


def build_subkey_lattice_certificate(p2519: dict[str, Any]) -> dict[str, Any]:
    lattice = subkey_rank_lattice()
    candidates = candidate_pair_audit()
    normal_form = source_obligation_normal_form()
    return {
        "frontier_atom_under_attack": "beta_eta_numeric_source decomposed into endpoint-anchor source subkeys",
        "p2519_endpoint_anchor_acceptance_inherited": p2519.get("conditional_beta_eta_target_recovered_from_anchors") is True,
        "certificate_type": "endpoint-anchor subkey lattice and rank certificate",
        "subkey_rank_lattice": lattice,
        "candidate_pair_audit": candidates,
        "source_obligation_normal_form": normal_form,
        "only_both_endpoint_subkeys_accept_numeric_target": lattice["only_both_subkeys_accept"] and candidates["only_strict_pair_accepted_by_both"],
        "left_subkey_alone_is_proper_obstruction": lattice["left_only_fixes_beta_not_eta"],
        "right_subkey_alone_is_proper_obstruction": lattice["right_only_fixes_neither_beta_nor_eta"],
        "numeric_source_refinement_exported": True,
        "beta_normalization_left_anchor_source_exported": False,
        "right_endpoint_value_anchor_source_exported": False,
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
## P2520/S1470 endpoint-anchor subkey lattice certificate

`P2520/S1470` refines the P2519 conditional endpoint-anchor acceptance result into a rank/subkey lattice.  The numeric anchor target decomposes into two independent endpoint source obligations: the left beta-normalization anchor `y(0)=0`, and the right endpoint value anchor `y(log 11)=(4/5)log 11`.  The affine constraint matrix has rank/nullity profile: no anchors leaves `(b,a)` two-dimensional, the left anchor fixes `b=0` but leaves slope free, the right anchor leaves a one-parameter intercept/slope tradeoff, and only both anchors give rank 2 and the unique target `b=0, a=4/5`, hence `beta=1, eta=9/5`.

The finite candidate-pair audit confirms the same proper-subkey obstruction on an explicit intercept/slope grid: only `(log beta, delta)=(0,4/5)` passes both anchors.  This exports a refined source-obligation normal form only; it does not source either endpoint subkey, does not export `beta_eta_numeric_source`, does not source the `m=2` operator signature, and exports no strict damping source closure, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, physical-value generator, or ToE closure.
"""
    lag_section = """
## P2520/S1470 endpoint-anchor subkey guard

`P2520/S1470` records that the P2519 endpoint anchors are themselves two independent source subkeys.  The left anchor supplies beta normalization but not the slope, while the right anchor alone permits an intercept/slope tradeoff; only both subkeys pin `beta=1, eta=9/5`.  Because neither subkey is sourced here, no nonlinear compression-flow source or role-bearing `L_total` term is licensed.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2520/S1470 endpoint-anchor subkey lattice certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2520/S1470 endpoint-anchor subkey guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2519 = theorem(sources["P2519_ENDPOINT_ANCHOR_ACCEPTANCE"], "strict_damping_biharmonic_endpoint_anchor_acceptance_certificate")
    cert = build_subkey_lattice_certificate(p2519)
    lattice = cert["subkey_rank_lattice"]
    candidates = cert["candidate_pair_audit"]
    theorem_export = {
        "theorem_name": "P2520_T1_strict_damping_endpoint_anchor_subkey_lattice_certificate",
        "audited_chain": ["P2519/S1469", "P2518/S1468", "P2516/S1466"],
        "strict_damping_endpoint_anchor_subkey_lattice_certificate": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2519_endpoint_anchor_acceptance_inherited": cert["p2519_endpoint_anchor_acceptance_inherited"],
        "subkey_lattice_row_count": lattice["row_count"],
        "accepted_row_count": lattice["accepted_row_count"],
        "proper_subkey_obstruction_count": lattice["proper_subkey_obstruction_count"],
        "only_both_subkeys_accept": lattice["only_both_subkeys_accept"],
        "left_subkey_alone_is_proper_obstruction": cert["left_subkey_alone_is_proper_obstruction"],
        "right_subkey_alone_is_proper_obstruction": cert["right_subkey_alone_is_proper_obstruction"],
        "candidate_grid_row_count": candidates["grid_row_count"],
        "candidate_grid_only_strict_pair_accepted_by_both": candidates["only_strict_pair_accepted_by_both"],
        "numeric_source_refinement_exported": cert["numeric_source_refinement_exported"],
        "beta_normalization_left_anchor_source_exported": False,
        "right_endpoint_value_anchor_source_exported": False,
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
            "P2520 exports a source-obligation refinement, not either endpoint-anchor source theorem.",
            "The left beta-normalization subkey alone and the right endpoint-value subkey alone are both proper-subkey obstructions.",
            "The m=2 operator signature source remains separate and unsourced by this packet.",
            "No damping bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Find a strict nadsoliton mechanism that sources both endpoint subkeys, or prove that one endpoint subkey follows from an already strict internal invariant without axiomatic insertion.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2519_inherited": theorem_export["p2519_endpoint_anchor_acceptance_inherited"],
        "lattice_shape_ok": theorem_export["subkey_lattice_row_count"] == 4,
        "unique_accepting_row": theorem_export["accepted_row_count"] == 1 and theorem_export["only_both_subkeys_accept"],
        "proper_subkeys_blocked": theorem_export["proper_subkey_obstruction_count"] == 3 and theorem_export["left_subkey_alone_is_proper_obstruction"] and theorem_export["right_subkey_alone_is_proper_obstruction"],
        "candidate_grid_confirms_unique_pair": theorem_export["candidate_grid_only_strict_pair_accepted_by_both"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "beta_normalization_left_anchor_source_exported",
            "right_endpoint_value_anchor_source_exported",
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
        "packet_id": "P2520",
        "stage_id": "S1470",
        "status": "STRICT_DAMPING_ENDPOINT_ANCHOR_SUBKEY_LATTICE_CERTIFICATE_SOURCE_OBLIGATION_REFINEMENT_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_endpoint_anchor_subkey_lattice_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_endpoint_anchor_subkey_lattice_certificate"]["theorem_export"]
    lines = [
        "# P2520/S1470 strict damping endpoint-anchor subkey lattice certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2519 endpoint-anchor acceptance inherited: `{t['p2519_endpoint_anchor_acceptance_inherited']}`.",
        f"- Subkey lattice row count: `{t['subkey_lattice_row_count']}`.",
        f"- Accepted row count: `{t['accepted_row_count']}`.",
        f"- Proper subkey obstruction count: `{t['proper_subkey_obstruction_count']}`.",
        f"- Only both subkeys accept: `{t['only_both_subkeys_accept']}`.",
        f"- Left subkey alone is obstruction: `{t['left_subkey_alone_is_proper_obstruction']}`.",
        f"- Right subkey alone is obstruction: `{t['right_subkey_alone_is_proper_obstruction']}`.",
        f"- Candidate grid rows: `{t['candidate_grid_row_count']}`.",
        f"- Candidate grid only strict pair accepted by both: `{t['candidate_grid_only_strict_pair_accepted_by_both']}`.",
        f"- Numeric source refinement exported: `{t['numeric_source_refinement_exported']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet refines the source obligation only. It does not source either endpoint subkey, does not source the m=2 operator signature, and exports no bridge completion, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_endpoint_anchor_subkey_lattice_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_endpoint_anchor_subkey_lattice_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
