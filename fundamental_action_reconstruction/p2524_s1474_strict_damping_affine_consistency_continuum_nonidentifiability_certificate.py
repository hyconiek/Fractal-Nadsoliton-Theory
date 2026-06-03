#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from fractions import Fraction
from itertools import combinations
from statistics import mean, pstdev
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
OUT = GEN / "p2524_s1474_strict_damping_affine_consistency_continuum_nonidentifiability_certificate.json"
MD = GEN / "p2524_s1474_strict_damping_affine_consistency_continuum_nonidentifiability_certificate.md"

SOURCE_FILES = {
    "P2523_PAIRWISE_SECANT_CONSISTENCY": GEN / "p2523_s1473_strict_damping_pairwise_secant_consistency_certificate.json",
}

NODE_DOMAIN = list(range(1, 12))
INTERCEPT_CANDIDATES = [Fraction(-1, 1), Fraction(-1, 2), Fraction(0, 1), Fraction(1, 2), Fraction(1, 1)]
SLOPE_CANDIDATES = [Fraction(-1, 1), Fraction(0, 1), Fraction(1, 2), Fraction(4, 5), Fraction(1, 1), Fraction(9, 5), Fraction(2, 1)]
STRICT_INTERCEPT = Fraction(0, 1)
STRICT_DELTA = Fraction(4, 5)
STRICT_ETA = Fraction(9, 5)
TOL = 1e-14


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
        "new_packet": "P2524|S1474|affine consistency continuum nonidentifiability|secant consistency nonidentifiability|affine-consistency continuum|basis-independent nonidentifiability|secant continuum",
        "precursor_packets": "P2523|S1473|pairwise secant consistency|P2522|two-node anchor basis equivalence|P2521|single-node anchor equivalence",
        "continuum_language": "affine continuum|secant continuum|intercept slope continuum|constant secant|triangle cocycle zero|basis-independent consistency",
        "source_blockers": "beta_eta_numeric_source|node data source|anchor basis source|strict source theorem|proper subset obstruction",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def frac_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def ell(d: int) -> float:
    return math.log(d)


def y_affine(d: int, intercept: Fraction, slope: Fraction) -> float:
    return float(intercept) + float(slope) * ell(d)


def affine_consistency_rows() -> list[dict[str, Any]]:
    rows = []
    for intercept in INTERCEPT_CANDIDATES:
        for slope in SLOPE_CANDIDATES:
            secants = []
            for d_i, d_j in combinations(NODE_DOMAIN, 2):
                secants.append((y_affine(d_j, intercept, slope) - y_affine(d_i, intercept, slope)) / (ell(d_j) - ell(d_i)))
            triangle_cocycles = []
            for d_i, d_j, d_k in combinations(NODE_DOMAIN, 3):
                y_i = y_affine(d_i, intercept, slope)
                y_j = y_affine(d_j, intercept, slope)
                y_k = y_affine(d_k, intercept, slope)
                triangle_cocycles.append((y_j - y_i) + (y_k - y_j) - (y_k - y_i))
            secant_spread = max(secants) - min(secants)
            max_abs_cocycle = max(abs(value) for value in triangle_cocycles)
            rows.append({
                "intercept_log_beta": frac_text(intercept),
                "slope_delta": frac_text(slope),
                "eta_if_slope_delta": frac_text(Fraction(1, 1) + slope),
                "is_strict_numeric_target": intercept == STRICT_INTERCEPT and slope == STRICT_DELTA,
                "pair_count": len(secants),
                "triangle_count": len(triangle_cocycles),
                "secant_min": min(secants),
                "secant_max": max(secants),
                "secant_mean": mean(secants),
                "secant_population_std": pstdev(secants),
                "secant_spread": secant_spread,
                "all_secants_equal_candidate_slope": secant_spread < TOL and all(abs(value - float(slope)) < TOL for value in secants),
                "max_abs_triangle_additive_cocycle": max_abs_cocycle,
                "all_triangle_additive_cocycles_zero": max_abs_cocycle < TOL,
                "basis_independent_affine_consistency_accepts": secant_spread < TOL and max_abs_cocycle < TOL,
            })
    return rows


def endpoint_anchor_filter(rows: list[dict[str, Any]]) -> dict[str, Any]:
    accepted_by_left = [row for row in rows if row["intercept_log_beta"] == "0"]
    accepted_by_slope = [row for row in rows if row["slope_delta"] == "4/5"]
    accepted_by_both = [row for row in rows if row["intercept_log_beta"] == "0" and row["slope_delta"] == "4/5"]
    return {
        "affine_consistency_accepting_rows": sum(1 for row in rows if row["basis_independent_affine_consistency_accepts"]),
        "left_normalization_filter_count": len(accepted_by_left),
        "slope_value_filter_count": len(accepted_by_slope),
        "both_numeric_filters_count": len(accepted_by_both),
        "both_numeric_filters_unique_strict_target": len(accepted_by_both) == 1 and accepted_by_both[0]["is_strict_numeric_target"],
        "left_normalization_alone_leaves_slope_continuum_on_grid": len(accepted_by_left) == len(SLOPE_CANDIDATES),
        "slope_value_alone_leaves_intercept_continuum_on_grid": len(accepted_by_slope) == len(INTERCEPT_CANDIDATES),
    }


def build_continuum_certificate(p2523: dict[str, Any]) -> dict[str, Any]:
    rows = affine_consistency_rows()
    filters = endpoint_anchor_filter(rows)
    accepted = [row for row in rows if row["basis_independent_affine_consistency_accepts"]]
    strict_rows = [row for row in rows if row["is_strict_numeric_target"]]
    return {
        "frontier_atom_under_attack": "basis-independent affine consistency cannot source beta_eta numeric values",
        "p2523_basis_independent_consistency_inherited": p2523.get("basis_independent_affine_consistency_exported") is True,
        "certificate_type": "affine-consistency continuum nonidentifiability certificate",
        "symbolic_statement": "The pairwise-secant/triangle-cocycle consistency predicates characterize affine node data y_d=b+a log d. They do not identify b=0 or a=4/5: every affine intercept/slope pair has constant secants and zero additive triangle cocycles. Numeric beta/eta selection therefore requires independent value/source data beyond affine consistency.",
        "affine_consistency_rows": rows,
        "endpoint_anchor_filter_audit": filters,
        "candidate_grid_row_count": len(rows),
        "affine_consistency_accepting_row_count": len(accepted),
        "strict_numeric_target_row_count": len(strict_rows),
        "all_candidate_grid_rows_pass_affine_consistency": len(accepted) == len(rows),
        "strict_target_is_one_member_of_affine_consistency_continuum": len(strict_rows) == 1 and len(accepted) > 1,
        "basis_independent_affine_consistency_nonidentifiability_exported": True,
        "node_data_source_exported": False,
        "left_normalization_source_exported": False,
        "slope_value_source_exported": False,
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
## P2524/S1474 affine-consistency continuum nonidentifiability certificate

`P2524/S1474` is the converse guard for P2523.  The P2523 pairwise-secant and triangle-cocycle predicates certify that supplied node data are affine in `ell=log d`, but affine consistency alone leaves the continuum `y_d=b+a log d`.  The finite grid audit checks 35 intercept/slope candidates: every affine candidate has constant pairwise secants and zero additive triangle cocycles, while the strict target `(b,a)=(0,4/5)` is only one accepted member.

Thus basis-independent affine consistency is a necessary consistency check but not a numeric source theorem.  It does not export `beta_eta_numeric_source`, does not source the `m=2` operator signature, and exports no strict damping source closure, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, physical-value generator, or ToE closure.
"""
    lag_section = """
## P2524/S1474 affine-consistency continuum guard

`P2524/S1474` records that constant pairwise secants and zero triangle cocycles only prove affine node-data form, not the strict values `b=0` and `a=4/5`.  A role-bearing damping source must still supply the numeric intercept/slope data, so no nonlinear compression-flow source or role-bearing `L_total` term is licensed.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2524/S1474 affine-consistency continuum nonidentifiability certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2524/S1474 affine-consistency continuum guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2523 = theorem(sources["P2523_PAIRWISE_SECANT_CONSISTENCY"], "strict_damping_pairwise_secant_consistency_certificate")
    cert = build_continuum_certificate(p2523)
    filters = cert["endpoint_anchor_filter_audit"]
    theorem_export = {
        "theorem_name": "P2524_T1_strict_damping_affine_consistency_continuum_nonidentifiability_certificate",
        "audited_chain": ["P2523/S1473", "P2522/S1472", "P2521/S1471"],
        "strict_damping_affine_consistency_continuum_nonidentifiability_certificate": cert,
        "frontier_atom_under_attack": cert["frontier_atom_under_attack"],
        "p2523_basis_independent_consistency_inherited": cert["p2523_basis_independent_consistency_inherited"],
        "candidate_grid_row_count": cert["candidate_grid_row_count"],
        "affine_consistency_accepting_row_count": cert["affine_consistency_accepting_row_count"],
        "strict_numeric_target_row_count": cert["strict_numeric_target_row_count"],
        "all_candidate_grid_rows_pass_affine_consistency": cert["all_candidate_grid_rows_pass_affine_consistency"],
        "strict_target_is_one_member_of_affine_consistency_continuum": cert["strict_target_is_one_member_of_affine_consistency_continuum"],
        "left_normalization_alone_leaves_slope_continuum_on_grid": filters["left_normalization_alone_leaves_slope_continuum_on_grid"],
        "slope_value_alone_leaves_intercept_continuum_on_grid": filters["slope_value_alone_leaves_intercept_continuum_on_grid"],
        "both_numeric_filters_unique_strict_target": filters["both_numeric_filters_unique_strict_target"],
        "basis_independent_affine_consistency_nonidentifiability_exported": cert["basis_independent_affine_consistency_nonidentifiability_exported"],
        "node_data_source_exported": False,
        "left_normalization_source_exported": False,
        "slope_value_source_exported": False,
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
            "P2524 exports nonidentifiability of affine-consistency predicates, not a source for beta/eta numeric values.",
            "Constant secants and zero triangle cocycles leave an intercept/slope continuum.",
            "The m=2 operator signature source remains separate and unsourced by this packet.",
            "No damping bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
        "next_honest_step": "Search for strict dynamics that selects b=0 and a=4/5 inside the affine-consistency continuum, rather than checking affine consistency after numeric data are supplied.",
    }
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2523_inherited": theorem_export["p2523_basis_independent_consistency_inherited"],
        "grid_shape_ok": theorem_export["candidate_grid_row_count"] == len(INTERCEPT_CANDIDATES) * len(SLOPE_CANDIDATES),
        "all_grid_rows_pass_affine_consistency": theorem_export["all_candidate_grid_rows_pass_affine_consistency"],
        "strict_target_not_unique_under_affine_consistency": theorem_export["strict_target_is_one_member_of_affine_consistency_continuum"],
        "proper_filters_leave_continua": theorem_export["left_normalization_alone_leaves_slope_continuum_on_grid"] and theorem_export["slope_value_alone_leaves_intercept_continuum_on_grid"],
        "negative_controls_preserved": not any(theorem_export[key] for key in [
            "node_data_source_exported",
            "left_normalization_source_exported",
            "slope_value_source_exported",
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
        "packet_id": "P2524",
        "stage_id": "S1474",
        "status": "STRICT_DAMPING_AFFINE_CONSISTENCY_CONTINUUM_NONIDENTIFIABILITY_CERTIFICATE_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_affine_consistency_continuum_nonidentifiability_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_affine_consistency_continuum_nonidentifiability_certificate"]["theorem_export"]
    lines = [
        "# P2524/S1474 strict damping affine-consistency continuum nonidentifiability certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2523 basis-independent consistency inherited: `{t['p2523_basis_independent_consistency_inherited']}`.",
        f"- Candidate grid row count: `{t['candidate_grid_row_count']}`.",
        f"- Affine-consistency accepting row count: `{t['affine_consistency_accepting_row_count']}`.",
        f"- Strict numeric target row count: `{t['strict_numeric_target_row_count']}`.",
        f"- All grid rows pass affine consistency: `{t['all_candidate_grid_rows_pass_affine_consistency']}`.",
        f"- Strict target is one continuum member: `{t['strict_target_is_one_member_of_affine_consistency_continuum']}`.",
        f"- Both numeric filters unique strict target: `{t['both_numeric_filters_unique_strict_target']}`.",
        f"- Strict damping source exported: `{t['strict_damping_beta_eta_source_exported']}`.",
        "",
        "## Negative controls",
        "",
        "This packet exports only nonidentifiability of affine-consistency predicates. It does not source node data, left normalization, slope value, m=2 operator signature, bridge completion, role-transfer theorem, selector/QW-2191 closure, role-bearing L_total, physical-value generator, or ToE closure.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_affine_consistency_continuum_nonidentifiability_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload["strict_damping_affine_consistency_continuum_nonidentifiability_certificate"]["theorem_export"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
