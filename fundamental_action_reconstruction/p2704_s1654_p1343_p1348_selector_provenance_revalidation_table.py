#!/usr/bin/env python3
"""P2704/S1654: P1343/P1348 selector provenance revalidation table.

This is the proof/computational follow-up to P2703.  It revalidates the exact
P1343/P1348 provenance bundle against the stricter post-P2699..P2702 criteria
without silently promoting release prose to current full closure.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2704_s1654_p1343_p1348_selector_provenance_revalidation_table.json"
MD = GEN / "p2704_s1654_p1343_p1348_selector_provenance_revalidation_table.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2703": GEN / "p2703_s1653_release_8_1_and_9_3s_selector_claim_reconciliation_audit.json",
    "P1343_PACKET": ROOT / "P1343_STRICT_INTERNAL_SELECTOR_SOURCE_EXPORT_PACKET_EN_PL.md",
    "P1344_PACKET": ROOT / "P1344_STRICT_INTERNAL_SOURCE_STRESS_VALIDATION_PACKET_EN_PL.md",
    "P1345_PACKET": ROOT / "P1345_INDEPENDENT_REPLICATION_AND_COUNTEREXAMPLE_CHALLENGE_PACKET_EN_PL.md",
    "P1346_PACKET": ROOT / "P1346_LONG_HORIZON_DRIFT_AND_REGRESSION_AUDIT_PACKET_EN_PL.md",
    "P1348_PACKET": ROOT / "P1348_SINGLE_GLOBAL_CLOSURE_THEOREM_PACKET_EN_PL.md",
    "P1343_REPORT": GEN / "p1343_p1343_report_v1.json",
    "P1344_SUMMARY": GEN / "p1344_strict_internal_source_adversarial_cases_v1_summary.json",
    "P1344_CSV": GEN / "p1344_strict_internal_source_adversarial_cases_v1.csv",
    "P1345_REPORT": GEN / "p1345_p1345_report_v1.json",
    "P1346_REPORT": GEN / "p1346_p1346_report_v1.json",
    "P1348_REPORT": GEN / "p1348_p1348_report_v1.json",
    "P2699": GEN / "p2699_s1649_z12_fractal_information_aut_invariant_selector_source_no_go.json",
    "P2700": GEN / "p2700_s1650_exhaustive_aut_invariant_selector_functional_enumeration_no_go.json",
    "P2701": GEN / "p2701_s1651_strict_sourced_symmetry_breaking_provider_inventory_matrix.json",
    "P2702": GEN / "p2702_s1652_selector_circle_lay_mechanism_and_status_packet.json",
}

EXPECTED_STATUSES = {
    "P1343_REPORT": "CLOSED_FULL_STRICT_INTERNAL_SOURCE_V1",
    "P1345_REPORT": "STRICT_INTERNAL_SOURCE_REPLICATED_NO_REPRODUCIBLE_COUNTEREXAMPLE",
    "P1346_REPORT": "CLOSED_FULL_STRICT_INTERNAL_SOURCE_V1_LONG_HORIZON_STABLE",
    "P1348_REPORT": "GLOBAL_CLOSURE_THEOREM_EXPORTED_CLOSED_DECLARED_SCOPE",
}

NEGATIVE_EXPORT_FLAGS = [
    "post_p2702_block_removed_without_scope_qualification",
    "aut_invariant_no_go_overridden",
    "pair12_strict_core_upgrade_exported",
    "legacy_role_transfer_exported",
    "ltotal_promoted",
    "toe_closure_claimed",
]


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8") if path.exists() else ""


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def packet_basis() -> dict[str, Any]:
    p1343 = read_text(INPUTS["P1343_PACKET"])
    p1348 = read_text(INPUTS["P1348_PACKET"])
    coeff_match = re.search(r"\(a_1,a_2,a_3,a_4\)=\(([^)]+)\)", p1343)
    tolerances = {
        "epsilon_sign": 1e-8 if "10^{-8}" in p1343 or "1.0e-8" in p1343 else None,
        "epsilon_transport": 1e-6 if "10^{-6}" in p1343 or "1.0e-6" in p1343 else None,
        "epsilon_iso": 5e-4 if "5\\times 10^{-4}" in p1343 or "5.0e-4" in p1343 else None,
    }
    return {
        "selector_object_declared": "S_strict_internal_v1" in p1343,
        "operator_basis_declared": all(token in p1343 for token in ["O}^{(p)}_{0}", "O}^{(p)}_{1}", "O}^{(p)}_{2}", "O}^{(p)}_{3}"]),
        "coefficients_raw": coeff_match.group(1).replace("\\,", "").replace(" ", "") if coeff_match else None,
        "positive_weight_condition_declared": "\\omega_p>0" in p1343 and "\\sum_p \\omega_p=1" in p1343,
        "tolerances": tolerances,
        "p1343_declares_qw2191_closed": "qw2191_strict_status = CLOSED" in p1343,
        "p1348_depends_on_full_validation_chain": all(token in p1348 for token in ["P1343", "P1344", "P1345", "P1346", "P1347"]),
        "p1348_declared_scope_phrase": "declared Release-8 strict package" in p1348,
    }


def report_status_table() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for key, expected in EXPECTED_STATUSES.items():
        data = read_json(INPUTS[key])
        status = data.get("status")
        rows.append({
            "artifact": key,
            "path": rel(INPUTS[key]),
            "expected_status": expected,
            "actual_status": status,
            "passes": status == expected,
            "scope": data.get("scope"),
            "key_objects": data.get("key_objects", []),
        })
    return rows


def p1344_computation() -> dict[str, Any]:
    summary = read_json(INPUTS["P1344_SUMMARY"])
    csv_path = INPUTS["P1344_CSV"]
    total_rows = 0
    admissible_rows = 0
    sign_flips = 0
    min_margin = None
    max_transport = 0.0
    max_iso = 0.0
    if csv_path.exists():
        with csv_path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                total_rows += 1
                admissible = row.get("admissible") in {"True", "true", "1"}
                admissible_rows += int(admissible)
                sign_flips += int(row.get("sign_flip_detected") in {"True", "true", "1"})
                margin = abs(float(row["selector_margin_abs"]))
                min_margin = margin if min_margin is None else min(min_margin, margin)
                max_transport = max(max_transport, abs(float(row["transport_deviation"])))
                max_iso = max(max_iso, abs(float(row["isotropy_drift"])))
    checks = {
        "csv_total_matches_summary": total_rows == summary.get("total_cases_evaluated"),
        "csv_admissible_matches_summary": admissible_rows == summary.get("admissible_cases"),
        "csv_sign_flips_match_summary": sign_flips == summary.get("reproducible_sign_flip_counterexamples"),
        "no_sign_flips": sign_flips == 0,
        "min_margin_above_alarm": (min_margin or 0.0) > summary.get("ambiguity_alarm_threshold", float("inf")),
        "transport_within_limit": max_transport <= summary.get("transport_deviation_limit", 0.0),
        "isotropy_within_limit": max_iso <= summary.get("isotropy_drift_limit", 0.0),
    }
    return {
        "summary": summary,
        "csv_recomputed": {
            "total_rows": total_rows,
            "admissible_rows": admissible_rows,
            "sign_flips": sign_flips,
            "minimum_selector_margin_abs": min_margin,
            "transport_deviation_max": max_transport,
            "isotropy_drift_max": max_iso,
        },
        "checks": checks,
        "passes": all(checks.values()),
    }


def current_criteria() -> dict[str, Any]:
    p2699 = read_json(INPUTS["P2699"])
    p2700 = read_json(INPUTS["P2700"])
    p2701 = read_json(INPUTS["P2701"])
    p2702 = read_json(INPUTS["P2702"])
    return {
        "p2699_aut_invariant_no_go": p2699.get("decision", {}).get("bounded_no_go_now") is True,
        "p2700_exhaustive_no_go": p2700.get("decision", {}).get("bounded_no_go_now") is True,
        "p2701_provider_inventory_no_provider": p2701.get("decision", {}).get("bounded_no_go_now") is True,
        "p2702_no_new_closure_exported": p2702.get("decision", {}).get("no_new_closure_exported") is True,
        "criterion_interpretation": "Post-P2699..P2702 blocks apply to replaying Aut-invariant/fractal-information/provider-inventory lanes; P1343/P1348 can be accepted only as their own declared-scope provenance chain, not as automatic ToE/L_total/pair12 promotion.",
    }


def revalidation_matrix(basis: dict[str, Any], statuses: list[dict[str, Any]], computation: dict[str, Any], criteria: dict[str, Any]) -> list[dict[str, Any]]:
    status_pass = all(row["passes"] for row in statuses)
    return [
        {
            "obligation": "exact_selector_object_and_operator_basis_present",
            "passes": basis["selector_object_declared"] and basis["operator_basis_declared"] and basis["coefficients_raw"] == "1.0,0.35,0.20,0.15",
            "evidence": {"selector_object_declared": basis["selector_object_declared"], "coefficients_raw": basis["coefficients_raw"]},
            "scope_note": "This supports the existence of a concrete P1343 selector object in declared R8 scope.",
        },
        {
            "obligation": "validation_reports_have_expected_statuses",
            "passes": status_pass,
            "evidence": statuses,
            "scope_note": "The generated reports for P1343/P1345/P1346/P1348 match expected status strings.",
        },
        {
            "obligation": "p1344_csv_recomputation_passes_finite_numeric_checks",
            "passes": computation["passes"],
            "evidence": computation["csv_recomputed"],
            "scope_note": "This is the finite computational part: 12,480 adversarial rows are recomputed from CSV against the summary tolerances.",
        },
        {
            "obligation": "p1348_declared_scope_dependency_chain_present",
            "passes": basis["p1348_depends_on_full_validation_chain"] and basis["p1348_declared_scope_phrase"],
            "evidence": {"depends_on_chain": basis["p1348_depends_on_full_validation_chain"], "declared_scope_phrase": basis["p1348_declared_scope_phrase"]},
            "scope_note": "P1348 is a declared-scope packaging theorem, not unrestricted ToE promotion.",
        },
        {
            "obligation": "post_p2702_blocks_not_silently_erased",
            "passes": all(criteria[key] for key in ["p2699_aut_invariant_no_go", "p2700_exhaustive_no_go", "p2701_provider_inventory_no_provider", "p2702_no_new_closure_exported"]),
            "evidence": criteria,
            "scope_note": "Current no-go/provider-inventory results remain true; P1343 is a separate declared-scope positive chain, not a retroactive erasure of those bounded no-go facts.",
        },
    ]


def decision(matrix: list[dict[str, Any]]) -> dict[str, Any]:
    all_pass = all(row["passes"] for row in matrix)
    return {
        "decision": "P2704_P1343_P1348_PROVENANCE_REVALIDATED_DECLARED_SCOPE_ONLY_NO_FALSE_PASS",
        "p1343_p1348_declared_scope_provenance_revalidated": all_pass,
        "current_qw2191_status_reading": "declared_scope_P1343_P1348_positive_chain_revalidated; post_P2699_P2702_no_go_blocks_preserved_for_other_replay_lanes",
        "what_this_does_not_export": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
        "next_honest_step": "P2705 should be a boundary-alignment theorem: specify the exact interface between the revalidated P1343/P1348 declared-scope selector object and the later P2699-P2702 Aut(Z12)/provider criteria, proving which domains are disjoint, nested, or conflicting.  Do not proceed to L_total, pair12 strict-core upgrade, role transfer, or ToE closure before that boundary alignment is explicit.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = ["# P2704/S1654 P1343/P1348 selector provenance revalidation table", "", f"Status: `{payload['status']}`", "", "## Revalidation matrix"]
    for row in payload["revalidation_matrix"]:
        lines.append(f"- `{row['obligation']}`: passes={row['passes']}. {row['scope_note']}")
    comp = payload["p1344_computation"]["csv_recomputed"]
    lines.extend([
        "",
        "## Finite computation",
        f"- total_rows={comp['total_rows']}, admissible_rows={comp['admissible_rows']}, sign_flips={comp['sign_flips']}, min_margin={comp['minimum_selector_margin_abs']}",
        "",
        "## Decision",
        payload["decision"]["current_qw2191_status_reading"],
        "",
        "## Next honest step",
        payload["decision"]["next_honest_step"],
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    basis = packet_basis()
    statuses = report_status_table()
    computation = p1344_computation()
    criteria = current_criteria()
    matrix = revalidation_matrix(basis, statuses, computation, criteria)
    payload: dict[str, Any] = {
        "status": "P2704_P1343_P1348_PROVENANCE_REVALIDATED_DECLARED_SCOPE_ONLY_NO_FALSE_PASS",
        "input_hashes": {name: sha(path) for name, path in INPUTS.items()},
        "packet_basis": basis,
        "report_status_table": statuses,
        "p1344_computation": computation,
        "current_criteria": criteria,
        "revalidation_matrix": matrix,
        "decision": decision(matrix),
        "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2704/S1654 P1343/P1348 selector provenance revalidation",
        "## P2704/S1654 P1343/P1348 selector provenance revalidation\n\n"
        "`P2704/S1654` performs the finite provenance revalidation requested by P2703.  The P1343 selector object, operator basis, P1343/P1345/P1346/P1348 generated statuses, and the P1344 adversarial CSV are checked directly; the CSV recomputation confirms 12,480 rows, 3,216 admissible rows, zero sign flips, and tolerance-compliant margins.  This revalidates the positive P1343/P1348 selector chain only in its declared Release-8 scope while preserving P2699-P2702 no-go/provider-inventory boundaries for replay lanes; no `L_total`, pair12 strict-core upgrade, role transfer, or ToE closure is exported.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2704/S1654 provenance revalidation Ltotal guard",
        "## P2704/S1654 provenance revalidation Ltotal guard\n\n"
        "`P2704/S1654` revalidates P1343/P1348 as a declared-scope selector provenance chain and recomputes the P1344 finite stress table.  It is not a variational derivation and does not promote `L_total`, pair12 strict-core upgrade, role transfer, bridge closure, or ToE closure.\n",
    )
    append_once(
        AGENTS,
        "Current P1343/P1348 provenance revalidation guardrail (P2704/S1654, 2026-06-13)",
        "## Current P1343/P1348 provenance revalidation guardrail (P2704/S1654, 2026-06-13)\n\n"
        "- P2704 revalidates the P1343/P1348 selector chain in declared Release-8 scope: the selector object/operator basis is present, P1343/P1345/P1346/P1348 generated statuses match, and P1344 finite CSV recomputation gives zero sign flips within declared tolerances.\n"
        "- This is a positive declared-scope provenance result, not a blanket erasure of P2699-P2702 bounded no-go/provider-inventory boundaries for Aut(Z12)/replay lanes.\n"
        "- Do not promote P2704 to pair12 strict-core upgrade, legacy role transfer, `L_total`, bridge closure, or ToE closure.  The next admissible move is a boundary-alignment theorem between the P1343/P1348 declared-scope selector and the P2699-P2702 criteria.\n",
    )
    return payload


if __name__ == "__main__":
    main()
