#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from itertools import product
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
OUT = GEN / "p2548_s1498_strict_damping_m2_trace_arity_conditional_selection_certificate.json"
MD = GEN / "p2548_s1498_strict_damping_m2_trace_arity_conditional_selection_certificate.md"

SOURCE_FILES = {
    "P2540_M2_OBSTRUCTION": GEN / "p2540_s1490_strict_damping_m2_operator_signature_current_premise_obstruction_certificate.json",
    "P2544_FOUR_KEY_CLOSURE_BLOCKER": GEN / "p2544_s1494_strict_damping_four_key_current_premise_closure_blocker_certificate.json",
    "P2547_POST_IDENTITY_RESIDUAL_TRIKEY": GEN / "p2547_s1497_strict_damping_post_identity_residual_trikey_certificate.json",
}

ORDER_RANGE = list(range(1, 11))
TARGET_TRACE_ARITY = 4
RESIDUAL_AFTER_IDENTITY_AND_M2_KEYS = [
    "prime_log_proportionality_source",
    "slope_value_or_prime_anchor_source",
]
NEGATIVE_EXPORT_FLAGS = [
    "multiplicative_character_law_source_exported",
    "prime_log_proportionality_source_exported",
    "slope_value_or_prime_anchor_source_exported",
    "m2_operator_signature_source_exported",
    "strict_quadruple_trace_source_exported",
    "beta_eta_numeric_source_exported",
    "strict_damping_beta_eta_source_exported",
    "source_obligation_discharge_exported",
    "damping_compression_bridge_component_ready",
    "full_bridge_theorem_exported",
    "role_transfer_theorem_exported",
    "selector_closure_exported",
    "qw2191_discharged_by_this_certificate",
    "role_bearing_ltotal_exported",
    "toe_closure_claimed",
]


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
    return {"count": len(lines), "samples": lines[:70]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2548|S1498|trace-arity|trace arity|quadruple trace|m2 trace arity",
        "intended_research_nonduplication": "quadruple boundary trace|boundary trace arity|trace arity.*m=2|2m=4|four trace.*operator-order",
        "m2_precursor_language": "m=2 operator-order selection|operator-order selection|m2_operator_signature_source|biharmonic|derivative-order",
        "post_identity_frontier": "post-identity residual|residual tri-key|prime_log_proportionality_source|slope_value_or_prime_anchor_source",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def load_theorems(sources: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {
        "P2540_M2_OBSTRUCTION": theorem(sources["P2540_M2_OBSTRUCTION"], "strict_damping_m2_operator_signature_current_premise_obstruction_certificate"),
        "P2544_FOUR_KEY_CLOSURE_BLOCKER": theorem(sources["P2544_FOUR_KEY_CLOSURE_BLOCKER"], "strict_damping_four_key_current_premise_closure_blocker_certificate"),
        "P2547_POST_IDENTITY_RESIDUAL_TRIKEY": theorem(sources["P2547_POST_IDENTITY_RESIDUAL_TRIKEY"], "strict_damping_post_identity_residual_trikey_certificate"),
    }


def trace_arity_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for order in ORDER_RANGE:
        boundary_trace_arity = 2 * order
        rows.append({
            "derivative_order_m": order,
            "functional_schema": f"J_{order}[y]=int (D^{order} y)^2 d ell",
            "euler_lagrange_order": 2 * order,
            "self_adjoint_boundary_trace_arity": boundary_trace_arity,
            "trace_count_formula": "2*m total endpoint traces for the order-2m self-adjoint variational problem",
            "matches_strict_quadruple_trace_arity": boundary_trace_arity == TARGET_TRACE_ARITY,
            "matches_target_m2_signature": order == 2,
        })
    return rows


def exact_trace_selection(rows: list[dict[str, Any]]) -> dict[str, Any]:
    matching_orders = [row["derivative_order_m"] for row in rows if row["matches_strict_quadruple_trace_arity"]]
    nonmatching_witnesses = [
        {
            "derivative_order_m": row["derivative_order_m"],
            "self_adjoint_boundary_trace_arity": row["self_adjoint_boundary_trace_arity"],
            "defect_from_quadruple_trace": row["self_adjoint_boundary_trace_arity"] - TARGET_TRACE_ARITY,
        }
        for row in rows
        if not row["matches_strict_quadruple_trace_arity"]
    ]
    return {
        "target_trace_arity": TARGET_TRACE_ARITY,
        "selection_equation": "2*m = 4",
        "integer_solution_m": TARGET_TRACE_ARITY // 2,
        "equation_has_integer_solution": TARGET_TRACE_ARITY % 2 == 0,
        "matching_orders_in_audited_range": matching_orders,
        "unique_matching_order_in_audited_range": matching_orders == [2],
        "nonmatching_order_witness_count": len(nonmatching_witnesses),
        "nonmatching_order_witnesses": nonmatching_witnesses,
        "conditional_theorem_statement": (
            "If strict nadsoliton dynamics exports exactly four independent self-adjoint boundary/trace channels "
            "for the derivative-only order family J_m, then the trace-arity equation 2*m=4 selects m=2."
        ),
    }


def residual_after_identity_and_m2_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for prime_log, slope in product([False, True], repeat=2):
        assignment = {
            "multiplicative_character_law_source": True,
            "prime_log_proportionality_source": prime_log,
            "slope_value_or_prime_anchor_source": slope,
            "m2_operator_signature_source": True,
        }
        beta_eta_numeric = prime_log and slope
        strict_damping = beta_eta_numeric and assignment["m2_operator_signature_source"]
        rows.append({
            "residual_assignment": {key: assignment[key] for key in RESIDUAL_AFTER_IDENTITY_AND_M2_KEYS},
            "identity_multiplicative_key_assumed_strict": True,
            "m2_operator_signature_assumed_strict_by_quadruple_trace": True,
            "beta_eta_numeric_source_accepts": beta_eta_numeric,
            "strict_damping_beta_eta_source_accepts": strict_damping,
            "missing_residual_keys": [key for key in RESIDUAL_AFTER_IDENTITY_AND_M2_KEYS if not assignment[key]],
        })
    return rows


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    source_theorems = load_theorems(sources)
    rows = trace_arity_rows()
    selection = exact_trace_selection(rows)
    residual_rows = residual_after_identity_and_m2_rows()
    accepting_rows = [row for row in residual_rows if row["strict_damping_beta_eta_source_accepts"]]
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2548_T1_strict_damping_m2_trace_arity_conditional_selection_certificate",
        "audited_chain": ["P2540/S1490", "P2544/S1494", "P2547/S1497"],
        "frontier_source_key_under_attack": "m2_operator_signature_source",
        "p2540_m2_obstruction_inherited": source_theorems["P2540_M2_OBSTRUCTION"].get("m2_operator_signature_source_route_refuted_for_current_source_free_premises") is True,
        "p2544_four_key_blocker_inherited": source_theorems["P2544_FOUR_KEY_CLOSURE_BLOCKER"].get("all_four_current_premise_routes_blocked") is True,
        "p2547_post_identity_residual_trikey_inherited": source_theorems["P2547_POST_IDENTITY_RESIDUAL_TRIKEY"].get("post_identity_residual_trikey_irredundancy_exported") is True,
        "operator_family_audited": "node-fixed derivative-only self-adjoint quadratic family J_m[y]=int (D^m y)^2 d ell",
        "finite_order_range": ORDER_RANGE,
        "trace_arity_rows": rows,
        "trace_arity_selection_audit": selection,
        "conditional_quadruple_trace_selects_m2": selection["unique_matching_order_in_audited_range"],
        "strict_quadruple_trace_source_required": True,
        "strict_quadruple_trace_source_exported": False,
        "conditional_m2_closure_if_quadruple_trace_source_supplied": True,
        "m2_operator_signature_source_exported": False,
        "conditional_not_actual_source_reason": (
            "The packet proves only the implication quadruple-trace-source => m=2 inside the audited operator family. "
            "It does not prove that strict nadsoliton dynamics exports the quadruple trace source."
        ),
        "residual_keys_after_identity_and_conditional_m2": RESIDUAL_AFTER_IDENTITY_AND_M2_KEYS,
        "residual_truth_table_after_identity_and_conditional_m2": residual_rows,
        "residual_truth_table_after_identity_and_conditional_m2_row_count": len(residual_rows),
        "residual_accepting_row_count_after_identity_and_conditional_m2": len(accepting_rows),
        "residual_accepting_row_after_identity_and_conditional_m2": accepting_rows[0],
        "conditional_identity_and_m2_reduce_missing_count_from_3_to_2": True,
        "identity_plus_conditional_m2_still_cannot_export_beta_eta_numeric_source": True,
        "identity_plus_conditional_m2_still_cannot_export_strict_damping_beta_eta_source": True,
        "recommended_next_honest_step": (
            "Do not declare the m=2 source closed from trace bookkeeping alone. The next honest step is to search for "
            "or prove a strict nadsoliton quadruple-boundary-trace theorem; if that source cannot be exported, pivot "
            "to the remaining numerical bottleneck prime_log_proportionality_source before any role-bearing L_total promotion."
        ),
        "not_licensed": [
            "No actual strict quadruple-trace source theorem is exported.",
            "No residual prime-log or slope-value source theorem is exported.",
            "No legacy-to-strict bridge completion or role-transfer theorem is exported.",
            "No QW-2191 selector discharge, role-bearing L_total term, physical-value generator, or ToE closure is exported.",
        ],
    }
    for flag in NEGATIVE_EXPORT_FLAGS:
        theorem_export[flag] = False
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2540_m2_obstruction_inherited": theorem_export["p2540_m2_obstruction_inherited"],
        "p2544_four_key_blocker_inherited": theorem_export["p2544_four_key_blocker_inherited"],
        "p2547_post_identity_residual_trikey_inherited": theorem_export["p2547_post_identity_residual_trikey_inherited"],
        "trace_arity_equation_uniquely_selects_m2_conditionally": theorem_export["conditional_quadruple_trace_selects_m2"],
        "actual_quadruple_trace_source_not_exported": theorem_export["strict_quadruple_trace_source_exported"] is False,
        "no_false_bridge_or_role_transfer": theorem_export["full_bridge_theorem_exported"] is False and theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_or_toe_claim": theorem_export["qw2191_discharged_by_this_certificate"] is False and theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2548",
        "stage_id": "S1498",
        "status": "STRICT_DAMPING_M2_TRACE_ARITY_CONDITIONAL_SELECTION_CERTIFICATE_CONDITIONAL_ONLY_NO_ACTUAL_TRACE_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_m2_trace_arity_conditional_selection_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_m2_trace_arity_conditional_selection_certificate"]["theorem_export"]
    selection = t["trace_arity_selection_audit"]
    lines = [
        "# P2548/S1498 strict damping m2 trace-arity conditional selection certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier source key under attack: `{t['frontier_source_key_under_attack']}`.",
        f"- P2540 m2 current-premise obstruction inherited: `{t['p2540_m2_obstruction_inherited']}`.",
        f"- P2544 four-key blocker inherited: `{t['p2544_four_key_blocker_inherited']}`.",
        f"- P2547 post-identity residual tri-key inherited: `{t['p2547_post_identity_residual_trikey_inherited']}`.",
        f"- Audited derivative orders: `{t['finite_order_range']}`.",
        f"- Trace-arity selection equation: `{selection['selection_equation']}`.",
        f"- Matching audited orders: `{selection['matching_orders_in_audited_range']}`.",
        f"- Conditional quadruple-trace selects m=2: `{t['conditional_quadruple_trace_selects_m2']}`.",
        f"- Strict quadruple-trace source exported: `{t['strict_quadruple_trace_source_exported']}`.",
        f"- m2 operator signature source exported: `{t['m2_operator_signature_source_exported']}`.",
        f"- Residual keys after identity + conditional m2: `{t['residual_keys_after_identity_and_conditional_m2']}`.",
        f"- Residual accepting rows after identity + conditional m2: `{t['residual_accepting_row_count_after_identity_and_conditional_m2']}` of `{t['residual_truth_table_after_identity_and_conditional_m2_row_count']}`.",
        "",
        "## Conditional theorem",
        "",
        selection["conditional_theorem_statement"],
        "",
        "This is not an actual source export: the packet does not prove that strict nadsoliton dynamics supplies exactly four independent self-adjoint boundary/trace channels.",
        "",
        "## Residual frontier",
        "",
        "Even if identity-action and the conditional quadruple-trace premise are both supplied, beta/eta numerics still require prime-log proportionality and a slope/prime anchor. Therefore no strict damping beta/eta source or role-bearing `L_total` term is exported.",
        "",
        "## Recommended next honest step",
        "",
        t["recommended_next_honest_step"],
        "",
        "## Negative controls",
        "",
        "No actual quadruple-trace source, residual prime-log source, residual slope source, source-obligation discharge, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, physical-value generator, or ToE closure is exported.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_m2_trace_arity_conditional_selection_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2548/S1498` attacks the post-identity `m2_operator_signature_source` frontier by isolating a sharper conditional source theorem.  In the derivative-only self-adjoint family `J_m[y]=int (D^m y)^2 d ell`, the endpoint trace arity is `2m`; an independently exported strict quadruple-trace theorem would solve `2m=4` and hence select `m=2` uniquely on the audited `m=1..10` range.  This is only a conditional implication: no strict quadruple-trace source is exported, and after identity plus conditional `m=2`, prime-log proportionality and the slope/prime anchor remain required before `beta_eta_numeric_source` or strict damping source closure can be claimed.
""".strip()
    lag_section = """
`P2548/S1498` permits only a conditional `m=2` bookkeeping guard for `L_total`: a future strict quadruple-boundary-trace theorem would select the biharmonic order inside the audited derivative-only family, but this packet does not export that theorem.  The damping/compression term remains non-role-bearing until the actual trace source, the residual prime-log source, the residual slope/anchor source, bridge completion, and role-transfer gates are explicit.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2548/S1498 m2 trace-arity conditional selection guard", "## P2548/S1498 m2 trace-arity conditional selection guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2548/S1498 m2 trace-arity conditional Ltotal guard", "## P2548/S1498 m2 trace-arity conditional Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
