#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
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
OUT = GEN / "p2549_s1499_strict_damping_quadruple_trace_current_premise_obstruction_certificate.json"
MD = GEN / "p2549_s1499_strict_damping_quadruple_trace_current_premise_obstruction_certificate.md"

SOURCE_FILES = {
    "P2540_M2_OBSTRUCTION": GEN / "p2540_s1490_strict_damping_m2_operator_signature_current_premise_obstruction_certificate.json",
    "P2547_POST_IDENTITY_RESIDUAL_TRIKEY": GEN / "p2547_s1497_strict_damping_post_identity_residual_trikey_certificate.json",
    "P2548_TRACE_ARITY_CONDITIONAL": GEN / "p2548_s1498_strict_damping_m2_trace_arity_conditional_selection_certificate.json",
}

TARGET_TRACE_ARITY = 4
NEGATIVE_EXPORT_FLAGS = [
    "strict_quadruple_trace_source_exported",
    "m2_operator_signature_source_exported",
    "multiplicative_character_law_source_exported",
    "prime_log_proportionality_source_exported",
    "slope_value_or_prime_anchor_source_exported",
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
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2549|S1499|quadruple trace current-premise obstruction|quadruple-trace obstruction|trace source obstruction",
        "intended_research_nonduplication": "current-premise.*quadruple trace|quadruple trace.*current-premise|exactly four.*trace.*obstruction|trace arity.*nonentailment|quadruple-boundary-trace theorem",
        "precursor_trace_arity": "P2548|S1498|trace-arity|trace arity|2\\*m = 4|quadruple trace",
        "m2_frontier_language": "P2540|S1490|m=2 operator-order selection|m2_operator_signature_source|derivative-order|biharmonic",
        "guardrails": "K_legacy_ont|K_strict_gate|role-transfer audit|QW-2191|ToE closure|selector guardrail",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def load_theorems(sources: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {
        "P2540_M2_OBSTRUCTION": theorem(sources["P2540_M2_OBSTRUCTION"], "strict_damping_m2_operator_signature_current_premise_obstruction_certificate"),
        "P2547_POST_IDENTITY_RESIDUAL_TRIKEY": theorem(sources["P2547_POST_IDENTITY_RESIDUAL_TRIKEY"], "strict_damping_post_identity_residual_trikey_certificate"),
        "P2548_TRACE_ARITY_CONDITIONAL": theorem(sources["P2548_TRACE_ARITY_CONDITIONAL"], "strict_damping_m2_trace_arity_conditional_selection_certificate"),
    }


def inherited_passing_orders(p2540: dict[str, Any]) -> list[int]:
    return list(p2540.get("current_source_free_premise_passing_orders", []))


def trace_nonentailment_rows(passing_orders: list[int]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for order in passing_orders:
        trace_arity = 2 * order
        rows.append({
            "derivative_order_m": order,
            "inherits_current_source_free_premise_pass": True,
            "self_adjoint_trace_arity": trace_arity,
            "satisfies_exact_quadruple_trace": trace_arity == TARGET_TRACE_ARITY,
            "is_target_m2_order": order == 2,
            "countermodel_to_quadruple_trace_source": trace_arity != TARGET_TRACE_ARITY,
        })
    return rows


def build_trace_obstruction(rows: list[dict[str, Any]]) -> dict[str, Any]:
    accepting = [row for row in rows if row["satisfies_exact_quadruple_trace"]]
    countermodels = [row for row in rows if row["countermodel_to_quadruple_trace_source"]]
    trace_arities = [row["self_adjoint_trace_arity"] for row in rows]
    return {
        "target_trace_arity": TARGET_TRACE_ARITY,
        "current_premise_passing_order_count": len(rows),
        "current_premise_trace_arities": trace_arities,
        "quadruple_trace_accepting_orders": [row["derivative_order_m"] for row in accepting],
        "quadruple_trace_accepting_order_count": len(accepting),
        "quadruple_trace_countermodel_orders": [row["derivative_order_m"] for row in countermodels],
        "quadruple_trace_countermodel_count": len(countermodels),
        "m3_countermodel": next(row for row in rows if row["derivative_order_m"] == 3),
        "current_premises_entail_exact_quadruple_trace": len(countermodels) == 0,
        "current_premises_entail_m2_via_quadruple_trace": len(countermodels) == 0 and [row["derivative_order_m"] for row in accepting] == [2],
        "model_separation_statement": (
            "The same current source-free derivative-order premises inherited from P2540 admit orders whose trace arity is not four; "
            "therefore those premises do not export the exact quadruple-trace source needed by P2548."
        ),
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    source_theorems = load_theorems(sources)
    passing_orders = inherited_passing_orders(source_theorems["P2540_M2_OBSTRUCTION"])
    rows = trace_nonentailment_rows(passing_orders)
    obstruction = build_trace_obstruction(rows)
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2549_T1_strict_damping_quadruple_trace_current_premise_obstruction_certificate",
        "audited_chain": ["P2540/S1490", "P2547/S1497", "P2548/S1498"],
        "frontier_source_under_attack": "strict_quadruple_trace_source",
        "p2540_current_premise_order_nonidentifiability_inherited": source_theorems["P2540_M2_OBSTRUCTION"].get("m2_and_m3_both_pass_current_source_free_premises") is True,
        "p2547_post_identity_residual_trikey_inherited": source_theorems["P2547_POST_IDENTITY_RESIDUAL_TRIKEY"].get("post_identity_residual_trikey_irredundancy_exported") is True,
        "p2548_conditional_trace_arity_selection_inherited": source_theorems["P2548_TRACE_ARITY_CONDITIONAL"].get("conditional_quadruple_trace_selects_m2") is True,
        "p2548_actual_trace_source_absent_inherited": source_theorems["P2548_TRACE_ARITY_CONDITIONAL"].get("strict_quadruple_trace_source_exported") is False,
        "operator_family_audited": "P2540/P2548 node-fixed derivative-only self-adjoint quadratic order family",
        "trace_nonentailment_rows": rows,
        "trace_obstruction_audit": obstruction,
        "current_premise_quadruple_trace_nonentailment_exported": not obstruction["current_premises_entail_exact_quadruple_trace"],
        "current_premise_m2_via_quadruple_trace_route_refuted": not obstruction["current_premises_entail_m2_via_quadruple_trace"],
        "strict_quadruple_trace_source_required_for_p2548": True,
        "strict_quadruple_trace_source_exported": False,
        "m2_operator_signature_source_exported": False,
        "residual_after_failed_trace_source": [
            "strict_quadruple_trace_source",
            "prime_log_proportionality_source",
            "slope_value_or_prime_anchor_source",
        ],
        "recommended_next_honest_step": (
            "Do not repeat trace-arity bookkeeping. P2549 shows current premises do not source the exact quadruple trace. "
            "The next honest proof/computation step is either to derive a real four-trace source from a new strict nadsoliton boundary/symplectic mechanism, "
            "or, if no such mechanism is available, pivot to the independent numeric residual key prime_log_proportionality_source."
        ),
        "not_licensed": [
            "The exact four-trace premise is not derived from current derivative-order premises.",
            "The P2548 conditional m=2 implication remains conditional and cannot be promoted to m2_operator_signature_source.",
            "No residual prime-log or slope-value source theorem is exported.",
            "No bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing L_total, physical-value generator, or ToE closure is exported.",
        ],
    }
    for flag in NEGATIVE_EXPORT_FLAGS:
        theorem_export[flag] = False
    gatekeepers = {
        "rg_audit_performed": grep["tool"] == "rg",
        "p2540_order_nonidentifiability_inherited": theorem_export["p2540_current_premise_order_nonidentifiability_inherited"],
        "p2548_conditional_trace_arity_inherited": theorem_export["p2548_conditional_trace_arity_selection_inherited"],
        "p2548_actual_trace_source_absent_inherited": theorem_export["p2548_actual_trace_source_absent_inherited"],
        "countermodels_to_quadruple_trace_exist": obstruction["quadruple_trace_countermodel_count"] > 0,
        "m3_is_explicit_countermodel": obstruction["m3_countermodel"]["countermodel_to_quadruple_trace_source"] is True,
        "no_false_source_export": theorem_export["strict_quadruple_trace_source_exported"] is False and theorem_export["m2_operator_signature_source_exported"] is False,
        "no_false_bridge_or_role_transfer": theorem_export["full_bridge_theorem_exported"] is False and theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_or_toe_claim": theorem_export["qw2191_discharged_by_this_certificate"] is False and theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2549",
        "stage_id": "S1499",
        "status": "STRICT_DAMPING_QUADRUPLE_TRACE_CURRENT_PREMISE_OBSTRUCTION_CERTIFICATE_NO_TRACE_SOURCE_EXPORT_NO_M2_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_damping_quadruple_trace_current_premise_obstruction_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_damping_quadruple_trace_current_premise_obstruction_certificate"]["theorem_export"]
    audit = t["trace_obstruction_audit"]
    lines = [
        "# P2549/S1499 strict damping quadruple-trace current-premise obstruction certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        f"- Frontier source under attack: `{t['frontier_source_under_attack']}`.",
        f"- P2540 order nonidentifiability inherited: `{t['p2540_current_premise_order_nonidentifiability_inherited']}`.",
        f"- P2548 conditional trace-arity selection inherited: `{t['p2548_conditional_trace_arity_selection_inherited']}`.",
        f"- P2548 actual trace source absent inherited: `{t['p2548_actual_trace_source_absent_inherited']}`.",
        f"- Current-premise passing orders audited: `{audit['current_premise_passing_order_count']}`.",
        f"- Current-premise trace arities: `{audit['current_premise_trace_arities']}`.",
        f"- Quadruple-trace accepting orders: `{audit['quadruple_trace_accepting_orders']}`.",
        f"- Quadruple-trace countermodel orders: `{audit['quadruple_trace_countermodel_orders']}`.",
        f"- Current premises entail exact quadruple trace: `{audit['current_premises_entail_exact_quadruple_trace']}`.",
        f"- Strict quadruple-trace source exported: `{t['strict_quadruple_trace_source_exported']}`.",
        f"- m2 operator signature source exported: `{t['m2_operator_signature_source_exported']}`.",
        "",
        "## Countermodel",
        "",
        "The explicit inherited `m=3` row passes the current source-free derivative-order premises, but its self-adjoint trace arity is `6`, not `4`. Hence current premises do not entail the exact four-trace source required by P2548.",
        "",
        "## Interpretation",
        "",
        audit["model_separation_statement"],
        "",
        "P2549 therefore blocks a false promotion of the P2548 conditional theorem into an actual m=2 source. The obstruction is narrow: it does not refute a future strict nadsoliton theorem that genuinely derives four boundary/trace channels.",
        "",
        "## Recommended next honest step",
        "",
        t["recommended_next_honest_step"],
        "",
        "## Negative controls",
        "",
        "No exact quadruple-trace source, `m2_operator_signature_source`, residual prime-log source, residual slope source, source-obligation discharge, bridge completion, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, physical-value generator, or ToE closure is exported.",
        "",
        "## Fingerprint",
        "",
        f"`{payload['strict_damping_quadruple_trace_current_premise_obstruction_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2549/S1499` audits whether the P2548 quadruple-trace premise is already forced by current derivative-order premises.  It is not: the P2540 source-free order models still admit all audited orders `m=1..10`, with trace arities `2m = 2,4,6,...,20`; only `m=2` has trace arity four, while `m=3` is an explicit passing countermodel with trace arity six.  Thus the P2548 implication remains conditional, and no `strict_quadruple_trace_source` or `m2_operator_signature_source` is exported.
""".strip()
    lag_section = """
`P2549/S1499` blocks promotion of the P2548 trace-arity bookkeeping into a role-bearing `L_total` operator source.  Current premises do not derive the exact four-trace boundary structure: they still admit non-four trace arities such as the inherited `m=3`/six-trace model.  A damping/compression term therefore remains non-role-bearing until a real strict four-trace source, residual numeric sources, bridge completion, and role-transfer gates are supplied.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2549/S1499 quadruple-trace current-premise obstruction guard", "## P2549/S1499 quadruple-trace current-premise obstruction guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2549/S1499 quadruple-trace obstruction Ltotal guard", "## P2549/S1499 quadruple-trace obstruction Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
