#!/usr/bin/env python3
"""Scratch probe: strict-completion closure-plan dependency certificate.

This probe is deliberately not another local GF(2) gadget.  It converts the
current repo frontier statements into a finite, auditable dependency plan for
closing the remaining bridge/theory gap without violating guardrails.

It reads the current strict-completion chain report plus selected frontier
packets and emits a dependency matrix, topological work order, and hard blockers.
The certificate is a planning/audit object only: it does not close the bridge,
does not discharge QW-2191, and does not claim ToE closure.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FAR = ROOT / "fundamental_action_reconstruction"
OUT_JSON = HERE / "bridge_strict_completion_closure_plan_dependency_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_closure_plan_dependency_certificate_report.md"
CHAIN_REPORT = HERE / "bridge_strict_completion_certificate_chain_integrity_report.json"
SOURCE_TEXTS = {
    "K1_kernel_split_guardrail": FAR / "K1_LEGACY_ONTOLOGICAL_KERNEL_VS_STRICT_GATE_KERNEL_SPLIT_NOTE.md",
    "K2_strict_gate_chain": FAR / "K2_STRICT_GATE_KERNEL_DERIVATION_CHAIN_NOTE.md",
    "F3_frontier_classification": FAR / "F3_CURRENT_FAR_FRONTIER_KERNEL_ARTIFACT_SENSITIVITY_CLASSIFICATION_PACKET.md",
    "S2_strategic_reorientation": FAR / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
    "P1445_next_honest_step": FAR / "P1445_STRICT_NEXT_HONEST_STEP_F_NADSOLITON_TO_LSM_PLUS_LGR_NO_LEGACY_BRIDGE_PACKET_PL.md",
    "N679_selector_boundary": FAR / "N679_CURRENT_STRICT_T172_STRICT_CORE_SELECTOR_CLOSURE_FRONTIER_BOUNDARY_THEOREM.md",
}

PLAN_STEPS = [
    "finite_strict_completion_ledger",
    "legacy_to_strict_completion_bridge_guardrail",
    "strict_local_selector_margin_monotonicity_witness",
    "selector_compatibility_certificate_scope_limited",
    "strict_orientation_or_global_qw2191_source",
    "strict_F_nadsoliton_to_LSM_plus_LGR_bridge",
    "strict_core_ToE_closure",
]

DEPENDENCIES = {
    "finite_strict_completion_ledger": [],
    "legacy_to_strict_completion_bridge_guardrail": ["finite_strict_completion_ledger"],
    "strict_local_selector_margin_monotonicity_witness": ["legacy_to_strict_completion_bridge_guardrail"],
    "selector_compatibility_certificate_scope_limited": ["strict_local_selector_margin_monotonicity_witness"],
    "strict_orientation_or_global_qw2191_source": ["selector_compatibility_certificate_scope_limited"],
    "strict_F_nadsoliton_to_LSM_plus_LGR_bridge": ["strict_orientation_or_global_qw2191_source"],
    "strict_core_ToE_closure": ["strict_F_nadsoliton_to_LSM_plus_LGR_bridge"],
}


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    if not path.exists():
        raise FileNotFoundError(f"missing source text: {path}")
    return path.read_text(encoding="utf-8")


def contains_all(text: str, needles: list[str]) -> bool:
    return all(needle in text for needle in needles)


def dependency_matrix() -> list[dict[str, Any]]:
    rows = []
    for dependent in PLAN_STEPS:
        rows.append({
            "step": dependent,
            "depends_on": DEPENDENCIES[dependent],
            "dependency_bits_in_plan_order": [1 if prerequisite in DEPENDENCIES[dependent] else 0 for prerequisite in PLAN_STEPS],
        })
    return rows


def transitive_prerequisites(step: str) -> set[str]:
    result: set[str] = set()
    stack = list(DEPENDENCIES[step])
    while stack:
        current = stack.pop()
        if current in result:
            continue
        result.add(current)
        stack.extend(DEPENDENCIES[current])
    return result


def build_payload() -> dict[str, Any]:
    chain = load_json(CHAIN_REPORT)
    texts = {name: load_text(path) for name, path in SOURCE_TEXTS.items()}

    source_keyword_checks = {
        "S2_requires_legacy_strict_completion_bridge_priority": contains_all(texts["S2_strategic_reorientation"], ["legacy -> strict` completion bridge", "role-transfer audit after bridge completion", "QW-2191"]),
        "P1445_recommends_S1_local_selector_margin": contains_all(texts["P1445_next_honest_step"], ["strict-local selector margin monotonicity witness", "bez claimu global closure"]),
        "N679_keeps_qw2191_and_orientation_open": contains_all(texts["N679_selector_boundary"], ["kernel-alone/global QW-2191 discharge", "directed/sign-sensitive physical orientation datum"]),
        "K1_blocks_legacy_strict_identity": contains_all(texts["K1_kernel_split_guardrail"], ["K_legacy_ont", "K_strict_gate", "not yet rigorously identified"]),
        "K2_blocks_direct_legacy_derivation_reading": contains_all(texts["K2_strict_gate_chain"], ["not directly derived", "bridge theorem"]),
        "F3_prefers_kernel_split_robust_frontier": contains_all(texts["F3_frontier_classification"], ["kernel-split-robust", "QW-2191"]),
    }

    base_ledger_certified = (
        chain["chain_summary"]["exact_APD_completion_certified"]
        and chain["chain_summary"]["phase_gf2_linear_system_certified"]
        and chain["chain_summary"]["positive_factor_sign_separation_certified"]
        and not chain["chain_summary"]["strict_dynamic_derivation_exported"]
        and not chain["chain_summary"]["bridge_theorem_exported"]
    )

    step_status = {
        "finite_strict_completion_ledger": "certified" if base_ledger_certified else "blocked",
        "legacy_to_strict_completion_bridge_guardrail": "next_recommended_guardrail_updated",
        "strict_local_selector_margin_monotonicity_witness": "blocked_until_bridge_guardrail_accepted",
        "selector_compatibility_certificate_scope_limited": "blocked_until_S1_PASS",
        "strict_orientation_or_global_qw2191_source": "open_hard_blocker",
        "strict_F_nadsoliton_to_LSM_plus_LGR_bridge": "blocked_by_selector_orientation_source",
        "strict_core_ToE_closure": "not_claimed_blocked_by_bridge_and_selector",
    }
    closure_order = [step for step in PLAN_STEPS if step_status[step] != "certified"]
    transitive_rows = [
        {
            "step": step,
            "transitive_prerequisites": sorted(transitive_prerequisites(step), key=PLAN_STEPS.index),
            "status": step_status[step],
        }
        for step in PLAN_STEPS
    ]

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_CLOSURE_PLAN_DEPENDENCY_CERTIFICATE__REPO_FRONTIER_MATRIX_NO_CLOSURE_CLAIM",
        "status": "strict-completion-closure-plan-dependency-matrix-exported-no-false-pass",
        "source_reports": {
            "strict_completion_chain_integrity_report": str(CHAIN_REPORT.relative_to(ROOT)),
            **{name: str(path.relative_to(ROOT)) for name, path in SOURCE_TEXTS.items()},
        },
        "grep_disambiguation": {
            "searched_terms": [
                "closure plan",
                "domknięcie mostu",
                "bridge plan",
                "strict-local selector margin monotonicity witness",
                "coordinate syndrome generator basis",
                "kernel-alone/global QW-2191 discharge",
                "F_nadsoliton => L_SM + L_GR",
                "strict-core ToE closure",
            ],
            "finding": "Repo grep finds many blocker/frontier packets and the existing coordinate-syndrome certificates, but no single strict-completion closure-plan dependency matrix tying the current ledger to S1 selector-margin, QW-2191/orientation, F_nadsoliton=>L_SM+L_GR, and ToE closure before this report.",
        },
        "source_keyword_checks": source_keyword_checks,
        "dependency_certificate": {
            "plan_steps_in_topological_order": PLAN_STEPS,
            "dependency_matrix_rows": dependency_matrix(),
            "transitive_prerequisite_rows": transitive_rows,
            "next_open_steps_in_order": closure_order,
            "recommended_next_step": "legacy_to_strict_completion_bridge_guardrail",
            "recommended_next_step_source": "S2_strategic_reorientation",
            "next_selector_step_after_bridge_guardrail": "strict_local_selector_margin_monotonicity_witness",
        },
        "closure_plan_summary": {
            "chain_ledger_currently_cross_consistent": base_ledger_certified,
            "all_source_keyword_checks_pass": all(source_keyword_checks.values()),
            "dependency_matrix_is_triangular_in_plan_order": all(
                all(PLAN_STEPS.index(prereq) < PLAN_STEPS.index(row["step"]) for prereq in row["depends_on"])
                for row in dependency_matrix()
            ),
            "recommended_next_step_is_legacy_strict_bridge_guardrail": closure_order[0] == "legacy_to_strict_completion_bridge_guardrail",
            "S1_selector_margin_remains_next_selector_subproblem": "strict_local_selector_margin_monotonicity_witness" in closure_order,
            "qw2191_or_orientation_remains_hard_blocker": step_status["strict_orientation_or_global_qw2191_source"] == "open_hard_blocker",
            "toe_closure_not_claimed": step_status["strict_core_ToE_closure"] != "certified",
            "legacy_strict_identity_not_used": source_keyword_checks["K1_blocks_legacy_strict_identity"] and source_keyword_checks["K2_blocks_direct_legacy_derivation_reading"],
            "role_transfer_audit_required": source_keyword_checks["S2_requires_legacy_strict_completion_bridge_priority"],
        },
        "closure_plan_for_bridge_and_theory": [
            {
                "order": 1,
                "step": "Accept the restored legacy -> strict completion-bridge guardrail and keep raw identity blocked.",
                "pass_condition": "K_legacy_ont is treated as an intermediate incomplete bridge kernel; K_strict_gate is the completed/enriched continuation only through explicit completion evidence.",
                "failure_action": "do not proceed to role transfer or selector closure; repair the bridge guardrail first.",
            },
            {
                "order": 2,
                "step": "Finish local strict-only selector-margin witness S1 on a kernel-split-robust provider class.",
                "pass_condition": "monotonic selector-margin growth on a fixed perturbation grid plus replay stability; scope remains local only.",
                "failure_action": "export obstruction packet and choose a new provider class or blocker cut.",
            },
            {
                "order": 3,
                "step": "If S1 passes, build a scope-limited selector compatibility certificate.",
                "pass_condition": "explicitly state selector premise/scope and show compatibility with current finite completion ledger.",
                "failure_action": "do not claim strict-core selector closure; return to provider-class search.",
            },
            {
                "order": 4,
                "step": "Attack the remaining QW-2191/orientation source boundary.",
                "pass_condition": "export either an internal strict selector/orientation source or a theorem showing exactly what extra premise is unavoidable.",
                "failure_action": "mark bridge/theory closure blocked by selector/orientation source.",
            },
            {
                "order": 5,
                "step": "Only after selector/orientation source exists, attempt strict F_nadsoliton => L_SM + L_GR bridge.",
                "pass_condition": "derive the field/action target without legacy role transfer and with all pass-scope premises explicit.",
                "failure_action": "export a non-bridge theorem or narrowed missing-object theorem.",
            },
            {
                "order": 6,
                "step": "Only after the strict bridge exists, evaluate strict-core ToE closure.",
                "pass_condition": "all needed strict-side sources, selector, orientation, SM and GR action pieces are exported theorem-level objects.",
                "failure_action": "state the residual blocker explicitly; no ToE closure claim.",
            },
        ],
        "blocker_context": {
            "resolved_locally": [
                "The current finite strict-completion ledger is cross-consistent.",
                "The restored legacy->strict bridge guardrail is represented in the finite matrix before S1 selector work.",
                "The plan explicitly separates local S1 progress from global closure claims.",
            ],
            "still_open": [
                "legacy->strict completion bridge still lacks a full map and role-transfer audit",
                "strict-local selector margin monotonicity witness S1 not exported here",
                "kernel-alone/global QW-2191 discharge remains open",
                "directed/sign-sensitive physical orientation datum remains open",
                "strict F_nadsoliton => L_SM + L_GR bridge remains open",
                "strict ToE closure remains open",
            ],
        },
        "proof_certificate": {
            "grep_step": "Repo grep disambiguates this from prior component-quotient and selector-boundary reports: no prior single dependency matrix for bridge/theory closure was exported.",
            "ledger_step": "The chain-integrity report is loaded and its non-closure finite strict-completion ledger flags must be certified before the closure plan treats the ledger as input.",
            "matrix_step": "The dependency matrix is triangular in plan order: each non-initial step depends only on earlier steps.",
            "next_step": "The first uncertified step is the restored legacy->strict completion-bridge guardrail from S2; S1 remains the next selector subproblem after that guardrail is accepted.",
            "blocker_step": "N679 keeps kernel-alone/global QW-2191 discharge and directed/sign-sensitive orientation open; therefore bridge and ToE closure remain blocked.",
            "theoretical_limit": "This report is a plan/dependency certificate; it does not complete the legacy->strict bridge, prove S1, discharge QW-2191, build F_nadsoliton=>L_SM+L_GR, or close ToE.",
        },
        "hard_limits": [
            "K_strict_gate remains the current live/full operational kernel.",
            "No unqualified identity K_legacy_ont == K_strict_gate is claimed or used.",
            "No legacy physical-role transfer to K_strict_gate is claimed.",
            "No strict-local selector margin witness is claimed in this report.",
            "No QW-2191 selector discharge is claimed.",
            "No strict F_nadsoliton => L_SM + L_GR bridge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    summary = payload["closure_plan_summary"]
    lines = [
        "# Strict-completion closure-plan dependency certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This report exports a finite dependency matrix and closure plan extracted",
        "from the current repo frontier. It is not a closure theorem.",
        "",
        "## Summary",
        "",
    ]
    for key, value in summary.items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Plan", ""])
    for row in payload["closure_plan_for_bridge_and_theory"]:
        lines.append(f"{row['order']}. {row['step']} Pass: {row['pass_condition']}")
    lines.extend(["", "## Hard limits", ""])
    for limit in payload["hard_limits"]:
        lines.append(f"- {limit}")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(payload)
    print(json.dumps(payload["closure_plan_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
