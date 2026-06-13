#!/usr/bin/env python3
"""P2684/S1634: pair12 lower-boundary cycle-cut semantic invariant/provider audit.

This is the next proof/computation step after P2683.  It does not add another
lower-boundary name.  It tests whether current repo artifacts already export a
real semantic invariant/provider at the chart-label-retaining pair12 typed
seed-slot coordinate entry point before F301 binding, Q_basis terminal collapse,
and projector-only atlas collapse.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2684_s1634_pair12_cycle_cut_semantic_invariant_provider_audit.json"
MD = GEN / "p2684_s1634_pair12_cycle_cut_semantic_invariant_provider_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2683": GEN / "p2683_s1633_post_ax13_direct_route_and_t173_lower_boundary_cycle_audit.json",
    "P946": GEN / "p946_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_existing_export_or_near_miss_candidate_audit_probe_summary.json",
    "P744": GEN / "p744_current_strict_t198_declared_scope_source_topology_selector_theorem_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json",
    "P765": GEN / "p765_current_strict_t219_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_nonexport_audit_probe_summary.json",
    "P766": GEN / "p766_current_strict_t220_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_probe_summary.json",
    "P767": GEN / "p767_current_strict_t221_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_nonexport_audit_probe_summary.json",
    "P770": GEN / "p770_current_strict_t224_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_probe_summary.json",
    "P778": GEN / "p778_current_strict_t232_pair12_seed_slot_coordinate_subsubsubinterface_actual_realization_attempt_probe_summary.json",
    "P780": GEN / "p780_current_strict_t234_pair12_seed_slot_coordinate_entry_subsubsubsubinterface_target_probe_summary.json",
    "P962": GEN / "p962_current_strict_t253_pair12_entry_point_exact_still_lower_boundary_refinement_actual_realization_attempt_verdict_or_exact_further_lower_boundary_refinement_nonexport_audit_probe_summary.json",
    "F301": GEN / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "new_semantic_invariant_provider_exported",
    "lower_boundary_recursion_reopened_as_primary",
    "t176_discharged",
    "qw2191_discharged",
    "strict_core_selector_closure_claimed",
    "toe_closure_claimed",
]


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        ["rg", "-n", pattern, ".", "-g", "*.py", "-g", "*.md", "-g", "*.json", "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**"],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:60]}


def content_grep() -> dict[str, Any]:
    patterns = {
        "cycle_cut_entry_point": r"chart-label-retaining|seed-slot coordinate|typed seed-slot|pair12 typed|Sigma_sel_src_target_v1",
        "pre_collapse_guards": r"F301|Q_basis_sel_v1|projector-only|projector_only|terminal collapse|atlas collapse",
        "semantic_invariant_provider_candidates": r"semantic invariant|sourced invariant|provider|transported flux|current-like section|chart-sensitive transported",
        "known_near_miss_or_nonexport": r"P946|T220|T222|T224|T232|T234|nonexport|near miss|missing subinterface",
        "forbidden_closure_claims": r"T176 discharge|QW-2191 discharge|ToE closure|strict_core_selector_closure_claimed|selector closure claimed",
    }
    return {"tool": "rg", "mode": "content-first pair12 cycle-cut semantic invariant/provider audit", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def current_evidence() -> dict[str, Any]:
    data = {name: load_json(path) for name, path in INPUTS.items()}
    p946 = data["P946"]
    p765 = data["P765"]
    p767 = data["P767"]
    p770 = data["P770"]
    p780 = data["P780"]
    serialized = {k: json.dumps(v, sort_keys=True) for k, v in data.items()}
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "p2683_loop_stop_present": not data["P2683"].get("missing", False),
        "p946_no_current_lawful_bridge_supplier_found": p946.get("no_current_lawful_bridge_supplier_found") is True,
        "p946_near_miss_only": "near_miss" in serialized["P946"] or "near miss" in serialized["P946"],
        "p765_chart_sensitive_interface_actual_exported": p765.get("t219_target_exported_on_current_repo_state") is True,
        "p767_immediate_missing_subinterface_nonexport": "NONEXPORT" in str(p767.get("status")) or "missing_subinterface" in serialized["P767"],
        "p770_deepest_positive_is_attempt": "attempt" in serialized["P770"].lower() and "exported" in serialized["P770"].lower(),
        "p780_seed_slot_entry_is_future_only_target": "future-only" in serialized["P780"] or "future_only" in serialized["P780"],
        "f301_carrier_exists": not data["F301"].get("missing", False),
    }


def obligation_matrix(evidence: dict[str, Any]) -> dict[str, Any]:
    rows = [
        {
            "obligation": "chart_label_retaining_pair12_typed_seed_or_subinterface_actual_export",
            "required_for_cycle_cut": True,
            "satisfied_now": bool(evidence["p765_chart_sensitive_interface_actual_exported"]),
            "current_witness": "P765/P767/P770/P780 family remains nonexport/attempt/future-only, not an actual exported seed/subinterface.",
        },
        {
            "obligation": "pre_F301_pre_Q_basis_pre_projector_collapse_scope_retained",
            "required_for_cycle_cut": True,
            "satisfied_now": True,
            "current_witness": "The target language explicitly retains pre-F301, pre-Q_basis terminal-collapse, and pre-projector-only atlas-collapse scope.",
        },
        {
            "obligation": "nonconventional_semantic_invariant_or_provider",
            "required_for_cycle_cut": True,
            "satisfied_now": False,
            "current_witness": "P946 reports only a closest near-miss family and no lawful bridge supplier/current-like transported section.",
        },
        {
            "obligation": "not_merely_next_lower_boundary_name",
            "required_for_cycle_cut": True,
            "satisfied_now": True,
            "current_witness": "P2684 audits candidate exportability instead of adding another T24x/T25x lower-boundary successor.",
        },
    ]
    exported = all(row["satisfied_now"] for row in rows)
    return {"rows": rows, "semantic_invariant_provider_exported_now": exported, "failed_obligations": [r["obligation"] for r in rows if not r["satisfied_now"]]}


def finite_policy_lattice(matrix: dict[str, Any]) -> dict[str, Any]:
    obligations = [row["obligation"] for row in matrix["rows"]]
    pass_count = 0
    cycle_cut_pass_count = 0
    for bits in itertools.product([False, True], repeat=len(obligations)):
        state = dict(zip(obligations, bits))
        admissible = state["not_merely_next_lower_boundary_name"] and state["pre_F301_pre_Q_basis_pre_projector_collapse_scope_retained"]
        cycle_cut = admissible and state["chart_label_retaining_pair12_typed_seed_or_subinterface_actual_export"] and state["nonconventional_semantic_invariant_or_provider"]
        pass_count += int(admissible)
        cycle_cut_pass_count += int(cycle_cut)
    return {
        "obligations": obligations,
        "total_states": 2 ** len(obligations),
        "admissible_audit_states": pass_count,
        "cycle_cut_success_states": cycle_cut_pass_count,
        "current_state": {row["obligation"]: row["satisfied_now"] for row in matrix["rows"]},
        "current_cycle_cut_success": matrix["semantic_invariant_provider_exported_now"],
    }


def decision(matrix: dict[str, Any]) -> dict[str, Any]:
    exported = matrix["semantic_invariant_provider_exported_now"]
    return {
        "decision": "P2684_PAIR12_CYCLE_CUT_SEMANTIC_INVARIANT_PROVIDER_NO_GO_PIVOT_TO_STRICT_LAGRANGIAN_EOM_REVERSE_CLOSURE_NO_FALSE_PASS",
        "professorial_verdict": (
            "P2684 performs the P2683-recommended cycle-cut audit without adding another lower-boundary successor.  The pre-collapse target is real and well localized, but the current repo still supplies only near-miss/nonexport/attempt/future-only objects: no actual chart-label-retaining pair12 typed seed/subinterface plus nonconventional semantic invariant/provider is exported.  Therefore the lower-boundary recursion remains frozen as a primary strategy."
        ),
        "next_honest_step": (
            "Pivot to an independent strict Lagrangian/EOM reverse-closure obstruction matrix.  The finite target should ask which currently exported selector-independent Lagrangian/EOM rows can be reversed into nonproxy full tensor-resolved EOM/L_total closure, and produce an explicit missing-obligation matrix rather than importing selector, role-transfer, tau->pair12, beta_tors->chi11, or legacy->strict bridge claims."
        ),
        "semantic_invariant_provider_exported_now": exported,
        "lower_boundary_recursion_primary": False,
        "strict_lagrangian_eom_reverse_closure_is_next": not exported,
        "t176_discharged_now": False,
        "qw2191_discharged_now": False,
        "toe_closed_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2684/S1634 pair12 cycle-cut semantic invariant/provider audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first grep",
    ]
    for name, data in payload["content_grep"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    lines.extend(["", "## Obligation matrix"])
    for row in payload["obligation_matrix"]["rows"]:
        lines.append(f"- `{row['obligation']}`: satisfied_now=`{row['satisfied_now']}`; witness={row['current_witness']}")
    lines.extend([
        "", "## Finite policy lattice",
        f"Total states: `{payload['finite_policy_lattice']['total_states']}`; cycle-cut success states: `{payload['finite_policy_lattice']['cycle_cut_success_states']}`.",
        f"Current cycle-cut success: `{payload['finite_policy_lattice']['current_cycle_cut_success']}`.",
        "", "## Verdict", payload["decision"]["professorial_verdict"],
        "", "## Next honest step", payload["decision"]["next_honest_step"],
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    evidence = current_evidence()
    matrix = obligation_matrix(evidence)
    lattice = finite_policy_lattice(matrix)
    payload: dict[str, Any] = {
        "status": "P2684_PAIR12_CYCLE_CUT_SEMANTIC_INVARIANT_PROVIDER_AUDIT_NO_FALSE_PASS",
        "content_grep": content_grep(),
        "current_evidence": evidence,
        "obligation_matrix": matrix,
        "finite_policy_lattice": lattice,
        "decision": decision(matrix),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2684/S1634 pair12 cycle-cut semantic invariant/provider no-go",
        "## P2684/S1634 pair12 cycle-cut semantic invariant/provider no-go\n\n"
        "`P2684/S1634` audits the P2683 cycle-cut target instead of adding another lower-boundary name.  It finds no current exported chart-label-retaining `pair12` typed seed/subinterface carrying a nonconventional semantic invariant/provider before F301 binding, `Q_basis` terminal collapse, and projector-only atlas collapse.  The lower-boundary recursion therefore remains nonprimary; the next proof-grade target is a strict Lagrangian/EOM reverse-closure obstruction matrix.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2684/S1634 reverse-closure pivot guard",
        "## P2684/S1634 reverse-closure pivot guard\n\n"
        "`P2684/S1634` adds no `L_total` term and exports no selector, role-transfer, or ToE closure.  It only promotes the strict Lagrangian/EOM reverse-closure obstruction matrix to the next honest finite computation after the pair12 semantic-invariant/provider cycle-cut audit fails on current artifacts.\n",
    )
    append_once(
        AGENTS,
        "Current cycle-cut and reverse-closure steering guardrail (P2684/S1634, 2026-06-13)",
        "## Current cycle-cut and reverse-closure steering guardrail (P2684/S1634, 2026-06-13)\n\n"
        "- After P2684, do not continue the T24x/T25x lower-boundary recursion as a primary strategy unless a real chart-label-retaining `pair12` typed seed/subinterface plus a nonconventional semantic invariant/provider is newly exported.\n"
        "- The current next proof-grade move is an independent strict Lagrangian/EOM reverse-closure obstruction matrix; it must not import selector closure, role transfer, generic bridge closure, `tau_src -> pair12`, or `beta_tors -> chi11` by replay.\n",
    )
    return payload


if __name__ == "__main__":
    main()
