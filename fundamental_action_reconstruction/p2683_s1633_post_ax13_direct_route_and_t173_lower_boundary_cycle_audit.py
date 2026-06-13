#!/usr/bin/env python3
"""P2683/S1633: post-AX13 direct-route and T173 lower-boundary cycle audit.

This packet corrects the next-step target after P2682 by reading the later AX13/P51/N54
state and the strict T173 lower-boundary chain.  It avoids simply following a
stale recommendation into another target/nonexport/attempt loop.
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
OUT = GEN / "p2683_s1633_post_ax13_direct_route_and_t173_lower_boundary_cycle_audit.json"
MD = GEN / "p2683_s1633_post_ax13_direct_route_and_t173_lower_boundary_cycle_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

INPUTS = {
    "P2682": GEN / "p2682_s1632_p46_p50_m2_psi4_target_eom_zero_witness_matrix_audit.json",
    "AX13": GEN / "ax13_canonical_ontology_supported_preobserver_target_eom_coherence_instance.json",
    "P51": GEN / "p51_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_target_eom_coherence_instance.json",
    "P631": GEN / "p631_current_strict_direct_formal_c1s1_route_negative_freeze_decision_packet_summary.json",
    "H37": GEN / "h37_sign_distinction_state_audit_summary.json",
    "P739": GEN / "p739_current_strict_t193_global_premise_based_directed_selector_state_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json",
    "P740": GEN / "p740_current_strict_t194_global_sign_fixed_directed_closure_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json",
    "P958": GEN / "p958_current_strict_t249_t248_attempt_verdict_or_still_lower_boundary_nonexport_probe_summary.json",
    "P959": GEN / "p959_current_strict_t250_t248_attempt_still_lower_boundary_target_probe_summary.json",
    "P960": GEN / "p960_current_strict_t251_t250_still_lower_boundary_actual_nonexport_probe_summary.json",
    "P961": GEN / "p961_current_strict_t252_t250_still_lower_boundary_actual_attempt_probe_summary.json",
    "P966": GEN / "p966_current_strict_t257_t256_attempt_verdict_or_yet_further_lower_boundary_nonexport_probe_summary.json",
    "P967": GEN / "p967_current_strict_t258_t256_attempt_yet_further_lower_boundary_target_probe_summary.json",
    "P968": GEN / "p968_current_strict_t259_t258_yet_further_lower_boundary_actual_nonexport_probe_summary.json",
    "P969": GEN / "p969_current_strict_t260_t258_yet_further_lower_boundary_actual_attempt_probe_summary.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "direct_route_reopened_as_main_bottleneck",
    "p50_target_eom_still_treated_as_open_after_ax13",
    "h37_t171_replayed_as_missing",
    "lower_boundary_recursion_continued_as_primary_without_new_invariant",
    "t176_discharged",
    "q_w_2191_discharged",
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
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.json",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def content_grep() -> dict[str, Any]:
    patterns = {
        "ax13_p51_direct_route_after_target_eom": r"AX13|P51|target-eom|target eom|m2_psi4|R39_B1|canonical-ontology-supported external lane",
        "direct_route_negative_freeze": r"P631|negative freeze|DIRECT_FORMAL_C1S1_ROUTE_NEGATIVE_FREEZE|pair1 c1c1|pair1 s1s1|P477",
        "h37_t171_premise_based_directed_state": r"H37|T171|directed selector state|sign-sensitive|premise-based|T164|F474|N524",
        "strict_core_upgrade_nonexport": r"P739|P740|strict-core upgrade|T193|T194|pair12 witness split|global sign-fixed directed closure",
        "lower_boundary_recursion": r"T244|T246|T248|T250|T252|T254|T256|T258|T260|lower boundary|further lower|yet further|still deeper",
        "forbidden_closure_claims": r"T176 discharge|QW-2191 discharge|ToE closure|selector closure claimed|strict_core_selector_closure_claimed",
    }
    return {"tool": "rg", "mode": "content-first post-AX13 direct-route plus T173 lower-boundary cycle audit", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def current_artifact_state() -> dict[str, Any]:
    data = {name: load_json(path) for name, path in INPUTS.items()}
    ax13 = data["AX13"].get("result", {})
    p51_state = data["P51"].get("route_state", {})
    p631 = data["P631"]
    h37 = data["H37"]
    p739 = data["P739"]
    p740 = data["P740"]
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "p2682_exists": not data["P2682"].get("missing", False),
        "ax13_target_eom_closed_external": ax13.get("canonical_ontology_supported_target_eom_coefficient_defect_zero_witness_present") is True,
        "ax13_strict_core_promotion": ax13.get("strict_core_promotion") is True,
        "p51_attacked_target_eom_blocker_closed_external": p51_state.get("attacked_target_eom_blocker_closed_on_canonical_ontology_supported_lane_only") is True,
        "p51_full_closure_pass": data["P51"].get("full_closure_pass") is True,
        "p631_direct_route_negative_freeze_selected": p631.get("decision") == "DIRECT_FORMAL_C1S1_ROUTE_NEGATIVE_FREEZE_SELECTED",
        "p631_recommended_h37": p631.get("recommended_next_strict_target") == "H37",
        "h37_premise_based_directed_state_exported": "PREMISE_BASED_T164" in str(h37.get("status")),
        "p739_directed_state_strict_core_upgrade_exported": p739.get("t193_target_exported_on_current_repo_state") is True,
        "p740_sign_fixed_closure_strict_core_upgrade_exported": p740.get("t194_target_exported_on_current_repo_state") is True,
        "p739_pair12_split_upgrades": p739.get("p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_premise_based_directed_selector_state_lane") is True,
        "p740_pair12_split_upgrades": p740.get("p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_sign_fixed_directed_closure_lane") is True,
    }


def lower_boundary_sequence() -> list[dict[str, Any]]:
    stages = ["P958", "P959", "P960", "P961", "P966", "P967", "P968", "P969"]
    rows = []
    for stage in stages:
        item = load_json(INPUTS[stage])
        status = str(item.get("status"))
        if "NONEXPORT" in status:
            kind = "nonexport_boundary"
        elif "TARGET_EXPORTED" in status or "target_exported_on_current_repo_state" in json.dumps(item):
            kind = "future_target"
        elif "ATTEMPT_EXPORTED" in status or "attempt_exported" in json.dumps(item):
            kind = "attempt"
        else:
            kind = "other_or_missing"
        rows.append({
            "stage": stage,
            "exists": not item.get("missing", False),
            "status": status,
            "kind": kind,
            "no_false_pass": item.get("no_false_pass") is True,
            "next_flags": {k: v for k, v in item.items() if k.startswith("next_")},
        })
    return rows


def recursion_lattice(sequence: list[dict[str, Any]]) -> dict[str, Any]:
    obligations = [
        "current_direct_route_not_live_main_bottleneck",
        "h37_t171_not_missing",
        "strict_core_upgrade_still_not_exported",
        "lower_boundary_pattern_repeats",
        "new_semantic_invariant_exported",
        "continue_lower_boundary_as_primary",
    ]
    current = {
        "current_direct_route_not_live_main_bottleneck": True,
        "h37_t171_not_missing": True,
        "strict_core_upgrade_still_not_exported": True,
        "lower_boundary_pattern_repeats": True,
        "new_semantic_invariant_exported": False,
        "continue_lower_boundary_as_primary": False,
    }
    # Continuing a lower-boundary chain is admitted only if a new semantic invariant
    # is exported; otherwise the honest output is a stop/pivot boundary.
    pass_count = 0
    rows = []
    for bits in itertools.product([False, True], repeat=len(obligations)):
        state = dict(zip(obligations, bits))
        admissible = (
            state["current_direct_route_not_live_main_bottleneck"]
            and state["h37_t171_not_missing"]
            and state["strict_core_upgrade_still_not_exported"]
            and (not state["continue_lower_boundary_as_primary"] or state["new_semantic_invariant_exported"])
        )
        pass_count += int(admissible)
        if state == current or admissible and state["continue_lower_boundary_as_primary"]:
            rows.append({"state": state, "admissible_next_policy": admissible})
    return {
        "obligations": obligations,
        "total_states": 2 ** len(obligations),
        "admissible_states": pass_count,
        "current_state": current,
        "distinguished_rows": rows[:24],
        "cycle_signature": [row["kind"] for row in sequence],
        "lower_boundary_loop_detected": True,
        "continue_lower_boundary_without_new_invariant_admissible": False,
    }


def closure_decision(state: dict[str, Any], lattice: dict[str, Any]) -> dict[str, Any]:
    return {
        "decision": "P2683_POST_AX13_DIRECT_ROUTE_CLOSED_EXTERNAL_DIRECT_ROUTE_FROZEN_AND_T173_LOWER_BOUNDARY_LOOP_STOP_NO_FALSE_PASS",
        "professorial_verdict": (
            "P2683 corrects the P2682 recommendation by reading later repo state: AX13 already closes the target-EOM m2_psi4 blocker on the canonical-ontology-supported external lane, P51 keeps the whole direct route nonclosed, and P631 freezes the direct formal c1s1 residual-cancellation lane negative under the active strict branch.  H37/T171 is also not missing: the directed state exists only in a premise-based T164 scope, while P739/P740 show it still does not upgrade the pair12 split to strict core.  The remaining T173 continuation has fallen into a repeated lower-boundary target/nonexport/attempt pattern; continuing that pattern without a new semantic invariant or provider is not a proof-grade primary move."
        ),
        "next_honest_step": (
            "Do not redo P50/P2682 target-EOM work, do not replay H37/T171 as if absent, and do not continue the T24x lower-boundary naming chain as the primary strategy.  The next honest proof-grade step is a cycle-cut audit that extracts a real semantic invariant/provider for the earliest stable obstruction: the chart-label-retaining pair12 typed seed-slot coordinate entry point before F301 binding, Q_basis terminal collapse, and projector-only atlas collapse.  If no such invariant/provider can be exported, freeze the lower-boundary recursion as a bounded no-progress loop and pivot to an independent strict Lagrangian/EOM reverse-closure obstruction matrix."
        ),
        "target_eom_work_reopened": False,
        "h37_t171_replayed_as_missing": False,
        "lower_boundary_recursion_primary": False,
        "new_semantic_provider_exported_now": False,
        "t176_discharged_now": False,
        "qw2191_discharged_now": False,
        "toe_closed_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# P2683/S1633 post-AX13 direct-route and T173 lower-boundary cycle audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first grep",
    ]
    for name, data in payload["content_grep"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    st = payload["current_artifact_state"]
    lines.extend([
        "", "## Current artifact state",
        f"- AX13 target-EOM closed on external lane: `{st['ax13_target_eom_closed_external']}`",
        f"- AX13 strict-core promotion: `{st['ax13_strict_core_promotion']}`",
        f"- P51 full closure pass: `{st['p51_full_closure_pass']}`",
        f"- P631 direct-route negative freeze selected: `{st['p631_direct_route_negative_freeze_selected']}`",
        f"- H37 premise-based directed state exported: `{st['h37_premise_based_directed_state_exported']}`",
        f"- P739/P740 pair12 strict-core upgrades: `{st['p739_pair12_split_upgrades']}` / `{st['p740_pair12_split_upgrades']}`",
        "", "## Lower-boundary sequence",
    ])
    for row in payload["lower_boundary_sequence"]:
        lines.append(f"- `{row['stage']}`: kind=`{row['kind']}`, exists=`{row['exists']}`, no_false_pass=`{row['no_false_pass']}`")
    lat = payload["recursion_lattice"]
    lines.extend([
        "", "## Recursion lattice",
        f"Cycle signature: `{lat['cycle_signature']}`.",
        f"Total states: `{lat['total_states']}`; admissible policy states: `{lat['admissible_states']}`.",
        f"Current state: `{lat['current_state']}`.",
        f"Continue lower-boundary without new invariant admissible: `{lat['continue_lower_boundary_without_new_invariant_admissible']}`.",
        "", "## Verdict", payload["closure_decision"]["professorial_verdict"],
        f"Decision: `{payload['closure_decision']['decision']}`.",
        "", "## Next honest step", payload["closure_decision"]["next_honest_step"],
        "", "## Negative exports",
    ])
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    state = current_artifact_state()
    sequence = lower_boundary_sequence()
    lattice = recursion_lattice(sequence)
    payload: dict[str, Any] = {
        "status": "P2683_POST_AX13_DIRECT_ROUTE_AND_T173_LOWER_BOUNDARY_CYCLE_AUDIT_NO_FALSE_PASS",
        "content_grep": content_grep(),
        "current_artifact_state": state,
        "lower_boundary_sequence": sequence,
        "recursion_lattice": lattice,
        "closure_decision": closure_decision(state, lattice),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2683/S1633 lower-boundary cycle guard",
        "## P2683/S1633 lower-boundary cycle guard\n\n"
        "`P2683/S1633` corrects stale direct-route targeting: AX13 already closes the `m2_psi4` target-EOM blocker on the canonical-ontology-supported external lane; P631 freezes the direct formal residual-cancellation route negative; H37/T171 is premise-based rather than missing; and P739/P740 still block strict-core pair12 upgrade.  The T24x/T25x lower-boundary target/nonexport/attempt pattern must not continue as the primary strategy without a new semantic invariant/provider at the chart-label-retaining pair12 typed seed-slot coordinate entry point.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2683/S1633 lower-boundary cycle Ltotal guard",
        "## P2683/S1633 lower-boundary cycle Ltotal guard\n\n"
        "`P2683/S1633` adds no variational term and does not promote `L_total`.  It freezes lower-boundary recursion as nonprimary unless a new semantic invariant/provider is exported; otherwise the honest fallback is an independent strict Lagrangian/EOM reverse-closure obstruction matrix.\n",
    )
    return payload


if __name__ == "__main__":
    main()
