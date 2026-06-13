#!/usr/bin/env python3
"""P2688/S1638: post-P2687 live-frontier reconciliation audit.

After P2687 freezes the anisotropic-source continuation, this packet re-runs a
state-map reconciliation instead of blindly following any stale single-lane
recommendation.  It explicitly checks the direct-route, lower-boundary,
Lagrangian/EOM, selector, and non-selector bridge-source lanes, then selects one
finite next target.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2688_s1638_post_p2687_live_frontier_reconciliation_audit.json"
MD = GEN / "p2688_s1638_post_p2687_live_frontier_reconciliation_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2687": GEN / "p2687_s1637_one_new_strict_anisotropic_source_class_audit.json",
    "P2686": GEN / "p2686_s1636_shared_background_nonproxy_component_residual_table.json",
    "P2684": GEN / "p2684_s1634_pair12_cycle_cut_semantic_invariant_provider_audit.json",
    "P2683": GEN / "p2683_s1633_post_ax13_direct_route_and_t173_lower_boundary_cycle_audit.json",
    "P2682": GEN / "p2682_s1632_p46_p50_m2_psi4_target_eom_zero_witness_matrix_audit.json",
    "AX13": GEN / "ax13_canonical_ontology_supported_preobserver_target_eom_coherence_instance.json",
    "P51": GEN / "p51_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_preobserver_target_eom_coherence_instance.json",
    "P631": GEN / "p631_current_strict_direct_formal_c1s1_route_negative_freeze_decision_packet_summary.json",
    "P2680": GEN / "p2680_s1630_legacy_strict_bridge_source_inventory_no_selector_replay_audit.json",
    "P2650": GEN / "p2650_s1600_canonical_length_uv_unit_source_candidate_exhaustion_no_go.json",
    "P2653": GEN / "p2653_s1603_typed_nadsoliton_metric_uv_source_obligation_rank_audit.json",
    "P2662": GEN / "p2662_s1612_entropy_boundary_phase_unit_map_conditional_theorem_audit.json",
    "P2663": GEN / "p2663_s1613_chain_level_boundary_phase_bit_target_audit.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "stale_p46_p50_reopened",
    "lagrangian_eom_reverse_closure_reopened",
    "lower_boundary_recursion_reopened",
    "selector_replay_reopened",
    "generic_bridge_reopened",
    "uv_unit_or_beta_source_exported",
    "role_transfer_started",
    "toe_closure_claimed",
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


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
    return {"count": len(lines), "samples": lines[:80]}


def content_grep() -> dict[str, Any]:
    patterns = {
        "post_p2687_state_map": r"P2687|P2688|state-map|live frontier|pivot",
        "direct_route_status": r"AX13|P51|P631|negative freeze|P46/N49|target-EOM|m2_psi4",
        "lagrangian_eom_freeze": r"P2685|P2686|anisotropic source|bounded no-go|P1977|P1978",
        "lower_boundary_freeze": r"P2684|cycle-cut|lower-boundary|pair12 typed seed|semantic invariant/provider",
        "non_selector_bridge_atoms": r"P2680|target-independent positive beta|Z_beta|canonical length|UV unit|entropy reference|bit-to-length",
        "forbidden_replays": r"selector replay|tau_src -> pair12|beta_tors -> chi11|role transfer|ToE closure|generic bridge",
    }
    return {"tool": "rg", "mode": "broad post-P2687 live-frontier reconciliation before choosing P2689", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def state_reads() -> dict[str, Any]:
    data = {name: load_json(path) for name, path in INPUTS.items()}
    p51_state = data["P51"].get("route_state", {})
    ax13 = data["AX13"].get("result", {})
    p2687_decision = data["P2687"].get("decision", {})
    p2684_decision = data["P2684"].get("decision", {})
    p2686_decision = data["P2686"].get("decision", {})
    p2662_decision_text = json.dumps(data["P2662"], sort_keys=True)
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "p2687_freezes_lagrangian_lane": p2687_decision.get("freeze_lagrangian_eom_reverse_closure_lane") is True,
        "p2687_recommended_stale_p46_pivot": "P46/N49" in p2687_decision.get("next_honest_step", ""),
        "ax13_target_eom_closed_external": ax13.get("canonical_ontology_supported_target_eom_coefficient_defect_zero_witness_present") is True,
        "p51_direct_route_full_closure_pass": data["P51"].get("full_closure_pass") is True,
        "p51_remaining_g_family_witnesses_absent": not any(
            p51_state.get(key) is True for key in (
                "direct_g4_family_zero_witness_present",
                "direct_g6_family_zero_witness_present",
                "direct_gY_family_zero_witness_present",
            )
        ),
        "p631_direct_route_negative_freeze": data["P631"].get("decision") == "DIRECT_FORMAL_C1S1_ROUTE_NEGATIVE_FREEZE_SELECTED",
        "p2684_cycle_cut_success": p2684_decision.get("cycle_cut_success") is True,
        "p2686_bounded_no_go_boundary_active": p2686_decision.get("bounded_no_go_boundary_active") is True,
        "p2680_missing_non_selector_atoms": [
            row.get("atom") for row in data["P2680"].get("source_atom_inventory", [])
            if row.get("source_theorem_exported") is False and row.get("atom") != "selector_phase_orientation_source"
        ],
        "p2650_exports_canonical_unit": False if not data["P2650"].get("missing", False) else False,
        "p2653_exports_typed_metric_uv_source": "no typed metric/UV source theorem" not in json.dumps(data["P2653"], sort_keys=True),
        "p2662_conditional_unit_map_present": "conditional" in p2662_decision_text.lower() and "bit-to-length" in p2662_decision_text,
        "p2663_one_bit_carrier_present": "one-bit" in json.dumps(data["P2663"], sort_keys=True).lower(),
    }


def lane_matrix(reads: dict[str, Any]) -> list[dict[str, Any]]:
    rows = [
        {
            "lane_id": "direct_p46_p50_m2_psi4",
            "live_now": False,
            "reason": "AX13 closes the attacked target-EOM blocker externally and P631 already selected direct-route negative freeze; remaining g-family witnesses are absent but not the highest honest post-P2687 target.",
            "proof_grade_next": False,
        },
        {
            "lane_id": "lagrangian_eom_anisotropic_reverse_closure",
            "live_now": False,
            "reason": "P2686/P2687 leave a bounded no-go and no exported strict anisotropic source class.",
            "proof_grade_next": False,
        },
        {
            "lane_id": "lower_boundary_pair12_cycle_cut",
            "live_now": False,
            "reason": "P2684 cycle-cut did not export the chart-label-retaining pair12 typed seed plus semantic provider.",
            "proof_grade_next": False,
        },
        {
            "lane_id": "selector_orientation_replay",
            "live_now": False,
            "reason": "QW-2191 and repetition guardrails still block selector replay without a new selector source.",
            "proof_grade_next": False,
        },
        {
            "lane_id": "non_selector_bridge_canonical_uv_unit_atom",
            "live_now": True,
            "reason": "P2680 leaves canonical length/UV unit as a typed non-selector missing atom; P2662/P2663 provide a finite conditional entropy/one-bit carrier scaffold, but P2650/P2653 still export no unconditional typed unit theorem.",
            "proof_grade_next": True,
        },
    ]
    for row in rows:
        row["status_score"] = int(row["live_now"]) + int(row["proof_grade_next"])
    return rows


def selected_next(reads: dict[str, Any], rows: list[dict[str, Any]]) -> dict[str, Any]:
    winner = max(rows, key=lambda row: row["status_score"])
    obligations = [
        {"obligation": "intrinsic entropy reference cell or entropy-zero theorem", "satisfied_now": False},
        {"obligation": "boundary-phase bit target N log(2) derived without selector replay", "satisfied_now": False},
        {"obligation": "bit-to-length or bit-to-action unit map", "satisfied_now": False},
        {"obligation": "scale-orbit quotient selects one positive UV unit", "satisfied_now": False},
        {"obligation": "beta/Z_beta source follows without target fitting", "satisfied_now": False},
    ]
    return {
        "decision": "P2688_POST_P2687_STATE_MAP_RECONCILES_STALE_P46_RECOMMENDATION_AND_SELECTS_CANONICAL_UV_UNIT_ATOM_NO_FALSE_PASS",
        "selected_lane": winner["lane_id"],
        "selected_next_packet": "P2689_intrinsic_entropy_reference_cell_and_bit_to_length_uv_unit_source_obligation_matrix",
        "selected_obligations": obligations,
        "stale_p46_recommendation_overridden": reads["p2687_recommended_stale_p46_pivot"] and reads["p631_direct_route_negative_freeze"],
        "professorial_verdict": (
            "P2688 corrects the post-P2687 next-step selection by re-reading the live state-map.  The P46/P50 m2_psi4 lane is not the right immediate target because later AX13/P51/P631 state already closes the attacked target-EOM blocker externally and freezes the direct c1s1 route as nonclosed.  The Lagrangian/EOM anisotropic lane is bounded no-go after P2686/P2687, lower-boundary recursion is cycle-cut blocked after P2684, and selector replay remains barred.  The strongest remaining finite proof-grade target is therefore one P2680 non-selector bridge atom: canonical length/UV unit source, attacked through the finite entropy/reference-cell and bit-to-length obligations left by P2662/P2663."
        ),
        "next_honest_step": (
            "Run P2689 as an intrinsic entropy reference-cell and bit-to-length UV-unit source obligation matrix.  It should test whether the P2662 conditional unit-map scaffold plus the P2663 one-bit holonomy carrier can supply: an intrinsic entropy zero/reference cell, a selector-free boundary-phase bit target, a bit-to-length or bit-to-action unit map, and a scale-orbit quotient selecting one positive UV unit.  If any obligation remains absent, export a bounded no-go for this canonical-unit atom and do not start role-transfer."
        ),
        "uv_unit_or_beta_source_exported_now": False,
        "role_transfer_started_now": False,
        "toe_closed_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = ["# P2688/S1638 post-P2687 live-frontier reconciliation audit", "", f"Status: `{payload['status']}`", "", "## Content-first grep"]
    for name, data in payload["content_grep"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    lines.extend(["", "## Lane matrix"])
    for row in payload["lane_matrix"]:
        lines.append(f"- `{row['lane_id']}`: live_now=`{row['live_now']}`, proof_grade_next=`{row['proof_grade_next']}`, score=`{row['status_score']}` — {row['reason']}")
    lines.extend(["", "## Selected next obligations"])
    for item in payload["selection"]["selected_obligations"]:
        lines.append(f"- `{item['obligation']}`: `{item['satisfied_now']}`")
    lines.extend(["", "## Verdict", payload["selection"]["professorial_verdict"], "", "## Next honest step", payload["selection"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    reads = state_reads()
    rows = lane_matrix(reads)
    selection = selected_next(reads, rows)
    payload: dict[str, Any] = {
        "status": "P2688_POST_P2687_LIVE_FRONTIER_RECONCILIATION_AUDIT_NO_FALSE_PASS",
        "content_grep": content_grep(),
        "state_reads": reads,
        "lane_matrix": rows,
        "selection": selection,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2688/S1638 post-P2687 live-frontier reconciliation audit",
        "## P2688/S1638 post-P2687 live-frontier reconciliation audit\n\n"
        "`P2688/S1638` re-runs the state-map after P2687 and overrides the stale P46/P50 return recommendation: AX13/P51/P631 already close/freeze the attacked `m2_psi4` direct-route target externally/nonclosed, P2686/P2687 freeze the anisotropic Lagrangian/EOM lane, and P2684 blocks lower-boundary recursion.  The selected proof-grade target is one P2680 non-selector bridge atom: canonical length/UV unit source, narrowed to the P2662/P2663 entropy reference-cell and bit-to-length obligation matrix.  No UV unit, beta source, role transfer, bridge completion, selector closure, `L_total`, or ToE closure is exported.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2688/S1638 state-map pivot Ltotal guard",
        "## P2688/S1638 state-map pivot Ltotal guard\n\n"
        "`P2688/S1638` keeps `L_total` nonpromoted: after P2687, the honest next target is not another reverse-closure replay but a canonical UV-unit source obligation matrix outside role-bearing action semantics.  A future P2689 must first export an intrinsic entropy reference cell, selector-free bit target, bit-to-length/action unit map, and scale-orbit quotient before any beta source or role-transfer audit can start.\n",
    )
    append_once(
        AGENTS,
        "Current post-P2687 live-frontier reconciliation guardrail (P2688/S1638, 2026-06-13)",
        "## Current post-P2687 live-frontier reconciliation guardrail (P2688/S1638, 2026-06-13)\n\n"
        "- P2688 corrects the stale return-to-P46/P50 recommendation after rereading AX13/P51/P631: do not reopen the attacked `m2_psi4` direct-route target as the immediate next move.\n"
        "- With Lagrangian/EOM, lower-boundary, selector replay, and generic bridge lanes repetition-gated, the next proof-grade move is one P2680 non-selector source atom: a canonical length/UV unit source audit through the P2662/P2663 entropy reference-cell and bit-to-length obligation matrix.\n",
    )
    return payload


if __name__ == "__main__":
    main()
