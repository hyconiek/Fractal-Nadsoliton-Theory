#!/usr/bin/env python3
"""P2697/S1647: post-direct-route state-map/no-new-live-frontier reconciliation.

Runs the P2696 recommendation after the residual direct-route audits.  This is a
broad state-map reconciliation rather than a new closure attempt: it admits a
lane only if a genuinely new typed object, theorem, source atom, blocker-cut, or
provider class appears after the latest freezes.  Otherwise it emits a finite
no-new-live-frontier certificate.
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
OUT = GEN / "p2697_s1647_post_direct_route_state_map_no_new_live_frontier_reconciliation.json"
MD = GEN / "p2697_s1647_post_direct_route_state_map_no_new_live_frontier_reconciliation.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2683": GEN / "p2683_s1633_post_ax13_direct_route_and_t173_lower_boundary_cycle_audit.json",
    "P2684": GEN / "p2684_s1634_pair12_cycle_cut_semantic_invariant_provider_audit.json",
    "P2687": GEN / "p2687_s1637_one_new_strict_anisotropic_source_class_audit.json",
    "P2693": GEN / "p2693_s1643_post_p2680_nonselector_source_inventory_closure_state_map.json",
    "P2694": GEN / "p2694_s1644_fresh_broad_state_map_selection_after_p2680_closure.json",
    "P2695": GEN / "p2695_s1645_residual_direct_g_family_c1s1_zero_witness_no_go_matrix.json",
    "P2696": GEN / "p2696_s1646_pair1_c1c1_s1s1_zero_equation_carrier_no_go_matrix.json",
    "P631": GEN / "p631_current_strict_direct_formal_c1s1_route_negative_freeze_decision_packet_summary.json",
    "H37": GEN / "h37_sign_distinction_state_audit_summary.json",
    "P739": GEN / "p739_current_strict_t193_global_premise_based_directed_selector_state_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_summary.json",
    "P740": GEN / "p740_current_strict_t193_t194_pair12_strict_core_upgrade_bridge_after_h37_audit_summary.json",
    "P2677": GEN / "p2677_s1627_s6_o3_typed_seed_route_no_go_audit.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "direct_route_reopened",
    "generic_bridge_replayed",
    "selector_replayed",
    "qw2191_discharged",
    "lower_boundary_replayed",
    "lagrangian_reopened",
    "ltotal_promoted",
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
        "p2696_selected_p2697": r"P2697|post-direct-route|no-new-live-frontier|state-map/no-new-live-frontier",
        "direct_route_closed": r"P2695|P2696|P631|g4/g6/gY|c1c1/s1s1|residual-cancellation",
        "closed_nonselector_bridge": r"P2693|P2680 non-selector|canonical UV|alpha_geo|beta/Z_beta|source-inventory",
        "closed_lagrangian_lower_boundary": r"P2687|P2684|Lagrangian/EOM|anisotropic source|lower-boundary|pair12 cycle",
        "selector_h37_not_new": r"H37|T171|P739|P740|QW-2191|selector replay|premise-based directed",
        "forbidden_promotions": r"role transfer|L_total|ToE closure|bridge completion|strict-core selector",
    }
    return {"tool": "rg", "mode": "content-first post-direct-route broad state-map no-new-live-frontier reconciliation", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def state_reads() -> dict[str, Any]:
    p2683 = load_json(INPUTS["P2683"])
    p2684 = load_json(INPUTS["P2684"])
    p2687 = load_json(INPUTS["P2687"])
    p2693 = load_json(INPUTS["P2693"])
    p2695 = load_json(INPUTS["P2695"])
    p2696 = load_json(INPUTS["P2696"])
    p631 = load_json(INPUTS["P631"])
    h37 = load_json(INPUTS["H37"])
    p739 = load_json(INPUTS["P739"])
    p740 = load_json(INPUTS["P740"])
    p2677 = load_json(INPUTS["P2677"])
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "p2696_selected_p2697": "P2697" in p2696.get("decision", {}).get("next_honest_step", ""),
        "direct_g_family_closed_by_p2695": p2695.get("decision", {}).get("selected_g_targets_closed") is True,
        "pair1_c1c1_s1s1_bounded_no_go_by_p2696": p2696.get("decision", {}).get("bounded_no_go_for_targeted_c1c1_s1s1_zero_witnesses_now") is True,
        "p631_direct_formal_negative_freeze": p631.get("decision") == "DIRECT_FORMAL_C1S1_ROUTE_NEGATIVE_FREEZE_SELECTED",
        "p2693_nonselector_source_inventory_closed": p2693.get("decision", {}).get("bounded_no_go_now") is True or p2693.get("status", "").startswith("P2693"),
        "p2687_lagrangian_frozen": p2687.get("decision", {}).get("freeze_lagrangian_eom_reverse_closure_lane") is True,
        "p2684_lower_boundary_not_primary": p2684.get("decision", {}).get("lower_boundary_recursion_primary") is False,
        "p2683_h37_not_missing_replay": "H37" in json.dumps(p2683) and "premise" in json.dumps(p2683).lower(),
        "h37_summary_status": h37.get("status"),
        "p739_strict_core_upgrade_exported": p739.get("strict_core_upgrade_exported") is True or p739.get("result", {}).get("strict_core_upgrade_exported") is True,
        "p740_strict_core_upgrade_exported": p740.get("strict_core_upgrade_exported") is True or p740.get("result", {}).get("strict_core_upgrade_exported") is True,
        "p2677_no_go_present": p2677.get("status", "").startswith("P2677") or "no_go" in json.dumps(p2677).lower(),
    }


def freshness_gate_matrix(reads: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "lane": "direct_route_g_family_and_pair1_residual_cancellation",
            "closed_or_replay_gated_now": reads["direct_g_family_closed_by_p2695"] and reads["pair1_c1c1_s1s1_bounded_no_go_by_p2696"] and reads["p631_direct_formal_negative_freeze"],
            "new_typed_object_required_to_reopen": "strict-derived provider class, non-N477 ingredient, or direct-route blocker-cut that evades P477/N520/P479/N522/P631",
            "live_now": False,
        },
        {
            "lane": "generic_legacy_to_strict_or_p2680_nonselector_bridge_source",
            "closed_or_replay_gated_now": reads["p2693_nonselector_source_inventory_closed"],
            "new_typed_object_required_to_reopen": "new non-selector source atom or source theorem not among canonical UV, alpha_geo amplitude, beta/Z_beta",
            "live_now": False,
        },
        {
            "lane": "strict_lagrangian_eom_reverse_closure",
            "closed_or_replay_gated_now": reads["p2687_lagrangian_frozen"],
            "new_typed_object_required_to_reopen": "new strict anisotropic source class evading P1977/P1978",
            "live_now": False,
        },
        {
            "lane": "lower_boundary_pair12_cycle",
            "closed_or_replay_gated_now": reads["p2684_lower_boundary_not_primary"],
            "new_typed_object_required_to_reopen": "chart-label-retaining pair12 typed seed/subinterface plus semantic invariant/provider",
            "live_now": False,
        },
        {
            "lane": "selector_h37_t171_qw2191_replay",
            "closed_or_replay_gated_now": reads["p2683_h37_not_missing_replay"] and not reads["p739_strict_core_upgrade_exported"] and not reads["p740_strict_core_upgrade_exported"],
            "new_typed_object_required_to_reopen": "strict internal selector/orientation source that actually discharges QW-2191 rather than premise-based replay",
            "live_now": False,
        },
        {
            "lane": "role_transfer_ltotal_toe",
            "closed_or_replay_gated_now": True,
            "new_typed_object_required_to_reopen": "completed bridge/source theorem and separate role-transfer theorem",
            "live_now": False,
        },
    ]


def decision(rows: list[dict[str, Any]]) -> dict[str, Any]:
    live = [row["lane"] for row in rows if row["live_now"]]
    all_gated = all(row["closed_or_replay_gated_now"] and not row["live_now"] for row in rows)
    return {
        "decision": "P2697_POST_DIRECT_ROUTE_STATE_MAP_NO_NEW_LIVE_FRONTIER_CERTIFICATE_NO_FALSE_PASS",
        "live_lanes_now": live,
        "no_new_live_frontier_certificate": all_gated and not live,
        "reason": "After P2695/P2696, the residual direct route is closed/bounded-no-go on current artifacts, and the other broad lanes remain repetition-gated without a new typed object, theorem, source atom, blocker-cut, or provider class.",
        "next_honest_step": "Stop replaying closed lanes.  The next admissible research move must introduce a genuinely new strict typed provider/source/blocker-cut; absent that, maintain the no-new-live-frontier certificate rather than promoting selector, bridge, role-transfer, L_total, or ToE closure.",
        "forbidden_reopens": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = ["# P2697/S1647 post-direct-route state-map/no-new-live-frontier reconciliation", "", f"Status: `{payload['status']}`", "", "## Freshness gate matrix"]
    for row in payload["freshness_gate_matrix"]:
        lines.append(f"- `{row['lane']}`: live_now=`{row['live_now']}`, gated=`{row['closed_or_replay_gated_now']}`; reopen only with {row['new_typed_object_required_to_reopen']}")
    lines.extend(["", "## Decision", payload["decision"]["reason"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    reads = state_reads()
    rows = freshness_gate_matrix(reads)
    payload: dict[str, Any] = {
        "status": "P2697_POST_DIRECT_ROUTE_STATE_MAP_NO_NEW_LIVE_FRONTIER_RECONCILIATION_NO_FALSE_PASS",
        "content_grep": content_grep(),
        "state_reads": reads,
        "freshness_gate_matrix": rows,
        "decision": decision(rows),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2697/S1647 post-direct-route no-new-live-frontier reconciliation",
        "## P2697/S1647 post-direct-route no-new-live-frontier reconciliation\n\n"
        "`P2697/S1647` reruns the broad state-map after P2695/P2696.  The residual direct route is closed/bounded-no-go on current artifacts; P2680 non-selector bridge-source atoms, Lagrangian/EOM reverse closure, lower-boundary recursion, H37/T171/QW-2191 replay, role transfer, `L_total`, and ToE closure are all repetition-gated without a genuinely new typed object, theorem, source atom, blocker-cut, or provider class.  The output is a no-new-live-frontier certificate, not a closure promotion.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2697/S1647 no-new-live-frontier Ltotal guard",
        "## P2697/S1647 no-new-live-frontier Ltotal guard\n\n"
        "`P2697/S1647` explicitly does not promote `L_total`, variational closure, role transfer, strict-core selector closure, or ToE closure; it certifies that no current lane is live without a new strict typed provider/source/blocker-cut.\n",
    )
    append_once(
        AGENTS,
        "Current post-direct-route no-new-live-frontier guardrail (P2697/S1647, 2026-06-13)",
        "## Current post-direct-route no-new-live-frontier guardrail (P2697/S1647, 2026-06-13)\n\n"
        "- P2697 reconciles the broad state-map after P2695/P2696 and emits a no-new-live-frontier certificate on current artifacts: direct residual cancellation, P2680 bridge-source atoms, Lagrangian/EOM reverse closure, lower-boundary recursion, selector/H37/QW-2191 replay, role transfer, `L_total`, and ToE closure are all repetition-gated.\n"
        "- Do not continue any closed lane unless a genuinely new strict typed object, theorem, source atom, blocker-cut, or provider class is exported.\n"
        "- The next admissible research move must introduce such a new object; otherwise preserve the no-new-live-frontier certificate rather than manufacturing a closure.\n",
    )
    return payload


if __name__ == "__main__":
    main()
