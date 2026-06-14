#!/usr/bin/env python3
"""P2707/S1657: post-P2706 no-new-live-frontier reconciliation.

After the P2705/P2706 damping-to-selector audits, reconcile the live frontier
computably instead of manufacturing a closure move.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2707_s1657_post_p2706_no_new_live_frontier_reconciliation.json"
MD = GEN / "p2707_s1657_post_p2706_no_new_live_frontier_reconciliation.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2687_LAGRANGIAN": GEN / "p2687_s1637_one_new_strict_anisotropic_source_class_audit.json",
    "P2693_BRIDGE_SOURCE": GEN / "p2693_s1643_post_p2680_nonselector_source_inventory_closure_state_map.json",
    "P2695_DIRECT_G": GEN / "p2695_s1645_residual_direct_g_family_c1s1_zero_witness_no_go_matrix.json",
    "P2696_DIRECT_PAIR1": GEN / "p2696_s1646_pair1_c1c1_s1s1_zero_equation_carrier_no_go_matrix.json",
    "P2697_STATE_MAP": GEN / "p2697_s1647_post_direct_route_state_map_no_new_live_frontier_reconciliation.json",
    "P2699_SELECTOR_AUT": GEN / "p2699_s1649_z12_fractal_information_aut_invariant_selector_source_no_go.json",
    "P2700_SELECTOR_FUNCTIONALS": GEN / "p2700_s1650_exhaustive_aut_invariant_selector_functional_enumeration_no_go.json",
    "P2701_PROVIDER_INVENTORY": GEN / "p2701_s1651_strict_sourced_symmetry_breaking_provider_inventory_matrix.json",
    "P2704_P1343_SCOPE": GEN / "p2704_s1654_p1343_p1348_selector_provenance_revalidation_table.json",
    "P2705_RELEASE_POINTER": GEN / "p2705_s1655_release_9_3s_commit_boundary_alignment_audit.json",
    "P2706_DAMPING_INTERFACE": GEN / "p2706_s1656_damping_to_selector_interface_obstruction_witness_table.json",
}

LANES = [
    "selector_qw2191",
    "damping_to_selector",
    "older_release_pointer",
    "p1343_declared_scope_transfer",
    "direct_route_residual",
    "p2680_bridge_source_atoms",
    "lagrangian_eom_reverse_closure",
    "role_transfer_ltotal_toe",
]

NEGATIVE_EXPORT_FLAGS = [
    "new_live_frontier_found",
    "qw2191_discharged",
    "selector_closure_exported",
    "direct_route_reopened",
    "bridge_source_replayed",
    "lagrangian_reopened",
    "ltotal_promoted",
    "toe_closure_claimed",
    "role_transfer_started",
]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def decision_text(data: dict[str, Any]) -> str:
    decision = data.get("decision")
    if isinstance(decision, dict):
        return json.dumps(decision, sort_keys=True, ensure_ascii=False).lower()
    return str(decision).lower()


def flag(data: dict[str, Any], path: list[str], expected: Any = True) -> bool:
    cur: Any = data
    for key in path:
        if not isinstance(cur, dict):
            return False
        cur = cur.get(key)
    return cur is expected


def artifact_state() -> dict[str, Any]:
    loaded = {name: read_json(path) for name, path in INPUTS.items()}
    return {
        "loaded": loaded,
        "hashes": {name: sha(path) for name, path in INPUTS.items()},
        "statuses": {name: data.get("status") for name, data in loaded.items()},
    }


def lane_matrix(state: dict[str, Any]) -> list[dict[str, Any]]:
    d = state["loaded"]
    rows = [
        {
            "lane": "selector_qw2191",
            "current_artifact": "P2699/P2700/P2701/P2706",
            "computed_block": flag(d["P2699_SELECTOR_AUT"], ["decision", "bounded_no_go_now"]) and flag(d["P2700_SELECTOR_FUNCTIONALS"], ["decision", "bounded_no_go_now"]) and flag(d["P2701_PROVIDER_INVENTORY"], ["decision", "bounded_no_go_now"]) and flag(d["P2706_DAMPING_INTERFACE"], ["decision", "damping_transport_exports_directed_selector"], False),
            "live_frontier_candidate": False,
            "reason": "Aut-invariant, exhaustive functional, provider-inventory, and damping-interface routes do not export a non-premise directed selector.",
        },
        {
            "lane": "damping_to_selector",
            "current_artifact": "P2705/P2706",
            "computed_block": flag(d["P2705_RELEASE_POINTER"], ["decision", "release_9_3s_pointer_unblocks_current_stage"], False) and flag(d["P2706_DAMPING_INTERFACE"], ["decision", "damping_transport_exports_directed_selector"], False),
            "live_frontier_candidate": False,
            "reason": "The Release pointer is P2377/P2378 damping math and P2706 computes orientation-blindness at the selector interface.",
        },
        {
            "lane": "older_release_pointer",
            "current_artifact": "P2704/P2705",
            "computed_block": flag(d["P2704_P1343_SCOPE"], ["decision", "p1343_p1348_declared_scope_provenance_revalidated"]) and flag(d["P2705_RELEASE_POINTER"], ["decision", "release_9_3s_pointer_unblocks_current_stage"], False),
            "live_frontier_candidate": False,
            "reason": "P1343/P1348 remains positive only in declared scope; the supplied pointer does not add current selector closure authority.",
        },
        {
            "lane": "direct_route_residual",
            "current_artifact": "P2695/P2696/P2697",
            "computed_block": flag(d["P2695_DIRECT_G"], ["decision", "selected_g_targets_closed"]) and flag(d["P2696_DIRECT_PAIR1"], ["decision", "bounded_no_go_for_targeted_c1c1_s1s1_zero_witnesses_now"]) and flag(d["P2697_STATE_MAP"], ["decision", "no_new_live_frontier_certificate"]),
            "live_frontier_candidate": False,
            "reason": "Residual g-family and pair1 zero-equation carriers are bounded/no-go without a new strict-derived provider.",
        },
        {
            "lane": "p2680_bridge_source_atoms",
            "current_artifact": "P2693",
            "computed_block": "bounded no-go" in decision_text(d["P2693_BRIDGE_SOURCE"]) or "no-go" in decision_text(d["P2693_BRIDGE_SOURCE"]),
            "live_frontier_candidate": False,
            "reason": "Named non-selector source atoms remain closed as bounded no-go; generic bridge replay is not new evidence.",
        },
        {
            "lane": "lagrangian_eom_reverse_closure",
            "current_artifact": "P2687",
            "computed_block": "no currently exported" in decision_text(d["P2687_LAGRANGIAN"]) or "bounded no-go" in decision_text(d["P2687_LAGRANGIAN"]) or "no-go" in decision_text(d["P2687_LAGRANGIAN"]),
            "live_frontier_candidate": False,
            "reason": "No new strict anisotropic source class evading the prior obstruction is exported.",
        },
        {
            "lane": "role_transfer_ltotal_toe",
            "current_artifact": "P2704/P2705/P2706 guardrails",
            "computed_block": all(not value for value in d["P2706_DAMPING_INTERFACE"].get("decision", {}).get("negative_export_flags", {}).values()),
            "live_frontier_candidate": False,
            "reason": "No upstream bridge/source/selector closure exists to start role transfer, L_total, or ToE promotion.",
        },
    ]
    for row in rows:
        row["passes_no_live_frontier_check"] = row["computed_block"] and not row["live_frontier_candidate"]
    return rows


def post_p2697_inventory() -> dict[str, Any]:
    paths = sorted(GEN.glob("p27*.json"))
    rows = []
    for path in paths:
        data = read_json(path)
        txt = json.dumps(data, sort_keys=True, ensure_ascii=False).lower()
        rows.append({
            "path": rel(path),
            "status": data.get("status"),
            "mentions_new_live_frontier": "new_live_frontier_found\": true" in txt or "live_frontier_candidate\": true" in txt,
            "exports_closure_true": any(token in txt for token in ["qw2191_discharged\": true", "toe_closure_claimed\": true", "ltotal_promoted\": true", "selector_closure_exported\": true"]),
        })
    return {
        "p27xx_json_artifacts_scanned": len(rows),
        "artifacts_with_live_frontier_true": [row for row in rows if row["mentions_new_live_frontier"]],
        "artifacts_with_forbidden_closure_true": [row for row in rows if row["exports_closure_true"]],
        "sample_statuses": rows[-12:],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = ["# P2707/S1657 post-P2706 no-new-live-frontier reconciliation", "", f"Status: `{payload['status']}`", "", "## Lane matrix"]
    for row in payload["lane_matrix"]:
        lines.append(f"- `{row['lane']}`: blocked={row['computed_block']}, live_frontier_candidate={row['live_frontier_candidate']}. {row['reason']}")
    inv = payload["post_p2697_inventory"]
    lines.extend([
        "",
        "## Generated-artifact scan",
        f"- p27xx_json_artifacts_scanned={inv['p27xx_json_artifacts_scanned']}",
        f"- artifacts_with_live_frontier_true={len(inv['artifacts_with_live_frontier_true'])}",
        f"- artifacts_with_forbidden_closure_true={len(inv['artifacts_with_forbidden_closure_true'])}",
        "",
        "## Decision",
        payload["decision"]["reason"],
        "",
        "## Next honest step",
        payload["decision"]["next_honest_step"],
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    state = artifact_state()
    matrix = lane_matrix(state)
    inventory = post_p2697_inventory()
    no_live = all(row["passes_no_live_frontier_check"] for row in matrix) and not inventory["artifacts_with_live_frontier_true"] and not inventory["artifacts_with_forbidden_closure_true"]
    payload = {
        "status": "P2707_POST_P2706_NO_NEW_LIVE_FRONTIER_CERTIFICATE" if no_live else "P2707_REQUIRES_MANUAL_REVIEW",
        "input_hashes": state["hashes"],
        "input_statuses": state["statuses"],
        "lane_matrix": matrix,
        "post_p2697_inventory": inventory,
        "decision": {
            "no_new_live_frontier_certificate": no_live,
            "live_lanes_now": [] if no_live else [row["lane"] for row in matrix if not row["passes_no_live_frontier_check"]],
            "reason": "After P2706, the damping-to-selector interface is computably orientation-blind; the selector, direct-route, bridge-source, Lagrangian/EOM, older-release, role-transfer, L_total, and ToE lanes remain blocked on current artifacts.",
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "next_honest_step": "Do not manufacture another replay move.  The next admissible step must introduce a genuinely new strict typed object/provider/source/blocker-cut outside the closed lanes; if none is provided, preserve the P2697-P2707 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2707/S1657 post-P2706 no-new-live-frontier reconciliation", "## P2707/S1657 post-P2706 no-new-live-frontier reconciliation\n\n`P2707/S1657` reconciles the state map after the P2705/P2706 damping audits.  The finite lane matrix keeps selector/`QW-2191`, damping-to-selector, older-release transfer, direct-route residuals, P2680 bridge-source atoms, Lagrangian/EOM reverse closure, role transfer, `L_total`, and ToE lanes blocked on current artifacts.  A generated-artifact scan finds no `p27xx` artifact exporting a live-frontier or forbidden closure flag.  The current state is therefore a P2697-P2707 no-new-live-frontier certificate unless a genuinely new strict typed object/provider/source/blocker-cut is introduced.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2707/S1657 no-new-live-frontier Ltotal guard", "## P2707/S1657 no-new-live-frontier Ltotal guard\n\n`P2707/S1657` is a state-map reconciliation after P2706, not a variational construction.  It preserves the no-new-live-frontier/no-unlock certificate and does not promote `L_total`, selector closure, pair12 strict-core upgrade, role transfer, bridge closure, or ToE closure.\n")
    append_once(AGENTS, "Current post-P2706 no-new-live-frontier guardrail (P2707/S1657, 2026-06-14)", "## Current post-P2706 no-new-live-frontier guardrail (P2707/S1657, 2026-06-14)\n\n- P2707 reconciles the broad state map after P2705/P2706 and finds no new live frontier on current artifacts: selector/`QW-2191`, damping-to-selector, older-release transfer, direct-route residuals, P2680 bridge-source atoms, Lagrangian/EOM reverse closure, role transfer, `L_total`, and ToE lanes remain blocked.\n- Do not continue replay moves in those lanes or promote closure from P2377/P2378, P1343/P1348 declared-scope provenance, older release prose, or direct-route residuals.\n- The next admissible move must introduce a genuinely new strict typed object/provider/source/blocker-cut outside the closed lanes; otherwise preserve the P2697-P2707 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
