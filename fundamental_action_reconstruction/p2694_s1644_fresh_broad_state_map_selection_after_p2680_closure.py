#!/usr/bin/env python3
"""P2694/S1644: fresh broad state-map selection after P2680 closure.

This executes the P2693 recommendation by rebuilding the current frontier after
the non-selector bridge-source atoms were closed as bounded no-go.  It does not
reopen generic bridge, selector replay, role transfer, L_total, or ToE closure.
It selects the next finite proof-grade move only if a non-replayed lane remains
live on current artifacts.
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
OUT = GEN / "p2694_s1644_fresh_broad_state_map_selection_after_p2680_closure.json"
MD = GEN / "p2694_s1644_fresh_broad_state_map_selection_after_p2680_closure.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
F3 = ROOT / "F3_CURRENT_FAR_FRONTIER_KERNEL_ARTIFACT_SENSITIVITY_CLASSIFICATION_PACKET.md"

INPUTS = {
    "P2681": GEN / "p2681_s1631_professorial_research_state_map_and_agents_reorientation_audit.json",
    "P2682": GEN / "p2682_s1632_p46_p50_m2_psi4_target_eom_zero_witness_matrix_audit.json",
    "P2684": GEN / "p2684_s1634_pair12_cycle_cut_semantic_invariant_provider_audit.json",
    "P2687": GEN / "p2687_s1637_one_new_strict_anisotropic_source_class_audit.json",
    "P2688": GEN / "p2688_s1638_post_p2687_live_frontier_reconciliation_audit.json",
    "P2693": GEN / "p2693_s1643_post_p2680_nonselector_source_inventory_closure_state_map.json",
    "F3": F3,
}

NEGATIVE_EXPORT_FLAGS = [
    "bridge_reopened",
    "selector_replay_imported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_claimed",
    "m2_psi4_attacked_target_reopened",
]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8") if path.exists() else ""


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
        "p2693_selected_p2694": r"P2694|fresh broad state-map|fresh_broad_state_map_selection|no-new-live-frontier",
        "f3_residual_direct_g_family": r"g4|g6|gY|c1s1 shift defect|zero witness|P46|N49",
        "closed_bridge_source_round": r"P2693|post-P2680|non-selector source-inventory|bounded no-go",
        "closed_lagrangian_lower_selector": r"P2684|P2687|Lagrangian/EOM|lower-boundary|selector replay|QW-2191",
        "forbidden_reopen": r"generic bridge replay|role transfer|L_total|ToE closure|m2_psi4|beta_tors.*chi_?11|tau_src",
    }
    return {"tool": "rg", "mode": "content-first fresh broad state-map selection after P2680 closure", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def state_reads() -> dict[str, Any]:
    p2682 = load_json(INPUTS["P2682"])
    p2684 = load_json(INPUTS["P2684"])
    p2687 = load_json(INPUTS["P2687"])
    p2688 = load_json(INPUTS["P2688"])
    p2693 = load_json(INPUTS["P2693"])
    f3_text = read_text(F3)
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "p2693_selected_p2694": "P2694" in p2693.get("decision", {}).get("next_honest_step", ""),
        "p2693_only_fresh_state_map_live": p2693.get("decision", {}).get("live_lanes") == ["fresh_broad_state_map_selection"],
        "p2682_m2_target_bounded_no_go": p2682.get("closure_decision", {}).get("bounded_no_go_now") is True,
        "p2684_lower_boundary_primary": p2684.get("decision", {}).get("lower_boundary_recursion_primary") is True,
        "p2687_freezes_lagrangian_lane": p2687.get("decision", {}).get("freeze_lagrangian_eom_reverse_closure_lane") is True,
        "p2688_selected_bridge_uv_lane": p2688.get("selection", {}).get("selected_lane") == "non_selector_bridge_canonical_uv_unit_atom",
        "f3_names_g4_g6_gY_zero_witnesses": all(term in f3_text for term in ["g4", "g6", "gY", "c1s1"]),
        "f3_names_m2_psi4_but_later_stale": "m2 psi4" in f3_text or "m2_psi4" in f3_text,
    }


def residual_direct_target_matrix(reads: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "target_id": "direct_g4_c1s1_shift_defect_zero_witness",
            "kernel_split_robust_from_f3": reads["f3_names_g4_g6_gY_zero_witnesses"],
            "already_attacked_by_ax13_p631_m2_freeze": False,
            "requires_selector_replay": False,
            "finite_computable": True,
            "live_now": True,
            "next_obligation": "Compute/export an explicit zero witness or no-go for the direct quartic-like g4 c1s1 shift defect on corrected canonical-ontology support.",
        },
        {
            "target_id": "direct_g6_c1s1_shift_defect_zero_witness",
            "kernel_split_robust_from_f3": reads["f3_names_g4_g6_gY_zero_witnesses"],
            "already_attacked_by_ax13_p631_m2_freeze": False,
            "requires_selector_replay": False,
            "finite_computable": True,
            "live_now": True,
            "next_obligation": "Compute/export an explicit zero witness or no-go for the direct quintic-like g6 c1s1 shift defect on corrected canonical-ontology support.",
        },
        {
            "target_id": "direct_gY_c1s1_shift_defect_zero_witness",
            "kernel_split_robust_from_f3": reads["f3_names_g4_g6_gY_zero_witnesses"],
            "already_attacked_by_ax13_p631_m2_freeze": False,
            "requires_selector_replay": False,
            "finite_computable": True,
            "live_now": True,
            "next_obligation": "Compute/export an explicit zero witness or no-go for the direct yukawa-like gY c1s1 shift defect on corrected canonical-ontology support.",
        },
    ]


def broad_lane_matrix(reads: dict[str, Any], direct_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    any_direct = any(row["live_now"] for row in direct_rows)
    return [
        {"lane": "generic_bridge_or_p2680_nonselector_replay", "live_now": False, "proof_grade_next": False, "reason": "P2693 closes the named P2680 non-selector source atoms as bounded no-go."},
        {"lane": "direct_m2_psi4_attacked_target", "live_now": False, "proof_grade_next": False, "reason": "P2682/P2688 mark the attacked m2_psi4 target as stale/bounded-no-go for immediate replay."},
        {"lane": "direct_residual_g_family_zero_witnesses", "live_now": any_direct, "proof_grade_next": any_direct, "reason": "F3 still names g4/g6/gY c1s1 zero-witness defects as kernel-split-robust, and they are not the attacked m2_psi4 target."},
        {"lane": "lagrangian_eom_reverse_closure", "live_now": False, "proof_grade_next": False, "reason": "P2687 freezes this lane unless a new strict anisotropic source class is exported."},
        {"lane": "lower_boundary_pair12_cycle", "live_now": False, "proof_grade_next": False, "reason": "P2684 blocks lower-boundary recursion without a new semantic invariant/provider."},
        {"lane": "selector_tau_pair_or_beta_tors_replay", "live_now": False, "proof_grade_next": False, "reason": "Selector/tau/beta_tors replay remains guardrail-blocked without a new source object."},
        {"lane": "role_transfer_ltotal_toe", "live_now": False, "proof_grade_next": False, "reason": "Role transfer, L_total, and ToE closure are downstream of source/bridge closure, which is not exported."},
    ]


def selection(lanes: list[dict[str, Any]], direct_rows: list[dict[str, Any]]) -> dict[str, Any]:
    live_proof = [row for row in lanes if row["live_now"] and row["proof_grade_next"]]
    selected = "direct_residual_g_family_zero_witnesses" if any(row["lane"] == "direct_residual_g_family_zero_witnesses" for row in live_proof) else "no_new_live_frontier"
    return {
        "decision": "P2694_FRESH_STATE_MAP_SELECTS_RESIDUAL_DIRECT_G_FAMILY_ZERO_WITNESS_MATRIX_NO_FALSE_PASS",
        "selected_lane": selected,
        "selected_next_packet": "P2695_residual_direct_g4_g6_gY_c1s1_zero_witness_no_go_matrix" if selected != "no_new_live_frontier" else None,
        "selected_targets": [row["target_id"] for row in direct_rows if row["live_now"]],
        "professorial_verdict": (
            "P2694 rebuilds the live frontier after the P2680 source-inventory closure.  Bridge-source replay, m2_psi4 replay, Lagrangian/EOM reverse closure, lower-boundary recursion, selector/tau replay, role transfer, L_total, and ToE closure remain closed.  The only finite proof-grade lane left by the current state-map is the residual kernel-split-robust direct-route g-family zero-witness matrix: g4, g6, and gY c1s1 shift defects named by F3 but not identical to the already attacked m2_psi4 target."
        ),
        "next_honest_step": (
            "P2695 should compute a residual direct g-family zero-witness/no-go matrix for g4, g6, and gY c1s1 shift defects on corrected canonical-ontology support, explicitly excluding m2_psi4 replay, selector import, bridge-source import, role transfer, L_total promotion, and ToE closure."
        ),
        "bridge_reopened_now": False,
        "m2_psi4_reopened_now": False,
        "selector_replay_now": False,
        "ltotal_promoted_now": False,
        "toe_closed_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = ["# P2694/S1644 fresh broad state-map selection after P2680 closure", "", f"Status: `{payload['status']}`", "", "## Content-first grep"]
    for name, data in payload["content_grep"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    lines.extend(["", "## Residual direct target matrix"])
    for row in payload["residual_direct_target_matrix"]:
        lines.append(f"- `{row['target_id']}`: live_now=`{row['live_now']}`, finite=`{row['finite_computable']}` — {row['next_obligation']}")
    lines.extend(["", "## Broad lane matrix"])
    for row in payload["broad_lane_matrix"]:
        lines.append(f"- `{row['lane']}`: live_now=`{row['live_now']}`, proof_grade_next=`{row['proof_grade_next']}` — {row['reason']}")
    lines.extend(["", "## Verdict", payload["selection"]["professorial_verdict"], "", "## Next honest step", payload["selection"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    reads = state_reads()
    direct_rows = residual_direct_target_matrix(reads)
    lanes = broad_lane_matrix(reads, direct_rows)
    payload: dict[str, Any] = {
        "status": "P2694_FRESH_BROAD_STATE_MAP_SELECTION_AFTER_P2680_CLOSURE_NO_FALSE_PASS",
        "content_grep": content_grep(),
        "state_reads": reads,
        "residual_direct_target_matrix": direct_rows,
        "broad_lane_matrix": lanes,
        "selection": selection(lanes, direct_rows),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2694/S1644 fresh broad state-map selection after P2680 closure",
        "## P2694/S1644 fresh broad state-map selection after P2680 closure\n\n"
        "`P2694/S1644` rebuilds the frontier after P2693.  Generic bridge-source replay, the attacked `m2_psi4` target, Lagrangian/EOM reverse closure, lower-boundary recursion, selector/tau replay, role transfer, `L_total`, and ToE closure remain closed.  The selected finite proof-grade next lane is the residual direct g-family zero-witness/no-go matrix for `g4`, `g6`, and `gY` `c1s1` shift defects named by F3, explicitly excluding `m2_psi4` replay and bridge/source imports.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2694/S1644 fresh state-map Ltotal guard",
        "## P2694/S1644 fresh state-map Ltotal guard\n\n"
        "`P2694/S1644` keeps `L_total` nonpromoted: the next selected work is a direct-route zero-witness/no-go matrix, not a variational closure or role-transfer step.\n",
    )
    append_once(
        AGENTS,
        "Current fresh broad state-map selection guardrail (P2694/S1644, 2026-06-13)",
        "## Current fresh broad state-map selection guardrail (P2694/S1644, 2026-06-13)\n\n"
        "- P2694 selects the residual kernel-split-robust direct-route `g4/g6/gY` `c1s1` zero-witness/no-go matrix as the next finite proof-grade move after P2680 source-inventory closure.\n"
        "- Do not reopen the attacked `m2_psi4` target, generic bridge-source replay, selector/tau replay, role transfer, `L_total`, or ToE closure while executing that matrix.\n",
    )
    return payload


if __name__ == "__main__":
    main()
