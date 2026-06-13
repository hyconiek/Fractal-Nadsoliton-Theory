#!/usr/bin/env python3
"""P2693/S1643: post-P2680 non-selector source-inventory closure/state-map.

This executes the post-P2692 recommendation: reconcile the P2680 non-selector
bridge-source inventory after the canonical UV-unit, boundary-sector provider,
alpha_geo amplitude, and beta/Z_beta source audits.  The goal is to compute the
current bridge-source bit state and decide whether any non-selector bridge atom
remains honestly live without importing generic bridge replay, selector replay,
role transfer, L_total, or ToE closure.
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
OUT = GEN / "p2693_s1643_post_p2680_nonselector_source_inventory_closure_state_map.json"
MD = GEN / "p2693_s1643_post_p2680_nonselector_source_inventory_closure_state_map.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

INPUTS = {
    "P2680": GEN / "p2680_s1630_legacy_strict_bridge_source_inventory_no_selector_replay_audit.json",
    "P2688": GEN / "p2688_s1638_post_p2687_live_frontier_reconciliation_audit.json",
    "P2689": GEN / "p2689_s1639_entropy_reference_cell_bit_to_length_uv_unit_obligation_matrix.json",
    "P2690": GEN / "p2690_s1640_selector_free_nonexact_boundary_phase_sector_provider_audit.json",
    "P2691": GEN / "p2691_s1641_alpha_geo_role_safe_amplitude_source_audit.json",
    "P2692": GEN / "p2692_s1642_target_independent_positive_beta_zbeta_source_audit.json",
}

NEGATIVE_EXPORT_FLAGS = [
    "p2680_nonselector_inventory_closed_as_pass",
    "bridge_completion_claimed",
    "role_transfer_started",
    "selector_replay_imported",
    "canonical_uv_unit_replayed",
    "beta_tors_chi11_imported",
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
        "p2692_selected_p2693": r"P2693|post-P2680|source-inventory closure|state-map reconciliation",
        "p2680_nonselector_atoms": r"amplitude_role_safe_source|canonical_length_or_uv_unit_source|target_independent_positive_beta_or_z_beta_source",
        "uv_unit_route_freeze": r"P2689|P2690|canonical UV|boundary-phase sector|entropy/UV-unit route|bounded no-go",
        "amplitude_beta_freeze": r"P2691|P2692|alpha_geo amplitude|beta/Z_beta|bounded no-go",
        "forbidden_replay_imports": r"generic bridge replay|bridge completion|selector replay|role transfer|canonical UV-unit replay|beta_tors.*chi_?11|L_total|ToE closure",
    }
    return {"tool": "rg", "mode": "content-first post-P2680 non-selector source-inventory closure/state-map", "patterns": {k: rg_count(v) for k, v in patterns.items()}}


def state_reads() -> dict[str, Any]:
    data = {name: load_json(path) for name, path in INPUTS.items()}
    p2680_atoms = {row.get("atom"): row for row in data["P2680"].get("source_atom_inventory", [])}
    p2689_decision = data["P2689"].get("decision", {})
    p2690_decision = data["P2690"].get("decision", {})
    p2691_decision = data["P2691"].get("decision", {})
    p2692_decision = data["P2692"].get("decision", {})
    return {
        "hashes": {name: sha256_file(path) for name, path in INPUTS.items()},
        "p2692_selected_p2693": "P2693" in p2692_decision.get("next_honest_step", ""),
        "p2680_amplitude_source_exported": p2680_atoms.get("amplitude_role_safe_source", {}).get("source_theorem_exported") is True,
        "p2680_canonical_uv_source_exported": p2680_atoms.get("canonical_length_or_uv_unit_source", {}).get("source_theorem_exported") is True,
        "p2680_positive_beta_source_exported": p2680_atoms.get("target_independent_positive_beta_or_z_beta_source", {}).get("source_theorem_exported") is True,
        "p2689_uv_unit_or_beta_source_exported": p2689_decision.get("uv_unit_or_beta_source_exported_now") is True,
        "p2690_freezes_entropy_uv_route": p2690_decision.get("freeze_entropy_uv_unit_route") is True,
        "p2691_bounds_alpha_atom": p2691_decision.get("bounded_no_go_for_current_alpha_geo_amplitude_atom") is True,
        "p2692_bounds_beta_atom": p2692_decision.get("bounded_no_go_for_current_beta_zbeta_atom") is True,
        "p2692_beta_source_exported": p2692_decision.get("target_independent_positive_beta_or_zbeta_source_exported_now") is True,
    }


def nonselector_atom_status(reads: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "atom": "canonical_length_or_uv_unit_source",
            "component": "damping/compression",
            "audited_by": ["P2689", "P2690"],
            "formal_near_miss_present": True,
            "source_exported_now": reads["p2689_uv_unit_or_beta_source_exported"],
            "bounded_no_go_now": reads["p2690_freezes_entropy_uv_route"],
            "blocking_reason": "Conditional entropy scale selection and one-bit carrier exist, but no selector-free nonexact sector provider or bit-to-length/action map is exported.",
        },
        {
            "atom": "amplitude_role_safe_source",
            "component": "amplitude/normalization",
            "audited_by": ["P2691"],
            "formal_near_miss_present": True,
            "source_exported_now": reads["p2680_amplitude_source_exported"],
            "bounded_no_go_now": reads["p2691_bounds_alpha_atom"],
            "blocking_reason": "Strict alpha_geo and exact scalar-shape normalization exist, but no amplitude absorption bridge source, physical-role safety theorem, or APD/Lagrangian source is exported.",
        },
        {
            "atom": "target_independent_positive_beta_or_z_beta_source",
            "component": "damping/compression",
            "audited_by": ["P2692"],
            "formal_near_miss_present": True,
            "source_exported_now": reads["p2692_beta_source_exported"],
            "bounded_no_go_now": reads["p2692_bounds_beta_atom"],
            "blocking_reason": "Positive beta orbit and target-dependent inversion exist, but no target-independent Z_beta/beta source theorem or unit/conservation identity is exported.",
        },
    ]


def bridge_source_lattice(rows: list[dict[str, Any]]) -> dict[str, Any]:
    obligations = [row["atom"] for row in rows]
    current_state = {row["atom"]: row["source_exported_now"] for row in rows}
    lattice_rows = []
    passing_states = 0
    for bits in itertools.product([False, True], repeat=len(obligations)):
        state = dict(zip(obligations, bits))
        passes = all(state.values())
        passing_states += int(passes)
        hamming_from_current = sum(state[key] != current_state[key] for key in obligations)
        if passes or hamming_from_current <= 1 or state == current_state:
            lattice_rows.append({"state": state, "passes_nonselector_bridge_source_gate": passes, "hamming_from_current": hamming_from_current})
    missing = [key for key, value in current_state.items() if not value]
    return {
        "obligations": obligations,
        "total_states": 2 ** len(obligations),
        "passing_states": passing_states,
        "current_state": current_state,
        "missing_current_obligations": missing,
        "hamming_distance_to_nonselector_bridge_source_pass": len(missing),
        "distinguished_rows": lattice_rows,
        "nonselector_bridge_source_gate_passes_now": all(current_state.values()),
    }


def lane_reconciliation(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    all_bridge_atoms_frozen = all(row["bounded_no_go_now"] and not row["source_exported_now"] for row in rows)
    return [
        {
            "lane": "generic_legacy_to_strict_bridge_completion",
            "live_now": False,
            "reason": "P2680 converted the generic bridge into finite atoms; P2689-P2692 leave the missing non-selector atoms bounded no-go.",
        },
        {
            "lane": "p2680_nonselector_bridge_source_atoms",
            "live_now": not all_bridge_atoms_frozen,
            "reason": "All currently named non-selector source atoms are either formal near-misses or bounded no-go; a new typed atom would be needed.",
        },
        {
            "lane": "selector_or_tau_pair_replay",
            "live_now": False,
            "reason": "QW-2191, selector replay, tau_src->pair12, and beta_tors->chi11 remain explicitly forbidden without a new source object.",
        },
        {
            "lane": "role_transfer_or_ltotal_promotion",
            "live_now": False,
            "reason": "Role transfer and L_total promotion are downstream of bridge/source closure, which is not exported.",
        },
        {
            "lane": "fresh_broad_state_map_selection",
            "live_now": True,
            "reason": "The admissible next move is to rebuild the live frontier from current state-map evidence and require a new typed object/theorem/source atom before reopening closed lanes.",
        },
    ]


def decision(lattice: dict[str, Any], lanes: list[dict[str, Any]]) -> dict[str, Any]:
    live_lanes = [row["lane"] for row in lanes if row["live_now"]]
    return {
        "decision": "P2693_POST_P2680_NONSELECTOR_SOURCE_INVENTORY_CLOSED_AS_BOUNDED_NO_GO_REQUIRES_FRESH_STATE_MAP_NO_FALSE_PASS",
        "nonselector_bridge_source_gate_passes_now": lattice["nonselector_bridge_source_gate_passes_now"],
        "missing_nonselector_source_atoms": lattice["missing_current_obligations"],
        "hamming_distance_to_nonselector_bridge_source_pass": lattice["hamming_distance_to_nonselector_bridge_source_pass"],
        "live_lanes": live_lanes,
        "professorial_verdict": (
            "P2693 closes the current P2680 non-selector source-inventory round without pretending to solve the bridge.  The finite lattice has three missing source obligations and only the all-true state would pass; the current bit vector is all false for canonical UV unit, role-safe amplitude source, and target-independent positive beta/Z_beta source.  Prior audits supply meaningful near-misses, but every named non-selector atom is bounded no-go on current artifacts.  Thus generic bridge completion, role transfer, selector replay, L_total promotion, and ToE closure remain inadmissible."
        ),
        "next_honest_step": (
            "P2694 should run a fresh broad state-map selection audit across non-bridge lanes and admit a next proof-grade move only if it introduces a genuinely new typed object, theorem, source atom, blocker-cut, or provider class; otherwise the correct output is a no-new-live-frontier certificate rather than replay."
        ),
        "bridge_completion_claimed_now": False,
        "role_transfer_started_now": False,
        "ltotal_promoted_now": False,
        "toe_closed_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = ["# P2693/S1643 post-P2680 non-selector source-inventory closure/state-map", "", f"Status: `{payload['status']}`", "", "## Content-first grep"]
    for name, data in payload["content_grep"]["patterns"].items():
        lines.append(f"- `{name}`: `{data['count']}` hits")
    lines.extend(["", "## Non-selector atom status"])
    for row in payload["nonselector_atom_status"]:
        lines.append(f"- `{row['atom']}`: exported=`{row['source_exported_now']}`, bounded_no_go=`{row['bounded_no_go_now']}` — {row['blocking_reason']}")
    lat = payload["bridge_source_lattice"]
    lines.extend(["", "## Bridge-source lattice", f"Total states: `{lat['total_states']}`; passing states: `{lat['passing_states']}`; current hamming distance to pass: `{lat['hamming_distance_to_nonselector_bridge_source_pass']}`."])
    lines.extend(["", "## Lane reconciliation"])
    for row in payload["lane_reconciliation"]:
        lines.append(f"- `{row['lane']}`: live_now=`{row['live_now']}` — {row['reason']}")
    lines.extend(["", "## Verdict", payload["decision"]["professorial_verdict"], "", "## Next honest step", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    reads = state_reads()
    atom_rows = nonselector_atom_status(reads)
    lattice = bridge_source_lattice(atom_rows)
    lanes = lane_reconciliation(atom_rows)
    payload: dict[str, Any] = {
        "status": "P2693_POST_P2680_NONSELECTOR_SOURCE_INVENTORY_CLOSURE_STATE_MAP_NO_FALSE_PASS",
        "content_grep": content_grep(),
        "state_reads": reads,
        "nonselector_atom_status": atom_rows,
        "bridge_source_lattice": lattice,
        "lane_reconciliation": lanes,
        "decision": decision(lattice, lanes),
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2693/S1643 post-P2680 non-selector source-inventory closure/state-map",
        "## P2693/S1643 post-P2680 non-selector source-inventory closure/state-map\n\n"
        "`P2693/S1643` reconciles the P2680 non-selector source inventory after P2689-P2692.  The finite bridge-source lattice has three currently missing obligations: canonical length/UV unit source, role-safe amplitude source, and target-independent positive `beta/Z_beta` source.  All three are bounded no-go on current artifacts despite real near-misses.  Therefore the generic legacy-to-strict bridge lane is closed as a replay path until a genuinely new typed object/theorem/source atom appears; no bridge completion, role transfer, selector replay, role-bearing `L_total`, or ToE closure is exported.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2693/S1643 post-P2680 source-inventory Ltotal guard",
        "## P2693/S1643 post-P2680 source-inventory Ltotal guard\n\n"
        "`P2693/S1643` keeps `L_total` nonpromoted because every currently named P2680 non-selector source atom remains unexported.  Variational promotion requires a new source theorem, not a replay of bridge-source near-misses.\n",
    )
    append_once(
        AGENTS,
        "Current post-P2680 non-selector source-inventory closure guardrail (P2693/S1643, 2026-06-13)",
        "## Current post-P2680 non-selector source-inventory closure guardrail (P2693/S1643, 2026-06-13)\n\n"
        "- P2693 closes the currently named P2680 non-selector source atoms as bounded no-go on current artifacts: canonical UV unit, role-safe `alpha_geo` amplitude source, and target-independent positive `beta/Z_beta` source all remain unexported.\n"
        "- Do not continue generic legacy->strict bridge-source replay.  The next move must be a fresh broad state-map selection audit and may reopen a lane only with a genuinely new typed object, theorem, source atom, blocker-cut, or provider class; otherwise emit a no-new-live-frontier certificate.\n",
    )
    return payload


if __name__ == "__main__":
    main()
